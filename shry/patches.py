# Copyright (c) SHRY Development Team.
# Distributed under the terms of the MIT License.

"""
Runtime patches for pymatgen classes.
"""

from __future__ import annotations

import collections
import itertools
import re
from typing import Tuple

import numpy as np
from monty.fractions import gcd_float
from pymatgen.core import Composition, Species
from pymatgen.core.composition import CompositionError, reduce_formula
from pymatgen.core.periodic_table import get_el_sp
from pymatgen.io.cif import CifParser, str2float
from pymatgen.io.vasp.inputs import Poscar

from . import const

_APPLIED = False

# Patched extra functionalities and bug fixes on top of Pymatgen's classes.

def get_integer_formula_and_factor(
    self,
    max_denominator: int = int(1 / const.DEFAULT_ATOL),
    iupac_ordering: bool = False,
) -> Tuple[str, float]:
    """The default Composition groups together different ox states which is not ideal..."""
    el_amt = self.as_dict()
    g = gcd_float(list(el_amt.values()), 1 / max_denominator)

    d = {k: round(v / g) for k, v in el_amt.items()}
    (formula, factor) = reduce_formula(d, iupac_ordering=iupac_ordering)
    if formula in Composition.special_formulas:
        formula = Composition.special_formulas[formula]
        factor /= 2
    return formula, factor * g


def to_int_dict(self):
    """
    Returns:
        Dict with element symbol and integer amount
    """
    _, factor = self.get_integer_formula_and_factor(max_denominator=int(1 / const.DEFAULT_ATOL))
    int_dict = {e: int(a) for e, a in (self / factor).as_dict().items()}

    # be safe: Composition groups together different ox states which is not ideal...
    for x, y in zip(int_dict.values(), self.as_dict().values()):
        if not np.isclose(x * factor, y, atol=const.DEFAULT_ATOL):
            raise ValueError(
                "Composition (Occupancy) is not rational! Please try to increase significant digits "
                "e.g., 1/3 = 0.3333 -> 1/3 = 0.3333333333333."
            )

    return int_dict


@property
def inted_composition(self):
    """
    Return Composition instance with integer formula
    """
    _, factor = self.get_integer_formula_and_factor(max_denominator=int(1 / const.DEFAULT_ATOL))
    int_comp = self / factor

    # be safe
    int_dict = {e: int(a) for e, a in int_comp.as_dict().items()}
    if not all(np.isclose(x * factor, y, atol=const.DEFAULT_ATOL) for x, y in zip(int_dict.values(), self.as_dict().values())):
        raise ValueError(
            "Composition (Occupancy) is not rational! Please try to increase significant digits "
            "e.g., 1/3 = 0.3333 -> 1/3 = 0.3333333333333."
        )

    return int_comp


def formula_double_format_tol(afloat, ignore_ones=True, tol: float = const.DEFAULT_ATOL * 10):
    """
    This function is used to make pretty formulas by formatting the amounts.
    Instead of Li1.0 Fe1.0 P1.0 O4.0, you get LiFePO4.

    Args:
        afloat (float): a float
        ignore_ones (bool): if true, floats of 1 are ignored.
        tol (float): Tolerance to round to nearest int. i.e. 2.0000000001 -> 2

    Returns:
        A string representation of the float for formulas.
    """
    if ignore_ones and afloat == 1:
        return ""
    if abs(afloat - round(afloat)) < tol:
        return round(afloat)
    return round(afloat, 8)


@property
def formula(self) -> str:
    """
    Returns a formula string, with elements sorted by electronegativity,
    e.g., Li4 Fe4 P4 O16.
    """
    sym_amt = self.get_el_amt_dict()
    syms = sorted(sym_amt, key=lambda sym: get_el_sp(sym).X)
    formula = [f"{s}{formula_double_format_tol(sym_amt[s], False)}" for s in syms]
    return " ".join(formula)


def from_string(composition_string) -> Composition:
    """
    Workaround when working with strings including oxication states
    """

    def composition_builder():
        components = re.findall(r"[A-z][a-z]*[0-9.\+\-]*[0-9.]*", composition_string)
        amount = re.compile(r"(?<=[\+\-A-Za-z])[0-9.]+(?![\+\-])")
        amount_part = [amount.findall(x) for x in components]
        amount_part = [x[0] if x else "1" for x in amount_part]
        species_part = [x.strip(amount) for x, amount in zip(components, amount_part)]
        species_part = [x + "0" if not re.search(r"[0-9\+\-]+", x) else x for x in species_part]
        amount_part = [float(x) for x in amount_part]

        return Composition({Species.from_string(species): amount for species, amount in zip(species_part, amount_part)})

    try:
        return Composition(composition_string)
    except (CompositionError, ValueError, IndexError):
        # - CompositionError: get_sym_dict() error
        #   when "+" oxidation is used.
        # - ValueError: Wrong parse by get_sym_dict()
        #   if the oxidation negative.
        # - IndexError: Sometimes appear when the string gets more complex
        return composition_builder()


@property
def site_symbols(self):
    """
    On Poscar: sometimes we would like to use a separate pseudopotential
    for each oxidation states.

    Write oxidation states, if, and only if, for the given element,
    there are more than 1 state.
    """
    return [a[0] for a in itertools.groupby(self.syms)]


@property
def syms(self):
    """
    Replaced self.syms in some functions of Poscar
    """
    e_oxstates = collections.defaultdict(set)
    s_symbols = collections.defaultdict(list)

    for site in self.structure:
        e_oxstates[site.specie.symbol].add(str(site.specie))
    for e, oes in e_oxstates.items():
        if len(oes) == 1:
            s_symbols[list(oes)[0]] = e
        else:
            for oe in oes:
                s_symbols[oe] = oe

    return [s_symbols[str(site.specie)] for site in self.structure]


@property
def natoms(self):
    """
    Sequence of number of sites of each type associated with the Poscar.
    Similar to 7th line in vasp 5+ POSCAR or the 6th line in vasp 4 POSCAR.
    """
    return [len(tuple(a[1])) for a in itertools.groupby(self.syms)]


@staticmethod
def parse_oxi_states(data):
    """
    Parse oxidation states from data dictionary
    """
    ox_state_regex = re.compile(r"\d?[\+,\-]?$")

    try:
        oxi_states = {
            data["_atom_type_symbol"][i]: str2float(data["_atom_type_oxidation_number"][i])
            for i in range(len(data["_atom_type_symbol"]))
        }
        # attempt to strip oxidation state from _atom_type_symbol
        # in case the label does not contain an oxidation state
        for i, symbol in enumerate(data["_atom_type_symbol"]):
            oxi_states[ox_state_regex.sub("", symbol)] = str2float(data["_atom_type_oxidation_number"][i])
    except (ValueError, KeyError):
        # Some CIF (including pymatgen's output) are like this.
        try:
            oxi_states = dict()
            for i, symbol in enumerate(data["_atom_site_type_symbol"]):
                _symbol = ox_state_regex.sub("", symbol)

                parsed_oxi_state = ox_state_regex.search(symbol).group(0)
                if not parsed_oxi_state:
                    oxi_states[_symbol] = None
                    continue

                sign = re.search(r"[-+]", parsed_oxi_state)
                if sign is None:
                    sign = ""
                else:
                    sign = sign.group(0)
                parsed_oxi_state = parsed_oxi_state.replace("+", "").replace("-", "")
                parsed_oxi_state = str2float(sign + parsed_oxi_state)
                oxi_states[ox_state_regex.sub("", symbol)] = parsed_oxi_state
        except (ValueError, KeyError):
            oxi_states = None
    return oxi_states


def apply_pymatgen_patches() -> None:
    """
    Apply runtime patches to pymatgen classes.
    """
    global _APPLIED
    if _APPLIED:
        return

    Composition.to_int_dict = to_int_dict
    Composition.get_integer_formula_and_factor = get_integer_formula_and_factor
    Composition.inted_composition = inted_composition
    Composition.formula = formula
    Composition.from_string = from_string

    Poscar.site_symbols = site_symbols
    Poscar.natoms = natoms
    Poscar.syms = syms

    CifParser.parse_oxi_states = parse_oxi_states

    _APPLIED = True
