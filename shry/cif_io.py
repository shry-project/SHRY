# Copyright (c) SHRY Development Team.
# Distributed under the terms of the MIT License.

"""
CIF I/O using PyCifRW instead of pymatgen.io.cif
Supports magnetic CIF files for future development.
"""

import io
import re
import warnings
from typing import Dict, List, Optional, Tuple

import numpy as np
from CifFile import CifFile as PyCifFile, ReadCif
from pymatgen.core import Element, Lattice, PeriodicSite, Species, Structure
from pymatgen.core.composition import Composition
from pymatgen.core.operations import SymmOp
from pymatgen.util.coord import in_coord_list_pbc


def parse_cif_to_structure(
    cif_path: str,
    primitive: bool = False,
    sort: bool = False,
    merge_tol: float = 1e-3,
    keep_oxidation_states: bool = False,
) -> Structure:
    """
    Parse a CIF file using PyCifRW and convert to pymatgen Structure.

    Args:
        cif_path: Path to CIF file
        primitive: Whether to convert to primitive cell
        sort: Whether to sort sites
        merge_tol: Tolerance for merging sites
        keep_oxidation_states: Whether to keep oxidation states when present

    Returns:
        pymatgen.Structure object
    """
    # Use permissive mode to handle CIF format errors
    cf = ReadCif(cif_path, permissive=True)

    # Get first data block
    block_name = list(cf.keys())[0]
    block = cf[block_name]

    return _parse_cif_block_to_structure(
        block,
        primitive=primitive,
        sort=sort,
        merge_tol=merge_tol,
        keep_oxidation_states=keep_oxidation_states,
    )


def parse_cif_string_to_structure(
    cif_text: str,
    primitive: bool = False,
    sort: bool = False,
    merge_tol: float = 1e-3,
    keep_oxidation_states: bool = False,
) -> Structure:
    """
    Parse a CIF string using PyCifRW and convert to pymatgen Structure.

    Args:
        cif_text: CIF contents as a string
        primitive: Whether to convert to primitive cell
        sort: Whether to sort sites
        merge_tol: Tolerance for merging sites
        keep_oxidation_states: Whether to keep oxidation states when present

    Returns:
        pymatgen.Structure object
    """
    cf = ReadCif(io.StringIO(cif_text), permissive=True)
    block_name = list(cf.keys())[0]
    block = cf[block_name]

    return _parse_cif_block_to_structure(
        block,
        primitive=primitive,
        sort=sort,
        merge_tol=merge_tol,
        keep_oxidation_states=keep_oxidation_states,
    )


def _parse_cif_block_to_structure(
    block,
    primitive: bool,
    sort: bool,
    merge_tol: float,
    keep_oxidation_states: bool,
) -> Structure:
    """
    Parse a single CIF block into a pymatgen Structure.
    """
    # Parse lattice parameters
    a = _get_cif_value(block, "_cell_length_a")
    b = _get_cif_value(block, "_cell_length_b")
    c = _get_cif_value(block, "_cell_length_c")
    alpha = _get_cif_value(block, "_cell_angle_alpha")
    beta = _get_cif_value(block, "_cell_angle_beta")
    gamma = _get_cif_value(block, "_cell_angle_gamma")

    lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)

    # Get symmetry operations
    cif_dict = {k: block.get(k) for k in block.keys()}
    symops = get_symops_from_cif(cif_dict)

    # Parse atomic sites (asymmetric unit)
    asym_sites = []

    # Get site data
    labels = _get_loop_values(block, "_atom_site_label")
    type_symbols = _get_loop_values(block, "_atom_site_type_symbol")
    fract_x = _get_loop_values(block, "_atom_site_fract_x")
    fract_y = _get_loop_values(block, "_atom_site_fract_y")
    fract_z = _get_loop_values(block, "_atom_site_fract_z")

    # Optional: occupancy
    try:
        occupancies = _get_loop_values(block, "_atom_site_occupancy")
    except KeyError:
        occupancies = [1.0] * len(labels)

    oxi_states = _get_oxidation_states(block) if keep_oxidation_states else None

    # Group sites by coordinates to handle partial occupancy.
    # NOTE: Use a rounded coordinate tuple only as a dictionary key; keep the original
    # (unrounded) coordinates for subsequent symmetry detection.
    coord_sites = {}  # {coord_key: {"coords": np.ndarray, "species_occs": [(species, occ), ...]}}

    for i, label in enumerate(labels):
        # Parse species
        symbol = type_symbols[i] if type_symbols else label
        species = _parse_species(symbol, oxi_states=oxi_states)
        if species is None:
            continue

        coords = np.array(
            [
                _parse_cif_float(fract_x[i]),
                _parse_cif_float(fract_y[i]),
                _parse_cif_float(fract_z[i]),
            ]
        )

        occ = _parse_cif_float(occupancies[i])

        # Use tuple of coordinates as key (rounded to 4 decimals to match tolerance)
        coord_key = tuple(np.round(coords, decimals=4))

        if coord_key not in coord_sites:
            coord_sites[coord_key] = {"coords": coords, "species_occs": []}
        coord_sites[coord_key]["species_occs"].append((species, occ))

    # Create asymmetric unit sites with merged species (partial occupancy)
    for _coord_key, payload in coord_sites.items():
        coords = payload["coords"]
        species_occs = payload["species_occs"]

        if len(species_occs) == 1:
            # Single species at this site
            species, occ = species_occs[0]
            asym_sites.append((species, coords, occ))
        else:
            # Multiple species at same site (partial occupancy)
            comp_dict = {sp: occ for sp, occ in species_occs}
            composition = Composition(comp_dict)
            asym_sites.append((composition, coords, 1.0))

    # Apply symmetry operations to generate all sites in unit cell
    all_coords = []
    all_sites = []

    for species, coords, occ in asym_sites:
        for symop in symops:
            new_coords = symop.operate(coords)
            # Wrap to [0, 1)
            new_coords = new_coords % 1.0

            # Check if this site already exists using in_coord_list_pbc
            if not in_coord_list_pbc(all_coords, new_coords, atol=1e-4):
                all_coords.append(new_coords)
                all_sites.append(PeriodicSite(species, new_coords, lattice, properties={"occupancy": occ}))

    structure = Structure.from_sites(all_sites)

    # Post-processing
    if merge_tol > 1e-3:
        structure.merge_sites(merge_tol, mode="sum")

    if primitive:
        structure = structure.get_primitive_structure()

    if sort:
        structure = structure.get_sorted_structure()

    return structure


def get_cif_dict(cif_path: str) -> dict:
    """
    Get CIF data as dictionary (similar to pymatgen CifParser.as_dict()).

    Args:
        cif_path: Path to CIF file

    Returns:
        Dictionary with CIF block names as keys and CIF data as values
    """
    cf = ReadCif(cif_path, permissive=True)
    result = {}

    for block_name in cf.keys():
        block = cf[block_name]
        block_dict = {}

        # Convert all CIF items to dictionary
        for key in block.keys():
            value = block[key]
            # Store as-is (lists stay lists, single values stay single)
            block_dict[key] = value

        result[block_name] = block_dict

    return result


def get_cif_block(cif_path: str, block_index: int = 0):
    """
    Get a specific CIF block from file.

    Args:
        cif_path: Path to CIF file
        block_index: Index of block to return (default: 0 for first block)

    Returns:
        CifFile block object
    """
    cf = ReadCif(cif_path, permissive=True)
    block_names = list(cf.keys())
    if block_index >= len(block_names):
        raise IndexError(f"Block index {block_index} out of range (only {len(block_names)} blocks)")
    return cf[block_names[block_index]]


def get_symops_from_cif(cif_dict: dict) -> List[SymmOp]:
    """
    Extract symmetry operations from CIF dictionary.
    Compatible with pymatgen CifParser.get_symops() behavior.

    Args:
        cif_dict: Dictionary from get_cif_dict() for a specific block

    Returns:
        List of SymmOp objects
    """
    symops = []

    # Try to get symmetry operations from CIF
    # Look for _symmetry_equiv_pos_as_xyz or _space_group_symop_operation_xyz
    symm_key = None
    for key in ["_symmetry_equiv_pos_as_xyz", "_space_group_symop_operation_xyz"]:
        if key in cif_dict:
            symm_key = key
            break

    if symm_key is None:
        # No symmetry operations found, return identity
        return [SymmOp.from_xyz_str("x,y,z")]

    symm_strs = cif_dict[symm_key]
    if not isinstance(symm_strs, list):
        symm_strs = [symm_strs]

    for s in symm_strs:
        # Clean up the string
        s = s.strip().replace("'", "").replace('"', "")
        try:
            symops.append(SymmOp.from_xyz_str(s))
        except Exception as e:
            warnings.warn(f"Could not parse symmetry operation: {s}, error: {e}")

    if not symops:
        # If parsing failed, return identity
        return [SymmOp.from_xyz_str("x,y,z")]

    return symops


def structure_to_cif_string(structure: Structure, symprec: Optional[float] = None) -> str:
    """
    Convert pymatgen Structure to CIF string using PyCifRW.

    Args:
        structure: pymatgen Structure object
        symprec: Symmetry precision (not used yet, for future implementation)

    Returns:
        CIF file content as string
    """
    cf = PyCifFile()

    # Create data block
    block_name = structure.formula.replace(" ", "")
    cf.NewBlock(block_name)
    block = cf[block_name]

    # Add lattice parameters
    lattice = structure.lattice
    block["_cell_length_a"] = f"{lattice.a:.6f}"
    block["_cell_length_b"] = f"{lattice.b:.6f}"
    block["_cell_length_c"] = f"{lattice.c:.6f}"
    block["_cell_angle_alpha"] = f"{lattice.alpha:.6f}"
    block["_cell_angle_beta"] = f"{lattice.beta:.6f}"
    block["_cell_angle_gamma"] = f"{lattice.gamma:.6f}"

    # Add symmetry info (P1 for now)
    block["_symmetry_space_group_name_H-M"] = "P 1"
    block["_symmetry_Int_Tables_number"] = "1"

    # Add atomic sites
    labels = []
    type_symbols = []
    fract_xs = []
    fract_ys = []
    fract_zs = []
    occupancies = []

    element_count = {}
    for site in structure.sites:
        # Generate unique label
        species_string = site.species_string
        if species_string not in element_count:
            element_count[species_string] = 0
        element_count[species_string] += 1
        label = f"{species_string}{element_count[species_string]}"

        labels.append(label)
        type_symbols.append(species_string)
        fract_xs.append(f"{site.frac_coords[0]:.6f}")
        fract_ys.append(f"{site.frac_coords[1]:.6f}")
        fract_zs.append(f"{site.frac_coords[2]:.6f}")

        occ = site.properties.get("occupancy", 1.0)
        occupancies.append(f"{occ:.4f}")

    # Assign loop data directly (PyCifRW automatically creates loops from list values)
    block["_atom_site_label"] = labels
    block["_atom_site_type_symbol"] = type_symbols
    block["_atom_site_fract_x"] = fract_xs
    block["_atom_site_fract_y"] = fract_ys
    block["_atom_site_fract_z"] = fract_zs
    block["_atom_site_occupancy"] = occupancies

    return str(cf)


def _get_cif_value(block, key: str) -> float:
    """Get a single numeric value from CIF block."""
    value = block.get(key)
    if value is None:
        raise KeyError(f"Key {key} not found in CIF block")
    return _parse_cif_float(value)


def _get_loop_values(block, key: str) -> List:
    """Get loop values from CIF block."""
    value = block.get(key)
    if value is None:
        raise KeyError(f"Key {key} not found in CIF block")

    # Handle both single values and lists
    if isinstance(value, (list, tuple)):
        return list(value)
    else:
        return [value]


def _parse_cif_float(value) -> float:
    """Parse a CIF numeric value (may have uncertainty in parentheses)."""
    if isinstance(value, (int, float)):
        return float(value)

    # Remove uncertainty: 1.234(5) -> 1.234
    value_str = str(value).strip()
    value_str = re.sub(r"\([0-9]+\)", "", value_str)

    try:
        return float(value_str)
    except ValueError:
        return 0.0


def _get_oxidation_states(block) -> Optional[Dict[str, float]]:
    """
    Extract oxidation states from CIF block if available.
    """
    try:
        type_symbols = _get_loop_values(block, "_atom_type_symbol")
        oxi_numbers = _get_loop_values(block, "_atom_type_oxidation_number")
    except KeyError:
        return None

    oxi_states = {}
    for symbol, oxi in zip(type_symbols, oxi_numbers, strict=False):
        try:
            oxi_states[str(symbol).strip()] = str2float(oxi)
        except Exception:
            continue
    return oxi_states or None


def _parse_species(symbol: str, oxi_states: Optional[Dict[str, float]] = None):
    """
    Parse a CIF type symbol into a pymatgen Species, optionally keeping oxidation states.
    """
    symbol = str(symbol).strip()
    if not symbol:
        return None

    # If oxidation states are explicitly provided, keep them.
    if oxi_states and symbol in oxi_states:
        base_symbol = re.sub(r"[0-9+-]+", "", symbol)
        try:
            return Species(base_symbol, oxi_states[symbol])
        except Exception:
            pass

    # Otherwise, behave like typical CIF readers: treat this as a plain element.
    # Always strip charge / oxidation decorations and digits.
    cleaned = re.sub(r"[0-9+-]+", "", symbol)
    element_match = re.match(r"([A-Z][a-z]?)", cleaned)
    if element_match:
        try:
            return Element(element_match.group(1))
        except Exception:
            return None

    warnings.warn(f"Could not parse species: {symbol}, skipping")
    return None


def str2float(text):
    """
    Remove uncertainty brackets from strings and return float.
    Compatible with pymatgen's str2float function.

    Args:
        text: String potentially containing uncertainty notation

    Returns:
        Floating point number
    """
    return _parse_cif_float(text)
