# Copyright (c) SHRY Development Team.
# Distributed under the terms of the MIT License.

"""
CIF I/O using PyCifRW instead of pymatgen.io.cif
Supports magnetic CIF files for future development.
"""

import contextlib
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


def _sanitize_block_name(name: str) -> str:
    """Make a safe CIF data block name (letters/digits/underscore only)."""
    s = re.sub(r"\s+", "", str(name))
    s = re.sub(r"[^0-9A-Za-z_]", "_", s)
    s = s.strip("_")
    return s or "shry"


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
    # For magCIF, coordinate generation should use the magnetic symmetry operations
    # (time reversal does not affect coordinates). Keep the parsing logic separate:
    # - non-magnetic CIF: get_symops_from_cif
    # - magCIF (MSG): get_magnetic_symops_from_cif
    if "_space_group_symop_magn_operation.xyz" in cif_dict or "_space_group_symop_magn_operation_xyz" in cif_dict:
        symops = [op for op, _tr in _get_magnetic_symops_from_cif(cif_dict)]
    else:
        symops = _get_symops_from_cif(cif_dict)

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


def _get_cif_dict(cif_path: str) -> dict:
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


def _get_cif_block(cif_path: str, block_index: int = 0):
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


def _get_symops_from_cif(cif_dict: dict) -> List[SymmOp]:
    """
    Extract symmetry operations from CIF dictionary.
    Compatible with pymatgen CifParser.get_symops() behavior.

    Note:
        This function intentionally handles *non-magnetic* CIF symmetry operations.
        For magCIF (MSG) operations (including time reversal), use
        :func:`get_magnetic_symops_from_cif`.

    Args:
        cif_dict: Dictionary from get_cif_dict() for a specific block

    Returns:
        List of SymmOp objects
    """
    symops = []

    def _as_list(v):
        if isinstance(v, (list, tuple)):
            return list(v)
        return [v]

    def _clean_xyz(xyz: str) -> str:
        """Clean CIF symop strings."""
        return str(xyz).strip().replace("'", "").replace('"', "")

    # Try to get symmetry operations from CIF (non-magnetic).
    # Standard CIF tags:
    #   - _symmetry_equiv_pos_as_xyz
    #   - _space_group_symop_operation_xyz (and some variants)
    ops_key = None
    for key in [
        "_symmetry_equiv_pos_as_xyz",
        "_space_group_symop_operation_xyz",
        "_space_group_symop_operation.xyz",
    ]:
        if key in cif_dict:
            ops_key = key
            break

    if ops_key is None:
        return [SymmOp.from_xyz_str("x,y,z")]

    ops = []
    for s in _as_list(cif_dict[ops_key]):
        s = _clean_xyz(s)
        try:
            ops.append(SymmOp.from_xyz_str(s))
        except Exception as e:
            warnings.warn(f"Could not parse symmetry operation: {s}, error: {e}")

    symops.extend(ops)

    if not symops:
        return [SymmOp.from_xyz_str("x,y,z")]

    return symops


def _get_magnetic_symops_from_cif(cif_dict: dict) -> list[tuple[SymmOp, int]]:
    """Extract magnetic symmetry operations (including time reversal) from a magCIF dict.

    magCIF operations may include a trailing time-reversal term, e.g. "x,y,z,+1".
    Returns a list of (SymmOp, time_reversal) pairs, where time_reversal is +1 or -1.

    Centering operations (if present) are composed similarly, with time reversal multiplied.
    """

    def _as_list(v):
        if isinstance(v, (list, tuple)):
            return list(v)
        return [v]

    def _parse_op_and_tr(xyz: str) -> tuple[SymmOp, int] | None:
        s = str(xyz).strip().replace("'", "").replace('"', "")
        parts = [p.strip() for p in s.split(",") if p.strip()]
        if len(parts) < 3:
            return None
        tr = 1
        if len(parts) >= 4:
            try:
                tr = int(str(parts[3]).replace("+", ""))
                tr = 1 if tr >= 0 else -1
            except Exception:
                tr = 1
        try:
            op = SymmOp.from_xyz_str(",".join(parts[:3]))
        except Exception:
            return None
        return op, tr

    ops_key = None
    for key in [
        "_space_group_symop_magn_operation.xyz",
        "_space_group_symop_magn_operation_xyz",
    ]:
        if key in cif_dict:
            ops_key = key
            break

    if ops_key is None:
        # Fallback to non-magnetic symops (time reversal +1)
        return [(op, 1) for op in _get_symops_from_cif(cif_dict)]

    ops: list[tuple[SymmOp, int]] = []
    for raw in _as_list(cif_dict[ops_key]):
        parsed = _parse_op_and_tr(raw)
        if parsed is not None:
            ops.append(parsed)

    centers_key = None
    for key in [
        "_space_group_symop_magn_centering.xyz",
        "_space_group_symop_magn_centering_xyz",
    ]:
        if key in cif_dict:
            centers_key = key
            break

    if centers_key is None:
        return ops or [(SymmOp.from_xyz_str("x,y,z"), 1)]

    centers: list[tuple[SymmOp, int]] = []
    for raw in _as_list(cif_dict[centers_key]):
        parsed = _parse_op_and_tr(raw)
        if parsed is not None:
            centers.append(parsed)

    if not centers:
        return ops or [(SymmOp.from_xyz_str("x,y,z"), 1)]

    composed: list[tuple[SymmOp, int]] = []
    seen: set[tuple] = set()
    for op, tr_op in ops or [(SymmOp.from_xyz_str("x,y,z"), 1)]:
        for cen, tr_cen in centers:
            c = op * cen
            tr = int(tr_op) * int(tr_cen)
            tr = 1 if tr >= 0 else -1
            key = (tuple(np.round(c.affine_matrix, decimals=12).flatten()), tr)
            if key in seen:
                continue
            seen.add(key)
            composed.append((c, tr))
    return composed


def _structure_to_cif_string(structure: Structure, symprec: Optional[float] = None) -> str:
    """
    Convert pymatgen Structure to CIF string using PyCifRW.

    Args:
        structure: pymatgen Structure object
        symprec: Symmetry precision (not used yet, for future implementation)

    Returns:
        CIF file content as string
    """
    # NOTE: We write CIF manually (classic loop_ form) to avoid PyCifRW emitting
    # CIF2 list syntax like "_atom_site_label [ ... ]" which some tools (e.g. VESTA)
    # cannot parse reliably.

    del symprec  # reserved for future

    block_name = _sanitize_block_name(structure.formula)

    def _label_to_symbol(lbl: str) -> str:
        m = re.match(r"([A-Za-z]+)", str(lbl).strip())
        return m.group(1) if m else str(lbl).strip()

    rows: list[tuple[str, str, float, float, float, float]] = []
    element_count: dict[str, int] = {}

    def _next_label(sym: str) -> str:
        """Generate CIF label as <symbol><index>, e.g. Cu1, Cu2, ..."""
        sym = str(sym).strip() or "X"
        element_count[sym] = element_count.get(sym, 0) + 1
        return f"{sym}{element_count[sym]}"

    for site in structure.sites:
        occ_map = site.properties.get("_atom_site_occupancy_by_label") or {}
        label_prop = site.properties.get("_atom_site_label")
        if isinstance(label_prop, str):
            site_labels = [label_prop]
        elif isinstance(label_prop, (list, tuple)):
            site_labels = [x for x in label_prop if isinstance(x, str)]
        else:
            site_labels = []

        x, y, z = (float(site.frac_coords[0]), float(site.frac_coords[1]), float(site.frac_coords[2]))

        # Expand by per-label occupancy when available.
        if occ_map and len(occ_map) > 1:
            for lbl in site_labels or list(occ_map.keys()):
                occ = float(occ_map.get(lbl, 0.0))
                if occ <= 0.0:
                    continue
                sym = _label_to_symbol(lbl)
                rows.append((_next_label(sym), sym, x, y, z, occ))
            continue

        # Otherwise expand by pymatgen's species for disordered sites.
        if not site.is_ordered:
            try:
                items = list(site.species.items())
            except Exception:
                items = []
            if items:
                for sp, occ in items:
                    occ_f = float(occ)
                    if occ_f <= 0.0:
                        continue
                    sym = getattr(sp, "symbol", str(sp))
                    rows.append((_next_label(sym), sym, x, y, z, occ_f))
                continue

        # Ordered site.
        sym = site.specie.symbol if hasattr(site, "specie") else site.species_string
        occ = float(site.properties.get("occupancy", 1.0))
        rows.append((_next_label(sym), sym, x, y, z, occ))

    lines: list[str] = []
    lines.append(f"data_{block_name}")
    lat = structure.lattice
    lines.append(f"_cell_length_a    {lat.a:.6f}")
    lines.append(f"_cell_length_b    {lat.b:.6f}")
    lines.append(f"_cell_length_c    {lat.c:.6f}")
    lines.append(f"_cell_angle_alpha {lat.alpha:.6f}")
    lines.append(f"_cell_angle_beta  {lat.beta:.6f}")
    lines.append(f"_cell_angle_gamma {lat.gamma:.6f}")
    lines.append("_symmetry_space_group_name_H-M 'P 1'")
    lines.append("_symmetry_Int_Tables_number 1")
    lines.append("loop_")
    lines.append("  _symmetry_equiv_pos_as_xyz")
    lines.append("  'x,y,z'")
    lines.append("loop_")
    lines.append("  _atom_site_label")
    lines.append("  _atom_site_type_symbol")
    lines.append("  _atom_site_fract_x")
    lines.append("  _atom_site_fract_y")
    lines.append("  _atom_site_fract_z")
    lines.append("  _atom_site_occupancy")
    for lbl, sym, x, y, z, occ in rows:
        lines.append(f"  {lbl} {sym} {x:.6f} {y:.6f} {z:.6f} {occ:.6f}")
    return "\n".join(lines) + "\n"


def _structure_to_mcif_string(structure: Structure, symprec: Optional[float] = None) -> str:
    """Convert a (possibly magnetic) pymatgen Structure to a minimal magCIF string.

    - Preserves magnetic moments when present via a magCIF loop:
      _atom_site_moment.label + crystalaxis components.
    - For disordered/mixed sites, emits one _atom_site_ row per label using
      site.properties['_atom_site_occupancy_by_label'] when available.

    Notes:
        This writer intentionally outputs P1 symmetry for simplicity. SHRY's
        enumeration and symmetry handling is performed via spglib; the output
        here focuses on retaining magnetic moment data.
    """

    # Manual magCIF writer using classic loop_ sections.
    del symprec  # reserved for future

    block_name = _sanitize_block_name(structure.formula)

    def _label_to_symbol(lbl: str) -> str:
        m = re.match(r"([A-Za-z]+)", str(lbl).strip())
        return m.group(1) if m else str(lbl).strip()

    site_rows: list[tuple[str, str, float, float, float, float]] = []
    moment_rows: list[tuple[str, float, float, float]] = []
    element_count: dict[str, int] = {}

    def _next_label(sym: str) -> str:
        """Generate magCIF label as <symbol><index>, e.g. Cu1, Cu2, ..."""
        sym = str(sym).strip() or "X"
        element_count[sym] = element_count.get(sym, 0) + 1
        return f"{sym}{element_count[sym]}"

    for site in structure.sites:
        occ_map = site.properties.get("_atom_site_occupancy_by_label") or {}
        label_prop = site.properties.get("_atom_site_label")
        if isinstance(label_prop, str):
            site_labels = [label_prop]
        elif isinstance(label_prop, (list, tuple)):
            site_labels = [x for x in label_prop if isinstance(x, str)]
        else:
            site_labels = []

        mom_map = site.properties.get("_atom_site_moment_by_label") or {}
        x, y, z = (float(site.frac_coords[0]), float(site.frac_coords[1]), float(site.frac_coords[2]))

        # Expand by per-label occupancy (preferred).
        if occ_map and len(occ_map) > 1:
            for lbl in site_labels or list(occ_map.keys()):
                occ = float(occ_map.get(lbl, 0.0))
                if occ <= 0.0:
                    continue
                sym = _label_to_symbol(lbl)
                out_lbl = _next_label(sym)
                site_rows.append((out_lbl, sym, x, y, z, occ))
                vec = mom_map.get(lbl)
                if vec is not None:
                    v = np.array(vec, dtype=float)
                    if v.shape == (3,) and not np.allclose(v, 0.0):
                        moment_rows.append((out_lbl, float(v[0]), float(v[1]), float(v[2])))
            continue

        # Otherwise expand by species if disordered.
        if not site.is_ordered:
            try:
                items = list(site.species.items())
            except Exception:
                items = []
            if items:
                for sp, occ in items:
                    occ_f = float(occ)
                    if occ_f <= 0.0:
                        continue
                    sym = getattr(sp, "symbol", str(sp))
                    lbl = _next_label(sym)
                    site_rows.append((lbl, sym, x, y, z, occ_f))
                continue

        # Ordered.
        sym = site.specie.symbol if hasattr(site, "specie") else site.species_string
        base_lbl = site_labels[0] if len(site_labels) == 1 else sym
        lbl = _next_label(sym)
        occ = float(site.properties.get("occupancy", 1.0))
        site_rows.append((str(lbl), sym, x, y, z, occ))

        vec = mom_map.get(base_lbl)
        if vec is None:
            vec = site.properties.get("magmom")
        if vec is not None:
            v = np.array(vec, dtype=float)
            if v.shape == ():
                v = np.array([float(v), 0.0, 0.0], dtype=float)
            if v.shape == (3,) and not np.allclose(v, 0.0):
                moment_rows.append((str(lbl), float(v[0]), float(v[1]), float(v[2])))

    lines: list[str] = []
    # VESTA supports CIF1.1 and magCIF 1.0 (CIF1.1-based). Avoid CIF2.0 headers.
    lines.append(f"data_{block_name}")
    lat = structure.lattice
    lines.append(f"_cell_length_a    {lat.a:.6f}")
    lines.append(f"_cell_length_b    {lat.b:.6f}")
    lines.append(f"_cell_length_c    {lat.c:.6f}")
    lines.append(f"_cell_angle_alpha {lat.alpha:.6f}")
    lines.append(f"_cell_angle_beta  {lat.beta:.6f}")
    lines.append(f"_cell_angle_gamma {lat.gamma:.6f}")
    lines.append("_symmetry_space_group_name_H-M 'P 1'")
    lines.append("_symmetry_Int_Tables_number 1")
    lines.append("loop_")
    lines.append("  _symmetry_equiv_pos_as_xyz")
    lines.append("  'x,y,z'")
    lines.append("loop_")
    lines.append("  _atom_site_label")
    lines.append("  _atom_site_type_symbol")
    lines.append("  _atom_site_fract_x")
    lines.append("  _atom_site_fract_y")
    lines.append("  _atom_site_fract_z")
    lines.append("  _atom_site_occupancy")
    for lbl, sym, x, y, z, occ in site_rows:
        lines.append(f"  {lbl} {sym} {x:.6f} {y:.6f} {z:.6f} {occ:.6f}")

    if moment_rows:
        lines.append("loop_")
        lines.append("  _atom_site_moment.label")
        lines.append("  _atom_site_moment.crystalaxis_x")
        lines.append("  _atom_site_moment.crystalaxis_y")
        lines.append("  _atom_site_moment.crystalaxis_z")
        for lbl, mx, my, mz in moment_rows:
            lines.append(f"  {lbl} {mx:.6f} {my:.6f} {mz:.6f}")

    return "\n".join(lines) + "\n"


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

    # Common CIF missing values
    if value_str in {"?", ".", ""}:
        return 0.0

    # Fractions occasionally appear (e.g. symmetry-related values), handle "3/4".
    if "/" in value_str:
        parts = value_str.split("/")
        if len(parts) == 2:
            try:
                num = float(parts[0])
                den = float(parts[1])
                if den != 0:
                    return num / den
            except ValueError:
                pass

    # Try direct float first.
    try:
        return float(value_str)
    except ValueError:
        pass

    # Bilbao magCIF sometimes uses tokens like "-0.5." (note the extra trailing dot).
    # More generally, extract the leading numeric token and parse that.
    m = re.match(r"^[+-]?(?:\d+\.\d*|\d*\.\d+|\d+)(?:[eE][+-]?\d+)?", value_str)
    if m:
        try:
            return float(m.group(0))
        except ValueError:
            return 0.0

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
            oxi_states[str(symbol).strip()] = _str2float(oxi)
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


def _str2float(text):
    """
    Remove uncertainty brackets from strings and return float.
    Compatible with pymatgen's str2float function.

    Args:
        text: String potentially containing uncertainty notation

    Returns:
        Floating point number
    """
    return _parse_cif_float(text)
