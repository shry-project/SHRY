# Copyright (c) SHRY Development Team.
# Distributed under the terms of the MIT License.
"""Test for SHRY structure generation.

1. Check count consistency: number of generated structures should match Substitutor.count()
2. Check benchmark consistency: count matches benchmark_SG_all.xls
3. Check internal redundancy: all generated structures should be symmetrically inequivalent

Configuration:
    Edit the constants at the top of this file to change test behavior:
    - BENCHMARK_COUNT: Number of benchmark cases to test (default: 10)
    - BENCHMARK_FILE: Path to benchmark Excel file
    - SORT_BY_COUNT: Sort by structure count (smallest first, default: True)
    - MAX_TEST_CASES: Maximum test cases for count/benchmark tests (default: 500)
    - MAX_REDUNDANCY_TEST_CASES: Maximum for redundancy tests (default: 100)
    - N_JOBS: Number of parallel jobs for redundancy checking (default: -1, use all cores)

Usage examples:
    # Run all tests except slow redundancy check (default configuration)
    pytest tests/test_redundancy_and_count.py -v -k "not no_redundancy"

    # Test only count consistency
    pytest tests/test_redundancy_and_count.py::test_count_consistency -v

    # Test only benchmark consistency
    pytest tests/test_redundancy_and_count.py::test_benchmark_count -v

    # Test redundancy (VERY SLOW)
    pytest tests/test_redundancy_and_count.py::test_no_redundancy -v

    # Skip slow tests marked with @pytest.mark.slow
    pytest tests/test_redundancy_and_count.py -v -m "not slow"

    # To change number of test cases, edit BENCHMARK_COUNT in this file
"""

import contextlib
import itertools
import os
import numpy as np
import pandas as pd
import pytest
import tqdm
from joblib import Parallel, delayed
from pymatgen.io.cif import CifFile, str2float
from pymatgen.util.coord import find_in_coord_list_pbc

from shry import LabeledStructure
from shry.core import NeedSupercellError, Substitutor, get_symmetry_operations_from_spglib

import warnings

warnings.simplefilter("ignore")

pytestmark = [
    pytest.mark.filterwarnings(
        "ignore:No Pauling electronegativity for .*:UserWarning",
    ),
    pytest.mark.filterwarnings(
        r"ignore:Set OLD_ERROR_HANDLING to false and catch the errors directly\.:DeprecationWarning",
    ),
]

# TOLERANCE
SHRY_TOLERANCE = 1.0e-2  # angstrom
SHRY_ANGLE_TOLERANCE = 5.0  # degree
SHRY_ATOL = 1.0e-2

# Display formatting
LINEWIDTH = 88
HLINE = "-" * LINEWIDTH
TQDM_CONF = {"ncols": LINEWIDTH}
DISABLE_PROGRESSBAR = False


# ============================================================================
# Test configuration (edit these values to change test behavior)
# ============================================================================

# Number of benchmark cases to test
BENCHMARK_COUNT = 10

# Path to benchmark Excel file (relative to repository root)
BENCHMARK_FILE = "tests/test_cifs/all_space_groups/benchmark_SG_all.xls"

# Sort by number of equivalent structures (smallest first)
SORT_BY_COUNT = True

# Maximum number of test cases for each test function
MAX_TEST_CASES = 10
MAX_REDUNDANCY_TEST_CASES = 1  # Redundancy test is very slow, so limit to 100

# Number of parallel jobs for redundancy checking (-1 = use all CPU cores)
N_JOBS = -1


# ============================================================================
# Pytest fixtures
# ============================================================================


@pytest.fixture(scope="module")
def benchmark_data(request):
    """Load benchmark data from Excel file."""
    df = pd.read_excel(BENCHMARK_FILE)

    # Sort if requested
    if SORT_BY_COUNT:
        df = df.sort_values("Equivalent Structures").reset_index(drop=True)

    # Select N cases
    df_selected = df.head(BENCHMARK_COUNT)

    benchmark_dir = os.path.dirname(BENCHMARK_FILE)

    # Prepare test cases
    test_cases = []
    for _, row in df_selected.iterrows():
        cif_file_relative = row["File"].replace(".cif", "_partial.cif")
        cif_file = os.path.join(benchmark_dir, cif_file_relative)

        # Parse supercell
        supercell_str = row["Supercell"]
        if pd.isna(supercell_str) or supercell_str == "1x1x1":
            supercell_array = None
        else:
            supercell_array = np.array([int(x) for x in supercell_str.split("x")])

        # Get expected count
        expected_count = int(row["Equivalent Structures"]) if not pd.isna(row["Equivalent Structures"]) else None

        test_cases.append(
            {
                "cif_file": cif_file,
                "supercell_array": supercell_array,
                "expected_count": expected_count,
                "file_basename": os.path.basename(cif_file),
            }
        )

    return test_cases


def get_coordinate_and_symbols(cif_filename, round_digits=4):
    """Load coordinates and symbols from CIF file.

    Args:
        cif_filename: Path to CIF file
        round_digits: Number of digits to round coordinates

    Returns:
        Dictionary with 'symbols' and 'coords' keys
    """
    import re

    cif_file = CifFile.from_file(cif_filename)
    # Assume single fragment
    key = list(cif_file.data.keys())[0]
    cif_dict = cif_file.data[key].data

    symbols = [str(re.sub("[^a-zA-Z]+", "", i)) for i in cif_dict["_atom_site_type_symbol"]]
    x = [str2float(i) for i in cif_dict["_atom_site_fract_x"]]
    y = [str2float(i) for i in cif_dict["_atom_site_fract_y"]]
    z = [str2float(i) for i in cif_dict["_atom_site_fract_z"]]
    coords = [[float(x), float(y), float(z)] for x, y, z in zip(x, y, z)]
    coords = np.round(coords, round_digits)
    condition = coords == 1.0000000
    coords[condition] = 0.0000000

    return {"symbols": symbols, "coords": coords}


def check_structure_consistency(str_cif_1, str_cif_2, symmops, round_digits=4):
    """Check if two structures are symmetrically equivalent.

    Args:
        str_cif_1: Path to first CIF file
        str_cif_2: Path to second CIF file
        symmops: Symmetry operations
        round_digits: Number of digits to round coordinates

    Returns:
        True if structures are equivalent, False otherwise
    """
    str_ref_1 = get_coordinate_and_symbols(str_cif_1, round_digits)
    str_ref_2 = get_coordinate_and_symbols(str_cif_2, round_digits)

    perms = []
    for _, symmop in enumerate(symmops):
        coords = str_ref_1["coords"]
        o_coords = symmop.operate_multi(coords)
        atol = 1e-8
        step = 10
        is_up = True

        for _ in range(10):
            matches = [find_in_coord_list_pbc(coords, coord, atol=atol) for coord in o_coords]

            match_lengths = [len(x) for x in matches]
            if all(x == 1 for x in match_lengths):
                break
            elif any(not x for x in match_lengths):
                # Too strict
                if not is_up:
                    step /= 2
                    is_up = True
                atol *= step
            elif any(x > 1 for x in match_lengths):
                # Too loose
                if is_up:
                    step /= 2
                    is_up = False
                atol /= step
        else:
            raise RuntimeError("Failed to build symmetry list.")

        indices = [x[0] for x in matches]
        perms.append(indices)

    # Check matching
    ref_coords = np.array(str_ref_2["coords"])
    ref_symbols = str_ref_2["symbols"]

    match_symbols = []
    for _, perm in enumerate(perms):
        coords = np.array(str_ref_1["coords"])
        symbols = str_ref_1["symbols"]
        perm_coords = coords  # Do not apply perm to coords (same as check_full_serial.py)
        perm_symbols = [symbols[i] for i in perm]

        # Lexsort ref1
        lex_perm_coords = np.array([perm_coords[i] for i in np.lexsort(perm_coords.T)])
        lex_perm_symbols = [perm_symbols[i] for i in np.lexsort(perm_coords.T)]

        # Lexsort ref2
        lex_ref_coords = np.array([ref_coords[i] for i in np.lexsort(ref_coords.T)])
        lex_ref_symbols = [ref_symbols[i] for i in np.lexsort(ref_coords.T)]

        coords_sum = np.sum(np.abs(np.array(lex_ref_coords) - np.array(lex_perm_coords)))
        coords_sum_each = np.abs(
            [np.array(ref_coord) - np.array(perm_coord) for ref_coord, perm_coord in zip(lex_ref_coords, lex_perm_coords)]
        )

        stol_coord = 1e0
        atol_coord = 1e-2

        if not coords_sum < stol_coord:
            continue
        if not any([np.sum(coord) < atol_coord for coord in coords_sum_each]):
            continue

        match_symbols.append(lex_perm_symbols == lex_ref_symbols)

    return any(match_symbols)


def _create_substitutor(cif_file, supercell_array=None):
    """Helper function to create Substitutor from CIF file.

    Args:
        cif_file: Path to CIF file
        supercell_array: Optional supercell expansion array

    Returns:
        Tuple of (substitutor, symmops_structure, error_dict or None)
    """
    # Load structure
    try:
        base_structure = LabeledStructure.from_file(cif_file)
    except Exception as e:
        return None, None, {"status": "error_loading", "error": str(e)}

    # Apply supercell if specified (keep a separate copy for symmops)
    structure_for_substitution = base_structure.copy()
    try:
        if supercell_array is not None:
            structure_for_substitution *= supercell_array
    except Exception as e:
        return None, structure_for_substitution, {"status": "error_supercell", "error": str(e)}

    structure_for_symmops = structure_for_substitution.copy()

    # Create Substitutor
    try:
        s = Substitutor(
            structure_for_substitution, symprec=SHRY_TOLERANCE, angle_tolerance=SHRY_ANGLE_TOLERANCE, atol=SHRY_ATOL
        )
        return s, structure_for_symmops, None
    except NeedSupercellError as e:
        return None, structure_for_symmops, {"status": "need_supercell", "error": str(e)}
    except Exception as e:
        return None, structure_for_symmops, {"status": "error_substitutor", "error": str(e)}


def _generate_all_structures(substitutor):
    """Helper function to generate all structures from Substitutor.

    Args:
        substitutor: Substitutor instance

    Returns:
        List of CIF text strings, or None if error
    """
    try:
        configs = []
        for packet in substitutor.quantities(["cifwriter"]):
            # Store CIF text to avoid shared mutable writer instances.
            configs.append(str(packet["cifwriter"]))
        return configs
    except Exception as e:
        return None


# ============================================================================
# Pytest test functions
# ============================================================================


@pytest.mark.parametrize("test_case", range(MAX_TEST_CASES), indirect=True)
def test_count_consistency(test_case):
    """Test that Substitutor.count() matches the actual number of generated structures."""
    cif_file = test_case["cif_file"]
    supercell_array = test_case["supercell_array"]

    # Skip if file doesn't exist
    if not os.path.exists(cif_file):
        pytest.skip(f"File not found: {cif_file}")

    # Create Substitutor
    s, structure, error = _create_substitutor(cif_file, supercell_array)
    if error:
        pytest.skip(f"{error['status']}: {error['error']}")

    # Get expected count
    count_expected = s.count()

    # Generate all structures
    configs = _generate_all_structures(s)
    if configs is None:
        pytest.fail("Failed to generate structures")

    count_actual = len(configs)

    # Check consistency
    assert count_expected == count_actual, f"Count mismatch: Substitutor.count()={count_expected}, but generated {count_actual}"


@pytest.mark.parametrize("test_case", range(MAX_TEST_CASES), indirect=True)
def test_benchmark_count(test_case):
    """Test that the generated structure count matches the benchmark."""
    cif_file = test_case["cif_file"]
    supercell_array = test_case["supercell_array"]
    expected_count = test_case["expected_count"]

    # Skip if file doesn't exist or no benchmark count
    if not os.path.exists(cif_file):
        pytest.skip(f"File not found: {cif_file}")
    if expected_count is None:
        pytest.skip("No benchmark count available")

    # Create Substitutor
    s, structure, error = _create_substitutor(cif_file, supercell_array)
    if error:
        pytest.skip(f"{error['status']}: {error['error']}")

    # Generate all structures
    configs = _generate_all_structures(s)
    if configs is None:
        pytest.fail("Failed to generate structures")

    count_actual = len(configs)

    # Check against benchmark
    assert count_actual == expected_count, f"Benchmark mismatch: expected {expected_count}, got {count_actual}"


# @pytest.mark.skip(reason="Manual run only; enable explicitly when needed.")
@pytest.mark.parametrize("test_case", range(MAX_REDUNDANCY_TEST_CASES), indirect=True)
def test_no_redundancy(test_case):
    """Test that all generated structures are symmetrically inequivalent (no redundancy)."""
    cif_file = test_case["cif_file"]
    supercell_array = test_case["supercell_array"]
    expected_count = test_case["expected_count"]

    # Skip if file doesn't exist
    if not os.path.exists(cif_file):
        pytest.skip(f"File not found: {cif_file}")

    # Skip if expected count is too large (redundancy check is very slow)
    if expected_count is not None and expected_count > 1000:
        pytest.skip(f"Skipping redundancy check for large count ({expected_count} structures)")

    # Create Substitutor
    s, structure, error = _create_substitutor(cif_file, supercell_array)
    if error:
        pytest.skip(f"{error['status']}: {error['error']}")

    # Generate all structures
    configs = _generate_all_structures(s)
    if configs is None:
        pytest.fail("Failed to generate structures")

    # Check redundancy (inline for readability)
    import tempfile

    output_dir = None  # "./redundancy_tmp"
    temp_context = tempfile.TemporaryDirectory() if output_dir is None else contextlib.nullcontext()

    with temp_context as tmpdir:
        if output_dir is None:
            output_dir = tmpdir
        else:
            os.makedirs(output_dir, exist_ok=True)

        cif_files = []
        for i, cif_text in enumerate(configs):
            tmp_cif = os.path.join(output_dir, f"config_{i:05d}.cif")
            with open(tmp_cif, "w", encoding="utf-8") as handle:
                handle.write(cif_text)
            cif_files.append(tmp_cif)

        # Get symmetry operations (supercell, pre-substitution)
        symmops = get_symmetry_operations_from_spglib(
            structure,
            symprec=SHRY_TOLERANCE,
            angle_tolerance=SHRY_ANGLE_TOLERANCE,
        )

        all_combinations = list(itertools.combinations(range(len(cif_files)), 2))

        def check_pair(cif_set):
            str_cif_1 = cif_files[cif_set[0]]
            str_cif_2 = cif_files[cif_set[1]]
            match_flag = check_structure_consistency(str_cif_1, str_cif_2, symmops)
            return (match_flag, cif_set if match_flag else None)

        results_generator = Parallel(n_jobs=N_JOBS, return_as="generator")(
            delayed(check_pair)(cif_set) for cif_set in all_combinations
        )

        pbar = tqdm.tqdm(
            total=len(all_combinations),
            desc="Checking pairs",
            **TQDM_CONF,
            disable=DISABLE_PROGRESSBAR,
        )
        for result in results_generator:
            pbar.update(1)
            match_flag, redundant_pair = result
            if match_flag:
                pbar.close()
                assert False, f"Found redundant pair: {redundant_pair} (structure indices are symmetrically equivalent)"
        pbar.close()


@pytest.fixture
def test_case(request, benchmark_data):
    """Fixture to get test case from benchmark data."""
    test_index = request.param
    if test_index >= len(benchmark_data):
        pytest.skip(f"Test index {test_index} out of range (only {len(benchmark_data)} cases)")
    return benchmark_data[test_index]
