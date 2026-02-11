# Copyright (c) SHRY Development Team.
# Distributed under the terms of the MIT License.
"""Test core operations."""

# python modules
import pandas as pd
import pytest

# shry
from shry.core import (
    Substitutor,
)
from shry.main import LabeledStructure
from shry import const
from helper import chdir

# Tolerances
SHRY_TOLERANCE = const.DEFAULT_SYMPREC
SHRY_ANGLE_TOLERANCE = const.DEFAULT_ANGLE_TOLERANCE
SHRY_ATOL = const.DEFAULT_ATOL

pytestmark = [
    pytest.mark.filterwarnings(
        "ignore:No Pauling electronegativity for .*:UserWarning",
    ),
    pytest.mark.filterwarnings(
        r"ignore:Set OLD_ERROR_HANDLING to false and catch the errors directly\.:DeprecationWarning",
    ),
]


# @pytest.mark.skip(reason="Comprehensive but time consuming. It will be activated later.")
@chdir("./test_cifs/all_space_groups")
def test_benchmark():
    """Benchmark / the number of symmetry-inequivalent structures."""
    df = pd.read_excel("./benchmark_SG_all.xls")
    for _, zipped in enumerate(
        zip(
            df["Supercell"],
            df["File"],
            df["Substitutions"],
            df["Equivalent Structures"],
            df["Checked"],
            strict=True,
        )
    ):
        supercell, filename, substitution, equivalent, checked = zipped

        # cif_basename = os.path.basename(filename).replace(".cif", "")
        filename = filename.replace(".cif", "_partial.cif")
        print(f"filename={filename}")
        structure = LabeledStructure.from_file(filename)
        supercell_size = list(map(int, supercell.split("x")))
        print(supercell_size)
        structure *= supercell_size
        s = Substitutor(structure, symprec=SHRY_TOLERANCE, angle_tolerance=SHRY_ANGLE_TOLERANCE, atol=SHRY_ATOL)
        count_obtained = s.count()
        count_ref = equivalent
        print(f"count_obtained={count_obtained}")
        print(f"count_ref={count_ref}")
        assert count_obtained == count_ref
