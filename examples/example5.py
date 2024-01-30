# Copyright (c) SHRY Development Team.
# Distributed under the terms of the MIT License.

"""
Equivalent enumlib operations in SHRY
"""

from pymatgen.core import Structure
import shry
from shry import Substitutor

shry.const.DISABLE_PROGRESSBAR = True

# PbSnTe structure
cif_file = "PbSnTe.cif"
structure = Structure.from_file(cif_file)
structure *= (2, 2, 2)

# Generate ASE atoms instances with shry
s = Substitutor(structure)
# Shry uses generator; below is to put the Structures into a list
shry_ase_atoms = [x for x in s.ase_atoms_writers()]
shry_num_structs = s.count()
print(
    f"SHRY (group equivalent sites) resulted in {shry_num_structs} structures"
)
print(shry_ase_atoms[0])
