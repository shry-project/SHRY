# Copyright (c) SHRY Development Team.
# Distributed under the terms of the MIT License.

"""
List substituted structures with the highest symmetry (lowest weight).

In SHRY, the multiplicity (weight) of each inequivalent structure is
inversely proportional to the number of possible symmetry operations.
So the smallest weight corresponds to the highest symmetry.
"""

import os
from heapq import heappush, heapreplace

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from shry import _ScriptHelper


cif_file = "PbSnTe.cif"  # Replace with your CIF if desired

helper = _ScriptHelper(
    structure_file=cif_file,
    from_species=("Sn",),
    to_species=("Sn0.25Sb0.75",),
    scaling_matrix=(2, 2, 2),
)

TOP_N = 5
WRITE_ALL_SORTED = False

# Heap containing TOP_N structures with the lowest weight (highest symmetry).
lowest = []

for i, packet in enumerate(helper.substitutor.quantities(("weight", "structure"))):
    weight = packet["weight"]
    structure = packet["structure"]

    if len(lowest) < TOP_N:
        heappush(lowest, (-weight, i, structure))
    elif weight < -lowest[0][0]:
        heapreplace(lowest, (-weight, i, structure))

output_dir = "output_symm_top"
os.makedirs(output_dir, exist_ok=True)

for rank, (neg_weight, idx, structure) in enumerate(sorted(lowest, key=lambda t: -t[0])):
    weight = -neg_weight
    sga = SpacegroupAnalyzer(structure)
    sg_number = sga.get_space_group_number()
    sg_symbol = sga.get_space_group_symbol().replace("/", "_")
    filename = os.path.join(
        output_dir,
        f"symm_rank{rank:03d}_{sg_number}_{sg_symbol}_w{weight}_idx{idx}.cif",
    )
    structure.to(fmt="cif", filename=filename)
    print(f"Saved {filename}")

if WRITE_ALL_SORTED:
    output_dir_all = "output_symm_all"
    os.makedirs(output_dir_all, exist_ok=True)

    structures = [
        (packet["weight"], i, packet["structure"])
        for i, packet in enumerate(helper.substitutor.quantities(("weight", "structure")))
    ]
    structures.sort(key=lambda t: t[0])

    for rank, (weight, idx, structure) in enumerate(structures):
        sga = SpacegroupAnalyzer(structure)
        sg_number = sga.get_space_group_number()
        sg_symbol = sga.get_space_group_symbol().replace("/", "_")
        filename = os.path.join(
            output_dir_all,
            f"symm_rank{rank:03d}_{sg_number}_{sg_symbol}_w{weight}_idx{idx}.cif",
        )
        structure.to(fmt="cif", filename=filename)
        print(f"Saved {filename}")
