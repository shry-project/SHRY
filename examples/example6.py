# Copyright (c) SHRY Development Team.
# Distributed under the terms of the MIT License.

"""
Rank substituted structures by Ewald energy and keep the lowest 100.

This script demonstrates how to iterate over SHRY-generated structures,
compute Ewald energies, and retain only the most stable ones using a
bounded min-heap. Oxidation states are guessed for the Ewald calculation.
"""

import os
from heapq import heappush, heapreplace

from pymatgen.analysis.ewald import EwaldSummation
from shry import _ScriptHelper


from shry import _ScriptHelper


def main():
    helper = _ScriptHelper(
        structure_file="PbSnTe.cif",  # Replace with your CIF if desired
        from_species=("Sn",),
        to_species=("Sn0.5Se0.5",),
        scaling_matrix=(2, 2, 1),  # enlarge to break more symmetry and get multiple patterns
        n_jobs=1,
    )

    # Heap of (negative_energy, index, structure) so the largest negative
    # (i.e., highest energy among kept items) sits at the root.
    lowest_100 = []

    for i, structure in enumerate(helper.substitutor.structure_writers()):
        try:
            # Ensure oxidation states are set explicitly for this system.
            structure.remove_oxidation_states()
            structure.add_oxidation_state_by_element(
                {
                    "Pb": +2,
                    "Sn": +2,
                    "Se": +4,
                    "Te": -2,
                }
            )
            energy = EwaldSummation(structure).total_energy
        except Exception as exc:  # noqa: BLE001 - narrow for clarity in real workflows
            print(f"Skipping structure {i}: {exc}")
            continue

        entry = (-energy, i, structure)
        if len(lowest_100) < 100:
            heappush(lowest_100, entry)
        elif energy < -lowest_100[0][0]:
            heapreplace(lowest_100, entry)

    output_dir = "output_ewald"
    os.makedirs(output_dir, exist_ok=True)

    for rank, (neg_energy, idx, structure) in enumerate(sorted(lowest_100, key=lambda t: -t[0])):
        energy = -neg_energy
        filename = os.path.join(output_dir, f"ewald_rank{rank:03d}_e{energy:.6f}_idx{idx}.cif")
        structure.to(fmt="cif", filename=filename)
        print(f"Saved {filename}")


if __name__ == "__main__":
    main()
