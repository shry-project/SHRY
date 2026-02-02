# Copyright (c) SHRY Development Team.
# Distributed under the terms of the MIT License.

"""
Example 8: Target specific Wyckoff labels using an existing CIF.

This uses SmFe12.cif (bundled in examples/) and programmatically relabels
multiple Wyckoff sites (e.g., Fe1 and Fe2) to a shared label "Target" so that
the substitution concentration is enforced across both positions together.
"""

from shry import LabeledStructure, Substitutor

cif_file = "SmFe12.cif"

# Load CIF with labels
structure = LabeledStructure.from_file(cif_file)

# Choose which labels to merge into the shared pool
target_labels = {"Fe1", "Fe2"}

# Relabel those sites to a shared label "Target"
for site in structure:
    labels = set(site.properties.get("_atom_site_label", ()))
    if labels & target_labels:
        site.properties["_atom_site_label"] = ("Target",)

# Replace the shared label with a disordered composition
structure.replace_species({"Target": "Fe0.5Ti0.5"})

substitutor = Substitutor(structure)

# Enumerate and report counts
count = substitutor.count()
print(f"Unique structures: {count}")

# Example: write all structures with weights
# for i, packet in enumerate(substitutor.quantities(("cifwriter", "weight"))):
#     cifwriter = packet["cifwriter"]
#     weight = packet["weight"]
#     filename = f"cif_i{i}w{weight}.cif"
#     cifwriter.write_file(filename=filename)
