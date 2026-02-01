# Examples

What each script does, summarized per example. All required CIF inputs are bundled in this directory.

## Example 1: Basic enumeration and CIF output
- File: [examples/example1.py](examples/example1.py)
- What it does: Load `SmFe7Ti.cif`, use `Substitutor` to generate all symmetry-inequivalent structures, and write CIFs with their multiplicity (weight).
- Key point: Uses `Substitutor.quantities(("cifwriter", "weight"))` and saves to `output/` as `cif_i{index}w{weight}.cif`.

## Example 2: Comparison with enumlib
- File: [examples/example2.py](examples/example2.py)
- What it does: Expand `PbSnTe.cif` to a 2×2×2 supercell and compare enumeration results between Pymatgen's `EnumlibAdaptor` and SHRY.
- Comparison modes:
  - SHRY default (group equivalent sites).
  - SHRY with `groupby=lambda x: x.species` (group by species, similar to enumlib behavior).
- Output: Prints the number of structures from each approach to stdout.

## Example 3: LabeledStructure and label-based substitutions
- File: [examples/example3.py](examples/example3.py)
- What it does: Load `SmFe12.cif` as a `LabeledStructure`, replace labels `Fe1`→`Fe3Ti`, `Fe2`→`Fe7Ti`, `Fe3`→`Ti`, then expand to a non-diagonal supercell `((2,0,1),(0,1,0),(1,0,1))`. Generates all structures, writes CIFs and weights to `output/`, and prints configuration letters (e.g., `aba...`).
- Key point: Keeps `_atom_site_label`, so substitutions remain label-consistent even after supercell expansion.

## Example 4a: Caching pattern generation and saving it
- File: [examples/example4a.py](examples/example4a.py)
- What it does: Load `SmFe12.cif` as a `LabeledStructure`, replace `Fe1`→`Fe7Ti`, run `Substitutor(cache=True)`, write CIFs, and pickle the internal `pattern_makers` to `pg.pkl`.
- Benefit: Reuse the cached pattern makers to skip expensive recomputation for symmetry-equivalent systems later.

## Example 4b: Reusing cache and applying to different compositions
- File: [examples/example4b.py](examples/example4b.py)
- What it does: Load `pg.pkl` from Example 4a. Starting from `SmFe12.cif`, evaluate two substitution cases (`Fe`→`Fe3Ti` and `Fe`→`FeTi3`) in sequence. Attach the loaded `pattern_makers` to a `Substitutor`, generate CIFs for `run1` and `run2`, then overwrite `pg.pkl` with the updated cache.
- Key point: When symmetry is the same, you can rapidly iterate different substitution patterns using the cached pattern makers.

## Example 5: Generating ASE Atoms
- File: [examples/example5.py](examples/example5.py)
- What it does: Expand `PbSnTe.cif` to 2×2×2, enumerate all structures with `Substitutor`, and obtain them as ASE `Atoms` objects.
- Output: Prints the number of structures to stdout and shows one sample `Atoms` instance; useful as a bridge to downstream ASE workflows.

## Example 6: Ranking by Ewald energy
- File: [examples/example6.py](examples/example6.py)
- What it does: Enumerates substitutions, guesses oxidation states, computes Ewald energies, and writes the 100 lowest-energy structures to output_ewald/.
- Note: If oxidation states cannot be guessed, Ewald energies may be zero or skipped.

## Example 7: Highest symmetry structures (lowest weight)
- File: [examples/example7.py](examples/example7.py)
- What it does: Enumerates substitutions and writes the top-N structures with the smallest weight (highest symmetry) to output_symm_top/. Optionally writes all structures sorted by weight.

## Example 8: Merge multiple Wyckoff labels for shared concentration
- File: [examples/example8.py](examples/example8.py)
- What it does: Loads `SmFe12.cif`, relabels multiple Wyckoff labels (e.g., `Fe1` and `Fe2`) into a shared label (`Target`), then applies a single disordered composition (e.g., `Fe0.5Ti0.5`) across that pooled site set. This enables specifying a global substitution concentration across multiple Wyckoff positions.
