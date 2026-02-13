(changelog)=

# Change Log

## v1.2.0b1

*   **magCIF / MSG support (experimental)**: `.mcif` inputs now parse magCIF symmetry operations and `_atom_site_moment.*` into `Structure` site properties, and switch symmetry handling to spglib's magnetic symmetry (MSG). `.cif` inputs continue to use conventional space-group symmetry (SG).

*   **Parallel tree search**: Tree search now supports dynamic parallel execution using Python `multiprocessing`.
    *   Added `n_jobs` parameter to `Substitutor` and `PatternMaker` classes (default: `-1`, uses all CPU cores)
    *   Parallelization is applied at subtree level for both `_invarless_search` and `_invar_search`
    *   **Adaptive depth expansion**: subtree seeds are expanded until enough splits are prepared (`target = max(n_cores * 8, 32)`) or depth limit is reached (`max_depth <= 10`, and bounded by `stop`)
    *   **Dynamic load balancing**: worker tasks are scheduled with `multiprocessing.Pool(...).imap_unordered(..., chunksize=1)`
    *   **Progress reporting**: parallel mode reports subtree processing progress (`Processing subtrees`), while sequential mode reports pattern generation progress (`Making patterns`)
    *   Logging includes expansion and mode details, e.g. `Expanding to depth ...` and `Parallel search (dynamic): ...`


### Internal Changes

*   **CIF reader/writer backend update**: Replaced CIF read/write path from `pymatgen.io.cif` to `PyCifRW` (`CifFile`).

*   **Refactored symmetry detection**: Changed from pymatgen's internal `SpacegroupAnalyzer` to direct `spglib` calls for improved transparency and control over symmetry operations.
    *   Added helper functions: `structure_to_spglib_cell()`, `get_symmetry_operations_from_spglib()`, and `get_symmetrized_structure_from_spglib()`
    *   Removed dependency on `PatchedSpacegroupAnalyzer` class
    *   Symmetry information is now obtained directly from spglib's dataset


## v1.2.0a1

- Initial release of **SHRY** documentation.

### Key Features

*   **Canonical Augmentation**: Efficient generation of substitution patterns.
*   **Supercell support**: Easy supercell expansion with various matrix specifications.
*   **Pymatgen integration**: Uses Pymatgen structures for IO and analysis.

### Examples

*   **Example 6**: Ranking substitutions by Ewald energy.
*   **Example 7**: Enumerates highest symmetry structures. Output filenames now include space group number and symbol.
*   **Example 8**: New example demonstrating how to merge multiple Wyckoff labels into a shared label to enforce a global substitution concentration.

### Dependency Updates

*   Python version support updated to `>=3.10.0` (dropped support for 3.8 and 3.9).
*   Updated core dependencies:
    *   `numpy >= 2.0.0`
    *   `pymatgen >= 2025.4.10`
    *   `scipy >= 1.14.0`
    *   `spglib >= 2.7.0`

### Known Limitations

*   The public API is currently unstable and subject to change.
*   Major refactoring is currently in progress.
