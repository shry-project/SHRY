(changelog)=

# Change Log

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

*   Multi-core processing is not yet supported.
*   The public API is currently unstable and subject to change.
*   Major refactoring is currently in progress.
