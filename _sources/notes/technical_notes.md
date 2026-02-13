(technical_notes)=

# Technical Notes

This section provides an in-depth look at the theoretical background and implementation details of SHRY.
For a formal mathematical description and extensive benchmarks, please refer to the original paper:
> **Shry: Application of Canonical Augmentation to the Atomic Substitution Problem**
> Genki I. Prayogo, Andrea Tirelli, Keishu Utimula, Kenta Hongo, Ryo Maezono, and Kousuke Nakano
> *Journal of Chemical Information and Modeling* **2022**, 62 (12), 2909–2915
> DOI: [10.1021/acs.jcim.2c00389](https://doi.org/10.1021/acs.jcim.2c00389)

## Overview

```{image} ../_static/Ce8Pd24Sb.pdf
:align: center
```
*Source: J. Chem. Inf. Model. 2022, 62, 2909−2915 (CC-BY-NC-ND 4.0)*

We employ a technique known as **Canonical Augmentation**. Essentially, this is a method for exhaustively exploring 'coloring' patterns. As illustrated in the figure (see below), the problem of atomic substitution in a crystal structure is mapped to a crystal site coloring problem. We search for unique coloring configurations by incrementally increasing the number of colored sites within a tree structure. Crucially, strictly local information is used to decide whether to prune or extend a branch from a given node, avoiding the need for global checks. This strategy prevents the explosion of computational complexity and memory usage. To enable these local decisions, a 'canonical' (standardized) representation of substitution information is adopted.

```{image} ../_static/tree_example.pdf
:align: center
```
*Source: J. Chem. Inf. Model. 2022, 62, 2909−2915 (CC-BY-NC-ND 4.0)*

The core problem SHRY solves is the **Atomic Substitution Problem**: given a crystal structure, how to generate all symmetry-inequivalent structures obtained by substituting a subset of atoms with different species.

In a naive approach, one might generate all possible permutations ($\binom{N}{k}$ combinations) and then group them by symmetry. However, this suffers from **combinatorial explosion**. For example, substituting atoms in a large supercell can easily result in trillions of candidates, making post-generation filtering impossible.

SHRY addresses this by using **Canonical Augmentation**, an algorithm that generates *only* the symmetry-inequivalent structures directly, bypassing the generation of redundant candidates.

## Method: Canonical Augmentation

Canonical Augmentation, originally introduced by Brendan McKay (1998), is a method for generating non-isomorphic combinatorial objects.

### The Algorithm
The problem is mapped to a **vertex coloring problem** on a graph where vertices represent atomic sites and colors represent atomic species.

1.  **Search Tree**: The generation process is visualized as traversing a search tree. The root is the unsubstituted structure. Each edge represents identifying one more atom to be substituted.
2.  **Canonical Reprensentation**: Every structure $X$ is assigned a "canonical parent" $p(X)$. The canonical parent is defined strictly based on the structure's properties (e.g., using a lexicographically minimal representation under symmetry operations).
3.  **Generation Rule**: A structure $X$ is accepted if and only if it was generated from its canonical parent $p(X)$. If $X$ is generated from a different parent $Y \neq p(X)$, it is discarded. This ensures that each unique structure is generated exactly once.
4.  **Isomorphic Rejection**: By checking local symmetry criteria, entire branches of the search tree leading to redundant structures can be pruned early.

Unlike **Orderly Generation** (Read, 1978), which requires every intermediate node to be canonical, Canonical Augmentation only requires the search tree to be traversed in a canonical manner. This makes it particularly flexible and efficient for this class of problems.

### Invariant and Optimization
A key optimization in SHRY is the use of **invariants**. To check if two structures are equivalent, or to find a canonical parent, SHRY computes invariants (such as the sum of distance matrices or other symmetry-invariant properties). If invariants differ, full isomorphism testing (which is expensive) can be avoided.

## Count Estimation: Polya's Enumeration Theorem
Before performing the actual enumeration, SHRY can estimate the total number of unique structures using **Polya's Enumeration Theorem**. This theorem relates the number of distinct colorings of a set under a group action to the cycle index of the group.

$$
|Z| = \frac{1}{|G|} \sum_{g \in G} \prod_{k=1}^n c_k^{j_k(g)}
$$

where $|Z|$ is the number of unique structures, $|G|$ is the order of the symmetry group, and $c_k$ terms relate to the cycle structure of the permutation $g$.

This estimate allows users to decide whether a full enumeration is computationally feasible before starting the run.

## Implementation Details

SHRY is implemented in pure **Python 3**.

### Core Libraries
*   **Pymatgen**: Used for general crystal structure handling (`Structure` objects), input/output (CIF files), and defining compositions.
*   **Spglib**: Used for detecting the symmetry group of the initial input structure.
*   **Numpy**: Used for high-performance matrix operations (e.g., coordinate transformations, distance matrices).

### Workflow
1.  **Initialization**:
    *   Read input CIF.
    *   Detect Space Group operation $G$ using `spglib`.
    *   Expand to supercell if requested (diagonal or non-diagonal expansions supported).
2.  **Pattern Generation (`shry.core.PatternMaker`)**:
    *   The `PatternMaker` class implements the Canonical Augmentation algorithm.
    *   It recursively successfully builds substitution patterns, pruning branches using symmetry operations derived from $G$.
3.  **Structure Construction**:
    *   For each unique pattern found, a new `Structure` object is created.
4.  **Post-Processing**:
    *   **Symmetrization**: (Optional) Use `spglib` to standardized the output structure.
    *   **Ewald Summation**: (Optional) Rank structures by electrostatic energy using Ewald summation (`pymatgen.analysis.ewald`). This is useful for finding typically stable configurations (low electrostatic energy) among thousands of candidates.

## Validation Test

To verify the reliability of SHRY, rigorous validation tests were conducted as described in *J. Chem. Inf. Model. 2022, 62, 2909−2915*. The validation process involves two main strategies: **Intra-software check** and **Inter-software check**.

### Intra-software and Inter-software Checks

```{image} ../_static/intra_and_inter_software_checks.pdf
:align: center
```
*Source: J. Chem. Inf. Model. 2022, 62, 2909−2915 (CC-BY-NC-ND 4.0)*

1.  **Intra-software check**: This validates the internal consistency of the algorithm. We confirm that the set of unique structures generated does not change depending on the choice of the initial structure. For example, applying a symmetry operation to the input structure or using a supercell with a different origin should result in an identical set of unique substituted structures (modulo symmetry). This ensures that the canonicalization process is robust and independent of the arbitrary choice of the starting representation.

2.  **Inter-software check**: This validates the correctness of the results by comparing them with other established software packages. We compared SHRY's output with:
    *   **Supercell** (Okhotnikov et al., 2016)

    The comparison involves verifying that the total count of unique structures matches exactly across different tools for various complex test cases (e.g., disordered perovskites, alloys). In all benchmarked cases, SHRY produced the exact same number of unique isomers as the reference software.

These tests confirm that SHRY correctly enumerates symmetry-inequivalent structures without missing any candidates or generating duplicates.

### Symmetry Detection

SHRY's symmetry detection relies on `spglib`, a crystal symmetry identification tool widely recognized for its numerical robustness. `spglib` implements both the ordinary space groups (230) and magnetic space groups (1651) [as known as Shubnikov groups], and SHRY supports both.

