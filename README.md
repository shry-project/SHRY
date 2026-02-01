# README

![logo](https://github.com/shry-project/SHRY/blob/main/logo/logo.jpg?raw=true)

[![license](https://img.shields.io/github/license/shry-project/SHRY)](https://github.com/shry-project/SHRY/blob/main/LICENSE)
[![DL](https://img.shields.io/pypi/dm/SHRY)](https://pypi.org/project/SHRY/)
[![release](https://img.shields.io/github/release/shry-project/SHRY/all.svg)](https://github.com/shry-project/SHRY/releases)
[![PYPI_version](https://badge.fury.io/py/SHRY.svg)](https://pypi.org/project/SHRY/)
[![Python_version](https://img.shields.io/pypi/pyversions/SHRY)](https://pypi.org/project/SHRY/)
[![workflows](https://github.com/shry-project/SHRY/actions/workflows/shry-pytest.yml/badge.svg)](https://github.com/shry-project/SHRY/actions/workflows/shry-pytest.yml)
[![fork](https://img.shields.io/github/forks/shry-project/SHRY?style=social)](https://github.com/shry-project/SHRY/forks)
[![stars](https://img.shields.io/github/stars/shry-project/SHRY?style=social)](https://github.com/shry-project/SHRY/stargazers)

SHRY (**S**uite for **H**igh-th**r**oughput generation of models with atomic substitutions implemented by p**y**thon) is a tool for generating unique ordered structures corresponding to a given disordered structure.

## How to cite

Please cite the following paper:

- [SHRY: Application of Canonical Augmentation to the Atomic Substitution Problem](https://doi.org/10.1021/acs.jcim.2c00389),
  G.I. Prayogo*, A. Tirelli, K. Utimula, K. Hongo, R. Maezono, and K. Nakano*,
  *J. Chem. Inf. Model.*, 62, 2909-2915 (2022),
  DOI: [10.1021/acs.jcim.2c00389](https://doi.org/10.1021/acs.jcim.2c00389)

## Installation

SHRY can be obtained from PyPI

```console
pip install shry
```


### Development

If you prefer to install from source, instead follow the procedures below.

```console
git clone https://github.com/shry-project/SHRY.git
cd shry
pip install -e .
```

## Quick use

### Preparation of an input file (a CIF file)

You can prepare a CIF file with partial occupations.

```text
# label element x y z occupation
Sm1 Sm 0.000 0.00 0.00 1.000
Fe1 Fe 0.250 0.25 0.25 0.400
Nb1 Nb 0.250 0.25 0.25 0.600
Fe2 Fe 0.278 0.50 0.00 1.000
```

`SHRY` will automatically stop if the total occupancy of a site is either less or more than 1.0. To simulate vacancies, create a pseudo atom with species `X`.

### Check total symmetry-inequivalent structures

You can readily check the number of total symmetry-inquivalent structures using the following command.

```console
shry --count-only STRUCTURE_CIF
```

This operation is based on Polya enumeration and takes much less time than a proper generation.

### Creating supercell

Sometimes a supercell is required to fit in finer concentrations. `SHRY` accepts either 3-digit  (diagonal) or 9-digit (non-diagonal) format to specify the supercell's scaling matrix. For example a 2x2x1 supercell can be specified by either

```console
shry -s 2 2 1 --count-only STRUCTURE_CIF
```

or

```console
shry -s 2 0 0 0 2 0 0 0 1 --count-only STRUCTURE_CIF
```

### Generating unique structures

Finally, you can generate symmetry-inequivalent structures using the following command:

```console
shry -s 2 2 1 STRUCTURE_CIF
```

The generated symmetry-inequivalent structures are saved in `sliceXX` directories.

### Additional information

For additional information, you can use the help command:

```console
shry -h
```

or you can refer to the documentation.

## Documentation

The documentation is available [here](https://shry.readthedocs.io/en/latest/).

## Contributing to the project

Please work on your **forked** repository, and send a pull request to the `main` branch of the original GitHub repository.

If you want to contribute to the project, report a bug, or ask for a new feature, please [raise an issue](https://github.com/shry-project/SHRY/issues).

## Branches

- `main`: main branch.
- `devel*`: development branches.
- `rc`: the latest stable version ready for deployment of the package.
- `rc-gh-pages`: the latest stable version ready for deployment of the documentation.

Every time a change is pushed to the `main` or `devel*` branch, the `GitHub` workflow launches the implemented unit and integration tests (`shry-pytest.yml`) for the `main` and `devel*` branches).

## How to deploy the package

Once the `main` branch is merged into the `rc` branch, the `GitHub` workflow launches the implemented unit and integration tests (`shry-pytest.yml`) and test a deployment using `test-PyPI`. Then, once a tag is attached to (the latest) commit in the `rc` branch, the `GitHub` workflow checks the tag format (PEP 440 with the starting v, e.g., v0.1.0b4, v0.1.1, v1.0) and deploy the package to `PyPI`.

## Formatting

Formatting rules are written in `pyproject.toml`.

## Pre-commit

Pre-commit (https://pre-commit.com/) is mainly used for applying the formatting rules automatically. Therefore, it is strongly encouraged to use it at or before git-commit. Pre-commit is set-up and used in the following way:

- Installed by `pip install pre-commit`, `conda install pre_commit` or see
  https://pre-commit.com/#install.
- pre-commit hook is installed by `pre-commit install`.
- pre-commit hook is run by `pre-commit run --all-files`.

Unless running pre-commit, pre-commit.ci may push the fix at PR by github action. In this case, the fix should be merged by the contributor's repository.

## VSCode setting
- Not strictly, but VSCode's `settings.json` may be written like below

  ```json
  "ruff.lint.args": [
      "--config=${workspaceFolder}/pyproject.toml",
  ],
  "[python]": {
      "editor.defaultFormatter": "charliermarsh.ruff",
      "editor.codeActionsOnSave": {
          "source.organizeImports": "explicit"
      }
  },
  ```


## How to run tests

Tests are written using `pytest`. To run tests, `pytest` has to be installed. The tests can be run by

```bash
% pytest -s -v
```


## Citation of `SHRY`

If you used `SHRY` in your research project, **please** cite the following articles. This indeed helps the `SHRY` project to continue:

- "Shry: Application of Canonical Augmentation to the Atomic Substitution Problem",

  [G.I. Prayogo*, A. Tirelli, K. Utimula, K. Hongo, R. Maezono, and K. Nakano*, J. Chem. Inf. Model. 62, 2909-2915 (2022)](https://doi.org/10.1021/acs.jcim.2c00389)

  ```
  @article{doi:10.1021/acs.jcim.2c00389,
    author = {Prayogo, Genki Imam and Tirelli, Andrea and Utimula, Keishu and Hongo, Kenta and Maezono, Ryo and Nakano, Kousuke},
    title = {SHRY: Application of Canonical Augmentation to the Atomic Substitution Problem},
    journal = {J. Chem. Inf. Model.},
    volume = {62},
    number = {12},
    pages = {2909-2915},
    year = {2022},
    doi = {10.1021/acs.jcim.2c00389},
  }
  ```
