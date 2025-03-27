README
==========

.. warning::
    Dear Users. If you are using `SHRY`, would you mind giving it a **star** from the button at the top right on the GitHub repo.? Maintaining and updating this software is not the primary job of the project owner. If there are users out there, we would like to continue the maintenance and updates, so we would like to gather some statistics. We cannot see how many people are using it from the developers' side. Your corporations would be greatly appreciated.

|logo|

|license| |DL| |release| |PYPI_version| |Python_version| |workflows| |fork| |stars|

.. |logo| image:: https://github.com/shry-project/SHRY/blob/main/logo/logo.jpg?raw=true
.. |license| image:: https://img.shields.io/github/license/shry-project/SHRY
.. |release| image:: https://img.shields.io/github/release/shry-project/SHRY/all.svg
.. |DL| image:: https://img.shields.io/pypi/dm/SHRY
.. |Python_version| image:: https://img.shields.io/pypi/pyversions/SHRY
.. |fork| image:: https://img.shields.io/github/forks/shry-project/SHRY?style=social
.. |stars| image:: https://img.shields.io/github/stars/shry-project/SHRY?style=social
.. |workflows| image:: https://github.com/shry-project/SHRY/actions/workflows/shry-pytest.yml/badge.svg
.. |PyPI_version| image:: https://badge.fury.io/py/SHRY.svg

SHRY (\ **S**\ uite for \ **H**\ igh-th\ **r**\ oughput generation of models
with atomic substitutions implemented by p\ **y**\ thon)
is a tool for generating unique ordered structures
corresponding to a given disordered structure.

How to cite
-------------
Please cite the following paper:

`SHRY: Application of Canonical Augmentation to the Atomic Substitution Problem <https://doi.org/10.1021/acs.jcim.2c00389>`_, G.I. Prayogo*, A. Tirelli, K. Utimula, K. Hongo, R. Maezono, and K. Nakano*, *J. Chem. Inf. Model.*, 62, 2909-2915 (2022), `DOI:10.1021/acs.jcim.2c00389 <https://doi.org/10.1021/acs.jcim.2c00389>`_

Installation
------------

SHRY can be obtained from PyPI

.. code-block:: console

    pip install shry

For Windows Users
^^^^^^^^^^^^^^^^^

For Windows users,
if you don't have Python already,
you can try, for example,
installing Python from the Microsoft store
following instructions on
`this page`_.

.. _`this page`: https://docs.microsoft.com/en-us/windows/python/beginners

Then install SHRY just like above
within PowerShell or your favourite terminal.

.. code-block:: console

    pip install shry

Development
^^^^^^^^^^^

If you prefer to install from source,
instead follow the procedures below.

.. code-block:: console

    git clone https://github.com/shry-project/SHRY.git
    cd shry
    pip install -e .

Quick use
---------

Preparation of an input file (a CIF file)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can prepare a CIF file with partial occupations.

.. code-block::

    # label element x y z occupation
    Sm1 Sm 0.000 0.00 0.00 1.000
    Fe1 Fe 0.250 0.25 0.25 0.400
    Nb1 Nb 0.250 0.25 0.25 0.600
    Fe2 Fe 0.278 0.50 0.00 1.000

``SHRY`` will automatically stop if the total occupancy of a site is
either less or more than 1.0. To simulate vacancies, create a pseudo
atom with species ``X``.

Check total symmetry-inequivalent structures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can readily check the number of total symmetry-inquivalent
structures using the following command.

.. code-block:: console

    shry --count-only STRUCTURE_CIF

This operation is based on Polya enumeration and takes much less time
than a proper generation.

Creating supercell
^^^^^^^^^^^^^^^^^^

Sometimes a supercell is required to fit in finer concentrations.
``SHRY`` accepts either 3-digit (diagonal) or 9-digit (non-diagonal)
format to specify the supercell's scaling matrix. For example a 2x2x1
supercell can be specified by either

.. code-block:: console

    shry -s 2 2 1 --count-only STRUCTURE_CIF # depending on your environment.
    shry -s='2 2 1' --count-only STRUCTURE_CIF # depending on your environment.

or

.. code-block:: console

    shry -s 2 0 0 0 2 0 0 0 1 --count-only STRUCTURE_CIF # depending on your environment.
    shry -s='2 0 0 0 2 0 0 0 1' --count-only STRUCTURE_CIF # depending on your environment.

Generating unique structures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Finally, you can generate symmetry-inequivalent structures using the
following command:

.. code-block:: console

    shry -s 2 2 1 STRUCTURE_CIF # depending on your environment.
    shry -s='2 2 1' STRUCTURE_CIF # depending on your environment.

The generated symmetry-inequivalent structures are saved in sliceXX
directories.

Additional information
^^^^^^^^^^^^^^^^^^^^^^

For additional information, you can use the help command:

.. code-block:: console

    shry -h

or you can refer to the documentation.

Documentation
-------------

The documentation is available `here <https://shry.readthedocs.io/en/latest/>`_.


Contributing to the project
---------------------------
Please work on your **forked** repository, and send a pull request to the devel branch of the original GitHub repository. 
After the pull request is approved and the devel branch is merged.

If you want to contribute to the project, report a bug, or ask for
a new feature, please `raise an issue <https://github.com/shry-project/SHRY/issues>`_.


How to release for maintainers
----------------------------------
Check the next-version version

.. code-block:: console

    # Confirm the version number via `setuptools-scm`
    python -m setuptools_scm
    e.g., 1.1.4.dev28+gceef293.d20221123 -> <next-version> = v1.1.4 or v1.1.4-alpha(for pre-release)

Add and push with the new tag

.. code-block:: console

    # Push with tag
    git tag <next-version>  # e.g., git tag v1.1.4  # Do not forget "v" before the version number!
    git push origin devel --tags  # or to the new branch # Do not push to the mater branch!!

Make a pull request for merging the devel branch to the main branch by hand. The implemented GitHub Action checks if the automatic deploy works using test-pyPI.

Finally, do a new release with a release note on GitHub. The new release triggers an implemented GitHub Action that automatically uploads the package to PyPI (if the commit is tagged correctly, e.g., v1.1.0 or v1.1.0-alpha).

