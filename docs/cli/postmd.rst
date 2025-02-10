PostMD CLI
==========

The ``postmd`` command provides a unified interface for various molecular dynamics tools and snippets.

Synopsis
--------

.. code-block:: bash

    postmd <command> [options]

Commands
--------

The following commands are available:

lmp
^^^
LAMMPS-related snippets and tools.

Options:
    --pre OPTION     Print pre-processing snippets. Available options:
                     
                     * all: Show all snippets
                     * msd: Mean Square Displacement
                     * vacf: Velocity Autocorrelation Function
                     * fric_coeff: Friction Coefficient
                     * 1d_profile: 1D Density Profile
                     * radial_profile: Radial Density Profile
                     * viscosity: Viscosity Calculation
                     * rdf: Radial Distribution Function
                     * dump: Trajectory Dump
                     * us: Umbrella Sampling

    --post OPTION    Print post-processing snippets. Available options:
                     
                     * all: Show all snippets
                     * msd: Mean Square Displacement Analysis
                     * us: Umbrella Sampling Analysis

plumed
^^^^^^
PLUMED-related snippets and tools.

Options:
    --pre OPTION     Print pre-processing snippets. Available options:
                     
                     * all: Show all snippets
                     * us: Umbrella Sampling
                     * metad: Metadynamics
                     * well-tempered metad: Well-tempered Metadynamics

    --post OPTION    Print post-processing snippets. Available options:
                     
                     * all: Show all snippets
                     * us: Umbrella Sampling Analysis
                     * metad: Metadynamics Analysis
                     * driver: PLUMED Driver Usage

cp2k
^^^^
CP2K-related snippets and tools.

Options:
    --pre OPTION     Print pre-processing snippets
    --post OPTION    Print post-processing snippets

Examples
--------

1. Show LAMMPS MSD pre-processing snippet:

.. code-block:: bash

    postmd lmp --pre msd

2. Show PLUMED umbrella sampling post-processing snippet:

.. code-block:: bash

    postmd plumed --post us

3. Show all CP2K pre-processing snippets:

.. code-block:: bash

    postmd cp2k --pre all

Notes
-----

* The snippets are meant to be used as templates and may need to be modified according to your specific needs.
* Some snippets may require additional software or dependencies to be installed.
* For detailed information about each snippet, please refer to the respective software's documentation.

See Also
--------

* `LAMMPS Documentation <https://docs.lammps.org/>`_
* `PLUMED Documentation <https://www.plumed.org/doc>`_
* `CP2K Documentation <https://manual.cp2k.org/>`_ 