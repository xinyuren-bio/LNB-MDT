Preparation Module
===================

This module provides tools for generating lipid nanobubble (LNB) systems and preparing them for molecular dynamics simulations using the Martini coarse-grained force field.

Overview
--------

The preparation module includes:

- **LNB System Generator**: Scripts to generate spherical lipid monolayer systems (lipid nanobubbles)
- **Index File Generator**: Utility to create GROMACS index files
- **Lipid Information**: Database of supported lipid types
- **Simulation Templates**: GROMACS parameter files for minimization, equilibration, and production runs

Components
----------

1. LNB Generator Scripts
~~~~~~~~~~~~~~~~~~~~~~~~

`lnb_gener_martini3.py` (Martini 3.0)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Main script for generating lipid nanobubble systems using Martini 3.0 force field.

**Key Parameters:**

**System Dimensions:**

- ``-d``: Distance between periodic images (nm)
- ``-r``: Radius of the membrane sphere (nm)
- ``-x``, ``-y``, ``-z``: Box dimensions (nm)
- ``-pbc``: Periodic boundary condition type (hexagonal, rectangular, square, cubic, optimal, or keep)

**Lipid Composition:**

- ``-u``: Lipid type and relative abundance (can be specified multiple times)
  
  - Format: ``LIPIDNAME:NUMBER`` (e.g., ``DPPC:50``, ``DAPC:30``, ``CHOL:20``)
  - Supported lipids: See ``lipidsInfo_martini3.py`` for complete list
- ``-a``: Area per lipid (nm², default: 0.60)
- ``-rand``: Random kick size for lipid placement (default: 0.1)

**Solvent and Gas:**

- ``-sol``: Solvent type (e.g., ``W`` for water)
- ``-gas``: Gas type inside the nanobubble (e.g., ``O2``, ``N2``, ``CO2``, ``H2``)
- ``-gden``: Gas density
- ``-gasd``: Gas diameter
- ``-sold``: Solvent diameter (default: 0.5)
- ``-salt``: Salt concentration (e.g., ``0.15`` for 150 mM)

**Output:**

- ``-o``: Output GRO file path
- ``-p``: Output topology file path (optional)

**Example: Generate LNB with DPPC:DAPC:CHOL = 5:3:2 and O₂**

.. code-block:: bash

   python preparation/lnb_gener_martini3.py \
       -d 1 \
       -r 5 \
       -x 20 -y 20 -z 20 \
       -u DPPC:50 \
       -u DAPC:30 \
       -u CHOL:20 \
       -sol W \
       -salt 0.15 \
       -a 1 \
       -gas O2 \
       -gden 200 \
       -o cases_lnb/lnb.gro \
       -p cases_lnb/lnb.top

**Get Help:**

.. code-block:: bash

   python preparation/lnb_gener_martini3.py -h

`lnb_gener_martini2.py` (Martini 2.x)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Similar script for Martini 2.x force field. Usage is identical to the Martini 3.0 version.

2. Index File Generator
~~~~~~~~~~~~~~~~~~~~~~~

`ndx_generator.py`
^^^^^^^^^^^^^^^^^

Automatically generates GROMACS index files (``.ndx``) from GRO files with standard groups.

**Function Signature:**

.. code-block:: python

   from preparation.ndx_generator import generate_ndx

   generate_ndx(
       gro_path="system.gro",
       ndx_path="system.ndx",
       lipid_resnames=["DPPC", "DAPC", "CHOL"],
       gas_resnames=["O2", "N2"],
       water_resnames=["W"],
       ion_resnames=["NA", "CL"]
   )

**Generated Groups:**

- ``System``: All atoms
- ``IONS``: Ion residues
- ``Lipids``: Lipid residues
- ``Gas``: Gas residues
- ``LNB``: Lipids + gas residues
- ``W``: Water residues

**Note:** The LNB generator automatically calls this function when generating GRO files, so manual generation is usually not necessary.

3. Simulation Files
~~~~~~~~~~~~~~~~~~~

The ``files/`` directory contains GROMACS parameter files:

**Force Field Files (Martini 3.0):**

- ``martini_v3.0.0.itp``: Main force field file
- ``martini_v3.0.0_phospholipids.itp``: Phospholipid definitions
- ``martini_v3.0.0_gas.itp``: Gas molecule definitions
- ``martini_v3.0.0_ions.itp``: Ion definitions
- ``martini_v3.0.0_solvents.itp``: Solvent definitions
- ``martini_v3.0.0_ffbonded.itp``: Bonded interactions

**Simulation Parameter Files:**

- ``minimization.mdp``: Energy minimization parameters
- ``equilibration_1.mdp`` to ``equilibration_4.mdp``: Equilibration steps
- ``production.mdp``: Production run parameters
- ``pipeline.sh``: Automated simulation pipeline script

Complete Workflow Example
--------------------------

Step 1: Generate LNB System
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   python preparation/lnb_gener_martini3.py \
       -d 1 \
       -r 5 \
       -x 20 -y 20 -z 20 \
       -u DPPC:50 \
       -u DAPC:30 \
       -u CHOL:20 \
       -sol W \
       -salt 0.15 \
       -a 1 \
       -gas O2 \
       -gden 200 \
       -o system.gro \
       -p system.top

This will generate:

- ``system.gro``: Structure file
- ``system.top``: Topology file
- ``system.ndx``: Index file (auto-generated)

Step 2: Prepare Topology File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Edit ``system.top`` to include the necessary force field includes:

.. code-block:: text

   #include "preparation/files/martini_v3.0.0.itp"
   #include "preparation/files/martini_v3.0.0_phospholipids.itp"
   #include "preparation/files/martini_v3.0.0_gas.itp"
   #include "preparation/files/martini_v3.0.0_ions.itp"
   #include "preparation/files/martini_v3.0.0_solvents.itp"

   [ system ]
   LNB System

   [ molecules ]
   ; name  number
   DPPC    50
   DAPC    30
   CHOL    20
   O2      200
   W       5000
   NA      150
   CL      150

Step 3: Run Simulations
~~~~~~~~~~~~~~~~~~~~~~~~

Use the provided pipeline script or run manually:

**Using Pipeline Script:**

.. code-block:: bash

   cd preparation/files
   chmod +x pipeline.sh
   ./pipeline.sh

Advanced Options
---------------

Custom Lipid Compositions
~~~~~~~~~~~~~~~~~~~~~~~~~~

You can specify multiple lipid types with different abundances:

.. code-block:: bash

   python preparation/lnb_gener_martini3.py \
       -u DPPC:50 \
       -u DAPC:30 \
       -u CHOL:20 \
       -u PA:10 \
       -o system.gro

Multiple Gas Types
~~~~~~~~~~~~~~~~~~

Support for multiple gas molecules:

.. code-block:: bash

   python preparation/lnb_gener_martini3.py \
       -gas O2 \
       -gden 200 \
       -gas N2 \
       -gden 100 \
       -o system.gro

Custom Box Dimensions
~~~~~~~~~~~~~~~~~~~~~

Adjust simulation box size based on your system:

.. code-block:: bash

   python preparation/lnb_gener_martini3.py \
       -x 30 -y 30 -z 30 \
       -pbc cubic \
       -o system.gro

Troubleshooting
---------------

Common Issues
~~~~~~~~~~~~~

**Error: "Lipid type not found"**

- **Solution**: Check that the lipid name matches exactly with entries in ``lipidsInfo_martini3.py``
- Use ``python preparation/lnb_gener_martini3.py -h`` to see supported options

**Error: "Box dimensions too small"**

- **Solution**: Ensure box dimensions are at least 2-3 times the membrane radius
- Increase ``-x``, ``-y``, ``-z`` values

**Error: "Gas density too high"**

- **Solution**: Reduce ``-gden`` value or increase box size
- Typical gas densities: 50-300 molecules per nm³

**Index file missing atoms**

- **Solution**: This was fixed in recent versions. Ensure you're using the latest ``ndx_generator.py``
- The generator now correctly handles gro files with empty lines after box vectors

Tips
----

1. **Start with default parameters**: Use default values first, then adjust based on your specific needs
2. **Check lipid compatibility**: Ensure all lipid types are compatible with Martini 3.0 force field
3. **Validate topology**: Always check the generated ``.top`` file for correct molecule counts
4. **Test with small systems**: Generate a small test system first before creating large production systems
5. **Monitor disk space**: Large systems with many water molecules can generate very large GRO files

References
----------

The LNB generator is based on the ``insane.py`` script by:

- Djurre H. de Jong et al. J. Chem. Theory Comput. 2013, 9, 687-697.

Modified and extended for lipid nanobubble generation by Yuan He and Yuxuan Wang in Dr. Xubo Lin's group at Beihang University.

