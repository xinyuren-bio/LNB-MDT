Preparation Module
===================

This module provides tools for generating lipid nanobubble (LNB) systems and preparing them for molecular dynamics simulations using the Martini coarse-grained force field.

Overview
--------

The preparation module includes:

- **LNB System Generator**: Command-line script to generate complete folder structures with all necessary files
- **Lipid Information**: Database of supported lipid types
- **Simulation Templates**: GROMACS parameter files for minimization, equilibration, and production runs

LNB System Generator
--------------------

`lnb_generate_with_folder.py`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Recommended script for command-line usage.** This script generates a complete folder structure automatically with all necessary files for GROMACS simulations.

**Key Features:**

- Automatically creates folder structure with all necessary files
- Generates ``system.gro``, ``system.top``, ``system.ndx`` with standard names
- Copies all ``.itp`` files to ``toppar/`` subfolder
- Copies simulation parameter files (``.mdp``)
- Creates ``parameter.txt`` with all generation parameters
- Automatically generates index file (``.ndx``) if not created by the underlying generator

**Usage:**

The ``-o`` parameter specifies the **folder name** (not file path). The script will create a folder named ``LNB_MDT_<folder_name>``.

**Example: Generate LNB with complete folder structure**

.. code-block:: bash

   python preparation/lnb_generate_with_folder.py \
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
       -o my_lnb_system

This will create a folder ``LNB_MDT_my_lnb_system/`` containing:

- ``system.gro`` - Structure file
- ``system.top`` - Topology file (already configured with proper includes)
- ``system.ndx`` - Index file with standard groups (System, lnb_layer, lnb_gas, solute)
- ``parameter.txt`` - Generation parameters and command history
- ``toppar/`` - Directory with all ``.itp`` files
- ``*.mdp`` files - Simulation parameter files (minimization, equilibration, production)

**Get Help:**

.. code-block:: bash

   python preparation/lnb_generate_with_folder.py -h

For detailed parameter descriptions, see:

.. code-block:: bash

   python preparation/lnb_gener_martini3.py -h

**Note:** All parameters are the same as ``lnb_gener_martini3.py``, except that ``-o`` specifies a folder name instead of a file path.

Complete Workflow Example
--------------------------

Step 1: Generate LNB System with Complete Folder Structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   python preparation/lnb_generate_with_folder.py \
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
       -o my_lnb_system

This will automatically create a folder ``LNB_MDT_my_lnb_system/`` with:

- ``system.gro`` - Structure file
- ``system.top`` - Topology file (already configured)
- ``system.ndx`` - Index file
- ``parameter.txt`` - Generation parameters
- ``toppar/`` - All ``.itp`` files
- ``*.mdp`` files - Simulation parameters

Step 2: Run Simulations
~~~~~~~~~~~~~~~~~~~~~~~~

The folder is ready for simulation. Navigate to the folder and run:

.. code-block:: bash

   cd LNB_MDT_my_lnb_system
   # Use the provided pipeline script or run GROMACS commands manually

Key Parameters
--------------

System Dimensions
~~~~~~~~~~~~~~~~~

- ``-d``: Distance between periodic images (nm)
- ``-r``: Radius of the membrane sphere (nm)
- ``-x``, ``-y``, ``-z``: Box dimensions (nm)
- ``-pbc``: Periodic boundary condition type (hexagonal, rectangular, square, cubic, optimal, or keep)

Lipid Composition
~~~~~~~~~~~~~~~~~

- ``-u``: Lipid type and relative abundance (can be specified multiple times)
  
  - Format: ``LIPIDNAME:NUMBER`` (e.g., ``DPPC:50``, ``DAPC:30``, ``CHOL:20``)
  - Supported lipids: See ``lipidsInfo_martini3.py`` for complete list
- ``-a``: Area per lipid (nm², default: 0.60)
- ``-rand``: Random kick size for lipid placement (default: 0.1)

Solvent and Gas
~~~~~~~~~~~~~~~

- ``-sol``: Solvent type (e.g., ``W`` for water)
- ``-gas``: Gas type inside the nanobubble (e.g., ``O2``, ``N2``, ``CO2``, ``H2``)
- ``-gden``: Gas density
- ``-gasd``: Gas diameter
- ``-sold``: Solvent diameter (default: 0.5)
- ``-salt``: Salt concentration (e.g., ``0.15`` for 150 mM)

Output
~~~~~~

- ``-o``: Output folder name (not file path) - **This is different from the underlying generator**

Advanced Options
----------------

Custom Lipid Compositions
~~~~~~~~~~~~~~~~~~~~~~~~~~

You can specify multiple lipid types with different abundances:

.. code-block:: bash

   python preparation/lnb_generate_with_folder.py \
       -u DPPC:50 \
       -u DAPC:30 \
       -u CHOL:20 \
       -u PA:10 \
       -o my_system

Multiple Gas Types
~~~~~~~~~~~~~~~~~~

Support for multiple gas molecules:

.. code-block:: bash

   python preparation/lnb_generate_with_folder.py \
       -gas O2 \
       -gden 200 \
       -gas N2 \
       -gden 100 \
       -o my_system

Custom Box Dimensions
~~~~~~~~~~~~~~~~~~~~~

Adjust simulation box size based on your system:

.. code-block:: bash

   python preparation/lnb_generate_with_folder.py \
       -x 30 -y 30 -z 30 \
       -pbc cubic \
       -o my_system

Troubleshooting
---------------

Common Issues
~~~~~~~~~~~~~

**Error: "Lipid type not found"**

- **Solution**: Check that the lipid name matches exactly with entries in ``lipidsInfo_martini3.py``
- Use ``python preparation/lnb_generate_with_folder.py -h`` to see supported options

**Error: "Box dimensions too small"**

- **Solution**: Ensure box dimensions are at least 2-3 times the membrane radius
- Increase ``-x``, ``-y``, ``-z`` values

**Error: "Gas density too high"**

- **Solution**: Reduce ``-gden`` value or increase box size
- Typical gas densities: 50-300 molecules per nm³

**Missing NDX file:**

- **Solution**: The script automatically generates the NDX file if it's not created by the underlying generator. If you encounter issues, check that MDAnalysis is properly installed.

Tips
~~~~

1. **Start with default parameters**: Use default values first, then adjust based on your specific needs
2. **Check lipid compatibility**: Ensure all lipid types are compatible with Martini 3.0 force field
3. **Validate topology**: Always check the generated ``system.top`` file for correct molecule counts
4. **Test with small systems**: Generate a small test system first before creating large production systems
5. **Monitor disk space**: Large systems with many water molecules can generate very large GRO files

References
----------

The LNB generator is based on the ``insane.py`` script by:

- Djurre H. de Jong et al. J. Chem. Theory Comput. 2013, 9, 687-697.

Modified and extended for lipid nanobubble generation by Yuan He and Yuxuan Wang in Dr. Xubo Lin's group at Beihang University.
