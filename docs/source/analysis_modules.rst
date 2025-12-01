Analysis Modules Details
========================

This directory contains various analysis modules for lipid nanobubble molecular dynamics simulations. All modules can be executed from the command line and support both single-frame (GRO only) and multi-frame (GRO + XTC) trajectory analysis.

Prerequisites
-------------

Before running the analysis modules, ensure you have:

- Activated the LNB-MDT conda environment: ``conda activate LNB-MDT``
- Example data files in ``cases_lnb/`` directory (see ``cases_lnb/README.md`` for details)

General Usage Pattern
---------------------

All analysis modules follow a similar command-line interface pattern. **Run from the project root directory:**

.. code-block:: bash

   python analysis/<module_name>.py [OPTIONS]

Common options available across most modules:

- ``--gro-file, -g``: Path to GRO topology file (default: ``cases_lnb/lnb.gro``)
- ``--xtc-file, -x``: Path to XTC trajectory file (optional for some modules)
- ``--output-csv, -o``: Output CSV file path
- ``--residues, -r``: Residue groups dictionary string
- ``--parallel, -p``: Enable parallel processing
- ``--n-jobs, -j``: Number of parallel jobs (-1 for all cores)
- ``--start-frame, -s``: Starting frame (0-indexed)
- ``--stop-frame, -e``: Stopping frame (exclusive)
- ``--step-frame, -t``: Frame step size
- ``--verbose, -v``: Enable verbose output

Parameter Reference Table
--------------------------

.. list-table:: Parameter Reference Table
   :header-rows: 1
   :widths: 8 15 30 20 27

   * - Short Parameter
     - Long Parameter
     - Description
     - Default Value
     - Example
   * - ``-g``
     - ``--gro-file``
     - GRO/TPR file path
     - ``cases_lnb/lnb.gro``
     - ``-g cases_lnb/lnb.gro``
   * - ``-x``
     - ``--xtc-file``
     - XTC/TRR file path
     - -
     - ``-x cases_lnb/lnb.xtc``
   * - ``-o``
     - ``--output-csv``
     - Output CSV file path
     - -
     - ``-o results.csv``
   * - ``-r``
     - ``--residues``
     - Residue group definition
     - -
     - ``-r "{'DPPC': ['PO4']}"``
   * - ``-gas``
     - ``--gas-group``
     - Gas group definition
     - -
     - ``-gas O2:O2``
   * - ``-rad``
     - ``--radius``
     - Radius (Å) for density analysis
     - ``50.0``
     - ``-rad 50.0``
   * - ``-max-rad``
     - ``--max-radius``
     - Maximum radius (Å) for multi-radius density analysis
     - ``50.0``
     - ``-max-rad 50.0``
   * - ``-segments``
     - ``--number-segments``
     - Number of radius segments
     - ``5``
     - ``-segments 5``
   * - ``-mw``
     - ``--MW``
     - Molecular weight (g/mol)
     - ``14.0``
     - ``-mw 32.0``
   * - ``-p``
     - ``--parallel``
     - Enable parallel processing
     - ``False``
     - ``-p``
   * - ``-j``
     - ``--n-jobs``
     - Number of parallel jobs
     - ``-1`` (all cores)
     - ``-j 4``
   * - ``-s``
     - ``--start-frame``
     - Start frame (0-indexed)
     - ``0``
     - ``-s 0``
   * - ``-e``
     - ``--stop-frame``
     - Stop frame (exclusive)
     - All frames
     - ``-e 100``
   * - ``-t``
     - ``--step-frame``
     - Frame step size
     - ``1``
     - ``-t 5``
   * - ``-v``
     - ``--verbose``
     - Verbose output
     - ``False``
     - ``-v``
   * - ``-k``
     - ``--k-value``
     - K value for Voronoi/Sz calculation
     - ``20`` (area), ``15`` (sz)
     - ``-k 20``

Available Analysis Modules
~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Area Analysis (``area.py``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Calculates the area per lipid using Voronoi tessellation.

**Basic usage:**

.. code-block:: bash

   python analysis/area.py \
       --gro-file cases_lnb/lnb.gro \
       --xtc-file cases_lnb/lnb.xtc \
       --output-csv cases_lnb/area_results.csv \
       --residues "{'DPPC': ['PO4'], 'DAPC': ['PO4'], 'CHOL': ['ROH']}" \
       --k-value 20 \
       --max-normal-angle 140

**With parallel processing:**

.. code-block:: bash

   python analysis/area.py \
       -g cases_lnb/lnb.gro \
       -x cases_lnb/lnb.xtc \
       -o cases_lnb/area_results.csv \
       -r "{'DPPC': ['PO4'], 'DAPC': ['PO4'], 'CHOL': ['ROH']}" \
       -k 20 \
       --parallel \
       -j 4

**Single frame analysis (GRO only):**

.. code-block:: bash

   python analysis/area.py \
       -g cases_lnb/lnb.gro \
       -o cases_lnb/area_single_frame.csv \
       -r "{'DPPC': ['PO4'], 'DAPC': ['PO4'], 'CHOL': ['ROH']}"

**Key parameters:**

- ``--k-value, -k``: K value for Voronoi tessellation (default: 20)
- ``--max-normal-angle``: Maximum normal angle in degrees (default: 140)

**Note:** To generate plots, use the Python API (see "Plotting and Visualization" section below).

2. Height Analysis (``height.py``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Calculates the height of lipid molecules relative to a reference plane.

**Basic usage:**

.. code-block:: bash

   python analysis/height.py \
       --gro-file cases_lnb/lnb.gro \
       --xtc-file cases_lnb/lnb.xtc \
       --output-csv cases_lnb/height_results.csv \
       --residues "{'DPPC': (['PO4'], ['C4B', 'C4A']), 'DAPC': (['PO4'], ['C4B', 'C4A']), 'CHOL': (['ROH'], ['R5'])}" \
       --k-value 20

**Note:** The residues argument for height analysis uses tuples: ``(head_atoms, tail_atoms)`` for each lipid type.

**With frame selection:**

.. code-block:: bash

   python analysis/height.py \
       -g cases_lnb/lnb.gro \
       -x cases_lnb/lnb.xtc \
       -o cases_lnb/height_results.csv \
       -r "{'DPPC': (['PO4'], ['C4B', 'C4A']), 'DAPC': (['PO4'], ['C4B', 'C4A']), 'CHOL': (['ROH'], ['R5'])}" \
       -s 0 \
       -e 100 \
       -t 5

**Note:** To generate plots, use the Python API (see "Plotting and Visualization" section below).

3. Cluster Analysis (``cluster.py``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Analyzes lipid clustering based on distance cutoff.

**Basic usage:**

.. code-block:: bash

   python analysis/cluster.py \
       --gro-file cases_lnb/lnb.gro \
       --xtc-file cases_lnb/lnb.xtc \
       --output-csv cases_lnb/cluster_results.csv \
       --residues "{'DPPC': ['PO4'], 'DAPC': ['PO4'], 'CHOL': ['ROH']}" \
       --cutoff 8.0

**With custom cutoff and parallel processing:**

.. code-block:: bash

   python analysis/cluster.py \
       -g cases_lnb/lnb.gro \
       -x cases_lnb/lnb.xtc \
       -o cases_lnb/cluster_results.csv \
       -r "{'DPPC': ['PO4'], 'DAPC': ['PO4'], 'CHOL': ['ROH']}" \
       --cutoff 10.0 \
       --parallel \
       -j 8

**Key parameters:**

- ``--cutoff``: Cutoff distance for clustering in Angstroms (default: 8.0)

**Note:** To generate plots, use the Python API (see "Plotting and Visualization" section below).

4. Anisotropy Analysis (``anisotropy.py``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Calculates the asphericity/anisotropy of lipid molecules.

**Basic usage:**

.. code-block:: bash

   python analysis/anisotropy.py \
       --gro-file cases_lnb/lnb.gro \
       --xtc-file cases_lnb/lnb.xtc \
       --output-csv cases_lnb/anisotropy_results.csv \
       --residues "{'DPPC': ['PO4'], 'DAPC': ['PO4'], 'CHOL': ['ROH']}"

**With plotting:**

.. code-block:: bash

   python analysis/anisotropy.py \
       -g cases_lnb/lnb.gro \
       -x cases_lnb/lnb.xtc \
       -o cases_lnb/anisotropy_results.csv \
       -r "{'DPPC': ['PO4'], 'DAPC': ['PO4'], 'CHOL': ['ROH']}" \
       --plot both \
       --plot-output cases_lnb/anisotropy_plot

**Single frame analysis:**

.. code-block:: bash

   python analysis/anisotropy.py \
       -g cases_lnb/lnb.gro \
       -o cases_lnb/anisotropy_single_frame.csv \
       -r "{'DPPC': ['PO4'], 'DAPC': ['PO4'], 'CHOL': ['ROH']}"

**Key parameters:**

- ``--plot``: Plot type - 'line', 'bar', 'both', or 'none' (default: 'none')
- ``--plot-output``: Output file path for plots

**Note:** Supports multiple atoms per residue. When multiple atoms are provided, their center of geometry will be calculated.

5. Gyration Analysis (``gyration.py``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Calculates the radius of gyration for lipid molecules.

**Basic usage:**

.. code-block:: bash

   python analysis/gyration.py \
       --gro-file cases_lnb/lnb.gro \
       --xtc-file cases_lnb/lnb.xtc \
       --output-csv cases_lnb/gyration_results.csv \
       --residues "{'DPPC': ['PO4'], 'DAPC': ['PO4'], 'CHOL': ['ROH']}"

**With parallel processing:**

.. code-block:: bash

   python analysis/gyration.py \
       -g cases_lnb/lnb.gro \
       -x cases_lnb/lnb.xtc \
       -o cases_lnb/gyration_results.csv \
       -r "{'DPPC': ['PO4'], 'DAPC': ['PO4'], 'CHOL': ['ROH']}" \
       --parallel

**Note:** To generate plots, use the Python API (see "Plotting and Visualization" section below).

6. Sz Order Parameter Analysis (``sz.py``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Calculates the Sz order parameter for lipid acyl chains.

**Basic usage:**

.. code-block:: bash

   python analysis/sz.py \
       --gro-file cases_lnb/lnb.gro \
       --xtc-file cases_lnb/lnb.xtc \
       --output-csv cases_lnb/sz_results.csv \
       --residues "{'DPPC': ['PO4'], 'DAPC': ['PO4']}" \
       --chain both \
       --k-value 15

**Analyze specific chain:**

.. code-block:: bash

   python analysis/sz.py \
       -g cases_lnb/lnb.gro \
       -x cases_lnb/lnb.xtc \
       -o cases_lnb/sz_sn1_results.csv \
       -r "{'DPPC': ['PO4'], 'DAPC': ['PO4']}" \
       --chain sn1 \
       -k 15

**Key parameters:**

- ``--chain``: Chain type to analyze - 'sn1', 'sn2', or 'both' (default: 'both')
- ``--k-value, -k``: K value for Sz calculation (default: 15)

**Note:** To generate plots, use the Python API (see "Plotting and Visualization" section below).

7. Density Analysis (``density.py``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Performs density analysis with two methods: single radius (frame) or multi-radius analysis.

**Single radius analysis (frame method):**

.. code-block:: bash

   python analysis/density.py \
       --method frame \
       --gro-file cases_lnb/lnb.gro \
       --xtc-file cases_lnb/lnb.xtc \
       --output-csv cases_lnb/density_frame_results.csv \
       --residues "DPPC:PO4,CHOL:ROH" \
       --gas-group "O2:O2" \
       --radius 50.0 \
       --MW 32.0

**Multi-radius analysis:**

.. code-block:: bash

   python analysis/density.py \
       --method radius \
       --gro-file cases_lnb/lnb.gro \
       --xtc-file cases_lnb/lnb.xtc \
       --output-csv cases_lnb/density_radius_results.csv \
       --residues "DPPC:PO4,CHOL:ROH" \
       --gas-group "O2:O2" \
       --max-radius 50.0 \
       --number-segments 5 \
       --MW 32.0

**Key parameters:**

- ``--method, -M``: Analysis method - 'frame' or 'radius' (default: 'frame')
- ``--residues, -r``: Residue groups in format 'RESNAME:ATOMNAME,RESNAME:ATOMNAME'
- ``--gas-group, -gas``: Gas group in format 'GASNAME:ATOMNAME'
- ``--radius, -rad``: Radius for single radius method in Angstroms (default: 50.0)
- ``--max-radius, -max-rad``: Maximum radius for multi-radius method (default: 50.0)
- ``--number-segments, -segments``: Number of radius segments for multi-radius method (default: 5)
- ``--MW, -mw``: Molecular weight of gas molecules (default: 14.0, use 32.0 for O₂)

**Note:** To generate plots, use the Python API (see "Plotting and Visualization" section below).

Common Workflow Examples
------------------------

Analyzing a subset of frames:

.. code-block:: bash

   python analysis/area.py \
       -g cases_lnb/lnb.gro \
       -x cases_lnb/lnb.xtc \
       -o cases_lnb/area_subset.csv \
       -r "{'DPPC': ['PO4'], 'DAPC': ['PO4'], 'CHOL': ['ROH']}" \
       -s 0 \
       -e 50 \
       -t 2

Getting help for any module:

.. code-block:: bash

   python analysis/<module_name>.py --help

Plotting and Visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All analysis modules support plotting functionality. There are two ways to generate plots:

Method 1: Command-line Plotting (Anisotropy Module Only)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``anisotropy.py`` module supports direct plotting via command-line arguments:

.. code-block:: bash

   python analysis/anisotropy.py \
       -g cases_lnb/lnb.gro \
       -x cases_lnb/lnb.xtc \
       -o cases_lnb/anisotropy_results.csv \
       -r "{'DPPC': ['PO4'], 'DAPC': ['PO4'], 'CHOL': ['ROH']}" \
       --plot both \
       --plot-output cases_lnb/anisotropy_plot

**Plot options:**

- ``--plot line``: Generate line chart
- ``--plot bar``: Generate bar chart
- ``--plot both``: Generate both line and bar charts
- ``--plot-output``: Output file path (without extension for 'both', will create ``_line.png`` and ``_bar.png``)

Method 2: Python Script Plotting (All Modules)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For other modules, you can create a Python script to generate plots after running the analysis. See :doc:`figure_details` for details on using the ``make_figure.py`` tool to generate plots from CSV files.

**Available Plot Types by Module:**

.. list-table:: Plot Types by Module
   :header-rows: 1
   :widths: 20 15 15 15 15 20

   * - Module
     - Line Chart
     - Bar Chart
     - Scatter
     - Heatmap
     - Notes
   * - Area
     - ✅
     - ✅
     - ✅
     - ❌
     - Use make_figure.py
   * - Height
     - ✅
     - ✅
     - ✅
     - ❌
     - Use make_figure.py
   * - Cluster
     - ✅
     - ✅
     - ❌
     - ❌
     - Use make_figure.py
   * - Anisotropy
     - ✅
     - ✅
     - ❌
     - ❌
     - Command-line or make_figure.py
   * - Gyration
     - ✅
     - ✅
     - ❌
     - ❌
     - Use make_figure.py
   * - Sz
     - ✅
     - ✅
     - ✅
     - ❌
     - Use make_figure.py
   * - Density (single)
     - ✅
     - ✅
     - ❌
     - ❌
     - Use make_figure.py
   * - Density (multi)
     - ✅
     - ❌
     - ❌
     - ✅
     - Use make_figure.py

**Note:** For detailed information on generating figures, see :doc:`figure_details`.

Output Format
~~~~~~~~~~~~~~

All modules output CSV files with the following general structure:

- **Time column**: Simulation time (if trajectory is provided)
- **Frame column**: Frame number (if trajectory is provided)
- **Residue-specific columns**: Analysis results for each residue type specified
- **Additional columns**: Module-specific metrics

Notes
~~~~~

1. **Residue dictionary format**: Most modules require a dictionary string for residues. Use single quotes inside double quotes: ``"{'DPPC': ['PO4'], 'CHOL': ['ROH']}"``

2. **Relative paths**: All paths in the examples are relative to the project root directory. Make sure to run commands from the project root directory.

3. **Parallel processing**: Enable with ``--parallel`` flag. Use ``--n-jobs`` to control the number of parallel workers.

4. **Frame indexing**: Frame numbers are 0-indexed. Use ``--start-frame`` and ``--stop-frame`` to select a subset of frames.

5. **Single frame analysis**: Some modules (area, anisotropy, sz) support analyzing only the GRO file without XTC trajectory.

6. **Example data**: The example files in ``cases_lnb/`` use Martini 3.0 force field with DPPC:DAPC:CHOL = 5:3:2 composition and O₂ gas molecules.

Troubleshooting
~~~~~~~~~~~~~~~~

**Import errors**: Ensure the LNB-MDT environment is activated: ``conda activate LNB-MDT``

**File not found**: Check that file paths are correct relative to your current working directory

**Memory issues**: Use ``--step-frame`` to analyze fewer frames, or reduce ``--n-jobs`` for parallel processing

**Residue parsing errors**: Ensure the residues dictionary string uses proper Python syntax with matching quotes
