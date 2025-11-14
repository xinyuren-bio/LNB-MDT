Analysis Modules Details
========================

LNB-MDT provides rich molecular dynamics analysis modules, each targeting specific physical properties for analysis.

Command Line Parameters Details
===============================

LNB-MDT provides rich command-line parameters with support for simplified input formats and short parameter aliases.

Parameter Reference Table
-------------------------

All command-line parameters have corresponding short aliases to make the command line more concise:

.. list-table:: Parameter Reference Table
   :header-rows: 1
   :widths: 8 15 25 20 32

   * - Short Parameter
     - Long Parameter
     - Description
     - Default Value
     - Example
   * - ``-g``
     - ``--gro-file``
     - GRO/TPR file path
     - -
     - ``-g cases/lnb.gro``
   * - ``-x``
     - ``--xtc-file``
     - XTC/TRR file path
     - -
     - ``-x cases/md.xtc``
   * - ``-o``
     - ``--output-csv``
     - Output CSV file path
     - ``cases/csv/results.csv``
     - ``-o results.csv``
   * - ``-r``
     - ``--residues``
     - Residue group definition
     - ``DPPC:PO4,CHOL:ROH``
     - ``-r DPPC:PO4``
   * - ``-gas``
     - ``--gas-group``
     - Gas group definition
     - ``N2:N2``
     - ``-gas N2:N2``
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
     - ``-mw 14.0``
   * - ``-p``
     - ``--parallel``
     - Enable parallel processing
     - ``False``
     - ``-p``
   * - ``-j``
     - ``--n-jobs``
     - Number of parallel jobs
     - ``2``
     - ``-j 4``
   * - ``-s``
     - ``--start-frame``
     - Start frame
     - ``0``
     - ``-s 0``
   * - ``-e``
     - ``--stop-frame``
     - Stop frame
     - ``All frames``
     - ``-e 100``
   * - ``-t``
     - ``--step-frame``
     - Frame step
     - ``1``
     - ``-t 5``
   * - ``-v``
     - ``--verbose``
     - Verbose output
     - ``False``
     - ``-v``
   * - ``-k``
     - ``--k-value``
     - k value
     - ``20``
     - ``-k 20``
   * - ``-M``
     - ``--method``
     - Calculation method
     - ``mean``
     - ``-M mean``
   * - ``-T``
     - ``--threshold``
     - Threshold
     - ``0.5``
     - ``-T 0.5``
   * - ``-P``
     - ``--plot-type``
     - Plot type
     - ``all``
     - ``-P line``
   * - ``-d``
     - ``--plot-dir``
     - Plot directory
     - ``plots/``
     - ``-d plots/``

Simplified Format Description
-----------------------------

The residues and gas-group parameters now support more intuitive input formats:

Detailed Module Description
----------------------------

Area Analysis (area.py)
~~~~~~~~~~~~~~~~~~~~~~~

**Function Description**
Uses Voronoi tessellation method to calculate area distribution of lipid molecules.

**Algorithm Principle**
- Construct Voronoi diagram
- Calculate Voronoi area for each molecule
- Analyze area distribution and changes

**Key Parameters**

k-value *number* (``-k``)
    k-value for Voronoi tessellation, default value is ``20``

max-normal-angle *angle*
    Maximum normal angle in degrees, default value is ``140``

residues *residue-definition* (``-r``)
    Residue group definition, specifying molecular types and atoms to analyze. Supports simplified format like ``DPPC:PO4,CHOL:ROH`` or dictionary format ``"{'DPPC': ['PO4']}"``

xtc-file *file-path* (``-x``)
    XTC file path (optional). If not provided, only GRO file will be analyzed (single frame)

**Usage Example**

.. code-block:: bash

   python analysis/area.py \
     -g cases/lnb.gro \
     -x cases/md.xtc \
     -o results/area_results.csv \
     -r DPPC:PO4,CHOL:ROH \
     -k 20 \
     --max-normal-angle 140 \
     -p \
     -v

**Output Results**
- Voronoi area for each molecule
- Area distribution statistics
- Can be used to analyze membrane density and packing

Height Analysis (height.py)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Function Description**
Analyzes height distribution and membrane thickness of lipid molecules.

**Algorithm Principle**
- Calculate molecular positions in Z direction
- Analyze height distribution
- Calculate membrane thickness and surface roughness

**Key Parameters**

k-value *number* (``-k``)
    k-value for height calculation, default value is ``20``

residues *residue-definition* (``-r``)
    Residue group definition, specifying molecular types and atoms to analyze. Supports multiple atom groups using tuple format: ``"{'DPPC': (['PO4'], ['C4B', 'C4A']), 'CHOL':(['ROH'], ['R5'])}"``

**Usage Example**

.. code-block:: bash

   python analysis/height.py \
     -g cases/lnb.gro \
     -x cases/md.xtc \
     -o results/height_results.csv \
     -r "{'DPPC': (['PO4'], ['C4B', 'C4A']), 'CHOL':(['ROH'], ['R5'])}" \
     -k 20 \
     -p \
     -v

**Output Results**
- Height values for each molecule
- Height distribution statistics
- Membrane thickness analysis

Cluster Analysis (cluster.py)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Function Description**
Analyzes aggregation behavior and clustering patterns of lipid molecules.

**Algorithm Principle**
- Distance-based clustering algorithm
- Identify molecular aggregations
- Analyze cluster size and distribution

**Key Parameters**

cutoff *distance*
    Clustering cutoff distance in Angstroms, default value is ``8.0``

residues *residue-definition* (``-r``)
    Residue group definition, specifying molecular types and atoms to analyze. Supports simplified format like ``DPPC:PO4,CHOL:ROH``

**Usage Example**

.. code-block:: bash

   python analysis/cluster.py \
     -g cases/lnb.gro \
     -x cases/md.xtc \
     -o results/cluster_results.csv \
     -r DPPC:PO4,CHOL:ROH \
     --cutoff 8.0 \
     -p \
     -v

**Output Results**
- Cluster size distribution
- Cluster count statistics
- Aggregation behavior analysis

Anisotropy Analysis (anisotropy.py)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Function Description**
Calculates anisotropy parameters of molecular orientation.

**Algorithm Principle**
- Calculate molecular orientation vectors
- Analyze orientation distribution
- Calculate anisotropy parameters

**Key Parameters**

residues *residue-definition* (``-r``)
    Residue group definition, specifying molecular types and atoms to analyze. Supports multiple atoms per residue: ``"{'DPPC': ['PO4', 'C1', 'C2'], 'CHOL': ['ROH']}"``. When multiple atoms are provided, their center of geometry will be calculated.

xtc-file *file-path* (``-x``)
    XTC file path (optional). If not provided, only GRO file will be analyzed (single frame)

**Usage Example**

.. code-block:: bash

   python analysis/anisotropy.py \
     -g cases/lnb.gro \
     -x cases/md.xtc \
     -o results/anisotropy_results.csv \
     -r DPPC:PO4,CHOL:ROH \
     -p \
     -v

**Output Results**
- Anisotropy parameters
- Orientation distribution statistics
- Molecular alignment analysis

Gyration Analysis (gyration.py)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Function Description**
Calculates molecular radius of gyration, reflecting molecular compactness.

**Algorithm Principle**
- Calculate molecular center of mass
- Calculate radius of gyration
- Analyze molecular shape changes

**Key Parameters**

residues *residue-definition* (``-r``)
    Residue group definition, specifying molecular types and atoms to analyze. Supports simplified format like ``DPPC:PO4,CHOL:ROH`` or dictionary format ``"{'DPPC': ['PO4']}"``

**Usage Example**

.. code-block:: bash

   python analysis/gyration.py \
     -g cases/lnb.gro \
     -x cases/md.xtc \
     -o results/gyration_results.csv \
     -r DPPC:PO4,CHOL:ROH \
     -p \
     -v

**Output Results**
- Radius of gyration values
- Shape change analysis
- Molecular compactness

Sz Order Parameter Analysis (sz.py)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Function Description**
Calculates Sz order parameter of lipid chains, reflecting chain ordering degree.

**Algorithm Principle**
- Calculate chain orientation vectors
- Calculate Sz order parameter
- Analyze chain ordering

**Key Parameters**

chain *chain-type*
    Chain type, options are ``sn1``, ``sn2`` or ``both``, default value is ``both``

k-value *number* (``-k``)
    k-value for Sz calculation, default value is ``15``

residues *residue-definition* (``-r``)
    Residue group definition, specifying molecular types and atoms to analyze. Supports simplified format like ``DPPC:PO4,DUPC:PO4``

xtc-file *file-path* (``-x``)
    XTC file path (optional). If not provided, only GRO file will be used

**Usage Example**

.. code-block:: bash

   python analysis/sz.py \
     -g cases/lnb.gro \
     -x cases/md.xtc \
     -o results/sz_results.csv \
     -r DPPC:PO4,DUPC:PO4 \
     --chain sn1 \
     -k 15 \
     -p \
     -v

**Output Results**
- Sz order parameter values
- Chain ordering analysis
- Phase transition behavior

Density Analysis (density.py)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Function Description**
Analyzes gas density distribution in lipid nanobubbles. Supports two analysis methods: single-radius density analysis (frame method) and multi-radius density analysis (radius method).

**Algorithm Principle**
- Calculate gas molecule density within specified radius
- Analyze density changes over time (frame method)
- Analyze density distribution across different radii (radius method)
- Support for multiple gas types

**Key Parameters**

method *analysis-type* (``-M``)
    Analysis method, options are ``frame`` (single radius) or ``radius`` (multi-radius), default value is ``frame``

residues *residue-definition* (``-r``)
    Residue group definition for lipid molecules. Format: ``DPPC:PO4,CHOL:ROH``

gas-group *gas-definition* (``-gas``)
    Gas group definition. Format: ``N2:N2`` or ``N2:N2,O2:O2`` for multiple gases

radius *distance* (``-rad``)
    Radius for density calculation in Angstroms (used for single radius method), default value is ``50.0``

max-radius *distance* (``-max-rad``)
    Maximum radius for density calculation in Angstroms (used for multi-radius method), default value is ``50.0``

number-segments *number* (``-segments``)
    Number of radius segments for multi-radius analysis, default value is ``5``

MW *molecular-weight* (``-mw``)
    Molecular weight of gas molecules in g/mol, default value is ``14.0``

**Usage Example**

Single-radius analysis (frame method):
.. code-block:: bash

   python analysis/density.py \
     -g cases/lnb.gro \
     -x cases/md.xtc \
     -o results/density_frame.csv \
     -r DPPC:PO4,CHOL:ROH \
     -gas N2:N2 \
     -rad 50.0 \
     -mw 14.0 \
     -M frame \
     -p \
     -v

Multi-radius analysis (radius method):
.. code-block:: bash

   python analysis/density.py \
     -g cases/lnb.gro \
     -x cases/md.xtc \
     -o results/density_radius.csv \
     -r DPPC:PO4,CHOL:ROH \
     -gas N2:N2 \
     -max-rad 50.0 \
     -segments 5 \
     -mw 14.0 \
     -M radius \
     -p \
     -v

**Output Results**
- CSV file containing density values for each frame (frame method)
- CSV file containing density values for each radius segment (radius method)
- Can be used to analyze gas distribution and bubble dynamics
- Supports visualization of analysis results

Parameter Optimization Recommendations
----------------------------------------

k-value Selection
~~~~~~~~~~~~~~~~~

.. raw:: html

   <div style="background-color: #e3f2fd; padding: 15px; border-radius: 8px; border-left: 4px solid #2196f3;">

**k-value Selection Principles:**

- **Small systems**: k = 10-15
- **Medium systems**: k = 15-25  
- **Large systems**: k = 25-35
- **High density**: Increase k-value
- **Low density**: Decrease k-value

**Optimization Method:**
Use the k-value optimizer from the machine learning module:

.. code:: python

   from machine_learning import KValueOptimizer
   optimizer = KValueOptimizer('area')
   best_k = optimizer.optimize()

   </div>

Cutoff Distance Selection
~~~~~~~~~~~~~~~~~~~~~~~~

.. raw:: html

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800;">

**Cutoff Distance Selection:**

- **Cluster analysis**: 8-12 Å
- **Interactions**: 5-8 Å
- **Long-range interactions**: 12-20 Å

**Selection Criteria:**
- Molecular size
- Interaction strength
- System density


Parallel Processing Optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px; border-left: 4px solid #4caf50;">

**Parallel Processing Recommendations:**

- **CPU cores**: Use `--n-jobs -1` for automatic detection
- **Memory considerations**: Reduce parallel jobs for large systems
- **I/O limitations**: SSD drives allow more parallel jobs

**Performance Optimization:**
- Use SSD storage for trajectory files
- Increase system memory
- Optimize network file systems
