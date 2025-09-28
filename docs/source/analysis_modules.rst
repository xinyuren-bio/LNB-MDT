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
     - GRO file path
     - -
     - ``-g cases/lnb.gro``
   * - ``-x``
     - ``--xtc-file``
     - XTC file path
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
   * - ``-a``
     - ``--gas-group``
     - Gas group definition
     - ``N2:N2``
     - ``-a N2:N2``
   * - ``-m``
     - ``--MW``
     - Molecular weight (g/mol)
     - ``14``
     - ``-m 14``
   * - ``-R``
     - ``--radius``
     - Radius (Å)
     - ``50``
     - ``-R 50``
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

Basic Format
~~~~~~~~~~~~

**Simple format (recommended):**
.. code-block:: bash

   # Basic format: RESIDUE:ATOM
   -r DPPC:PO4,CHOL:ROH
   -a N2:N2
   
   # Multiple residues/gases
   -r DPPC:PO4,DUPC:PO4,CHOL:ROH
   -a N2:N2,O2:O2

**Multi-atom format:**
.. code-block:: bash

   # Multi-atom: RESIDUE:ATOM1+ATOM2
   -r DPPC:PO4+GLY,CHOL:ROH
   -r DPPC:PO4+GLY+CH2,CHOL:ROH

**Name-only format:**
.. code-block:: bash

   # Only residue/gas name (atom name same as name)
   -r DPPC
   -a N2

Traditional Format
~~~~~~~~~~~~~~~~~~

**Dictionary string format (still supported):**
.. code-block:: bash

   # Traditional dictionary format
   -r "{'DPPC': ['PO4'], 'CHOL': ['ROH']}"
   -a "{'N2': ['N2']}"

Module Overview
---------------

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; margin: 20px 0;">

   <div style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 20px; border-radius: 10px;">
   <h3 style="margin-top: 0;">📐 PCA Analysis</h3>
   <p>Principal Component Analysis for studying molecular conformational changes</p>
   <ul style="margin-bottom: 0;">
   <li>Dimensionality reduction</li>
   <li>Conformational changes</li>
   <li>Motion patterns</li>
   </ul>
   </div>

   <div style="background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%); color: white; padding: 20px; border-radius: 10px;">
   <h3 style="margin-top: 0;">📏 Area Analysis</h3>
   <p>APL area calculation</p>
   <ul style="margin-bottom: 0;">
   <li>Molecular area</li>
   <li>Density distribution</li>
   <li>Packing efficiency</li>
   </ul>
   </div>

   <div style="background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%); color: white; padding: 20px; border-radius: 10px;">
   <h3 style="margin-top: 0;">🌊 Curvature Analysis</h3>
   <p>Membrane curvature calculation (mean/Gaussian)</p>
   <ul style="margin-bottom: 0;">
   <li>Mean curvature</li>
   <li>Gaussian curvature</li>
   <li>Membrane deformation</li>
   </ul>
   </div>

   <div style="background: linear-gradient(135deg, #43e97b 0%, #38f9d7 100%); color: white; padding: 20px; border-radius: 10px;">
   <h3 style="margin-top: 0;">📊 Height Analysis</h3>
   <p>Molecular height distribution analysis</p>
   <ul style="margin-bottom: 0;">
   <li>Z-coordinate distribution</li>
   <li>Membrane thickness</li>
   <li>Surface roughness</li>
   </ul>
   </div>

   <div style="background: linear-gradient(135deg, #fa709a 0%, #fee140 100%); color: white; padding: 20px; border-radius: 10px;">
   <h3 style="margin-top: 0;">🔗 Cluster Analysis</h3>
   <p>Molecular aggregation behavior analysis</p>
   <ul style="margin-bottom: 0;">
   <li>Aggregation patterns</li>
   <li>Cluster size</li>
   <li>Interactions</li>
   </ul>
   </div>

   <div style="background: linear-gradient(135deg, #a8edea 0%, #fed6e3 100%); color: #333; padding: 20px; border-radius: 10px;">
   <h3 style="margin-top: 0;">🎯 Anisotropy Analysis</h3>
   <p>Molecular orientation anisotropy calculation</p>
   <ul style="margin-bottom: 0;">
   <li>Orientation distribution</li>
   <li>Order parameter</li>
   <li>Molecular alignment</li>
   </ul>
   </div>

   </div>

Detailed Module Description
----------------------------

PCA Analysis (pca.py)
~~~~~~~~~~~~~~~~~~~~~

**Function Description**
Principal Component Analysis is used to study conformational changes and motion patterns of lipid molecules.

**Algorithm Principle**
- Perform principal component analysis on molecular coordinates
- Extract main motion patterns
- Reduce dimensionality to principal component space

**Key Parameters**

residues *residue-definition*
    Residue group definition, specifying molecular types and atoms to analyze. Supports simplified format like ``DPPC:PO4,CHOL:ROH``

n_components *number*
    Number of principal components, default value is ``3``

start_frame *frame-number*
    Start frame, default value is ``0``

stop_frame *frame-number*
    Stop frame, default value is ``-1`` (means analyze to the end)

step_frame *frame-step*
    Frame step, default value is ``1``

**Usage Example**

.. code-block:: bash

   python analysis/pca.py \
     -g cases/lnb.gro \
     -x cases/md.xtc \
     -o results/pca_results.csv \
     -r DPPC:PO4,CHOL:ROH \
     --n-components 3 \
     -p \
     -v

**Output Results**
- CSV file contains principal component values for each frame
- Can be used to analyze molecular conformational change trends
- Supports visualization of analysis results

Area Analysis (area.py)
~~~~~~~~~~~~~~~~~~~~~~~

**Function Description**
Uses Voronoi tessellation method to calculate area distribution of lipid molecules.

**Algorithm Principle**
- Construct Voronoi diagram
- Calculate Voronoi area for each molecule
- Analyze area distribution and changes

**Key Parameters**

k-value *number*
    k-value for Voronoi tessellation, default value is ``20``

max-normal-angle *angle*
    Maximum normal angle, default value is ``140`` degrees

residues *residue-definition*
    Residue group definition, specifying molecular types and atoms to analyze

**Usage Example**

.. code-block:: bash

   python analysis/area.py \
     -g cases/lnb.gro \
     -x cases/md.xtc \
     -o results/area_results.csv \
     -r DPPC:PO4 \
     -k 20 \
     --max-normal-angle 140 \
     -p \
     -v

**Output Results**
- Voronoi area for each molecule
- Area distribution statistics
- Can be used to analyze membrane density and packing

Curvature Analysis (curvature.py)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Function Description**
Calculates mean curvature and Gaussian curvature of lipid membranes.

**Algorithm Principle**
- Based on local surface fitting
- Calculate curvature tensor
- Extract mean curvature and Gaussian curvature

**Key Parameters**

method *curvature-type*
    Curvature type, options are ``mean`` or ``gaussian``, default value is ``mean``

k-value *number*
    k-value for curvature calculation, default value is ``20``

residues *residue-definition*
    Residue group definition, specifying molecular types and atoms to analyze

**Usage Example**

.. code-block:: bash

   python analysis/curvature.py \
     -g cases/lnb.gro \
     -x cases/md.xtc \
     -o results/curvature_results.csv \
     -r DPPC:PO4 \
     -k 20 \
     -M mean \
     -p \
     -v

**Output Results**
- Curvature values for each molecule
- Curvature distribution statistics
- Can be used to analyze membrane deformation and stability

高度分析 (height.py)
~~~~~~~~~~~~~~~~~~~~~

**功能描述**
分析脂质分子的高度分布和膜厚度。

**算法原理**
- 计算分子在Z方向的位置
- 分析高度分布
- 计算膜厚度和表面粗糙度

**关键参数**

k-value *number*
    高度计算的k值，默认值为 ``20``

residues *residue-definition*
    残基组定义，指定要分析的分子类型和原子，支持多组原子

**使用示例**

.. code-block:: bash

   python analysis/height.py \
     -g cases/lnb.gro \
     -x cases/md.xtc \
     -o results/height_results.csv \
     -r DPPC:PO4,CHOL:ROH \
     -k 20 \
     -p \
     -v

**输出结果**
- 每个分子的高度值
- 高度分布统计
- 膜厚度分析

聚类分析 (cluster.py)
~~~~~~~~~~~~~~~~~~~~~~

**功能描述**
分析脂质分子的聚集行为和聚类模式。

**算法原理**
- 基于距离的聚类算法
- 识别分子聚集
- 分析聚类大小和分布

**关键参数**

cutoff *distance*
    聚类截止距离，默认值为 ``8.0`` 埃

residues *residue-definition*
    残基组定义，指定要分析的分子类型和原子

**使用示例**

.. code-block:: bash

   python analysis/cluster.py \
     -g cases/lnb.gro \
     -x cases/md.xtc \
     -o results/cluster_results.csv \
     -r DPPC:PO4,CHOL:ROH \
     --cutoff 8.0 \
     -p \
     -v

**输出结果**
- 聚类大小分布
- 聚类数量统计
- 聚集行为分析

各向异性分析 (anisotropy.py)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**功能描述**
计算分子取向的各向异性参数。

**算法原理**
- 计算分子取向向量
- 分析取向分布
- 计算各向异性参数

**关键参数**

residues *residue-definition*
    残基组定义，指定要分析的分子类型和原子

**使用示例**

.. code-block:: bash

   python analysis/anisotropy.py \
     -g cases/lnb.gro \
     -x cases/md.xtc \
     -o results/anisotropy_results.csv \
     -r DPPC:PO4,CHOL:ROH \
     -p \
     -v

**输出结果**
- 各向异性参数
- 取向分布统计
- 分子排列分析

回转半径分析 (gyration.py)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**功能描述**
计算分子的回转半径，反映分子的紧凑程度。

**算法原理**
- 计算分子质心
- 计算回转半径
- 分析分子形状变化

**关键参数**

residues *residue-definition*
    残基组定义，指定要分析的分子类型和原子

**使用示例**

.. code-block:: bash

   python analysis/gyration.py \
     -g cases/lnb.gro \
     -x cases/md.xtc \
     -o results/gyration_results.csv \
     -r DPPC:PO4,CHOL:ROH \
     -p \
     -v

**输出结果**
- 回转半径值
- 形状变化分析
- 分子紧凑度

Sz序参数分析 (sz.py)
~~~~~~~~~~~~~~~~~~~~~

**功能描述**
计算脂质链的Sz序参数，反映链的有序程度。

**算法原理**
- 计算链取向向量
- 计算Sz序参数
- 分析链有序性

**关键参数**

chain *chain-type*
    链类型，可选值为 ``sn1``、``sn2`` 或 ``both``

k-value *number*
    Sz计算的k值，默认值为 ``15``

residues *residue-definition*
    残基组定义，指定要分析的分子类型和原子

**使用示例**

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

**输出结果**
- Sz序参数值
- 链有序性分析
- 相变行为

N-聚类分析 (n_cluster.py)
~~~~~~~~~~~~~~~~~~~~~~~~~~

**功能描述**
统计聚类数量，分析聚集模式。

**算法原理**
- 基于距离的聚类
- 统计聚类数量
- 分析聚集模式

**关键参数**

cutoff *distance*
    聚类截止距离，默认值为 ``12.0`` 埃

n-cutoff *number*
    最小聚类大小阈值，默认值为 ``10``

residues *residue-definition*
    残基组定义，指定要分析的分子类型和原子

**使用示例**

.. code-block:: bash

   python analysis/n_cluster.py \
     -g cases/lnb.gro \
     -x cases/md.xtc \
     -o results/ncluster_results.csv \
     -r DAPC:GL1+GL2,DPPC:PO4 \
     --cutoff 12.0 \
     --n-cutoff 10 \
     -p \
     -v

**输出结果**
- 聚类数量统计
- 聚集模式分析
- 相互作用强度

径向分布分析 (rad.py)
~~~~~~~~~~~~~~~~~~~~~~

**功能描述**
计算径向分布函数，分析分子间的距离分布。

**算法原理**
- 计算分子间距离
- 构建径向分布函数
- 分析相互作用

**关键参数**

n-circle *number*
    径向分析的同心圆数量，默认值为 ``50``

residues *residue-definition*
    残基组定义，指定要分析的分子类型和原子

**使用示例**

.. code-block:: bash

   python analysis/rad.py \
     -g cases/lnb.gro \
     --output-excel results/radial_distribution.xlsx \
     -r DPPC:NC3,CHOL:ROH \
     --n-circle 50

**输出结果**
- Excel文件包含径向分布数据
- 距离分布统计
- 相互作用分析

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
