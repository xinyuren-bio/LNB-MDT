Quick Start Guide
=================

This guide will help you get started with LNB-MDT for lipid nanobubble analysis in 5 minutes.

Configure VMD Path
------------------

First-time use of LNB-MDT requires configuring the VMD path. VMD is used for molecular visualization and trajectory analysis.

1. **Edit Configuration File**
   
   Open the `config.ini` file in the project root directory and modify `vmd_path` to your actual VMD installation path:

.. code:: text

   # Common path examples
   Windows: C:/Program Files/VMD/vmd.exe
   macOS:   /Applications/VMD.app/Contents/vmd/vmd_MACOSXARM64
   Linux:   /usr/local/bin/vmd

2. **Save and Restart**
   
   Save the configuration file and restart the program.

Starting the Program
--------------------

Graphical Interface Launch
~~~~~~~~~~~~~~~~~~~~~~~~~~

Using the graphical interface is the easiest way to start using LNB-MDT:

.. code:: python

   # Activate environment
   conda activate LNB-MDT
   
   # Start main program
   python main.py

After launching, you will see the LNB-MDT main interface with the following functional modules:

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin: 20px 0;">

   <div style="background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%); color: white; padding: 10px; border-radius: 8px; text-align: center;">
   <h4 style="margin-center: 0;">🧬 Generation Module</h4>
   </div>

   <div style="background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%); color: white; padding: 10px; border-radius: 8px; text-align: center;">
   <h4 style="margin-center: 0;">📊 Analysis Module</h4>
   </div>

   <div style="background: linear-gradient(135deg, #43e97b 0%, #38f9d7 100%); color: white; padding: 10px; border-radius: 8px; text-align: center;">
   <h4 style="margin-center: 0;">📈 Figure Module</h4>
   </div>

   <div style="background: linear-gradient(135deg, #fa709a 0%, #fee140 100%); color: white; padding: 10px; border-radius: 8px; text-align: center;">
   <h4 style="margin-center: 0;">🔧 VMD Module</h4>
   </div>

   </div>

Command Line Launch
~~~~~~~~~~~~~~~~~~~

For batch processing and automated analysis, you can use command-line tools:

.. code:: python

   # Activate environment
   conda activate LNB-MDT
   
   # View help information
   python analysis/pca.py --help

Basic Analysis Workflow
-----------------------

Prepare Data Files
~~~~~~~~~~~~~~~~~~

LNB-MDT requires the following files for analysis:

- **GRO file**: Molecular topology structure file
- **XTC file**: Molecular dynamics trajectory file

The project includes sample data files:
- `cases/lnb.gro` - Sample topology file  
- `cases/md.xtc` - Sample trajectory file

Select Analysis Type
~~~~~~~~~~~~~~~~~~~~

LNB-MDT provides comprehensive analysis capabilities for lipid nanobubbles. For detailed information about each analysis type, see the main Features section on the homepage.


运行分析
~~~~~~~~

图形界面运行
^^^^^^^^^^^^

1. 在界面中加载GRO和XTC文件
2. 选择分析类型
3. 配置参数
4. 点击"运行"按钮
5. 查看结果

命令行运行
^^^^^^^^^^

LNB-MDT支持简化的命令行参数输入，让您更轻松地使用命令行工具：

**短参数别名对照表**

+-------------------------------+-------------------------------+----------------------------+
| 短参数                        | 长参数                        | 说明                       |
+===============================+===============================+============================+
| ``-g``                        | ``--gro-file``                | GRO文件路径                |
+-------------------------------+-------------------------------+----------------------------+
| ``-x``                        | ``--xtc-file``                | XTC文件路径                |
+-------------------------------+-------------------------------+----------------------------+
| ``-o``                        | ``--output-csv``              | 输出CSV文件路径            |
+-------------------------------+-------------------------------+----------------------------+
| ``-r``                        | ``--residues``                | 残基组定义                 |
+-------------------------------+-------------------------------+----------------------------+
| ``-a``                        | ``--gas-group``               | 气体组定义                 |
+-------------------------------+-------------------------------+----------------------------+
| ``-m``                        | ``--MW``                      | 分子量 (g/mol)             |
+-------------------------------+-------------------------------+----------------------------+
| ``-R``                        | ``--radius``                  | 半径 (Å)                   |
+-------------------------------+-------------------------------+----------------------------+
| ``-p``                        | ``--parallel``                | 启用并行处理               |
+-------------------------------+-------------------------------+----------------------------+
| ``-j``                        | ``--n-jobs``                  | 并行任务数                 |
+-------------------------------+-------------------------------+----------------------------+
| ``-s``                        | ``--start-frame``             | 起始帧                     |
+-------------------------------+-------------------------------+----------------------------+
| ``-e``                        | ``--stop-frame``              | 结束帧                     |
+-------------------------------+-------------------------------+----------------------------+
| ``-t``                        | ``--step-frame``              | 帧步长                     |
+-------------------------------+-------------------------------+----------------------------+
| ``-v``                        | ``--verbose``                 | 详细输出                   |
+-------------------------------+-------------------------------+----------------------------+

**Simplified residues and gas-group formats**

.. code-block:: python

   # Simple format (recommended)
   -r DPPC:PO4,CHOL:ROH
   -a N2:N2
   
   # Multi-atom format
   -r DPPC:PO4+GLY,CHOL:ROH
   
   # Traditional dictionary format (still supported)
   -r "{'DPPC': ['PO4'], 'CHOL': ['ROH']}"

**Traditional approach (still supported):**

.. code-block:: python

   # PCA analysis example
   python analysis/pca.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/pca_results.csv \
     --residues "{'DPPC': ['PO4']}" \
     --parallel \
     --verbose

**New simplified approach (recommended):**

.. code-block:: python

   # Using short parameters and simple format
   python analysis/pca.py \
     -g cases/lnb.gro \
     -x cases/md.xtc \
     -o results/pca_results.csv \
     -r DPPC:PO4 \
     -p \
     -v

Viewing Results
~~~~~~~~~~~~~~~

After analysis completion, LNB-MDT generates the following outputs:

- **CSV files**: Numerical data containing analysis results
- **Charts**: Visualization of analysis results  
- **Logs**: Detailed information about the analysis process

Result interpretation:

- View numerical results in CSV files
- Use chart module to visualize data
- Combine with VMD for molecular visualization

Practical Examples
-------------------

PCA Analysis
~~~~~~~~~~~~

Analyze conformational changes of lipid molecules:

**Traditional approach:**

.. code-block:: python

   python analysis/pca.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/pca_test.csv \
     --residues "{'DPPC': ['PO4'], 'CHOL': ['ROH']}" \
     --start-frame 0 \
     --stop-frame 100 \
     --parallel \
     --verbose

**Simplified approach:**

.. code-block:: python

   python analysis/pca.py \
     -g cases/lnb.gro \
     -x cases/md.xtc \
     -o results/pca_test.csv \
     -r DPPC:PO4,CHOL:ROH \
     -s 0 \
     -e 100 \
     -p \
     -v

Area Analysis
~~~~~~~~~~~~~

Calculate Voronoi tessellation area of lipid molecules:

**Traditional approach:**

.. code-block:: python

   python analysis/area.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/area_test.csv \
     --residues "{'DPPC': ['PO4']}" \
     --k-value 20 \
     --max-normal-angle 140 \
     --parallel \
     --verbose

**Simplified approach:**

.. code-block:: python

   python analysis/area.py \
     -g cases/lnb.gro \
     -x cases/md.xtc \
     -o results/area_test.csv \
     -r DPPC:PO4 \
     -k 20 \
     --max-normal-angle 140 \
     -p \
     -v

Curvature Analysis
~~~~~~~~~~~~~~~~~~

Calculate curvature properties of lipid membranes:

**Traditional approach:**

.. code-block:: python

   python analysis/curvature.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/curvature_test.csv \
     --residues "{'DPPC': ['PO4']}" \
     --k-value 20 \
     --method mean \
     --parallel \
     --verbose

**Simplified approach:**

.. code-block:: python

   python analysis/curvature.py \
     -g cases/lnb.gro \
     -x cases/md.xtc \
     -o results/curvature_test.csv \
     -r DPPC:PO4 \
     -k 20 \
     -M mean \
     -p \
     -v

Density Analysis
~~~~~~~~~~~~~~~~

Analyze gas density changes over time in bubbles:

**Simplified approach (recommended):**

.. code-block:: python

   python analysis/densitywithframe.py \
     -g cases/lnb.gro \
     -x cases/md.xtc \
     -o results/density_test.csv \
     -r DPPC:PO4,CHOL:ROH \
     -a N2:N2 \
     -m 14 \
     -R 50 \
     -p \
     -v


VMD Integration
---------------

LNB-MDT supports seamless integration with VMD for molecular visualization and trajectory analysis.

VMD Path Configuration
~~~~~~~~~~~~~~~~~~~~~~

First-time use requires configuring the VMD path:

1. **Find VMD Installation Path**

.. code:: text

   Windows: Usually at C:/Program Files/VMD/vmd.exe
   macOS:   Usually at /Applications/VMD.app/Contents/vmd/vmd_MACOSXARM64
   Linux:   Usually at /usr/local/bin/vmd

2. **Edit Configuration File**
   
   Open the `config.ini` file in the project root directory and modify `vmd_path` to your actual VMD installation path:

.. code:: ini

   [VMD]
   vmd_path = /Applications/VMD.app/Contents/vmd/vmd_MACOSXARM64

3. **Verify Configuration**
   
   Save the configuration file and restart the LNB-MDT program.

Starting VMD
~~~~~~~~~~~~

Graphical interface launch:

1. Click the "Start VMD" button
2. Wait for VMD to start
3. Drag CSV files to the VMD window
4. Select molecules for visualization

Command line launch:

.. code:: python

   # Start VMD
   python -c "from modules.vmd_control import VMDTcp; vmd = VMDTcp(); vmd.start()"

Visualization Operations
~~~~~~~~~~~~~~~~~~~~~~~~

Operation steps:

1. Load analysis results in LNB-MDT
2. Select frames and molecules to visualize
3. VMD automatically jumps to the corresponding frame
4. Highlight selected molecules
5. Adjust visualization parameters

Next Steps
----------

Congratulations! You have successfully completed the LNB-MDT quick start!

What you can do next:

- Learn advanced usage of :doc:`analysis_modules`  
- Check :doc:`api_reference` for API details
