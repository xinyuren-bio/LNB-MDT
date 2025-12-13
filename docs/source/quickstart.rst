Quick Start Guide
=================

This guide will help you get started with LNB-MDT for lipid nanobubble analysis in 5 minutes.

Configure VMD Path
------------------

First-time use of LNB-MDT requires configuring the VMD path. VMD is used for molecular visualization and trajectory analysis.

1. **Configure VMD Path Using Command Line (Recommended)**
   
   Use the built-in VMD configuration command:

.. code:: bash

   LNB-MDT VMD --path /path/to/vmd
   
   # Example for macOS
   LNB-MDT VMD --path /Applications/VMD.app/Contents/MacOS/VMD
   
   # Example for Linux
   LNB-MDT VMD --path /usr/local/bin/vmd
   
   # Example for Windows
   LNB-MDT VMD --path "C:/Program Files/VMD/vmd.exe"

2. **Or Edit Configuration File Manually**
   
   Open the `config.ini` file in the project root directory and modify `vmd_path` to your actual VMD installation path:

.. code:: text

   # Common path examples
   Windows: C:/Program Files/VMD/vmd.exe
   macOS:   /Applications/VMD.app/Contents/vmd/vmd_MACOSXARM64
   Linux:   /usr/local/bin/vmd

3. **Verify Configuration**
   
   Check if VMD path is correctly configured:

.. code:: bash

   LNB-MDT VMD

Starting the Program
--------------------

Graphical Interface Launch
~~~~~~~~~~~~~~~~~~~~~~~~~~

Using the graphical interface is the easiest way to start using LNB-MDT:

.. code:: bash

   # Launch GUI (after installation)
   LNB-MDT UI
   
   # Or simply (UI is the default)
   LNB-MDT

After launching, you will see the LNB-MDT main interface with four functional modules: Generation, Analysis, Figure, and VMD modules. For detailed information about these modules, see the main Features section on the homepage.

Command Line Launch
~~~~~~~~~~~~~~~~~~~

For batch processing and automated analysis, you can use command-line tools:

.. code:: bash

   # View help information
   LNB-MDT --help
   
   # View help for specific analysis
   LNB-MDT AREA --help

Basic Analysis Workflow
-----------------------

Prepare Data Files
~~~~~~~~~~~~~~~~~~

LNB-MDT requires the following files for analysis:

- **GRO/TPR file**: Molecular coordinate structure file
- **XTC/TRR file**: Molecular dynamics trajectory file

The project includes example data files:
- `cases/lnb.gro` - Example GRO file  
- `cases/md.xtc` - Example XTC file

Select Analysis Type
~~~~~~~~~~~~~~~~~~~~

LNB-MDT provides comprehensive analysis capabilities for lipid nanobubbles. For detailed information about each analysis type, see the main Features section on the homepage.


Running Analysis
~~~~~~~~~~~~~~~~

Graphical Interface Execution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Load GRO and XTC files in the interface
2. Select analysis type
3. Configure parameters
4. Click the "Run" button
5. View results

Command Line Execution
^^^^^^^^^^^^^^^^^^^^^^

LNB-MDT supports simplified command-line parameter input, making it easier to use command-line tools:

**Short Parameter Alias Reference Table**

+-------------------------------+-------------------------------+----------------------------+
| Short Parameter               | Long Parameter                | Description                |
+===============================+===============================+============================+
| ``-g``                        | ``--gro-file``                | GRO file path               |
+-------------------------------+-------------------------------+----------------------------+
| ``-x``                        | ``--xtc-file``                | XTC file path               |
+-------------------------------+-------------------------------+----------------------------+
| ``-o``                        | ``--output-csv``              | Output CSV file path        |
+-------------------------------+-------------------------------+----------------------------+
| ``-r``                        | ``--residues``                | Residue group definition    |
+-------------------------------+-------------------------------+----------------------------+
| ``-a``                        | ``--gas-group``               | Gas group definition        |
+-------------------------------+-------------------------------+----------------------------+
| ``-m``                        | ``--MW``                      | Molecular weight (g/mol)    |
+-------------------------------+-------------------------------+----------------------------+
| ``-R``                        | ``--radius``                  | Radius (Å)                  |
+-------------------------------+-------------------------------+----------------------------+
| ``-p``                        | ``--parallel``                | Enable parallel processing  |
+-------------------------------+-------------------------------+----------------------------+
| ``-j``                        | ``--n-jobs``                  | Number of parallel jobs     |
+-------------------------------+-------------------------------+----------------------------+
| ``-s``                        | ``--start-frame``             | Start frame                 |
+-------------------------------+-------------------------------+----------------------------+
| ``-e``                        | ``--stop-frame``              | Stop frame                  |
+-------------------------------+-------------------------------+----------------------------+
| ``-t``                        | ``--step-frame``              | Frame step                  |
+-------------------------------+-------------------------------+----------------------------+
| ``-v``                        | ``--verbose``                 | Verbose output              |
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

**Using LNB-MDT command (recommended):**

.. code-block:: bash

   # Gyration analysis example with full parameters
   LNB-MDT GYRATION \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/gyration_results.csv \
     --residues "{'DPPC': ['PO4']}" \
     --parallel \
     --verbose

**Simplified approach with short parameters:**

.. code-block:: bash

   # Using short parameters and simple format
   LNB-MDT GYRATION \
     -g cases/lnb.gro \
     -x cases/md.xtc \
     -o results/gyration_results.csv \
     -r DPPC:PO4 \
     -p \
     -v
   
   # Or use test mode for quick testing
   LNB-MDT GYRATION -test

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

Next Steps
----------

Congratulations! You have successfully completed the LNB-MDT quick start!

What you can do next:

- Learn advanced usage of :doc:`analysis_modules`  
- Check :doc:`api_reference` for API details
