LNB-MDT: Lipid NanoBubble Molecular Dynamics Toolbox
====================================================

**LNB-MDT** (Lipid NanoBubble Molecular Dynamics Toolbox) is a comprehensive toolbox designed for molecular dynamics analysis of lipid nanobubbles.

Features
========

Graphical User Interface Modules
--------------------------------

LNB-MDT provides three main functional modules:

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin: 20px 0;">

   <div style="background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%); color: white; padding: 10px; border-radius: 8px; text-align: center;">
   <h4 style="margin-center: 0;">Preparation Module</h4>
   </div>

   <div style="background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%); color: white; padding: 10px; border-radius: 8px; text-align: center;">
   <h4 style="margin-center: 0;">Analysis Module</h4>
   </div>

   <div style="background: linear-gradient(135deg, #43e97b 0%, #38f9d7 100%); color: white; padding: 10px; border-radius: 8px; text-align: center;">
   <h4 style="margin-center: 0;">Visualization Module</h4>
   </div>

   </div>

Molecular Dynamics Analysis
----------------------------

LNB-MDT provides various analysis types:

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 15px; margin: 20px 0;">

   <div style="background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%); color: white; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0;">Area Per Lipid</h4>
   <p style="margin-bottom: 0;">Lipid-level; LNB-level</p>
   </div>

   <div style="background: linear-gradient(135deg, #43e97b 0%, #38f9d7 100%); color: white; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0;">Anisotropy</h4>
   <p style="margin-bottom: 0;">LNB-level</p>
   </div>

   <div style="background: linear-gradient(135deg, #fa709a 0%, #fee140 100%); color: white; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0;">Cluster Analysis</h4>
   <p style="margin-bottom: 0;">LNB-level</p>
   </div>

   <div style="background: linear-gradient(135deg, #a8edea 0%, #fed6e3 100%); color: #333; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0;">Order Parameter</h4>
   <p style="margin-bottom: 0;">Lipid-level; LNB-level</p>
   </div>

   <div style="background: linear-gradient(135deg, #ffeaa7 0%, #fab1a0 100%); color: #333; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0;">Gyration</h4>
   <p style="margin-bottom: 0;">LNB-level</p>
   </div>

   <div style="background: linear-gradient(135deg, #a29bfe 0%, #6c5ce7 100%); color: white; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0;">Density</h4>
   <p style="margin-bottom: 0;">LNB-level</p>
   </div>

   </div>

Graphical User Interface
------------------------

* Modern Qt6 interface design
* Intuitive data visualization
* Drag-and-drop file operations
* VMD integration support

High-Performance Computing
--------------------------

* Parallel processing support
* Configurable computation parameters
* Optimized algorithm implementation

Simplified Command Line Interface
---------------------------------

* **Short Parameter Aliases**: All parameters have short aliases (e.g., ``-g`` for ``--gro-file``)
* **Simplified Input Formats**: Support intuitive formats like ``DPPC:PO4,CHOL:ROH``
* **Backward Compatibility**: Traditional formats are still fully supported

Contents
========

.. toctree::
   :maxdepth: 2

   installation
   quickstart
   analysis_modules

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`