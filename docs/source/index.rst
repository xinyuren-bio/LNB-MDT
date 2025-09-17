LNB-MDT: Lipid NanoBubble Molecular Dynamics Toolbox
====================================================

**LNB-MDT** (Lipid NanoBubble Molecular Dynamics Toolbox) is a comprehensive toolbox designed for molecular dynamics analysis of lipid nanobubbles.

Features
========

Molecular Dynamics Analysis
----------------------------

* **PCA Analysis**: Principal Component Analysis for analyzing lipid molecular conformational changes
* **Area Analysis**: Voronoi tessellation area calculation  
* **Curvature Analysis**: Mean and Gaussian curvature calculation
* **Height Analysis**: Lipid molecular height distribution analysis
* **Cluster Analysis**: Lipid molecular aggregation behavior analysis
* **Anisotropy Analysis**: Molecular orientation anisotropy calculation
* **Gyration Analysis**: Molecular radius of gyration calculation
* **Sz Order Parameter Analysis**: Lipid chain order parameter analysis
* **N-Cluster Analysis**: Cluster count statistics
* **Radial Distribution Analysis**: Radial distribution function calculation

Machine Learning Integration
----------------------------

* **Automatic Parameter Optimization**: Bayesian optimization for finding optimal analysis parameters
* **Anomaly Detection**: Identify unusual patterns in molecular dynamics trajectories
* **Property Prediction**: Predict molecular properties using machine learning models
* **Feature Engineering**: Advanced feature extraction from trajectory data
* **Model Evaluation**: Comprehensive model performance assessment and visualization

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
   parameter_input_guide
   analysis_modules
   machine_learning
   api_reference

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`