Installation
============

Conda
-----

The easiest way to install **LNB-MDT** is through the `conda-forge
<https://anaconda.org/conda-forge>`__ channel of `Conda
<https://docs.conda.io/en/latest/index.html>`__::

    conda config --add channels conda-forge
    conda create -n LNB-MDT -c conda-forge python=3.11
    conda activate LNB-MDT
    git clone https://github.com/xinyuren-bio/LNB-MDT.git
    cd LNB-MDT
    pip install -r requirements.txt

This will install **LNB-MDT** along with all of its dependencies into a new virtual environment.

If you do not already have Conda installed on your machine, we recommend
downloading and installing `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`__
--- a lightweight version of Conda.

PyPI
----

It's also possible to install **LNB-MDT** from the `Python Package
Index <https://pypi.org/>`__. You can do this using `pip`::

    pip install LNB-MDT

Alternatively, you can also install the in-development version with::

    pip install https://github.com/xinyuren-bio/LNB-MDT/archive/main.zip

Dependencies
------------

**LNB-MDT** uses `MDAnalysis <https://www.mdanalysis.org/>`__ to carry out all analysis
calculations, and `NumPy <https://numpy.org/>`__ for numerical computations.

As mentioned above, the simplest way to install these packages,
along with **LNB-MDT**, is with `Conda <https://docs.conda.io/en/latest/index.html>`__.
However, it is also possible to install MDAnalysis and NumPy using pip, or from source. See
the `MDAnalysis <https://userguide.mdanalysis.org/stable/installation.html>`_ and
`NumPy <https://numpy.org/install/>`_
installation instructions for further information.

System Requirements
-------------------

Operating System Support
~~~~~~~~~~~~~~~~~~~~~~~~

**LNB-MDT** supports the following operating systems:

- **Windows**: Windows 10/11 (64-bit)
- **macOS**: macOS 10.15 (Catalina) or later
- **Linux**: Ubuntu 18.04+, CentOS 7+, or other mainstream distributions

Hardware Requirements
~~~~~~~~~~~~~~~~~~~~~

**Minimum Requirements:**
- CPU: Dual-core processor
- Memory: 8GB RAM
- Storage: 2GB available space

**Recommended Requirements:**
- CPU: Quad-core or more processors
- Memory: 16GB+ RAM
- Storage: 5GB+ available space
- GPU: CUDA-compatible graphics card (optional, for acceleration)

Installation Methods
--------------------

Method 1: Using Installation Scripts (Recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is the simplest installation method. The installation scripts will automatically handle all dependencies.

Windows Installation
^^^^^^^^^^^^^^^^^^^^

1. **Download and install Git** (if not already installed)
   - Visit https://git-scm.com/download/win
   - Download and run the installer

2. **Clone repository and run installation script**
   
   .. code-block:: cmd

      git clone https://github.com/xinyuren-bio/LNB-MDT.git
      cd LNB-MDT
      install.bat

3. **Wait for installation to complete**
   - The script will automatically create a conda environment
   - Install all required Python packages
   - Verify installation success

macOS/Linux Installation
^^^^^^^^^^^^^^^^^^^^^^^^^^

1. **Ensure Git is installed**
   
   .. code-block:: bash

      # macOS (using Homebrew)
      brew install git
      
      # Ubuntu/Debian
      sudo apt update && sudo apt install git

2. **Clone repository and run installation script**
   
   .. code-block:: bash

      git clone https://github.com/xinyuren-bio/LNB-MDT.git
      cd LNB-MDT
      chmod +x install.sh
      ./install.sh

3. **Wait for installation to complete**

Method 2: Manual Installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you need more control or encounter issues with the installation scripts, you can install manually.

Step 1: Install Conda
^^^^^^^^^^^^^^^^^^^^^

Choose a Conda distribution:

- **Miniconda**: Lightweight, contains only conda and Python
- **Anaconda**: Full distribution with many scientific packages

Download links:
- Miniconda: https://docs.conda.io/en/latest/miniconda.html
- Anaconda: https://www.anaconda.com/products/distribution

Step 2: Create Virtual Environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   # Create new conda environment
   conda create -n LNB-MDT python=3.11 -y
   
   # Activate environment
   conda activate LNB-MDT

Step 3: Clone Repository
^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   git clone https://github.com/xinyuren-bio/LNB-MDT.git
   cd LNB-MDT

Step 4: Install Dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   # Install basic dependencies
   pip install -r requirements.txt
   
   # Install machine learning dependencies (optional)
   pip install scikit-learn scipy matplotlib seaborn joblib

Step 5: Verify Installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   # Check Python version
   python --version
   
   # Check key dependencies
   python -c "import MDAnalysis, numpy, pandas, PySide6; print('All dependencies installed successfully!')"
   
   # Test main program
   python main.py --version

Getting Help
------------

If you encounter issues during installation:

1. **Check log files**: Review logs generated by installation scripts
2. **Check system requirements**: Ensure all system requirements are met
3. **Search known issues**: Check GitHub Issues
4. **Contact support**: Send email to zy2310205@buaa.edu.cn

.. raw:: html

   <div style="background-color: #e3f2fd; padding: 20px; border-radius: 8px; margin: 20px 0; text-align: center;">
   <h3 style="color: #1976d2; margin-top: 0;">🎉 Installation Complete!</h3>
   <p>Congratulations! You have successfully installed LNB-MDT. You can now start using this powerful molecular dynamics analysis toolbox.</p>
   <p><strong>Next step:</strong> Check out the <a href="quickstart.html">Quick Start Guide</a> to learn basic usage.</p>
   </div>
