Installation Guide
==================

Conda Installation
------------------

The easiest way to install LNB-MDT is through the Conda environment manager:

.. code-block:: bash

   conda config --add channels conda-forge
   conda create -n LNB-MDT -c conda-forge python=3.11
   conda activate LNB-MDT
   git clone https://github.com/xinyuren-bio/LNB-MDT.git
   cd LNB-MDT
   pip install -r requirements.txt

This will install LNB-MDT and all its dependencies in a new virtual environment.

If you don't have Conda installed yet, we recommend downloading and installing Miniconda, which is the lightweight version of Conda.

PyPI Installation
-----------------

You can also install LNB-MDT from the Python Package Index:

.. code-block:: bash

   pip install LNB-MDT

Or install the development version:

.. code-block:: bash

   pip install https://github.com/xinyuren-bio/LNB-MDT/archive/main.zip

Dependencies
------------

LNB-MDT uses MDAnalysis for all analysis calculations and NumPy for numerical computations.

As mentioned above, the easiest way to install these packages along with LNB-MDT is using Conda. However, you can also install MDAnalysis and NumPy using pip or from source.

System Requirements
-------------------

- **Windows**: Windows 10/11 (64-bit)
- **macOS**: macOS 10.15 (Catalina) or later
- **Linux**: Ubuntu 18.04+, CentOS 7+, or other mainstream distributions

Installation Complete
---------------------

Congratulations! You have successfully installed LNB-MDT. You can now start using this powerful molecular dynamics analysis toolbox.

Next step: Check out :doc:`quickstart` to learn the basics.
