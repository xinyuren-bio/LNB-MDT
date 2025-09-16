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

.. raw:: html

   <div style="background-color: #e3f2fd; padding: 20px; border-radius: 8px; margin: 20px 0; text-align: center;">
   <h3 style="color: #1976d2; margin-top: 0;">🎉 Installation Complete!</h3>
   <p>Congratulations! You have successfully installed LNB-MDT. You can now start using this powerful molecular dynamics analysis toolbox.</p>
   <p><strong>Next step:</strong> Check out the <a href="quickstart.html">Quick Start Guide</a> to learn basic usage.</p>
   </div>
