# LNB-MDT v1.0

![LNB-MDT Logo](LNB-MDT.jpg)

**LNB-MDT** (Lipid NanoBubble Molecular Dynamics Toolkit) is a comprehensive toolkit designed for molecular dynamics simulations of lipid nanobubbles.

## Installation Guide

### Method 1: Using Installation Scripts (Recommended)

**Linux/macOS:**
```bash
git clone https://github.com/xinyuren-bio/LNB-MDT.git
cd LNB-MDT
./install.sh
```

**Windows:**
```cmd
git clone https://github.com/xinyuren-bio/LNB-MDT.git
cd LNB-MDT
install.bat
```

### Method 2: Manual Installation

#### 1. Install Conda

If you haven't installed Conda yet, please visit the following links to download the version suitable for your operating system:

- **Miniconda**: [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html)
- **Anaconda**: [https://www.anaconda.com/products/distribution](https://www.anaconda.com/products/distribution)

#### 2. Clone the Project

```bash
git clone https://github.com/xinyuren-bio/LNB-MDT.git
cd LNB-MDT
```

#### 3. Create Virtual Environment

```bash
# Create new conda environment
conda create -n LNB-MDT python=3.11

# Activate environment
conda activate LNB-MDT

# Install dependencies
pip install -r requirements.txt
```

#### 4. Verify Installation

```bash
# Check Python version
python --version

# Check key dependencies
python -c "import MDAnalysis, numpy, pandas, PySide6; print('All dependencies installed successfully!')"
```

## Documentation

For detailed documentation, including installation guide, quick start, user guide, and command line tools, please visit:

📚 **[Read the Docs - LNB-MDT Documentation](https://lnb-mdt.readthedocs.io/en/latest/)**

## File Structure

```
LNB-MDT/
├── main.py                 # Main program entry
├── requirements.txt        # Python dependencies
├── analysis/              # Analysis modules
│   ├── area.py           # Area analysis
│   ├── height.py         # Height analysis
│   ├── cluster.py        # Cluster analysis
│   ├── anisotropy.py     # Anisotropy analysis
│   ├── gyration.py       # Gyration analysis
│   ├── sz.py             # Sz order parameter analysis
│   └── density.py        # Density analysis (time and radius)
├── preparation/            # Preparation module
├── cases/                 # Example data
    ├── lnb.gro           # Example topology file
    ├── md.xtc            # Example trajectory file
    └── csv/              # Results file directory
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

**LNB-MDT v1.0** - Making lipid nanobubble simulations simpler and more efficient!
