# LNB-MDT v1.0

![LNB-MDT Logo](LNB-MDT.jpg)

**LNB-MDT** (Lipid NanoBubble Molecular Dynamics Toolkit) is a comprehensive toolkit designed for molecular dynamics simulations of lipid nanobubbles.

<details>
<summary>Table of Contents</summary>

- [Features](#features)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Quick Start](#quick-start)
- [User Guide](#user-guide)
- [Command Line Tools](#command-line-tools)
- [File Structure](#file-structure)
- [FAQ](#faq)
- [Contributing](#contributing)
- [License](#license)

</details>

## Features

### Molecular Dynamics Analysis
- **PCA Analysis**: Principal Component Analysis for analyzing lipid molecular conformational changes
- **Area Analysis**: Voronoi tessellation area calculation
- **Curvature Analysis**: Mean and Gaussian curvature calculation
- **Height Analysis**: Lipid molecular height distribution analysis
- **Cluster Analysis**: Lipid molecular aggregation behavior analysis
- **Anisotropy Analysis**: Molecular orientation anisotropy calculation
- **Gyration Analysis**: Molecular radius of gyration calculation
- **Sz Order Parameter Analysis**: Lipid chain order parameter analysis
- **N-Cluster Analysis**: Cluster count statistics
- **Radial Distribution Analysis**: Radial distribution function calculation


### Graphical User Interface
- Modern Qt6 interface design
- Intuitive data visualization
- Drag-and-drop file operations
- VMD integration support

### High-Performance Computing
- Parallel processing support
- Configurable computation parameters
- Optimized algorithm implementation

## System Requirements

### Operating System
- **Windows**: Windows 10/11 (64-bit)
- **macOS**: macOS 10.15 or higher
- **Linux**: Ubuntu 18.04+ or other mainstream distributions

### Software Dependencies
- **Python**: 3.11 or higher
- **Conda**: Miniconda or Anaconda
- **VMD**: 1.9.4 or higher (optional, for visualization)

### Hardware Requirements
- **Memory**: Minimum 8GB RAM, recommended 16GB+
- **Storage**: At least 2GB available space
- **Processor**: Multi-core processor supporting parallel computation

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

## Quick Start

### Launch Graphical Interface

```bash
# Activate environment
conda activate LNB-MDT

# Start main program
python main.py
```

### Use Command Line Tools

```bash
# PCA analysis example
python analysis/pca.py --gro-file cases/lnb.gro --xtc-file cases/md.xtc --residues "{'DPPC': ['PO4']}" --parallel --verbose

# Area analysis example
python analysis/area.py --gro-file cases/lnb.gro --xtc-file cases/md.xtc --residues "{'DPPC': ['PO4']}" --k-value 20 --parallel --verbose

```

## User Guide

### Graphical Interface Usage

1. **Launch Program**: Run `python main.py`
2. **Load Data**: Load GRO and XTC files through the interface
3. **Select Analysis**: Choose analysis type from the left menu
4. **Configure Parameters**: Set analysis parameters
5. **Run Analysis**: Click the run button to start analysis
6. **View Results**: Check analysis results in the results panel

### Main Function Modules

#### Generation Module
- Lipid nanobubble structure generation
- Custom parameter configuration

#### Analysis Module
- Multiple physical property analysis
- Real-time progress display
- Result visualization

#### Visualization Module
- Multiple chart types
- Custom styles
- Interactive with VMD

### VMD Integration

1. **Load Trajectory**: Drag and drop CSV files to Visualization Module
2. **Start VMD**: Click "Start VMD" in the interface
3. **Visualization**: Use VMD for molecular visualization

## Command Line Tools

LNB-MDT provides a complete command-line interface supporting batch processing and automated analysis.

### Common Parameters

All analysis modules support the following common parameters:

```bash
--gro-file GRO_FILE      # GRO file path (topology file)
--xtc-file XTC_FILE      # XTC file path (trajectory file)
--output-csv OUTPUT_CSV  # Output CSV file path
--residues RESIDUES      # Residue group dictionary string
--parallel               # Enable parallel processing
--n-jobs N_JOBS          # Number of parallel jobs (-1 means use all CPU cores)
--start-frame START_FRAME # Starting frame for analysis (0-indexed)
--stop-frame STOP_FRAME  # Ending frame for analysis (exclusive)
--step-frame STEP_FRAME  # Frame step size
--verbose                # Enable verbose output
```

### Analysis Module Examples

#### PCA Analysis
```bash
python analysis/pca.py \
  --gro-file cases/lnb.gro \
  --xtc-file cases/md.xtc \
  --output-csv results/pca_results.csv \
  --residues "{'DPPC': ['PO4'], 'CHOL': ['ROH']}" \
  --parallel \
  --n-jobs 4 \
  --verbose
```

#### Area Analysis
```bash
python analysis/area.py \
  --gro-file cases/lnb.gro \
  --xtc-file cases/md.xtc \
  --output-csv results/area_results.csv \
  --residues "{'DPPC': ['PO4']}" \
  --k-value 20 \
  --max-normal-angle 140 \
  --parallel \
  --verbose
```

#### Curvature Analysis
```bash
python analysis/curvature.py \
  --gro-file cases/lnb.gro \
  --xtc-file cases/md.xtc \
  --output-csv results/curvature_results.csv \
  --residues "{'DPPC': ['PO4']}" \
  --k-value 20 \
  --method mean \
  --parallel \
  --verbose
```

For detailed usage instructions, please refer to [analysis/README_COMMAND_LINE.md](analysis/README_COMMAND_LINE.md)


## File Structure

```
LNB-MDT/
├── main.py                 # Main program entry
├── requirements.txt        # Python dependencies
├── analysis/              # Analysis modules
│   ├── pca.py            # PCA analysis
│   ├── area.py           # Area analysis
│   ├── curvature.py      # Curvature analysis
│   ├── height.py         # Height analysis
│   ├── cluster.py        # Cluster analysis
│   ├── anisotropy.py     # Anisotropy analysis
│   ├── gyration.py       # Gyration analysis
│   ├── sz.py             # Sz order parameter analysis
│   ├── n_cluster.py      # N-cluster analysis
│   ├── rad.py            # Radial distribution analysis
│   └── README_COMMAND_LINE.md  # Command line usage guide
├── generation/            # Generation module
├── cases/                 # Example data
    ├── lnb.gro           # Example topology file
    ├── md.xtc            # Example trajectory file
    └── csv/              # Results file directory
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

**LNB-MDT v1.0** - Making lipid nanobubble simulations simpler and more efficient!
