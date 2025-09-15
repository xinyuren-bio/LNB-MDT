# 🚀 LNB-MDT Quick Start Guide

This guide will help you get started with LNB-MDT in 5 minutes.

## 📋 Prerequisites

- ✅ Miniconda or Anaconda installed
- ✅ Internet connection (for downloading dependencies)

## ⚡ Quick Installation

### Method 1: Using Installation Scripts (Recommended)

**Linux/macOS:**
```bash
./install.sh
```

**Windows:**
```cmd
install.bat
```

### Method 2: Manual Installation

```bash
# 1. Create environment
conda create -n LNB-MDT python=3.11 -y

# 2. Activate environment
conda activate LNB-MDT

# 3. Install dependencies
pip install -r requirements.txt
```

## 🎯 Quick Usage

### 1. Launch Graphical Interface

```bash
conda activate LNB-MDT
python main.py
```

### 2. Use Command Line Tools

```bash
# Activate environment
conda activate LNB-MDT

# PCA analysis (using example data)
python analysis/pca.py \
  --gro-file cases/lnb.gro \
  --xtc-file cases/md.xtc \
  --output-csv results/pca_test.csv \
  --residues "{'DPPC': ['PO4']}" \
  --start-frame 0 \
  --stop-frame 10 \
  --parallel \
  --verbose
```

## 📊 Example Analysis

### Basic Analysis Workflow

1. **Prepare Data Files**
   - `*.gro` - Topology file
   - `*.xtc` - Trajectory file

2. **Select Analysis Type**
   - PCA Analysis: Molecular conformational changes
   - Area Analysis: Voronoi tessellation area
   - Curvature Analysis: Membrane curvature calculation
   - Cluster Analysis: Molecular aggregation behavior

3. **Configure Parameters**
   - Residue groups: Specify molecular types to analyze
   - Frame range: Select time range for analysis
   - Computation parameters: k-value, cutoff distance, etc.

4. **Run Analysis**
   - Serial processing: Suitable for small datasets
   - Parallel processing: Suitable for large datasets

### Common Analysis Commands

```bash
# Area analysis
python analysis/area.py \
  --gro-file your_system.gro \
  --xtc-file your_trajectory.xtc \
  --residues "{'DPPC': ['PO4'], 'CHOL': ['ROH']}" \
  --k-value 20 \
  --parallel

# Curvature analysis
python analysis/curvature.py \
  --gro-file your_system.gro \
  --xtc-file your_trajectory.xtc \
  --residues "{'DPPC': ['PO4']}" \
  --method mean \
  --parallel

# Cluster analysis
python analysis/cluster.py \
  --gro-file your_system.gro \
  --xtc-file your_trajectory.xtc \
  --residues "{'DPPC': ['PO4']}" \
  --cutoff 8.0 \
  --parallel

### 🤖 Machine Learning Features

```bash
# Test ML module functionality
python test_ml_module.py

# Run ML demos
python ml_demo.py

# Parameter optimization example
python -c "
from ml import AnalysisParameterOptimizer
optimizer = AnalysisParameterOptimizer('area')
print('Parameter optimizer ready!')
"

# Anomaly detection example
python -c "
from ml import MDAnomalyDetector
detector = MDAnomalyDetector(method='isolation_forest')
print('Anomaly detector ready!')
"
```
```

## 🔧 Parameter Description

### Common Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--gro-file` | GRO topology file path | `cases/lnb.gro` |
| `--xtc-file` | XTC trajectory file path | `cases/md.xtc` |
| `--residues` | Residue group dictionary | `"{'DPPC': ['PO4']}"` |
| `--parallel` | Enable parallel processing | No parameter |
| `--n-jobs` | Number of parallel jobs | `4` or `-1` |
| `--start-frame` | Starting frame | `0` |
| `--stop-frame` | Ending frame | `1000` |
| `--verbose` | Verbose output | No parameter |

### Specific Parameters

| Analysis Type | Specific Parameter | Description |
|---------------|-------------------|-------------|
| Area Analysis | `--k-value` | K-value for Voronoi tessellation |
| Area Analysis | `--max-normal-angle` | Maximum normal angle |
| Curvature Analysis | `--method` | Curvature type (mean/gaussian) |
| Height Analysis | `--k-value` | K-value for height calculation |
| Cluster Analysis | `--cutoff` | Cluster cutoff distance |
| Sz Analysis | `--chain` | Chain type (sn1/sn2/both) |

## 📁 Output Files

After analysis completion, the following files will be generated:

```
results/
├── pca_results.csv          # PCA analysis results
├── pca_results_frames.csv   # Frame information
├── area_results.csv         # Area analysis results
├── area_results_frames.csv  # Frame information
└── ...
```

### CSV File Format

```csv
# Created by LNB-MDT v1.0
# PCA Analysis
# TYPE:Bubble
# Parameters:{'DPPC': ['PO4']}
Frames,Values
0,0.787
1,0.801
2,0.800
...
```

## 🚨 Common Issues

### Q: Program fails to start
**A:** Check if conda environment is properly activated:
```bash
conda activate LNB-MDT
python main.py
```

### Q: Analysis is very slow
**A:** Use parallel processing:
```bash
python analysis/pca.py --parallel --n-jobs 4
```

### Q: Insufficient memory
**A:** Process large trajectories in segments:
```bash
python analysis/pca.py --start-frame 0 --stop-frame 1000
python analysis/pca.py --start-frame 1000 --stop-frame 2000
```

### Q: Parameter format error
**A:** Check residues parameter format:
```bash
# Correct format
--residues "{'DPPC': ['PO4'], 'CHOL': ['ROH']}"

# Incorrect format
--residues {'DPPC': ['PO4']}  # Missing quotes
```

## 📚 Next Steps

- 📖 View complete documentation: [README.md](README.md)
- 🖥️ Command line detailed guide: [analysis/README_COMMAND_LINE.md](analysis/README_COMMAND_LINE.md)
- 🧪 Run more examples: [cases/](cases/)

## 🆘 Get Help

- View help information: `python analysis/pca.py --help`
- Check version: `python main.py --version`
- Report issues: Create GitHub Issue

---

**🎉 Congratulations! You have successfully started using LNB-MDT!**
