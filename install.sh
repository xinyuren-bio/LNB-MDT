#!/bin/bash

# LNB-MDT Installation Script
# For Linux and macOS

echo "🚀 LNB-MDT v1.0 Installation Script"
echo "================================"

# Check if conda is installed
if ! command -v conda &> /dev/null; then
    echo "❌ Error: conda not found. Please install Miniconda or Anaconda first."
    echo "Download: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

echo "✅ Detected conda: $(conda --version)"

# Check Python version
python_version=$(python --version 2>&1 | awk '{print $2}' | cut -d. -f1,2)
required_version="3.11"

if [ "$(printf '%s\n' "$required_version" "$python_version" | sort -V | head -n1)" != "$required_version" ]; then
    echo "⚠️  Warning: Current Python version is $python_version, recommended Python $required_version or higher"
fi

# Create conda environment
echo "📦 Creating conda environment..."
if conda env list | grep -q "LNB-MDT"; then
    echo "⚠️  Environment LNB-MDT already exists, do you want to recreate it? (y/N)"
    read -r response
    if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]; then
        conda env remove -n LNB-MDT -y
        conda create -n LNB-MDT python=3.11 -y
    else
        echo "Using existing environment"
    fi
else
    conda create -n LNB-MDT python=3.11 -y
fi

# Activate environment
echo "🔧 Activating environment and installing dependencies..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate LNB-MDT

# Install dependencies
echo "📥 Installing Python dependencies..."
pip install -r requirements.txt

# Verify installation
echo "🔍 Verifying installation..."
python -c "
import sys
import importlib

required_packages = ['MDAnalysis', 'numpy', 'pandas', 'PySide6', 'scipy', 'matplotlib']

print('Python version:', sys.version)
print('Checking dependencies...')

for package in required_packages:
    try:
        importlib.import_module(package)
        print(f'✅ {package}: Installed')
    except ImportError:
        print(f'❌ {package}: Not installed')

print('\\n🎉 Installation complete!')
print('\\nUsage:')
print('1. Activate environment: conda activate LNB-MDT')
print('2. Start program: python main.py')
print('3. Command line analysis: python analysis/pca.py --help')
"

echo ""
echo "🎯 Installation complete!"
echo "================================"
echo "📋 Next steps:"
echo "1. Activate environment: conda activate LNB-MDT"
echo "2. Start GUI: python main.py"
echo "3. View command line help: python analysis/pca.py --help"
echo "4. View detailed documentation: cat README.md"
echo ""
echo "📚 For more information, see README.md file"
