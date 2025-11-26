@echo off
chcp 65001 >nul
setlocal enabledelayedexpansion

echo 🚀 LNB-MDT v1.0 Installation Script
echo ================================

REM Check if conda is installed
where conda >nul 2>&1
if %errorlevel% neq 0 (
    echo ❌ Error: conda not found. Please install Miniconda or Anaconda first.
    echo Download: https://docs.conda.io/en/latest/miniconda.html
    pause
    exit /b 1
)

echo ✅ Detected conda

REM Check Python version
for /f "tokens=2" %%i in ('python --version 2^>^&1') do set python_version=%%i
echo Current Python version: %python_version%

REM Create conda environment
echo 📦 Creating conda environment...
conda env list | findstr "LNB-MDT" >nul
if %errorlevel% equ 0 (
    echo ⚠️  Environment LNB-MDT already exists, do you want to recreate it? (Y/N)
    set /p response=
    if /i "!response!"=="Y" (
        conda env remove -n LNB-MDT -y
        conda create -n LNB-MDT python=3.11 -y
    ) else (
        echo Using existing environment
    )
) else (
    conda create -n LNB-MDT python=3.11 -y
)

REM Activate environment and install dependencies
echo 🔧 Activating environment and installing dependencies...
call conda activate LNB-MDT

echo 📥 Installing Python dependencies...
pip install -r requirements.txt

REM Verify installation
echo 🔍 Verifying installation...
python -c "import sys; import importlib; required_packages = ['MDAnalysis', 'numpy', 'pandas', 'PySide6', 'scipy', 'matplotlib', 'seaborn']; print('Python version:', sys.version); print('Checking dependencies...'); [print(f'✅ {package}: Installed') if importlib.util.find_spec(package) else print(f'❌ {package}: Not installed') for package in required_packages]; print('\n🎉 Installation complete!'); print('\nUsage:'); print('1. Activate environment: conda activate LNB-MDT'); print('2. Start program: python main.py'); print('3. Command line analysis: python analysis/pca.py --help')"

echo.
echo 🎯 Installation complete!
echo ================================
echo 📋 Next steps:
echo 1. Activate environment: conda activate LNB-MDT
echo 2. Start GUI: python main.py
echo 3. View command line help: python analysis/pca.py --help
echo 4. View detailed documentation: type README.md
echo.
echo 📚 For more information, see README.md file
pause
