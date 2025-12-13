from setuptools import setup, find_packages

# Read requirements from requirements.txt
try:
    with open("requirements.txt", "r", encoding="utf-8") as f:
        requirements = [line.strip() for line in f if line.strip() and not line.startswith("#")]
except FileNotFoundError:
    requirements = []

# Read README for long description
try:
    with open("README.md", "r", encoding="utf-8") as f:
        long_description = f.read()
except FileNotFoundError:
    long_description = ""

setup(
    name="lnb-mdt",
    version="1.0.1",
    description="Lipid Nanobubble Molecular Dynamics Toolkit",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="RenXinYU",
    author_email="liqianyanqin01@163.com",  # 请根据实际情况修改
    url="https://github.com/xinyuren-bio/LNB-MDT",
    project_urls={
        "Bug Reports": "https://github.com/xinyuren-bio/LNB-MDT/issues",
        "Source": "https://github.com/xinyuren-bio/LNB-MDT",
        "Documentation": "https://lnb-mdt.readthedocs.io/",
    },
    packages=find_packages(),
    install_requires=requirements,
    python_requires=">=3.7",
    
    # Entry point for command-line interface
    entry_points={
        'console_scripts': [
            'LNB-MDT=LNB_MDT.main:main',
        ],
    },
    
    # Include package data files
    package_data={
        'LNB_MDT': [
            '*.ico',
            '*.tcl',
            '*.ini',
            '*.ui',
            '*.qrc',
            '*.mdp',
            '*.itp',
            '*.png',
            '*.qss',
            'preparation/files/*.mdp',
            'preparation/files/*.itp',
            'preparation/files/*.sh',
            'themes/*.qss',
            'images/**/*',
            'cases_lnb/*.gro',
            'cases_lnb/*.xtc',
        ],
    },
    
    # Ensure non-Python files are included
    include_package_data=True,
    
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],
    keywords="molecular-dynamics, lipid, nanobubble, simulation, analysis, mda-analysis",
)
