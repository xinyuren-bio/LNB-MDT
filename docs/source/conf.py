# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'LNB-MDT'
copyright = '2025, renxinyu'
author = 'renxinyu'
release = 'v.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',  # 核心：从代码 docstrings 自动生成文档
    'sphinx.ext.napoleon', # 支持 Google 和 NumPy 风格的 docstrings
    'myst_parser'          # 允许你使用 Markdown (.md) 文件代替 .rst
]

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

import os
import sys

sys.path.insert(0, os.path.abspath('../../'))

