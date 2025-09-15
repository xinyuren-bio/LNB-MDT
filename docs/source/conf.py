# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'LNB-MDT'
copyright = '2025, XinyuRen'
author = 'XinyuRen'
release = 'v1.0'
version = '1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',      # 从代码 docstrings 自动生成文档
    'sphinx.ext.napoleon',     # 支持 Google 和 NumPy 风格的 docstrings
    'sphinx.ext.viewcode',     # 添加源代码链接
    'sphinx.ext.githubpages',  # GitHub Pages 支持
    'sphinx.ext.intersphinx',  # 链接到其他项目的文档
    'sphinx.ext.todo',         # 支持 todo 指令
    'sphinx.ext.coverage',     # 文档覆盖率检查
    'myst_parser',             # 支持 Markdown 文件
    'sphinx_copybutton',       # 复制代码按钮
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# 支持中文搜索
html_search_language = 'zh'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# Read the Docs 主题配置
html_theme_options = {
    'analytics_id': '',  # 如果需要，添加 Google Analytics ID
    'logo_only': False,
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': True,
    'vcs_pageview_mode': '',
    'style_nav_header_background': '#2980B9',
    # 侧边栏配置
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': False,  # 隐藏隐藏的页面
    'titles_only': False
}

# 自定义CSS
html_css_files = [
    'custom.css',
]

# 自定义JavaScript
html_js_files = [
    'custom.js',
]

# 网站图标
html_favicon = '../../icon.ico'

# 网站标题
html_title = 'LNB-MDT Documentation'

# 网站短标题
html_short_title = 'LNB-MDT'

# -- Extension configuration -------------------------------------------------

# Napoleon 配置
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = False
napoleon_type_aliases = None
napoleon_attr_annotations = True

# MyST 配置
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_admonition",
    "html_image",
    "linkify",
    "replacements",
    "smartquotes",
    "substitution",
]

# Intersphinx 配置
intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'pandas': ('https://pandas.pydata.org/docs/', None),
    'matplotlib': ('https://matplotlib.org/stable/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
    'sklearn': ('https://scikit-learn.org/stable/', None),
}

# Todo 配置
todo_include_todos = True

# 自动文档配置
autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'special-members': '__init__',
    'undoc-members': True,
    'exclude-members': '__weakref__'
}

# 代码高亮主题
pygments_style = 'sphinx'

# 添加路径
import os
import sys
sys.path.insert(0, os.path.abspath('../../'))

# 自定义配置
def setup(app):
    app.add_css_file('custom.css')
    app.add_js_file('custom.js')

