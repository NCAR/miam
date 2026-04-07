# Configuration file for the Sphinx documentation builder.
import os
import re

conf_dir = os.path.dirname(os.path.abspath(__file__))

# -- Project information -----------------------------------------------------

project = 'MIAM'
copyright = '2026, University Corporation for Atmospheric Research'
author = 'NCAR'

regex = r'project\(.*VERSION\s+(\d+\.\d+\.\d+)'
version = '0.0.0'
# Read the version from the cmake files (use absolute path from conf.py location)
cmake_file = os.path.join(conf_dir, '..', '..', 'CMakeLists.txt')
with open(cmake_file, 'r') as f:
    content = f.read()
    match = re.search(regex, content)
    if match:
        version = match.group(1)
release = version

# -- General configuration ---------------------------------------------------

extensions = [
    'breathe',
    'sphinx.ext.mathjax',
    'sphinx_copybutton',
    'sphinx_design',
]

# -- Breathe configuration ---------------------------------------------------

this_dir = os.path.dirname(os.path.abspath(__file__))
breathe_projects_dir = os.path.abspath(
    os.path.join(this_dir, "..", "..", "build", "docs", "doxygen", "xml")
)
breathe_projects = {"miam": breathe_projects_dir}
breathe_default_project = 'miam'

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_book_theme'
html_theme_options = {
    'repository_url': 'https://github.com/NCAR/miam',
    'use_repository_button': True,
    'use_issues_button': True,
    'use_edit_page_button': False,
}

# -- General ---------------------------------------------------

exclude_patterns = ['_build']
templates_path = []
