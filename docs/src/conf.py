# Configuration file for the Sphinx documentation builder.

project = 'MIAM'
copyright = '2026, University Corporation for Atmospheric Research'
author = 'NCAR'

extensions = [
    'breathe',
    'sphinx.ext.mathjax',
    'sphinx_copybutton',
    'sphinx_design',
]

# Breathe configuration
breathe_default_project = 'miam'

# Theme
html_theme = 'sphinx_book_theme'
html_static_path = ['_static']
html_theme_options = {
    'repository_url': 'https://github.com/NCAR/miam',
    'use_repository_button': True,
    'use_issues_button': True,
    'use_edit_page_button': False,
}

# General
exclude_patterns = ['_build']
templates_path = []
