# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'AIGCT'
copyright = '2025, Claudio Fratarcangeli, Ian Lee'
author = 'Claudio Fratarcangeli, Ian Lee'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["sphinx.ext.todo", "sphinx.ext.viewcode", "sphinx.ext.autodoc",
              # 'sphinx.ext.autosummary',  # Create neat summary tables
              'autoapi.extension'
]
# autosummary_generate = True  # Turn on sphinx.ext.autosummary

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']

# Changes made by CF

import sphinx_rtd_theme

html_theme = 'sphinx_rtd_theme'

autodoc_mock_imports = ['bs4', 'requests']

autoapi_dirs = ['../aigct']
autoapi_template_dir = "_templates/autoapi"

import os
import sys

sys.path.insert(0, os.path.abspath('..'))
