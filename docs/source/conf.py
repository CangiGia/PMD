# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
current_file_path = os.path.abspath(__file__)
current_directory = os.path.dirname(current_file_path)
target_path = os.path.abspath(os.path.join(current_directory, '../../src'))
sys.path.insert(0, target_path)

project = 'PMD: Planar Multi-Body Dynamics Open Source Simulation Software'
copyright = '2025, Giacomo Cangi'
author = 'Giacomo Cangi'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon', 'myst_parser']
templates_path = ['_templates']
exclude_patterns = []
language = 'sphinx-quickstart'
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'classic'
html_static_path = ['_static']