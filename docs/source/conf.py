# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
# current_file_path = os.path.abspath(__file__)
# current_directory = os.path.dirname(current_file_path)
# target_path = os.path.abspath(os.path.join(current_directory, '../../src'))
# sys.path.insert(0, target_path)

current_file_path = os.path.abspath(__file__)
target_directory = os.path.abspath(os.path.join(os.path.dirname(current_file_path), '..', '..'))
sys.path.insert(0, target_directory)

project = 'PMD: Planar Multi-Body Dynamics Open Source Software'
copyright = '2025, Giacomo Cangi'
author = 'Giacomo Cangi'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon', 'myst_parser', 'sphinx_design']
myst_enable_extensions = ["colon_fence"]
templates_path = ['_templates']
exclude_patterns = []
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

napoleon_google_docstring = True                    # use google style docstrings
napoleon_numpy_docstring = False                    # do not use numpy style docstrings
napoleon_include_init_with_doc = True               # include the docstring of the __init__ method
napoleon_include_private_with_doc = False           # include the docstring of private methods
napoleon_include_special_with_doc = False           # include the docstring of special methods
napoleon_use_admonition_for_examples = False        # Use `Example` instead of `.. code-block:: python`
napoleon_use_rtype = False                          # include the return type in the docstring

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'sphinx_rtd_theme'
html_theme = 'sphinx_book_theme'
# html_theme = 'pydata_sphinx_theme'
html_static_path = ['_static']
html_css_files = ['custom.css']