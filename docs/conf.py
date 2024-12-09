# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys

# sys.path.insert(0,os.path.abspath(".."))  # must be .., ".." dir is D:\OneDrive\Github\MyRepos\PostMD




project = 'postmd Documentation'
copyright = '2024, Shusong Zhang'
author = 'Shusong Zhang'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon", # Google style docstring
    "sphinx_multiversion",
    "sphinx_copybutton",   # add copy button to code block
    "myst_parser",         # markdown parser
    # "autodoc2",
    "autoapi.extension"
]

autoapi_dirs = ['../postmd']
def skip_submodules(app, what, name, obj, skip, options):
    if what == "module":
        skip = True
    return skip


def setup(sphinx):
    sphinx.connect("autoapi-skip-member", skip_submodules)

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

# myst_parser configuration
myst_enable_extensions = [
    "amsmath",
    "attrs_inline",
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "html_admonition",
    "html_image",
    # "linkify",
    "replacements",
    "smartquotes",
    "strikethrough",
    "substitution",
    "tasklist",
]


# autodoc2_module_all_regexes = [
#     r"postmd\..*",
# ]

templates_path = ['_templates']
html_sidebars = {
    '**': [
        'versioning.html',
    ],
}

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']


autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    # 'special-members': '__init__',
    'undoc-members': True,
    'exclude-members': '__weakref__',
    # 'private-members': True, # open the docstring of _xxx method
}


autosummary_generate = True
# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
add_module_names = True

# Both the class’ and the __init__ method’s docstring are concatenated and inserted.
autoclass_content = "both"

# This value controls the docstrings inheritance. If set to True the docstring for classes or methods, if not explicitly set, is inherited from parents.
autodoc_inherit_docstrings = True



# from sphinx.builders.html import StandaloneHTMLBuilder
# StandaloneHTMLBuilder.supported_image_types = [
#     'image/svg+xml',
#     'image/gif',
#     'image/png',
#     'image/jpeg'
# ]