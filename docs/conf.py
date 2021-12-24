# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
from pathlib import Path
import subprocess
import json

currdir = Path(__file__).resolve().parent # PROJECT_SOURCE_DIR/docs
PROJECT_SOURCE_DIR = currdir.parent

def obtain_cpp_files():
    script_path = PROJECT_SOURCE_DIR.joinpath('UnionFindPy/cpp/build_utils/cpp_files.py')

    if not script_path.exists():
        print('The project directory structure is corrupted.')
        sys.exit(1)

    p = subprocess.Popen([script_path], shell=False, stdin=None,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
    out, _ = p.communicate()
    parsed = json.loads(out)

    file_list = []
    for item in parsed:
        file_list.append('../' + str(Path(item['name']).relative_to(PROJECT_SOURCE_DIR)))
    return file_list

CPP_FILES = obtain_cpp_files()

sys.path.insert(0, PROJECT_SOURCE_DIR)
sys.path.insert(0, os.path.abspath('_ext'))


# -- Project information -----------------------------------------------------

project = 'UnionFind'
copyright = '2021, Chae-Yeun Park'
author = 'Chae-Yeun Park'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'breathe',
    'exhale',
    'sphinx.ext.autodoc'
]

# Setup the breathe extension
breathe_projects = {
    "UnionFindPy": "./_doxygen/xml"
}
breathe_default_project = "UnionFindPy"

# Setup the exhale extension
exhale_args = {
    # These arguments are required
    "containmentFolder":     "./api",
    "rootFileName":          "library_root.rst",
    "doxygenStripFromPath":  "..",
    # Heavily encouraged optional argument (see docs)
    "rootFileTitle":         "C++ Library API",
    # Suggested optional arguments
    "createTreeView":        True,
    # TIP: if using the sphinx-bootstrap-theme, you need
    # "treeViewIsBootstrap": True,
    "exhaleExecutesDoxygen": True,
    "exhaleDoxygenStdin":    "INPUT = " + ' '.join(CPP_FILES)
}

# Tell sphinx what the primary language being documented is.
primary_domain = 'cpp'

# Tell sphinx what the pygments highlight language should be.
highlight_language = 'cpp'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'UnionFindPy/cpp']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_context = {
    "display_github": True, # Integrate GitHub
    "github_user": "chaeyeunpark", # Username
    "github_repo": "UnionFind", # Repo name
    "github_version": "main", # Version
    "conf_py_path": "/docs/", # Path in the checkout to the docs root
}
