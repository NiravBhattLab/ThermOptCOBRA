import os
import sys
sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------

project = 'ThermOptCobra'
copyright = '2024, Pavan'
author = 'Pavan and Nirav'

# The full version, including alpha/beta/rc tags
release = '0.1'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.napoleon',
              'sphinx.ext.viewcode',
              'sphinx_rtd_theme',
              'sphinx.ext.mathjax',
              'sphinxcontrib.matlab',
              'sphinx.ext.linkcode'
]

def linkcode_resolve(domain, info):
    filename = info['module'].replace('.', '/')
    return 'https://github.com/NiravBhattLab/ThermOptCOBRA/blob/master/'+filename+'/'+info['fullname']+'.m'

napoleon_google_docstring = True
napoleon_custom_sections = [('INPUTS','params_style'),('INPUT','params_style'),
                            ('OUTPUTS','params_style'),('OUTPUT','params_style'),
                            ('REQUIRED INPUTS','params_style'),('REQUIRED INPUT','params_style'),
                            'Authors',('OPTIONAL INPUTS','params_style'),
                            ('OPTIONAL INPUT','params_style'),('USAGE','params_style')]

# Matlab module related configurations
matlab_src_dir = os.path.abspath(os.path.join('..'))
primary_domain = 'mat'

# remove path in function names
add_module_names = False
# You can specify multiple suffix as a list of string:
#
source_suffix = ['.rst']

# The encoding of source files.
#
source_encoding = 'utf-8'

root_doc = 'index'
# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'
# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
html_favicon = 'favicon.png'

# -- Options for EPUB output
epub_show_urls = 'footnote'
