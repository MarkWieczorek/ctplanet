import datetime
import pycrust

# Project information
year = datetime.date.today().year
project = 'pycrust'
copyright = "2018-{}, The pyCrust Developers".format(year)
author = 'Mark A. Wieczorek'
version = pycrust.__version__.split(sep='-')[0]  # use version of last tag


# General configuration
extensions = ['sphinx_rtd_theme',
              'sphinx.ext.autodoc',
              'sphinx.ext.napoleon']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# Options for HTML output
html_theme = 'sphinx_rtd_theme'
html_theme_options = {"display_version": True}
html_static_path = ['_static']
