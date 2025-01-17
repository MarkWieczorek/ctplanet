import datetime
import ctplanet
from packaging.version import Version

# Project information
year = datetime.date.today().year
project = 'ctplanet'
copyright = "{}".format(year)
author = 'the ctplanet developers'
v = Version(ctplanet.__version__)
version = f'v{v.major}.{v.minor}'

# General configuration
extensions = [
    'sphinx_copybutton',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx_toolbox.collapse',
    'myst_parser',
]
myst_enable_extensions = [
    'attrs_block',
    'colon_fence',
]

# Autosummary pages will be generated by sphinx-autogen instead of sphinx-build
autosummary_generate = True

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
source_suffix = ['.rst', '.md']

# Options for HTML output
html_theme = 'sphinx_book_theme'
html_logo = "_static/logo.png"
html_theme_options = {
    "repository_branch": "master",
    "repository_url": "https://github.com/MarkWieczorek/ctplanet",
    "use_repository_button": True,
    "use_sidenotes": False,
    "use_edit_page_button": False,
    "use_source_button": False,
    "use_issues_button": False,
    "use_download_button": False,
    "logo": {
        "text": f'<span class="project-version">{version}</span>'
    },
    "home_page_in_toc": False,
}
html_static_path = ['_static']
html_show_sourcelink = False
html_show_sphinx = False
html_show_copyright = True
html_title = "ctplanet"
html_short_title = "ctplanet"
html_last_updated_fmt = "%b %d, %Y"
html_css_files = ["custom.css"]
html_context = {}
