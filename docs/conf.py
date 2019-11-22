# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('..'))
import sphinx_bootstrap_theme
import sphinx_gallery
from sphinx_gallery.sorting import FileNameSortKey
# -- Project information -----------------------------------------------------
# sys.path.append(os.path.pardir)


project = 'CANA'
copyright = '2019, M. De Pra'
author = 'M. De Pra'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
#    'sphinx.ext.autodoc',
#    'sphinx.ext.autosummary',
#    'sphinx.ext.coverage',
#    'sphinx.ext.mathjax',
#    'sphinx.ext.doctest',
#    'sphinx.ext.viewcode',
#    'sphinx.ext.extlinks',
#    'matplotlib.sphinxext.plot_directive',
    'sphinx_gallery.gen_gallery']

plot_formats = [("png", 90)]
plot_html_show_formats = False
plot_html_show_source_link = False
#min_reported_time = 0
#if 'SOURCE_DATE_EPOCH' in os.environ:
#    min_reported_time = sys.maxint if sys.version_info[0] == 2 else sys.maxsize

#image_scrapers = ('matplotlib', 'cana')

sphinx_gallery_conf = {
    'doc_module': ('sphinx_gallery', 'numpy'),
    'examples_dirs': '../cookbook',
    'gallery_dirs': 'gallery',
    'backreferences_dir': False,
#    'image_scrapers': image_scrapers,
#    # specify the order of examples to be according to filename
#    'within_subsection_order': FileNameSortKey,
#    'min_reported_time': min_reported_time,
#    'show_memory': True,
#    # capture raw HTML or, if not present, __repr__ of last expression in
#    # each code block
#    'capture_repr': ('_repr_html_', '__repr__'),
}


# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The master toctree document.
master_doc = 'index'

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'bootstrap'
html_theme_path = sphinx_bootstrap_theme.get_html_theme_path()
html_favicon = u'_static/Cana_favicon.png'
#html_logo = '_static/Cana_favicon.png'

html_theme_options = {
    'bootswatch_theme': "sandstone",
    'navbar_title': 'CANA',
    'navbar_site_name': "Site",
    'navbar_links': [("About", "about"),
                     ("Install", "install"),
                     ("Cookbook", "gallery/index"),
                     ("Docs", "docs"),
                     ("Contribute", "contribute"),
                     ("Cite", "cite")],
    # Render the next and previous page links in navbar. (Default: true)
    'navbar_sidebarrel': False,
    # Render the current pages TOC in the navbar. (Default: true)
    'navbar_pagenav': False,
    # Tab name for the current pages TOC. (Default: "Page")
    'navbar_pagenav_name': "This page",
    # Global TOC depth for "site" navbar tab. (Default: 1)
    # Switching to -1 shows all levels.
    'globaltoc_depth': 1,
    # Include hidden TOCs in Site navbar?
    # Note: If this is "false", you cannot have mixed ``:hidden:`` and
    # non-hidden ``toctree`` directives in the same page, or else the build
    # will break.
    # Values: "true" (default) or "false"
    'globaltoc_includehidden': "false",
    # HTML navbar class (Default: "navbar") to attach to <div> element.
    # For black navbar, do "navbar navbar-inverse"
    'navbar_class': "navbar navbar-default",
    # Fix navigation bar to top of page?
    # Values: "true" (default) or "false"
    'navbar_fixed_top': "false",
    # Location of link to source.
    # Options are "nav" (default), "footer" or anything else to exclude.
    'source_link_position': "footer",
    'bootstrap_version': "3",

}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_css_files = [
    '"style.css',
]
def setup(app):
	app.add_stylesheet("style.css")
