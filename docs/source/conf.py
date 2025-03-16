# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'template-intergrated'
copyright = '2025, brucekomike'
author = 'brucekomike'
release = 'v0.2.0'
language='en'
# 'en' 'jp' 'zh_CN'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["sphinx_inline_tabs",
"sphinxext.opengraph",
'sphinx_copybutton',
'myst_parser',
'sphinx.ext.githubpages',
]
myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "html_admonition",
    "html_image",
    "linkify",
    "replacements",
    "smartquotes",
    "strikethrough",
    "substitution",
    "tasklist",
]
source_suffix = {
   '.rst': 'restructuredtext',
   '.txt': 'markdown',
   '.md': 'markdown',
   '': 'markdown',
}

templates_path = ['_templates']
exclude_patterns = ['_build', 'build', 'Thumbs.db',  
  '.gitignore', '.gitattributes', '.git', 'Makefile', 
  '*.py', '*.bat', '*.sh', 'LICENSE'
  'requirements.txt', '*venv', 
  '.DS_Store', '._*',
]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_static_path = ['_static']
html_title = "furo template"
html_last_updated_fmt = ''
html_css_files = [
        "footer.css",
]
html_theme_options = {
    "source_repository": "https://github.com/brucekomike/furo-embeded",
    "source_branch": "main",
    "source_directory": "docs/source/",
    "navigation_with_keys": True,
}
highlight_language = 'text'
#html_logo = '_static/logo.svg'
ogp_site_url = 'https://brucekomike.github.io/furo-embeded/'
#ogp_image = '_static/logo.svg'
#ogp_image_alt = 'site logo'
ogp_site_name = 'furo template'
ogp_use_first_image = True
