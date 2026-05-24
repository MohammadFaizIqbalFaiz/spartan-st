import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

project = "Spartan"
author = "Mohammad Faiz Iqbal Faiz"
release = "0.2.0"
version = "0.2.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "myst_parser",
    "nbsphinx",
]

autosummary_generate = True
autodoc_typehints = "description"
napoleon_google_docstring = True
napoleon_numpy_docstring = True

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_logo = "_static/Spartan_logo.png"

html_context = {
    "display_github": True,
    "github_user": "MohammadFaizIqbalFaiz",
    "github_repo": "spartan-st",
    "github_version": "main",
    "conf_py_path": "/docs/",
}

html_css_files = [
    "custom.css",
]


myst_enable_extensions = ["dollarmath", "amsmath"]
nbsphinx_execute = "never"
autodoc_mock_imports = [
    "scanpy",
    "squidpy",
    "spatialdata",
    "spatialdata_plot",
    "anndata",
    "igraph",
    "leidenalg",
    "sklearn",
    "scipy",
    "statsmodels",
    "joblib",
    "matplotlib",
    "pandas",
    "numpy",
]
