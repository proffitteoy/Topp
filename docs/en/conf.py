from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

project = "Topp"
author = "Topp contributors"
copyright = "2026, Topp contributors"
version = "1.0"
release = "1.0.1"
language = "en"

extensions = [
    "myst_parser",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx_copybutton",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}
master_doc = "index"
exclude_patterns = ["_build"]

myst_enable_extensions = [
    "colon_fence",
    "dollarmath",
    "fieldlist",
    "substitution",
]
myst_dmath_double_inline = True
myst_heading_anchors = 3

html_theme = "furo"
html_title = "Topp 1.0 documentation"
html_static_path: list[str] = []
html_theme_options = {
    "source_repository": "https://github.com/proffitteoy/Topp/",
    "source_branch": "main",
    "source_directory": "docs/en/",
    "light_css_variables": {
        "color-brand-primary": "#1f5f79",
        "color-brand-content": "#1f5f79",
    },
    "dark_css_variables": {
        "color-brand-primary": "#6cb6d4",
        "color-brand-content": "#6cb6d4",
    },
}

copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True
