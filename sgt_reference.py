"""
Links to sections of the SGTool reference (docs/reference.md) on GitHub, for the
Ctrl+click help of the QGIS dialog (sgt_help.py) and the ArcGIS Pro GUI.

No QGIS or ArcGIS dependency. The docs are not in the plugin zip (export-ignore),
so the guide is read online, at the heading's anchor.
"""

import re

REFERENCE_URL = "https://github.com/swaxi/SGTool/blob/main/docs/reference.md"

HINT = "Ctrl+click: open this in the SGTool reference"


def anchor(heading):
    """GitHub's anchor for a heading: lower case, punctuation (including symbols
    such as ●) dropped, spaces as hyphens."""
    text = re.sub(r"[^\w\- ]", "", heading.strip().lower())
    return text.replace(" ", "-")


def section_url(heading):
    """The URL of a heading of docs/reference.md on GitHub."""
    return REFERENCE_URL + "#" + anchor(heading)
