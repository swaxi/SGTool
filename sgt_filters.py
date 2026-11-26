# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2024-2026 Mark Jessell
"""
The shared filter table is calcs/filter_table.py, so that it is part of the pysgtool package too.
This module keeps the name the QGIS plugin and the ArcGIS Pro tools import it by.
"""

try:  # as part of the QGIS plugin package
    from .calcs import filter_table as _module
except ImportError:  # from the ArcGIS Pro GUI, where calcs is a top-level package
    from calcs import filter_table as _module

globals().update({k: v for k, v in vars(_module).items() if not k.startswith("__")})
