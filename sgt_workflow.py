# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2024-2026 Mark Jessell
"""
The replay of a processing history is calcs/replay.py, so that it is part of the pysgtool package too.
This module keeps the name the QGIS plugin and the ArcGIS Pro tools import it by.
"""

try:  # as part of the QGIS plugin package
    from .calcs import replay as _module
except ImportError:  # from the ArcGIS Pro GUI, where calcs is a top-level package
    from calcs import replay as _module

globals().update({k: v for k, v in vars(_module).items() if not k.startswith("__")})
