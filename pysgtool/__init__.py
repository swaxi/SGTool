"""
pysgtool: the SGTool potential field and radiometric processing, without QGIS or ArcGIS.

    import pysgtool as sg

    mag = sg.read("mag.tif")                       # or sg.Grid(array, dx, dy, crs="EPSG:28350")
    rtp = mag.reduction_to_pole(date="2015-06-30") # or inclination=, declination=
    rtp.derivative(direction="z", power=1).save()  # mag_RTP_d1z.tif, beside mag.tif

The calculations are those of the SGTool plugin (the same code), and the files saved
carry the same processing history. sg.filters() lists the filters.
"""

from .calcs.raster_io import backend as io_backend, set_backend as set_io_backend
from .calcs.sgt_cancel import OperationCancelled
from .calcs.filter_table import FilterError, _plugin_version
from .grid import Grid, Stack
from ._filters import filters
from .io import Points, read, read_points, read_stack
from .gridding import grid_points, normalise
from .analysis import (
    efimov, history, igrf, kd, multifractal_spectrum, plot_wtmm, replay, ternary_gradient, ud, wtmm,
)
from .results import EulerSolutions, Features, Worms
from . import petrophysics

__version__ = _plugin_version()

__all__ = [
    "Grid", "Stack", "Points", "EulerSolutions", "Features", "Worms",
    "read", "read_stack", "read_points", "grid_points", "normalise", "filters",
    "efimov", "kd", "ud", "ternary_gradient", "igrf",
    "wtmm", "plot_wtmm", "multifractal_spectrum", "history", "replay", "petrophysics",
    "FilterError", "OperationCancelled", "io_backend", "set_io_backend",
]
