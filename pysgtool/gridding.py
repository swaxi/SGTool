# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2024-2026 Mark Jessell
"""Gridding scattered points, and normalising several grids together."""

import os

import numpy as np

from .calcs import raster_io, sgt_metadata
from .calcs.mba_suggest import suggest_mba_parameters
from .calcs.saga_mba_gridding import mba_gridding
from .calcs.filter_table import _plugin_version
from .grid import Grid


def grid_points(x, y, z, cell_size=None, levels=None, epsilon=1e-4, ignore_below=None, crs=None,
                name="points_bspline", field="z", source=None, callback=None):
    """Grid scattered points with a multilevel B-spline (after SAGA), on the same cells as
    the plugin. cell_size and levels: chosen from the spacing of the points if not given.
    epsilon: threshold error. ignore_below: leave out smaller values."""
    x, y, z = (np.asarray(v, dtype=float).ravel() for v in (x, y, z))
    if not x.size == y.size == z.size:
        raise ValueError("x, y and z must have the same length")
    keep = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    if ignore_below is not None:
        keep &= z >= ignore_below
    x, y, z = x[keep], y[keep], z[keep]
    if x.size < 3:
        raise ValueError("Need at least 3 points with a value")
    wkt = raster_io.crs_to_wkt(crs)
    geographic = raster_io.is_geographic(wkt)
    auto = cell_size is None
    if auto or levels is None:
        chosen = suggest_mba_parameters(x, y, geographic=geographic, cell_size=cell_size)
        cell_size = chosen["cell_size"]
        levels = chosen["max_levels"] if levels is None else levels
    if cell_size <= 0:
        raise ValueError("cell_size must be greater than 0")
    # the lower-left node is snapped down to a multiple of the cell size, as in the plugin
    xmin = float(np.floor(x.min() / cell_size) * cell_size)
    ymin = float(np.floor(y.min() / cell_size) * cell_size)
    nx = 2 + int((x.max() - xmin) / cell_size)
    ny = 2 + int((y.max() - ymin) / cell_size)
    values = mba_gridding(x, y, z, cellsize=cell_size, xmin=xmin, ymin=ymin, nx=nx, ny=ny,
                          epsilon=epsilon, level_max=int(levels), refinement=False,
                          callback=callback)
    gt = (xmin - cell_size / 2, cell_size, 0.0, ymin + (ny - 0.5) * cell_size, 0.0, -cell_size)
    grid = Grid(np.flipud(values), gt=gt, crs=wkt, geographic=geographic, name=name)
    if source:
        grid.path = os.path.join(os.path.dirname(source), name + ".tif")
    record = {"data_field": field, "cell_size": cell_size, "epsilon": epsilon,
              "max_levels": int(levels), "points_used": int(x.size)}
    if auto:
        record["settings_chosen_automatically"] = True
    if ignore_below is not None:
        record["ignore_values_below"] = ignore_below
    grid._provenance = sgt_metadata._build(
        grid.path or f"{name}.tif", source, "Multilevel B-spline gridding", record,
        _plugin_version(), file_stats=False)
    return grid


def normalise(grids, order=1):
    """Remove the gradient (order 1 or 2) of each grid, centre it on 0 and scale it to the
    standard deviation of the first, as the plugin's Normalise Grids does for a folder."""
    from .calcs.GeophysicalProcessor import GeophysicalProcessor

    if order not in (1, 2):
        raise ValueError("order is 1 or 2")
    out, std = [], None
    for g in grids:
        # float32, as the plugin: the same numbers as its files
        data = g.values.astype(np.float32)
        mask = np.isnan(data)
        values, std = GeophysicalProcessor(g.dx, g.dy, 0).normalise_array(
            data, mask, order == 1, std)
        if values is None:
            raise ValueError(f"{g.name} has no variation to normalise")
        values = np.where(mask, np.nan, values)
        out.append(g._derive(values, "Normalise grids", {"polynomial_order": order}, "_norm"))
    return out
