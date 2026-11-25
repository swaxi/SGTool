"""
Henkel plots and other scatter plots, heat maps and histograms (need matplotlib).

    d, s = petrophysics.pair(density, susceptibility)
    fig = petrophysics.heatmap(d.values, s.values, "Henkel plot", "Density (g/cc)",
                               "Susceptibility (SI)")
    petrophysics.save_figure(fig, "henkel", ["png"])

heatmap, scatter, histogram and the *_by_category versions are those of the plugin.
"""

import numpy as np

from .calcs import raster_io
from .calcs.Petrophysics import (  # noqa: F401
    category_label, group_distributions, heatmap, heatmap_by_category, histogram,
    histogram_by_category, index_categories, save_figure, scatter, to_number,
)


def _ravelled(plot):
    """The plots take columns of points: let them take grids too."""
    import functools

    @functools.wraps(plot)
    def call(*args, **kwargs):
        args = [a.ravel() if isinstance(a, np.ndarray) and a.ndim == 2 else a for a in args]
        return plot(*args, **kwargs)

    return call


heatmap, scatter, histogram = _ravelled(heatmap), _ravelled(scatter), _ravelled(histogram)
heatmap_by_category = _ravelled(heatmap_by_category)
histogram_by_category = _ravelled(histogram_by_category)


def pair(density, susceptibility):
    """Two Grids on the cells of the coarser one, each of its cells taking the nearest
    cell of the finer grid."""
    same = not (density.crs and susceptibility.crs) or raster_io.same_crs(
        density.crs, susceptibility.crs)
    area = susceptibility.dx * susceptibility.dy
    if not same:
        area = raster_io.cell_area(susceptibility.gt, susceptibility.shape,
                                   susceptibility.crs, density.crs)
    if area > density.dx * density.dy:
        return density.resample_to(susceptibility, "nearest"), susceptibility
    return density, susceptibility.resample_to(density, "nearest")


def categories(grid, path, field):
    """(index of the polygon category at each cell of a grid, category names) from a field
    of a polygon file, -1 where there is no polygon: for the *_by_category plots."""
    import shapely

    v = raster_io.read_vector(path)
    if field not in v["fields"]:
        raise ValueError(f"{path} has no field {field}: {', '.join(v['fields'])}")
    if v["wkt"] and grid.crs and not raster_io.same_crs(v["wkt"], grid.crs):
        raise ValueError("The polygons and the grid have different coordinate systems")
    rows = [(w, category_label(t)) for w, t in zip(v["wkb"], v["fields"][field]) if w is not None]
    index, names = index_categories([label for _w, label in rows])
    polygons = [(shapely.from_wkb(bytes(w)).wkt, int(i)) for (w, _label), i in zip(rows, index)]
    return raster_io.rasterize(polygons, grid.gt, grid.shape, grid.crs), names
