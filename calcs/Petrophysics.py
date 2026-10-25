"""
Henkel plots (density against log susceptibility) for apparent petrophysics
grids and for measured point data, with no QGIS dependency.
"""

import math

import numpy as np
from osgeo import gdal, ogr, osr

# Sequential ramp (light -> dark blue) for the heat maps, and the fixed
# categorical order for category colours, each paired with a marker shape so
# identity never rests on colour alone. Categories past the eighth share "Other".
HEAT_RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
CATEGORY_COLOURS = [
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100",
    "#e87ba4", "#008300", "#6250d6", "#e34948",
]
CATEGORY_MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]
OTHER_COLOUR = "#8a8a85"
TEXT_SECONDARY = "#5f5e5a"
GRID_COLOUR = "#e4e3df"


def log_shift(values):
    """The amount added to susceptibilities so all are > 0 for a log axis:
    -min + 1% of the range when the minimum is <= 0, otherwise 0."""
    vmin = float(np.nanmin(values))
    if vmin > 0:
        return 0.0
    span = float(np.nanmax(values)) - vmin
    return -vmin + (0.01 * span if span > 0 else 1e-6)


def _read_warped(path, wkt, bounds, width, height):
    """Band 1 of a raster resampled (nearest neighbour) onto the target grid,
    NaN where it has no data or does not cover the target."""
    ds = gdal.Warp(
        "", path, format="MEM", dstSRS=wkt, outputBounds=bounds,
        width=width, height=height, resampleAlg="near",
        outputType=gdal.GDT_Float64, dstNodata=np.nan,
    )
    if ds is None:
        raise ValueError(f"Could not resample {path}")
    arr = ds.GetRasterBand(1).ReadAsArray().astype(float)
    ds = None
    return arr


def _cell_area_in(path, wkt):
    """Cell area of a raster once expressed in the CRS given by wkt."""
    vrt = gdal.Warp("", path, format="VRT", dstSRS=wkt)
    gt = vrt.GetGeoTransform()
    vrt = None
    return abs(gt[1] * gt[5])


def paired_grids(density_path, sus_path):
    """Apparent density and susceptibility on one common grid: that of the
    coarser of the two, each of its cells taking the nearest cell of the finer
    grid (reprojected if the CRSs differ).

    Returns (density, susceptibility, geotransform, wkt); cells where either
    grid has no data are NaN.
    """
    d_ds = gdal.Open(density_path)
    s_ds = gdal.Open(sus_path)
    if d_ds is None or s_ds is None:
        raise ValueError("Could not open the apparent density or susceptibility grid")
    d_gt = d_ds.GetGeoTransform()
    d_area = abs(d_gt[1] * d_gt[5])
    s_area = _cell_area_in(sus_path, d_ds.GetProjection())
    target = s_ds if s_area > d_area else d_ds
    gt = target.GetGeoTransform()
    wkt = target.GetProjection()
    w, h = target.RasterXSize, target.RasterYSize
    xs = (gt[0], gt[0] + w * gt[1])
    ys = (gt[3], gt[3] + h * gt[5])
    bounds = (min(xs), min(ys), max(xs), max(ys))
    d_ds = s_ds = target = None
    density = _read_warped(density_path, wkt, bounds, w, h)
    sus = _read_warped(sus_path, wkt, bounds, w, h)
    out_gt = (bounds[0], abs(gt[1]), 0.0, bounds[3], 0.0, -abs(gt[5]))
    return density, sus, out_gt, wkt


def rasterize_categories(polygons, wkt, gt, shape):
    """Category index of each cell (-1 outside every polygon).

    polygons: (geometry WKT in the target CRS, category index) pairs.
    """
    srs = osr.SpatialReference()
    srs.ImportFromWkt(wkt)
    # vector "MEM" from GDAL 3.11, which deprecates "Memory"
    driver = ogr.GetDriverByName("MEM") or ogr.GetDriverByName("Memory")
    src = driver.CreateDataSource("")
    layer = src.CreateLayer("categories", srs, ogr.wkbUnknown)
    layer.CreateField(ogr.FieldDefn("cid", ogr.OFTInteger))
    for geom_wkt, cid in polygons:
        feature = ogr.Feature(layer.GetLayerDefn())
        feature.SetField("cid", int(cid))
        feature.SetGeometry(ogr.CreateGeometryFromWkt(geom_wkt))
        layer.CreateFeature(feature)
        feature = None
    rows, cols = shape
    ds = gdal.GetDriverByName("MEM").Create("", cols, rows, 1, gdal.GDT_Int32)
    ds.SetGeoTransform(gt)
    ds.SetProjection(wkt)
    band = ds.GetRasterBand(1)
    band.Fill(-1)
    gdal.RasterizeLayer(ds, [1], layer, options=["ATTRIBUTE=cid"])
    out = band.ReadAsArray()
    ds = src = None
    return out


# ---------------------------------------------------------------- plotting
def _style_axes(ax, x_label, y_label):
    ax.set_yscale("log")
    ax.set_xlabel(x_label, color=TEXT_SECONDARY)
    ax.set_ylabel(y_label, color=TEXT_SECONDARY)
    ax.tick_params(colors=TEXT_SECONDARY, labelsize=8)
    ax.grid(True, which="major", color=GRID_COLOUR, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def _sus_label(name, shift):
    return f"{name} (SI)" if shift == 0 else f"{name} (SI) + {shift:.3g}"


def _limits(density, sus, robust):
    """Shared axis limits: the 0.1-99.9 percentile range for grids (FFT edge
    effects), the full range for measured points."""
    lo, hi = (0.1, 99.9) if robust else (0, 100)
    x0, x1 = np.percentile(density, [lo, hi])
    y0, y1 = np.percentile(sus, [lo, hi])
    if x1 <= x0:
        x0, x1 = x0 - 0.01, x1 + 0.01
    if y1 <= y0:
        y0, y1 = y0 / 2, y1 * 2
    return (x0, x1), (y0, y1)


def _bin_edges(xlim, ylim, bins):
    return (
        np.linspace(xlim[0], xlim[1], bins + 1),
        np.logspace(math.log10(ylim[0]), math.log10(ylim[1]), bins + 1),
    )


def _heat(ax, density, sus, edges, norm, weights=None):
    from matplotlib.colors import LinearSegmentedColormap

    counts, _, _ = np.histogram2d(density, sus, bins=edges, weights=weights)
    cmap = LinearSegmentedColormap.from_list("sgt_heat", HEAT_RAMP)
    return ax.pcolormesh(
        edges[0], edges[1], np.ma.masked_less_equal(counts.T, 0),
        cmap=cmap, norm=norm, shading="flat",
    )


def heat_bins(n):
    """Bins per axis for a heat map of n points."""
    return int(np.clip(math.sqrt(n) / 2, 20, 200))


def henkel_heatmap(density, sus, title, sus_name="Apparent susceptibility",
                   density_name="Apparent density", robust=True, bins=200):
    """One Henkel heat map (count of points per bin, log colour scale)."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    ok = np.isfinite(density) & np.isfinite(sus)
    density, sus = density[ok], sus[ok]
    if density.size == 0:
        raise ValueError("No locations have both a density and a susceptibility value")
    shift = log_shift(sus)
    sus = sus + shift
    xlim, ylim = _limits(density, sus, robust)
    edges = _bin_edges(xlim, ylim, bins)

    fig, ax = plt.subplots(figsize=(7.5, 6))
    mesh = _heat(ax, density, sus, edges, LogNorm())
    _style_axes(ax, f"{density_name} (g/cc)", _sus_label(sus_name, shift))
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_title(f"{title}\n{density.size:,} points", fontsize=11)
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label("Points per bin", color=TEXT_SECONDARY)
    cbar.ax.tick_params(colors=TEXT_SECONDARY, labelsize=8)
    fig.tight_layout()
    return fig


def henkel_category_grid(density, sus, cats, names, title, bins=100):
    """Henkel heat maps of each category that has data, in a grid of
    ceil(sqrt(N)) columns with shared axes. Each panel shows the fraction of its
    category's points per bin, on one colour scale, so panels can be compared
    whatever their size."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    ok = np.isfinite(density) & np.isfinite(sus) & (cats >= 0)
    density, sus, cats = density[ok], sus[ok], cats[ok]
    if density.size == 0:
        raise ValueError("No polygon covers locations with both grids defined")
    # only categories with data get a panel: renumber them 0..N-1
    present = np.unique(cats)
    names = [names[c] for c in present]
    cats = np.searchsorted(present, cats)
    shift = log_shift(sus)
    sus = sus + shift
    xlim, ylim = _limits(density, sus, True)
    edges = _bin_edges(xlim, ylim, bins)

    # one colour scale for all panels, from the smallest to the largest fraction
    counts = [int(np.count_nonzero(cats == cid)) for cid in range(len(names))]
    fractions = [
        np.histogram2d(density[cats == cid], sus[cats == cid], bins=edges)[0] / n
        for cid, n in enumerate(counts) if n
    ]
    vmax = max(float(f.max()) for f in fractions)
    vmin = min(float(f[f > 0].min()) for f in fractions if (f > 0).any())
    norm = LogNorm(vmin=vmin, vmax=max(vmax, vmin * 1.01))

    # ceil(sqrt(N)) columns, and only as many rows as are needed
    side = math.ceil(math.sqrt(len(names)))
    rows = math.ceil(len(names) / side)
    fig, axes = plt.subplots(
        rows, side, figsize=(3.2 * side + 1.2, 3.0 * rows + 1.0),
        sharex=True, sharey=True, squeeze=False,
    )
    mesh = None
    for i, ax in enumerate(axes.flat):
        if i >= len(names):
            ax.set_visible(False)
            # the panel above is the last in its column: give it the x labels
            axes.flat[i - side].xaxis.set_tick_params(labelbottom=True)
            continue
        sel = cats == i
        n = counts[i]
        mesh = _heat(ax, density[sel], sus[sel], edges, norm,
                     weights=np.full(n, 1.0 / n))
        _style_axes(ax, "", "")
        ax.set_title(f"{names[i]}  (n={n:,})", fontsize=9)
    axes.flat[0].set_xlim(xlim)
    axes.flat[0].set_ylim(ylim)
    fig.supxlabel("Apparent density (g/cc)", color=TEXT_SECONDARY)
    fig.supylabel(_sus_label("Apparent susceptibility", shift), color=TEXT_SECONDARY)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout(rect=(0, 0, 0.9, 1))
    if mesh is not None:
        cax = fig.add_axes((0.92, 0.15, 0.015, 0.7))
        cbar = fig.colorbar(mesh, cax=cax)
        cbar.set_label("Fraction of category's pixels per bin", color=TEXT_SECONDARY)
        cbar.ax.tick_params(colors=TEXT_SECONDARY, labelsize=8)
    return fig


def henkel_points(density, sus, title, categories=None):
    """Henkel plot of measured samples, coloured by category if given."""
    import matplotlib.pyplot as plt

    ok = np.isfinite(density) & np.isfinite(sus)
    density, sus = density[ok], sus[ok]
    if categories is not None:
        categories = np.asarray(categories, dtype=object)[ok]
    if density.size == 0:
        raise ValueError("No samples have both a density and a susceptibility value")
    shift = log_shift(sus)
    sus = sus + shift
    xlim, ylim = _limits(density, sus, False)
    pad = 0.03 * (xlim[1] - xlim[0])

    fig, ax = plt.subplots(figsize=(8.5 if categories is not None else 7, 6))
    style = dict(s=24, linewidths=0.4, edgecolors="white", alpha=0.85)
    if categories is None:
        ax.scatter(density, sus, color=CATEGORY_COLOURS[0], marker="o", **style)
    else:
        names = sorted({str(c) for c in categories}, key=_natural_key)
        labels = np.array([str(c) for c in categories], dtype=object)
        shown = names[:len(CATEGORY_COLOURS)]
        if len(names) > len(CATEGORY_COLOURS):
            # more categories than distinguishable colours: the largest keep one
            counts = {n: int(np.count_nonzero(labels == n)) for n in names}
            shown = sorted(
                sorted(names, key=lambda n: -counts[n])[:len(CATEGORY_COLOURS) - 1],
                key=_natural_key,
            )
        for i, name in enumerate(shown):
            sel = labels == name
            ax.scatter(density[sel], sus[sel], color=CATEGORY_COLOURS[i],
                       marker=CATEGORY_MARKERS[i], label=f"{name} ({sel.sum()})", **style)
        rest = ~np.isin(labels, shown)
        if rest.any():
            ax.scatter(density[rest], sus[rest], color=OTHER_COLOUR, marker=".",
                       label=f"Other ({rest.sum()})", **style)
        ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False,
                  fontsize=8, labelcolor=TEXT_SECONDARY)
    _style_axes(ax, "Density (g/cc)", _sus_label("Susceptibility", shift))
    ax.set_xlim(xlim[0] - pad, xlim[1] + pad)
    ax.set_ylim(ylim[0] / 1.3, ylim[1] * 1.3)
    ax.set_title(f"{title}\n{density.size:,} samples", fontsize=11)
    fig.tight_layout()
    return fig


NO_CATEGORY = "(none)"


def _natural_key(text):
    """Sort '2' before '10', numbers before words, samples without a category last."""
    if text == NO_CATEGORY:
        return (2, 0.0, "")
    try:
        return (0, float(text), "")
    except ValueError:
        return (1, 0.0, text.lower())


def save_figure(fig, base_path, formats):
    """Save fig as base_path.<fmt> for each format ('png', 'pdf'); returns the paths."""
    paths = []
    for fmt in formats:
        path = f"{base_path}.{fmt}"
        fig.savefig(path, dpi=200 if fmt == "png" else None, bbox_inches="tight")
        paths.append(path)
    return paths
