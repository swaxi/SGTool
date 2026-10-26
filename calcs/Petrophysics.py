"""
Henkel plots (density against log susceptibility) and other scatter plots, heat
maps and histograms for apparent petrophysics grids and measured point data,
with no QGIS dependency.
"""

import math

import numpy as np
from osgeo import gdal, ogr, osr

# Colour map for the heat maps (viridis: its changing hue shows the variation
# in the densest part more clearly than a one-hue ramp), and the fixed
# categorical order for category colours, combined with marker shapes so
# identity never rests on colour alone (8 x 8 = 64 symbols).
HEAT_CMAP = "viridis"
CATEGORY_COLOURS = [
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100",
    "#e87ba4", "#008300", "#6250d6", "#e34948",
]
CATEGORY_MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]
OTHER_COLOUR = "#8a8a85"
TEXT_SECONDARY = "#5f5e5a"
GRID_COLOUR = "#e4e3df"


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

    Either path may be None, when the other grid is read on its own grid (and
    None is returned in its place).

    Returns (density, susceptibility, geotransform, wkt); cells where either
    grid has no data are NaN.
    """
    d_ds = gdal.Open(density_path) if density_path else None
    s_ds = gdal.Open(sus_path) if sus_path else None
    if (density_path and d_ds is None) or (sus_path and s_ds is None):
        raise ValueError("Could not open the apparent density or susceptibility grid")
    if d_ds is None or s_ds is None:
        target = d_ds or s_ds
    else:
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
    density = _read_warped(density_path, wkt, bounds, w, h) if density_path else None
    sus = _read_warped(sus_path, wkt, bounds, w, h) if sus_path else None
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
# Two-variable plots (heat maps and scatter plots) and one-variable histograms.
# Each axis is linear or log; on a log axis values <= 0 can't be shown and are
# left out (their number is given in the title). Henkel plots are the case
# x = density (linear), y = susceptibility (log).
def valid_points(x, y, x_log, y_log, *others):
    """The points with both values defined (and > 0 on a log axis). Returns
    the kept arrays (x, y, *others) and how many were left out for being <= 0."""
    ok = np.isfinite(x) & np.isfinite(y)
    keep = ok.copy()
    if x_log:
        keep &= x > 0
    if y_log:
        keep &= y > 0
    dropped = int(np.count_nonzero(ok & ~keep))
    kept = [x[keep], y[keep]] + [o[keep] for o in others]
    if kept[0].size == 0:
        raise ValueError(
            "No points have both values" + (" above 0 on the log axes" if x_log or y_log else "")
        )
    return kept, dropped


def valid_values(values, log, *others):
    """The defined values (and > 0 for a log axis), with *others, and how many
    were left out for being <= 0."""
    ok = np.isfinite(values)
    keep = ok & (values > 0) if log else ok
    dropped = int(np.count_nonzero(ok & ~keep))
    kept = [values[keep]] + [o[keep] for o in others]
    if kept[0].size == 0:
        raise ValueError("No values" + (" above 0" if log else "") + " to plot")
    return kept, dropped


def _count_text(n, dropped, what="points"):
    text = f"{n:,} {what}"
    if dropped:
        text += f" ({dropped:,} with values <= 0 on a log axis left out)"
    return text


def _style_axes(ax, x_label, y_label, x_log=False, y_log=False, grid="both"):
    if x_log:
        ax.set_xscale("log")
    if y_log:
        ax.set_yscale("log")
    ax.set_xlabel(x_label, color=TEXT_SECONDARY)
    ax.set_ylabel(y_label, color=TEXT_SECONDARY)
    ax.tick_params(colors=TEXT_SECONDARY, labelsize=8)
    ax.grid(True, axis=grid, which="major", color=GRID_COLOUR, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def _axis_range(values, log, robust):
    """Axis range: the 0.1-99.9 percentiles for grids (FFT edge effects), the
    full range for measured samples."""
    lo, hi = np.percentile(values, [0.1, 99.9] if robust else [0, 100])
    if hi <= lo:
        lo, hi = (lo / 2, hi * 2) if log else (lo - 0.01, hi + 0.01)
    return lo, hi


def _edges(lims, log, bins):
    if log:
        return np.logspace(math.log10(lims[0]), math.log10(lims[1]), bins + 1)
    return np.linspace(lims[0], lims[1], bins + 1)


def _padded(lims, log):
    """Axis limits with a margin so points at the edges aren't clipped."""
    if log:
        return lims[0] / 1.3, lims[1] * 1.3
    pad = 0.03 * (lims[1] - lims[0])
    return lims[0] - pad, lims[1] + pad


def _heat(ax, x, y, edges, norm, weights=None):
    counts, _, _ = np.histogram2d(x, y, bins=edges, weights=weights)
    return ax.pcolormesh(
        edges[0], edges[1], np.ma.masked_less_equal(counts.T, 0),
        cmap=HEAT_CMAP, norm=norm, shading="flat",
    )


def heat_bins(n):
    """Bins per axis for a heat map of n points."""
    return int(np.clip(math.sqrt(n) / 2, 20, 200))


# The QGIS plugin shows plots in pyplot windows. The ArcGIS Pro GUI turns this
# off: it builds figures in a worker thread and shows them in its own Tk windows,
# and pyplot figures belong to the thread that made them.
USE_PYPLOT = True


def _subplots(*args, **kwargs):
    """A figure and its axes, made with pyplot or (USE_PYPLOT off) without."""
    if USE_PYPLOT:
        import matplotlib.pyplot as plt

        return plt.subplots(*args, **kwargs)
    from matplotlib.figure import Figure

    fig = Figure(figsize=kwargs.pop("figsize", None))
    return fig, fig.subplots(*args, **kwargs)


MAX_COLUMNS = 8  # panels per row in a grid of categories


def _panel_scale(columns):
    """Panels (and their titles) shrink as the number of columns grows."""
    return 1.0 if columns <= 3 else 0.75 if columns <= 5 else 0.62


def _grid_of_panels(n, width, height):
    """ceil(sqrt(n)) columns (at most MAX_COLUMNS) and only as many rows as
    needed, with shared axes; panels smaller when there are many columns."""
    side = min(math.ceil(math.sqrt(n)), MAX_COLUMNS)
    rows = math.ceil(n / side)
    scale = _panel_scale(side)
    fig, axes = _subplots(
        rows, side, figsize=(width * scale * side + 1.2, height * scale * rows + 1.0),
        sharex=True, sharey=True, squeeze=False,
    )
    for i in range(n, rows * side):
        axes.flat[i].set_visible(False)
        # the panel above is the last in its column: give it the x labels
        axes.flat[i - side].xaxis.set_tick_params(labelbottom=True)
    return fig, axes


def _panel_title(ax, name, n, columns, group=None):
    """The category name (after its group letter, if grouped) wrapped onto at most
    two lines, with its count beneath, in a font that shrinks with the number of
    columns, so titles don't overlap."""
    import textwrap

    width = 30 if columns <= 3 else 24 if columns <= 5 else 18
    if group is not None:
        name = f"{group} | {name}"
    lines = textwrap.wrap(name, width) or [name]
    if len(lines) > 2:
        lines = [lines[0], textwrap.shorten(" ".join(lines[1:]), width, placeholder="...")]
    ax.set_title("\n".join(lines) + f"\n(n={n:,})",
                 fontsize=9 if columns <= 3 else 8 if columns <= 5 else 7, linespacing=1.1)


# ------------------------------------------------ grouping similar categories
GROUP_MIN_CELLS = 30      # grid cells a category needs to be grouped
GROUP_MIN_SAMPLES = 10    # measured samples a category needs to be grouped
MAX_GROUPS = 8
UNGROUPED = "-"


def _js_distance(p, q):
    """Jensen-Shannon distance (base 2: 0 identical, 1 no overlap) of two
    distributions that each sum to 1."""
    m = (p + q) / 2

    def kl(a, b):
        nz = a > 0
        return float(np.sum(a[nz] * np.log2(a[nz] / b[nz])))

    return math.sqrt(max(0.0, (kl(p, m) + kl(q, m)) / 2))


def group_distributions(distributions, counts, n_groups=0, min_count=GROUP_MIN_CELLS):
    """Group categories whose distributions are alike.

    distributions: for each category its fraction of points per bin (1-D or 2-D,
    on bins shared by all categories). They are lightly smoothed, compared with
    the Jensen-Shannon distance and grouped by average-linkage clustering.
    n_groups: 0 chooses the number from the largest jump in the clustering (one
    group if all are alike), otherwise 2 to MAX_GROUPS.
    Returns a label per category: "A", "B" ... in order of the group's mean
    position along the first (x) axis, or UNGROUPED for categories with fewer
    than min_count points.
    """
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.ndimage import gaussian_filter

    labels = [UNGROUPED] * len(distributions)
    use = [i for i, n in enumerate(counts) if n >= min_count]
    if len(use) < 2:
        for i in use:
            labels[i] = "A"
        return labels
    smooth = []
    for i in use:
        s = gaussian_filter(np.asarray(distributions[i], dtype=float), sigma=1.0, mode="constant")
        smooth.append(s / s.sum() if s.sum() > 0 else s)
    n = len(use)
    pairs = [_js_distance(smooth[a], smooth[b]) for a in range(n) for b in range(a + 1, n)]
    tree = linkage(np.array(pairs), method="average")
    heights = tree[:, 2]
    if n_groups:
        k = int(n_groups)
    elif heights[-1] < 0.25:      # even the least alike are close: one group
        k = 1
    elif n == 2:
        k = 2
    else:
        # cut in the largest jump between merges: below merge j + 1, n - j - 1 groups
        k = n - int(np.argmax(np.diff(heights))) - 1
    k = int(np.clip(k, 1, min(MAX_GROUPS, n)))
    cluster = fcluster(tree, k, criterion="maxclust")

    # letters in order of the groups' mean position along the x axis
    position = np.arange(smooth[0].shape[0])
    centres = {}
    for c, s in zip(cluster, smooth):
        profile = s.sum(axis=tuple(range(1, s.ndim))) if s.ndim > 1 else s
        centres.setdefault(c, []).append(float((position * profile).sum()))
    order = sorted(centres, key=lambda c: np.mean(centres[c]))
    letter = {c: "ABCDEFGH"[j] for j, c in enumerate(order)}
    for i, c in zip(use, cluster):
        labels[i] = letter[c]
    return labels


def group_colour(label):
    return OTHER_COLOUR if label == UNGROUPED else CATEGORY_COLOURS["ABCDEFGH".index(label)]


def _group_order(names, labels):
    """Panel order: by group, then by name; ungrouped last."""
    rank = lambda i: (labels[i] == UNGROUPED, labels[i], _natural_key(names[i]))
    return sorted(range(len(names)), key=rank)


def _frame(ax, label):
    """Outline a panel in its group's colour."""
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_edgecolor(group_colour(label))
        spine.set_linewidth(1.6)


def _group_note(labels, what, min_count):
    groups = sorted({g for g in labels if g != UNGROUPED})
    text = (f"Groups {groups[0]}-{groups[-1]}" if len(groups) > 1 else f"Group {groups[0]}") \
        if groups else "No groups"
    text += (": categories with similar distributions (Jensen-Shannon distance on the "
             "plotted bins, average-linkage clustering), lettered along the x axis.")
    if UNGROUPED in labels:
        text += f"\n{UNGROUPED} = fewer than {min_count} {what}, not grouped."
    return text


def _bottom_notes(fig, x_label, notes, right=1.0):
    """Notes under the x axis label (a figure-wide label, if given), with room
    made for them."""
    lines = sum(n.count("\n") + 1 for n in notes)
    bottom = 0.02 + 0.022 * lines
    fig.tight_layout(rect=(0, bottom + (0.03 if x_label else 0.0), right, 1))
    if x_label:
        fig.supxlabel(x_label, color=TEXT_SECONDARY, y=bottom, va="bottom")
    fig.text(0.5, 0.01, "\n".join(notes), ha="center", va="bottom", fontsize=8,
             color=TEXT_SECONDARY)


def _by_category(cats, names):
    """Renumber categories 0..N-1, keeping only those with values."""
    present = np.unique(cats)
    return [names[c] for c in present], np.searchsorted(present, cats)


def heatmap(x, y, title, x_label, y_label, x_log=False, y_log=True,
            robust=True, bins=200, what="points"):
    """One heat map of y against x (points per bin, log colour scale)."""
    from matplotlib.colors import LogNorm

    (x, y), dropped = valid_points(x, y, x_log, y_log)
    xlim, ylim = _axis_range(x, x_log, robust), _axis_range(y, y_log, robust)
    edges = (_edges(xlim, x_log, bins), _edges(ylim, y_log, bins))

    fig, ax = _subplots(figsize=(7.5, 6))
    mesh = _heat(ax, x, y, edges, LogNorm())
    _style_axes(ax, x_label, y_label, x_log, y_log)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_title(f"{title}\n{_count_text(x.size, dropped, what)}", fontsize=11)
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label(f"{what.capitalize()} per bin", color=TEXT_SECONDARY)
    cbar.ax.tick_params(colors=TEXT_SECONDARY, labelsize=8)
    fig.tight_layout()
    return fig


def _position_colours(edges, x_log, y_log):
    """RGB of each bin from where it is on the plot: red to green along x, red
    to blue up y (low/low red, high x green, high y blue, high/high cyan)."""
    def fraction(e, log):
        c = np.log10(e) if log else e
        mid = (c[:-1] + c[1:]) / 2
        return (mid - c[0]) / (c[-1] - c[0])

    u, v = np.meshgrid(fraction(edges[0], x_log), fraction(edges[1], y_log))
    return np.stack([(1 - u) * (1 - v), u, v], axis=-1)  # (y bins, x bins, 3)


def _position_heat(ax, x, y, edges, norm, rgb, weights):
    """A heat map coloured by bin position, with the number of points per bin
    shown by opacity (log scale; empty bins transparent)."""
    counts = np.histogram2d(x, y, bins=edges, weights=weights)[0].T
    filled = counts > 0
    level = np.clip(np.asarray(norm(np.where(filled, counts, norm.vmin))), 0, 1)
    alpha = np.where(filled, 0.15 + 0.85 * level, 0.0)
    mesh = ax.pcolormesh(edges[0], edges[1], np.zeros_like(counts), shading="flat")
    mesh.set_array(None)  # colours set directly, not from a colour map
    mesh.set_facecolor(np.dstack([rgb, alpha]).reshape(-1, 4))
    return mesh


def heatmap_by_category(x, y, cats, names, title, x_label, y_label, x_log=False,
                        y_log=True, robust=True, bins=100, what="cells",
                        position_colours=False, grouping=None, min_count=GROUP_MIN_CELLS):
    """Heat maps of each category that has data, in a grid of panels sharing
    their axes and colour scale. Each shows the fraction of its category's
    points per bin, so panels compare whatever their size.
    cats: index into names of each point (-1 for none).
    position_colours: colour each bin by its position (x red to green, y red to
    blue) and show the fraction per bin by opacity instead of a colour scale.
    grouping: None for no grouping, 0 to group similar categories automatically,
    or the number of groups (see group_distributions)."""
    from matplotlib.colors import LogNorm

    inside = cats >= 0
    if not inside.any():
        raise ValueError("No polygon covers the grids")
    (x, y, cats), dropped = valid_points(x[inside], y[inside], x_log, y_log, cats[inside])
    names, cats = _by_category(cats, names)
    xlim, ylim = _axis_range(x, x_log, robust), _axis_range(y, y_log, robust)
    edges = (_edges(xlim, x_log, bins), _edges(ylim, y_log, bins))

    # one colour scale for all panels, from the smallest to the largest fraction
    counts = [int(np.count_nonzero(cats == i)) for i in range(len(names))]
    fractions = [
        np.histogram2d(x[cats == i], y[cats == i], bins=edges)[0] / n
        for i, n in enumerate(counts)
    ]
    vmax = max(float(f.max()) for f in fractions)
    vmin = min(float(f[f > 0].min()) for f in fractions if (f > 0).any())
    norm = LogNorm(vmin=vmin, vmax=max(vmax, vmin * 1.01))

    labels = (group_distributions(fractions, counts, grouping, min_count)
              if grouping is not None else None)
    order = _group_order(names, labels) if labels else range(len(names))

    rgb = _position_colours(edges, x_log, y_log) if position_colours else None
    fig, axes = _grid_of_panels(len(names), 3.2, 3.0)
    for slot, i in enumerate(order):
        ax = axes.flat[slot]
        sel = cats == i
        weights = np.full(counts[i], 1.0 / counts[i])
        if position_colours:
            mesh = _position_heat(ax, x[sel], y[sel], edges, norm, rgb, weights)
        else:
            mesh = _heat(ax, x[sel], y[sel], edges, norm, weights=weights)
        _style_axes(ax, "", "", x_log, y_log)
        _panel_title(ax, names[i], counts[i], axes.shape[1], labels[i] if labels else None)
        if labels:
            _frame(ax, labels[i])
    axes.flat[0].set_xlim(xlim)
    axes.flat[0].set_ylim(ylim)
    fig.supxlabel(x_label, color=TEXT_SECONDARY)
    fig.supylabel(y_label, color=TEXT_SECONDARY)
    if dropped:
        title += f"\n{dropped:,} {what} with values <= 0 on a log axis left out"
    fig.suptitle(title, fontsize=11)
    notes = []
    if labels:
        notes.append(_group_note(labels, what, min_count))
    if position_colours:
        notes.append(
            "Colour shows position: red to green along the x axis, red to blue up the "
            f"y axis.\nOpacity shows the fraction of the category's {what} per bin "
            "(log scale): solid where most crowded, faint where sparse.")
        _bottom_notes(fig, x_label, notes)
        return fig
    if notes:
        _bottom_notes(fig, x_label, notes, right=0.9)
    else:
        fig.tight_layout(rect=(0, 0, 0.9, 1))
    cax = fig.add_axes((0.92, 0.15, 0.015, 0.7))
    cbar = fig.colorbar(mesh, cax=cax)
    cbar.set_label(f"Fraction of category's {what} per bin", color=TEXT_SECONDARY)
    cbar.ax.tick_params(colors=TEXT_SECONDARY, labelsize=8)
    return fig


def scatter(x, y, title, x_label, y_label, x_log=False, y_log=True,
            categories=None, legend_title=None, what="samples", grouping=None,
            min_count=GROUP_MIN_SAMPLES):
    """Every point of y against x, with a symbol per category if given.
    grouping: None for no grouping, 0 to group categories with similar
    distributions automatically, or the number of groups; the group letter is
    shown in the legend (see group_distributions)."""
    if categories is None:
        (x, y), dropped = valid_points(x, y, x_log, y_log)
    else:
        (x, y, categories), dropped = valid_points(
            x, y, x_log, y_log, np.asarray(categories, dtype=object)
        )
    xlim, ylim = _axis_range(x, x_log, False), _axis_range(y, y_log, False)

    style = dict(s=24, linewidths=0.4, edgecolors="white", alpha=0.85)
    if categories is None:
        fig, ax = _subplots(figsize=(7, 6))
        ax.scatter(x, y, color=CATEGORY_COLOURS[0], marker="o", **style)
    else:
        # Each category gets a colour + marker pair: 8 x 8 = 64 symbols. They go
        # to the categories with the most samples, in rank order (the 8 largest
        # are circles in the 8 colours, the next 8 squares, ...); any further
        # categories are grouped as Other.
        labels = np.array([str(c) for c in categories], dtype=object)
        names, counts = np.unique(labels, return_counts=True)
        ranked = [str(n) for n in names[np.lexsort((names, -counts))]]
        n_symbols = len(CATEGORY_COLOURS) * len(CATEGORY_MARKERS)
        if len(ranked) > n_symbols:
            ranked = ranked[:n_symbols - 1]  # the last symbol is Other
        symbol = {
            name: (CATEGORY_COLOURS[i % len(CATEGORY_COLOURS)],
                   CATEGORY_MARKERS[i // len(CATEGORY_COLOURS)])
            for i, name in enumerate(ranked)
        }
        group = {}
        if grouping is not None:
            # each category's distribution on a coarse grid of the plotted space
            edges = (_edges(xlim, x_log, 20), _edges(ylim, y_log, 20))
            sizes = [int(np.count_nonzero(labels == n)) for n in names]
            fractions = [np.histogram2d(x[labels == n], y[labels == n], bins=edges)[0] / s
                         for n, s in zip(names, sizes)]
            group = dict(zip((str(n) for n in names),
                             group_distributions(fractions, sizes, grouping, min_count)))
        entries = len(symbol) + (len(names) > len(symbol))
        columns = math.ceil(entries / 24)
        fig, ax = _subplots(figsize=(7 + 1.7 * columns, 6))
        # legend in name order, or by group then name
        key = ((lambda n: (group[n] == UNGROUPED, group[n], _natural_key(n))) if group
               else _natural_key)
        for name in sorted(symbol, key=key):
            sel = labels == name
            colour, marker = symbol[name]
            text = f"{name} ({sel.sum()})"
            ax.scatter(x[sel], y[sel], color=colour, marker=marker,
                       label=f"{group[name]} | {text}" if group else text, **style)
        rest = ~np.isin(labels, list(symbol))
        if rest.any():
            ax.scatter(x[rest], y[rest], color=OTHER_COLOUR, marker=".",
                       label=f"Other ({rest.sum()})", **style)
        legend = ax.legend(
            loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False,
            fontsize=7 if entries > 24 else 8, labelcolor=TEXT_SECONDARY,
            ncol=columns, columnspacing=1.0, handletextpad=0.3,
            title=legend_title, title_fontsize=8,
        )
        legend.get_title().set_color(TEXT_SECONDARY)
    _style_axes(ax, x_label, y_label, x_log, y_log)
    ax.set_xlim(_padded(xlim, x_log))
    ax.set_ylim(_padded(ylim, y_log))
    ax.set_title(f"{title}\n{_count_text(x.size, dropped, what)}", fontsize=11)
    if categories is not None and group:
        _bottom_notes(fig, None, [_group_note(list(group.values()), what, min_count)])
    else:
        fig.tight_layout()
    return fig


MEDIAN_COLOUR = "#d62728"


def _round_up(value):
    """Round up to the next half step of the value's leading digit:
    47 -> 50, 0.137 -> 0.15, 0.083 -> 0.085, 120 -> 150."""
    if value <= 0:
        return 1.0
    step = 10 ** math.floor(math.log10(value)) / 2
    return math.ceil(value / step - 1e-9) * step


def _add_stats(ax, values, fontsize=8):
    """Median as a red vertical line, and median, mean and standard deviation
    as text in the top right corner."""
    median = float(np.median(values))
    ax.axvline(median, color=MEDIAN_COLOUR, linewidth=1.5)
    ax.text(
        0.98, 0.97,
        f"Median {median:.4g}\nMean {float(np.mean(values)):.4g}\nSD {float(np.std(values)):.4g}",
        transform=ax.transAxes, ha="right", va="top", fontsize=fontsize,
        color=TEXT_SECONDARY, linespacing=1.3,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=GRID_COLOUR, alpha=0.85),
    )


def histogram(values, title, label, log=False, robust=True, what="points"):
    """Histogram of one variable, with linear or log bins."""
    (values,), dropped = valid_values(values, log)
    lims = _axis_range(values, log, robust)
    edges = _edges(lims, log, int(np.clip(math.sqrt(values.size), 10, 100)))
    fig, ax = _subplots(figsize=(7.5, 5))
    counts, _, _ = ax.hist(values, bins=edges, color=CATEGORY_COLOURS[0],
                           histtype="stepfilled", linewidth=0)
    _add_stats(ax, values, fontsize=9)
    _style_axes(ax, label, "Count", x_log=log, grid="y")
    ax.set_xlim(edges[0], edges[-1])
    ax.set_ylim(0, _round_up(counts.max()))
    ax.set_title(f"{title}\n{_count_text(values.size, dropped, what)}", fontsize=11)
    fig.tight_layout()
    return fig


def index_categories(labels):
    """Category labels -> (index of each label, sorted category names)."""
    names = sorted({str(label) for label in labels}, key=_natural_key)
    index = {name: i for i, name in enumerate(names)}
    return np.array([index[str(label)] for label in labels], dtype=int), names


def histogram_by_category(values, cats, names, title, label, log=False,
                          robust=True, what="cells", grouping=None,
                          min_count=GROUP_MIN_CELLS):
    """Histograms of one variable for each category that has values, in a grid
    of panels sharing their bins and axes. Each shows the fraction of its
    category's values per bin, so categories compare whatever their size.
    cats: index into names of each value (-1 for none).
    grouping: None for no grouping, 0 to group similar categories automatically,
    or the number of groups (see group_distributions)."""
    inside = cats >= 0
    if not inside.any():
        raise ValueError("No polygon covers the grid")
    (values, cats), dropped = valid_values(values[inside], log, cats[inside])
    names, cats = _by_category(cats, names)
    lims = _axis_range(values, log, robust)
    edges = _edges(lims, log, int(np.clip(math.sqrt(values.size), 10, 100)))

    sizes = [int(np.count_nonzero(cats == i)) for i in range(len(names))]
    labels = None
    if grouping is not None:
        fractions = [np.histogram(values[cats == i], bins=edges)[0] / n
                     for i, n in enumerate(sizes)]
        labels = group_distributions(fractions, sizes, grouping, min_count)
    order = _group_order(names, labels) if labels else range(len(names))

    fig, axes = _grid_of_panels(len(names), 3.2, 2.6)
    top = 0.0
    for slot, i in enumerate(order):
        ax = axes.flat[slot]
        v = values[cats == i]
        colour = group_colour(labels[i]) if labels else CATEGORY_COLOURS[0]
        counts, _, _ = ax.hist(v, bins=edges, weights=np.full(v.size, 1.0 / v.size),
                               color=colour, histtype="stepfilled", linewidth=0)
        top = max(top, float(counts.max()))
        _add_stats(ax, v, fontsize=7)
        _style_axes(ax, "", "", x_log=log, grid="y")
        _panel_title(ax, names[i], v.size, axes.shape[1], labels[i] if labels else None)
        if labels:
            _frame(ax, labels[i])
    axes.flat[0].set_xlim(edges[0], edges[-1])
    axes.flat[0].set_ylim(0, _round_up(top))  # shared by all panels
    fig.supxlabel(label, color=TEXT_SECONDARY)
    fig.supylabel(f"Fraction of category's {what}", color=TEXT_SECONDARY)
    if dropped:
        title += f"\n{dropped:,} {what} with values <= 0 on a log axis left out"
    fig.suptitle(title, fontsize=11)
    if labels:
        _bottom_notes(fig, label, [_group_note(labels, what, min_count)])
    else:
        fig.tight_layout()
    return fig


NO_CATEGORY = "(none)"


def to_number(value):
    """A numeric or text attribute as a number (text such as ' 2.67' or '2,67'
    too); NaN for NULL or text that isn't a number."""
    if value is None or str(value) == "NULL":
        return np.nan
    try:
        return float(value)
    except (TypeError, ValueError):
        try:
            return float(str(value).strip().replace(",", "."))
        except ValueError:
            return np.nan


def category_label(value):
    """A category attribute as text, NO_CATEGORY for NULL or blank."""
    if value is None or str(value) == "NULL" or not str(value).strip():
        return NO_CATEGORY
    return str(value)


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
