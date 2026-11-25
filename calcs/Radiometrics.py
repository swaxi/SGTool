"""
Radiometric (gamma-ray spectrometry) ratios and transforms, no QGIS dependency.

    efimov_f            Efimov F-parameter, K * eU / eTh
    thorium_normalised  Kd / Ud: deviation of K or eU from the background
                        predicted from eTh (Saunders et al. 1987; Pires 1995)
    clr, ilr            centred / isometric log-ratio transforms of a band stack
    ternary_gradient    multichannel (vector Sobel / Di Zenzo) gradient of K, eU, eTh

and read_band_on_frame, which reads one band of any raster resampled onto a
target grid frame, so K, eU and eTh may come from bands of one multiband
raster or from separate rasters on different grids.
"""

import numpy as np
from scipy.ndimage import distance_transform_edt, sobel

try:
    from osgeo import gdal
except ImportError:  # array calculations only
    from .no_gdal import gdal


# ---------------------------------------------------------------- reading
def north_up_frame(path):
    """(geotransform, cols, rows, wkt) of a raster, as a north-up frame (the
    plugin flips south-up grids to north-up when it reads them)."""
    ds = gdal.Open(path)
    if ds is None:
        raise ValueError(f"Could not open {path}")
    gt = list(ds.GetGeoTransform())
    cols, rows = ds.RasterXSize, ds.RasterYSize
    if gt[5] > 0:
        gt = [gt[0], gt[1], 0.0, gt[3] + rows * gt[5], 0.0, -gt[5]]
    return tuple(gt), cols, rows, ds.GetProjection()


def read_band_on_frame(path, band, gt, cols, rows, wkt):
    """Band `band` (1-based) of the raster at `path` on the north-up frame
    (gt, cols, rows, wkt), as float with NaN for no data. Read directly when
    the raster already is on that frame, otherwise resampled with gdal.Warp
    (block average when the frame is coarser, as for a zoomed-out preview,
    bilinear otherwise)."""
    src = gdal.Open(path)
    if src is None:
        raise ValueError(f"Could not open {path}")
    if band < 1 or band > src.RasterCount:
        raise ValueError(f"{path} has no band {band}")
    sgt = src.GetGeoTransform()
    if (
        src.RasterXSize == cols
        and src.RasterYSize == rows
        and np.allclose(sgt, gt, rtol=0, atol=1e-9 * max(1.0, abs(gt[1])))
    ):
        b = src.GetRasterBand(band)
        arr = b.ReadAsArray().astype(float)
        nd = b.GetNoDataValue()
        if nd is not None:
            arr[arr == nd] = np.nan
        return arr

    vrt = gdal.Translate("", src, format="VRT", bandList=[band])
    coarser = abs(gt[1]) > 1.5 * abs(sgt[1])
    kwargs = dict(
        format="MEM",
        outputBounds=(gt[0], gt[3] + rows * gt[5], gt[0] + cols * gt[1], gt[3]),
        width=cols,
        height=rows,
        resampleAlg="average" if coarser else "bilinear",
        outputType=gdal.GDT_Float32,
        dstNodata=np.nan,
    )
    if wkt and src.GetProjection():
        kwargs["srcSRS"] = src.GetProjection()
        kwargs["dstSRS"] = wkt
    out = gdal.Warp("", vrt, **kwargs)
    if out is None:
        raise ValueError(f"Could not resample {path} onto the grid")
    return out.GetRasterBand(1).ReadAsArray().astype(float)


def band_names(path):
    """Band descriptions of a raster ('' where there is none)."""
    ds = gdal.Open(path)
    if ds is None:
        return []
    return [ds.GetRasterBand(i + 1).GetDescription() or "" for i in range(ds.RasterCount)]


# ---------------------------------------------------------------- ratios
def efimov_f(k, u, th):
    """Efimov (1978) F-parameter, F = K * eU / eTh (equivalently K / (eTh/eU)).
    High where K and U are enriched relative to Th, as in potassic and
    uraniferous hydrothermal alteration. NaN where K, eU or eTh <= 0
    (a negative concentration, e.g. noisy eU, would make F negative, which
    has no physical meaning)."""
    with np.errstate(invalid="ignore", divide="ignore"):
        f = k * u / th
    f[~((k > 0) & (u > 0) & (th > 0))] = np.nan
    return f


def thorium_normalised(x, th, method="ratio", min_fraction=0.1):
    """Thorium-normalised deviation (Saunders et al. 1987; Pires 1995):

        X_s = predicted background of X (K or eU) from eTh
        X_d = (X - X_s) / X_s

    method 'ratio': X_s = (mean X / mean eTh) * eTh (the usual form)
    method 'regression': X_s from the least-squares line X = a + b * eTh
    fitted to the whole grid.

    Positive values are K (Kd) or U (Ud) in excess of what the thorium
    content predicts for the rock type, i.e. candidate alteration.

    Where X_s is below `min_fraction` of the mean of X (very low eTh, or
    where a regression line nears zero) the ratio explodes, so those cells
    are left blank."""
    ok = np.isfinite(x) & np.isfinite(th)
    if ok.sum() < 3:
        return np.full(x.shape, np.nan)
    if method == "regression":
        b, a = np.polyfit(th[ok], x[ok], 1)
        xs = a + b * th
    else:
        mean_th = th[ok].mean()
        if mean_th == 0:
            return np.full(x.shape, np.nan)
        xs = x[ok].mean() / mean_th * th
    floor = max(min_fraction * abs(x[ok].mean()), 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        xd = (x - xs) / xs
        xd[~(xs > floor)] = np.nan
    return xd


# ---------------------------------------------------------------- log-ratios
def _positive(stack):
    """Log-ratios need strictly positive parts: values <= 0 in each band are
    replaced by half the smallest positive value of that band. Returns the
    stack and the number of cells replaced."""
    out = np.array(stack, dtype=float)
    replaced = 0
    for i in range(out.shape[0]):
        b = out[i]
        bad = np.isfinite(b) & (b <= 0)
        if bad.any():
            pos = b[np.isfinite(b) & (b > 0)]
            floor = pos.min() / 2.0 if pos.size else 1e-6
            b[bad] = floor
            replaced += int(bad.sum())
    return out, replaced


def clr(stack):
    """Centred log-ratio (Aitchison 1986): clr_i = ln(x_i) - mean_j ln(x_j).
    `stack` is (bands, rows, cols). Returns (clr stack, cells replaced)."""
    pos, replaced = _positive(stack)
    logs = np.log(pos)
    return logs - logs.mean(axis=0, keepdims=True), replaced


def ilr(stack):
    """Isometric log-ratio, pivot coordinates (Egozcue et al. 2003; Fiserova
    & Hron 2011): for D parts, D - 1 coordinates

        ilr_i = sqrt((D - i) / (D - i + 1)) * ln( x_i / gmean(x_{i+1} ... x_D) )

    so the first coordinate contrasts band 1 with all the others.
    Returns (ilr stack, cells replaced)."""
    pos, replaced = _positive(stack)
    logs = np.log(pos)
    d = logs.shape[0]
    if d < 2:
        raise ValueError("ILR needs at least two bands")
    out = np.empty((d - 1,) + logs.shape[1:])
    for i in range(d - 1):  # i is 0-based, so D - i parts remain
        rest = logs[i + 1:].mean(axis=0)
        out[i] = np.sqrt((d - i - 1) / (d - i)) * (logs[i] - rest)
    return out, replaced


# ---------------------------------------------------------------- gradient
def _fill_nearest(grid):
    mask = ~np.isfinite(grid)
    if not mask.any() or mask.all():
        return np.where(mask, 0.0, grid)
    idx = distance_transform_edt(mask, return_distances=False, return_indices=True)
    return grid[idx[0], idx[1]]


def ternary_gradient(stack, scaling="zscore", dx=1.0, dy=1.0):
    """Multichannel gradient magnitude of a band stack (Di Zenzo 1986):
    Sobel derivatives of every band are summed into a 2 x 2 structure tensor
    and the square root of its largest eigenvalue is the rate of change of
    the whole composition, so boundaries show where any radioelement changes,
    even if the total count does not.

    scaling 'zscore': each band is standardised (mean 0, std 1) first, so
    each element weighs the same; 'clr': the bands are CLR-transformed first,
    so the gradient follows changes in the element ratios only.
    Gradients are per unit distance of the grid (dx, dy)."""
    stack = np.asarray(stack, dtype=float)
    nan = ~np.all(np.isfinite(stack), axis=0)
    if scaling == "clr":
        bands, _ = clr(stack)
    else:
        bands = np.empty_like(stack)
        for i in range(stack.shape[0]):
            b = stack[i]
            sd = np.nanstd(b)
            bands[i] = (b - np.nanmean(b)) / (sd if sd > 0 else 1.0)
    sxx = np.zeros(stack.shape[1:])
    syy = np.zeros(stack.shape[1:])
    sxy = np.zeros(stack.shape[1:])
    for b in bands:
        f = _fill_nearest(b)
        gx = sobel(f, axis=1) / (8.0 * dx)
        gy = -sobel(f, axis=0) / (8.0 * dy)  # rows run south: north-positive
        sxx += gx * gx
        syy += gy * gy
        sxy += gx * gy
    lam = 0.5 * (sxx + syy + np.sqrt((sxx - syy) ** 2 + 4.0 * sxy * sxy))
    out = np.sqrt(lam)
    out[nan] = np.nan
    return out


# ---------------------------------------------------------------- shared table
# The K / eU / eTh calculations as the QGIS dialog, the QGIS Processing
# algorithms and the ArcGIS Pro GUI and toolbox all run and name them:
# id -> (operation recorded in the provenance, output file suffix)
RADIOMETRIC_OPERATIONS = {
    "efimov": ("Efimov F-parameter (K*eU/eTh)", "_Efimov"),
    "kd": ("Thorium-normalised K deviation (Kd)", "_Kd"),
    "ud": ("Thorium-normalised U deviation (Ud)", "_Ud"),
    "ternary_gradient": ("Ternary (multichannel Sobel) gradient", "_TernGrad"),
}
LOG_RATIO_OPERATIONS = {
    "CLR": ("Centred log-ratio", "_CLR"),
    "ILR": ("Isometric log-ratio", "_ILR"),
}


def compute_radiometric(kind, stack, background="ratio", scaling="zscore", dx=1.0, dy=1.0):
    """One of RADIOMETRIC_OPERATIONS from a (3, rows, cols) K, eU, eTh stack.
    background ('ratio' or 'regression') is for Kd / Ud, scaling ('zscore' or
    'clr') for the ternary gradient."""
    k, u, th = stack[0], stack[1], stack[2]
    if kind == "efimov":
        return efimov_f(k, u, th)
    if kind == "kd":
        return thorium_normalised(k, th, background)
    if kind == "ud":
        return thorium_normalised(u, th, background)
    if kind == "ternary_gradient":
        return ternary_gradient(stack, scaling, dx=dx, dy=dy)
    raise ValueError(f"Unknown radiometric calculation: {kind}")


def log_ratio(kind, stack):
    """CLR or ILR of a band stack: (result stack, cells <= 0 replaced)."""
    return (ilr if kind == "ILR" else clr)(stack)


def log_ratio_band_names(kind, names):
    """Band descriptions of a log-ratio result, from the input band names."""
    if kind == "ILR":
        return [f"ILR{i + 1} {names[i]} vs rest" for i in range(len(names) - 1)]
    return [f"CLR {n}" for n in names]


_BAND_HINTS = {
    "K": ("k", "pot", "potassium", "k_pct", "kpct"),
    "U": ("u", "eu", "ura", "uranium", "u_ppm", "eu_ppm"),
    "Th": ("th", "eth", "tho", "thorium", "th_ppm", "eth_ppm"),
}


def guess_band(names, key):
    """1-based band of K ('K'), eU ('U') or eTh ('Th') among band descriptions:
    by name if one matches (K_pct, eU_ppm, Thorium ...), else bands 1, 2, 3 of
    a grid of 3 or more bands, else band 1."""
    import re

    for i, name in enumerate(names):
        low = (name or "").lower()
        tokens = set(re.split(r"[^a-z0-9_]+", low)) | set(re.split(r"[^a-z0-9]+", low))
        if tokens & set(_BAND_HINTS[key]):
            return i + 1
    order = {"K": 1, "Th": 2, "U": 3}[key]  # ternary order K, Th, U
    return order if len(names) >= 3 else 1


def read_stack(sources, gt, cols, rows, wkt):
    """[(path, band), ...] read onto one north-up frame, as a (bands, rows, cols)
    array. Raises ValueError if the same band of the same grid is given twice
    (every ratio would then be constant)."""
    if len(set(sources)) < len(sources):
        raise ValueError(
            "K, eU and eTh point at the same band of the same grid: choose three "
            "different bands (or grids)"
        )
    return np.stack([read_band_on_frame(p, b, gt, cols, rows, wkt) for p, b in sources])


def write_multiband(path, stack, gt, wkt, descriptions=None):
    """A (bands, rows, cols) stack as a float32 GeoTIFF (NaN no data)."""
    rows, cols = stack.shape[1:]
    ds = gdal.GetDriverByName("GTiff").Create(
        path, cols, rows, stack.shape[0], gdal.GDT_Float32,
        options=["COMPRESS=DEFLATE", "TILED=YES"],
    )
    if ds is None:
        raise ValueError(f"Could not create {path}")
    ds.SetGeoTransform(gt)
    if wkt:
        ds.SetProjection(wkt)
    for i in range(stack.shape[0]):
        band = ds.GetRasterBand(i + 1)
        band.SetNoDataValue(np.nan)
        if descriptions:
            band.SetDescription(descriptions[i])
        band.WriteArray(stack[i].astype(np.float32))
    ds.FlushCache()
    ds = None
