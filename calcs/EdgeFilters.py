"""
Windowed edge detectors and edge-preserving smoothers (no QGIS dependency).

    nstd                  Normalised standard deviation (Cooper & Cowan 2008)
    kuwahara              Kuwahara filter: mean of the least-variance quadrant
    snn                   Symmetric nearest neighbour (Harwood et al. 1987)
    nav_edge              Normalised anisotropy variance (Zhang et al. 2014)

All work on north-up 2D float arrays with NaN for no data, and return a grid
of the same shape with NaN where the input was NaN. Window sizes are in cells.
"""

import numpy as np
from scipy.ndimage import correlate, distance_transform_edt
from scipy.signal import fftconvolve


# ---------------------------------------------------------------- helpers
def _check(callback, fraction):
    if callback is not None:
        callback(fraction)


def _fill_nearest(grid):
    """NaNs replaced by the nearest valid value (for filters that need a full grid)."""
    mask = ~np.isfinite(grid)
    if not mask.any() or mask.all():
        return np.where(mask, 0.0, grid), mask
    idx = distance_transform_edt(mask, return_distances=False, return_indices=True)
    filled = grid[idx[0], idx[1]]
    return filled, mask


def _box_sum(x, k):
    """Sum of every k x k box of x (summed-area table); result[i, j] is the box
    whose top-left cell is x[i, j], shape (rows - k + 1, cols - k + 1)."""
    c = np.zeros((x.shape[0] + 1, x.shape[1] + 1))
    c[1:, 1:] = x.cumsum(0).cumsum(1)
    return c[k:, k:] - c[:-k, k:] - c[k:, :-k] + c[:-k, :-k]


def _box_moments(grid, k, pad):
    """NaN-aware count, sum and sum of squares of every k x k box of `grid`
    padded by `pad` cells (padding counts as no data). The grid is centred on
    its mean first so the summed-area tables keep their precision."""
    valid = np.isfinite(grid)
    offset = np.nanmean(grid) if valid.any() else 0.0
    a = np.where(valid, grid - offset, 0.0)
    a = np.pad(a, pad)
    v = np.pad(valid.astype(float), pad)
    return _box_sum(v, k), _box_sum(a, k), _box_sum(a * a, k), offset


def windowed_std(grid, window):
    """NaN-aware standard deviation in a square window centred on each cell."""
    window = max(1, int(window)) | 1
    r = window // 2
    n, s1, s2, _ = _box_moments(grid, window, r)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = s1 / n
        var = s2 / n - mean * mean
    std = np.sqrt(np.clip(var, 0.0, None))
    std[n < 2] = np.nan
    return std


# ---------------------------------------------------------------- NSTD
def nstd(fx, fy, fz, window=5):
    """Normalised standard deviation (Cooper & Cowan 2008, Geophysics 73(3),
    H1-H4):

        NSTD = s(dz) / (s(dx) + s(dy) + s(dz))

    with s the standard deviation of each derivative in a moving window.
    Edges are the maxima; the ratio balances strong and weak anomalies."""
    sx = windowed_std(fx, window)
    sy = windowed_std(fy, window)
    sz = windowed_std(fz, window)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = sz / (sx + sy + sz)
    out[~np.isfinite(fz)] = np.nan
    return out


# ---------------------------------------------------------------- Kuwahara
def _least_variance_mean(grid, sub, offsets, callback=None):
    """For every cell, the mean of the sub x sub window (among those whose
    top-left corner sits at the given offsets from cell - (sub - 1)) with the
    least variance. NaNs are ignored inside each window."""
    pad = sub - 1
    n, s1, s2, offset = _box_moments(grid, sub, pad)
    rows, cols = grid.shape
    best_var = np.full(grid.shape, np.inf)
    best_mean = np.full(grid.shape, np.nan)
    total = len(offsets)
    for k, (i, j) in enumerate(offsets):
        _check(callback, k / total)
        nn = n[i:i + rows, j:j + cols]
        ss1 = s1[i:i + rows, j:j + cols]
        ss2 = s2[i:i + rows, j:j + cols]
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = ss1 / nn
            var = ss2 / nn - mean * mean
        var[nn < 1] = np.inf
        better = var < best_var
        best_var[better] = var[better]
        best_mean[better] = mean[better]
    out = best_mean + offset
    out[~np.isfinite(grid)] = np.nan
    return out


def kuwahara(grid, window=5, callback=None):
    """Kuwahara filter: the window (odd size 2a+1) is split into four
    overlapping (a+1) x (a+1) quadrants that share the centre cell; the cell
    takes the mean of the quadrant with the lowest variance. Smooths the
    inside of domains while keeping steps sharp."""
    window = max(3, int(window) | 1)
    a = window // 2
    sub = a + 1
    # top-left corners: quadrants up-left, up-right, down-left, down-right
    offsets = [(0, 0), (0, a), (a, 0), (a, a)]
    return _least_variance_mean(grid, sub, offsets, callback)



# ---------------------------------------------------------------- SNN
def snn(grid, window=3, callback=None):
    """Symmetric nearest neighbour (Harwood et al. 1987): for each pair of
    cells placed symmetrically about the centre of the window, keep the one
    closer in value to the centre cell, and average the kept values. Removes
    spikes and line noise without widening narrow features."""
    window = max(3, int(window) | 1)
    r = window // 2
    rows, cols = grid.shape
    p = np.pad(grid, r, constant_values=np.nan)
    centre = grid
    total = np.zeros(grid.shape)
    count = np.zeros(grid.shape)
    pairs = [
        (dy, dx)
        for dy in range(0, r + 1)
        for dx in range(-r, r + 1)
        if dy > 0 or dx > 0
    ]
    for k, (dy, dx) in enumerate(pairs):
        _check(callback, k / len(pairs))
        a = p[r + dy:r + dy + rows, r + dx:r + dx + cols]
        b = p[r - dy:r - dy + rows, r - dx:r - dx + cols]
        with np.errstate(invalid="ignore"):
            take_a = np.abs(a - centre) <= np.abs(b - centre)
        take_a = np.where(np.isnan(b), True, np.where(np.isnan(a), False, take_a))
        v = np.where(take_a, a, b)
        ok = np.isfinite(v)
        total[ok] += v[ok]
        count[ok] += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(count > 0, total / count, centre)
    out[~np.isfinite(grid)] = np.nan
    return out


# ---------------------------------------------------------------- NAV-Edge
def _corr(padded, kernel):
    """Correlation of the padded grid with an odd kernel, 'valid' part only
    (so the result has the unpadded shape). Direct for small kernels, FFT for
    large ones."""
    r = kernel.shape[0] // 2
    if kernel.shape[0] <= 11:
        out = correlate(padded, kernel, mode="constant", cval=0.0)
        return out[r:-r, r:-r] if r else out
    return fftconvolve(padded, kernel[::-1, ::-1], mode="valid")


def nav_edge(grid, window=7, cof=2.0, n_angles=18, callback=None):
    """NAV-Edge, normalised anisotropy variance (Zhang, Ravat, Marangoni & Hu
    2014, Geophysics 79(3), J43-J53), equations 10, 11, 14 and 15.

    In an M x M window (M = `window`) with Gaussian standard deviation
    sigma = (M - 1)/6:

        NAV = sum[(f - fbar) Q] / ( sqrt(sum (f - fbar)^2) * sqrt(sum Q^2) )

    where fbar is the window mean and Q an anisotropic Gaussian narrowed by
    the anisotropy scale `cof` across the local gradient direction (which is
    found from Gaussian-derivative filters, eq. 15) and so stretched along
    strike. cof = 1 is isotropic. NAV is dimensionless, positive over the
    source and negative outside it; edges are the zero crossings.

    The direction is quantised into `n_angles` bins over 180 degrees (one
    kernel per bin); each cell takes the response of its own bin.
    """
    M = max(3, int(window) | 1)
    r = M // 2
    sigma = (M - 1) / 6.0
    cof = max(float(cof), 1e-6)

    filled, mask = _fill_nearest(np.asarray(grid, dtype=float))
    fp = np.pad(filled, r, mode="reflect")

    v, u = np.mgrid[-r:r + 1, -r:r + 1].astype(float)  # v: row offset, u: column offset
    inside = (u * u + v * v) <= r * r + 1e-9  # points beyond 3 sigma are discarded

    # gradient direction from Gaussian-derivative filters (eq. 15)
    G = np.exp(-(u * u + v * v) / (2 * sigma * sigma)) * inside
    gc = _corr(fp, -u / sigma**2 * G)  # d/d column
    gr = _corr(fp, -v / sigma**2 * G)  # d/d row
    phi = np.mod(np.arctan2(gr, gc), np.pi)
    bins = np.minimum((phi / np.pi * n_angles).astype(int), n_angles - 1)
    _check(callback, 0.05)

    # window mean and sqrt(sum (f - fbar)^2) over the full M x M square
    s1 = _box_sum(fp, M)
    s2 = _box_sum(fp * fp, M)
    fbar = s1 / (M * M)
    spread = np.sqrt(np.clip(s2 - M * M * fbar * fbar, 0.0, None))

    out = np.zeros(filled.shape)
    for b in range(n_angles):
        _check(callback, 0.05 + 0.95 * b / n_angles)
        sel = bins == b
        if not sel.any():
            continue
        t = (b + 0.5) * np.pi / n_angles  # gradient direction of this bin
        along_grad = u * np.cos(t) + v * np.sin(t)
        along_strike = -u * np.sin(t) + v * np.cos(t)
        Q = np.exp(
            -(along_strike**2 + cof**2 * along_grad**2) / (2 * sigma * sigma)
        ) * inside
        Q /= Q.sum()
        num = _corr(fp, Q) - fbar * Q.sum()
        den = spread * np.sqrt((Q * Q).sum())
        with np.errstate(invalid="ignore", divide="ignore"):
            nav = np.where(den > 0, num / den, 0.0)
        out[sel] = nav[sel]
    out[mask] = np.nan
    return out
