"""
Moving-window statistics on a grid, computed for many windows at once.

The windows of a NaN-padded grid are gathered in blocks of rows into an array
holding the values of each window along its last axis, so a statistic over
every window is one vectorised NumPy call rather than a Python call per pixel.
"""

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view


def iter_windows(padded, size, out_shape, footprint=None, max_values=4000000):
    """
    Yield (r0, r1, values) for blocks of output rows r0:r1, where
    values[r, c, :] holds the cells of the size x size window whose top-left
    corner is padded[r0 + r, c] (only the cells where footprint is True, if
    given).
    """
    rows, cols = out_shape
    padded = padded[: rows + size - 1, : cols + size - 1]
    cells = None if footprint is None else np.flatnonzero(np.ravel(footprint))
    n_cells = size * size if cells is None else cells.size
    step = max(1, max_values // max(1, cols * n_cells))
    for r0 in range(0, rows, step):
        r1 = min(r0 + step, rows)
        windows = sliding_window_view(padded[r0 : r1 + size - 1], (size, size))
        values = windows.reshape(r1 - r0, cols, size * size)
        if cells is not None:
            values = values[..., cells]
        yield r0, r1, values


def nanmedian_last(values):
    """
    Median of the non-NaN values along the last axis (NaN where there are
    none), as np.nanmedian gives, from one sort of all the windows.
    """
    count = np.count_nonzero(~np.isnan(values), axis=-1)
    ordered = np.sort(values, axis=-1)  # NaNs sort last
    lo = np.take_along_axis(ordered, np.maximum((count - 1) // 2, 0)[..., None], -1)
    hi = np.take_along_axis(ordered, (count // 2)[..., None], -1)
    lo, hi = lo[..., 0], hi[..., 0]
    median = np.where(count % 2 == 1, hi, (lo + hi) / 2.0)
    median[count == 0] = np.nan
    return median


def masked_stat(values, valid, count, stat_type):
    """
    One statistic of the valid values along the last axis, matching NumPy
    (mean, variance, std, median, min, max) and scipy.stats with its defaults
    (biased skewness, Fisher kurtosis, NaN for a near-constant window).
    """
    if stat_type == "median":
        return nanmedian_last(np.where(valid, values, np.nan))
    if stat_type == "min":
        return np.where(valid, values, np.inf).min(axis=-1)
    if stat_type == "max":
        return np.where(valid, values, -np.inf).max(axis=-1)

    with np.errstate(all="ignore"):
        mean = np.where(valid, values, 0.0).sum(axis=-1) / count
        if stat_type == "mean":
            return mean
        dev = np.where(valid, values - mean[..., None], 0.0)
        m2 = (dev**2).sum(axis=-1) / count
        if stat_type == "variance":
            return m2
        if stat_type == "std":
            return np.sqrt(m2)
        zero = m2 <= (np.finfo(float).eps * mean) ** 2
        if stat_type == "skewness":
            m3 = (dev**2 * dev).sum(axis=-1) / count
            return np.where(zero, np.nan, m3 / m2**1.5)
        if stat_type == "kurtosis":
            m4 = ((dev**2) ** 2).sum(axis=-1) / count
            return np.where(zero, np.nan, m4 / m2**2.0) - 3
    return np.full(values.shape[:-1], np.nan)
