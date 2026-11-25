"""Slope of a grid, as QGIS's Slope algorithm (native:slope) gives it (numpy only)."""

import numpy as np


def slope_degrees(values, dx, dy, z_factor=1.0):
    """
    Slope in degrees by Horn's 3 x 3 method. values: 2-D array, row 0 = north;
    dx, dy: cell sizes (positive). As in QGIS, the outer row and column of cells
    and every cell next to a gap (NaN) are NaN. Agrees with QGIS to its float32
    rounding (about 1e-6 degrees).
    """
    z = np.asarray(values, dtype=float) * z_factor
    p = np.pad(z, 1, constant_values=np.nan)
    ny, nx = z.shape

    def w(r, c):  # the neighbour r rows down and c columns right of each cell
        return p[1 + r:1 + r + ny, 1 + c:1 + c + nx]

    gx = ((w(-1, 1) - w(-1, -1)) + 2 * (w(0, 1) - w(0, -1)) + (w(1, 1) - w(1, -1))) / (8 * dx)
    gy = ((w(1, -1) - w(-1, -1)) + 2 * (w(1, 0) - w(-1, 0)) + (w(1, 1) - w(-1, 1))) / (8 * dy)
    with np.errstate(invalid="ignore"):
        slope = np.degrees(np.arctan(np.sqrt(gx * gx + gy * gy)))
    # Horn's method leaves out the cell itself: as in QGIS, a single empty cell
    # with data all round gets a slope (its neighbours, next to it, get none)
    slope[0, :] = slope[-1, :] = np.nan
    slope[:, 0] = slope[:, -1] = np.nan
    return slope
