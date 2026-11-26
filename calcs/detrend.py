# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2024-2026 Mark Jessell
"""
The regional plane removed from a grid before Fourier filtering.

numpy only, so both the FFT filters (GeophysicalProcessor) and the variable
continuation (TaylorContinuation) can use it.
"""

import numpy as np

MAX_FIT_CELLS = 250_000  # larger grids are fitted on a regular subsample


def robust_plane(filled, iterations=200, tolerance=1e-7):
    """
    Plane a + b*row + c*col through a gap-filled grid (no NaN), fitted over
    the whole rectangle by least absolute deviations (iteratively reweighted
    least squares).

    The fit is over the whole rectangle, not just the cells with data: a plane
    fitted to an irregular survey alone and extrapolated to the grid's corners
    leaves large values at the edges, which the padding turns into broad false
    anomalies. It is robust (L1) so that it follows the background rather
    than being pulled up or down by a few strong anomalies, as a least-squares
    plane is: the edges of the residual then sit near zero, which is what the
    padding expects.

    The iterations go on until the plane at the grid's edges moves by less than
    tolerance times the grid's standard deviation (or up to `iterations`): run
    to convergence, the plane no longer depends on rounding, so it is the same
    with any numpy version.

    Returns (plane grid, (a, b, c)).
    """
    ny, nx = filled.shape
    step = max(1, int(np.ceil(np.sqrt(filled.size / MAX_FIT_CELLS))))
    rows, cols = np.mgrid[0:ny:step, 0:nx:step]
    y = filled[::step, ::step].ravel().astype(float)
    # rows and columns centred on the grid: a well-conditioned 3 x 3 system
    r0, c0 = (ny - 1) / 2, (nx - 1) / 2
    r, c = rows.ravel() - r0, cols.ravel() - c0
    A = np.column_stack([np.ones(y.size), r, c])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    scale = np.std(y) + 1e-30
    floor = 1e-6 * scale
    rr, rc, cc, ry, cy = r * r, r * c, c * c, r * y, c * y
    for _ in range(int(iterations)):
        # weighted least squares with weights 1/|residual|, by its normal equations
        w = 1.0 / np.maximum(np.abs(y - coef[0] - coef[1] * r - coef[2] * c), floor)
        sw, swr, swc, swrc = w.sum(), w @ r, w @ c, w @ rc
        M = np.array([[sw, swr, swc], [swr, w @ rr, swrc], [swc, swrc, w @ cc]])
        previous, coef = coef, np.linalg.solve(M, [w @ y, w @ ry, w @ cy])
        if np.max(np.abs(coef - previous) * [1, ny, nx]) <= tolerance * scale:
            break
    a, b, c = coef[0] - coef[1] * r0 - coef[2] * c0, coef[1], coef[2]
    rows, cols = np.indices(filled.shape)
    return a + b * rows + c * cols, (a, b, c)
