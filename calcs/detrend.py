"""
The regional plane removed from a grid before Fourier filtering.

numpy only, so both the FFT filters (GeophysicalProcessor) and the variable
continuation (TaylorContinuation) can use it.
"""

import numpy as np

MAX_FIT_CELLS = 250_000  # larger grids are fitted on a regular subsample


def robust_plane(filled, iterations=20):
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

    Returns (plane grid, (a, b, c)).
    """
    ny, nx = filled.shape
    step = max(1, int(np.ceil(np.sqrt(filled.size / MAX_FIT_CELLS))))
    rows, cols = np.mgrid[0:ny:step, 0:nx:step]
    y = filled[::step, ::step].ravel().astype(float)
    A = np.column_stack([np.ones(y.size), rows.ravel(), cols.ravel()])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    floor = 1e-6 * (np.std(y) + 1e-30)
    for _ in range(int(iterations)):
        sw = 1.0 / np.sqrt(np.maximum(np.abs(y - A @ coef), floor))
        coef, *_ = np.linalg.lstsq(A * sw[:, None], y * sw, rcond=None)
    r, c = np.indices(filled.shape)
    return coef[0] + coef[1] * r + coef[2] * c, tuple(coef)
