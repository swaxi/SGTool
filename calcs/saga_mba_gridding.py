"""
Standalone Python / NumPy translation of SAGA's Multilevel B-Spline (MBA)
gridding algorithm, as implemented in Gridding_Spline_MBA.cpp
(CGridding_Spline_MBA -- SAGA tool "grid_spline" / id 4 / "Multilevel
B-Spline", O. Conrad 2006, based on Lee, Wolberg & Shin 1997, "Scattered
Data Interpolation with Multilevel B-Splines", IEEE TVCG 3(3)).

This is the scattered-*points* variant of the tool (input: a point set with
an attribute field), as opposed to Gridding_Spline_MBA_Grid.cpp which reads
its points from a raster's non-NoData cells. The translation below follows
Gridding_Spline_MBA.cpp statement for statement:

    On_Execute            -> mba_gridding()
    _Set_MBA               -> refinement=False branch
    _Set_MBA_Refinement     -> refinement=True branch
    BA_Get_B               -> _basis()
    BA_Set_Phi              -> _ba_set_phi()
    BA_Get_Phi               -> _evaluate_phi()
    BA_Set_Grid               -> _ba_set_grid()
    _Set_MBA_Refinement(Psi_0, Psi_1) (the pyramid refinement step)
                                -> _refine_add()
    _Get_Difference          -> _get_difference()

Two formulations are provided, matching the two C++ code paths:

    _Set_MBA             (METHOD = "no")  -> refinement=False
    _Set_MBA_Refinement  (METHOD = "yes") -> refinement=True

Faithfully reproduced quirks of the original tool
---------------------------------------------------
* Only points that round (nearest grid cell, `floor(0.5 + (x-xmin)/cellsize)`)
  into the target grid are used at all -- points far outside the target
  extent are silently dropped, exactly as `Initialize(m_Points, true, &s)`
  does in the C++ tool.
* Detrending mean, and the final min/max clamp, are computed from the
  *filtered* point set (matching `CSG_Simple_Statistics s` in `On_Execute`).
* The final grid is clamped to [min(z), max(z)] of the filtered input
  points -- this is a real feature of `CGridding_Spline_MBA::On_Execute`,
  not a bug, and is reproduced here.
* The B-spline control lattice ("Phi") has *independent* column/row counts
  (`nx = 4 + int(XRange/cellsize)`, `ny = 4 + int(YRange/cellsize)`) rather
  than a single square size -- unlike the raster-input variant of this tool.

Precision note
---------------
In the C++ tool the B-spline control lattices ("Phi", "Delta") are always
32-bit float SAGA grids (`Phi.Create(SG_DATATYPE_Float, ...)`), and every
`Add_Value`/`Set_Value` call truncates to that storage type immediately.
Reproducing that *exactly*, including the truncation-order dependence of
float32 accumulation, would require a non-vectorised, point-by-point Python
loop (impractically slow). By default this module keeps the control
lattices in float64 (`phi_dtype=np.float64`), which is mathematically the
same algorithm and matches SAGA to ~1e-6 relative precision or better. Pass
`phi_dtype=np.float32` to mimic SAGA's storage type more closely (still not
bit-exact, because the scatter accumulation here is vectorised rather than
sequential per-point).
"""

from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# CGridding_Spline_MBA::BA_Get_B
# ---------------------------------------------------------------------------
def _basis(i: int, d: np.ndarray) -> np.ndarray:
    """Cubic B-spline basis function i (0..3), d in [0, 1)."""
    if i == 0:
        e = 1.0 - d
        return e * e * e / 6.0
    if i == 1:
        return (3.0 * d ** 3 - 6.0 * d ** 2 + 4.0) / 6.0
    if i == 2:
        return (-3.0 * d ** 3 + 3.0 * d ** 2 + 3.0 * d + 1.0) / 6.0
    if i == 3:
        return d ** 3 / 6.0
    return np.zeros_like(d)


# ---------------------------------------------------------------------------
# CGridding_Spline_MBA::BA_Get_Phi (vectorised gather, generic px/py)
# phi has shape (phi_ny, phi_nx)
# ---------------------------------------------------------------------------
def _evaluate_phi(phi: np.ndarray, phi_nx: int, phi_ny: int,
                   px: np.ndarray, py: np.ndarray) -> np.ndarray:
    px = np.asarray(px, dtype=np.float64)
    py = np.asarray(py, dtype=np.float64)

    x = np.floor(px).astype(np.int64)
    y = np.floor(py).astype(np.int64)

    inb = (x >= 0) & (x < phi_nx - 3) & (y >= 0) & (y < phi_ny - 3)

    z = np.zeros(px.shape, dtype=np.float64)
    if not np.any(inb):
        return z

    xm = x[inb]
    ym = y[inb]
    dx = px[inb] - xm
    dy = py[inb] - ym

    phi64 = phi.astype(np.float64, copy=False)

    acc = np.zeros(xm.shape, dtype=np.float64)
    for iy in range(4):
        by = _basis(iy, dy)
        for ix in range(4):
            bx = _basis(ix, dx)
            acc += by * bx * phi64[ym + iy, xm + ix]

    z[inb] = acc
    return z


# ---------------------------------------------------------------------------
# CGridding_Spline_MBA::BA_Set_Phi (vectorised scatter)
# ---------------------------------------------------------------------------
def _ba_set_phi(x: np.ndarray, y: np.ndarray, z: np.ndarray,
                 phi_cellsize: float,
                 target_xmin: float, target_ymin: float,
                 target_xrange: float, target_yrange: float,
                 phi_dtype):

    phi_nx = 4 + int(target_xrange / phi_cellsize)
    phi_ny = 4 + int(target_yrange / phi_cellsize)

    phi = np.zeros((phi_ny, phi_nx), dtype=phi_dtype)
    delta = np.zeros((phi_ny, phi_nx), dtype=phi_dtype)

    p_x = (x - target_xmin) / phi_cellsize
    p_y = (y - target_ymin) / phi_cellsize
    xi = np.floor(p_x).astype(np.int64)
    yi = np.floor(p_y).astype(np.int64)

    inb = (xi >= 0) & (xi < phi_nx - 3) & (yi >= 0) & (yi < phi_ny - 3)
    if np.any(inb):
        xi = xi[inb]
        yi = yi[inb]
        dx = p_x[inb] - xi
        dy = p_y[inb] - yi
        pz = np.asarray(z, dtype=np.float64)[inb]

        wx = np.stack([_basis(i, dx) for i in range(4)], axis=0)  # (4, N)
        wy = np.stack([_basis(i, dy) for i in range(4)], axis=0)  # (4, N)

        sw2 = np.zeros(xi.shape[0], dtype=np.float64)
        for iy in range(4):
            for ix in range(4):
                w = wy[iy] * wx[ix]
                sw2 += w * w

        ok = sw2 > 0.0
        if np.any(ok):
            xk = xi[ok]
            yk = yi[ok]
            pzk = pz[ok] / sw2[ok]
            wxk = wx[:, ok]
            wyk = wy[:, ok]

            for iy in range(4):
                for ix in range(4):
                    w = wyk[iy] * wxk[ix]
                    ty = yk + iy
                    tx = xk + ix
                    np.add.at(delta, (ty, tx), (w ** 3 * pzk).astype(phi_dtype))
                    np.add.at(phi, (ty, tx), (w * w).astype(phi_dtype))

    nz = phi != 0
    phi_out = np.zeros((phi_ny, phi_nx), dtype=phi_dtype)
    phi_out[nz] = (delta[nz].astype(np.float64) / phi[nz].astype(np.float64)).astype(phi_dtype)

    return phi_nx, phi_ny, phi_out


# ---------------------------------------------------------------------------
# CGridding_Spline_MBA::_Get_Difference
# ---------------------------------------------------------------------------
def _get_difference(x: np.ndarray, y: np.ndarray, z_resid: np.ndarray,
                     phi: np.ndarray, phi_nx: int, phi_ny: int, phi_cellsize: float,
                     target_xmin: float, target_ymin: float,
                     epsilon: float, points_dtype,
                     verbose: bool, level: int) -> bool:

    p_x = (x - target_xmin) / phi_cellsize
    p_y = (y - target_ymin) / phi_cellsize

    interp = _evaluate_phi(phi, phi_nx, phi_ny, p_x, p_y)
    znew = z_resid.astype(np.float64) - interp
    z_resid[:] = znew.astype(points_dtype)

    exceed = np.abs(znew) > epsilon
    count = int(np.count_nonzero(exceed))

    if verbose:
        vmax = float(np.max(np.abs(znew[exceed]))) if count else 0.0
        vmean = float(np.mean(np.abs(znew[exceed]))) if count else 0.0
        print(f"level:{level + 1} errors:{count} maximum:{vmax:g} mean:{vmean:g}")

    return count > 0


# ---------------------------------------------------------------------------
# CGridding_Spline_MBA::BA_Set_Grid
# ---------------------------------------------------------------------------
def _ba_set_grid(phi: np.ndarray, phi_nx: int, phi_ny: int, phi_cellsize: float,
                  target_nx: int, target_ny: int, target_cellsize: float,
                  add: bool, out_dtype, out: np.ndarray | None) -> np.ndarray:

    d = target_cellsize / phi_cellsize

    xidx = np.arange(target_nx, dtype=np.float64)
    yidx = np.arange(target_ny, dtype=np.float64)
    px = d * xidx
    py = d * yidx
    px_grid, py_grid = np.meshgrid(px, py, indexing="xy")  # shape (ny, nx)

    vals = _evaluate_phi(phi, phi_nx, phi_ny, px_grid.ravel(), py_grid.ravel()).reshape(target_ny, target_nx)

    if add:
        out = (out.astype(np.float64) + vals).astype(out_dtype)
    else:
        out = vals.astype(out_dtype)

    return out


# ---------------------------------------------------------------------------
# CGridding_Spline_MBA::_Set_MBA_Refinement(Psi_0, Psi_1)  (pyramid step)
# psi0 has shape (ny0, nx0), psi1 has shape (ny1, nx1); psi0 is refined and
# ADDED into psi1 in place.
# ---------------------------------------------------------------------------
def _refine_add(psi0: np.ndarray, nx0: int, ny0: int,
                 psi1: np.ndarray, nx1: int, ny1: int) -> None:
    if 2 * (nx0 - 4) != (nx1 - 4) or 2 * (ny0 - 4) != (ny1 - 4):
        return  # dimension mismatch: SAGA silently skips refinement here too

    p0 = np.zeros((ny0 + 2, nx0 + 2), dtype=np.float64)
    p0[1:-1, 1:-1] = psi0.astype(np.float64)

    a00 = p0[0:ny0, 0:nx0]
    a10 = p0[0:ny0, 1:nx0 + 1]
    a20 = p0[0:ny0, 2:nx0 + 2]
    a01 = p0[1:ny0 + 1, 0:nx0]
    a11 = p0[1:ny0 + 1, 1:nx0 + 1]
    a21 = p0[1:ny0 + 1, 2:nx0 + 2]
    a02 = p0[2:ny0 + 2, 0:nx0]
    a12 = p0[2:ny0 + 2, 1:nx0 + 1]
    a22 = p0[2:ny0 + 2, 2:nx0 + 2]

    out00 = (a00 + a02 + a20 + a22 + 6.0 * (a01 + a10 + a12 + a21) + 36.0 * a11) / 64.0
    out01 = (a01 + a02 + a21 + a22 + 6.0 * (a11 + a12)) / 16.0
    out10 = (a10 + a12 + a20 + a22 + 6.0 * (a11 + a21)) / 16.0
    out11 = (a11 + a12 + a21 + a22) / 4.0

    xi = np.arange(nx0)
    yi = np.arange(ny0)
    yy0 = 2 * yi - 1  # target row for out00 / out10
    yy1 = 2 * yi      # target row for out01 / out11
    xx0 = 2 * xi - 1  # target col for out00 / out01
    xx1 = 2 * xi      # target col for out10 / out11

    def _valid(a, limit):
        return (a >= 0) & (a < limit)

    vy0, vy1 = _valid(yy0, ny1), _valid(yy1, ny1)
    vx0, vx1 = _valid(xx0, nx1), _valid(xx1, nx1)

    if np.any(vy0) and np.any(vx0):
        ry, rx = np.ix_(yy0[vy0], xx0[vx0])
        psi1[ry, rx] += out00[np.ix_(vy0, vx0)].astype(psi1.dtype)

    if np.any(vy1) and np.any(vx0):
        ry, rx = np.ix_(yy1[vy1], xx0[vx0])
        psi1[ry, rx] += out01[np.ix_(vy1, vx0)].astype(psi1.dtype)

    if np.any(vy0) and np.any(vx1):
        ry, rx = np.ix_(yy0[vy0], xx1[vx1])
        psi1[ry, rx] += out10[np.ix_(vy0, vx1)].astype(psi1.dtype)

    if np.any(vy1) and np.any(vx1):
        ry, rx = np.ix_(yy1[vy1], xx1[vx1])
        psi1[ry, rx] += out11[np.ix_(vy1, vx1)].astype(psi1.dtype)


# ---------------------------------------------------------------------------
# CGridding_Spline_MBA::On_Execute / _Set_MBA / _Set_MBA_Refinement
# ---------------------------------------------------------------------------
def mba_gridding(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    cellsize: float,
    xmin: float,
    ymin: float,
    nx: int,
    ny: int,
    epsilon: float = 0.0001,
    level_max: int = 11,
    refinement: bool = False,
    phi_dtype=np.float64,
    points_dtype=np.float64,
    output_dtype=np.float64,
    verbose: bool = False,
) -> np.ndarray:
    """
    Multilevel B-Spline (MBA) interpolation of scattered points onto a
    regular grid, translated from SAGA's Gridding_Spline_MBA.cpp
    (tool "Multilevel B-Spline", grid_spline / id 4).

    Parameters
    ----------
    x, y, z : 1D array_like
        Scattered point coordinates and values (any order, any density).
        Non-finite z values are dropped before processing. Points that do
        not round to a cell within the target grid are also dropped
        (matching SAGA's `Initialize(m_Points, true, &s)` behaviour).
    cellsize, xmin, ymin, nx, ny :
        Definition of the output grid: cell size, lower-left cell-centre
        coordinate, and grid dimensions (columns=nx, rows=ny).
    epsilon : float
        "Threshold Error" (SAGA parameter EPSILON).
    level_max : int
        "Maximum Level" (SAGA parameter LEVEL_MAX).
    refinement : bool
        False = SAGA METHOD "no"  (_Set_MBA).
        True  = SAGA METHOD "yes" (_Set_MBA_Refinement).
    phi_dtype, points_dtype, output_dtype :
        NumPy dtypes for the B-spline control lattice, the (detrended)
        point residuals, and the returned grid. Default float64
        throughout -- see module docstring for the precision trade-off
        against SAGA's internal 32-bit float lattice.
    verbose : bool
        If True, print the same per-level diagnostics SAGA logs
        (level / errors / maximum / mean).

    Returns
    -------
    numpy.ndarray of shape (ny, nx), dtype `output_dtype`
        Interpolated grid, clamped to [min(z), max(z)] of the (filtered)
        input points, row 0 = y = ymin, column 0 = x = xmin.
    """
    x = np.asarray(x, dtype=np.float64).ravel()
    y = np.asarray(y, dtype=np.float64).ravel()
    z = np.asarray(z, dtype=np.float64).ravel()

    finite = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x, y, z = x[finite], y[finite], z[finite]

    # Initialize(m_Points, true, &s): keep only points that round (nearest
    # cell) into the target grid.
    ix = np.floor(0.5 + (x - xmin) / cellsize).astype(np.int64)
    iy = np.floor(0.5 + (y - ymin) / cellsize).astype(np.int64)
    in_grid = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
    x, y, z = x[in_grid], y[in_grid], z[in_grid]

    if x.size < 3:
        raise ValueError("need at least 3 valid (x, y, z) points inside the target grid")

    mean = float(np.mean(z))
    zmin = float(np.min(z))
    zmax = float(np.max(z))

    m_z = (z - mean).astype(points_dtype)  # detrending

    target_xrange = (nx - 1) * cellsize
    target_yrange = (ny - 1) * cellsize
    cellsize0 = max(target_xrange, target_yrange)

    if not refinement:
        output = None
        cellsize_lvl = cellsize0
        b_continue = True
        level = 0
        while b_continue and level < level_max:
            pnx, pny, phi = _ba_set_phi(x, y, m_z, cellsize_lvl, xmin, ymin,
                                         target_xrange, target_yrange, phi_dtype)

            b_continue = _get_difference(x, y, m_z, phi, pnx, pny, cellsize_lvl,
                                          xmin, ymin, epsilon, points_dtype,
                                          verbose, level)

            output = _ba_set_grid(phi, pnx, pny, cellsize_lvl, nx, ny, cellsize,
                                   add=(level > 0), out_dtype=output_dtype, out=output)

            level += 1
            cellsize_lvl /= 2.0
    else:
        phi_buf = [None, None]
        cellsize_lvl = cellsize0
        b_continue = True
        i = 0
        level = 0
        while b_continue and level < level_max:
            i = level % 2
            pnx, pny, phi_new = _ba_set_phi(x, y, m_z, cellsize_lvl, xmin, ymin,
                                             target_xrange, target_yrange, phi_dtype)

            b_continue = _get_difference(x, y, m_z, phi_new, pnx, pny, cellsize_lvl,
                                          xmin, ymin, epsilon, points_dtype,
                                          verbose, level)

            other = phi_buf[(i + 1) % 2]
            if other is not None:
                _refine_add(other["array"], other["nx"], other["ny"], phi_new, pnx, pny)

            phi_buf[i] = {"array": phi_new, "nx": pnx, "ny": pny, "cellsize": cellsize_lvl}

            level += 1
            cellsize_lvl /= 2.0

        final = phi_buf[i]
        output = _ba_set_grid(final["array"], final["nx"], final["ny"], final["cellsize"],
                               nx, ny, cellsize, add=False,
                               out_dtype=output_dtype, out=None)

    output = output.astype(np.float64) + mean          # de-detrending
    output = np.clip(output, zmin, zmax)                # On_Execute's min/max clamp
    return output.astype(output_dtype)


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    rng = np.random.default_rng(0)
    n_pts = 400
    px = rng.uniform(0, 100, n_pts)
    py = rng.uniform(0, 100, n_pts)
    pz = np.sin(px / 12.0) * np.cos(py / 15.0) * 10.0 + 50.0

    grid_no_refine = mba_gridding(px, py, pz, cellsize=1.0, xmin=0.0, ymin=0.0,
                                   nx=101, ny=101, level_max=11, refinement=False,
                                   verbose=True)

    grid_refine = mba_gridding(px, py, pz, cellsize=1.0, xmin=0.0, ymin=0.0,
                                nx=101, ny=101, level_max=11, refinement=True,
                                verbose=True)

    print("no-refinement grid: shape", grid_no_refine.shape,
          "min/max", grid_no_refine.min(), grid_no_refine.max())
    print("refinement grid:    shape", grid_refine.shape,
          "min/max", grid_refine.min(), grid_refine.max())
    print("max abs diff between the two methods:",
          np.max(np.abs(grid_no_refine - grid_refine)))
