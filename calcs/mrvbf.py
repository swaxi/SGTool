"""
Standalone Python / NumPy translation of SAGA's Multiresolution Index of
Valley Bottom Flatness tool (mrvbf.cpp, CMRVBF; ta_morphometry / tool id 8;
O. Conrad 2006, based on Gallant & Dowling 2003, "A multiresolution index
of valley bottom flatness for mapping depositional areas", Water Resources
Research 39/12).

Input is a DEM as a 2D NumPy array (row 0 = y = ymin, "row increases
northward", matching SAGA's raster convention).

Translated pieces (C++ -> Python):

    On_Execute        -> mrvbf()
    Get_Smoothed       -> _smoothed()               (11x11 weighted mean)
    Get_Slopes          -> _slope_percent()          (Zevenbergen & Thorne)
    Get_Percentiles      -> _percentile()             (disk-neighbourhood rank)
    Get_Values             -> inlined in the level>=3 loop of mrvbf()
                              (uses _smoothed / _slope_percent / _nearest_resample
                              / _percentile exactly as CMRVBF::Get_Values does)
    Get_Flatness              -> _flatness()
    Get_MRVBF                  -> _blend()
    Get_Classified               -> _classify()
    Get_Transformation             -> _transform()
    CSG_Grid::Get_Value(...,Bicubic_2) -> _sample_bspline()  (this is the
        resampling method used by Get_Flatness -- SAGA's default Get_Value()
        interpolation is the *uniform cubic B-spline* ("BSpline"), not
        bilinear or an interpolating bicubic; it is reproduced here
        including its exact behaviour at grid/NoData boundaries: a point is
        only sampled if its nearest source cell is valid, and any NoData
        cell inside the 4x4 support window is filled by the same iterative
        "average of surrounding data" diffusion SAGA's
        _Get_ValAtPos_Fill4x4Submatrix() performs (up to 16 passes, looking
        one ring of true source data beyond the 4x4 window itself).

Everything here mirrors the C++ control flow and arithmetic statement for
statement; the vectorisation (NumPy arrays instead of per-pixel loops) does
not change any numeric result.

Precision note
---------------
All arithmetic is carried out in float64; SAGA stores its intermediate
grids (DEM copies, Slopes, Percentiles, Smoothed, CF/VF/RF) as 32-bit float
SAGA grids, truncating on every write. As with the companion MBA
translation in ../grid_spline/mba_gridding.py, this is not reproduced
bit-for-bit (would require per-pixel truncation matching the exact write
order); expect agreement to roughly float32 relative precision (~1e-6) when
compared against the compiled SAGA tool.
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import convolve


# ---------------------------------------------------------------------------
# generic helpers
# ---------------------------------------------------------------------------
def _shift(arr: np.ndarray, valid: np.ndarray, dy: int, dx: int):
    """result[y, x] = (arr[y+dy, x+dx], valid[y+dy, x+dx]) if in bounds,
    else (0, False). Mirrors CSG_Grid::is_InGrid(x+ix, y+iy) with the
    default bCheckNoData=True."""
    ny, nx = arr.shape
    out_val = np.zeros_like(arr)
    out_ok = np.zeros(valid.shape, dtype=bool)

    y0, y1 = max(0, -dy), ny - max(0, dy)
    x0, x1 = max(0, -dx), nx - max(0, dx)
    if y0 < y1 and x0 < x1:
        ys0, ys1 = y0 + dy, y1 + dy
        xs0, xs1 = x0 + dx, x1 + dx
        out_val[y0:y1, x0:x1] = arr[ys0:ys1, xs0:xs1]
        out_ok[y0:y1, x0:x1] = valid[ys0:ys1, xs0:xs1]
    return out_val, out_ok


def _bspline_basis(i: int, d: np.ndarray) -> np.ndarray:
    """Uniform cubic B-spline basis function i (0..3), d in [0, 1). This is
    the same blending family SAGA's _Get_ValAtPos_BiCubic_2 computes via a
    different (but numerically equivalent) summation formula."""
    if i == 0:
        e = 1.0 - d
        return e * e * e / 6.0
    if i == 1:
        return (3.0 * d ** 3 - 6.0 * d ** 2 + 4.0) / 6.0
    if i == 2:
        return (-3.0 * d ** 3 + 3.0 * d ** 2 + 3.0 * d + 1.0) / 6.0
    return d ** 3 / 6.0


def _transform(x: np.ndarray, t: float, p: float) -> np.ndarray:
    """CMRVBF::Get_Transformation: 1 / (1 + (x/t)^p)."""
    return 1.0 / (1.0 + np.power(x / t, p))


# ---------------------------------------------------------------------------
# CMRVBF::Get_Smoothed
# ---------------------------------------------------------------------------
def _smoothed(dem: np.ndarray, valid: np.ndarray, radius: int = 5):
    iy, ix = np.mgrid[-radius:radius + 1, -radius:radius + 1]
    d = np.sqrt(ix.astype(np.float64) ** 2 + iy.astype(np.float64) ** 2)
    kernel = 4.3565 * np.exp(-(d / 3.0) ** 2)

    data0 = np.where(valid, dem, 0.0)
    w = valid.astype(np.float64)

    num = convolve(data0, kernel, mode="constant", cval=0.0)
    den = convolve(w, kernel, mode="constant", cval=0.0)

    out_valid = den > 0.0
    out = np.divide(num, den, out=np.zeros_like(num), where=out_valid)
    return out, out_valid


# ---------------------------------------------------------------------------
# CMRVBF::Get_Slopes (CSG_Grid::Get_Gradient, Zevenbergen & Thorne 1986)
# ---------------------------------------------------------------------------
def _slope_percent(dem: np.ndarray, valid: np.ndarray, cellsize: float):
    z = dem

    N_val, N_ok = _shift(dem, valid, 1, 0)
    S_val, S_ok = _shift(dem, valid, -1, 0)
    E_val, E_ok = _shift(dem, valid, 0, 1)
    W_val, W_ok = _shift(dem, valid, 0, -1)

    dz0 = np.where(N_ok, N_val - z, np.where(S_ok, z - S_val, 0.0))
    dz2 = np.where(S_ok, S_val - z, np.where(N_ok, z - N_val, 0.0))
    dz1 = np.where(E_ok, E_val - z, np.where(W_ok, z - W_val, 0.0))
    dz3 = np.where(W_ok, W_val - z, np.where(E_ok, z - E_val, 0.0))

    G = (dz0 - dz2) / (2.0 * cellsize)
    H = (dz1 - dz3) / (2.0 * cellsize)

    slope = np.arctan(np.sqrt(G * G + H * H))
    slope_pct = 100.0 * np.tan(slope)

    return slope_pct, valid.copy()


# ---------------------------------------------------------------------------
# CMRVBF::Get_Percentile / Get_Percentiles (CSG_Grid_Radius disk neighbourhood)
# ---------------------------------------------------------------------------
def _percentile(dem: np.ndarray, valid: np.ndarray, radius: int):
    offsets = [(dx, dy) for dy in range(-radius, radius + 1)
               for dx in range(-radius, radius + 1)
               if dx * dx + dy * dy < radius * radius]

    n_points = np.zeros(dem.shape, dtype=np.int32)
    n_lower = np.zeros(dem.shape, dtype=np.int32)

    for dx, dy in offsets:
        v, ok = _shift(dem, valid, dy, dx)
        n_points += ok
        n_lower += ok & (v < dem)

    good = valid & (n_points > 1)
    pct = np.zeros(dem.shape, dtype=np.float64)
    pct[good] = n_lower[good] / (n_points[good] - 1.0)
    return pct, good


# ---------------------------------------------------------------------------
# CSG_Grid::Get_Value(x, y, ..., Bicubic_2)  ("BSpline" resampling), used by
# Get_Flatness to sample Slopes/Percentiles at the master grid's node
# positions. Includes the exact NoData behaviour of
# _Get_ValAtPos_Fill4x4Submatrix.
# ---------------------------------------------------------------------------
def _sample_bspline(src: np.ndarray, src_valid: np.ndarray,
                     cellsize: float, xmin: float, ymin: float,
                     px: np.ndarray, py: np.ndarray):
    ny, nx = src.shape
    px = np.asarray(px, dtype=np.float64)
    py = np.asarray(py, dtype=np.float64)
    m_total = px.size

    out_val = np.zeros(m_total, dtype=np.float64)
    out_ok = np.zeros(m_total, dtype=bool)

    # m_System.Get_Extent(true).Contains(x, y)  (cell-based extent, +/-0.5 cell)
    xlo = xmin - 0.5 * cellsize
    xhi = xmin + (nx - 1) * cellsize + 0.5 * cellsize
    ylo = ymin - 0.5 * cellsize
    yhi = ymin + (ny - 1) * cellsize + 0.5 * cellsize
    in_extent = (px >= xlo) & (px <= xhi) & (py >= ylo) & (py <= yhi)

    fx = (px - xmin) / cellsize
    fy = (py - ymin) / cellsize
    ix = np.floor(fx).astype(np.int64)
    iy = np.floor(fy).astype(np.int64)
    dx = fx - ix
    dy = fy - iy

    # is_InGrid(ix + round(dx), iy + round(dy)) gate
    rx = ix + np.floor(0.5 + dx).astype(np.int64)
    ry = iy + np.floor(0.5 + dy).astype(np.int64)
    nb = (rx >= 0) & (rx < nx) & (ry >= 0) & (ry < ny)
    nearest_ok = np.zeros(m_total, dtype=bool)
    if np.any(nb):
        rxc, ryc = rx[nb], ry[nb]
        nearest_ok[nb] = src_valid[ryc, rxc]

    candidate = in_extent & nearest_ok
    idx = np.nonzero(candidate)[0]
    if idx.size == 0:
        return out_val, out_ok

    ixc, iyc = ix[idx], iy[idx]
    dxc, dyc = dx[idx], dy[idx]
    m = idx.size

    # raw 6x6 window: local r (0..5) <-> global offset r-2  (covers one extra
    # ring beyond the 4x4 support window needed by the gap-fill diffusion)
    raw_val = np.zeros((6, 6, m), dtype=np.float64)
    raw_ok = np.zeros((6, 6, m), dtype=bool)
    for rj in range(6):
        gy = iyc - 2 + rj
        gy_in = (gy >= 0) & (gy < ny)
        gy_c = np.clip(gy, 0, ny - 1)
        for rc in range(6):
            gx = ixc - 2 + rc
            gx_in = (gx >= 0) & (gx < nx)
            m_in = gy_in & gx_in
            if not np.any(m_in):
                continue
            gx_c = np.clip(gx, 0, nx - 1)
            vv = src[gy_c, gx_c]
            ok = m_in & src_valid[gy_c, gx_c]
            raw_val[rj, rc, m_in] = vv[m_in]
            raw_ok[rj, rc] = ok

    # inner 4x4 support window == raw6[1:5, 1:5]
    v = raw_val[1:5, 1:5, :].copy()
    ok4 = raw_ok[1:5, 1:5, :].copy()

    n_missing = (~ok4).sum(axis=(0, 1))
    need_fill = (n_missing > 0) & (n_missing < 16)

    for _ in range(16):
        if not np.any(need_fill):
            break
        v_snap = v.copy()
        ok_snap = ok4.copy()

        for iy4 in range(4):
            for ix4 in range(4):
                missing_here = (~ok_snap[iy4, ix4]) & need_fill
                if not np.any(missing_here):
                    continue
                s = np.zeros(m, dtype=np.float64)
                n = np.zeros(m, dtype=np.int64)
                for jy in range(iy4 - 1, iy4 + 2):
                    ry6 = jy + 1
                    for jx in range(ix4 - 1, ix4 + 2):
                        rx6 = jx + 1
                        rawv = raw_val[ry6, rx6]
                        rawok = raw_ok[ry6, rx6]
                        in_inner = 0 <= jx < 4 and 0 <= jy < 4
                        if in_inner:
                            snap_ok = ok_snap[jy, jx] & (~rawok)
                            val = np.where(rawok, rawv, v_snap[jy, jx])
                            take = missing_here & (rawok | snap_ok)
                        else:
                            val = rawv
                            take = missing_here & rawok
                        s = np.where(take, s + val, s)
                        n = np.where(take, n + 1, n)

                fillable = missing_here & (n > 0)
                if np.any(fillable):
                    v[iy4, ix4] = np.where(fillable, s / np.maximum(n, 1), v[iy4, ix4])
                    ok4[iy4, ix4] = ok4[iy4, ix4] | fillable

        n_missing = (~ok4).sum(axis=(0, 1))
        need_fill = (n_missing > 0) & (n_missing < 16)

    fully_ok = ok4.all(axis=(0, 1))

    Rx = np.stack([_bspline_basis(i, dxc) for i in range(4)], axis=0)  # (4, m)
    Ry = np.stack([_bspline_basis(i, dyc) for i in range(4)], axis=0)  # (4, m)

    val = np.zeros(m, dtype=np.float64)
    for iy4 in range(4):
        for ix4 in range(4):
            val += v[iy4, ix4] * Rx[ix4] * Ry[iy4]

    sel = idx[fully_ok]
    out_val[sel] = val[fully_ok]
    out_ok[sel] = True
    return out_val, out_ok


# ---------------------------------------------------------------------------
# CSG_Grid::_Assign_Interpolated with NearestNeighbour resampling, used by
# CMRVBF::Get_Values to down-sample the smoothed DEM onto the next
# (coarser) working resolution.
# ---------------------------------------------------------------------------
def _nearest_resample(src: np.ndarray, src_valid: np.ndarray, src_cs: float,
                       src_xmin: float, src_ymin: float,
                       dst_nx: int, dst_ny: int, dst_cs: float,
                       dst_xmin: float, dst_ymin: float):
    ny_s, nx_s = src.shape

    xg = dst_xmin + np.arange(dst_nx) * dst_cs
    yg = dst_ymin + np.arange(dst_ny) * dst_cs
    fx = (xg - src_xmin) / src_cs
    fy = (yg - src_ymin) / src_cs
    ixg = np.floor(fx)
    dxg = fx - ixg
    iyg = np.floor(fy)
    dyg = fy - iyg

    nxg = (ixg + np.floor(0.5 + dxg)).astype(np.int64)
    nyg = (iyg + np.floor(0.5 + dyg)).astype(np.int64)

    NX, NY = np.meshgrid(nxg, nyg)  # shape (dst_ny, dst_nx)
    inb = (NX >= 0) & (NX < nx_s) & (NY >= 0) & (NY < ny_s)

    NXc = np.clip(NX, 0, nx_s - 1)
    NYc = np.clip(NY, 0, ny_s - 1)

    vals = src[NYc, NXc]
    ok = inb & src_valid[NYc, NXc]

    out = np.where(ok, vals, 0.0)
    return out, ok


# ---------------------------------------------------------------------------
# CMRVBF::Get_Flatness
# ---------------------------------------------------------------------------
def _flatness(slopes, slopes_valid, s_cs, s_xmin, s_ymin,
              pctl, pctl_valid, p_cs, p_xmin, p_ymin,
              CF, XPf, YPf, out_shape, t_slope, p_slope, t_pctl_v, t_pctl_r, p_pctl):

    sv, sok = _sample_bspline(slopes, slopes_valid, s_cs, s_xmin, s_ymin, XPf, YPf)
    pv, pok = _sample_bspline(pctl, pctl_valid, p_cs, p_xmin, p_ymin, XPf, YPf)

    ok = (sok & pok).reshape(out_shape)
    sv = sv.reshape(out_shape)
    pv = pv.reshape(out_shape)

    CF_new = np.where(ok, CF * _transform(sv, t_slope, p_slope), CF)

    vf_raw = CF_new * _transform(pv, t_pctl_v, p_pctl)
    rf_raw = CF_new * _transform(1.0 - pv, t_pctl_r, p_pctl)

    VF = np.where(ok, 1.0 - _transform(vf_raw, 0.3, 4.0), 0.0)
    RF = np.where(ok, 1.0 - _transform(rf_raw, 0.3, 4.0), 0.0)

    return VF, RF, ok, CF_new


# ---------------------------------------------------------------------------
# CMRVBF::Get_MRVBF
# ---------------------------------------------------------------------------
def _blend(level: int, MRVBF, MRRTF, mask_valid, VF, RF, level_ok):
    t = 0.4
    p = np.log((level - 0.5) / 0.1) / np.log(1.5)

    do = mask_valid & level_ok

    w = 1.0 - _transform(VF, t, p)
    new_mrvbf = np.where(do, w * (level - 1 + VF) + (1.0 - w) * MRVBF, MRVBF)

    w2 = 1.0 - _transform(RF, t, p)
    new_mrrtf = np.where(do, w2 * (level - 1 + RF) + (1.0 - w2) * MRRTF, MRRTF)

    return new_mrvbf, new_mrrtf


# ---------------------------------------------------------------------------
# CMRVBF::Get_Classified
# ---------------------------------------------------------------------------
def _classify(a: np.ndarray) -> np.ndarray:
    valid = ~np.isnan(a)
    conds = [a < 0.5, a < 1.5, a < 2.5, a < 3.5, a < 4.5, a < 5.5]
    choices = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    out = np.select(conds, choices, default=6.0)
    return np.where(valid, out, np.nan)


# ---------------------------------------------------------------------------
# CMRVBF::On_Execute
# ---------------------------------------------------------------------------
def mrvbf(
    dem: np.ndarray,
    cellsize: float,
    xmin: float = 0.0,
    ymin: float = 0.0,
    nodata_mask: np.ndarray | None = None,
    t_slope: float = 16.0,
    t_pctl_v: float = 0.40,
    t_pctl_r: float = 0.35,
    p_slope: float = 4.0,
    p_pctl: float = 3.0,
    max_res: float = 100.0,
    classify: bool = False,
    verbose: bool = False,
):
    """
    Multiresolution Index of Valley Bottom Flatness (and its complementary
    Ridge Top Flatness), translated from SAGA's mrvbf.cpp (CMRVBF,
    ta_morphometry tool id 8).

    Parameters
    ----------
    dem : 2D array (ny, nx)
        Elevation grid; row 0 = y = ymin (row index increases northward).
    cellsize, xmin, ymin :
        DEM grid geometry (SAGA's Get_System()).
    nodata_mask : 2D bool array, optional
        True where dem is NoData. If omitted, non-finite dem cells are
        treated as NoData.
    t_slope, t_pctl_v, t_pctl_r, p_slope, p_pctl, max_res, classify :
        Same meaning as the SAGA tool's T_SLOPE / T_PCTL_V / T_PCTL_R /
        P_SLOPE / P_PCTL / MAX_RES / CLASSIFY parameters.
    verbose : bool
        Print per-level diagnostics (level, resolution, threshold slope).

    Returns
    -------
    (mrvbf, mrrtf) : tuple of 2D float64 arrays (ny, nx)
        NaN where the DEM was NoData (or where sampling failed at every
        level, matching SAGA's NoData propagation).
    """
    dem = np.asarray(dem, dtype=np.float64)
    ny, nx = dem.shape

    if nodata_mask is not None:
        valid0 = ~np.asarray(nodata_mask, dtype=bool)
    else:
        valid0 = np.isfinite(dem)

    xrange0 = (nx - 1) * cellsize
    yrange0 = (ny - 1) * cellsize
    diag = np.sqrt(xrange0 ** 2 + yrange0 ** 2)
    max_resolution = (max_res / 100.0) * diag

    xs = xmin + np.arange(nx) * cellsize
    ys = ymin + np.arange(ny) * cellsize
    XP, YP = np.meshgrid(xs, ys)  # (ny, nx)
    XPf, YPf = XP.ravel(), YP.ravel()

    CF = np.ones((ny, nx), dtype=np.float64)

    dem_cur, dem_cur_valid = dem.copy(), valid0.copy()
    dem_xmin, dem_ymin = xmin, ymin
    dem_cs, dem_w, dem_h = cellsize, nx, ny

    # ---- Level 1 --------------------------------------------------------
    level = 1
    T_Slope = t_slope
    resolution = cellsize

    if verbose:
        print(f"step: {level}, resolution: {resolution:.2f}, threshold slope {T_Slope:.2f}")

    Slopes, slopes_valid = _slope_percent(dem_cur, dem_cur_valid, dem_cs)
    Percentiles, pctl_valid = _percentile(dem_cur, dem_cur_valid, 3)

    VF, RF, ok, CF = _flatness(Slopes, slopes_valid, dem_cs, dem_xmin, dem_ymin,
                                Percentiles, pctl_valid, dem_cs, dem_xmin, dem_ymin,
                                CF, XPf, YPf, (ny, nx), T_Slope, p_slope, t_pctl_v, t_pctl_r, p_pctl)

    MRVBF, MRRTF = VF, RF          # Get_Flatness wrote directly into pMRVBF/pMRRTF at level 1
    out_valid = ok.copy()           # fixed for the rest of the run

    # ---- Level 2 --------------------------------------------------------
    T_Slope /= 2.0
    level += 1

    if verbose:
        print(f"step: {level}, resolution: {resolution:.2f}, threshold slope {T_Slope:.2f}")

    Percentiles, pctl_valid = _percentile(dem_cur, dem_cur_valid, 6)

    VF, RF, ok, CF = _flatness(Slopes, slopes_valid, dem_cs, dem_xmin, dem_ymin,
                                Percentiles, pctl_valid, dem_cs, dem_xmin, dem_ymin,
                                CF, XPf, YPf, (ny, nx), T_Slope, p_slope, t_pctl_v, t_pctl_r, p_pctl)

    MRVBF, MRRTF = _blend(level, MRVBF, MRRTF, out_valid, VF, RF, ok)

    # ---- Level 3+ ---------------------------------------------------------
    while resolution < max_resolution:
        resolution *= 3.0
        T_Slope /= 2.0
        level += 1

        if verbose:
            print(f"step: {level}, resolution: {resolution:.2f}, threshold slope {T_Slope:.2f}")

        prev_cs, prev_xmin, prev_ymin = dem_cs, dem_xmin, dem_ymin

        Smoothed, smoothed_valid = _smoothed(dem_cur, dem_cur_valid, 5)
        Slopes, slopes_valid = _slope_percent(Smoothed, smoothed_valid, prev_cs)

        new_xrange = (dem_w - 1) * prev_cs
        new_yrange = (dem_h - 1) * prev_cs
        new_nx = 2 + int(new_xrange / resolution)
        new_ny = 2 + int(new_yrange / resolution)

        dem_cur, dem_cur_valid = _nearest_resample(
            Smoothed, smoothed_valid, prev_cs, prev_xmin, prev_ymin,
            new_nx, new_ny, resolution, dem_xmin, dem_ymin)
        dem_cs, dem_w, dem_h = resolution, new_nx, new_ny

        Percentiles, pctl_valid = _percentile(dem_cur, dem_cur_valid, 6)

        VF, RF, ok, CF = _flatness(Slopes, slopes_valid, prev_cs, prev_xmin, prev_ymin,
                                    Percentiles, pctl_valid, dem_cs, dem_xmin, dem_ymin,
                                    CF, XPf, YPf, (ny, nx), T_Slope, p_slope, t_pctl_v, t_pctl_r, p_pctl)

        MRVBF, MRRTF = _blend(level, MRVBF, MRRTF, out_valid, VF, RF, ok)

    # ---- finalize -----------------------------------------------------
    MRVBF_out = np.where(out_valid, MRVBF, np.nan)
    MRRTF_out = np.where(out_valid, MRRTF, np.nan)

    if classify:
        MRVBF_out = _classify(MRVBF_out)
        MRRTF_out = _classify(MRRTF_out)

    return MRVBF_out, MRRTF_out


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    rng = np.random.default_rng(1)
    ny, nx = 60, 60
    yy, xx = np.mgrid[0:ny, 0:nx]
    dem = (
        20.0 * np.exp(-((xx - 20) ** 2 + (yy - 15) ** 2) / 200.0)
        + 0.05 * xx
        + 0.02 * yy
        + rng.normal(0, 0.2, (ny, nx))
    )

    mv, mr = mrvbf(dem, cellsize=10.0, xmin=0.0, ymin=0.0, verbose=True)
    print("MRVBF range:", np.nanmin(mv), np.nanmax(mv))
    print("MRRTF range:", np.nanmin(mr), np.nanmax(mr))
    print("NaN count MRVBF:", np.isnan(mv).sum(), "/", mv.size)
