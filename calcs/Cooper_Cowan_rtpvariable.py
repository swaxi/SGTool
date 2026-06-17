r"""
RTP with variable geomagnetic parameters
Uses a Taylor series expansion in the space domain
GRJ Cooper July 2002 - School of Geosciences, University of the Witwatersrand
Python conversion: follows line-by-line logic of original MATLAB (rtpvariable.m)

Any mistakes due to conversion by Mark Jessell, UWA 2026-06.  Please report any issues.

Run with QGIS Python 3.12 (has scipy + GDAL):
  "<QGIS>\apps\Python312\python.exe" Cooper_Cowan_rtpvariable.py
"""

import numpy as np
import matplotlib.pyplot as plt
from osgeo import gdal, osr


# =============================================================================
# NaN helpers (ported from GeophysicalProcessor.fill_nan / restore_nan)
# =============================================================================

def _fill_nan(data):
    """Replace NaN values via 2-D nearest-neighbour interpolation."""
    from scipy.interpolate import griddata
    nan_mask = np.isnan(data)
    if not np.any(nan_mask) or np.all(nan_mask):
        return data.copy(), nan_mask
    filled = data.copy()
    rows, cols = data.shape
    y, x = np.mgrid[0:rows, 0:cols]
    valid = ~nan_mask
    try:
        filled[nan_mask] = griddata(
            np.column_stack((y[valid].ravel(), x[valid].ravel())),
            data[valid].ravel(),
            np.column_stack((y[nan_mask].ravel(), x[nan_mask].ravel())),
            method='nearest',
            fill_value=np.nanmean(data),
        )
    except Exception:
        filled[nan_mask] = np.nanmean(data)
    return filled, nan_mask


def _restore_nan(data, nan_mask):
    result = data.copy()
    result[nan_mask] = np.nan
    return result


# =============================================================================
# Main function
# =============================================================================

def rtpvariable(magc, nr, nc, inc_corners, dec_corners, incv=None, decv=None,
                inc_center=None, dec_center=None):
    """
    magc        : input data as a numpy array, shape (nr, nc)
    nr, nc      : number of rows and columns (must match magc.shape)
    inc_corners : inclination (degrees) at the four grid corners in order:
                    [row=nr col=1, row=nr col=nc, row=1 col=1, row=1 col=nc]
                  i.e. [NW, NE, SW, SE] when displayed with axis xy (row nr = north)
    dec_corners : declination/azimuth (degrees) at the same four corners,
                  same order as inc_corners
    inc_center  : optional inclination (degrees) at grid centre for a 5-point fit
    dec_center  : optional declination (degrees) at grid centre for a 5-point fit
    incv        : optional full-grid inclination array (nr x nc). When supplied,
                  inc_corners and inc_center are ignored for inclination
    decv        : optional full-grid declination array (nr x nc). When supplied,
                  dec_corners and dec_center are ignored for declination
    """
    dtr = np.pi / 180

    # Preserve valid-data mean before filling (used as neutral fill for derivatives)
    data_mean = float(np.nanmean(magc))

    # Fill NaN with nearest-neighbour for the base RTP — this gives a smooth
    # transition at the survey boundary and good base-RTP quality.
    magc, nan_mask = _fill_nan(magc)

    # For the Taylor derivative computations use a constant (mean) fill instead.
    # Nearest-neighbour fill creates blocky step patterns inside the NaN region;
    # those steps are amplified by numerical differentiation and produce a spurious
    # smooth gradient in the correction terms.  A constant fill has zero derivatives,
    # so it contributes nothing to drtp1_inc / drtp1_dec.
    if np.any(nan_mask):
        magc_deriv = magc.copy()
        magc_deriv[nan_mask] = data_mean
    else:
        magc_deriv = magc

    xi, yi = np.meshgrid(np.arange(1, nc + 1), np.arange(1, nr + 1))

    # Bump basis: zero at all four corners, 1 at grid centre — used when inc_center/dec_center given
    r_mid = (nr + 1) / 2.0
    c_mid = (nc + 1) / 2.0
    norm = ((nr - 1) / 2.0) ** 2 * ((nc - 1) / 2.0) ** 2
    bump = (yi - 1) * (yi - nr) * (xi - 1) * (xi - nc) / norm

    corner_mat = np.array([[1, nr, 1,  nr * 1 ],
                            [1, nr, nc, nr * nc],
                            [1, 1,  1,  1  * 1 ],
                            [1, 1,  nc, 1  * nc]], dtype=float)

    if incv is None:
        print('Generating grid of inclination data from values at grid corners')
        p = np.linalg.solve(corner_mat, np.array(inc_corners, dtype=float))
        incv = p[0] + p[1] * yi + p[2] * xi + p[3] * yi * xi
        if inc_center is not None:
            bilinear_at_centre = p[0] + p[1]*r_mid + p[2]*c_mid + p[3]*r_mid*c_mid
            delta_inc = inc_center - bilinear_at_centre
            print(f'  inc bilinear at centre: {bilinear_at_centre:.2f}°, IGRF centre: {inc_center:.2f}°, delta: {delta_inc:.2f}°')
            incv = incv + delta_inc * bump

    if decv is None:
        print('Generating grid of declination data from values at grid corners')
        p = np.linalg.solve(corner_mat, np.array(dec_corners, dtype=float))
        decv = p[0] + p[1] * yi + p[2] * xi + p[3] * yi * xi
        if dec_center is not None:
            bilinear_at_centre = p[0] + p[1]*r_mid + p[2]*c_mid + p[3]*r_mid*c_mid
            delta_dec = dec_center - bilinear_at_centre
            print(f'  dec bilinear at centre: {bilinear_at_centre:.2f}°, IGRF centre: {dec_center:.2f}°, delta: {delta_dec:.2f}°')
            decv = decv + delta_dec * bump

    incv_deg = incv.copy()   # preserve degrees for optional export
    decv_deg = decv.copy()
    incv = incv * dtr;  inc = incv[nr // 2, nc // 2]   # inclination at centre of grid
    decv = decv * dtr;  azimuth = decv[nr // 2, nc // 2]   # declination at centre of grid

    nmax = max([nr, nc]);  npts = 2 ** int(np.ceil(np.log2(nmax)))
    cdiff = int(np.floor((npts - nc) / 2));  rdiff = int(np.floor((npts - nr) / 2))
    print('Computing rtp')
    rtpdata = rtp(magc, npts, nr, nc, rdiff, cdiff, azimuth, inc)
    incv = incv - inc
    decv = decv - azimuth

    print('Computing inc dx1')
    drtp1_inc = taylortp1(magc_deriv, npts, nr, nc, rdiff, cdiff, azimuth, inc)
    print('Computing inc dx2')
    drtp2_inc = taylortp1(drtp1_inc, npts, nr, nc, rdiff, cdiff, azimuth, inc)
    print('Computing inc dx3')
    drtp3_inc = taylortp1(drtp2_inc, npts, nr, nc, rdiff, cdiff, azimuth, inc)

    print('Computing dec dx1')
    drtp1_dec = taylordec1(magc_deriv, npts, nr, nc, rdiff, cdiff, azimuth, inc)
    print('Computing dec dx2')
    drtp2_dec = taylordec1(drtp1_dec, npts, nr, nc, rdiff, cdiff, azimuth, inc)
    print('Computing dec dx3')
    drtp3_dec = taylordec1(drtp2_dec, npts, nr, nc, rdiff, cdiff, azimuth, inc)

    print('Computing rtpvar')
    varrtp = (rtpdata
              + drtp1_inc * incv + 0.5 * drtp2_inc * incv ** 2 + 0.1666666 * drtp3_inc * incv ** 3
              + drtp1_dec * decv + 0.5 * drtp2_dec * decv ** 2 + 0.1666666 * drtp3_dec * decv ** 3)

    if np.any(nan_mask):
        rtpdata = _restore_nan(rtpdata, nan_mask)
        varrtp  = _restore_nan(varrtp,  nan_mask)

    return rtpdata, varrtp, incv_deg, decv_deg


# =============================================================================
# Helper: RTP derivatives via Taylor series
# =============================================================================

def taylortp1(magc, npts, nr, nc, rdiff, cdiff, azimuth, inc):
    dtheta = 0.01
    rtpdata1 = rtp(magc, npts, nr, nc, rdiff, cdiff, azimuth, inc)
    rtpdata2 = rtp(magc, npts, nr, nc, rdiff, cdiff, azimuth, inc + dtheta)
    drtp = (rtpdata2 - rtpdata1) / dtheta
    return drtp


def taylordec1(magc, npts, nr, nc, rdiff, cdiff, azimuth, inc):
    dtheta = 0.01
    rtpdata1 = rtp(magc, npts, nr, nc, rdiff, cdiff, azimuth, inc)
    rtpdata2 = rtp(magc, npts, nr, nc, rdiff, cdiff, azimuth + dtheta, inc)
    drtp = (rtpdata2 - rtpdata1) / dtheta
    return drtp


# =============================================================================
# Helper: core RTP computation
# =============================================================================

def rtp(magc, npts, nr, nc, rdiff, cdiff, azimuth, inc):
    gf = taper2d(magc, npts, nc, nr, cdiff, rdiff)
    f = np.fft.fft2(gf, s=(npts, npts))
    fx = f.copy()

    # Direction cosines
    L = np.cos(azimuth) * np.cos(inc)
    M = np.sin(azimuth) * np.cos(inc)
    N = np.sin(inc)
    N2 = N * N

    wn = 2.0 * np.pi / (npts - 1)
    half = npts // 2
    ii = np.arange(1, half + 1)
    jj = np.arange(1, half + 1)
    U = ii[:, None] * wn        # (half, 1) — row wavenumber
    V = jj[None, :] * wn        # (1, half) — col wavenumber
    temp = U * U + V * V        # (half, half)
    temp_15 = temp ** 1.5

    def _apply(block, t1):
        t2 = N2 * temp - t1 * t1
        base = t2 * t2 + 4.0 * N2 * t1 * t1 * temp
        Rp = temp * t2 / base
        Ip = -2.0 * N * t1 * temp_15 / base
        re, im = block.real, block.imag
        return (Rp * re - Ip * im) + 1j * (Rp * im + Ip * re)

    # Four symmetric wavenumber quadrants — disjoint, cover the full spectrum
    fx[np.ix_(ii - 1,    jj - 1)]    = _apply(f[np.ix_(ii - 1,    jj - 1)],     L * U + M * V)
    fx[np.ix_(npts - ii, jj - 1)]    = _apply(f[np.ix_(npts - ii, jj - 1)],    -L * U + M * V)
    fx[np.ix_(ii - 1,    npts - jj)] = _apply(f[np.ix_(ii - 1,    npts - jj)],  L * U - M * V)
    fx[np.ix_(npts - ii, npts - jj)] = _apply(f[np.ix_(npts - ii, npts - jj)], -L * U - M * V)

    fxinv = np.fft.ifft2(fx);  rtpdata = fxinv[rdiff:rdiff + nr, cdiff:cdiff + nc].real
    return rtpdata


# =============================================================================
# Helper: 2D cosine taper
# =============================================================================

def taper2d(g, npts, nc, nr, cdiff, rdiff):
    """GRJ Cooper 1994"""
    gf = np.zeros((npts, npts));  gf[rdiff:rdiff + nr, cdiff:cdiff + nc] = g

    ii = np.arange(1, rdiff + 1)                                                # 1..rdiff
    wr = (1 + np.sin(-np.pi / 2 + (ii - 1) * np.pi / rdiff)) * 0.5            # row weights
    jj = np.arange(1, cdiff + 1)                                                # 1..cdiff
    wc = (1 + np.sin(-np.pi / 2 + (jj - 1) * np.pi / cdiff)) * 0.5            # col weights

    dc = np.arange(cdiff, cdiff + nc)   # data column indices (0-based)
    dr = np.arange(rdiff, rdiff + nr)   # data row indices (0-based)

    # Top and bottom row tapers over data columns
    gf[np.ix_(ii - 1,    dc)] = gf[np.ix_(2 * rdiff - ii,            dc)] * wr[:, None]
    gf[np.ix_(npts - ii, dc)] = gf[np.ix_(npts - 2 * rdiff + ii - 1, dc)] * wr[:, None]

    # Left and right column tapers over data rows
    gf[np.ix_(dr, jj - 1)]         = gf[np.ix_(dr, 2 * cdiff - jj)]             * wc[None, :]
    gf[np.ix_(dr, npts - jj - 1)]  = gf[np.ix_(dr, npts - 2 * cdiff + jj - 1)] * wc[None, :]

    # Corner tapers: left edge — reads boundary rows set by column tapers above
    lc = np.arange(cdiff)
    gf[np.ix_(ii - 1,        lc)] = gf[rdiff,            lc][None, :] * wr[:, None]
    gf[np.ix_(npts - ii - 1, lc)] = gf[npts - rdiff - 1, lc][None, :] * wr[:, None]

    # Corner tapers: right edge (incl. last data column per original loop range)
    rc = np.arange(cdiff + nc - 1, npts)
    gf[np.ix_(ii - 1,        rc)] = gf[rdiff,            rc][None, :] * wr[:, None]
    gf[np.ix_(npts - ii - 1, rc)] = gf[npts - rdiff - 1, rc][None, :] * wr[:, None]

    return gf


# =============================================================================
# GeoTIFF export
# =============================================================================

def save_geotiff(data, output_path, lat_top=-30.0, lat_bottom=-90.0,
                 lon_left=0.0, lon_right=None, crs_epsg=4326):
    """Save a 2D array as GeoTIFF.

    Coordinate origin from lines 30-33 of rtpvariable.m (EPSG:4326):
      - lat_top    (-30): inclination yval at row=nr corners (lines 30-31)
                          With 'axis xy', row nr appears at geographic NORTH.
      - lat_bottom (-90): inclination yval at row=1 corners (lines 32-33)
                          Row 1 appears at geographic SOUTH with 'axis xy'.
      - lon_left/right: NOT specified in lines 30-33. Default gives square pixels
                        (lon span = nc * lat_span / nr). Adjust as needed.

    data rows: row 0 = geographic south (MATLAB row 1, 'axis xy' bottom).
    GeoTIFF rows: row 0 = geographic north. Array is flipped vertically on write.
    """
    nr, nc = data.shape
    if lon_right is None:
        lat_span = lat_top - lat_bottom           # degrees (positive)
        deg_per_pixel = lat_span / nr
        lon_right = lon_left + nc * deg_per_pixel

    x_res = (lon_right - lon_left) / nc
    y_res = (lat_top - lat_bottom) / nr           # positive step size

    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(output_path, nc, nr, 1, gdal.GDT_Float32)

    # GeoTransform: (top_left_x, x_pixel_size, x_rotation, top_left_y, y_rotation, -y_pixel_size)
    ds.SetGeoTransform((lon_left, x_res, 0.0, lat_top, 0.0, -y_res))

    srs = osr.SpatialReference()
    srs.ImportFromEPSG(crs_epsg)
    ds.SetProjection(srs.ExportToWkt())

    # Flip vertically: MATLAB row 1 (Python index 0) is geographic south;
    # GeoTIFF row 0 must be geographic north.
    band = ds.GetRasterBand(1)
    band.WriteArray(np.flipud(data).astype(np.float32))
    band.SetNoDataValue(float('nan'))
    band.FlushCache()
    ds = None

    print(f'Saved: {output_path}')
    print(f'  Size : {nr} rows x {nc} cols')
    print(f'  Lat  : {lat_bottom} to {lat_top}  (south to north)')
    print(f'  Lon  : {lon_left} to {lon_right:.4f}')
    print(f'  CRS  : EPSG:{crs_epsg}')
    print(f'  NOTE : Longitude bounds are inferred (not in lines 30-33). Adjust lon_left/lon_right if needed.')


# =============================================================================
# Entry point
# =============================================================================

if __name__ == '__main__':
    import os
    import scipy.io as sio

    here = os.path.dirname(os.path.abspath(__file__))
    mat_dir = os.path.join(here, '..', '..', 'Dropbox', 'WAXI4', 'gis',
                           'SGTool_comparison', 'Cooper_DRTP')

    mat = sio.loadmat(os.path.join(mat_dir, 'rtpdemodata.mat'))
    magc = mat['magc']
    nr, nc = magc.shape

    # Corner values matching the original MATLAB demo
    inc_corners = (-30, -30, -90, -90)   # degrees
    dec_corners = (0, 0, 0, 0)           # degrees

    save_geotiff(magc, os.path.join(mat_dir, 'rtpdemodata.tif'))

    rtpdata, varrtp, incv_deg, decv_deg = rtpvariable(magc, nr, nc, inc_corners, dec_corners)

    save_geotiff(rtpdata,   os.path.join(mat_dir, 'rtpdata_simple.tif'))
    save_geotiff(varrtp,    os.path.join(mat_dir, 'rtpdata_differential.tif'))
    save_geotiff(incv_deg,  os.path.join(mat_dir, 'rtpdata_inc_field.tif'))
    save_geotiff(decv_deg,  os.path.join(mat_dir, 'rtpdata_dec_field.tif'))
