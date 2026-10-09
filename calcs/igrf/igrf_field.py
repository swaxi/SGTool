"""
Geomagnetic field (IGRF) values for a grid, with no QGIS dependency, so they
can be worked out in a background task or a Processing algorithm.

calc_igrf() is the plugin's calcIGRF() moved here unchanged; grid_field() adds
the step the plugin dialog does around it (grid centre and corners from the
extent and CRS code, converted to latitude/longitude).
"""

import os
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy import interpolate

from .igrf_utils import igrf_utils as IGRF

_SHC_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "SHC_files")
IGRF_GENERATION = "14"


def decimal_year(year, month, day):
    """A date as a decimal year (e.g. 2023.5 for mid-2023)."""
    date = datetime(year, month, day)
    start_of_year = datetime(year, 1, 1)
    end_of_year = datetime(year + 1, 1, 1)
    days_in_year = (end_of_year - start_of_year).days
    return year + (date - start_of_year).days / days_in_year


def calc_igrf(date, alt, lat, lon):
    """Inclination, declination (degrees) and total intensity (nT) of the
    International Geomagnetic Reference Field.

    date: decimal year; alt: altitude above sea level in km; lat, lon: degrees.
    """
    itype = 1
    d1 = d2 = d3 = None
    colat = 90 - lat
    iut = IGRF(d1, d2, d3)

    # Load in the file of coefficients
    igrf_file = str(Path(os.path.join(_SHC_DIR, "IGRF" + IGRF_GENERATION + ".SHC")))
    igrf = iut.load_shcfile(igrf_file, None)

    # Interpolate the geomagnetic coefficients to the desired date
    f = interpolate.interp1d(igrf.time, igrf.coeffs, fill_value="extrapolate")
    coeffs = f(date)

    # Main field B_r, B_theta and B_phi for the location
    Br, Bt, Bp = iut.synth_values(coeffs.T, alt, colat, lon, igrf.parameters["nmax"])

    # For the SV, find the 5 year period in which the date lies and compute
    # the SV within that period (IGRF has constant SV between each 5 year period)
    epoch = (date - 1900) // 5
    epoch_start = epoch * 5
    coeffs_sv = f(1900 + epoch_start + 1) - f(1900 + epoch_start)
    Brs, Bts, Bps = iut.synth_values(coeffs_sv.T, alt, colat, lon, igrf.parameters["nmax"])

    # Main field coefficients from the start of each five epoch
    coeffsm = f(1900 + epoch_start)
    Brm, Btm, Bpm = iut.synth_values(coeffsm.T, alt, colat, lon, igrf.parameters["nmax"])

    # Rearrange to X, Y, Z components
    X = -Bt
    Y = Bp
    Z = -Br
    dX = -Bts
    dZ = -Brs
    Xm = -Btm
    Zm = -Brm
    if itype == 1:
        alt, colat, sd, cd = iut.gg_to_geo(alt, colat)

    # Rotate back to geodetic coords if needed
    if itype == 1:
        t = X
        X = X * cd + Z * sd
        Z = Z * cd - t * sd
        t = dX
        dX = dX * cd + dZ * sd
        dZ = dZ * cd - t * sd
        t = Xm
        Xm = Xm * cd + Zm * sd
        Zm = Zm * cd - t * sd

    intensity = np.sqrt(X**2 + Y**2 + Z**2)
    dec, hoz, inc, eff = iut.xyz2dhif(X, Y, Z)
    return inc, dec, intensity


def _to_lonlat(authid):
    """Transformer from the grid's CRS (an 'EPSG:nnnn' code) to lon/lat."""
    from pyproj import CRS, Transformer

    return Transformer.from_crs(
        CRS.from_user_input(int(authid.split(":")[1])),
        CRS.from_user_input(4326),
        always_xy=True,
    )


def grid_field(extent, authid, date, alt=100.0, corners=False):
    """IGRF inclination and declination for a grid.

    extent: (xmin, ymin, xmax, ymax) in the grid's CRS; authid: e.g. 'EPSG:28350';
    date: decimal year; alt: km above sea level (the dialog uses 100).

    Returns (inc, dec) at the grid centre. With corners=True returns
    (inc_corners, dec_corners, inc_centre, dec_centre) with the corners in
    the order [NW, NE, SW, SE], as differential RTP expects.
    """
    xmin, ymin, xmax, ymax = extent
    proj = _to_lonlat(authid)
    lon_c, lat_c = proj.transform((xmin + xmax) / 2.0, (ymin + ymax) / 2.0)
    inc_c, dec_c, _ = calc_igrf(date, alt, lat_c, lon_c)
    if not corners:
        return inc_c, dec_c
    corners_xy = [(xmin, ymax), (xmax, ymax), (xmin, ymin), (xmax, ymin)]
    inc_corners, dec_corners = [], []
    for x, y in corners_xy:
        lon, lat = proj.transform(x, y)
        inc_v, dec_v, _ = calc_igrf(date, alt, lat, lon)
        inc_corners.append(inc_v)
        dec_corners.append(dec_v)
    return tuple(inc_corners), tuple(dec_corners), inc_c, dec_c
