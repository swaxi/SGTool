# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2024-2026 Mark Jessell
"""
Geomagnetic field (IGRF) values for a grid, with no QGIS dependency, so they
can be worked out in a background task or a Processing algorithm.

igrf_values() evaluates the field at many points at once (vectorised, in
chunks), and calc_igrf() at one point (or many); grid_igrf() gives the field
at every cell of a grid, as Remove IGRF and differential RTP use it, and
grid_field() at the grid's centre (and corners), as the IGRF button does.
Altitudes are in km above the WGS-84 ellipsoid; the grids use sea level.
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


_MODEL = None  # (igrf_utils instance, coefficient file, interpolator), loaded once


def _model():
    global _MODEL
    if _MODEL is None:
        iut = IGRF(None, None, None)
        igrf_file = str(Path(os.path.join(_SHC_DIR, "IGRF" + IGRF_GENERATION + ".SHC")))
        igrf = iut.load_shcfile(igrf_file, None)
        f = interpolate.interp1d(igrf.time, igrf.coeffs, fill_value="extrapolate")
        _MODEL = (iut, igrf, f)
    return _MODEL


def igrf_values(date, alt, lat, lon, chunk=100_000, callback=None):
    """Inclination, declination (degrees) and total intensity (nT) of the IGRF
    main field at many points at once (arrays of any matching shape).

    date: decimal year; alt: km above sea level (WGS-84); lat, lon: geodetic
    degrees. The points are converted to geocentric radius and colatitude, the
    field is evaluated there and rotated back to the geodetic frame. Points are
    done chunk at a time (memory); callback(fraction done) is called after each
    chunk and may raise to stop.
    """
    iut, igrf, f = _model()
    coeffs = f(date).T
    nmax = igrf.parameters["nmax"]
    lat, lon = np.broadcast_arrays(np.asarray(lat, dtype=float), np.asarray(lon, dtype=float))
    shape = lat.shape
    lat, lon = lat.ravel(), lon.ravel()
    inc = np.empty(lat.size)
    dec = np.empty(lat.size)
    F = np.empty(lat.size)
    for s in range(0, lat.size, chunk):
        e = min(s + chunk, lat.size)
        radius, theta, sd, cd = iut.gg_to_geo(np.full(e - s, float(alt)), 90.0 - lat[s:e])
        Br, Bt, Bp = iut.synth_values(coeffs, radius, theta, lon[s:e], nmax)
        X, Y, Z = -Bt, Bp, -Br
        X, Z = X * cd + Z * sd, Z * cd - X * sd  # geocentric -> geodetic
        dec[s:e], _hoz, inc[s:e], F[s:e] = iut.xyz2dhif(X, Y, Z)
        if callback is not None:
            callback(e / lat.size)
    return inc.reshape(shape), dec.reshape(shape), F.reshape(shape)


def calc_igrf(date, alt, lat, lon):
    """Inclination, declination (degrees) and total intensity (nT) of the
    International Geomagnetic Reference Field.

    date: decimal year; alt: altitude above sea level in km; lat, lon: degrees
    (numbers, or arrays for many points).
    """
    inc, dec, F = igrf_values(date, alt, lat, lon)
    if np.ndim(inc) == 0:
        return float(inc), float(dec), float(F)
    return inc, dec, F


def _to_lonlat(authid):
    """Transformer (with .transform(x, y) -> lon, lat) from the grid's CRS (an
    'EPSG:nnnn' code) to lon/lat. Uses pyproj, or GDAL where pyproj is not installed
    (the ArcGIS Pro Python environment)."""
    epsg = int(authid.split(":")[1])
    try:
        from pyproj import CRS, Transformer
    except ImportError:
        from osgeo import osr

        src, dst = osr.SpatialReference(), osr.SpatialReference()
        src.ImportFromEPSG(epsg)
        dst.ImportFromEPSG(4326)
        for ref in (src, dst):
            if hasattr(osr, "OAMS_TRADITIONAL_GIS_ORDER"):  # x = lon / easting, always
                ref.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        transform = osr.CoordinateTransformation(src, dst)

        class _OsrTransformer:
            @staticmethod
            def transform(x, y):
                lon, lat, _z = transform.TransformPoint(x, y)
                return lon, lat

        return _OsrTransformer()
    return Transformer.from_crs(
        CRS.from_user_input(epsg), CRS.from_user_input(4326), always_xy=True
    )


def grid_field(extent, authid, date, alt=0.0, corners=False):
    """IGRF inclination and declination for a grid.

    extent: (xmin, ymin, xmax, ymax) in the grid's CRS; authid: e.g. 'EPSG:28350';
    date: decimal year; alt: km above sea level.

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



def _lonlat_arrays(crs, x, y):
    """(lon, lat) arrays of points (x, y arrays) in a CRS given as 'EPSG:nnnn' or
    WKT. pyproj if it is installed, else GDAL (the ArcGIS Pro environment)."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    try:
        from pyproj import CRS, Transformer

        transformer = Transformer.from_crs(CRS.from_user_input(crs), CRS.from_user_input(4326),
                                           always_xy=True)
        lon, lat = transformer.transform(x, y)
        return np.asarray(lon, dtype=float), np.asarray(lat, dtype=float)
    except ImportError:
        pass
    from osgeo import osr

    src = osr.SpatialReference()
    if str(crs).upper().startswith("EPSG:"):
        src.ImportFromEPSG(int(str(crs).split(":")[1]))
    else:
        src.ImportFromWkt(str(crs))
    dst = osr.SpatialReference()
    dst.ImportFromEPSG(4326)
    for ref in (src, dst):
        if hasattr(osr, "OAMS_TRADITIONAL_GIS_ORDER"):  # x = lon / easting, always
            ref.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    out = np.array(osr.CoordinateTransformation(src, dst).TransformPoints(
        np.column_stack([x.ravel(), y.ravel()]).tolist()), dtype=float)
    return out[:, 0].reshape(x.shape), out[:, 1].reshape(x.shape)


def grid_igrf(gt, shape, crs, date, alt=0.0, callback=None, chunk=200_000):
    """IGRF inclination, declination (degrees) and total intensity (nT) at the
    centre of every cell of a north-up grid (each a (rows, cols) array).

    gt: GDAL geotransform of the north-up grid; shape: (rows, cols); crs: its
    CRS as 'EPSG:nnnn' or WKT; date: decimal year; alt: km above sea level.
    The cells are done a block of rows at a time (about `chunk` cells);
    callback(fraction done) is called after each block and may raise to stop.
    """
    rows, cols = int(shape[0]), int(shape[1])
    inc = np.empty((rows, cols))
    dec = np.empty((rows, cols))
    F = np.empty((rows, cols))
    x = gt[0] + (np.arange(cols) + 0.5) * gt[1]
    step = max(1, chunk // max(cols, 1))
    for r0 in range(0, rows, step):
        r1 = min(rows, r0 + step)
        y = gt[3] + (np.arange(r0, r1) + 0.5) * gt[5]
        X, Y = np.meshgrid(x, y)
        lon, lat = _lonlat_arrays(crs, X, Y)
        inc[r0:r1], dec[r0:r1], F[r0:r1] = igrf_values(date, alt, lat, lon)
        if callback is not None:
            callback(r1 / rows)
    return inc, dec, F


def looks_like_anomaly(data, field):
    """True if a grid meant to be total magnetic intensity sits far below the IGRF
    total field over it (its median under half the field's): it is then most
    likely an anomaly grid that already has the IGRF (or a regional) removed."""
    data = np.asarray(data, dtype=float)
    med = np.nanmedian(data) if np.isfinite(data).any() else np.nan
    return bool(np.isfinite(med) and abs(med) < 0.5 * float(np.nanmean(field)))
