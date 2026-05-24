"""
arcgis_utils.py  –  Shared raster I/O and path helpers for SGTool ArcGIS.

Tries arcpy first (inside ArcGIS Pro).  Falls back to osgeo.gdal when
arcpy is not available so the GUI can also run stand-alone or in QGIS.
"""

import os
import sys
import numpy as np

# ---------------------------------------------------------------------------
# Path setup  – make calcs importable as a proper package
# ---------------------------------------------------------------------------
_THIS_DIR  = os.path.dirname(os.path.abspath(__file__))
_SGTOOL    = os.path.dirname(_THIS_DIR)          # …/sgtool
_CALCS     = os.path.join(_SGTOOL, "calcs")
_WORMS     = os.path.join(_CALCS,  "worms")
_IGRF      = os.path.join(_CALCS,  "igrf")
_EULER     = os.path.join(_CALCS,  "euler")

for _p in (_SGTOOL, _CALCS, _WORMS, _IGRF, _EULER):
    if _p not in sys.path:
        sys.path.insert(0, _p)

NODATA = -9999.0


# ---------------------------------------------------------------------------
# Raster read
# ---------------------------------------------------------------------------
def raster_to_numpy(raster_path):
    """
    Read a single-band raster.
    Returns (array2d, nodata, ll_x, ll_y, cell_x, cell_y, spatial_ref)
        spatial_ref  – arcpy.SpatialReference when arcpy available, else WKT str
        array2d      – float64 with NaN replacing nodata
    """
    raster_path = str(raster_path)

    # --- arcpy path ---
    try:
        import arcpy
        desc   = arcpy.Describe(raster_path)
        raster = arcpy.Raster(raster_path)
        nodata = float(raster.noDataValue) if raster.noDataValue is not None else NODATA
        cx     = desc.meanCellWidth
        cy     = desc.meanCellHeight
        ll_x   = desc.extent.XMin
        ll_y   = desc.extent.YMin
        sr     = desc.spatialReference
        arr    = arcpy.RasterToNumPyArray(raster_path, nodata_to_value=np.nan).astype(np.float64)
        # Some rasters return int; nodata may already be filled – ensure NaN
        arr[arr == nodata] = np.nan
        return arr, nodata, ll_x, ll_y, cx, cy, sr
    except Exception:
        pass

    # --- GDAL path ---
    from osgeo import gdal
    ds   = gdal.Open(raster_path, gdal.GA_ReadOnly)
    if ds is None:
        raise IOError(f"Cannot open raster: {raster_path}")
    band   = ds.GetRasterBand(1)
    arr    = band.ReadAsArray().astype(np.float64)
    nodata = band.GetNoDataValue()
    if nodata is None:
        nodata = NODATA
    arr[arr == nodata] = np.nan
    gt   = ds.GetGeoTransform()
    cx   = abs(gt[1])
    cy   = abs(gt[5])
    ll_x = gt[0]
    ll_y = gt[3] + gt[5] * ds.RasterYSize
    wkt  = ds.GetProjection()
    ds   = None
    return arr, float(nodata), ll_x, ll_y, cx, cy, wkt


# ---------------------------------------------------------------------------
# Raster write
# ---------------------------------------------------------------------------
def numpy_to_raster(array, output_path, ll_x, ll_y,
                    cell_x, cell_y, spatial_ref=None,
                    nodata=NODATA):
    """
    Write a 2-D float array to a GeoTIFF.
    NaNs in *array* are replaced by *nodata* in the output.
    """
    output_path = str(output_path)
    if os.path.exists(output_path):
        try:
            os.remove(output_path)
        except OSError:
            pass
    out = np.where(np.isnan(array), nodata, array).astype(np.float32)

    # --- arcpy path ---
    try:
        import arcpy
        ll    = arcpy.Point(ll_x, ll_y)
        robj  = arcpy.NumPyArrayToRaster(out, ll, cell_x, cell_y, nodata)
        if spatial_ref is not None:
            arcpy.DefineProjection_management(robj, spatial_ref)
        robj.save(output_path)
        return
    except Exception:
        pass

    # --- GDAL path ---
    from osgeo import gdal, osr
    rows, cols = out.shape
    driver = gdal.GetDriverByName("GTiff")
    ds = driver.Create(output_path, cols, rows, 1, gdal.GDT_Float32)
    nrows = rows
    gt = (ll_x, cell_x, 0.0,
          ll_y + nrows * cell_y, 0.0, -cell_y)
    ds.SetGeoTransform(gt)
    if isinstance(spatial_ref, str) and spatial_ref:
        ds.SetProjection(spatial_ref)
    elif spatial_ref is not None:
        try:
            ds.SetProjection(spatial_ref.exportToString())
        except Exception:
            pass
    band = ds.GetRasterBand(1)
    band.WriteArray(out)
    band.SetNoDataValue(float(nodata))
    ds.FlushCache()
    ds = None


# ---------------------------------------------------------------------------
# Convenience: read → process → write wrapper
# ---------------------------------------------------------------------------
def run_raster_tool(input_raster, output_raster, func):
    """
    Read *input_raster*, call func(array, dx, dy) → result, write to *output_raster*.
    """
    arr, nodata, ll_x, ll_y, cx, cy, sr = raster_to_numpy(input_raster)
    result = func(arr, cx, cy)
    numpy_to_raster(result, output_raster, ll_x, ll_y, cx, cy, sr, nodata)


# ---------------------------------------------------------------------------
# GeophysicalProcessor factory  (handles relative import inside calcs/)
# ---------------------------------------------------------------------------
def make_processor(dx, dy, buffer_size=10):
    """Return a GeophysicalProcessor(dx, dy, buffer_size) instance."""
    from calcs.GeophysicalProcessor import GeophysicalProcessor
    return GeophysicalProcessor(dx, dy, buffer_size)


# ---------------------------------------------------------------------------
# Raster centre → geographic lat/lon  (header only, no pixel read)
# ---------------------------------------------------------------------------
def raster_center_latlon(raster_path):
    """
    Return (lat, lon) of the raster centre in decimal degrees, or None on failure.
    Reprojects to geographic WGS-84 if the raster is in a projected CRS.
    """
    raster_path = str(raster_path)

    # --- arcpy path ---
    try:
        import arcpy
        desc  = arcpy.Describe(raster_path)
        ext   = desc.extent
        cx    = (ext.XMin + ext.XMax) / 2.0
        cy    = (ext.YMin + ext.YMax) / 2.0
        sr    = desc.spatialReference
        if sr.type == "Projected":
            pt  = arcpy.PointGeometry(arcpy.Point(cx, cy), sr)
            geo = pt.projectAs(arcpy.SpatialReference(4326))
            return float(geo.centroid.Y), float(geo.centroid.X)
        return float(cy), float(cx)
    except Exception:
        pass

    # --- GDAL path ---
    try:
        from osgeo import gdal, osr
        ds  = gdal.Open(raster_path, gdal.GA_ReadOnly)
        if ds is None:
            return None
        gt   = ds.GetGeoTransform()
        cols = ds.RasterXSize
        rows = ds.RasterYSize
        cx   = gt[0] + gt[1] * cols / 2.0 + gt[2] * rows / 2.0
        cy   = gt[3] + gt[4] * cols / 2.0 + gt[5] * rows / 2.0
        wkt  = ds.GetProjection()
        ds   = None
        if not wkt:
            return float(cy), float(cx)
        src_srs = osr.SpatialReference()
        src_srs.ImportFromWkt(wkt)
        if not src_srs.IsGeographic():
            geo_srs = osr.SpatialReference()
            geo_srs.ImportFromEPSG(4326)
            geo_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
            ct  = osr.CoordinateTransformation(src_srs, geo_srs)
            lon, lat, _ = ct.TransformPoint(cx, cy)
            return float(lat), float(lon)
        return float(cy), float(cx)
    except Exception:
        return None


# ---------------------------------------------------------------------------
# IGRF helper
# ---------------------------------------------------------------------------
def calc_igrf(lat, lon, altitude_km, date_decimal):
    """
    Calculate IGRF magnetic field parameters.
    Returns (inclination_deg, declination_deg, intensity_nT) or None on failure.
    """
    try:
        from calcs.igrf import igrf_utils as iu

        shc_dir = os.path.join(_IGRF, "SHC_files")
        shc_files = sorted(
            [f for f in os.listdir(shc_dir) if f.endswith(".SHC")],
            reverse=True
        )
        if not shc_files:
            return None
        shc_path = os.path.join(shc_dir, shc_files[0])

        # load_shcfile returns a populated igrf_utils instance
        igrf = iu.igrf_utils(None, None, None).load_shcfile(shc_path)

        # Interpolate coefficients to the target decimal year
        # igrf.time: shape (N,)   igrf.coeffs: shape (n_coefficients, N)
        coeffs_t = np.array([
            np.interp(date_decimal, igrf.time, igrf.coeffs[i])
            for i in range(igrf.coeffs.shape[0])
        ])

        # Geodetic → geocentric; returns (radius_km, geocentric_colat, sd, cd)
        gdcolat = 90.0 - lat
        r, thc, sd, cd = igrf.gg_to_geo(float(altitude_km), float(gdcolat))

        # Field components in geocentric frame
        Br, Bt, Bp = igrf.synth_values(
            coeffs_t, r,
            np.float64(thc), np.float64(lon)
        )

        # Rotate geocentric → geodetic (BGS convention)
        X = -float(Bt) * float(cd) - float(Br) * float(sd)   # North
        Y =  float(Bp)                                          # East
        Z =  float(Bt) * float(sd) - float(Br) * float(cd)   # Down

        dec, _H, inc, F = igrf.xyz2dhif(X, Y, Z)
        return float(inc), float(dec), float(F)
    except Exception:
        return None
