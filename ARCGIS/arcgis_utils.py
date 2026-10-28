"""
arcgis_utils.py  –  Shared raster I/O and path helpers for SGTool ArcGIS.

Tries arcpy first (inside ArcGIS Pro).  Falls back to osgeo.gdal when
arcpy is not available so the GUI can also run stand-alone or in QGIS.
"""

import logging
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
    except Exception as e:
        print(f"SGTool: arcpy raster read unavailable, falling back to GDAL: {e}")

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

    # --- GDAL path (preferred: much lower overhead than arcpy.Raster.save(),
    #     which stays slow regardless of pyramid/statistics env settings
    #     since those only apply to geoprocessing tools, not Raster.save()).
    #     Output is always a plain .tif on disk, never a geodatabase raster,
    #     so GDAL is a safe primary path here. ---
    try:
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
            ds.SetProjection(spatial_ref.exportToString())
        band = ds.GetRasterBand(1)
        band.WriteArray(out)
        band.SetNoDataValue(float(nodata))
        ds.FlushCache()
        ds = None
        return
    except Exception as e:
        print(f"SGTool: GDAL raster write failed, falling back to arcpy: {e}")

    # --- arcpy fallback path ---
    import arcpy
    ll    = arcpy.Point(ll_x, ll_y)
    robj  = arcpy.NumPyArrayToRaster(out, ll, cell_x, cell_y, nodata)
    if spatial_ref is not None:
        arcpy.DefineProjection_management(robj, spatial_ref)
    robj.save(output_path)


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
    except Exception as e:
        print(f"SGTool: arcpy centre lat/lon unavailable, falling back to GDAL: {e}")

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


def get_noddy_mag_params(raster_path):
    """
    Read inclination, declination and intensity from the GeoTIFF metadata
    written when a Noddy grid is imported.
    Returns (inc, dec, intensity) as floats, or None if the metadata is absent.
    """
    try:
        from osgeo import gdal

        ds = gdal.Open(str(raster_path).split("|")[0], gdal.GA_ReadOnly)
        if ds is None:
            return None
        md = ds.GetMetadata()
        ds = None
        if not all(k in md for k in ("inclination", "declination", "intensity")):
            return None
        return (float(md["inclination"]), float(md["declination"]),
                float(md["intensity"]))
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Version, grid information and provenance (shared with the QGIS plugin)
# ---------------------------------------------------------------------------
def plugin_version():
    """SGTool version from metadata.txt (None if it cannot be read)."""
    try:
        with open(os.path.join(_SGTOOL, "metadata.txt"), encoding="utf-8") as f:
            for line in f:
                if line.startswith("version="):
                    return line.split("=", 1)[1].strip()
    except OSError:
        pass
    return None


def raster_header_info(path):
    """Size, cell size, extent and CRS of a raster, from its header only.

    Returns a dict (width, height, cell_x, cell_y, extent=(xmin, ymin, xmax, ymax),
    authid e.g. 'EPSG:28350' or '', crs_name, geographic) or None if the raster
    cannot be read."""
    try:
        from osgeo import gdal, osr

        ds = gdal.Open(str(path).split("|")[0], gdal.GA_ReadOnly)
        if ds is None:
            return None
        width, height = ds.RasterXSize, ds.RasterYSize
        gt = ds.GetGeoTransform()
        wkt = ds.GetProjection()
        ds = None
        x0, x1 = gt[0], gt[0] + gt[1] * width
        y0, y1 = gt[3], gt[3] + gt[5] * height
        authid, name, geographic = "", "", False
        if wkt:
            srs = osr.SpatialReference()
            srs.ImportFromWkt(wkt)
            name = srs.GetName() or ""
            geographic = bool(srs.IsGeographic())
            try:
                srs.AutoIdentifyEPSG()
            except Exception:
                logging.getLogger(__name__).debug("Ignored: ignored", exc_info=True)
            code = srs.GetAuthorityCode(None)
            if code and (srs.GetAuthorityName(None) or "").upper() == "EPSG":
                authid = f"EPSG:{code}"
            else:
                # the file's WKT may carry no code: take a confident database match
                try:
                    for match, confidence in srs.FindMatches() or []:
                        if (confidence >= 90 and match.GetAuthorityCode(None)
                                and (match.GetAuthorityName(None) or "").upper() == "EPSG"):
                            authid = f"EPSG:{match.GetAuthorityCode(None)}"
                            break
                except Exception:
                    logging.getLogger(__name__).debug("Ignored: ignored", exc_info=True)
        return {
            "width": width, "height": height,
            "cell_x": abs(gt[1]), "cell_y": abs(gt[5]),
            "extent": (min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1)),
            "authid": authid, "crs_name": name, "geographic": geographic,
        }
    except Exception:
        return None


def worms_reference(raster_path):
    """(layer, epsg) for GeophysicalProcessor.bsdwormer: a LayerShim with the
    raster's extent and EPSG code (it georeferences the padded grid) and the
    EPSG number (for the worms shapefile)."""
    from calcs.layer_shim import LayerShim

    info = raster_header_info(raster_path)
    if not info or not info["authid"]:
        raise ValueError(
            "Worms needs a grid whose coordinate system has an EPSG code: "
            f"none found for {raster_path}"
        )
    shim = LayerShim(info["extent"], info["authid"], info["geographic"])
    return shim, int(info["authid"].split(":")[1])


def grid_info_text(info):
    """e.g. "Units: m   CRS: EPSG:28350   Pixels: 120 x 120   Pixel size: 100 m"."""
    if not info:
        return ""
    unit = "deg" if info["geographic"] else "m"
    crs = info["authid"] or info["crs_name"] or "unknown"
    text = f"Units: {unit}   CRS: {crs}"
    w, h = info["width"], info["height"]
    if w > 0 and h > 0:
        dx, dy = info["cell_x"], info["cell_y"]
        size = (f"{dx:.6g}" if abs(dx - dy) <= 1e-9 * max(abs(dx), abs(dy))
                else f"{dx:.6g} x {dy:.6g}")
        text += f"   Pixels: {w} x {h}   Pixel size: {size} {unit}"
    return text


def write_provenance(output_path, source_path, operation, parameters=None):
    """Record how output_path was made, inside the GeoTIFF (see calcs/sgt_metadata).
    Never raises: provenance must not stop a calculation."""
    try:
        from calcs.sgt_metadata import write_sgt_metadata

        write_sgt_metadata(
            output_path, source_path, operation, parameters or {}, plugin_version()
        )
    except Exception as e:
        print(f"SGTool: could not write provenance for {output_path}: {e}")


def step_record(spec_id, **values):
    """(operation, parameters) for a filter in the form the replay understands:
    the filter table's own operation name and parameter names (settings not given
    take their defaults), so a history made here can be replayed on another grid."""
    from sgt_filters import FILTER_SPECS

    spec = next(s for s in FILTER_SPECS if s["id"] == spec_id)
    record = {}
    for p in spec["params"]:
        if p["name"] in values:
            record[p["name"]] = values[p["name"]]
        elif p["kind"] == "enum":
            record[p["name"]] = p["options"][p["default"]]
        else:
            record[p["name"]] = p["default"]
    for key, value in values.items():
        record.setdefault(key, value)  # extra items (such as the survey date)
    return spec["operation"], record


def step_suffix(step):
    """The file-name suffix the QGIS plugin gives the result of a (operation, parameters)
    step from step_record, e.g. "_RTP", "_d1z", "_UC_500", "_BP_50000_5000". None if the
    step is not in the shared filter table (the caller then names it)."""
    if not step:
        return None
    operation, values = step
    from sgt_filters import FILTER_SPECS, SUFFIXES

    for spec in FILTER_SPECS:
        if spec["operation"] == operation:
            try:
                return SUFFIXES[spec["id"]](values)
            except Exception:
                return None
    return None


# ---------------------------------------------------------------------------
# MRVBF / MRRTF (Gallant & Dowling 2003) and the slope and RGB that go with it
# ---------------------------------------------------------------------------
MRVBF_DEFAULTS = dict(t_slope=16.0, t_pctl_v=0.40, t_pctl_r=0.35,
                      p_slope=4.0, p_pctl=3.0, max_res=100.0)


def slope_degrees(arr, cx, cy):
    """Slope of a grid in degrees (3 x 3 Horn method); NaN where the cell or one of
    its neighbours is no data."""
    from scipy.ndimage import convolve

    a = np.asarray(arr, dtype=float)
    nan = ~np.isfinite(a)
    z = np.where(nan, np.nanmean(a) if (~nan).any() else 0.0, a)
    kx = np.array([[1, 0, -1], [2, 0, -2], [1, 0, -1]], dtype=float) / (8.0 * cx)
    ky = np.array([[1, 2, 1], [0, 0, 0], [-1, -2, -1]], dtype=float) / (8.0 * cy)
    slope = np.degrees(np.arctan(np.hypot(convolve(z, kx, mode="nearest"),
                                          convolve(z, ky, mode="nearest"))))
    slope[convolve(nan.astype(float), np.ones((3, 3)), mode="nearest") > 0] = np.nan
    return slope


def mrvbf_set(arr, ll_x, ll_y, cx, **params):
    """(MRVBF, MRRTF, slope in degrees) for a north-up DEM array (row 0 = north), as
    the QGIS plugin makes them. params override MRVBF_DEFAULTS."""
    from calcs.mrvbf import mrvbf as calc_mrvbf

    settings = dict(MRVBF_DEFAULTS, **params)
    dem = np.flipud(np.asarray(arr, dtype=float))  # mrvbf() wants row 0 = ymin
    valid_mask = ~np.isfinite(dem)
    mrvbf_su, mrrtf_su = calc_mrvbf(dem, float(cx), xmin=ll_x, ymin=ll_y,
                                    nodata_mask=valid_mask, **settings)
    return (np.flipud(mrvbf_su), np.flipud(mrrtf_su),
            slope_degrees(arr, cx, cx))


def numpy_to_raster_multiband(bands, output_path, ll_x, ll_y, cell_x, cell_y,
                              spatial_ref=None, nodata=NODATA):
    """Write equal-shaped 2-D arrays as the bands of one Float32 GeoTIFF."""
    from osgeo import gdal

    output_path = str(output_path)
    if os.path.exists(output_path):
        try:
            os.remove(output_path)
        except OSError:
            pass
    out = [np.where(np.isnan(b), nodata, b).astype(np.float32) for b in bands]
    rows, cols = out[0].shape
    ds = gdal.GetDriverByName("GTiff").Create(output_path, cols, rows, len(out),
                                              gdal.GDT_Float32)
    ds.SetGeoTransform((ll_x, cell_x, 0.0, ll_y + rows * cell_y, 0.0, -cell_y))
    if isinstance(spatial_ref, str) and spatial_ref:
        ds.SetProjection(spatial_ref)
    elif spatial_ref is not None:
        ds.SetProjection(spatial_ref.exportToString())
    for i, band_array in enumerate(out, start=1):
        band = ds.GetRasterBand(i)
        band.WriteArray(band_array)
        band.SetNoDataValue(float(nodata))
    ds.FlushCache()
    ds = None


def write_mrvbf_outputs(in_path, base, ext, arr, ll_x, ll_y, cx, cy, sr, nodata, **params):
    """Compute MRVBF / MRRTF / slope for in_path and write <base>_MRVBF, _MRRTF,
    _slope and the 3-band _MRVBF_RGB (R = MRRTF, G = MRVBF, B = slope) like QGIS,
    each with its provenance. Returns the paths written, in that order."""
    settings = dict(MRVBF_DEFAULTS, **params)
    mrvbf_grid, mrrtf_grid, slope = mrvbf_set(arr, ll_x, ll_y, cx, **settings)
    record = {k: settings[k] for k in ("t_slope", "t_pctl_v", "t_pctl_r",
                                       "p_slope", "p_pctl", "max_res")}
    written = []
    for suffix, grid, operation, rec in (
        ("_MRVBF", mrvbf_grid, "Multiresolution valley bottom flatness", record),
        ("_MRRTF", mrrtf_grid, "Multiresolution ridge top flatness", record),
        ("_slope", slope, "Slope (degrees)", {"z_factor": 1.0}),
    ):
        path = f"{base}{suffix}{ext}"
        numpy_to_raster(grid, path, ll_x, ll_y, cx, cy, sr, nodata)
        write_provenance(path, in_path, operation, rec)
        written.append(path)
    rgb = f"{base}_MRVBF_RGB{ext}"
    numpy_to_raster_multiband([mrrtf_grid, mrvbf_grid, slope], rgb, ll_x, ll_y, cx, cy, sr, nodata)
    write_provenance(rgb, in_path, "MRVBF/MRRTF/slope RGB composite",
                     dict(record, bands="R=MRRTF, G=MRVBF, B=slope"))
    written.append(rgb)
    return written


# ---------------------------------------------------------------------------
# Multilevel B-spline gridding of scattered points (native translation of SAGA)
# ---------------------------------------------------------------------------
def list_point_fields(path):
    """Attribute names of a point file readable by OGR (shapefile, GeoPackage, ...)."""
    from osgeo import ogr

    ds = ogr.Open(str(path))
    if ds is None:
        return []
    layer = ds.GetLayer(0)
    definition = layer.GetLayerDefn()
    names = [definition.GetFieldDefn(i).GetName() for i in range(definition.GetFieldCount())]
    ds = None
    return names


def read_points_ogr(path, field):
    """x, y and the numeric value of `field` for every point of an OGR point file, as
    arrays, plus the file's CRS as WKT. Values stored as text (as the SGTool importer
    writes them) are converted; points without a number are skipped."""
    from osgeo import ogr

    ds = ogr.Open(str(path))
    if ds is None:
        raise IOError(f"Cannot open points file: {path}")
    layer = ds.GetLayer(0)
    srs = layer.GetSpatialRef()
    wkt = srs.ExportToWkt() if srs is not None else ""
    xs, ys, zs = [], [], []
    for feat in layer:
        geom = feat.GetGeometryRef()
        if geom is None or geom.IsEmpty():
            continue
        try:
            z = float(feat.GetField(field))
        except (TypeError, ValueError, KeyError):
            continue
        if geom.GetGeometryCount() > 0:  # multipoint
            points = [geom.GetGeometryRef(i) for i in range(geom.GetGeometryCount())]
        else:
            points = [geom]
        for pt in points:
            xs.append(pt.GetX())
            ys.append(pt.GetY())
            zs.append(z)
    ds = None
    return np.array(xs), np.array(ys), np.array(zs), wkt


def read_attributes_ogr(path, fields):
    """The raw values of `fields` (None for a blank name) for every feature of an OGR
    file, as {field: list of values}, in feature order."""
    from osgeo import ogr

    ds = ogr.Open(str(path))
    if ds is None:
        raise IOError(f"Cannot open file: {path}")
    layer = ds.GetLayer(0)
    names = [f for f in fields if f]
    values = {f: [] for f in names}
    for feat in layer:
        for f in names:
            values[f].append(feat.GetField(f) if feat.IsFieldSet(f) else None)
    ds = None
    return values


def read_polygon_categories(path, field, target_wkt):
    """(category text, geometry WKT in the CRS of target_wkt) for every polygon of an OGR
    file with a value in `field`, reprojected if the CRSs differ."""
    from osgeo import ogr, osr

    ds = ogr.Open(str(path))
    if ds is None:
        raise IOError(f"Cannot open polygon file: {path}")
    layer = ds.GetLayer(0)
    src = layer.GetSpatialRef()
    transform = None
    if src is not None and target_wkt:
        dst = osr.SpatialReference()
        dst.ImportFromWkt(target_wkt)
        src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        dst.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        if not src.IsSame(dst):
            transform = osr.CoordinateTransformation(src, dst)
    shapes = []
    for feat in layer:
        value = feat.GetField(field) if feat.IsFieldSet(field) else None
        geom = feat.GetGeometryRef()
        if value is None or not str(value).strip() or geom is None or geom.IsEmpty():
            continue
        geom = geom.Clone()
        if transform is not None:
            geom.Transform(transform)
        shapes.append((str(value), geom.ExportToWkt()))
    ds = None
    return shapes


def suggest_bspline_settings(path, cell_size=None):
    """Cell size and Max Levels suggested for gridding the points in an OGR point file (see
    calcs/mba_suggest.py): the dict from suggest_mba_parameters plus "description"."""
    from osgeo import ogr
    from calcs.mba_suggest import describe, suggest_mba_parameters

    ds = ogr.Open(str(path))
    if ds is None:
        raise IOError(f"Cannot open points file: {path}")
    layer = ds.GetLayer(0)
    srs = layer.GetSpatialRef()
    geographic = bool(srs is not None and srs.IsGeographic())
    xs, ys = [], []
    for feat in layer:
        geom = feat.GetGeometryRef()
        if geom is None or geom.IsEmpty():
            continue
        if geom.GetGeometryCount() > 0:  # multipoint
            points = [geom.GetGeometryRef(i) for i in range(geom.GetGeometryCount())]
        else:
            points = [geom]
        for pt in points:
            xs.append(pt.GetX())
            ys.append(pt.GetY())
    ds = None
    suggestion = suggest_mba_parameters(np.array(xs), np.array(ys), geographic=geographic,
                                        cell_size=cell_size)
    suggestion["description"] = describe(suggestion)
    return suggestion


def grid_points_to_file(xs, ys, zs, out_path, source_path, field, cell_size,
                        epsilon=1e-4, level_max=11, ignore_below=None,
                        spatial_ref=None, callback=None):
    """Multilevel B-spline gridding of points into a GeoTIFF, with the QGIS plugin's
    grid geometry (lower-left node snapped down to a multiple of the cell size) and its
    provenance. Returns the number of points used."""
    from calcs.saga_mba_gridding import mba_gridding

    xs, ys, zs = (np.asarray(v, dtype=float) for v in (xs, ys, zs))
    if ignore_below is not None:
        keep = zs >= ignore_below
        xs, ys, zs = xs[keep], ys[keep], zs[keep]
    if xs.size < 3:
        raise ValueError("Need at least 3 valid points with a numeric value to grid")
    xmin = float(np.floor(xs.min() / cell_size) * cell_size)
    ymin = float(np.floor(ys.min() / cell_size) * cell_size)
    nx = 2 + int((xs.max() - xmin) / cell_size)
    ny = 2 + int((ys.max() - ymin) / cell_size)
    if nx < 2 or ny < 2:
        raise ValueError("Cell size too large for the point extent")
    grid = mba_gridding(xs, ys, zs, cellsize=cell_size, xmin=xmin, ymin=ymin, nx=nx, ny=ny,
                        epsilon=epsilon, level_max=level_max, refinement=False,
                        callback=callback)
    grid = np.flipud(grid)  # the gridder's row 0 is ymin; rasters are written north-up
    # cell (0, 0) is centred on (xmin, ymin): the raster origin is half a cell lower-left
    numpy_to_raster(grid, out_path, xmin - cell_size / 2.0, ymin - cell_size / 2.0,
                    cell_size, cell_size, spatial_ref, NODATA)
    record = {"data_field": field, "cell_size": cell_size, "epsilon": epsilon,
              "max_levels": level_max, "points_used": int(xs.size)}
    if ignore_below is not None:
        record["ignore_values_below"] = ignore_below
    write_provenance(out_path, source_path, "Multilevel B-spline gridding", record)
    return int(xs.size)


# ---------------------------------------------------------------------------
# The ArcGIS Pro map view, as a part of a grid (for the preview window)
# ---------------------------------------------------------------------------
def map_view_to_cells(view, raster_path, clip=True):
    """The extent of the ArcGIS Pro map view as a box of cells in a raster:
    (x0, y0, x1, y1), x = column, y = row (row 0 at the top), clipped to the raster unless
    clip is False (the whole map view, which may reach beyond the raster).

    view is what the toolbox tool writes: {"extent": [xmin, ymin, xmax, ymax], "epsg": code
    or 0, "wkt": text}, in the map's coordinate system. Returns None if the map view does
    not overlap the raster."""
    from osgeo import gdal, osr

    ds = gdal.Open(str(raster_path).split("|")[0], gdal.GA_ReadOnly)
    if ds is None:
        return None
    gt, width, height, wkt = ds.GetGeoTransform(), ds.RasterXSize, ds.RasterYSize, ds.GetProjection()
    ds = None
    xmin, ymin, xmax, ymax = (float(v) for v in view["extent"])
    corners = [(xmin, ymin), (xmin, ymax), (xmax, ymin), (xmax, ymax)]

    dst = osr.SpatialReference()
    if wkt:
        dst.ImportFromWkt(wkt)
    src = osr.SpatialReference()
    epsg = int(view.get("epsg") or 0)
    if epsg:
        src.ImportFromEPSG(epsg)
    elif view.get("wkt"):
        src.ImportFromWkt(view["wkt"])
        try:
            src.MorphFromESRI()
        except Exception:
            logging.getLogger(__name__).debug("Ignored: ignored", exc_info=True)
    if (epsg or view.get("wkt")) and wkt and not src.IsSame(dst):
        for ref in (src, dst):
            if hasattr(osr, "OAMS_TRADITIONAL_GIS_ORDER"):  # x = easting / longitude
                ref.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        transform = osr.CoordinateTransformation(src, dst)
        corners = [transform.TransformPoint(x, y)[:2] for x, y in corners]

    cols = [(x - gt[0]) / gt[1] for x, _y in corners]
    rows = [(y - gt[3]) / gt[5] for _x, y in corners]
    x0, x1 = max(0.0, min(cols)), min(float(width), max(cols))
    y0, y1 = max(0.0, min(rows)), min(float(height), max(rows))
    if x1 - x0 < 1.0 or y1 - y0 < 1.0:
        return None  # the map is not looking at this raster
    if not clip:
        return (min(cols), min(rows), max(cols), max(rows))
    return (x0, y0, x1, y1)
