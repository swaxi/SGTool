# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2024-2026 Mark Jessell
"""
Grid and vector files for the pysgtool package, with GDAL (osgeo) where it is
installed, as in QGIS, ArcGIS Pro and conda, or else with rasterio and pyogrio, which
pip can install everywhere (pip install pysgtool[io]).

Rasters are (bands, rows, cols) float arrays, north-up, NaN where there is no data,
with a GDAL geotransform and a WKT coordinate system ('' for none).

No QGIS or ArcGIS dependency. The plugins themselves always have GDAL and read and
write their grids as they did; nothing here is used there.
"""

import os

import numpy as np

_backend = None


def backend():
    """"gdal" or "rasterio": what reads and writes the files."""
    global _backend
    if _backend is None:
        wanted = os.environ.get("PYSGTOOL_IO", "").strip().lower()
        for name in ([wanted] if wanted else ["gdal", "rasterio"]):
            if _available(name):
                _backend = name
                break
        else:
            raise ImportError(
                "Reading and writing files needs GDAL or rasterio: pip install pysgtool[io]"
            )
    return _backend


def set_backend(name):
    """Use "gdal" or "rasterio" (None: choose again), instead of GDAL if it is installed."""
    global _backend
    if name is not None and not _available(name):
        raise ImportError(f"{name} is not installed")
    _backend = name


def _available(name):
    if name not in ("gdal", "rasterio"):
        raise ValueError(f"The file backend is 'gdal' or 'rasterio', not {name!r}")
    try:
        if name == "gdal":
            from osgeo import gdal  # noqa: F401
        else:
            import rasterio  # noqa: F401
        return True
    except ImportError:
        return False


def _gdal():
    from osgeo import gdal

    return gdal


# ---------------------------------------------------------------- rasters
def _north_up(stack, gt):
    if gt[5] > 0:  # south-up: flip to north-up
        stack = stack[:, ::-1, :]
        gt = (gt[0], gt[1], 0.0, gt[3] + stack.shape[1] * gt[5], 0.0, -gt[5])
    if gt[2] or gt[4]:
        raise ValueError("Rotated grids are not supported: resample the grid to north-up first")
    return np.array(stack), tuple(float(v) for v in gt)


def read(path, bands=None, masks=False):
    """Read a raster: a dict with "stack", "gt", "wkt", "names" (band descriptions) and
    "metadata" (text items). bands: 1-based band numbers (default all). masks: also
    leave out cells that an alpha band or mask marks as invalid."""
    path = str(path)
    if backend() == "gdal":
        gdal = _gdal()
        try:
            ds = gdal.Open(path)
        except RuntimeError as e:
            raise OSError(f"Could not open the raster {path}: {e}")
        if ds is None:
            raise OSError(f"Could not open the raster {path}")
        count, gt, wkt = ds.RasterCount, ds.GetGeoTransform(), ds.GetProjection()
        numbers = _bands(path, bands, count)
        stack, names = [], []
        for b in numbers:
            band = ds.GetRasterBand(b)
            arr = band.ReadAsArray().astype(float)
            nodata = band.GetNoDataValue()
            if nodata is not None:
                arr[arr == nodata] = np.nan
            if masks and band.GetMaskFlags() != gdal.GMF_ALL_VALID:
                arr[band.GetMaskBand().ReadAsArray() == 0] = np.nan
            stack.append(arr)
            names.append(band.GetDescription() or "")
        metadata = dict(ds.GetMetadata() or {})
        ds = band = None
    else:
        import rasterio

        try:
            ds = rasterio.open(path)
        except rasterio.errors.RasterioIOError as e:
            raise OSError(f"Could not open the raster {path}: {e}")
        with ds:
            numbers = _bands(path, bands, ds.count)
            gt, wkt = ds.transform.to_gdal(), _wkt_of(ds.crs)
            stack, names = [], []
            for b in numbers:
                arr = ds.read(b).astype(float)
                nodata = ds.nodatavals[b - 1]
                if nodata is not None:
                    arr[arr == nodata] = np.nan
                if masks:
                    arr[ds.read_masks(b) == 0] = np.nan
                stack.append(arr)
                names.append(ds.descriptions[b - 1] or "")
            metadata = dict(ds.tags())
    stack, gt = _north_up(np.stack(stack), gt)
    return {"stack": stack, "gt": gt, "wkt": wkt, "names": names, "metadata": metadata}


def _bands(path, bands, count):
    numbers = list(range(1, count + 1)) if bands is None else list(bands)
    for b in numbers:
        if not 1 <= b <= count:
            raise ValueError(f"{path} has no band {b}")
    return numbers


def read_grid(path):
    """Band 1 of a raster, its geotransform and WKT (as filter_table.read_grid)."""
    r = read(path, [1])
    return r["stack"][0], r["gt"], r["wkt"]


def write(path, stack, gt, wkt, names=None, metadata=None):
    """Write a (bands, rows, cols) stack as a float32 GeoTIFF with NaN no data."""
    path = str(path)
    stack = np.asarray(stack)
    n, rows, cols = stack.shape
    if backend() == "gdal":
        gdal = _gdal()
        options = ["COMPRESS=DEFLATE", "TILED=YES"] if n > 1 else []
        ds = gdal.GetDriverByName("GTiff").Create(path, cols, rows, n, gdal.GDT_Float32,
                                                  options=options)
        if ds is None:
            raise OSError(f"Could not create {path}")
        ds.SetGeoTransform(gt)
        if wkt:
            ds.SetProjection(wkt)
        for k, v in (metadata or {}).items():
            ds.SetMetadataItem(k, str(v))
        for i in range(n):
            band = ds.GetRasterBand(i + 1)
            band.SetNoDataValue(np.nan)
            if names and names[i]:
                band.SetDescription(names[i])
            band.WriteArray(stack[i].astype(np.float32))
        ds.FlushCache()
        band = ds = None
        return
    import rasterio
    from rasterio.transform import Affine

    options = {"compress": "deflate", "tiled": True} if n > 1 else {}
    with rasterio.open(
        path, "w", driver="GTiff", height=rows, width=cols, count=n, dtype="float32",
        crs=_rio_crs(wkt), transform=Affine.from_gdal(*gt), nodata=np.nan, **options,
    ) as ds:
        if metadata:
            ds.update_tags(**{k: str(v) for k, v in metadata.items()})
        for i in range(n):
            if names and names[i]:
                ds.set_band_description(i + 1, names[i])
            ds.write(stack[i].astype(np.float32), i + 1)


def resample(stack, gt, wkt, to_gt, to_shape, to_wkt, how=None):
    """A stack resampled (reprojected if need be) onto another north-up frame. how:
    "nearest", "bilinear" or "average"; by default average where the new frame is
    coarser, bilinear otherwise, as the plugin does."""
    stack = np.asarray(stack, dtype=float)
    rows, cols = to_shape
    if (stack.shape[1:] == (rows, cols)
            and np.allclose(gt, to_gt, rtol=0, atol=1e-9 * max(1.0, abs(to_gt[1])))
            and (not wkt or not to_wkt or same_crs(wkt, to_wkt))):
        return stack.copy()
    how = how or ("average" if abs(to_gt[1]) > 1.5 * abs(gt[1]) else "bilinear")
    reproject = bool(wkt and to_wkt)
    if backend() == "gdal":
        gdal = _gdal()
        src = gdal.GetDriverByName("MEM").Create("", stack.shape[2], stack.shape[1],
                                                 stack.shape[0], gdal.GDT_Float32)
        src.SetGeoTransform(gt)
        if wkt:
            src.SetProjection(wkt)
        for i in range(stack.shape[0]):
            band = src.GetRasterBand(i + 1)
            band.SetNoDataValue(np.nan)
            band.WriteArray(stack[i].astype(np.float32))
        kwargs = {"srcSRS": wkt, "dstSRS": to_wkt} if reproject else {}
        out = gdal.Warp(
            "", src, format="MEM", width=cols, height=rows,
            outputBounds=(to_gt[0], to_gt[3] + rows * to_gt[5], to_gt[0] + cols * to_gt[1], to_gt[3]),
            resampleAlg={"nearest": "near"}.get(how, how),
            outputType=gdal.GDT_Float32, dstNodata=np.nan, **kwargs,
        )
        if out is None:
            raise ValueError("Could not resample the grid")
        result = np.stack([out.GetRasterBand(i + 1).ReadAsArray() for i in range(stack.shape[0])])
        out = src = band = None
        return result.astype(float)
    from rasterio.enums import Resampling
    from rasterio.transform import Affine
    from rasterio.warp import reproject as warp

    out = np.full((stack.shape[0], rows, cols), np.nan, dtype=np.float32)
    # without both coordinate systems the two frames are taken to share one
    src_crs, dst_crs = (_rio_crs(wkt), _rio_crs(to_wkt)) if reproject else ("EPSG:3857",) * 2
    for i in range(stack.shape[0]):
        warp(
            stack[i].astype(np.float32), out[i], src_transform=Affine.from_gdal(*gt),
            src_crs=src_crs, dst_transform=Affine.from_gdal(*to_gt), dst_crs=dst_crs,
            src_nodata=np.nan, dst_nodata=np.nan, resampling=Resampling[how],
        )
    return out.astype(float)


def get_tag(path, domain, item):
    """A text item of a metadata domain of a raster, or None."""
    if backend() == "gdal":
        ds = _gdal().Open(str(path))
        text = ds.GetMetadataItem(item, domain) if ds is not None else None
        ds = None
        return text
    import rasterio

    with rasterio.open(str(path)) as ds:
        return ds.tags(ns=domain).get(item)


def set_tags(path, domains):
    """Store {domain: {item: text}} in an existing GeoTIFF ("" is the default domain)."""
    if backend() == "gdal":
        gdal = _gdal()
        ds = gdal.Open(str(path), gdal.GA_Update)
        if ds is None:
            raise OSError(f"Could not update {path}")
        for domain, items in domains.items():
            for k, v in items.items():
                ds.SetMetadataItem(k, v, domain)
        ds.FlushCache()
        ds = None
        return
    import rasterio

    with rasterio.open(str(path), "r+") as ds:
        for domain, items in domains.items():
            ds.update_tags(**(dict(items, ns=domain) if domain else items))


def is_geotiff(path):
    """True if the file is a GeoTIFF that can be opened."""
    if backend() == "gdal":
        ds = _gdal().Open(str(path))
        ok = ds is not None and ds.GetDriver().ShortName == "GTiff"
        ds = None
        return ok
    import rasterio

    with rasterio.open(str(path)) as ds:
        return ds.driver == "GTiff"


def rasterize(polygons, gt, shape, wkt=""):
    """The number of the polygon covering each cell centre (-1 outside them all, the
    later polygon where they overlap). polygons: (geometry WKT, number) pairs."""
    if backend() == "gdal":
        from .Petrophysics import rasterize_categories

        return rasterize_categories(polygons, wkt, gt, shape)
    from rasterio.features import rasterize as burn
    from rasterio.transform import Affine
    from shapely import wkt as shapely_wkt

    shapes = [(shapely_wkt.loads(g), int(i)) for g, i in polygons]
    if not shapes:
        return np.full(shape, -1, dtype=np.int32)
    return burn(shapes, out_shape=shape, fill=-1, transform=Affine.from_gdal(*gt), dtype="int32")


# ---------------------------------------------------------------- coordinate systems
def _wkt_of(crs):
    return crs.to_wkt() if crs else ""


def _rio_crs(wkt):
    from rasterio.crs import CRS

    return CRS.from_wkt(wkt) if wkt else None


def _osr(wkt):
    from osgeo import osr

    srs = osr.SpatialReference()
    srs.ImportFromWkt(wkt)
    if hasattr(osr, "OAMS_TRADITIONAL_GIS_ORDER"):
        srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    return srs


def _pyproj_crs(wkt):
    from pyproj import CRS

    return CRS.from_user_input(wkt)


def crs_to_wkt(crs):
    """WKT of a coordinate system given as an EPSG number, 'EPSG:nnnn', WKT or an object
    that has to_wkt() (rasterio, pyproj); '' for None."""
    if crs is None or crs == "":
        return ""
    if hasattr(crs, "to_wkt"):
        return crs.to_wkt()
    text = f"EPSG:{int(crs)}" if isinstance(crs, (int, np.integer)) else str(crs).strip()
    if not text.upper().startswith("EPSG:"):
        return text
    code = int(text.split(":")[1])
    try:
        from osgeo import osr

        srs = osr.SpatialReference()
        if srs.ImportFromEPSG(code) != 0:
            raise ValueError(f"Unknown coordinate system {text}")
        return srs.ExportToWkt()
    except ImportError:
        pass
    except RuntimeError as e:
        raise ValueError(f"Unknown coordinate system {text}: {e}")
    try:
        from pyproj import CRS
        from pyproj.enums import WktVersion
        from pyproj.exceptions import CRSError
    except ImportError:
        return text  # kept as the code; what needs more says it needs pyproj
    try:
        return CRS.from_epsg(code).to_wkt(WktVersion.WKT1_GDAL)
    except CRSError as e:
        raise ValueError(f"Unknown coordinate system {text}: {e}")


def epsg(wkt):
    """EPSG number of a coordinate system, or None."""
    if not wkt:
        return None
    if wkt.upper().startswith("EPSG:"):
        return int(wkt.split(":")[1])
    try:
        from osgeo import osr  # noqa: F401

        srs = _osr(wkt)
        code = srs.GetAuthorityCode(None)
        if code is None:
            try:
                srs.AutoIdentifyEPSG()
            except RuntimeError:
                return None
            code = srs.GetAuthorityCode(None)
        return int(code) if code and code.isdigit() else None
    except ImportError:
        return _pyproj_crs(wkt).to_epsg()


def is_geographic(wkt):
    """True if the coordinate system is in degrees (False if there is none)."""
    if not wkt:
        return False
    try:
        from osgeo import osr  # noqa: F401

        return bool(_osr(wkt).IsGeographic())
    except ImportError:
        return bool(_pyproj_crs(wkt).is_geographic)


def same_crs(a, b):
    """True if two WKT coordinate systems are the same."""
    if a == b:
        return True
    try:
        from osgeo import osr  # noqa: F401

        return bool(_osr(a).IsSame(_osr(b)))
    except ImportError:
        return _pyproj_crs(a).equals(_pyproj_crs(b), ignore_axis_order=True)


# ---------------------------------------------------------------- vectors
def read_vector(path, layer=None):
    """Read a vector file: a dict with "wkb" (geometries), "fields" {name: list} and
    "wkt" (coordinate system)."""
    path = str(path)
    if backend() == "gdal":
        from osgeo import ogr

        try:
            ds = ogr.Open(path)
        except RuntimeError as e:
            raise OSError(f"Could not open {path}: {e}")
        if ds is None:
            raise OSError(f"Could not open {path}")
        lyr = ds.GetLayerByName(layer) if layer else ds.GetLayer(0)
        if lyr is None:
            raise ValueError(f"{path} has no layer {layer}")
        defn = lyr.GetLayerDefn()
        names = [defn.GetFieldDefn(i).GetName() for i in range(defn.GetFieldCount())]
        fields, wkb = {n: [] for n in names}, []
        srs = lyr.GetSpatialRef()
        wkt = srs.ExportToWkt() if srs is not None else ""
        for feature in lyr:
            geom = feature.GetGeometryRef()
            wkb.append(bytes(geom.ExportToWkb()) if geom is not None else None)
            for i, n in enumerate(names):
                fields[n].append(feature.GetField(i))
        feature = geom = srs = defn = lyr = ds = None
        return {"wkb": wkb, "fields": fields, "wkt": wkt}
    from pyogrio.raw import read as read_raw

    try:
        meta, _fids, geometry, data = read_raw(path, layer=layer)
    except Exception as e:
        raise OSError(f"Could not open {path}: {e}")
    names = [str(n) for n in meta["fields"]]
    wkt = crs_to_wkt(meta["crs"]) if meta.get("crs") else ""
    return {"wkb": list(geometry), "wkt": wkt,
            "fields": {n: list(data[i]) for i, n in enumerate(names)}}


_DRIVERS = {".shp": "ESRI Shapefile", ".gpkg": "GPKG", ".geojson": "GeoJSON"}


def write_vector(path, wkb, geometry_type, fields, wkt):
    """Write features to a .shp, .gpkg or .geojson file. geometry_type: 'Point',
    'LineString Z', 'MultiPolygon' ...; fields: {name: values}, numbers or text."""
    path = str(path)
    ext = os.path.splitext(path)[1].lower()
    if ext not in _DRIVERS:
        raise ValueError(f"Vector files are written as .shp, .gpkg or .geojson, not {ext!r}")
    columns = {str(n): np.asarray(v) for n, v in fields.items()}
    for n, v in columns.items():
        if v.dtype.kind in "biu":
            columns[n] = v.astype(np.int64)
        elif v.dtype.kind == "f":
            columns[n] = v.astype(float)
        else:
            columns[n] = np.array(["" if t is None else str(t) for t in v], dtype=object)
    if backend() == "gdal":
        from osgeo import ogr

        types = {"Point": ogr.wkbPoint, "Point Z": ogr.wkbPoint25D,
                 "LineString Z": ogr.wkbLineString25D, "MultiPolygon": ogr.wkbMultiPolygon}
        driver = ogr.GetDriverByName(_DRIVERS[ext])
        if os.path.exists(path):
            driver.DeleteDataSource(path)
        ds = driver.CreateDataSource(path)
        if ds is None:
            raise OSError(f"Could not create {path}")
        name = os.path.splitext(os.path.basename(path))[0]
        lyr = ds.CreateLayer(name, _osr(wkt) if wkt else None, types[geometry_type])
        kinds = {"i": ogr.OFTInteger64, "f": ogr.OFTReal, "O": ogr.OFTString}
        for n, v in columns.items():
            lyr.CreateField(ogr.FieldDefn(n, kinds[v.dtype.kind]))
        values = [v.tolist() for v in columns.values()]
        for k, geometry in enumerate(wkb):
            feature = ogr.Feature(lyr.GetLayerDefn())
            for i, v in enumerate(values):
                feature.SetField(i, v[k])
            feature.SetGeometry(ogr.CreateGeometryFromWkb(geometry))
            lyr.CreateFeature(feature)
            feature = None
        lyr = ds = None
        return
    from pyogrio.raw import write as write_raw

    if os.path.exists(path):
        os.remove(path)
    write_raw(
        path, np.array(list(wkb), dtype=object), list(columns.values()), list(columns),
        driver=_DRIVERS[ext], geometry_type=geometry_type, crs=wkt or None,
    )


def cell_area(gt, shape, wkt, to_wkt):
    """Cell area of a raster once reprojected to another coordinate system."""
    rows, cols = shape
    if backend() == "gdal":
        gdal = _gdal()
        src = gdal.GetDriverByName("MEM").Create("", cols, rows, 1, gdal.GDT_Byte)
        src.SetGeoTransform(gt)
        src.SetProjection(wkt)
        out = gdal.Warp("", src, format="VRT", dstSRS=to_wkt).GetGeoTransform()
        src = None
        return abs(out[1] * out[5])
    from rasterio.warp import calculate_default_transform

    bounds = (gt[0], gt[3] + rows * gt[5], gt[0] + cols * gt[1], gt[3])
    transform, _w, _h = calculate_default_transform(_rio_crs(wkt), _rio_crs(to_wkt), cols, rows,
                                                    *bounds)
    return abs(transform.a * transform.e)
