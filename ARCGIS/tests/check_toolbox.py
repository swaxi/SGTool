"""
Checks of the ArcGIS Pro toolbox (SGTool.pyt): every tool builds its parameters, and
some are run.

Run it with ArcGIS Pro's own Python (it needs arcpy-era GDAL, numpy, scipy and tkinter):

    "C:\\Program Files\\ArcGIS\\Pro\\bin\\Python\\envs\\arcgispro-py3\\python.exe" ARCGIS\\tests\\check_toolbox.py

It builds small synthetic grids in a temporary folder, drives the code, prints PASS / FAIL per
check and a total, and exits with status 1 if anything failed. It needs arcpy (so ArcGIS Pro's Python), but not a running copy of ArcGIS Pro.
"""
import os, sys, tempfile, importlib.machinery, importlib.util, shutil
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = tempfile.mkdtemp(prefix="sgtool_arc_toolbox_")
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)
from osgeo import gdal, osr

results = []


def check(n, c, d=""):
    results.append((n, bool(c)))
    print(("PASS  " if c else "FAIL  ") + n + (f"   [{d}]" if d and not c else ""), flush=True)


rng = np.random.default_rng(2)
arr = (rng.normal(0, 5, (120, 120)) + np.repeat(rng.normal(0, 20, 30), 4)[:, None]).astype("float32")
src = os.path.join(DATA, "t.tif")
ds = gdal.GetDriverByName("GTiff").Create(src, 120, 120, 1, gdal.GDT_Float32)
ds.SetGeoTransform([500000, 100, 0, 6900000, 0, -100])
s = osr.SpatialReference()
s.ImportFromEPSG(28350)
ds.SetProjection(s.ExportToWkt())
b = ds.GetRasterBand(1)
b.SetNoDataValue(-9999.0)
b.WriteArray(arr)
ds = None

loader = importlib.machinery.SourceFileLoader("sgt_pyt", os.path.join(REPO, "ARCGIS", "SGTool.pyt"))
spec = importlib.util.spec_from_loader("sgt_pyt", loader)
mod = importlib.util.module_from_spec(spec)
loader.exec_module(mod)
tb = mod.Toolbox()
check("toolbox lists the new line-noise tool", mod.RemoveLineNoise in tb.tools)
broken = []
for cls_ in tb.tools:
    try:
        cls_().getParameterInfo()
    except Exception as e_:
        broken.append((cls_.__name__, str(e_)[:60]))
check("every toolbox tool builds its parameter list (%d tools)" % len(tb.tools), not broken, broken)


class Msg:
    def __init__(self):
        self.lines = []

    def addMessage(self, m):
        self.lines.append(m)

    addWarningMessage = addErrorMessage = addMessage


def run(tool_cls, values):
    tool = tool_cls()
    params = tool.getParameterInfo()
    for p in params:
        if p.name in values:
            p.value = values[p.name]
    msgs = Msg()
    tool.execute(params, msgs)
    return msgs


from calcs import sgt_metadata as sm

out1 = os.path.join(DATA, "t_noise.tif")
m = run(mod.RemoveLineNoise, {"in_raster": src, "out_raster": out1, "azimuth": 90.0,
                              "line_spacing_min": 400.0, "scale": 1.0, "wedge": 40.0})
check("line-noise tool wrote its output", os.path.exists(out1), m.lines)
root = sm.read_sgt_metadata(out1)
pr = {p.get("name"): p.text for p in root.iter("parameter")} if root is not None else {}
check("line-noise tool records replayable settings",
      root is not None and root.findtext("operation") == "Directional Cosine/Butterworth line-noise removal"
      and pr.get("WEDGE") in ("40.0", "40") and pr.get("AZIMUTH") in ("90.0", "90"), pr)
out2 = os.path.join(DATA, "t_mean.tif")
run(mod.MeanFilter, {"in_raster": src, "out_raster": out2})
r2 = sm.read_sgt_metadata(out2)
check("any toolbox tool now records provenance",
      r2 is not None and r2.findtext("operation") == mod.MeanFilter.label,
      None if r2 is None else r2.findtext("operation"))
check("the source is recorded",
      r2 is not None and os.path.basename(r2.findtext("sources/source/path") or "") == "t.tif")
yy, xx = np.mgrid[0:120, 0:120]
dem = os.path.join(DATA, "dem.tif")
z = (np.minimum(np.abs(xx - 60) * 6.0, 120.0) + 0.05 * yy).astype("float32")
ds = gdal.GetDriverByName("GTiff").Create(dem, 120, 120, 1, gdal.GDT_Float32)
ds.SetGeoTransform([500000, 25, 0, 6900000, 0, -25]); s2 = osr.SpatialReference(); s2.ImportFromEPSG(28350); ds.SetProjection(s2.ExportToWkt())
b = ds.GetRasterBand(1); b.SetNoDataValue(-9999.0); b.WriteArray(z); ds = None
outm = os.path.join(DATA, "dem_MRVBF.tif")
m = run(mod.MRVBFTool, {"in_raster": dem, "out_raster": outm})
check("the MRVBF tool writes MRVBF, MRRTF, slope and the RGB composite",
      all(os.path.exists(os.path.join(DATA, n)) for n in ("dem_MRVBF.tif", "dem_MRRTF.tif", "dem_slope.tif", "dem_MRVBF_RGB.tif")), m.lines)
check("it is listed in the toolbox", mod.MRVBFTool in tb.tools)
sys.stdout.flush()
from osgeo import ogr
shp = os.path.join(DATA, "bsp.shp")
dsv = ogr.GetDriverByName("ESRI Shapefile").CreateDataSource(shp)
sr_ = osr.SpatialReference(); sr_.ImportFromEPSG(28350)
lay = dsv.CreateLayer("bsp", sr_, ogr.wkbPoint); lay.CreateField(ogr.FieldDefn("num", ogr.OFTReal))
rg = np.random.default_rng(11)
for x, y in zip(500000 + rg.uniform(0, 2000, 400), 6890000 + rg.uniform(0, 2000, 400)):
    f = ogr.Feature(lay.GetLayerDefn()); g = ogr.Geometry(ogr.wkbPoint); g.AddPoint(float(x), float(y))
    f.SetGeometry(g); f.SetField("num", 0.02 * (x - 500000) + 0.01 * (y - 6890000)); lay.CreateFeature(f)
dsv = None
outb = os.path.join(DATA, "bsp_num_bspline.tif")
mm = run(mod.MultilevelBSplineGridding, {"in_points": shp, "z_field": "num", "cell_size": 50.0, "out_raster": outb})
rb_ = sm.read_sgt_metadata(outb)
check("the B-spline gridding tool writes a grid with provenance",
      os.path.exists(outb) and rb_ is not None and rb_.findtext("operation") == "Multilevel B-spline gridding", mm.lines)
dd = gdal.Open(outb); aa = dd.GetRasterBand(1).ReadAsArray().astype(float); gg = dd.GetGeoTransform(); dd = None
cx_ = gg[0] + (np.arange(aa.shape[1]) + 0.5) * gg[1]; cy_ = gg[3] + (np.arange(aa.shape[0]) + 0.5) * gg[5]
XX, YY = np.meshgrid(cx_, cy_); inner = (XX > 500300) & (XX < 501700) & (YY > 6890300) & (YY < 6891700)
check("the toolbox gridding follows the plane too", float(np.nanmax(np.abs((aa - (0.02 * (XX - 500000) + 0.01 * (YY - 6890000)))[inner]))) < 1.0)
check("it is listed in the toolbox", mod.MultilevelBSplineGridding in tb.tools)
outa = os.path.join(DATA, "bsp_auto_bspline.tif")
ma = run(mod.MultilevelBSplineGridding, {"in_points": shp, "z_field": "num", "cell_size": 9999.0, "max_levels": 2,
                                         "auto": True, "out_raster": outa})
from calcs.mba_suggest import suggest_mba_parameters as _suggest
_ds = ogr.Open(shp); _l = _ds.GetLayer(0)
_pts = np.array([(f.GetGeometryRef().GetX(), f.GetGeometryRef().GetY()) for f in _l]); _ds = None
_want = _suggest(_pts[:, 0], _pts[:, 1])
ra = sm.read_sgt_metadata(outa)
pa = {p_.get("name"): p_.text for p_ in ra.iter("parameter")} if ra is not None else {}
check("Auto in the B-spline tool ignores the cell size and levels given",
      os.path.exists(outa) and abs(float(pa.get("cell_size", 0)) - _want["cell_size"]) < 1e-9
      and pa.get("max_levels") == str(_want["max_levels"]), (pa, _want["cell_size"], _want["max_levels"], ma.lines))
# ---- the Cancel button of the geoprocessing pane stops the tools (arcpy.env.isCancelled is polled)
from calcs.sgt_cancel import OperationCancelled
_real_check = mod._check_cancel


def cancel_after(n):
    state = {"calls": 0}

    def f(fraction=None):
        state["calls"] += 1
        if state["calls"] > n:
            raise OperationCancelled()
    f.state = state
    return f


try:
    outc = os.path.join(DATA, "cancel_mean.tif")
    mod._check_cancel = cancel_after(0)
    mc = run(mod.MeanFilter, {"in_raster": src, "out_raster": outc})
    check("a cancelled filter writes nothing and says so",
          not os.path.exists(outc) and any("Cancelled" in l for l in mc.lines), mc.lines)

    # cancelled after the calculation, just before the result is saved: still nothing is written
    outc2 = os.path.join(DATA, "cancel_mean2.tif")
    mod._check_cancel = cancel_after(1)   # let _read through, stop at _write
    mc = run(mod.MeanFilter, {"in_raster": src, "out_raster": outc2})
    check("a result finished after Cancel is thrown away, not saved",
          not os.path.exists(outc2) and any("Cancelled" in l for l in mc.lines), mc.lines)

    outc3 = os.path.join(DATA, "cancel_bs.tif")
    mod._check_cancel = cancel_after(2)   # stops inside the gridding levels
    mc = run(mod.MultilevelBSplineGridding, {"in_points": shp, "z_field": "num", "cell_size": 50.0,
                                              "out_raster": outc3})
    check("B-spline gridding stops between levels when cancelled",
          not os.path.exists(outc3) and any("Cancelled" in l for l in mc.lines), mc.lines)

    outc4 = os.path.join(DATA, "cancel_mrvbf.tif")
    fn = cancel_after(3)
    mod._check_cancel = fn
    mc = run(mod.MRVBFTool, {"in_raster": dem, "out_raster": outc4})
    check("MRVBF stops part-way when cancelled (after %d polls)" % fn.state["calls"],
          not os.path.exists(outc4) and any("Cancelled" in l for l in mc.lines), mc.lines)
finally:
    mod._check_cancel = _real_check

# without a cancel the same tools still run (the poll is quiet outside a cancelled tool)
outc5 = os.path.join(DATA, "not_cancelled.tif")
mc = run(mod.MeanFilter, {"in_raster": src, "out_raster": outc5})
check("an uncancelled tool is not disturbed by the polling", os.path.exists(outc5))

# the calculation code reports progress and can be stopped from its callback
from calcs.mrvbf import mrvbf as _mrvbf
_dem = np.random.default_rng(4).normal(0, 1, (80, 90)).cumsum(0).cumsum(1)
_seen = []
_mrvbf(_dem, 25.0, callback=_seen.append)
check("mrvbf reports progress between 0 and 1, rising", len(_seen) >= 4 and all(0 <= f <= 1 for f in _seen)
      and _seen == sorted(_seen), _seen)


def _stop_at_once(fraction=None):
    raise OperationCancelled()


try:
    _mrvbf(_dem, 25.0, callback=_stop_at_once)
    check("mrvbf stops when its callback raises", False)
except OperationCancelled:
    check("mrvbf stops when its callback raises", True)

print(f"\n{sum(1 for _, ok in results if ok)}/{len(results)} PYT checks passed", flush=True)
os._exit(0 if all(ok for _, ok in results) else 1)
