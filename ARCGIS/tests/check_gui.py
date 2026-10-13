"""
Checks of the ArcGIS Pro GUI (SGTool_GUI.py), with the real tkinter main loop.

Run it with ArcGIS Pro's own Python (it needs arcpy-era GDAL, numpy, scipy and tkinter):

    "C:\\Program Files\\ArcGIS\\Pro\\bin\\Python\\envs\\arcgispro-py3\\python.exe" ARCGIS\\tests\\check_gui.py

It builds small synthetic grids in a temporary folder, drives the code, prints PASS / FAIL per
check and a total, and exits with status 1 if anything failed. The preview window checks open real (briefly visible) windows.
"""
import os, sys, tempfile, time, shutil, traceback, threading, datetime
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = tempfile.mkdtemp(prefix="sgtool_arc_gui_")
sys.path.insert(0, os.path.join(REPO, "ARCGIS"))
sys.path.insert(0, REPO)

from osgeo import gdal, osr

results = []


def check(name, cond, detail=""):
    results.append((name, bool(cond)))
    print(("PASS  " if cond else "FAIL  ") + name + (f"   [{detail}]" if detail and not cond else ""), flush=True)


def make_tif(name, size=100, origin=(500000, 6900000), cell=100.0, seed=1):
    path = os.path.join(DATA, name)
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:size, 0:size]
    arr = (100 * np.exp(-(((xx - 30) ** 2 + (yy - 40) ** 2) / 150.0)) + 0.2 * xx
           + rng.normal(0, 2, (size, size))).astype("float32")
    ds = gdal.GetDriverByName("GTiff").Create(path, size, size, 1, gdal.GDT_Float32)
    ds.SetGeoTransform([origin[0], cell, 0, origin[1], 0, -cell])
    srs = osr.SpatialReference(); srs.ImportFromEPSG(28350)
    ds.SetProjection(srs.ExportToWkt())
    b = ds.GetRasterBand(1); b.SetNoDataValue(-9999.0); b.WriteArray(arr)
    ds = None
    return path


def run_loop(cond, timeout):
    """Run the real Tk main loop (needed for callbacks from job threads) until cond()."""
    end = time.time() + timeout
    def tick():
        if cond() or time.time() > end:
            app.quit()
        else:
            app.after(50, tick)
    app.after(50, tick)
    app.mainloop()


def pump(sec=0.2):
    end = time.time() + sec
    while time.time() < end:
        app.update()
        time.sleep(0.01)


import SGTool_GUI as G
import arcgis_utils as au
read_points_ogr_ = au.read_points_ogr
app = G.SGToolApp()
app.withdraw()
pump(0.5)

src = make_tif("arcsrc.tif")
app.v_input.set(src)
pump(0.2)

# ---- grid information row
info = app.v_grid_info.get()
check("version and grid info shown", "Version" in info and "Units: m" in info and "CRS: EPSG:28350" in info
      and "Pixels: 100 x 100" in info and "Pixel size: 100 m" in info, info)

# ---- build every tab (they are built lazily)
for i in range(5):
    app.nb.select(i)
    app.update()
check("all tabs build", len(app._tab_built) == 5, app._tab_built)
check("Wedge entry exists with default 45", app.v_dc_wedge.get() == "45")

# ---- FFT filters: outputs carry provenance in replayable form
app.v_deriv.set(True); app.v_deriv_dir.set("z") if hasattr(app, "v_deriv_dir") else None
app.v_rtp.set(True); app.v_rtp_type.set("Pole"); app.v_inc.set("-60"); app.v_dec.set("0")
app.v_dirclean.set(True); app.v_dc_az.set("0"); app.v_dc_wl.set("500"); app.v_dc_wl_max.set(""); app.v_dc_wedge.set("40")
app._compute_output_path()
preview = app.v_output_preview.get()
app._do_apply_fft(src)
pump(0.3)
written = sorted(f for f in os.listdir(DATA) if f.startswith("arcsrc_") and f.endswith(".tif"))
check("outputs use the QGIS naming convention",
      {"arcsrc_RTP.tif", "arcsrc_DirC.tif", "arcsrc_d1z.tif"} <= set(written), written)
check("the Output (auto) line shows exactly the names written",
      all(name in preview for name in ("arcsrc_RTP.tif", "arcsrc_DirC.tif", "arcsrc_d1z.tif")), preview)
from calcs import sgt_metadata as sm
outs = sorted(f for f in os.listdir(DATA) if f.startswith("arcsrc_") and f.endswith(".tif"))
check("FFT filters wrote their outputs", {"arcsrc_RTP.tif", "arcsrc_DirC.tif"} <= set(outs) and any("dz" in o or "d" in o for o in outs), outs)

def params(path):
    return {p.get("name"): p.text for p in sm.read_sgt_metadata(path).iter("parameter")}

root = sm.read_sgt_metadata(os.path.join(DATA, "arcsrc_RTP.tif"))
check("provenance is embedded in the GeoTIFF", root is not None and root.findtext("operation") == "Reduction to pole",
      None if root is None else root.findtext("operation"))
check("RTP parameters recorded under the replay's names",
      params(os.path.join(DATA, "arcsrc_RTP.tif")).get("INCLINATION") in ("-60.0", "-60"), params(os.path.join(DATA, "arcsrc_RTP.tif")))
pd = params(os.path.join(DATA, "arcsrc_DirC.tif"))
check("line-noise settings include the wedge", pd.get("WEDGE") in ("40.0", "40") and pd.get("AZIMUTH") in ("0.0", "0"), pd)
check("no sidecar beside a GeoTIFF", not os.path.exists(os.path.join(DATA, "arcsrc_RTP.tif.sgt.xml")))

# ---- Conv filters
app.v_mean.set(True); app.v_mean_n.set("5")
app.v_ss_std.set(True); app.v_ss_win.set("5")
app._do_apply_conv_stats(src)
pump(0.3)
r = sm.read_sgt_metadata(os.path.join(DATA, "arcsrc_Mn.tif"))
check("mean filter provenance", r is not None and r.findtext("operation") == "Mean filter")
check("windowed statistic provenance",
      params(os.path.join(DATA, "arcsrc_SS_StdDev.tif")).get("STATISTIC") == "std")

# ---- metadata viewer helpers (no windows opened in the test)
app.v_input.set(os.path.join(DATA, "arcsrc_RTP.tif"))
found = app._selected_metadata()
check("selected grid's metadata is found", found is not None and found[0].endswith("arcsrc_RTP.tif"))
txt = sm.format_sgt_metadata_text(found[1])
check("history text lists the step", "Reduction to pole" in txt and "INCLINATION" in txt)
app._save_metadata_xml()
check("Save as XML writes <grid>.sgt.xml", os.path.exists(os.path.join(DATA, "arcsrc_RTP.tif.sgt.xml")))

# ---- build a chain: RTP -> derivative, then replay it on another grid
chain1 = os.path.join(DATA, "arcsrc_RTP.tif")
import sgt_workflow as wf
from calcs.workflow import history_steps
import arcgis_utils as au
op, rec = au.step_record("derivative", DIRECTION="z", POWER=1.0)
chain2 = os.path.join(DATA, "arcsrc_RTP_dz.tif")
shutil.copy(chain1, chain2)
au.write_provenance(chain2, chain1, op, rec)
# (a real derivative result is not needed for the history test)
origin, steps = history_steps(sm.read_sgt_metadata(chain2))
plan = wf.plan_steps(steps)
check("an ArcGIS-made history can be replayed", len(plan) == 2 and all(p["replayable"] for p in plan)
      and plan[0]["operation"] == "Reduction to pole" and plan[1]["operation"] == "Derivative", [p["operation"] for p in plan])

target = make_tif("arctarget.tif", size=80, origin=(620000, 6800000), seed=5)
tinfo = au.raster_header_info(target)
check("target header info", tinfo and tinfo["authid"] == "EPSG:28350" and tinfo["width"] == 80, tinfo)
check("IGRF can be calculated for the target (no pyproj needed)", wf.check_igrf_possible(tinfo) is None)

# the replay window: drive it as the user would, without showing it
app.v_input.set(chain2)
win = G.ReplayWindow(app, "arcsrc_RTP_dz.tif", plan)
win.withdraw()
win.v_target.set(target)
win.v_date.set("2021-03-01")
win._on_click  # exists
# untick Save for the first step, as a user would by clicking
win.save[0] = False
items = [dict(p, save=win.save[i]) for i, p in enumerate(plan) if win.apply[i]]
class T:
    @staticmethod
    def check(percent=None): app._check_cancel()
res = wf.run_replay(target, items, tinfo, when=datetime.date(2021, 3, 1), task=T())
names = [os.path.basename(p) for p in res["outputs"]]
check("replay on the target: the unsaved step is deleted, the last kept",
      names == ["arctarget_RTP_d1z.tif"] and not os.path.exists(os.path.join(DATA, "arctarget_RTP.tif")), names)
fin = res["final"]
steps_t = history_steps(sm.read_sgt_metadata(fin))[1]
check("the replayed grid has the full history", [s["operation"] for s in steps_t] == ["Reduction to pole", "Derivative"],
      [s["operation"] for s in steps_t])
first = [p.text for p in sm.read_sgt_metadata(fin).iter("parameter") if p.get("name") == "INCLINATION"]
check("RTP got a new inclination from IGRF", first and abs(float(first[0]) - (-60)) > 0.5, first)
win.destroy()

# ---- the replay window works end to end through the app (threaded, with Cancel available)
app.v_input.set(chain2)
win2 = G.ReplayWindow(app, "arcsrc_RTP_dz.tif", plan)
win2.withdraw()
win2.v_target.set(target)
win2.v_date.set("2021-03-01")
win2._apply()          # closes the window and starts the job
run_loop(lambda: app._active_token is None and "Applied" in app.v_status.get(), 90)
run_loop(lambda: False, 0.3)
check("replay through the window finished", app._active_token is None and "Applied 2 steps" in app.v_status.get(), app.v_status.get())
check("results are queued for ArcGIS Pro", os.path.exists(app._pending_file) and
      "arctarget_RTP_d1z.tif" in open(app._pending_file).read())

# ---- cancel ends a job at once; the abandoned work saves nothing
state = {}
def slow():
    time.sleep(1.5)            # an uninterruptible calculation
    app._check_cancel()        # the next checkpoint
    state["wrote"] = True
app._run_in_thread(slow)
run_loop(lambda: False, 0.2)
t0 = time.time()
check("Cancel is enabled while a job runs", str(app._cancel_btn["state"]) == "normal")
app._cancel_job()
check("cancelling frees the GUI at once", app._active_token is None and str(app._cancel_btn["state"]) == "disabled"
      and app.v_status.get() == "Cancelled." and time.time() - t0 < 0.5)
run_loop(lambda: False, 2.0)
check("the abandoned job stops at its checkpoint without saving", "wrote" not in state)
check("no error dialog from a cancelled job", "Error" not in app.v_status.get())

# ---- threshold writes provenance too
app.v_thresh_cond.set("above"); app.v_thresh_above.set("50"); app.v_thresh_below.set("-9e9")
out_t = os.path.join(DATA, "arcsrc_thresh.tif")
import importlib.util
check('shapely is not installed in this environment (so the threshold must not need it)', importlib.util.find_spec('shapely') is None)
app._do_threshold(src, out_t)
check("threshold provenance", sm.read_sgt_metadata(out_t) is not None and
      sm.read_sgt_metadata(out_t).findtext("operation") == "Threshold to NaN")

# ---- MRVBF / MRRTF / slope
def make_dem(name, epsg=28350, size=120, cell=25.0):
    path = os.path.join(DATA, name)
    yy, xx = np.mgrid[0:size, 0:size]
    # a flat-floored valley (|x-60| <= 8), walls rising 6 m per 25 m cell, then a plateau
    z = (np.clip((np.abs(xx - 60) - 8) * 6.0, 0.0, 120.0) + 0.02 * yy
         + np.random.default_rng(7).normal(0, 0.05, (size, size))).astype("float32")
    ds = gdal.GetDriverByName("GTiff").Create(path, size, size, 1, gdal.GDT_Float32)
    ds.SetGeoTransform([500000 if epsg == 28350 else 115.0, cell if epsg == 28350 else 0.001, 0,
                        6900000 if epsg == 28350 else -31.0, 0, -(cell if epsg == 28350 else 0.001)])
    srs = osr.SpatialReference(); srs.ImportFromEPSG(epsg); ds.SetProjection(srs.ExportToWkt())
    b = ds.GetRasterBand(1); b.SetNoDataValue(-9999.0); b.WriteArray(z); ds = None
    return path

dem = make_dem("dem.tif")
app.v_input.set(dem)
app.v_mrvbf.set(True)
app._compute_output_path()
prev = app.v_output_preview.get()
check("the Output (auto) line lists the four MRVBF results",
      "_MRVBF" in prev and "_MRRTF" in prev and "_slope" in prev, prev)
app._do_apply_conv_stats(dem)
pump(0.3)
names = {n: os.path.join(DATA, n) for n in ("dem_MRVBF.tif", "dem_MRRTF.tif", "dem_slope.tif", "dem_MRVBF_RGB.tif")}
check("MRVBF, MRRTF, slope and the RGB composite are written", all(os.path.exists(v) for v in names.values()),
      [n for n, v in names.items() if not os.path.exists(v)])
ds = gdal.Open(names["dem_MRVBF_RGB.tif"])
check("the composite has three bands", ds is not None and ds.RasterCount == 3)
ds = None
def band(path, i=1):
    d_ = gdal.Open(path); a = d_.GetRasterBand(i).ReadAsArray().astype(float)
    nd = d_.GetRasterBand(i).GetNoDataValue(); d_ = None
    if nd is not None: a[a == nd] = np.nan
    return a
v = band(names["dem_MRVBF.tif"]); sl = band(names["dem_slope.tif"])
rf = band(names["dem_MRRTF.tif"])
floor_v, plateau_v = np.nanmean(v[:, 54:66]), np.nanmean(v[:, 100:118])
check("MRVBF is higher on the valley floor than on the plateau",
      np.nanmax(v) > 1.0 and floor_v > plateau_v, (np.nanmax(v), floor_v, plateau_v))
check("MRRTF is higher on the plateau than on the valley floor",
      np.nanmean(rf[:, 100:118]) > np.nanmean(rf[:, 54:66]), (np.nanmean(rf[:, 100:118]), np.nanmean(rf[:, 54:66])))
check("slope is in degrees: about 13.5 on the walls (6 m per 25 m cell), near 0 on the floor",
      12.5 < np.nanmedian(sl[:, 70:80]) < 14.5 and np.nanmedian(sl[:, 54:66]) < 2.0,
      (np.nanmedian(sl[:, 70:80]), np.nanmedian(sl[:, 54:66])))
r = sm.read_sgt_metadata(names["dem_MRVBF.tif"])
check("MRVBF provenance records its settings",
      r is not None and r.findtext("operation") == "Multiresolution valley bottom flatness"
      and {p.get("name") for p in r.iter("parameter")} >= {"t_slope", "max_res"})
geo = make_dem("dem_geo.tif", epsg=4326)
check("a geographic DEM is detected (MRVBF is refused for it)", au.raster_header_info(geo)["geographic"])


def make_points(name, with_sentinels=False):
    """Points on a plane z = 0.02*(x-500000) + 0.01*(y-6890000); the value is stored as text
    (as the SGTool importer writes it) in 'val' and as a number in 'num'."""
    from osgeo import ogr
    path = os.path.join(DATA, name)
    drv = ogr.GetDriverByName("ESRI Shapefile")
    ds = drv.CreateDataSource(path)
    srs = osr.SpatialReference(); srs.ImportFromEPSG(28350)
    layer = ds.CreateLayer("pts", srs, ogr.wkbPoint)
    layer.CreateField(ogr.FieldDefn("val", ogr.OFTString))
    layer.CreateField(ogr.FieldDefn("num", ogr.OFTReal))
    rng = np.random.default_rng(11)
    n = 400
    xs = 500000 + rng.uniform(0, 2000, n); ys = 6890000 + rng.uniform(0, 2000, n)
    zs = 0.02 * (xs - 500000) + 0.01 * (ys - 6890000)
    if with_sentinels:
        zs[:8] = -99999.0
    for x, y, z in zip(xs, ys, zs):
        f = ogr.Feature(layer.GetLayerDefn()); g = ogr.Geometry(ogr.wkbPoint); g.AddPoint(float(x), float(y))
        f.SetGeometry(g); f.SetField("val", str(z)); f.SetField("num", float(z)); layer.CreateFeature(f)
    ds = None
    return path


def plane_error(tif):
    d_ = gdal.Open(tif); a = d_.GetRasterBand(1).ReadAsArray().astype(float); nd = d_.GetRasterBand(1).GetNoDataValue()
    gt = d_.GetGeoTransform(); d_ = None
    if nd is not None: a[a == nd] = np.nan
    rows, cols = a.shape
    xs = gt[0] + (np.arange(cols) + 0.5) * gt[1]; ys = gt[3] + (np.arange(rows) + 0.5) * gt[5]
    X, Y = np.meshgrid(xs, ys)
    truth = 0.02 * (X - 500000) + 0.01 * (Y - 6890000)
    inner = (X > 500300) & (X < 501700) & (Y > 6890300) & (Y < 6891700)
    return float(np.nanmax(np.abs((a - truth)[inner]))), a, gt

# ---- multilevel B-spline gridding
pts = make_points("bs_pts.shp")
app._set_bspline_points(pts)
check("point fields are listed after choosing the file", "val" in app._bs_field_combo["values"] and "num" in app._bs_field_combo["values"],
      app._bs_field_combo["values"])
app._do_bspline(pts, "val", 50.0, 1e-4, 11, None)    # text field, as the importer writes it
pump(0.3)
out_bs = os.path.join(DATA, "bs_pts_val_bspline.tif")
check("gridding writes <points>_<field>_bspline.tif (QGIS naming)", os.path.exists(out_bs))
err, arr_bs, gt_bs = plane_error(out_bs)
check("the gridded surface follows the plane away from the edges", err < 1.0, err)
check("cells are the requested size and the CRS is kept",
      abs(gt_bs[1] - 50.0) < 1e-9 and "EPSG:28350" == au.raster_header_info(out_bs)["authid"], gt_bs)
rb = sm.read_sgt_metadata(out_bs)
pbs = {p_.get("name"): p_.text for p_ in rb.iter("parameter")} if rb is not None else {}
check("provenance records the gridding settings",
      rb is not None and rb.findtext("operation") == "Multilevel B-spline gridding"
      and pbs.get("data_field") == "val" and pbs.get("points_used") == "400" and pbs.get("max_levels") == "11", pbs)
pts2 = make_points("bs_sent.shp", with_sentinels=True)
app._do_bspline(pts2, "num", 50.0, 1e-4, 11, -90000.0)
err2, _a, _g = plane_error(os.path.join(DATA, "bs_sent_num_bspline.tif"))
check("points below the ignore value are dropped before gridding", err2 < 1.0, err2)
pbs2 = {p_.get("name"): p_.text for p_ in sm.read_sgt_metadata(os.path.join(DATA, "bs_sent_num_bspline.tif")).iter("parameter")}
check("the ignore value and the points used are recorded",
      pbs2.get("points_used") == "392" and pbs2.get("ignore_values_below") in ("-90000.0", "-90000"), pbs2)
refused = False
try:
    au.grid_points_to_file([0, 1], [0, 1], [1, 2], os.path.join(DATA, "x.tif"), None, "z", 1.0)
except ValueError:
    refused = True
check("fewer than 3 points is refused with a clear error", refused)
ran_cancelled = {}
def slow_cb(fraction=None):
    raise G.OperationCancelled()
try:
    au.grid_points_to_file(*read_points_ogr_(pts, "num")[:3], os.path.join(DATA, "c.tif"), pts, "num", 50.0, callback=slow_cb)
except G.OperationCancelled:
    ran_cancelled["yes"] = True
check("a cancel raised from the gridder's progress callback stops it", ran_cancelled.get("yes") and not os.path.exists(os.path.join(DATA, "c.tif")))

# ---- external preview
import math
from sgt_preview_window import block_mean, stretch_to_grey
for v in (app.v_rtp, app.v_deriv, app.v_mean, app.v_dirclean, app.v_regrem):
    v.set(False)
pv = make_tif("prevsrc.tif", size=100, cell=100.0, seed=3)
app.v_input.set(pv)
app.nb.select(0); app.update()

# the calculation: whole grid <= 600 cells is the full-resolution result Apply would write
app.v_deriv.set(True); app.v_deriv_dir.set("z"); app.v_deriv_pow.set("1")
res = app._compute_preview(0, pv, "whole", (0, 0, 100, 100))
app._do_apply_fft(pv)                      # what Apply writes
pump(0.3)
written = band(os.path.join(DATA, "prevsrc_d1z.tif"))
check("whole-grid preview equals the filter Apply writes", np.allclose(res["filtered"], written, equal_nan=True, atol=1e-6),
      float(np.nanmax(np.abs(res["filtered"] - written))))
check("the preview reports the operation and the cell count", "Derivative" in res["text"] and "100 x 100" in res["text"], res["text"])
app.v_deriv.set(False)

# view mode: only what is in view is calculated, from the full-resolution cells
app.v_deriv.set(True)
res = app._compute_preview(0, pv, "view", (20.0, 30.0, 70.0, 80.0))
from calcs.GeophysicalProcessor import GeophysicalProcessor
arr_pv = au.raster_to_numpy(pv)[0]
crop = arr_pv[30:80, 20:70]
buf_ = min(int(app.v_buf.get() or 10), 50, 50)   # the GUI's FFT buffer, limited by the 50-cell crop
expect = GeophysicalProcessor(100.0, 100.0, buf_).compute_derivative(crop, direction="z", order=1.0, buffer_size=buf_)
check("view preview covers the view box", res["box"] == (20, 30, 70, 80) and res["filtered"].shape == (50, 50), (res["box"], res["filtered"].shape))
check("view preview is the filter of that crop at full detail", np.allclose(res["filtered"], expect, equal_nan=True, atol=1e-6))
check("the original crop is returned for the toggle", np.allclose(res["original"], crop, equal_nan=True))
app.v_deriv.set(False)

# big grids are block-averaged for the preview
big_pv = make_tif("prevbig.tif", size=700, cell=50.0, seed=4)
app.v_mean.set(True); app.v_mean_n.set("3")
res = app._compute_preview(1, big_pv, "whole", (0, 0, 700, 700))
check("a 700-cell grid is reduced 2 x for the whole-grid preview",
      res["filtered"].shape == (350, 350) and res["box"] == (0, 0, 700, 700) and "block-averaged 2" in res["text"],
      (res["filtered"].shape, res["text"]))
check("block_mean averages blocks and ignores no data",
      np.allclose(block_mean(np.array([[1.0, 3.0], [np.nan, 5.0]]), 2), [[3.0]]))
app.v_mean.set(False)

# sun shading is stretched min to max, the rest by standard deviation
app.nb.select(1); app.update()
app.v_sun.set(True)
res = app._compute_preview(1, pv, "whole", (0, 0, 100, 100))
check("sun shading is not stretched by standard deviation", res["std_clip"] is False)
app.v_sun.set(False)
rgb = stretch_to_grey(np.array([[0.0, 1.0], [np.nan, 0.5]]), std_clip=False)
check("the grey stretch maps min to 0, max to 255 and no data to light grey",
      rgb[0, 0, 0] == 0 and rgb[0, 1, 0] == 255 and rgb[1, 0, 0] == 217, rgb[:, :, 0].tolist())

# the window, the controls and following between tabs (needs the real main loop)
app.nb.select(0); app.update()
app.v_input.set(pv)
app.v_preview.set(True)
run_loop(lambda: app._preview_win is not None, 5)
win = app._preview_win
check("ticking Preview opens the preview window", win is not None and win.winfo_exists())
check("the Apply buttons read Apply Preview", all(b.cget("text") == "Apply Preview" for b in app._apply_btns.values()),
      [b.cget("text") for b in app._apply_btns.values()])
app.v_deriv.set(True)
app.v_rtp.set(True)
check("ticking a second filter unticks the first", app.v_rtp.get() and not app.v_deriv.get())
app.v_rtp.set(False)
app.v_deriv.set(True)
run_loop(lambda: win._result is not None, 40)
check("the window shows the filtered grid", win._result is not None and win._result["box"] == (0, 0, 100, 100) and "Derivative" in win.v_status.get(),
      win.v_status.get())
check("the window knows the grid size", (win.rows, win.cols) == (100, 100))
win.zoom(2.0, 100, 100)
check("zooming in narrows the view", win.view_box()[2] - win.view_box()[0] < 100, win.view_box())
view_before = win.view_box()
run_loop(lambda: win._result is not None and win._result["box"] != (0, 0, 100, 100), 40)
check("after zooming the filter is recalculated for the view",
      win._result["box"][2] - win._result["box"][0] < 100, win._result["box"])
win.zoom(0.1, 10, 10)
check("zooming out cannot go past the whole grid by much", win.cpp <= max(100 / win._canvas_size()[0], 100 / win._canvas_size()[1]) * 1.5 + 1e-9)
win.reset_view()
check("double-click / reset shows the whole grid", win.view_box() == (0.0, 0.0, 100.0, 100.0), win.view_box())

# following between tabs
app.nb.select(1); app.update()
check("changing to Conv + Stats unticks the FFT filter and moves the preview there",
      not app.v_deriv.get() and app.v_preview.get() and app._preview_tab == 1)
app.nb.select(2); app.update()
check("tabs without a preview just wait", app.v_preview.get() and app._preview_tab is None and app._preview_win is not None)
app.nb.select(0); app.update()
check("back on an FFT tab the preview is there again", app._preview_tab == 0 and app.v_preview.get())

# Apply Preview: saves the previewed filter at full resolution and ends the preview
app.v_thg.set(True)
app._apply_fft()
check("Apply Preview closes the preview but keeps the filter ticked for the run",
      not app.v_preview.get() and app._preview_win is None)
run_loop(lambda: app._active_token is None and os.path.exists(os.path.join(DATA, "prevsrc_THG.tif")), 60)
run_loop(lambda: False, 0.3)
check("Apply Preview wrote the full-resolution result", os.path.exists(os.path.join(DATA, "prevsrc_THG.tif")))
check("the filter is unticked afterwards and the buttons read Apply again",
      not app.v_thg.get() and app._apply_btns[0].cget("text") == "Apply FFT Processing", app._apply_btns[0].cget("text"))

# closing the window ends the preview and unticks the filters
app.v_preview.set(True); run_loop(lambda: app._preview_win is not None, 5)
app.v_agc.set(True)
app._preview_closed()
check("closing the preview window ends the preview and unticks the filters",
      not app.v_preview.get() and app._preview_win is None and not app.v_agc.get())

# ---- the preview starts where the ArcGIS Pro map is looking
import json, tempfile
# a private temp folder, so a real map-view file from a Pro session is never touched
_real_tmp = tempfile.tempdir
tempfile.tempdir = DATA
view_file = os.path.join(DATA, "sgtool_view.json")
if os.path.exists(view_file):
    os.remove(view_file)

# prevsrc.tif: 100 x 100 cells of 100 m from (500000, 6900000); cells (20,30)-(70,80)
view_same = {"extent": [502000.0, 6892000.0, 507000.0, 6897000.0], "epsg": 28350, "wkt": ""}
box = au.map_view_to_cells(view_same, pv)
check("a map view in the grid's own CRS becomes a box of cells", box is not None and
      all(abs(a - b) < 1e-6 for a, b in zip(box, (20, 30, 70, 80))), box)
srs_g = osr.SpatialReference(); srs_g.ImportFromEPSG(28350)
srs_ll = osr.SpatialReference(); srs_ll.ImportFromEPSG(4326)
for r_ in (srs_g, srs_ll):
    r_.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
to_ll = osr.CoordinateTransformation(srs_g, srs_ll)
ll_corners = [to_ll.TransformPoint(x, y)[:2] for x, y in ((502000, 6892000), (507000, 6897000))]
view_ll = {"extent": [ll_corners[0][0], ll_corners[0][1], ll_corners[1][0], ll_corners[1][1]], "epsg": 4326, "wkt": ""}
box_ll = au.map_view_to_cells(view_ll, pv)
check("a map in another CRS (lat/long) is transformed to the grid's cells",
      box_ll is not None and all(abs(a - b) < 1.5 for a, b in zip(box_ll, (20, 30, 70, 80))), box_ll)
srs_wkt = srs_g.ExportToWkt()
box_wkt = au.map_view_to_cells({"extent": [502000.0, 6892000.0, 507000.0, 6897000.0], "epsg": 0, "wkt": srs_wkt}, pv)
check("a map CRS given as text works too", box_wkt is not None and abs(box_wkt[0] - 20) < 1e-6, box_wkt)
partly = au.map_view_to_cells({"extent": [499000.0, 6890000.0, 503000.0, 6901000.0], "epsg": 28350, "wkt": ""}, pv)
check("a map view half off the grid is clipped to the grid", partly is not None and partly[0] == 0.0 and partly[1] == 0.0, partly)
check("a map view that misses the grid gives nothing",
      au.map_view_to_cells({"extent": [600000.0, 6000000.0, 601000.0, 6001000.0], "epsg": 28350, "wkt": ""}, pv) is None)

# the preview window opens on that view
with open(view_file, "w") as f:
    json.dump(view_same, f)
app.v_input.set(pv)
app.nb.select(0); app.update()
app.v_deriv.set(True)
app.v_preview.set(True)
run_loop(lambda: app._preview_win is not None and app._preview_win._result is not None, 60)
win = app._preview_win
vx0, vy0, vx1, vy1 = win._visible()
check("the preview opens centred on the map view",
      abs((vx0 + vx1) / 2 - 45) < 1.0 and abs((vy0 + vy1) / 2 - 55) < 1.0, ((vx0 + vx1) / 2, (vy0 + vy1) / 2))
check("and the whole map extent is visible", vx0 <= 20 + 1e-6 and vx1 >= 70 - 1e-6 and vy0 <= 30 + 1e-6 and vy1 >= 80 - 1e-6,
      (vx0, vy0, vx1, vy1))
check("the filter was calculated for that part of the grid, not the whole grid",
      win._result["box"][2] - win._result["box"][0] < 100 and win._result["box"][0] <= 20 + 1.5, win._result["box"])
win.reset_view()
check("double-click shows the whole grid again", np.allclose(win.view_box(), (0.0, 0.0, 100.0, 100.0), atol=0.01), win.view_box())
app._preview_match_map()
vx0, vy0, vx1, vy1 = win._visible()
check("the Match map view button goes back to the map's extent",
      abs((vx0 + vx1) / 2 - 45) < 1.0 and abs((vy0 + vy1) / 2 - 55) < 1.0, ((vx0 + vx1) / 2, (vy0 + vy1) / 2))
os.remove(view_file)
win.reset_view()
app._preview_match_map()
check("without a map view the button says why and leaves the view alone",
      "No ArcGIS Pro map view" in win.v_status.get() and np.allclose(win.view_box(), (0.0, 0.0, 100.0, 100.0), atol=0.01), win.v_status.get())
app.v_preview.set(False)
tempfile.tempdir = _real_tmp

# ---- Threshold to NaN preview, and the shape of the preview
for v in (app.v_deriv, app.v_rtp, app.v_mean, app.v_thresh_en):
    v.set(False)
app.nb.select(3); app.update()
app.v_input.set(pv)
app.v_thresh_cond.set("above"); app.v_thresh_above.set("40"); app.v_thresh_below.set("-9e9")
app.v_thresh_en.set(True)
res = app._compute_preview(3, pv, "whole", (0, 0, 100, 100))
check("the threshold preview marks the values above the threshold as no data",
      "Threshold to NaN" in res["text"] and np.array_equal(np.isnan(res["filtered"]), res["original"] > 40)
      and np.isnan(res["filtered"]).any(), res["text"])
out_thr = os.path.join(DATA, "prevsrc_Clean.tif")
if os.path.exists(out_thr):
    os.remove(out_thr)
app._do_threshold(pv, out_thr)
pump(0.3)
check("the threshold preview matches what Apply writes",
      np.array_equal(np.isnan(res["filtered"]), np.isnan(band(out_thr))))
app.v_thresh_en.set(False)

# following onto the Utils tab, and Apply Preview there
app.nb.select(0); app.update()
app.v_preview.set(True)
run_loop(lambda: app._preview_win is not None, 5)
app.nb.select(3); app.update()
check("Preview follows onto the Utils tab", app.v_preview.get() and app._preview_tab == 3
      and app._apply_btns[3].cget("text") == "Apply Preview", app._apply_btns[3].cget("text"))
app.v_thresh_en.set(True)
win = app._preview_win
run_loop(lambda: win._result is not None and "Threshold" in win.v_status.get(), 40)
check("the window shows the threshold", win._result is not None and "Threshold" in win.v_status.get(), win.v_status.get())
if os.path.exists(out_thr):
    os.remove(out_thr)
app._apply_threshold()
run_loop(lambda: app._active_token is None and os.path.exists(out_thr), 60)
run_loop(lambda: not app.v_thresh_en.get(), 5)
check("Apply Preview on the Utils tab writes the result, unticks the threshold and ends the preview",
      os.path.exists(out_thr) and not app.v_thresh_en.get() and not app.v_preview.get() and app._preview_win is None)
check("the button reads Apply Threshold again", app._apply_btns[3].cget("text") == "Apply Threshold \u2192 NaN", app._apply_btns[3].cget("text"))

# the preview window takes the shape of the ArcGIS Pro map canvas
_real_tmp = tempfile.tempdir
tempfile.tempdir = DATA
view_file = os.path.join(DATA, "sgtool_view.json")
wide = {"extent": [502000.0, 6893000.0, 507000.0, 6896000.0], "epsg": 28350, "wkt": ""}   # 5000 x 3000 m
with open(view_file, "w") as f:
    json.dump(wide, f)
app.nb.select(0); app.update()
app.v_input.set(pv)
app.v_deriv.set(True)
app.v_preview.set(True)
run_loop(lambda: app._preview_win is not None and app._preview_win._result is not None, 60)
win = app._preview_win
run_loop(lambda: False, 0.6)      # let the window manager resize the window
cw_, ch_ = win.canvas.winfo_width(), win.canvas.winfo_height()
check("the preview canvas has the map canvas's shape (5:3)", abs(cw_ / ch_ - 5 / 3) < 0.06, (cw_, ch_, cw_ / ch_))
vx0, vy0, vx1, vy1 = win._visible()
check("and shows exactly the map's extent (cells 20-70 by 40-70)",
      abs(vx0 - 20) < 1.5 and abs(vx1 - 70) < 1.5 and abs(vy0 - 40) < 1.5 and abs(vy1 - 70) < 1.5, (vx0, vy0, vx1, vy1))
tall = {"extent": [503000.0, 6892000.0, 505000.0, 6896000.0], "epsg": 28350, "wkt": ""}  # 2000 x 4000 m, 1:2
with open(view_file, "w") as f:
    json.dump(tall, f)
app._preview_match_map()
run_loop(lambda: False, 0.6)      # let the window manager resize the window
cw_, ch_ = win.canvas.winfo_width(), win.canvas.winfo_height()
vx0, vy0, vx1, vy1 = win._visible()
check("the tall map is also fitted exactly (cells 30-50 by 40-80)", abs(vx0 - 30) < 2.0 and abs(vx1 - 50) < 2.0 and abs(vy0 - 40) < 2.0 and abs(vy1 - 80) < 2.0, (vx0, vy0, vx1, vy1))
check("Match map view reshapes the window for a tall map (1:2)", abs(cw_ / ch_ - 0.5) < 0.04, (cw_, ch_, cw_ / ch_))
os.remove(view_file)
app.v_preview.set(False)

# cells that are not square are drawn true
ns = os.path.join(DATA, "prevnonsq.tif")
ds = gdal.GetDriverByName("GTiff").Create(ns, 80, 40, 1, gdal.GDT_Float32)
ds.SetGeoTransform([500000, 100, 0, 6900000, 0, -50])
srs_ = osr.SpatialReference(); srs_.ImportFromEPSG(28350); ds.SetProjection(srs_.ExportToWkt())
b_ = ds.GetRasterBand(1); b_.SetNoDataValue(-9999.0); b_.WriteArray(np.random.default_rng(2).normal(0, 1, (40, 80)).astype("float32")); ds = None
app.v_input.set(ns)
app.v_deriv.set(True)
app.v_preview.set(True)
run_loop(lambda: app._preview_win is not None and app._preview_win._result is not None, 60)
win = app._preview_win
check("the window knows the cells are twice as wide as high", abs(win.ratio - 0.5) < 1e-9, win.ratio)
cw_, ch_ = win._canvas_size()
shown_w = 80 / win.cpp; shown_h = 40 * win.ratio / win.cpp
check("the whole grid is drawn in its true proportions (8000 m by 2000 m = 4:1)", abs(shown_w / shown_h - 4.0) < 0.01, shown_w / shown_h)
check("fits the canvas", shown_w <= cw_ + 1 and shown_h <= ch_ + 1)
app.v_preview.set(False)
tempfile.tempdir = _real_tmp

passed = sum(1 for _, ok in results if ok)
print(f"\n{passed}/{len(results)} ARC checks passed", flush=True)
if any(not ok for _, ok in results):
    print("FAILED:", [n for n, ok in results if not ok])
os._exit(0 if all(ok for _, ok in results) else 1)
