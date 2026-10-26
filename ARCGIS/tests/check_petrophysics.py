"""
Checks of the Petrophys tab of the ArcGIS Pro GUI (SGTool_GUI.py) and the petrophysics
toolbox tools, with the real tkinter main loop.

Run it with ArcGIS Pro's own Python:

    "C:\\Program Files\\ArcGIS\\Pro\\bin\\Python\\envs\\arcgispro-py3\\python.exe" ARCGIS\\tests\\check_petrophysics.py

It builds synthetic grids, polygons and sample points in a temporary folder, drives the
tab, prints PASS / FAIL per check and a total, and exits with status 1 if anything failed.
Plot windows open briefly.
"""
import os, sys, tempfile, time, traceback
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = tempfile.mkdtemp(prefix="sgtool_arc_petro_")
sys.path.insert(0, os.path.join(REPO, "ARCGIS"))
sys.path.insert(0, REPO)

from osgeo import gdal, ogr, osr

results = []


def check(name, cond, detail=""):
    results.append((name, bool(cond)))
    print(("PASS  " if cond else "FAIL  ") + name + (f"   [{detail}]" if detail and not cond else ""), flush=True)


def make_tif(name, arr, origin=(500000, 6900000), cell=100.0, epsg=28350):
    path = os.path.join(DATA, name)
    rows, cols = arr.shape
    ds = gdal.GetDriverByName("GTiff").Create(path, cols, rows, 1, gdal.GDT_Float32)
    ds.SetGeoTransform([origin[0], cell, 0, origin[1], 0, -cell])
    srs = osr.SpatialReference(); srs.ImportFromEPSG(epsg)
    ds.SetProjection(srs.ExportToWkt())
    b = ds.GetRasterBand(1); b.SetNoDataValue(-9999.0); b.WriteArray(arr.astype("float32"))
    ds = None
    return path


def make_vector(name, geom_type, epsg, fields, features):
    """features: (geometry WKT, {field: value}) in the CRS of epsg."""
    path = os.path.join(DATA, name)
    srs = osr.SpatialReference(); srs.ImportFromEPSG(epsg)
    srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    ds = ogr.GetDriverByName("ESRI Shapefile").CreateDataSource(path)
    layer = ds.CreateLayer("data", srs, geom_type)
    for fname, ftype in fields:
        layer.CreateField(ogr.FieldDefn(fname, ftype))
    for wkt, values in features:
        feat = ogr.Feature(layer.GetLayerDefn())
        for k, v in values.items():
            if v is not None:
                feat.SetField(k, v)
        feat.SetGeometry(ogr.CreateGeometryFromWkt(wkt))
        layer.CreateFeature(feat)
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


# ---- synthetic data: two units, left dense and magnetic, right light and weakly magnetic
rng = np.random.default_rng(3)
yy, xx = np.mgrid[0:120, 0:120]
mag = np.where(xx < 60, 300.0, -100.0) + rng.normal(0, 20, xx.shape)      # nT, RTP anomaly
grav = np.where(xx < 60, -200.0, -500.0) + rng.normal(0, 5, xx.shape)     # um/s2, Bouguer
mag_tif = make_tif("mag_rtp.tif", mag)
grav_tif = make_tif("grav.tif", grav[::2, ::2], cell=200.0)                # coarser grid
x0, y0, x1, y1 = 500000, 6900000 - 12000, 512000, 6900000
# polygons in geographic coordinates (another CRS), two units
to_ll = osr.CoordinateTransformation(
    *[(lambda s: (s.ImportFromEPSG(e), s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER), s)[2])
      (osr.SpatialReference()) for e in (28350, 4326)])
def ll_box(a, b):
    pts = [to_ll.TransformPoint(x, y)[:2] for x, y in ((a, y0), (b, y0), (b, y1), (a, y1), (a, y0))]
    return "POLYGON((" + ",".join(f"{px} {py}" for px, py in pts) + "))"
poly_shp = make_vector("geology.shp", ogr.wkbPolygon, 4326, [("LITH", ogr.OFTString)],
                       [(ll_box(x0, 506000), {"LITH": "Basalt"}),
                        (ll_box(506000, x1), {"LITH": "Granite"})])
samples = []
for i in range(80):
    unit = ["BIF", "Basalt", "Granite"][i % 3]
    rho = 2.6 + 0.2 * (i % 3) + rng.normal(0, 0.03)
    samples.append((f"POINT({500100 + 100 * i} {6895000})",
                    {"RHO": f"{rho:.3f}".replace(".", ","),          # text, decimal comma
                     "CHI": float(10 ** rng.normal(-3 + i % 3, 0.3)),
                     "UNIT": unit if i % 20 else None}))
points_shp = make_vector("samples.shp", ogr.wkbPoint, 28350,
                         [("RHO", ogr.OFTString), ("CHI", ogr.OFTReal), ("UNIT", ogr.OFTString)],
                         samples)

# ---- the GUI, with message boxes captured instead of shown
import SGTool_GUI as G
boxes = []
for name in ("showwarning", "showerror", "showinfo"):
    setattr(G.messagebox, name, lambda *a, n=name, **k: boxes.append((n,) + a))

app = G.SGToolApp()
app.v_input.set(mag_tif)
app.v_preview.set(True)          # Preview on, then arrive on the Petrophys tab
app.nb.select(3); app.update()
check("the Petrophys tab is the fourth tab and is built", app.nb.tab(3, "text") == "Petrophys"
      and hasattr(app, "v_app_sus") and app.nb.tab(4, "text") == "Utils")
check("both calculations start unticked", not app.v_app_sus.get() and not app.v_app_dens.get())
run_loop(lambda: False, 2)  # long enough for a preview calculation to show
win = app._preview_win
check("arriving on the tab with Preview on calculates nothing",
      win is not None and win._result is None, None if win is None else win.v_status.get())
app.v_preview.set(False); app.update()
app.v_app_sus.set(True); app.v_app_dens.set(True)
check("ticking one unticks the other", app.v_app_dens.get() and not app.v_app_sus.get())
app.v_app_dens.set(False)
check("both can be unticked", not app.v_app_sus.get() and not app.v_app_dens.get())

# a Petrophys tick does not survive Preview choosing a filter on another tab
app.v_app_sus.set(True)
app.nb.select(0); app.update()
app.v_preview.set(True); app.update()
check("turning Preview on elsewhere unticks the Petrophys calculation", not app.v_app_sus.get())
app.nb.select(3); app.update()
app.v_app_dens.set(True); app.update()
app.nb.select(0); app.update()
app.v_thg.set(True); app.update()
check("choosing a filter on another tab while previewing unticks it",
      not app.v_app_dens.get() and app.v_thg.get())
app.v_preview.set(False); app.v_thg.set(False); app.update()
app.nb.select(3); app.update()

app.v_input.set(mag_tif)
app._petro_igrf()
field = float(app.v_app_field.get())
check("IGRF fills in the field strength for the grid's location", 50000 < field < 65000, field)

# ---- apparent susceptibility
app.v_app_sus.set(True)
out_sus = os.path.join(DATA, "mag_rtp_AppSus.tif")
app._apply_petro()
run_loop(lambda: app._active_token is None and os.path.exists(out_sus), 60)
check("Apply writes <grid>_AppSus.tif", os.path.exists(out_sus))
from calcs import sgt_metadata as sm
from calcs.workflow import history_steps
step = history_steps(sm.read_sgt_metadata(out_sus))[1][-1] if os.path.exists(out_sus) else {}
check("its provenance is the shared filter's operation and settings",
      step.get("operation") == "Apparent susceptibility"
      and abs(float(step["parameters"].get("FIELD", 0)) - field) < 1, step)
sus_arr = gdal.Open(out_sus).ReadAsArray() if os.path.exists(out_sus) else np.zeros(1)
check("the result's mean is the background susceptibility", abs(float(np.mean(sus_arr)) - 0.001) < 0.002,
      float(np.mean(sus_arr)))

# ---- apparent density, grid in um/s2
app.v_input.set(grav_tif)
app.v_app_dens.set(True)
app.v_app_units.set("um/s2 (gu)")
out_dens = os.path.join(DATA, "grav_AppDens.tif")
app._apply_petro()
run_loop(lambda: app._active_token is None and os.path.exists(out_dens), 60)
dens_arr = gdal.Open(out_dens).ReadAsArray() if os.path.exists(out_dens) else np.zeros(1)
check("Apply writes <grid>_AppDens.tif with a mean of the background density",
      os.path.exists(out_dens) and abs(float(np.mean(dens_arr)) - 2.67) < 0.01, float(np.mean(dens_arr)))
spread = float(np.percentile(dens_arr, 90) - np.percentile(dens_arr, 10))
check("the two units differ by a plausible density contrast (gu converted to mGal)",
      0.02 < spread < 1.0, spread)

check("Apply unticks the calculation", not app.v_app_sus.get() and not app.v_app_dens.get())

# ---- preview of the ticked calculation
app.v_app_dens.set(True)
res = app._compute_preview(3, grav_tif, "whole", None)
check("the preview runs the chosen calculation", "Apparent density" in res["text"], res["text"])

# ---- plots, saved as PNG beside the inputs
app.v_plot_png.set(True)
pngs = lambda: sorted(f for f in os.listdir(DATA) if f.endswith(".png"))
def plot(fn, what, expect):
    before = set(pngs())
    boxes.clear()
    fn()
    run_loop(lambda: app._active_token is None and set(pngs()) - before, 60)
    new = sorted(set(pngs()) - before)
    check(what, new == [expect] and not boxes, (new, boxes))

app.v_hk_sus.set(out_sus); app.v_hk_dens.set(out_dens)
plot(app._plot_apparent, "Plot Apparent: Henkel heat map of the two grids",
     "grav_AppDens_mag_rtp_AppSus_henkel.png")
app.v_hk_poly.set(poly_shp); app.v_hk_poly_field.set("LITH"); app.v_hk_position.set(True)
plot(app._plot_categories, "Plot per Category: polygons in another CRS, coloured by position",
     "grav_AppDens_mag_rtp_AppSus_henkel_by_LITH.png")
app.v_hk_sus.set("")
plot(app._plot_apparent, "one grid only: histogram", "grav_AppDens_hist.png")
plot(app._plot_categories, "one grid only, per category: histograms", "grav_AppDens_hist_by_LITH.png")

app.v_ms_points.set(points_shp)
app.v_ms_x.set("RHO"); app.v_ms_y.set("CHI"); app.v_ms_cat.set("UNIT")
plot(app._plot_samples, "Plot Samples: Y against X with categories (text numbers read)",
     "samples_CHI_vs_RHO.png")
app.v_ms_x.set("")
plot(app._plot_samples, "one field only, with a category: histograms", "samples_CHI_hist_by_UNIT.png")
windows = [w for w in app.winfo_children() if isinstance(w, G.tk.Toplevel)]
check("each plot opened in its own window", len(windows) >= 6, len(windows))

# ---- grouping similar categories
app.v_hk_sus.set(out_sus); app.v_hk_position.set(False)
app.v_hk_group.set(True); app.v_hk_groups.set("Auto")
os.remove(os.path.join(DATA, "grav_AppDens_mag_rtp_AppSus_henkel_by_LITH.png"))
plot(app._plot_categories, "Plot per Category with Group similar (auto)",
     "grav_AppDens_mag_rtp_AppSus_henkel_by_LITH.png")
app.v_ms_x.set("RHO"); app.v_ms_group.set(True); app.v_ms_groups.set("2")
os.remove(os.path.join(DATA, "samples_CHI_vs_RHO.png"))
plot(app._plot_samples, "Plot Samples with Group similar (2 groups)", "samples_CHI_vs_RHO.png")

# ---- replay planning of the new steps
from sgt_workflow import plan_steps
items = plan_steps(history_steps(sm.read_sgt_metadata(out_dens))[1])
check("an apparent density step can be replayed", items and items[-1]["replayable"]
      and items[-1]["values"]["UNITS"].startswith("um"), items[-1] if items else None)

# ---- toolbox: the new tools load
try:
    import importlib.machinery, importlib.util
    loader = importlib.machinery.SourceFileLoader("sgtool_pyt", os.path.join(REPO, "ARCGIS", "SGTool.pyt"))
    spec = importlib.util.spec_from_loader("sgtool_pyt", loader)
    pyt = importlib.util.module_from_spec(spec)
    loader.exec_module(pyt)
    tools = pyt.Toolbox().tools
    ok = pyt.ApparentSusceptibility in tools and pyt.ApparentDensity in tools
    params = [p.name for p in pyt.ApparentDensity().getParameterInfo()]
    check("the toolbox lists Apparent Susceptibility and Apparent Density", ok and "units" in params, params)
except Exception:
    check("the toolbox loads", False, traceback.format_exc().splitlines()[-1])

app.destroy()
failed = [n for n, ok in results if not ok]
print(f"\n{len(results) - len(failed)} / {len(results)} passed" + (f"; failed: {failed}" if failed else ""))
print("data in", DATA)
sys.exit(1 if failed else 0)
