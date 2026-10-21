"""
Checks of the ArcGIS Pro Clip All Layers tool (ARCGIS/arcgis_clip.py and the toolbox tool).

Run it with ArcGIS Pro's own Python:

    "C:\\Program Files\\ArcGIS\\Pro\\bin\\Python\\envs\\arcgispro-py3\\python.exe" ARCGIS\\tests\\check_clip.py

It builds a small project in a temporary folder (from Esri's own empty project template), clips
it, reopens the saved copy, and prints PASS / FAIL per check. It needs arcpy and an ArcGIS Pro
licence (Basic is enough), but not a running copy of ArcGIS Pro.
"""
import os, sys, shutil, tempfile, importlib.machinery, importlib.util
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, "ARCGIS"))
import arcpy
from osgeo import gdal, osr

DATA = tempfile.mkdtemp(prefix="sgtool_arc_clip_")
results = []


def check(n, c, d=""):
    results.append((n, bool(c)))
    print(("PASS  " if c else "FAIL  ") + n + (f"   [{d}]" if d and not c else ""), flush=True)


arcpy.env.overwriteOutput = True
SR = arcpy.SpatialReference(28350)

# ---- a project to start from: arcpy cannot make an .aprx from nothing, so any small project of
# yours is copied and its layers removed (the original is not touched). Set SGTOOL_TEST_APRX to
# the path of one.
TEMPLATE = os.environ.get("SGTOOL_TEST_APRX") or os.path.join(
    "C:" + os.sep, "Users", "00073294", "Dropbox", "cet_courses", "SG_courses_all", "#All_Labs",
    "Processing", "Prac_2_SGTool_examples.aprx")
if not os.path.exists(TEMPLATE):
    print("SKIPPED: set the environment variable SGTOOL_TEST_APRX to the path of a small .aprx "
          "project to run these checks")
    os._exit(0)
proj_path = os.path.join(DATA, "test_project.aprx")
arcpy.mp.ArcGISProject(TEMPLATE).saveACopy(proj_path)
aprx = arcpy.mp.ArcGISProject(proj_path)
the_map = aprx.listMaps()[0]
for lyr in list(the_map.listLayers()):
    try:
        the_map.removeLayer(lyr)
    except Exception:
        pass

# ---- the data: points, lines, two polygons in a geodatabase, and a grid
src_gdb = os.path.join(DATA, "source.gdb")
arcpy.management.CreateFileGDB(DATA, "source.gdb")


def make_fc(name, geom_type, geoms):
    fc = arcpy.management.CreateFeatureclass(src_gdb, name, geom_type, spatial_reference=SR).getOutput(0)
    arcpy.management.AddField(fc, "id", "LONG")
    with arcpy.da.InsertCursor(fc, ["SHAPE@", "id"]) as cur:
        for i, g in enumerate(geoms):
            cur.insertRow([g, i])
    return fc


xs_ = np.linspace(500000, 511000, 23)
ys_ = np.linspace(6889000, 6900000, 23)
pts = make_fc("my_points", "POINT", [arcpy.PointGeometry(arcpy.Point(float(x), float(y)), SR) for x in xs_ for y in ys_])
lines = make_fc("Faults", "POLYLINE", [arcpy.Polyline(arcpy.Array([arcpy.Point(500000, 6890000 + 600 * k),
                                                                   arcpy.Point(511000, 6891000 + 600 * k)]), SR)
                                       for k in range(10)])


def rect(x0, y0, x1, y1):
    return arcpy.Polygon(arcpy.Array([arcpy.Point(x0, y0), arcpy.Point(x0, y1), arcpy.Point(x1, y1),
                                      arcpy.Point(x1, y0), arcpy.Point(x0, y0)]), SR)


poly = make_fc("Clip_area", "POLYGON", [rect(502000, 6892000, 505000, 6895000), rect(507000, 6890000, 509000, 6892000)])

tif = os.path.join(DATA, "Magnetics.tif")
ds = gdal.GetDriverByName("GTiff").Create(tif, 100, 100, 1, gdal.GDT_Float32)
ds.SetGeoTransform([500000, 110, 0, 6900000, 0, -110])
s_ = osr.SpatialReference(); s_.ImportFromEPSG(28350)
ds.SetProjection(s_.ExportToWkt())
b_ = ds.GetRasterBand(1); b_.SetNoDataValue(-9999.0)
b_.WriteArray(np.random.default_rng(1).normal(0, 5, (100, 100)).astype("float32"))
ds = None

for path in (poly, lines, pts, tif):
    the_map.addDataFromPath(path)
by_name = {l.name: l for l in the_map.listLayers()}
check("the test map has the four layers", sorted(by_name) == ["Clip_area", "Faults", "Magnetics.tif", "my_points"], sorted(by_name))
# a symbol we can look for afterwards
sym = by_name["my_points"].symbology
sym.renderer.symbol.color = {"RGB": [255, 0, 0, 100]}
by_name["my_points"].symbology = sym
aprx.save()

# ---- the clip
out = os.path.join(DATA, "clipped")
from arcgis_clip import ClipJob
job = ClipJob(aprx, poly, out).run()
check("the vector layers and the grid were clipped, nothing failed",
      sorted(n for n, _ in job.clipped) == ["Clip_area", "Faults", "Magnetics.tif", "my_points"] and not job.errors,
      (job.clipped, job.errors))
gdb_path = os.path.join(out, "test_project_clip.gdb")
check("all the vectors are in one geodatabase called <project>_clip.gdb", os.path.isdir(gdb_path))
arcpy.env.workspace = gdb_path
fcs = sorted(arcpy.ListFeatureClasses())
check("feature classes are <name>_clip, the polygon (the mask used) saved as one of them",
      fcs == ["Clip_area_clip", "Faults_clip", "my_points_clip"], fcs)
n_expected = sum(1 for x in xs_ for y in ys_
                 if (502000 <= x <= 505000 and 6892000 <= y <= 6895000) or (507000 <= x <= 509000 and 6890000 <= y <= 6892000))
n_pts = int(arcpy.management.GetCount(os.path.join(gdb_path, "my_points_clip")).getOutput(0))
check("the clipped points are those inside the polygons", n_pts == n_expected, (n_pts, n_expected))
tifs = [f for f in os.listdir(out) if f.lower().endswith(".tif")]
check("the grid is written as <name>_clip.tif beside the geodatabase", tifs == ["Magnetics_clip.tif"], tifs)
ext = arcpy.Describe(os.path.join(out, tifs[0])).extent if tifs else None
check("the clipped grid is cropped to the polygons",
      ext is not None and ext.XMin >= 501890 and ext.XMax <= 509110 and ext.YMin >= 6889890 and ext.YMax <= 6895110,
      ext)
lyrx = sorted(f for f in os.listdir(out) if f.endswith(".lyrx"))
check("a .lyrx style file is saved for every layer", len(lyrx) == 4 and "Magnetics_clip.lyrx" in lyrx, lyrx)

# ---- the saved project
check("the project copy is saved beside them", job.aprx_path and os.path.exists(job.aprx_path)
      and os.path.basename(job.aprx_path) == "test_project_clip.aprx", job.aprx_path)
new = arcpy.mp.ArcGISProject(job.aprx_path)
nm = new.listMaps(the_map.name)[0]
names = sorted(l.name for l in nm.listLayers())
check("its layers are the clipped ones, renamed <oldname>_clip",
      names == ["Clip_area_clip", "Faults_clip", "Magnetics_clip", "my_points_clip"], names)
srcs = {l.name: l.dataSource for l in nm.listLayers()}
check("they point at the clipped data",
      all(os.path.normcase(os.path.dirname(v if v.lower().endswith(".tif") else os.path.dirname(v))).startswith(
          os.path.normcase(out)) for v in srcs.values()), srcs)
check("the original project is untouched",
      sorted(l.name for l in arcpy.mp.ArcGISProject(proj_path).listMaps()[0].listLayers()) == sorted(by_name))
pts_new = [l for l in nm.listLayers() if l.name == "my_points_clip"]
col = None
if pts_new and pts_new[0].symbology.renderer.type == "SimpleRenderer":
    col = pts_new[0].symbology.renderer.symbol.color
check("the point symbology (red) came through", col is not None and col["RGB"][:3] == [255, 0, 0], col)

# ---- cancel, and a failing layer
from calcs.sgt_cancel import OperationCancelled


def stop(i, n, name):
    if i >= 1:
        raise OperationCancelled()


try:
    ClipJob(aprx, poly, os.path.join(DATA, "cancelled"), check_cancel=stop).run()
    check("cancelling stops the clip", False)
except OperationCancelled:
    check("cancelling stops the clip", not os.path.exists(os.path.join(DATA, "cancelled", "test_project_clip.aprx")))

orig = ClipJob._clip_vector


def failing(self, layer, stem):
    if layer.name == "Faults":
        raise RuntimeError("simulated failure\nsecond line")
    return orig(self, layer, stem)


ClipJob._clip_vector = failing
try:
    job2 = ClipJob(aprx, poly, os.path.join(DATA, "failed")).run()
finally:
    ClipJob._clip_vector = orig
new2 = arcpy.mp.ArcGISProject(job2.aprx_path).listMaps(the_map.name)[0]
names2 = sorted(l.name for l in new2.listLayers())
check("a failing layer is reported and left out of the new project, the rest are clipped",
      [n for n, _ in job2.errors] == ["Faults"] and "Faults_clip" not in names2 and "Faults" not in names2
      and "my_points_clip" in names2, (job2.errors, names2))

# ---- the toolbox tool, with the open project standing in for "CURRENT"
loader = importlib.machinery.SourceFileLoader("sgt_pyt", os.path.join(REPO, "ARCGIS", "SGTool.pyt"))
spec = importlib.util.spec_from_loader("sgt_pyt", loader)
mod = importlib.util.module_from_spec(spec)
loader.exec_module(mod)
real_project = arcpy.mp.ArcGISProject
arcpy.mp.ArcGISProject = lambda p: real_project(proj_path if p == "CURRENT" else p)


class Msg:
    def __init__(self):
        self.lines = []

    def addMessage(self, m):
        self.lines.append(m)

    addWarningMessage = addErrorMessage = addMessage


try:
    tool = mod.ClipAllLayersToPolygon()
    params = tool.getParameterInfo()
    values = {"in_polygon": poly, "out_folder": os.path.join(DATA, "tool_out")}
    for p_ in params:
        if p_.name in values:
            p_.value = values[p_.name]
    msgs = Msg()
    tool.execute(params, msgs)
finally:
    arcpy.mp.ArcGISProject = real_project
check("the toolbox tool clips the project and saves its copy",
      os.path.exists(os.path.join(DATA, "tool_out", "test_project_clip.aprx"))
      and any("clipped" in m for m in msgs.lines), msgs.lines)
check("the tool is listed in the toolbox", mod.ClipAllLayersToPolygon in mod.Toolbox().tools)

print(f"\n{sum(1 for _, ok in results if ok)}/{len(results)} CLIP checks passed", flush=True)
os._exit(0 if all(ok for _, ok in results) else 1)
