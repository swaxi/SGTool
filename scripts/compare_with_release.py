#!/usr/bin/env python
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2024-2026 Mark Jessell
"""
Compare what the plugin code calculates now with what a release calculated.

Runs the functions the pysgtool work changed (since v0.4.8) from both versions
on the same synthetic inputs and reports, for each, whether the outputs are
identical, differ only by rounding, differ slightly, or differ:

  SAME       bit for bit
  CLOSE      max difference <= 1e-6 x the output's standard deviation (rounding)
  SMALL      max difference <= 1e-3 x the standard deviation
  DIFFERENT  anything larger, or different shapes, NaN cells, features, text
  ERROR      one version raised (the message says which)

The old version is taken from git (`git archive <ref>`), the new one from the
working tree, each into its own package so both can be imported at once. QGIS
itself is not needed: the qgis and Qt modules are replaced by stand-ins, and
the plugin methods that build QGIS layers are given small fakes that record
the features they add, which are then compared. GDAL is needed, so run it with
QGIS's Python:

  "C:\\Program Files\\QGIS 3.44.4\\bin\\python-qgis.bat" scripts\\compare_with_release.py
      [--ref 19e5093] [--grid my_grid.tif] [--out DIR]

--grid adds the Fourier-filter comparisons on one of your own grids (band 1).
The report is printed and also written to DIR/compare_report.md, with all
inputs and outputs of both versions kept under DIR for inspection.
"""
import argparse
import csv
import importlib
import importlib.abc
import importlib.machinery
import io
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import traceback
import types
import xml.etree.ElementTree as ET

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_REF = "19e5093"  # v0.4.8, the last release before the pysgtool work
EPSG = 28350
DX = 100.0
GT = [400000.0, DX, 0.0, 6500000.0, 0.0, -DX]
VERSIONS = ("old", "new")


# =========================================================================== QGIS stand-ins
STUB_ROOTS = {"qgis", "PyQt5", "PyQt6", "sip", "processing"}


class _StubMeta(type):
    def __getattr__(cls, name):
        if name.startswith("__"):
            raise AttributeError(name)
        return _StubMeta(name, (_Stub,), {})

    def __or__(cls, other):
        return cls

    __ror__ = __and__ = __rand__ = __or__

    def __invert__(cls):
        return cls

    def __int__(cls):  # e.g. [int(v) for v in QtCore.qVersion().split(".")] in resources.py
        return 0

    __index__ = __int__

    def __float__(cls):
        return 0.0


class _Stub(metaclass=_StubMeta):
    """Anything: every attribute, call, item or operator gives another stand-in."""

    def __init__(self, *args, **kwargs):
        pass

    def __getattr__(self, name):
        if name.startswith("__"):
            raise AttributeError(name)
        return _Stub()

    def __call__(self, *args, **kwargs):
        return _Stub()

    def __getitem__(self, key):
        return _Stub()

    def __iter__(self):  # unpacking, e.g. FORM_CLASS, _ = uic.loadUiType(...)
        return iter((_StubMeta("StubA", (_Stub,), {}), _StubMeta("StubB", (_Stub,), {})))

    def __len__(self):
        return 0

    def __bool__(self):
        return True

    def __or__(self, other):
        return self

    __ror__ = __and__ = __rand__ = __or__

    def __invert__(self):
        return self

    def __int__(self):
        return 0

    def __float__(self):
        return 0.0

    def __index__(self):
        return 0

    def __str__(self):
        return ""


class _StubModule(types.ModuleType):
    def __init__(self, name):
        super().__init__(name)
        self.__path__ = []
        self.__all__ = []
        self._cache = {}

    def __getattr__(self, name):
        if name.startswith("__") or name == "_cache":
            raise AttributeError(name)
        if name not in self._cache:
            self._cache[name] = _StubMeta(name, (_Stub,), {})
        return self._cache[name]


class _StubFinder(importlib.abc.MetaPathFinder, importlib.abc.Loader):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split(".")[0] in STUB_ROOTS:
            return importlib.machinery.ModuleSpec(fullname, self, is_package=True)
        return None

    def create_module(self, spec):
        return _StubModule(spec.name)

    def exec_module(self, module):
        pass


# =========================================================================== fakes that record
class FakeLayer:
    """A vector layer: records the fields and features the plugin adds."""

    def __init__(self, uri="", name="", provider="memory", log=None, features=None, fields=None):
        self.uri, self._name = uri, name
        self.log = log if log is not None else []
        self.added_fields, self.added = [], []
        self._features = features or []
        self._field_names = fields or []
        self.log.append(("layer", uri.split("?")[0] if "?" in uri else os.path.basename(uri), name))

    # writing
    def dataProvider(self):
        return self

    def addAttributes(self, fields):
        self.added_fields += [f.name() for f in fields]

    def addFeature(self, feature):
        self.log.append(("feature", self._name, feature.geom, list(feature.attrs)))
        return True

    def updateFields(self):
        pass

    def isValid(self):
        return True

    def id(self):
        return self._name

    def name(self):
        return self._name

    # reading (profiles)
    def fields(self):
        return FakeFields([FakeField(n) for n in self._field_names])

    def selectByExpression(self, expression):
        self.log.append(("select", expression))

    def selectedFeatures(self):
        return self._features


class FakeRaster:
    """QgsRasterLayer(path, name): records what was loaded."""

    def __init__(self, path, name, *args, log=None):
        self._name = name
        if log is not None:
            log.append(("raster layer", os.path.basename(path), name))

    def name(self):
        return self._name

    def id(self):
        return self._name


class FakeField:
    def __init__(self, name, *args):
        self._name = name

    def name(self):
        return self._name


class FakeFields(list):
    def append(self, field):
        super().append(field)

    def indexFromName(self, name):
        names = [f.name() for f in self]
        return names.index(name) if name in names else -1


class FakePoint:
    def __init__(self, x, y):
        self._x, self._y = x, y

    def x(self):
        return self._x

    def y(self):
        return self._y


class FakeGeometry:
    @staticmethod
    def fromPolylineXY(points):
        return ("line", tuple((p.x(), p.y()) for p in points))

    @staticmethod
    def fromPointXY(point):
        return ("point", (point.x(), point.y()))


class FakeFeature:
    def __init__(self, *args):
        self.geom, self.attrs, self._values = None, [], {}

    def setGeometry(self, geometry):
        self.geom = geometry

    def setAttributes(self, attributes):
        self.attrs = attributes

    def geometry(self):
        geom = self.geom

        class _G:
            def asPoint(self_inner):
                return FakePoint(*geom[1])
        return _G()

    def __getitem__(self, key):
        return self._values[key]


class FakeProject:
    """QgsProject.instance(): records added layers; serves the given layers by name."""

    def __init__(self, log, layers=None):
        self.log, self.layers = log, layers or {}

    def instance(self):
        return self

    def addMapLayer(self, layer):
        self.log.append(("add to project", layer._name if isinstance(layer, (FakeLayer, FakeRaster))
                         else type(layer).__name__))
        return layer

    def mapLayersByName(self, name):
        return [self.layers[name]] if name in self.layers else []

    def mapLayers(self):
        return {}

    def __getattr__(self, name):
        return _Stub()


# =========================================================================== versions
def find_git(given=None):
    """QGIS's python-qgis.bat resets PATH, so look in the usual places too."""
    candidates = [given, shutil.which("git"),
                  os.path.expandvars(r"%LOCALAPPDATA%\Programs\Git\cmd\git.exe"),
                  r"C:\Program Files\Git\cmd\git.exe", r"C:\Program Files (x86)\Git\cmd\git.exe"]
    for c in candidates:
        if c and os.path.exists(c):
            return c
    raise SystemExit("git not found: pass --git C:\\path\\to\\git.exe")


def build_versions(work, ref, git="git"):
    """work/old_sgt from git, work/new_sgt from the working tree."""
    os.makedirs(work, exist_ok=True)
    old = os.path.join(work, "old_sgt")
    shutil.rmtree(old, ignore_errors=True)
    tar = subprocess.run([git, "-C", REPO, "archive", "--format=tar", ref],
                         check=True, capture_output=True).stdout
    with tarfile.open(fileobj=io.BytesIO(tar)) as t:
        members = [m for m in t.getmembers()
                   if m.name.split("/")[0] not in {"ARCGIS", "docs", "tests", "pysgtool", "data", ".github"}]
        t.extractall(old, members=members)
    new = os.path.join(work, "new_sgt")
    shutil.rmtree(new, ignore_errors=True)
    shutil.copytree(REPO, new, ignore=shutil.ignore_patterns(
        ".git", "ARCGIS", "docs", "tests", "pysgtool", "pysgtool.egg-info", "data", "__pycache__",
        ".github", "*.tif", "*.zip"))
    sys.path.insert(0, work)
    sys.meta_path.insert(0, _StubFinder())


def mod(version, name):
    return importlib.import_module(f"{version}_sgt.{name}")


# =========================================================================== synthetic inputs
def synthetic_grid(ny=160, nx=200, seed=1):
    """Magnetic-looking anomalies on a regional plane, with noise (nT)."""
    rng = np.random.default_rng(seed)
    y, x = np.mgrid[0:ny, 0:nx].astype(float)
    g = 0.8 * x - 0.5 * y + 50.0
    for cy, cx, s, a in [(50, 60, 8, 300), (110, 140, 15, -200), (80, 100, 4, 150), (30, 170, 10, 120)]:
        r2 = (x - cx) ** 2 + (y - cy) ** 2
        g += a * np.exp(-r2 / (2 * s * s)) * (1 - (y - cy) / (2 * s))
    return g + rng.normal(0, 2.0, g.shape)


def with_holes(g):
    g = g.copy()
    ny, nx = g.shape
    y, x = np.mgrid[0:ny, 0:nx]
    g[x + y < 25] = np.nan
    g[((x - 40) / 12.0) ** 2 + ((y - 120) / 7.0) ** 2 < 1] = np.nan
    return g


def write_tif(path, bands, gt=GT, epsg=EPSG, nodata=None, dtype=None):
    from osgeo import gdal, osr
    bands = [np.asarray(b) for b in (bands if isinstance(bands, (list, tuple)) else [bands])]
    gtype = dtype if dtype is not None else gdal.GDT_Float32
    ds = gdal.GetDriverByName("GTiff").Create(path, bands[0].shape[1], bands[0].shape[0], len(bands), gtype)
    ds.SetGeoTransform(gt)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    ds.SetProjection(srs.ExportToWkt())
    for i, b in enumerate(bands, 1):
        band = ds.GetRasterBand(i)
        if nodata is not None:
            band.SetNoDataValue(nodata)
        band.WriteArray(b)
    ds = None
    return path


def write_noddy(path, values, code=333):
    ny, nx = values.shape
    with open(path, "w") as f:
        f.write(f"{code}\n")
        f.write(f"0 {nx} {ny} 1\n")
        f.write("-60.0 5.0 50000.0\n")
        f.write("0.0 0.0 5000.0\n")
        f.write(f"{nx * 50.0} {ny * 50.0} 0.0\n")
        f.write("50.0 80.0\n")
        f.write("# extra header line\n")
        f.write("# extra header line\n")
        for row in values:
            f.write(" ".join(f"{v:.4f}" for v in row) + "\n")
    return path


def write_xyz(path):
    rng = np.random.default_rng(3)
    lines = ["/ synthetic XYZ line data", "/ X Y MAG ALT"]
    for i, line_no in enumerate((10, 20, 30)):
        lines.append(f"LINE: {line_no}" if i != 1 else f"Line {line_no}")
        for k in range(25):
            x, y = 400000 + k * 100.0, 6500000 - i * 500.0
            lines.append(f"{x:.1f} {y:.1f} {100 + rng.normal(0, 30):.3f} {80 + k:.1f}")
        if i == 0:
            lines.append(f"{x + 100:.1f} {y:.1f} 2024-05-01 12:00:00 99.0 81.0")  # words among the numbers
    lines.append("TIE: 900")
    for k in range(10):
        lines.append(f"{400050 + k * 250.0:.1f} {6499000 + k * 100.0:.1f} {90.0 + k:.3f} {75.0:.1f}")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return path


def write_aseg(directory):
    dfn = os.path.join(directory, "survey.dfn")
    dat = os.path.join(directory, "survey.dat")
    with open(dfn, "w") as f:
        f.write("DEFN   ST=RECD,RT=COMM;RT:A4;COMMENTS:A76\n")
        f.write("DEFN 1 ST=RECD,RT=;LINE:I10:NAME=line\n")
        f.write("DEFN 2 ST=RECD,RT=;FID:F12.2:NAME=fid\n")
        f.write("DEFN 3 ST=RECD,RT=;EASTING:F12.2:UNIT=m,NAME=easting\n")
        f.write("DEFN 4 ST=RECD,RT=;NORTHING:F12.2:UNIT=m,NAME=northing\n")
        f.write("DEFN 5 ST=RECD,RT=;MAG:F10.2:UNIT=nT,NULL=-99999.00,NAME=mag\n")
        f.write("DEFN 6 ST=RECD,RT=;SPEC:4F8.1:NAME=spec\n")
        f.write("DEFN 7 ST=RECD,RT=;END DEFN\n")
    rng = np.random.default_rng(4)
    with open(dat, "w") as f:
        for line in (100, 200):
            for k in range(30):
                mag = -99999.0 if k == 7 else 1000 + rng.normal(0, 20)
                spec = " ".join(f"{v:8.1f}" for v in rng.uniform(0, 500, 4))
                f.write(f"{line:10d}{k:12.2f}{400000 + k * 50:12.2f}{6500000 + line:12.2f}{mag:10.2f}{spec}\n")
    return dfn, dat


# =========================================================================== comparing
def cmp_arrays(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.shape != b.shape:
        return "DIFFERENT", f"shape {a.shape} vs {b.shape}"
    na, nb = np.isnan(a), np.isnan(b)
    if (na != nb).any():
        return "DIFFERENT", f"NaN in different cells ({int((na != nb).sum())})"
    if np.array_equal(a, b, equal_nan=True):
        return "SAME", "identical"
    d = np.abs(a - b)[~na]
    scale = np.std(a[~na]) or np.max(np.abs(a[~na])) or 1.0
    rel = d.max() / scale
    status = "CLOSE" if rel <= 1e-6 else "SMALL" if rel <= 1e-3 else "DIFFERENT"
    return status, f"max |diff| {d.max():.3g} ({rel:.2g} x sd), rms {np.sqrt(np.mean(d ** 2)):.3g}"


WORST = ["SAME", "CLOSE", "SMALL", "DIFFERENT", "ERROR"]


def worst(*statuses):
    return max(statuses, key=WORST.index)


def cmp_values(a, b):
    """Nested lists/dicts/tuples/numbers/strings; numbers compared as arrays."""
    diffs = []

    def walk(x, y, where):
        if isinstance(x, dict) and isinstance(y, dict):
            for k in sorted(set(x) | set(y), key=str):
                if k not in x or k not in y:
                    diffs.append(("DIFFERENT", f"{where}/{k}: only in {'new' if k in y else 'old'}"))
                else:
                    walk(x[k], y[k], f"{where}/{k}")
        elif isinstance(x, (list, tuple)) and isinstance(y, (list, tuple)):
            if len(x) != len(y):
                diffs.append(("DIFFERENT", f"{where}: {len(x)} vs {len(y)} items"))
            for i, (p, q) in enumerate(zip(x, y)):
                walk(p, q, f"{where}[{i}]")
        elif isinstance(x, (int, float, np.number)) and isinstance(y, (int, float, np.number)) \
                and not isinstance(x, bool) and not isinstance(y, bool):
            s, _ = cmp_arrays([x], [y]) if x != y else ("SAME", "")
            if s != "SAME":
                rel = abs(x - y) / (abs(x) or 1.0)
                s = "CLOSE" if rel <= 1e-9 else "SMALL" if rel <= 1e-6 else "DIFFERENT"
                diffs.append((s, f"{where}: {x!r} vs {y!r}"))
        elif isinstance(x, np.ndarray) or isinstance(y, np.ndarray):
            s, d = cmp_arrays(x, y)
            if s != "SAME":
                diffs.append((s, f"{where}: {d}"))
        elif x != y:
            diffs.append(("DIFFERENT", f"{where}: {x!r} vs {y!r}"))

    walk(a, b, "")
    if not diffs:
        return "SAME", "identical"
    status = worst(*[s for s, _ in diffs])
    shown = [d for s, d in diffs if s == status][:3]
    more = f" (+{len(diffs) - len(shown)} more)" if len(diffs) > len(shown) else ""
    return status, "; ".join(shown) + more


def read_raster(path):
    from osgeo import gdal, osr
    ds = gdal.Open(path)
    if ds is None:
        raise FileNotFoundError(path)
    bands = np.array([ds.GetRasterBand(i + 1).ReadAsArray().astype(float) for i in range(ds.RasterCount)])
    nodata = ds.GetRasterBand(1).GetNoDataValue()
    if nodata is not None and not np.isnan(nodata):
        bands[bands == nodata] = np.nan
    srs = osr.SpatialReference(wkt=ds.GetProjection()) if ds.GetProjection() else None
    info = {"geotransform": list(ds.GetGeoTransform()),
            "epsg": srs.GetAuthorityCode(None) if srs else None,
            "metadata": {k: v for k, v in (ds.GetMetadata() or {}).items() if not k.startswith("SGTOOL")}}
    return bands, info


def cmp_rasters(pa, pb):
    a, ia = read_raster(pa)
    b, ib = read_raster(pb)
    s1, d1 = cmp_arrays(a, b)
    s2, d2 = cmp_values(ia, ib)
    status = worst(s1, s2)
    return status, d1 if s2 == "SAME" else f"{d1}; header: {d2}"


def read_vector(path):
    from osgeo import ogr
    ds = ogr.Open(path)
    if ds is None:
        raise FileNotFoundError(path)
    layer = ds.GetLayer(0)
    defn = layer.GetLayerDefn()
    names = [defn.GetFieldDefn(i).GetName() for i in range(defn.GetFieldCount())]
    feats = []
    for f in layer:
        g = f.GetGeometryRef()
        wkt = g.ExportToWkt() if g else ""
        coords = [float(v) for v in re.findall(r"-?\d+\.?\d*(?:[eE][-+]?\d+)?", wkt)]
        attrs = [f.GetField(n) for n in names]
        feats.append((wkt.split("(")[0].strip(), coords, attrs))
    feats.sort(key=lambda t: (t[0], [str(a) for a in t[2]], t[1][:4]))
    return {"fields": names, "count": len(feats),
            "types": [t for t, _, _ in feats], "attributes": [a for _, _, a in feats],
            "coords": np.array([c for _, cs, _ in feats for c in cs])}


def cmp_vectors(pa, pb):
    return cmp_values(read_vector(pa), read_vector(pb))


def cmp_csv(pa, pb):
    a = np.loadtxt(pa, delimiter=",", skiprows=1, ndmin=2)
    b = np.loadtxt(pb, delimiter=",", skiprows=1, ndmin=2)
    with open(pa) as f, open(pb) as g:
        header = "SAME" if f.readline() == g.readline() else "DIFFERENT"
    a = a[np.lexsort(a.T[::-1])] if a.size else a
    b = b[np.lexsort(b.T[::-1])] if b.size else b
    s, d = cmp_arrays(a, b)
    return worst(s, header), f"{len(a)} vs {len(b)} rows; {d}" + ("" if header == "SAME" else "; header differs")


VOLATILE = {"created", "modified", "mtime", "date_modified", "timestamp", "time"}


def canonical_xml(element, work):
    """Provenance XML as comparable text: no times, versions or version folders."""
    el = ET.fromstring(ET.tostring(element))
    for node in el.iter():
        for k in list(node.attrib):
            if k.lower() in VOLATILE or (node.tag == "software" and k == "version"):
                del node.attrib[k]
        if node.tag.lower() in VOLATILE:
            node.text = ""
    text = ET.tostring(el, encoding="unicode")
    for v in VERSIONS:
        for root in (os.path.join(work, v), os.path.join(work, v).replace("\\", "/")):
            text = text.replace(root, "<VERSION>")
    return re.sub(r"\s+", " ", text)


# =========================================================================== the comparisons
class Runner:
    def __init__(self, work):
        self.work = work
        self.rows = []

    def dir(self, version, *parts):
        d = os.path.join(self.work, version, *parts)
        os.makedirs(d, exist_ok=True)
        return d

    def record(self, area, test, status, detail):
        self.rows.append((area, test, status, detail))
        print(f"{status:9s} {area} | {test} | {detail}", flush=True)

    def check(self, area, test, run, compare):
        """run(version) -> output (anything); compare(old_output, new_output) -> (status, detail)."""
        outs = {}
        for v in VERSIONS:
            try:
                outs[v] = run(v)
            except Exception as e:  # noqa: BLE001  reported, not fatal
                tb = traceback.extract_tb(e.__traceback__)[-1]
                self.record(area, test, "ERROR",
                            f"{v} version raised {type(e).__name__}: {e} ({os.path.basename(tb.filename)}:{tb.lineno})")
                return
        try:
            status, detail = compare(outs["old"], outs["new"])
        except Exception as e:  # noqa: BLE001
            status, detail = "ERROR", f"comparison failed: {type(e).__name__}: {e}"
        self.record(area, test, status, detail)


def compare_fourier(r, grid, label, dx=DX):
    """Filters that remove the regional plane first (detrend.robust_plane)."""
    area = f"Fourier filters ({label})"
    holes = with_holes(grid) if label == "synthetic" else grid
    filled = np.where(np.isnan(holes), np.nanmedian(holes), holes)

    def plane(v):
        p, coef = mod(v, "calcs.detrend").robust_plane(filled)
        return p, np.array(coef)
    r.check(area, "detrend.robust_plane: plane grid", lambda v: plane(v)[0], cmp_arrays)

    cases = [
        ("upward continuation 500 m", lambda gp, g: gp.upward_continuation(g, 500.0)),
        ("downward continuation 50 m", lambda gp, g: gp.downward_continuation(g, 50.0)),
        ("RTP inc -60 dec 5", lambda gp, g: gp.reduction_to_pole(g, -60.0, 5.0)),
        ("RTE inc -60 dec 5", lambda gp, g: gp.reduction_to_equator(g, -60.0, 5.0)),
        ("RTP low latitude inc -8 dec 3, 3 iterations",
         lambda gp, g: gp.reduction_to_pole_low_latitude(g, -8.0, 3.0, iterations=3)),
        ("1st vertical derivative", lambda gp, g: gp.compute_derivative(g, "z", 1)),
        ("1st x derivative", lambda gp, g: gp.compute_derivative(g, "x", 1)),
        ("1st y derivative", lambda gp, g: gp.compute_derivative(g, "y", 1)),
        ("2nd vertical derivative", lambda gp, g: gp.compute_derivative(g, "z", 2)),
        ("tilt angle", lambda gp, g: gp.tilt_angle(g)),
        ("analytic signal", lambda gp, g: gp.analytic_signal(g)),
        ("total horizontal gradient", lambda gp, g: gp.total_horizontal_gradient(g)),
        ("low pass 2 km", lambda gp, g: gp.low_pass_filter(g, 2000.0)),
        ("high pass 2 km", lambda gp, g: gp.high_pass_filter(g, 2000.0)),
        ("band pass 1-5 km", lambda gp, g: gp.band_pass_filter(g, 5000.0, 1000.0)),
    ]
    def explained(run):
        """compare(old, new), plus whether the old code with the new plane fit gives the new result."""
        def compare(a, b):
            s, d = cmp_values(a, b)
            if s == "SAME":
                return s, d
            try:
                s2, d2 = cmp_values(with_new_plane(run), b)
                note = "all of it from the plane fit" if s2 == "SAME" else f"not only the plane fit: {d2}"
            except Exception as e:  # noqa: BLE001
                note = f"could not swap the plane fit in: {e}"
            return s, f"{d} [old code with the new plane fit: {note}]"
        return compare

    def taylor(v):
        tc = mod(v, "calcs.TaylorContinuation")
        ny, nx = holes.shape
        yy, xx = np.mgrid[0:ny, 0:nx]
        z_obs = 100.0 + 60.0 * np.sin(xx / 25.0) * np.cos(yy / 30.0)
        return tc.variable_continuation(holes.copy(), z_obs, 250.0, dx, dx)

    for name, f in cases:
        def run(v, f=f):
            GP = mod(v, "calcs.GeophysicalProcessor").GeophysicalProcessor
            return f(GP(dx, dx, 10), holes.copy())
        r.check(area, name, run, explained(run))
    r.check(area, "variable continuation (draped 40-160 m to 250 m)", taylor, explained(taylor))


def with_new_plane(run):
    """run("old") with the new detrend.robust_plane swapped into the old code."""
    modules = [mod("old", m) for m in ("calcs.detrend", "calcs.GeophysicalProcessor", "calcs.TaylorContinuation")]
    modules = [m for m in modules if hasattr(m, "robust_plane")]
    new = mod("new", "calcs.detrend").robust_plane
    saved = [m.robust_plane for m in modules]
    for m in modules:
        m.robust_plane = new
    try:
        return run("old")
    finally:
        for m, f in zip(modules, saved):
            m.robust_plane = f


class FakeRasterLayer:
    """What numpy_array_to_raster asks of a reference layer."""

    def __init__(self, shape, gt=GT, epsg=EPSG):
        self.shape, self.gt, self.epsg = shape, gt, epsg

    def dataProvider(self):
        return self

    def extent(self):
        ny, nx = self.shape
        gt = self.gt

        class _E:
            def width(s):
                return nx * gt[1]

            def height(s):
                return ny * -gt[5]

            def xMinimum(s):
                return gt[0]

            def yMaximum(s):
                return gt[3]
        return _E()

    def crs(self):
        epsg = self.epsg

        class _C:
            def authid(s):
                return f"EPSG:{epsg}"
        return _C()


def compare_worms(r, grid):
    area = "BSD worms"
    image = grid[:96, :128].copy()

    def pad(v, kind):
        W = mod(v, "calcs.worms.wormer").Wormer
        job = W()
        job.base_grid = image - image.mean()
        job.buildPaddedRaster((160, 192), rolloff_size=20, pad_type=kind)
        return job.padded_grid
    for kind in ("hann", "hamming", "cos"):  # the tapers the old FftUtils knew
        r.check(area, f"Wormer.buildPaddedRaster (FftUtils), taper {kind}", lambda v, k=kind: pad(v, k), cmp_arrays)

    def worms(v):
        d = r.dir(v, "worms")
        path = write_tif(os.path.join(d, "worm_grid.tif"), image)
        GP = mod(v, "calcs.GeophysicalProcessor").GeophysicalProcessor
        GP(DX, DX, 10).bsdwormer(image, FakeRasterLayer(image.shape), path, 4, 0.0, 400.0,
                                 True, EPSG, out_dir=d)
        return os.path.join(d, "worm_grid_worms.csv"), os.path.join(d, "worm_grid_worms.shp")
    outs = {}

    def run(v):
        outs[v] = worms(v)
        return outs[v]
    r.check(area, "bsdwormer: worm points (CSV), 4 levels", run, lambda a, b: cmp_csv(a[0], b[0]))
    if len(outs) == 2:
        r.check(area, "bsdwormer: worm polylines (shapefile)", lambda v: outs[v][1], cmp_vectors)


def compare_normalise(r, grid):
    area = "Normalise grids"
    for order, label in ((True, "first-order gradient"), (False, "second-order gradient")):
        def run(v, order=order, label=label):
            src = r.dir(v, "normalise_in")
            out = r.dir(v, f"normalise_out_{'1' if order else '2'}")
            for i, scale in enumerate((1.0, 0.1, 25.0)):
                g = with_holes(synthetic_grid(seed=10 + i)) * scale
                write_tif(os.path.join(src, f"grid{i}.tif"), g, nodata=-99999.0) if i != 1 else \
                    write_tif(os.path.join(src, f"grid{i}.tif"), np.where(np.isnan(g), -99999.0, g), nodata=-99999.0)
            GP = mod(v, "calcs.GeophysicalProcessor").GeophysicalProcessor
            GP(DX, DX, 10).normalise_geotiffs(src, out, order)
            return sorted(os.path.join(out, f) for f in os.listdir(out) if f.endswith(".tif"))

        def compare(a, b):
            if [os.path.basename(p) for p in a] != [os.path.basename(p) for p in b]:
                return "DIFFERENT", f"output files {[os.path.basename(p) for p in a]} vs {[os.path.basename(p) for p in b]}"
            res = [cmp_rasters(p, q) for p, q in zip(a, b)]
            return worst(*[s for s, _ in res]), " | ".join(f"{os.path.basename(p)}: {d}" for p, (_, d) in zip(a, res))
        r.check(area, f"normalise_geotiffs, {label}", run, compare)


def compare_boundary(r, grid):
    area = "Data boundary polygons"

    def run(v):
        d = r.dir(v, "boundary")
        path = write_tif(os.path.join(d, "holes.tif"), np.where(np.isnan(with_holes(grid)), -99999.0, grid),
                         nodata=-99999.0)
        out = os.path.join(d, "holes_boundary.shp")
        mod(v, "calcs.SG_Util").SG_Util(grid).create_data_boundary_lines(path, out)
        return out
    r.check(area, "SG_Util.create_data_boundary_lines (grid with a notch and a hole)", run, cmp_vectors)


def compare_aseg(r):
    area = "ASEG-GDF import"
    outs = {}

    def parse(v):
        d = r.dir(v, "aseg")
        dfn, dat = write_aseg(d)
        p = mod(v, "calcs.aseggdf2parser").AsegGdf2Parser()
        header, defs, epsg = p.parse_dfn_file(dfn)
        data = p.parse_dat_file(dat, header, defs)
        outs[v] = (p, d, data, defs)
        return {"header": header, "field_defs": defs, "epsg": epsg, "rows": data}
    r.check(area, "parse_dfn_file + parse_dat_file", parse, cmp_values)
    if len(outs) == 2:
        def export(v):
            p, d, data, defs = outs[v]
            path = os.path.join(d, "survey.shp")
            p.export_to_shapefile(data, defs, path, "survey", "easting", "northing", "28350")
            return path
        r.check(area, "export_to_shapefile", export, cmp_vectors)


def plugin(v, log=None, layers=None):
    """An SGTool object without QGIS: its __init__ is skipped, QGIS layer classes faked."""
    m = mod(v, "SGTool")
    log = [] if log is None else log
    m.QgsProject = FakeProject(log, layers)
    m.QgsVectorLayer = lambda *a, **k: FakeLayer(*a, log=log)
    m.QgsFields = FakeFields
    m.QgsField = FakeField
    m.QgsPointXY = FakePoint
    m.QgsGeometry = FakeGeometry
    m.QgsFeature = FakeFeature
    m.QgsRasterLayer = lambda path, name, *a: FakeRaster(path, name, log=log)
    obj = object.__new__(m.SGTool)
    obj.iface = _Stub()
    obj.dlg = _Stub()
    obj.dlg.mFieldComboBox_feature = types.SimpleNamespace(currentText=lambda: "1")  # the line chosen
    obj.save_layer_as_shapefile = lambda layer, path: False
    return obj, log


def compare_plugin_imports(r, grid):
    area = "Plugin imports (SGTool.py)"

    # Noddy grids
    for code, kind in ((333, "magnetic"), (444, "gravity")):
        def noddy(v, code=code):
            d = r.dir(v, f"noddy_{code}")
            src = write_noddy(os.path.join(d, f"model_{code}.mag" if code == 333 else f"model_{code}.grv"),
                              synthetic_grid(60, 80, seed=code))
            obj, log = plugin(v)
            obj.loadNoddyGrid(src)
            out = [os.path.join(d, f) for f in os.listdir(d) if f.endswith(".tif")]
            return out, log

        def compare(a, b):
            if len(a[0]) != 1 or len(b[0]) != 1:
                return "DIFFERENT", f"GeoTIFFs written: {len(a[0])} vs {len(b[0])}"
            s1, d1 = cmp_rasters(a[0][0], b[0][0])
            s2, d2 = cmp_values(a[1], b[1])
            return worst(s1, s2), d1 + ("" if s2 == "SAME" else f"; QGIS calls: {d2}")
        r.check(area, f"loadNoddyGrid ({kind}, file code {code})", noddy, compare)

    # XYZ and CSV headers
    def xyz_path(v):
        return write_xyz(os.path.join(r.dir(v, "xyz"), "lines.xyz"))
    r.check(area, "get_XYZ_header (values per row)", lambda v: plugin(v)[0].get_XYZ_header(xyz_path(v)), cmp_values)

    def plain_xyz(v):
        p = os.path.join(r.dir(v, "xyz"), "plain.xyz")
        with open(p, "w") as f:
            f.write("LINE: 5\n1.0 2.0 3.0 4.0\n1.5 2.5 3.5 4.5\n")
        return plugin(v)[0].get_XYZ_header(p)
    r.check(area, "get_XYZ_header (plain LINE: file)", plain_xyz, cmp_values)

    def header(v):
        p = os.path.join(r.dir(v, "csv"), "points.csv")
        with open(p, "w") as f:
            f.write("X,Y,,MAG,GRAV\n1,2,,3,4\n")
        return plugin(v)[0].read_header(p)
    r.check(area, "read_header (CSV)", header, cmp_values)

    for ties in (True, False):
        def xyz(v, ties=ties):
            obj, log = plugin(v)
            obj.import_XYZ(xyz_path(v), "28350", layer_name="lines", load_ties=ties)
            return log
        r.check(area, f"import_XYZ, tie lines {'loaded' if ties else 'left out'} (layers, features, attributes)",
                xyz, cmp_values)

    # RGB to grey
    def rgb(v):
        d = r.dir(v, "rgb")
        cmap_vals = (grid - grid.min()) / (grid.max() - grid.min())
        red = np.clip(255 * cmap_vals, 0, 255)
        green = np.clip(255 * (1 - np.abs(cmap_vals - 0.5) * 2), 0, 255)
        blue = np.clip(255 * (1 - cmap_vals), 0, 255)
        red[:5, :5] = green[:5, :5] = blue[:5, :5] = 255  # white: NaN
        red[-5:, -5:] = green[-5:, -5:] = blue[-5:, -5:] = 0  # black: NaN
        from osgeo import gdal
        path = write_tif(os.path.join(d, "image.tif"), [red, green, blue], dtype=gdal.GDT_Byte)
        obj, _ = plugin(v)
        ok, out = obj.convert_RGB_to_grey(path, "#0000ff, #00ff00, #ff0000", lut_min=-100.0, lut_max=400.0)
        if not ok:
            raise RuntimeError(f"convert_RGB_to_grey returned {ok}, {out}")
        return out
    r.check(area, "convert_RGB_to_grey (3-colour table, scaled -100 to 400)", rgb, cmp_rasters)

    # profiles
    rng = np.random.default_rng(7)
    xs = np.sort(rng.uniform(0, 5000, 40)) + 400000
    ys = 6500000 + 0.3 * (xs - 400000) + rng.normal(0, 20, 40)
    vals = 100 + 50 * np.sin((xs - 400000) / 700.0)
    order = rng.permutation(40)

    def profile_layer(log):
        feats = []
        for i in order:
            f = FakeFeature()
            f.geom = ("point", (xs[i], ys[i]))
            f._values = {"MAG": vals[i]}
            feats.append(f)
        return FakeLayer("Point", "pts", log=log, features=feats, fields=["LINE_ID", "MAG"])

    for spacing, num in (("auto", None), (150.0, None), (None, 25), (None, None)):
        def prof(v, spacing=spacing, num=num):
            log = []
            obj, _ = plugin(v, log, {"pts": profile_layer(log)})
            return obj.regularize_selected_points("pts", "MAG", spacing, num)
        label = f"spacing {spacing}" if spacing is not None else f"{num} points" if num else "original count"
        r.check(area, f"regularize_selected_points ({label})", prof, cmp_values)


def compare_replay(r):
    area = "Replay and filter table"

    def specs(v):
        table = mod(v, "calcs.filter_table") if v == "new" else mod(v, "sgt_filters")

        def plain(x):
            if isinstance(x, dict):
                return {k: plain(y) for k, y in x.items() if not callable(y)}
            if isinstance(x, (list, tuple)):
                return [plain(y) for y in x]
            return None if callable(x) else x
        return {"specs": plain(table.FILTER_SPECS), "suffixes": plain(table.SUFFIXES)}
    r.check(area, "FILTER_SPECS and SUFFIXES (everything but the functions)", specs, cmp_values)

    def plans(v, kind):
        wf = mod(v, "calcs.replay") if v == "new" else mod(v, "sgt_workflow")
        aniso = [o for o, s in wf._OPERATION_TO_SPEC.items() if s == "local_anisotropy"]
        if kind == "all":
            steps = [{"operation": op, "parameters": {}, "n_sources": 1, "software": "SGTool",
                      "output": "grid_out.tif"} for op in sorted(wf._OPERATION_TO_SPEC)]
            steps += [{"operation": op, "parameters": {}, "n_sources": 1, "output": "grid_AnisoOrient.tif"}
                      for op in aniso]
        else:  # the step recorded its output type, and the file name says the other one
            steps = [{"operation": op, "parameters": {"output_type": "Magnitude"}, "n_sources": 1,
                      "output": "grid_AnisoOrient.tif"} for op in aniso]
        return [{k: (str(x) if k == "date" else x) for k, x in item.items()} for item in wf.plan_steps(steps)]
    r.check(area, "plan_steps for every replayable operation, no parameters recorded",
            lambda v: plans(v, "all"), cmp_values)
    r.check(area, "plan_steps, local anisotropy with OUTPUT_TYPE recorded "
                  "(intended change: the recorded type now wins over the file name)",
            lambda v: plans(v, "recorded"), cmp_values)


def compare_provenance(r, grid):
    area = "Provenance"

    def chain(v):
        d = r.dir(v, "provenance")
        md = mod(v, "calcs.sgt_metadata")
        a = write_tif(os.path.join(d, "a.tif"), grid)
        b = write_tif(os.path.join(d, "b.tif"), grid + 1)
        c = write_tif(os.path.join(d, "c.tif"), grid + 2)
        md.write_sgt_metadata(b, a, "Upward continuation", {"height": 500}, "test")
        md.write_sgt_metadata(c, b, "Vertical derivative", {"order": 1}, "test")
        return canonical_xml(md.read_sgt_metadata(c), r.work)

    def compare(a, b):
        if a == b:
            return "SAME", "identical apart from times and versions"
        i = next(k for k in range(min(len(a), len(b))) if a[k] != b[k]) if a[:min(len(a), len(b))] != b[:min(len(a), len(b))] else min(len(a), len(b))
        return "DIFFERENT", f"first difference: ...{a[max(0, i - 60):i + 60]}... vs ...{b[max(0, i - 60):i + 60]}..."
    r.check(area, "write_sgt_metadata chain a -> b -> c, read back from c", chain, compare)


# =========================================================================== main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", default=DEFAULT_REF, help="git ref of the release to compare with")
    ap.add_argument("--grid", help="also run the Fourier comparisons on this grid (band 1)")
    ap.add_argument("--out", help="folder for inputs, outputs and the report")
    ap.add_argument("--git", help="path to git.exe, if it is not found")
    args = ap.parse_args()

    work = os.path.abspath(args.out or os.path.join(tempfile.gettempdir(), f"sgtool_compare_{args.ref}"))
    print(f"comparing {args.ref} with the working tree of {REPO}\nfiles in {work}\n", flush=True)
    build_versions(work, args.ref, find_git(args.git))
    import logging
    logging.disable(logging.WARNING)
    for v in VERSIONS:
        try:
            mod(v, "SGTool")
        except Exception as e:  # noqa: BLE001
            print(f"note: could not import the {v} SGTool.py ({type(e).__name__}: {e}); plugin imports will error")

    r = Runner(work)
    grid = synthetic_grid()
    compare_fourier(r, grid, "synthetic")
    if args.grid:
        g, info = read_raster(args.grid)
        compare_fourier(r, g[0], os.path.basename(args.grid), dx=abs(info["geotransform"][1]))
    compare_worms(r, grid)
    compare_normalise(r, grid)
    compare_boundary(r, grid)
    compare_aseg(r)
    compare_plugin_imports(r, grid)
    compare_replay(r)
    compare_provenance(r, grid)

    counts = {s: sum(1 for row in r.rows if row[2] == s) for s in WORST}
    summary = ", ".join(f"{n} {s}" for s, n in counts.items() if n)
    lines = [f"# SGTool: {args.ref} compared with the working tree", "", summary, "",
             "| Status | Area | Test | Detail |", "|---|---|---|---|"]
    lines += [f"| {s} | {a} | {t} | {d.replace('|', '/')} |" for a, t, s, d in r.rows]
    report = os.path.join(work, "compare_report.md")
    with open(report, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\n{summary}\nreport: {report}")
    return 0 if counts["DIFFERENT"] == 0 and counts["ERROR"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
