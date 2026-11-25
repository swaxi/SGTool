"""
Tests of the pysgtool package, on synthetic grids.

    pip install -e .[all,test]
    pytest tests

Tests of optional parts are skipped where what they need is not installed, so the
tests of the array calculations also run with numpy and scipy only. PYSGTOOL_IO=gdal or
rasterio chooses what reads and writes the files.
"""

import datetime
import importlib.util
import os
import subprocess
import sys

import numpy as np
import pytest

import pysgtool as sg
from pysgtool.calcs import raster_io
from pysgtool.calcs.layer_shim import LayerShim
from pysgtool.calcs.filter_table import FILTER_SPECS, GridContext


def has(*modules):
    return all(importlib.util.find_spec(m) is not None for m in modules)


HAS_IO = has("osgeo") or has("rasterio")
HAS_CRS = has("osgeo") or has("pyproj")
needs_io = pytest.mark.skipif(not HAS_IO, reason="needs GDAL or rasterio")
needs_crs = pytest.mark.skipif(not HAS_CRS, reason="needs GDAL or pyproj")
needs_sklearn = pytest.mark.skipif(not has("sklearn", "shapely"), reason="needs scikit-learn, shapely")
needs_mpl = pytest.mark.skipif(not has("matplotlib"), reason="needs matplotlib")
CRS = "EPSG:32750" if HAS_CRS else None
ORIGIN = (500000.0, 7000000.0)


def make(nan=True, crs=CRS, shape=(90, 110), seed=1):
    rng = np.random.default_rng(seed)
    y, x = np.mgrid[0:shape[0], 0:shape[1]] * 100.0
    a = (300 * np.exp(-((x - 4000) ** 2 + (y - 3500) ** 2) / 1.2e6)
         - 200 * np.exp(-((x - 7500) ** 2 + (y - 6000) ** 2) / 5e5)
         + 0.01 * x + rng.normal(0, 2, x.shape) + 400)
    if nan:
        a[10:16, 20:28] = np.nan
    return sg.Grid(a, 100, origin=ORIGIN, crs=crs, name="mag")


def settings(spec_id):
    if spec_id in ("reduction_to_pole", "reduction_to_equator", "reduction_to_pole_low_latitude"):
        return {"inclination": -60.0 if "low" not in spec_id else 5.0, "declination": 0.0}
    return {"field": 50000.0} if spec_id == "apparent_susceptibility" else {}


# ---------------------------------------------------------------- filters
@pytest.mark.parametrize("spec", FILTER_SPECS, ids=lambda s: s["id"])
def test_filter_is_the_plugins(spec):
    g = make(crs=None)
    before = g.values.copy()
    out = getattr(g, spec["id"])(**settings(spec["id"]))
    values = {p["name"]: p["options"][p["default"]] if p["kind"] == "enum" else p["default"]
              for p in spec["params"]}
    ctx = GridContext(g.values.copy(), g.gt, LayerShim(g.extent, "", False), 0)
    np.testing.assert_array_equal(out.values, np.asarray(spec["compute"](ctx, values), dtype=float))
    np.testing.assert_array_equal(g.values, before)  # its own grid is left alone
    assert out.gt == g.gt and out.history[-1]["operation"] == spec["operation"]


def test_settings_are_checked():
    g = make(crs=None)
    for call in (
        lambda: g.derivative(direction="q"), lambda: g.derivative(power=-1),
        lambda: g.derivative(power="two"), lambda: g.mean_filter(size=3.5),
        lambda: g.mean_filter(size=1), lambda: g.reduction_to_pole(),
        lambda: g.reduction_to_pole(inclination=-60), lambda: g.apparent_susceptibility(),
        lambda: g.sun_shading(relief="yes"), lambda: g.gaussian_filter(sigma=float("nan")),
    ):
        with pytest.raises(ValueError):
            call()
    with pytest.raises(TypeError):
        g.derivative(directoin="z")
    assert g.derivative("X", 2).name == "mag_d2x"


def test_filter_list():
    f = sg.filters()
    assert set(f) == {s["id"] for s in FILTER_SPECS}
    assert f["reduction_to_pole"]["settings"]["inclination"] is None
    assert "direction" in str(sg.Grid.derivative.__signature__)


def test_vertical_derivative_of_a_sine():
    """d/dz of a sinusoid of wavenumber k is k times it (z up, field above the sources)."""
    n, k = 256, 2 * np.pi / 3200.0
    x = np.arange(n) * 100.0
    g = sg.Grid(np.tile(np.sin(k * x), (n, 1)), 100)
    inner = (slice(80, 176),) * 2
    up = g.continuation("up", 500).values[inner]
    np.testing.assert_allclose(up, np.exp(-k * 500) * g.values[inner], atol=0.02)
    dz = g.derivative("z", 1).values[inner]
    np.testing.assert_allclose(np.abs(dz), np.abs(k * g.values[inner]), atol=0.03 * k)


# ---------------------------------------------------------------- Grid
def test_grid_geometry():
    g = make(crs=None)
    assert g.shape == (90, 110) and g.extent == (500000.0, 6991000.0, 511000.0, 7000000.0)
    assert g.x[0] == 500050 and g.y[0] == 6999950 and g.y[-1] == 6991050
    assert g.epsg is None and not g.geographic
    assert sg.Grid(np.ones((4, 5)), 1, geographic=True).geographic
    with pytest.raises(ValueError):
        sg.Grid(np.ones(5), 1)
    with pytest.raises(ValueError):
        sg.Grid(np.full((3, 3), np.nan), 1).mean_filter()
    with pytest.raises(ValueError):
        sg.Grid(np.ones((3, 3)), gt=(0, 1, 0, 0, 0, 1))
    with pytest.raises(ValueError):
        g.save()  # nowhere to save it


@needs_crs
def test_crs():
    g = make()
    assert g.epsg == 32750 and not g.geographic
    assert make(crs=32750).epsg == 32750
    assert make(crs="EPSG:4326").geographic
    with pytest.raises(ValueError):
        make(crs="EPSG:999999")
    with pytest.raises(sg.FilterError):
        make(crs="EPSG:4326").vertical_integration()
    with pytest.raises(ValueError):
        make(crs="EPSG:4326").euler_deconvolution()


def test_history_in_memory():
    g = make(crs=None)
    out = g.reduction_to_pole(-60, 2).local_anisotropy(output_type="Orientation").mean_filter(5)
    steps = out.history
    assert [s["operation"] for s in steps] == ["Reduction to pole", "Local anisotropy", "Mean filter"]
    assert steps[0]["parameters"]["DECLINATION"] == "2.0"
    assert g.history == [] and len(g.mean_filter().history) == 1
    from pysgtool.calcs.replay import plan_steps

    plan = plan_steps(steps)
    assert all(p["replayable"] for p in plan)
    assert plan[1]["values"]["OUTPUT_TYPE"] == "Orientation"


# ---------------------------------------------------------------- files
@needs_io
def test_save_and_read(tmp_path):
    g = make()
    path = tmp_path / "mag.tif"
    g.save(path)
    r = sg.read(path)
    np.testing.assert_allclose(r.values, g.values, rtol=1e-6)
    assert r.gt == g.gt and r.epsg == 32750 and r.name == "mag" and r.history == []
    out = r.reduction_to_pole(inclination=-60, declination=2).derivative("z", 1).save()
    assert os.path.basename(out.path) == "mag_RTP_d1z.tif"
    assert not (tmp_path / "mag_RTP.tif").exists()
    back = sg.read(out.path)
    assert [s["operation"] for s in back.history] == ["Reduction to pole", "Derivative"]
    assert os.path.basename(back.history[0]["source"]) == "mag.tif"
    assert [s["operation"] for s in sg.history(out.path)] == ["Reduction to pole", "Derivative"]
    # going on from a saved step, and over an earlier file
    assert len(back.tilt_angle().save().history) == 3
    assert len(sg.read(path).reduction_to_pole(-60, 2).derivative("z", 1).save().history) == 2
    with pytest.raises(ValueError):
        sg.read(path).save()  # never over the file it was read from, by default
    with pytest.raises(ValueError):
        out.save(tmp_path / "x.ers")
    with pytest.raises(FileNotFoundError):
        sg.read(tmp_path / "none.tif")
    with pytest.raises(ValueError):
        sg.read(path, band=2)


@needs_io
def test_south_up_and_nodata(tmp_path):
    g = make(nan=False)
    path = str(tmp_path / "south.tif")
    south = (g.gt[0], 100.0, 0.0, g.extent[1], 0.0, 100.0)
    a = g.values[::-1].astype(np.float32)
    a[0, 0] = -9999
    if raster_io.backend() == "gdal":
        from osgeo import gdal

        ds = gdal.GetDriverByName("GTiff").Create(path, 110, 90, 1, gdal.GDT_Float32)
        ds.SetGeoTransform(south)
        ds.GetRasterBand(1).SetNoDataValue(-9999)
        ds.GetRasterBand(1).WriteArray(a)
        ds = None
    else:
        import rasterio

        with rasterio.open(path, "w", driver="GTiff", width=110, height=90, count=1,
                           dtype="float32", nodata=-9999,
                           transform=rasterio.Affine.from_gdal(*south)) as ds:
            ds.write(a, 1)
    r = sg.read(path)
    assert r.gt == g.gt and np.isnan(r.values[-1, 0]) and np.isnan(r.values).sum() == 1
    np.testing.assert_allclose(r.values[:-1], g.values[:-1], rtol=1e-6)


@needs_io
def test_stack_files(tmp_path):
    g = make(nan=False)
    s = sg.Stack.from_grids([g, g.mean_filter(), g.gaussian_filter()])
    s.names = ["K_pct", "Th", "eU"]
    s.save(tmp_path / "s.tif")
    r = sg.read_stack(tmp_path / "s.tif")
    assert r.names == s.names and len(r) == 3
    np.testing.assert_allclose(r.values, s.values, rtol=1e-6)
    assert r.band("Th").name.endswith("_b2") and r.band(2).shape == g.shape
    assert sg.read_stack(tmp_path / "s.tif", [3, 1]).names == ["eU", "K_pct"]
    with pytest.raises(ValueError):
        r.band(4)


def write_geosoft(path, values, dx, dy, x0, y0, order):
    """A minimal uncompressed double Geosoft grid; values row 0 = south."""
    rows, cols = values.shape
    by_row = order == 1
    ne, nv = (cols, rows) if by_row else (rows, cols)
    de, dv = (dx, dy) if by_row else (dy, dx)
    h = bytearray(512)
    h[0:20] = np.array([8, 2, ne, nv, order], dtype="<i4").tobytes()
    h[20:60] = np.array([de, dv, x0, y0, 0], dtype="<f8").tobytes()
    h[60:76] = np.array([0, 1], dtype="<f8").tobytes()
    data = np.where(np.isnan(values), -1e32, values)
    with open(path, "wb") as f:
        f.write(bytes(h))
        f.write(data.astype("<f8").tobytes(order="C" if by_row else "F"))


@pytest.mark.parametrize("order", [1, -1])
def test_geosoft(tmp_path, order):
    values = np.arange(12.0).reshape(3, 4)  # row 0 = south
    path = str(tmp_path / "g.grd")
    write_geosoft(path, values, 50.0, 80.0, 1000.0, 2000.0, order)
    g = sg.read(path)
    np.testing.assert_array_equal(g.values, values[::-1])
    assert g.gt == (975.0, 50.0, 0.0, 2200.0, 0.0, -80.0)
    assert g.y[-1] == 2000 and g.x[0] == 1000


# ---------------------------------------------------------------- resampling
@needs_io
def test_resample():
    g = make(nan=False)
    frame = sg.Grid(np.zeros((30, 40)), 230, origin=(ORIGIN[0] + 130, ORIGIN[1] - 170), crs=CRS)
    r = g.resample_to(frame)
    assert r.shape == (30, 40) and r.gt == frame.gt
    plane = sg.Grid(np.add.outer(g.y, 2 * g.x), gt=g.gt, crs=CRS)
    p = plane.resample_to(sg.Grid(np.zeros((50, 50)), 60, origin=(ORIGIN[0] + 730, ORIGIN[1] - 870),
                                  crs=CRS))
    np.testing.assert_allclose(p.values, np.add.outer(p.y, 2 * p.x), atol=1.0)  # float32
    with pytest.raises(ValueError):
        g.resample_to(sg.Grid(np.zeros((5, 5)), 100, origin=(0, 0), crs=CRS))
    with pytest.raises(ValueError):
        g.resample_to(make(crs=None))


@needs_io
def test_reproject():
    x = np.linspace(116.0, 116.2, 81)
    g = sg.Grid(np.tile(x, (81, 1)), 0.0025, origin=(115.99875, -27.0), crs="EPSG:4326")
    to = sg.Grid(np.zeros((40, 40)), 200, origin=(403000, 7005000), crs="EPSG:32750")
    r = g.resample_to(to)
    from pysgtool.calcs.igrf.igrf_field import _lonlat_arrays

    X, Y = np.meshgrid(to.x, to.y)
    lon, _lat = _lonlat_arrays("EPSG:32750", X, Y)
    ok = np.isfinite(r.values)
    assert ok.sum() > 400
    np.testing.assert_allclose(r.values[ok], lon[ok], atol=3e-4)


# ---------------------------------------------------------------- field
def test_igrf_values():
    inc, dec, f = sg.igrf(-32.0, 116.0, 2020.0)  # Perth
    assert -67 < inc < -64 and -3 < dec < 0 and 57000 < f < 60000
    assert sg.igrf(-32.0, 116.0, "2020-01-01") == (inc, dec, f)
    assert np.shape(sg.igrf(np.zeros(3), np.arange(3.0), 2020.0)[0]) == (3,)


@needs_crs
def test_igrf_grids():
    g = make()
    date = datetime.date(2015, 7, 2)
    inc, dec, f = g.field_direction(date)
    gi, gd, gf = g.igrf("2015-07-02")
    assert abs(gi.values[45, 55] - inc) < 0.01 and abs(gf.values[45, 55] - f) < 1
    r = g.remove_igrf(date)
    np.testing.assert_allclose(r.values, g.values - gf.values)
    assert r.name == "mag_RMI"
    a = g.reduction_to_pole(date=date)
    np.testing.assert_array_equal(a.values, g.reduction_to_pole(inc, dec).values)
    assert a.history[0]["parameters"]["igrf_date"] == "2015-07-02"
    d = g.differential_rtp(date)
    assert np.isnan(d.values).sum() == np.isnan(g.values).sum()
    y, x = np.mgrid[0:90, 0:110] * 100.0  # anomalies on a zero background, away from the edges
    c = sg.Grid(300 * np.exp(-((x - 5000) ** 2 + (y - 4500) ** 2) / 1.2e6), gt=g.gt, crs=CRS)
    inner = (slice(15, 75), slice(15, 95))
    both = [m.values[inner].ravel() for m in (c.differential_rtp(date), c.reduction_to_pole(inc, dec))]
    assert np.corrcoef(*both)[0, 1] > 0.99
    with pytest.raises(ValueError):
        make(crs=None).remove_igrf(date)


def test_noddy_field():
    g = make(crs=None)
    g.metadata = {"inclination": "-50", "declination": "3", "intensity": "51000"}
    a = g.reduction_to_pole(date="2000-01-01")
    np.testing.assert_array_equal(a.values, g.reduction_to_pole(-50, 3).values)
    assert a.field_direction("2000-01-01") == (-50, 3, 51000)
    with pytest.raises(ValueError):
        g.remove_igrf("2000-01-01")


# ---------------------------------------------------------------- continuation
def test_variable_continuation_flat():
    g = make(nan=False, crs=None)
    up = g.variable_continuation(100, 400)
    inner = (slice(25, 65), slice(25, 85))
    ref = g.continuation("up", 300).values
    assert np.corrcoef(up.values[inner].ravel(), ref[inner].ravel())[0, 1] > 0.999
    assert up.history[0]["parameters"]["observed"] == "100 m" and "iterations" in up.info


@needs_io
def test_variable_continuation_surfaces():
    g = make(nan=False)
    dem = sg.Grid(np.add.outer(np.zeros(90), np.linspace(0, 100, 110)), gt=g.gt, crs=CRS, name="dem")
    out = g.variable_continuation((dem, 80), 400.0)
    assert out.history[0]["parameters"]["observed"] == "dem + 80 m"
    with pytest.warns(UserWarning):
        g.variable_continuation((dem, 80), 50.0)
    with pytest.raises(ValueError):
        g.variable_continuation((80, dem), 400)


# ---------------------------------------------------------------- Euler, worms
def test_euler(tmp_path):
    n, dx, z0 = 80, 50.0, 300.0
    c = (np.arange(n) - n / 2 + 0.5) * dx
    X, Y = np.meshgrid(c, c)
    # vertical field of a point mass: structural index 2
    f = 1e9 * z0 / (X ** 2 + Y ** 2 + z0 ** 2) ** 1.5
    g = sg.Grid(f, dx, origin=(1000.0, 5000.0), name="pt")
    s = g.euler_deconvolution(2, window=10, keep=0.02)
    assert len(s) == int(0.02 * 70 * 70)
    assert abs(np.median(s.x) - 3000) < dx and abs(np.median(s.y) - 3000) < dx
    assert abs(np.median(s.depth) - z0) < dx
    path = s.save(tmp_path / "e.csv")
    table = np.loadtxt(path, delimiter=",", skiprows=1)
    np.testing.assert_array_equal(table[:, 0], s.y)
    np.testing.assert_array_equal(table[:, 2], -s.depth)
    assert os.path.exists(path + ".sgt.xml") and sg.history(path)[-1]["operation"] == "Euler deconvolution"
    stats = g.euler_window_statistics(tmp_path / "st", keep=0.2, structural_indices=(0, 2))
    assert set(stats) == {0, 2} and (tmp_path / "st" / "pt_window_summary.txt").exists()
    with pytest.raises(ValueError):
        g.euler_deconvolution(2, window=100)


def test_worms_points(tmp_path):
    g = make(crs=None)
    w = g.worms(3, 200, 200)
    assert len(w) and set(w.z) == {200, 400, 600}
    xmin, ymin, xmax, ymax = g.extent
    assert w.x.min() >= xmin and w.x.max() <= xmax and w.y.min() >= ymin and w.y.max() <= ymax
    table = np.genfromtxt(w.save(tmp_path / "w.csv"), delimiter=",", skip_header=1)
    np.testing.assert_array_equal(table, w.table)


@pytest.mark.skipif(not has("osgeo"), reason="needs GDAL")
def test_worms_are_the_plugins(tmp_path):
    g = make()
    g.save(tmp_path / "g.tif")
    r = sg.read(tmp_path / "g.tif")
    from pysgtool.calcs.GeophysicalProcessor import GeophysicalProcessor

    GeophysicalProcessor(100.0, 100.0, 10).bsdwormer(
        r.values.copy(), LayerShim(r.extent, "EPSG:32750", False), r.path, 3, 200, 200, False, 32750)
    table = np.genfromtxt(tmp_path / "g_worms.csv", delimiter=",", skip_header=1)
    np.testing.assert_array_equal(table, r.worms(3, 200, 200).table)


@needs_io
@needs_sklearn
def test_worm_lines(tmp_path):
    lines = make().worms(3, 200, 200).polylines()
    path = lines.save(tmp_path / "w.shp")
    v = raster_io.read_vector(path)
    assert len(v["wkb"]) == len(lines) and set(v["fields"]["Z"]) <= {200, 400, 600}
    assert raster_io.epsg(v["wkt"]) == 32750
    import shapely

    first = shapely.from_wkb(bytes(v["wkb"][0]))
    assert first.has_z and first.coords[0][2] == v["fields"]["Z"][0]


# ---------------------------------------------------------------- others
def test_mrvbf():
    y, x = np.mgrid[0:60, 0:60] * 30.0
    dem = sg.Grid(np.abs(x - 900) * 0.2, 30, name="dem")
    valley, ridge = dem.mrvbf()
    assert valley.name == "dem_MRVBF" and valley.shape == (60, 60)
    assert np.nanmean(valley.values[:, 28:32]) > np.nanmean(valley.values[:, 10:20])


def test_power_spectrum():
    n, dx, k = 256, 50.0, 2 * np.pi / 800
    g = sg.Grid(np.tile(np.sin(k * np.arange(n) * dx), (n, 1)), dx)
    wavenumber, power = g.power_spectrum()
    assert abs(wavenumber[1:][np.argmax(power[1:])] - k) < 2 * np.pi / (n * dx) * 1.1


@needs_io
@pytest.mark.skipif(not has("shapely"), reason="needs shapely")
def test_boundary(tmp_path):
    g = make()
    b = g.data_boundary()
    area = sum(p.area for p in b.geometries)
    assert area == np.isfinite(g.values).sum() * 1e4
    for ext in (".shp", ".gpkg"):
        v = raster_io.read_vector(b.save(tmp_path / ("b" + ext)))
        assert len(v["wkb"]) == len(b)
    with pytest.raises(ValueError):
        b.save(tmp_path / "b.xyz")


def test_normalise():
    a, b = make(crs=None), make(crs=None, seed=2)
    x, y = sg.normalise([a, b])
    assert abs(np.nanmean(x.values)) < 1e-3 and abs(np.nanstd(x.values) - np.nanstd(y.values)) < 1e-2
    with pytest.raises(ValueError):
        sg.normalise([a], order=3)


def test_radiometrics():
    k, u, th = (make(crs=None, seed=s) for s in (1, 2, 3))
    np.testing.assert_allclose(sg.efimov(k, u, th).values, k.values * u.values / th.values)
    for f in (sg.kd(k, th), sg.kd(k, th, "regression"), sg.ud(u, th),
              sg.ternary_gradient(k, u, th), sg.ternary_gradient(k, u, th, "clr")):
        assert f.shape == k.shape and np.isfinite(f.values).any()
    with pytest.raises(ValueError):
        sg.kd(k, th, "mean")


def test_log_ratio():
    s = sg.Stack(np.stack([make(crs=None, seed=i).values for i in range(3)]), 100, names=list("abc"))
    c = s.log_ratio()
    np.testing.assert_allclose(np.nansum(c.values, axis=0), 0, atol=1e-9)
    assert c.names == ["CLR a", "CLR b", "CLR c"] and len(s.log_ratio("ilr")) == 2
    assert c.info["replaced"] == 0


@needs_sklearn
def test_components(tmp_path):
    s = sg.Stack(np.stack([make(seed=i).values for i in range(3)]), gt=make().gt, crs=CRS)
    p = s.pca(2)
    ratio = [float(v) for v in p.info["EXPLAINED_VARIANCE_RATIO"].split(",")]
    assert len(p) == 2 and 0 < ratio[1] <= ratio[0] < 1
    assert len(s.ica()) == 3
    if HAS_IO:
        p.save(tmp_path / "p.tif")
        r = sg.read_stack(tmp_path / "p.tif")
        assert r.names[0].startswith("PC1 (")
        assert raster_io.read(tmp_path / "p.tif")["metadata"]["LOADINGS"] == p.info["LOADINGS"]


# ---------------------------------------------------------------- points
def points():
    rng = np.random.default_rng(3)
    x, y = rng.uniform(0, 5000, 400), rng.uniform(0, 4000, 400)
    return x, y, 10 + 0.01 * x - 0.02 * y


def test_grid_points():
    x, y, z = points()
    g = sg.grid_points(x, y, z, cell_size=100, levels=11)
    assert g.x[0] % 100 == 0 and g.y[-1] % 100 == 0
    inner = (slice(5, -5), slice(5, -5))
    expected = 10 + 0.01 * g.x[None, :] - 0.02 * g.y[:, None]
    np.testing.assert_allclose(g.values[inner], np.broadcast_to(expected, g.shape)[inner], atol=3)
    auto = sg.grid_points(x, y, z)
    assert auto.history[0]["parameters"]["settings_chosen_automatically"] == "True"
    with pytest.raises(ValueError):
        sg.grid_points(x[:2], y[:2], z[:2])
    with pytest.raises(ValueError):
        sg.grid_points(x, y[:5], z)


@needs_io
@pytest.mark.skipif(not has("shapely"), reason="needs shapely")
def test_read_points(tmp_path):
    import shapely

    x, y, z = points()
    path = str(tmp_path / "p.gpkg")
    raster_io.write_vector(path, [shapely.Point(a, b).wkb for a, b in zip(x, y)], "Point",
                           {"z": z, "n": np.arange(400), "t": ["a"] * 400},
                           raster_io.crs_to_wkt("EPSG:32750"))
    p = sg.read_points(path)
    np.testing.assert_allclose(p.x, x)
    np.testing.assert_allclose(p.fields["z"], z)
    assert raster_io.epsg(p.crs) == 32750 and len(p) == 400
    g = p.grid("z", cell_size=100, levels=11)
    np.testing.assert_allclose(g.values, sg.grid_points(x, y, z, 100, 11).values, rtol=1e-6)
    assert os.path.basename(g.save().path) == "p_z_bspline.tif"
    with pytest.raises(ValueError):
        p.grid("nope")


def test_aseg_gdf2(tmp_path):
    (tmp_path / "s.dfn").write_text(
        "DEFN   ST=RECD,RT=COMM;RT:A4;COMMENTS:A76\n"
        "DEFN 1 ST=RECD,RT=;EAST:F12.2\n"
        "DEFN 2 ST=RECD,RT=;NORTH:F12.2\n"
        "DEFN 3 ST=RECD,RT=;MAG:F10.2:NULL=-9999.00\n"
        "DEFN 4 ST=RECD,RT=;END DEFN\n")
    (tmp_path / "s.dat").write_text("".join(
        f"{500000 + i * 10:12.2f}{7000000 + i:12.2f}{50000 + i:10.2f}\n" for i in range(20)))
    with pytest.raises(ValueError, match="EAST"):
        sg.read_points(tmp_path / "s.dat")
    p = sg.read_points(tmp_path / "s.dat", x="EAST", y="NORTH")
    assert len(p) == 20 and p.x[1] == 500010 and p.fields["MAG"][3] == 50003


# ---------------------------------------------------------------- WTMM, plots, replay
@pytest.mark.skipif(not has("pywt"), reason="needs PyWavelets")
def test_wtmm():
    rng = np.random.default_rng(5)
    profile = np.cumsum(rng.normal(size=256))
    r = sg.wtmm(profile)
    assert len(r["scales"]) == 15 and len(r["maxima_lines"])
    h, d = sg.multifractal_spectrum(r)
    assert len(h) == len(d)
    with pytest.raises(ValueError):
        sg.wtmm(np.r_[profile, np.nan])


@needs_io
@needs_mpl
def test_petrophysics(tmp_path):
    import matplotlib

    matplotlib.use("Agg")
    d, s = make(), sg.Grid(np.abs(make(seed=4).values[::2, ::2]) / 1e4, 200, origin=ORIGIN, crs=CRS)
    a, b = sg.petrophysics.pair(d, s)
    assert a.shape == b.shape == s.shape
    fig = sg.petrophysics.heatmap(a.values, b.values, "t", "d", "s", bins=40)
    files = sg.petrophysics.save_figure(fig, str(tmp_path / "h"), ["png"])
    assert os.path.getsize(files[0]) > 1000
    sg.petrophysics.histogram(a.values, "t", "d")
    if has("shapely"):
        import shapely

        path = str(tmp_path / "poly.shp")
        raster_io.write_vector(
            path, [shapely.MultiPolygon([shapely.box(500000, 6995000, 505000, 7000000)]).wkb],
            "MultiPolygon", {"UNIT": ["basalt"]}, d.crs)
        cats, names = sg.petrophysics.categories(b, path, "UNIT")
        assert names == ["basalt"] and (cats == 0).sum() == 25 * 25


@needs_io
def test_replay(tmp_path):
    make().save(tmp_path / "a.tif")
    make(seed=7).save(tmp_path / "b.tif")
    out = (sg.read(tmp_path / "a.tif").reduction_to_pole(-60, 2).continuation("up", 200)
           .derivative("z", 1).save())
    r = sg.replay(out.path, tmp_path / "b.tif", date="2015-07-02")
    assert [os.path.basename(p) for p in r["outputs"]] == [
        "b_RTP.tif", "b_RTP_UC.tif", "b_RTP_UC_d1z.tif"]
    b = sg.read(tmp_path / "b.tif")
    inc, dec, _ = b.field_direction("2015-07-02")
    b.reduction_to_pole(inc, dec).save(tmp_path / "c1.tif")
    sg.read(tmp_path / "c1.tif").continuation("up", 200).save(tmp_path / "c2.tif")
    mine = sg.read(tmp_path / "c2.tif").derivative("z", 1)
    np.testing.assert_array_equal(sg.read(r["final"]).values, mine.values.astype(np.float32))
    # straight from memory
    again = sg.replay(b.reduction_to_pole(-60, 0), tmp_path / "a.tif", date="2015-07-02")
    assert os.path.basename(again["final"]) == "a_RTP.tif"
    with pytest.raises(ValueError):
        sg.replay(tmp_path / "a.tif", tmp_path / "b.tif")


# ---------------------------------------------------------------- the package
def test_imports_are_light():
    """Importing pysgtool loads no GUI, GDAL binding, plotting or optional package."""
    code = ("import sys, pysgtool; bad = [m for m in ('qgis', 'PyQt5', 'arcpy', 'matplotlib', "
            "'sklearn', 'networkx', 'rasterio', 'pyogrio') if m in sys.modules]; "
            "sys.exit(' '.join(bad) if bad else 0)")
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


def test_quiet(capsys):
    g = make(crs=None)
    g.total_horizontal_gradient()
    g.worms(2, 200, 200)
    from pysgtool.calcs.Cooper_Cowan_rtpvariable import rtpvariable

    rtpvariable(g.values, 90, 110, (-60,) * 4, (0,) * 4)
    assert capsys.readouterr().out == ""


def test_version():
    assert sg.__version__ and sg.__version__[0].isdigit()
