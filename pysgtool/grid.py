"""
Grid and Stack: a raster in memory, with every SGTool calculation as a method.

Each method returns a new object and leaves its own alone. The processing history
goes along in memory, so a grid saved at the end of a chain carries every step, as a
grid made step by step in QGIS does, and can be inspected or replayed there.
"""

import copy
import datetime
import os
import xml.etree.ElementTree as ET  # nosec B405 - building trees only

import numpy as np

from .calcs import raster_io, sgt_metadata
from .calcs.layer_shim import LayerShim
from .calcs.filter_table import GridContext, _plugin_version, auto_output_path
from .calcs.workflow import history_steps


def as_date(date):
    """A datetime.date from a date, 'yyyy-mm-dd' or a decimal year (2015.5)."""
    if isinstance(date, datetime.datetime):
        return date.date()
    if isinstance(date, datetime.date):
        return date
    if isinstance(date, str):
        return datetime.date.fromisoformat(date.strip())
    year = int(date)
    start = datetime.date(year, 1, 1)
    days = (datetime.date(year + 1, 1, 1) - start).days
    return start + datetime.timedelta(days=int(round((float(date) - year) * days)))


class _Raster:
    """What Grid (2-D) and Stack (bands, rows, cols) share."""

    _ndim = 2

    def __init__(self, values, dx=None, dy=None, *, origin=(0.0, 0.0), gt=None, crs=None,
                 geographic=None, name="grid"):
        values = np.array(values, dtype=float)
        if values.ndim != self._ndim:
            raise ValueError(f"{type(self).__name__} needs a {self._ndim}-D array, "
                             f"not {values.ndim}-D")
        values[~np.isfinite(values)] = np.nan
        if gt is None:
            if dx is None:
                raise ValueError("Give the cell size dx (and dy), or a geotransform gt")
            dy = dx if dy is None else dy
            gt = (origin[0], abs(dx), 0.0, origin[1], 0.0, -abs(dy))
        gt = tuple(float(v) for v in gt)
        if gt[1] <= 0 or gt[5] >= 0 or gt[2] or gt[4]:
            raise ValueError("The geotransform must be north-up (row 0 = north) and not rotated")
        self.values = values
        self.gt = gt
        self.crs = raster_io.crs_to_wkt(crs)
        self._geographic = geographic
        self.name = name
        self.path = None  # where it is saved, or would be by default
        self.metadata = {}  # text items of the file it was read from
        self._on_disk = False
        self._provenance = None  # XML element, once something has been calculated

    # ------------------------------------------------------------ facts
    @property
    def shape(self):
        """(rows, cols)."""
        return self.values.shape[-2:]

    @property
    def dx(self):
        return self.gt[1]

    @property
    def dy(self):
        return -self.gt[5]

    @property
    def extent(self):
        """(xmin, ymin, xmax, ymax) of the outer cell edges."""
        rows, cols = self.shape
        return (self.gt[0], self.gt[3] - rows * self.dy, self.gt[0] + cols * self.dx, self.gt[3])

    @property
    def epsg(self):
        """EPSG number of the coordinate system, or None."""
        return raster_io.epsg(self.crs)

    @property
    def geographic(self):
        """True if the coordinates are in degrees."""
        if self._geographic is None:
            self._geographic = raster_io.is_geographic(self.crs)
        return self._geographic

    @property
    def x(self):
        """Cell-centre x of each column."""
        return self.gt[0] + (np.arange(self.shape[1]) + 0.5) * self.dx

    @property
    def y(self):
        """Cell-centre y of each row (row 0 is the northernmost)."""
        return self.gt[3] - (np.arange(self.shape[0]) + 0.5) * self.dy

    def __array__(self, dtype=None, copy=None):
        return np.array(self.values, dtype=dtype)

    def __repr__(self):
        shape = " x ".join(str(n) for n in self.values.shape)
        code = self.epsg if self.crs and not self._crs_unreadable() else None
        crs = f"EPSG:{code}" if code else ("with CRS" if self.crs else "no CRS")
        return f"<{type(self).__name__} {self.name!r}: {shape}, cell {self.dx:g} x {self.dy:g}, {crs}>"

    def _crs_unreadable(self):
        try:
            raster_io.epsg(self.crs)
            return False
        except Exception:
            return True

    # ------------------------------------------------------------ history
    @property
    def provenance(self):
        """The processing history as an XML element (None for data as first read)."""
        if self._provenance is None and self._on_disk:
            return sgt_metadata.read_sgt_metadata(self.path)
        return self._provenance

    @property
    def history(self):
        """The steps that made this, oldest first, as dicts (operation, parameters ...)."""
        root = self.provenance
        return [] if root is None else history_steps(root)[1]

    def _source(self):
        """(path to record as the source, provenance held in memory for it)."""
        path = self.path or f"{self.name}.tif"
        if self._on_disk:
            return path, {}
        return path, {str(path): self._provenance if self._provenance is not None else False}

    def _derive(self, values, operation, parameters, suffix, cls=None):
        """A new raster on the same cells, with its provenance."""
        out = (cls or type(self))(values, gt=self.gt, crs=self.crs, geographic=self._geographic,
                                  name=self.name + suffix)
        source, known = self._source()
        out.path = auto_output_path(self.path, suffix) if self.path else None
        out._provenance = sgt_metadata._build(
            out.path or f"{out.name}.tif", source, operation, parameters, _plugin_version(),
            file_stats=False, source_provenance=known,
        )
        if "inclination" in self.metadata:  # a Noddy model's own field
            out.metadata = {k: self.metadata[k] for k in ("inclination", "declination", "intensity")
                            if k in self.metadata}
        return out

    # ------------------------------------------------------------ saving
    def _band_names(self):
        return None

    def _file_metadata(self):
        return None

    def save(self, path=None):
        """Save as a GeoTIFF with the processing history inside, and return self.

        path: by default beside the file the first step was read from, with SGTool's
        usual name (mag.tif -> mag_RTP.tif). An earlier file of that name is replaced."""
        if path is None:
            if self._on_disk:
                raise ValueError(f"This is already the file {self.path}: give another path")
            if self.path is None:
                raise ValueError("Give a path: the data did not come from a file")
            path = self.path
        path = os.path.abspath(str(path))
        if os.path.splitext(path)[1].lower() not in (".tif", ".tiff"):
            raise ValueError("Grids are saved as GeoTIFF: give a path that ends in .tif")
        root = copy.deepcopy(self.provenance)
        for side in (".aux.xml", sgt_metadata.SIDECAR_SUFFIX):
            if os.path.exists(path + side):
                os.remove(path + side)
        stack = self.values if self.values.ndim == 3 else self.values[None]
        raster_io.write(path, stack, self.gt, self.crs, self._band_names(), self._file_metadata())
        if root is not None:
            root.find("output/path").text = path
            root.find("created").text = datetime.datetime.now().astimezone().isoformat(
                timespec="seconds")
            if not sgt_metadata._embed_in_geotiff(path, root):
                ET.ElementTree(root).write(sgt_metadata.sidecar_path(path), encoding="utf-8",
                                           xml_declaration=True)
        self._provenance = root
        self.path, self._on_disk = path, True
        self.name = os.path.splitext(os.path.basename(path))[0]
        return self

    # ------------------------------------------------------------ shared
    def resample_to(self, other, how=None):
        """This raster on the cells of another (reprojected if need be). how: 'nearest',
        'bilinear' or 'average'; by default average onto larger cells, else bilinear."""
        if bool(self.crs) != bool(other.crs):
            raise ValueError("One grid has a coordinate system and the other has none")
        stack = self.values if self.values.ndim == 3 else self.values[None]
        out = raster_io.resample(stack, self.gt, self.crs, other.gt, other.shape, other.crs, how)
        new = copy.copy(self)
        new.values = out if self.values.ndim == 3 else out[0]
        new.gt, new.crs, new._geographic = other.gt, other.crs, other._geographic
        if not np.isfinite(new.values).any():
            raise ValueError(f"{self.name} does not overlap {other.name}")
        new._provenance = self.provenance
        new._on_disk = False
        new.path = auto_output_path(self.path, "_resampled") if self.path else None
        return new

    def _require_projected(self, what):
        if self.geographic:
            raise ValueError(f"{what} needs a grid in a projected (metre-based) coordinate system")

    def _context(self, buffer=0):
        if not np.isfinite(self.values).any():
            raise ValueError("The grid has no data")
        authid = ""
        if self.crs:
            code = self.epsg if not self._crs_unreadable() else None
            authid = f"EPSG:{code}" if code else ""
        return GridContext(self.values.copy(), self.gt,
                           LayerShim(self.extent, authid, self.geographic), buffer)


class Grid(_Raster):
    """A single-band grid: a north-up array (row 0 = north, NaN = no data) with its
    position. The filters, with the names and settings of the QGIS Processing
    algorithms, are methods: see pysgtool.filters() or help(pysgtool.Grid).

    Grid(array, dx, dy=None, origin=(xmin, ymax), crs='EPSG:28350', name='mag'), or
    pysgtool.read(path).
    """

    # ------------------------------------------------------------ magnetic field
    def field_direction(self, date):
        """(inclination, declination, intensity in nT) of the IGRF at the grid's centre
        at sea level on a date, or of the model a Noddy grid was made from."""
        own = self._own_field()
        if own is not None:
            return own
        from .calcs.igrf.igrf_field import _lonlat_arrays, calc_igrf

        xmin, ymin, xmax, ymax = self.extent
        x, y = np.array([(xmin + xmax) / 2]), np.array([(ymin + ymax) / 2])
        lon, lat = _lonlat_arrays(self._crs_for_field(), x, y)
        return calc_igrf(_decimal(date), 0.0, float(lat[0]), float(lon[0]))

    def _own_field(self):
        try:
            return tuple(float(self.metadata[k]) for k in ("inclination", "declination", "intensity"))
        except (KeyError, ValueError):
            return None

    def _crs_for_field(self):
        if not self.crs:
            raise ValueError("The IGRF needs the grid's coordinate system: set crs")
        code = self.epsg
        return f"EPSG:{code}" if code else self.crs

    def _field_arrays(self, date, callback=None):
        from .calcs.igrf.igrf_field import grid_igrf

        return grid_igrf(self.gt, self.shape, self._crs_for_field(), _decimal(date),
                         callback=callback)

    def igrf(self, date, callback=None):
        """IGRF inclination, declination (degrees) and total intensity (nT) at every
        cell, at sea level, as three Grids."""
        arrays = self._field_arrays(date, callback)
        record = self._igrf_record(date)
        return tuple(
            self._derive(a, f"IGRF {what}", record, suffix, Grid)
            for a, what, suffix in zip(arrays, ("inclination", "declination", "total field"),
                                       ("_IGRF_inc", "_IGRF_dec", "_IGRF"))
        )

    @staticmethod
    def _igrf_record(date):
        from .calcs.igrf.igrf_field import IGRF_GENERATION

        return {"date": as_date(date).isoformat(), "altitude_km": 0.0,
                "model": "IGRF-" + IGRF_GENERATION}

    def remove_igrf(self, date, callback=None):
        """Total magnetic intensity minus the IGRF total field at every cell."""
        if self._own_field() is not None:
            raise ValueError("This grid carries its own field (a Noddy model): no IGRF to remove")
        field = self._field_arrays(date, callback)[2]
        return self._derive(self.values - field, "Remove IGRF (TMI - IGRF)",
                            self._igrf_record(date), "_RMI")

    def differential_rtp(self, date):
        """Reduction to the pole with the IGRF inclination and declination of every cell
        (Cooper & Cowan)."""
        from .calcs.Cooper_Cowan_rtpvariable import rtpvariable

        own = self._own_field()
        if own is not None:
            inc, dec = (np.full(self.shape, v) for v in own[:2])
        else:
            inc, dec, _ = self._field_arrays(date)
        rows, cols = self.shape
        _, rtp, _, _ = rtpvariable(np.flipud(self.values), rows, cols, None, None,
                                   incv=np.flipud(inc), decv=np.flipud(dec))
        return self._derive(np.flipud(rtp), "Differential reduction to pole",
                            {"date": as_date(date).isoformat()}, "_DRTP")

    # ------------------------------------------------------------ continuation
    def variable_continuation(self, observed, output, terms=6, smoothing=0.0, cancel_check=None):
        """Continue from an observed to an output surface, either of them uneven
        (Cordell & Grauch 1985). Each surface is a height (a number, m), a Grid of
        heights, or (DEM Grid, clearance in m).

        terms: of the Taylor series. smoothing: wavelength (m) below which downward
        continuation is damped (0 = automatic). Warnings about the surfaces are issued
        with the warnings module. Returns a Grid with the run's facts in .info."""
        import warnings

        from .calcs import TaylorContinuation as T

        def surface(s, what):
            if isinstance(s, Grid):
                mode, grid, value = "grid", s, 0.0
            elif isinstance(s, (tuple, list)):
                if len(s) != 2 or not isinstance(s[0], Grid):
                    raise ValueError(f"The {what} surface must be a height, a Grid or (DEM Grid, clearance)")
                mode, (grid, value) = "clearance", s
            else:
                mode, grid, value = "height", None, float(s)
            array = None if grid is None else grid.resample_to(self).values
            z, dem = T.surface_heights(mode, array, value)
            return (mode, z, dem, value), T.surface_text(mode, grid.name if grid else "", value)

        obs, obs_text = surface(observed, "observed")
        out, out_text = surface(output, "output")
        for text in T.check_surfaces(obs, out, np.isfinite(self.values)):
            warnings.warn(text, stacklevel=2)
        dx, dy = self.dx, self.dy
        if self.geographic:
            mx, my = T.metres_per_degree((self.extent[1] + self.extent[3]) / 2)
            dx, dy = dx * mx, dy * my
        values, info = T.variable_continuation(
            self.values, obs[1], out[1], dx, dy, order=terms, smoothing=smoothing,
            cancel_check=cancel_check)
        if info["diverged"]:
            warnings.warn(
                "the flattening stopped improving at a misfit of "
                f"{100 * info['relative_misfit']:.1f}% of the data; check the observed heights "
                "match the data", stacklevel=2)
        grid = self._derive(values, "Variable continuation (Taylor series)",
                            {"observed": obs_text, "output": out_text,
                             "smoothing": smoothing or "auto", "terms": terms}, "_VC")
        grid.info = info
        return grid

    # ------------------------------------------------------------ sources and edges
    def euler_deconvolution(self, structural_index=1, window=10, keep=0.1, callback=None):
        """Euler deconvolution for one structural index (0 contact, 1 line of poles,
        2 pole, 3 dipole). window: size in cells. keep: fraction of the solutions kept,
        those with the largest standard deviation of df/dz. Returns EulerSolutions."""
        from .results import EulerSolutions
        from .calcs.euler.euler_python import euler_deconv_grid

        self._require_projected("Euler deconvolution")
        if structural_index < 0 or not 0 < keep <= 1 or not 3 <= window < min(self.shape):
            raise ValueError("structural_index must be >= 0, keep in (0, 1] and window from 3 "
                             "to less than the grid's size")
        si = structural_index or 0.001  # 0.001 stands in for 0, as in the plugin
        table = euler_deconv_grid(self.values, self.extent, si, int(window), keep, callback=callback)
        return EulerSolutions(self, table, si, int(window), keep)

    def euler_window_statistics(self, folder, window=10, keep=0.1, structural_indices=(0, 1, 2, 3)):
        """Statistics of the Euler solutions in windows over the grid, written to folder as
        the plugin's <name>_window_stats_SI_<n>.txt and <name>_window_summary.txt.
        Returns {structural index: table}."""
        from .calcs.euler.estimates_statistics import window_stats

        solutions = [self.euler_deconvolution(si, window, keep) for si in structural_indices]
        xmin, ymin, xmax, ymax = self.extent
        os.makedirs(folder, exist_ok=True)
        tables = window_stats(
            [s.table[:, [1, 0, 2, 3, 4]] for s in solutions], [ymin, ymax, xmin, xmax],
            [s.structural_index for s in solutions], self.name, str(folder), self.shape,
            int(window), detailed_stats=True)
        return {si: tables[s.structural_index] for si, s in zip(structural_indices, solutions)}

    def worms(self, levels=10, base=1000.0, increment=1000.0, callback=None):
        """BSDWorms: multilevel wavelet edges of a gravity grid, or an RTP / RTE magnetic
        grid after vertical_integration(). base: upward continuation height of the first
        level, increment: between levels (map units). Returns Worms."""
        from .results import Worms
        from .calcs.GeophysicalProcessor import GeophysicalProcessor

        self._require_projected("Worms")
        processor = GeophysicalProcessor(self.dx, self.dy, 0)
        table = processor.worm_points(self.values, self.gt, int(levels), base, increment, callback)
        return Worms(self, table, {"levels": int(levels), "bottom_level": base,
                                   "level_increment": increment})

    # ------------------------------------------------------------ the rest
    def mrvbf(self, t_slope=16.0, t_pctl_v=0.40, t_pctl_r=0.35, p_slope=4.0, p_pctl=3.0,
              max_res=100.0, callback=None):
        """(MRVBF, MRRTF) of a DEM: the multiresolution index of valley bottom flatness
        and of ridge top flatness (Gallant & Dowling 2003), as two Grids."""
        from .calcs.mrvbf import mrvbf

        self._require_projected("MRVBF")
        record = dict(t_slope=t_slope, t_pctl_v=t_pctl_v, t_pctl_r=t_pctl_r, p_slope=p_slope,
                      p_pctl=p_pctl, max_res=max_res)
        dem = np.flipud(self.values)
        valley, ridge = mrvbf(dem, self.dx, xmin=self.extent[0], ymin=self.extent[1],
                              nodata_mask=~np.isfinite(dem), callback=callback, **record)
        return (
            self._derive(np.flipud(valley), "Multiresolution valley bottom flatness", record, "_MRVBF"),
            self._derive(np.flipud(ridge), "Multiresolution ridge top flatness", record, "_MRRTF"),
        )

    def power_spectrum(self):
        """(wavenumber in radians per map unit, power): the radially averaged power
        spectrum, with no data as 0 as in the plugin."""
        dock = self._spectrum_dock()
        _display, wavenumber, power, _limits = dock.compute()
        return wavenumber, power

    def plot_power_spectrum(self):
        """Show the grid and its radially averaged power spectrum (needs matplotlib)."""
        self._spectrum_dock().plot_grid_and_power_spectrum()

    def _spectrum_dock(self):
        from .calcs.PSplot import PowerSpectrumDock

        xmin, ymin, xmax, ymax = self.extent
        rows, cols = self.shape
        return PowerSpectrumDock(self.values.copy(), self.name, self.dx, self.dy,
                                 np.linspace(xmin, xmax, cols), np.linspace(ymin, ymax, rows))

    def data_boundary(self):
        """Outline of the cells with data (holes included) as polygons (needs shapely)."""
        from shapely.geometry import MultiPolygon, Polygon

        from .results import Features
        from .calcs.SG_Util import SG_Util

        mask = np.where(np.isnan(self.values), 0, 1)
        rings, exteriors, holes = SG_Util(self.values).boundary_rings(mask, self.gt)
        inside = {}
        for hole, outer in holes:
            inside.setdefault(outer, []).append(rings[hole])
        polygons = [MultiPolygon([Polygon(rings[e], inside.get(e, []))]) for e in exteriors]
        return Features(polygons, "MultiPolygon",
                        {"id": list(range(1, len(polygons) + 1)),
                         "type": ["Multipolygon boundary"] * len(polygons)},
                        self, "Data boundary polygon", {})


def _decimal(date):
    from .calcs.igrf.igrf_field import decimal_year

    if isinstance(date, (int, float)):
        return float(date)
    when = as_date(date)
    return decimal_year(when.year, when.month, when.day)


class Stack(_Raster):
    """Several bands on the same cells: (bands, rows, cols), with their names.
    Stack(array, dx, ...) , pysgtool.read_stack(path) or Stack.from_grids([k, th, u])."""

    _ndim = 3

    def __init__(self, values, dx=None, dy=None, *, names=None, **kwargs):
        super().__init__(values, dx, dy, **kwargs)
        n = self.values.shape[0]
        self.names = list(names) if names else [f"Band {i + 1}" for i in range(n)]
        if len(self.names) != n:
            raise ValueError(f"{len(self.names)} names for {n} bands")
        self.info = {}

    @classmethod
    def from_grids(cls, grids, name=None):
        """Grids as the bands of a Stack, on the cells of the first."""
        grids = list(grids)
        first = grids[0]
        values = np.stack([first.values] + [g.resample_to(first).values for g in grids[1:]])
        out = cls(values, gt=first.gt, crs=first.crs, geographic=first._geographic,
                  names=[g.name for g in grids], name=name or first.name + "_stack")
        if first.path:
            out.path = os.path.join(os.path.dirname(first.path), out.name + ".tif")
        out._provenance = sgt_metadata._build(
            out.path or f"{out.name}.tif", [g._source()[0] for g in grids], "Stack of grids",
            {"bands": ", ".join(out.names)}, _plugin_version(), file_stats=False,
            source_provenance={k: v for g in grids for k, v in g._source()[1].items()})
        return out

    def _band_names(self):
        return self.names

    def _file_metadata(self):
        return {k: v for k, v in self.info.items() if isinstance(v, str)}

    def __len__(self):
        return self.values.shape[0]

    def band(self, which):
        """One band as a Grid: its number (from 1) or its name."""
        i = self.names.index(which) if isinstance(which, str) else int(which) - 1
        if not 0 <= i < len(self):
            raise ValueError(f"There is no band {which}")
        return self._derive(self.values[i], "Band of a multiband grid",
                            {"band": i + 1, "name": self.names[i]}, f"_b{i + 1}", Grid)

    def _components(self, kind, n_components):
        from .calcs import PCAICA

        if len(self) < 2:
            raise ValueError("Component analysis needs at least two bands")
        values, model = getattr(PCAICA, kind + "_stack")(self.values, n_components)
        operation = {"pca": "Principal component analysis", "ica": "Independent component analysis"}
        out = self._derive(values, operation[kind], {"n_components": values.shape[0]},
                           f"_{kind.upper()}")
        return out, model

    def pca(self, n_components=None):
        """Principal components, one band each (needs scikit-learn). .info holds the
        explained variance (EXPLAINED_VARIANCE_RATIO) and the LOADINGS."""
        out, model = self._components("pca", n_components)
        ratio = model.explained_variance_ratio_
        out.names = [f"PC{i + 1} ({100 * r:.2f}%)" for i, r in enumerate(ratio)]
        out.info = {
            "EXPLAINED_VARIANCE": ",".join(f"{v:.6f}" for v in model.explained_variance_),
            "EXPLAINED_VARIANCE_RATIO": ",".join(f"{v:.6f}" for v in ratio),
            "LOADINGS": ";".join(",".join(f"{v:.6f}" for v in c) for c in model.components_),
        }
        return out

    def ica(self, n_components=None):
        """Independent components, one band each (needs scikit-learn)."""
        out, _model = self._components("ica", n_components)
        out.names = [f"IC{i + 1}" for i in range(len(out))]
        return out

    def log_ratio(self, kind="CLR"):
        """Centred ('CLR') or isometric ('ILR') log-ratio of every band. Values <= 0 are
        replaced by half the smallest positive value of their band: .info['replaced'] is
        how many."""
        from .calcs import Radiometrics as R

        kind = kind.upper()
        if kind not in R.LOG_RATIO_OPERATIONS:
            raise ValueError("kind is 'CLR' or 'ILR'")
        if len(self) < 2:
            raise ValueError("The log-ratio transform needs at least two bands")
        values, replaced = R.log_ratio(kind, self.values)
        operation, suffix = R.LOG_RATIO_OPERATIONS[kind]
        out = self._derive(values, operation, {"type": kind, "bands": ", ".join(self.names)}, suffix)
        out.names = R.log_ratio_band_names(kind, self.names)
        out.info = {"replaced": replaced}
        return out
