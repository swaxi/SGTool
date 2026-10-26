"""
The SGTool filters and the grid helpers they need, with no QGIS dependency.

FILTER_SPECS describes every filter (its settings and how to run it), and
GridContext carries a grid to a filter's compute function. The QGIS Processing
algorithms (sgt_processing.py), the "apply the same steps to another grid" replay
(sgt_workflow.py) and the ArcGIS Pro GUI all use this one table, so a filter gives
the same result wherever it is run.
"""

import os

import numpy as np
from osgeo import gdal

try:  # as part of the QGIS plugin package
    from .calcs.GeophysicalProcessor import GeophysicalProcessor
except ImportError:  # from the ArcGIS Pro GUI, where calcs is a top-level package
    from calcs.GeophysicalProcessor import GeophysicalProcessor


def _calc_class(module, name):
    """A calculation class, imported when first needed: some (SG_Util needs shapely)
    are not installed in every environment, and most filters do not use them."""
    import importlib

    try:
        return getattr(importlib.import_module(f".calcs.{module}", __package__), name)
    except (ImportError, TypeError, ValueError):
        return getattr(importlib.import_module(f"calcs.{module}"), name)


class FilterError(Exception):
    """A filter cannot be run on this grid (the message says why)."""


# ---------------------------------------------------------------- grid I/O
def _plugin_version():
    try:
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "metadata.txt")
        with open(path, encoding="utf-8") as f:
            for line in f:
                if line.startswith("version="):
                    return line.split("=", 1)[1].strip()
    except OSError:
        pass
    return None


def read_grid(path):
    """Band 1 of a raster as a north-up float array (NaN for no data), with its
    geotransform and projection, as the plugin reads grids."""
    ds = gdal.Open(path)
    if ds is None:
        raise FilterError(f"Could not open the raster {path}")
    band = ds.GetRasterBand(1)
    nodata = band.GetNoDataValue()
    arr = band.ReadAsArray().astype(float)
    if nodata is not None:
        arr[arr == nodata] = np.nan
    gt = ds.GetGeoTransform()
    projection = ds.GetProjection()
    ds = None
    if gt[5] > 0:  # south-up: flip to north-up
        arr = np.flipud(arr)
        gt = (gt[0], gt[1], 0.0, gt[3] + arr.shape[0] * gt[5], 0.0, -gt[5])
    return arr, gt, projection


def write_grid(path, arr, gt, projection):
    rows, cols = arr.shape
    ds = gdal.GetDriverByName("GTiff").Create(path, cols, rows, 1, gdal.GDT_Float32)
    if ds is None:
        raise FilterError(f"Could not create {path}")
    ds.SetGeoTransform(gt)
    if projection:
        ds.SetProjection(projection)
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(np.nan)
    band.WriteArray(arr.astype(np.float32))
    band.FlushCache()
    band = None
    ds = None


class GridContext:
    """Everything a calculation needs about the input grid, with the
    calculation helpers created on first use."""

    def __init__(self, arr, gt, layer, buffer):
        self.arr = arr
        self.gt = gt
        self.dx = abs(gt[1])
        self.dy = abs(gt[5])
        rows, cols = arr.shape
        self.buffer = int(buffer) if buffer and buffer > 0 else min(rows, cols, 5000)
        ext = layer.extent()
        self.geographic = layer.crs().isGeographic()
        self.centroid = (
            (ext.xMinimum() + ext.xMaximum()) / 2.0,
            (ext.yMinimum() + ext.yMaximum()) / 2.0,
        )
        self.extent = (ext.xMinimum(), ext.yMinimum(), ext.xMaximum(), ext.yMaximum())
        self._cache = {}

    def _get(self, key, make):
        if key not in self._cache:
            self._cache[key] = make()
        return self._cache[key]

    @property
    def processor(self):
        return self._get("p", lambda: GeophysicalProcessor(self.dx, self.dy, self.buffer))

    @property
    def convolution(self):
        return self._get("c", lambda: _calc_class("ConvolutionFilter", "ConvolutionFilter")(self.arr))

    @property
    def sg_util(self):
        return self._get("u", lambda: _calc_class("SG_Util", "SG_Util")(self.arr))

    @property
    def spatial(self):
        return self._get("s", lambda: _calc_class("SpatialStats", "SpatialStats")(self.arr))

    def require_projected(self, what):
        if self.geographic:
            raise FilterError(
                f"{what} needs a grid in a projected (metre-based) coordinate system"
            )

    def degrees_to_cells_scale(self):
        """Average metres per degree at the grid centre (geographic grids)."""
        _long, lat = self.centroid
        dx, dy = self.sg_util.arc_degree_to_meters(lat)
        return np.sqrt(dx**2.0 + dy**2.0) / 2


def _fmt(value):
    """A number as the plugin dialog would show it in a file name: 500.0 -> '500'."""
    return "%.10g" % float(value)


# File-name suffix for each filter, the same as the dialog gives its outputs
# (grid.tif -> grid_d1z.tif), worked out from the algorithm's settings.
SUFFIXES = {
    "directional_line_noise": lambda p: (
        "_DirC_noise" if p["OUTPUT_TYPE"].startswith("Noise") else "_DirC"
    ),
    "reduction_to_pole": lambda p: "_RTP",
    "reduction_to_equator": lambda p: "_RTE",
    "reduction_to_pole_low_latitude": lambda p: "_RTP_LL",
    "continuation": lambda p: ("_UC_" if p["DIRECTION"] == "up" else "_DC_") + _fmt(p["HEIGHT"]),
    "vertical_integration": lambda p: "_VI",
    "remove_regional": lambda p: "_RR_" + p["ORDER"][0] + "o",
    "band_pass": lambda p: (
        "_BP_" + _fmt(p["LOW_CUT"] if p["LOW_CUT"] > 0 else 1e-10) + "_" + _fmt(p["HIGH_CUT"])
    ),
    "high_low_pass": lambda p: ("_LP_" if p["TYPE"] == "Low" else "_HP_") + _fmt(p["CUTOFF"]),
    "automatic_gain_control": lambda p: "_AGC",
    "derivative": lambda p: "_d" + _fmt(p["POWER"]) + p["DIRECTION"],
    "tilt_angle": lambda p: "_TA",
    "analytic_signal": lambda p: "_AS",
    "total_horizontal_gradient": lambda p: "_THG",
    "mean_filter": lambda p: "_Mn",
    "median_filter": lambda p: "_Md",
    "gaussian_filter": lambda p: "_Gs",
    "directional_filter": lambda p: "_Dr",
    "sun_shading": lambda p: "_Sh",
    "windowed_statistic": lambda p: {
        "min": "_SS_Min", "max": "_SS_Max", "std": "_SS_StdDev",
        "variance": "_SS_Var", "skewness": "_SS_Skew", "kurtosis": "_SS_Kurt",
    }[p["STATISTIC"]],
    "local_anisotropy": lambda p: (
        "_SS_AnisoMag" if p["OUTPUT_TYPE"] == "Magnitude" else "_SS_AnisoOrient"
    ),
    "anisotropy_chain_length": lambda p: "_SS_ChainLen",
    "anisotropy_streamline_length": lambda p: "_SS_StreamLen",
    "threshold_to_nan": lambda p: "_Clean",
    "pca": lambda p: "_PCA",
    "ica": lambda p: "_ICA",
    "apparent_susceptibility": lambda p: "_AppSus",
    "apparent_density": lambda p: "_AppDens",
}


def auto_output_path(source_path, suffix, extension=".tif"):
    """Where the dialog would save a result: next to the source, named after it
    plus the processing step."""
    base, _ext = os.path.splitext(source_path)
    return base + suffix + extension


def _check_length(g, length):
    """Lengths in degrees for geographic grids (as the plugin dialog checks)."""
    if g.geographic and length > 100:
        raise FilterError(
            "Since this is a geographic projection, you need to specify lengths in degrees"
        )


# ---------------------------------------------------------------- the filters
# Each compute function takes (GridContext, values) and returns the new grid.
def _c_line_noise(g, p):
    return g.processor.line_noise_removal(
        g.arr,
        p["AZIMUTH"],
        p["LINE_SPACING_MIN"],
        p["LINE_SPACING_MAX"] or None,
        scale=p["SCALE"],
        buffer_size=g.buffer,
        direction_width=p["WEDGE"],
        return_noise=p["OUTPUT_TYPE"].startswith("Noise"),
    )


def _c_rtp(g, p):
    return g.processor.reduction_to_pole(
        g.arr, inclination=p["INCLINATION"], declination=p["DECLINATION"],
        buffer_size=g.buffer,
    )


def _c_rte(g, p):
    return g.processor.reduction_to_equator(
        g.arr, inclination=p["INCLINATION"], declination=p["DECLINATION"],
        buffer_size=g.buffer,
    )


def _c_rtp_ll(g, p):
    grid, _m = g.processor.reduction_to_pole_low_latitude(
        g.arr, inclination=p["INCLINATION"], declination=p["DECLINATION"],
        iterations=p["ITERATIONS"], buffer_size=g.buffer,
    )
    return grid


def _c_regional(g, p):
    data, nodata_value = g.processor.fix_extreme_values(g.arr.copy())
    data = data.astype(np.float32)
    mask = (data == nodata_value) | np.isnan(data)
    if p["ORDER"].startswith("1"):
        return g.processor.remove_gradient(data, mask)
    return g.processor.remove_2o_gradient(data, mask)


def _c_derivative(g, p):
    return g.processor.compute_derivative(
        g.arr, direction=p["DIRECTION"], order=p["POWER"], buffer_size=g.buffer
    )


def _c_tilt(g, p):
    # in degrees (-90 to 90), as the plugin dialog and the ArcGIS Pro version save it
    return np.degrees(g.processor.tilt_angle(g.arr, buffer_size=g.buffer))


def _c_analytic(g, p):
    return g.processor.analytic_signal(g.arr, buffer_size=g.buffer)


def _c_thg(g, p):
    return g.processor.total_hz_grad(g.arr, buffer_size=g.buffer)


def _c_vint(g, p):
    g.require_projected("Vertical integration")
    return g.processor.vertical_integration(
        g.arr, max_wavenumber=None, min_wavenumber=1e-4,
        buffer_size=g.buffer, buffer_method="mirror",
    )


def _petro_cut(g, cut):
    """High-cut wavelength for the apparent petrophysics filters: 0 means
    4 x the cell size, a negative value no high-cut."""
    if not cut:
        return 4 * max(g.dx, g.dy)
    return None if cut < 0 else cut


def _c_app_sus(g, p):
    g.require_projected("Apparent susceptibility")
    if p["FIELD"] <= 0:
        raise FilterError("Set the field strength (nT) for apparent susceptibility")
    return g.processor.apparent_susceptibility(
        g.arr, p["FIELD"], depth=p["DEPTH"], thickness=p["THICKNESS"] or None,
        background=p["BACKGROUND"], cutoff_wavelength=_petro_cut(g, p["HIGH_CUT"]),
        buffer_size=g.buffer,
    )


def _c_app_dens(g, p):
    g.require_projected("Apparent density")
    if p["THICKNESS"] <= 0:
        raise FilterError("The apparent density layer thickness must be greater than 0")
    to_mgal = 0.1 if p["UNITS"].startswith("um") else 1.0  # 1 mGal = 10 um/s^2
    return g.processor.apparent_density(
        g.arr * to_mgal, depth=p["DEPTH"], thickness=p["THICKNESS"],
        background=p["BACKGROUND"], cutoff_wavelength=_petro_cut(g, p["HIGH_CUT"]),
        buffer_size=g.buffer,
    )


def _c_continuation(g, p):
    height = p["HEIGHT"]
    if g.geographic:
        height = height / g.degrees_to_cells_scale()
    if p["DIRECTION"] == "up":
        return g.processor.upward_continuation(g.arr, height=height, buffer_size=g.buffer)
    return g.processor.downward_continuation(g.arr, height=height, buffer_size=g.buffer)


def _c_bandpass(g, p):
    low = p["LOW_CUT"] if p["LOW_CUT"] > 0 else 1e-10
    _check_length(g, low)
    _check_length(g, p["HIGH_CUT"])
    return g.processor.band_pass_filter(
        g.arr, low_cut=low, high_cut=p["HIGH_CUT"],
        high_transition_width=p["WIDTH"], low_transition_width=p["WIDTH"],
        buffer_size=g.buffer,
    )


def _c_highlow(g, p):
    _check_length(g, p["CUTOFF"])
    if p["TYPE"] == "Low":
        return g.processor.low_pass_filter(
            g.arr, cutoff_wavelength=p["CUTOFF"], transition_width=p["WIDTH"],
            buffer_size=g.buffer,
        )
    return g.processor.high_pass_filter(
        g.arr, cutoff_wavelength=p["CUTOFF"], transition_width=p["WIDTH"],
        buffer_size=g.buffer,
    )


def _c_agc(g, p):
    return g.processor.automatic_gain_control(g.arr, window_size=int(p["WINDOW"]))


def _c_mean(g, p):
    return g.convolution.mean_filter(int(p["SIZE"]))


def _c_median(g, p):
    return g.convolution.median_filter(int(p["SIZE"]))


def _c_gaussian(g, p):
    return g.convolution.gaussian_filter(p["SIGMA"])


def _c_directional(g, p):
    return g.convolution.directional_filter(p["DIRECTION"], n=3)


def _c_sunshade(g, p):
    if p["RELIEF"]:
        hz = g.degrees_to_cells_scale() if g.geographic else 1.0
        return g.convolution.sun_shading_filter_grass(
            g.arr, altitude=p["ZENITH"], azimuth=p["AZIMUTH"],
            resolution_ns=g.dy * hz, resolution_ew=g.dx * hz, scale=1.0, zscale=1.0,
        )
    return g.convolution.sun_shading_filter(
        g.arr, sun_alt=p["ZENITH"], sun_az=180 - p["AZIMUTH"]
    )


def _c_nan(g, p):
    return g.sg_util.Threshold2Nan(
        g.arr, condition=p["CONDITION"],
        above_threshold_value=p["ABOVE"], below_threshold_value=p["BELOW"],
    )


def _c_window_stat(g, p):
    return g.spatial.calculate_windowed_stats(
        window_size=int(p["WINDOW"]), stat_type=p["STATISTIC"]
    )


def _c_anisotropy(g, p):
    magnitude, orientation = g.spatial.calculate_local_anisotropy(
        g.spatial.grid, window_size=int(p["WINDOW"]), window_type=p["WINDOW_TYPE"]
    )
    return magnitude if p["OUTPUT_TYPE"] == "Magnitude" else orientation


def _anisotropy_maps(g, p):
    return g.spatial.calculate_local_anisotropy(
        g.spatial.grid, window_size=int(p["WINDOW"]), window_type=p["WINDOW_TYPE"]
    )


def _c_chain_length(g, p):
    magnitude, orientation = _anisotropy_maps(g, p)
    return g.spatial.calculate_anisotropy_chain_length(
        magnitude, orientation, aniso_threshold=p["ANISOTROPY_THRESHOLD"],
        angle_tolerance=p["ANGLE_TOLERANCE"], search_radius=int(p["SEARCH_RADIUS"]),
    )


def _c_streamline_length(g, p):
    magnitude, orientation = _anisotropy_maps(g, p)
    return g.spatial.calculate_anisotropy_streamline_length(
        magnitude, orientation, aniso_threshold=p["ANISOTROPY_THRESHOLD"],
        angle_tolerance=p["ANGLE_TOLERANCE"], max_steps=int(p["MAX_STEPS"]),
    )


def _dbl(name, label, default, minimum=None, maximum=None, optional=False):
    return dict(name=name, kind="double", label=label, default=default,
                minimum=minimum, maximum=maximum, optional=optional)


def _int(name, label, default, minimum=None, maximum=None):
    return dict(name=name, kind="int", label=label, default=default,
                minimum=minimum, maximum=maximum, optional=False)


def _opt(name, label, options, default=0):
    return dict(name=name, kind="enum", label=label, options=options, default=default)


def _bool(name, label, default):
    return dict(name=name, kind="bool", label=label, default=default)


GROUP_GRAV = ("grav_mag", "Grav/Mag filters")
GROUP_FREQ = ("frequency", "Frequency filters")
GROUP_GRAD = ("gradient", "Gradient filters")
GROUP_CONV = ("convolution", "Convolution filters")
GROUP_STATS = ("statistics", "Spatial statistics")
GROUP_UTIL = ("utilities", "Utilities")
GROUP_GRID = ("gridding", "Gridding")
GROUP_EULER = ("euler", "Euler deconvolution")
GROUP_MULTI = ("multivariate", "Multivariate analysis")
GROUP_PETRO = ("petrophysics", "Petrophysics")

_DIRECTIONS = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"]

# fft=True adds an "FFT buffer" parameter and runs the filter on a padded grid
FILTER_SPECS = [
    dict(
        id="directional_line_noise",
        name="Remove line noise (directional Cosine/Butterworth)",
        group=GROUP_FREQ, fft=True, compute=_c_line_noise,
        operation="Directional Cosine/Butterworth line-noise removal",
        help=(
            "Removes line-parallel acquisition noise. A zero-centred noise "
            "estimate is made from wavelengths between 2 x the smallest and 10 x "
            "the largest line spacing within a wedge (45 degrees by default) about the given "
            "azimuth, multiplied by the scale, and subtracted from the grid."
        ),
        params=[
            _dbl("AZIMUTH", "Azimuth of the noise wedge (degrees clockwise from north)", 0.0),
            _dbl("LINE_SPACING_MIN", "Smallest line spacing (map units)", 1000.0, minimum=0.0),
            _dbl("LINE_SPACING_MAX", "Largest line spacing (map units, 0 = same as smallest)",
                 0.0, minimum=0.0),
            _dbl("SCALE", "Scale applied to the noise estimate before subtracting", 1.0),
            _dbl("WEDGE", "Wedge half-width (degrees)", 45.0, minimum=1.0, maximum=90.0),
            _opt("OUTPUT_TYPE", "Output", ["Corrected grid", "Noise estimate"]),
        ],
    ),
    dict(
        id="reduction_to_pole", name="Reduction to the pole",
        group=GROUP_GRAV, fft=True, compute=_c_rtp, operation="Reduction to pole",
        help="Converts magnetic data to what it would be if measured at the magnetic pole.",
        params=[_dbl("INCLINATION", "Inclination (degrees)", -60.0),
                _dbl("DECLINATION", "Declination (degrees)", 0.0)],
    ),
    dict(
        id="reduction_to_pole_low_latitude", name="Reduction to the pole (low latitude)",
        group=GROUP_GRAV, fft=True, compute=_c_rtp_ll,
        operation="Low-latitude reduction to pole",
        help=(
            "Reduction to the pole for low magnetic latitudes by the Taylor series "
            "iterative method (Hao et al. 2018): the grid is reduced to the equator, "
            "then to the pole with a stable truncated operator applied iteratively. "
            "More iterations sharpen the result but bring back stripes along the "
            "declination; 0 chooses the number that best matches the analytic signal."
        ),
        params=[_dbl("INCLINATION", "Inclination (degrees)", 5.0),
                _dbl("DECLINATION", "Declination (degrees)", 0.0),
                _int("ITERATIONS", "Iterations (0 = automatic)", 0, minimum=0, maximum=1000)],
    ),
    dict(
        id="reduction_to_equator", name="Reduction to the equator",
        group=GROUP_GRAV, fft=True, compute=_c_rte, operation="Reduction to equator",
        help="Converts magnetic data to what it would be if measured at the magnetic equator.",
        params=[_dbl("INCLINATION", "Inclination (degrees)", -60.0),
                _dbl("DECLINATION", "Declination (degrees)", 0.0)],
    ),
    dict(
        id="continuation", name="Upward / downward continuation",
        group=GROUP_GRAV, fft=True, compute=_c_continuation, operation="Continuation",
        help="Continues the field up or down by the given height (map units; for "
             "geographic grids the height is roughly converted to degrees).",
        params=[_opt("DIRECTION", "Direction", ["up", "down"]),
                _dbl("HEIGHT", "Height", 500.0, minimum=0.0)],
    ),
    dict(
        id="vertical_integration", name="Vertical integration (pseudo-gravity)",
        group=GROUP_GRAV, fft=True, compute=_c_vint, operation="Vertical integration",
        help="Vertical integration; applied to an RTP/RTE grid it gives the pseudo-gravity.",
        params=[],
    ),
    dict(
        id="apparent_susceptibility", name="Apparent susceptibility",
        group=GROUP_PETRO, fft=True, compute=_c_app_sus,
        operation="Apparent susceptibility",
        help=(
            "Apparent susceptibility (SI) from a magnetic anomaly grid (nT, IGRF removed) "
            "reduced to the pole (or the equator, approximately), for a layer of vertical "
            "prisms magnetised by induction in the ambient field. The grid's absolute level "
            "is not used: the result is a variation about the background susceptibility."
        ),
        params=[
            _dbl("FIELD", "Ambient (IGRF) field strength (nT)", 50000.0, minimum=0.0),
            _dbl("DEPTH", "Depth from the sensor to the top of the layer (m)", 0.0, minimum=0.0),
            _dbl("THICKNESS", "Layer thickness (m, 0 = infinite)", 0.0, minimum=0.0),
            _dbl("BACKGROUND", "Background susceptibility (SI)", 0.001),
            _dbl("HIGH_CUT", "High-cut wavelength (m, 0 = 4 x the cell size, -1 = none)", 0.0),
        ],
    ),
    dict(
        id="apparent_density", name="Apparent density",
        group=GROUP_PETRO, fft=True, compute=_c_app_dens,
        operation="Apparent density",
        help=(
            "Apparent density (g/cc) from a Bouguer gravity grid, for a layer of vertical "
            "prisms. The grid's absolute level (datum, regional) is not used: the result is "
            "a variation about the background density."
        ),
        params=[
            _opt("UNITS", "Gravity grid units", ["mGal", "um/s2 (gu)"]),
            _dbl("DEPTH", "Depth from the sensor to the top of the layer (m)", 0.0, minimum=0.0),
            _dbl("THICKNESS", "Layer thickness (m)", 5000.0, minimum=0.0),
            _dbl("BACKGROUND", "Background density (g/cc)", 2.67),
            _dbl("HIGH_CUT", "High-cut wavelength (m, 0 = 4 x the cell size, -1 = none)", 0.0),
        ],
    ),
    dict(
        id="remove_regional", name="Remove regional trend",
        group=GROUP_FREQ, fft=False, compute=_c_regional, operation="Remove regional",
        help="Removes a 1st order (dipping plane) or 2nd order (parabolic) regional.",
        params=[_opt("ORDER", "Polynomial order", ["1st order", "2nd order"])],
    ),
    dict(
        id="band_pass", name="Band pass filter",
        group=GROUP_FREQ, fft=True, compute=_c_bandpass, operation="Band pass",
        help="Keeps wavelengths between the low and high cut (map units) with a cosine "
             "taper of the given transition width.",
        params=[_dbl("LOW_CUT", "Low cut wavelength", 50000.0),
                _dbl("HIGH_CUT", "High cut wavelength", 5000.0),
                _dbl("WIDTH", "Transition width", 5000.0, minimum=0.0)],
    ),
    dict(
        id="high_low_pass", name="High / low pass filter",
        group=GROUP_FREQ, fft=True, compute=_c_highlow, operation="High/Low pass",
        help="Removes wavelengths shorter (low pass) or longer (high pass) than the cutoff.",
        params=[_opt("TYPE", "Pass", ["Low", "High"]),
                _dbl("CUTOFF", "Cutoff wavelength", 5000.0),
                _dbl("WIDTH", "Transition width", 5000.0, minimum=0.0)],
    ),
    dict(
        id="automatic_gain_control", name="Automatic gain control",
        group=GROUP_FREQ, fft=False, compute=_c_agc, operation="Automatic gain control",
        help="Divides the grid by its RMS in a moving window.",
        params=[_int("WINDOW", "Window size (pixels)", 10, minimum=1)],
    ),
    dict(
        id="derivative", name="Derivative",
        group=GROUP_GRAD, fft=True, compute=_c_derivative, operation="Derivative",
        help="x, y or vertical (z) derivative of the given power.",
        params=[_opt("DIRECTION", "Direction", ["z", "x", "y"]),
                _dbl("POWER", "Power", 1.0, minimum=0.0)],
    ),
    dict(
        id="tilt_angle", name="Tilt angle",
        group=GROUP_GRAD, fft=True, compute=_c_tilt, operation="Tilt angle",
        help="Tilt angle of the field, in degrees (-90 to 90).", params=[],
    ),
    dict(
        id="analytic_signal", name="Analytic signal",
        group=GROUP_GRAD, fft=True, compute=_c_analytic, operation="Analytic signal",
        help="Total amplitude of the gradients.", params=[],
    ),
    dict(
        id="total_horizontal_gradient", name="Total horizontal gradient",
        group=GROUP_GRAD, fft=True, compute=_c_thg, operation="Total horizontal gradient",
        help="Total horizontal gradient of the field.", params=[],
    ),
    dict(
        id="mean_filter", name="Mean filter",
        group=GROUP_CONV, fft=False, compute=_c_mean, operation="Mean filter",
        help="Moving-window mean.",
        params=[_int("SIZE", "Filter size (pixels, odd number)", 3, minimum=3)],
    ),
    dict(
        id="median_filter", name="Median filter",
        group=GROUP_CONV, fft=False, compute=_c_median, operation="Median filter",
        help="Moving-window median.",
        params=[_int("SIZE", "Filter size (pixels, odd number)", 3, minimum=3)],
    ),
    dict(
        id="gaussian_filter", name="Gaussian filter",
        group=GROUP_CONV, fft=False, compute=_c_gaussian, operation="Gaussian filter",
        help="Gaussian smoothing.",
        params=[_dbl("SIGMA", "Sigma (pixels)", 1.0, minimum=0.0)],
    ),
    dict(
        id="directional_filter", name="Directional filter",
        group=GROUP_CONV, fft=False, compute=_c_directional, operation="Directional filter",
        help="3 x 3 directional (edge) filter.",
        params=[_opt("DIRECTION", "Direction", _DIRECTIONS)],
    ),
    dict(
        id="sun_shading", name="Sun shading",
        group=GROUP_CONV, fft=False, compute=_c_sunshade, operation="Sun shading",
        help="Shaded relief from the given sun azimuth and zenith (altitude) angles.",
        params=[_dbl("AZIMUTH", "Sun azimuth (degrees; 0 = north)", 45.0),
                _dbl("ZENITH", "Sun zenith / altitude angle (degrees)", 45.0),
                _bool("RELIEF", "Use the softer GRASS-like shading", True)],
    ),
    dict(
        id="windowed_statistic", name="Windowed statistic",
        group=GROUP_STATS, fft=False, compute=_c_window_stat, operation="Windowed statistic",
        help="A statistic of the values in a moving window.",
        params=[_opt("STATISTIC", "Statistic",
                     ["min", "max", "std", "variance", "skewness", "kurtosis"], default=3),
                _int("WINDOW", "Window size (pixels)", 5, minimum=3)],
    ),
    dict(
        id="local_anisotropy", name="Local anisotropy",
        group=GROUP_STATS, fft=False, compute=_c_anisotropy, operation="Local anisotropy",
        help="How directional the grid is in a moving window: the magnitude (0 = no preferred "
             "direction, 1 = strongly linear) or the dominant orientation (degrees, -90 to 90).",
        params=[_int("WINDOW", "Window size (pixels)", 7, minimum=3),
                _opt("WINDOW_TYPE", "Window", ["gaussian", "box"]),
                _opt("OUTPUT_TYPE", "Output", ["Magnitude", "Orientation"])],
    ),
    dict(
        id="anisotropy_chain_length", name="Anisotropy chain length",
        group=GROUP_STATS, fft=False, compute=_c_chain_length, operation="Anisotropy chain length",
        help="Scores each strongly directional cell by the size of the connected chain of cells "
             "with similar orientation that it belongs to.",
        params=[_int("WINDOW", "Window size (pixels)", 7, minimum=3),
                _opt("WINDOW_TYPE", "Window", ["gaussian", "box"]),
                _dbl("ANISOTROPY_THRESHOLD", "Anisotropy threshold (0-1)", 0.3, minimum=0.0, maximum=1.0),
                _dbl("ANGLE_TOLERANCE", "Orientation tolerance (degrees)", 22.5, minimum=0.0),
                _int("SEARCH_RADIUS", "Search radius (pixels)", 2, minimum=1)],
    ),
    dict(
        id="anisotropy_streamline_length", name="Anisotropy streamline length",
        group=GROUP_STATS, fft=False, compute=_c_streamline_length,
        operation="Anisotropy streamline length",
        help="Traces each strongly directional cell forwards and backwards along the orientation "
             "field and gives the length of the path (pixels).",
        params=[_int("WINDOW", "Window size (pixels)", 7, minimum=3),
                _opt("WINDOW_TYPE", "Window", ["gaussian", "box"]),
                _dbl("ANISOTROPY_THRESHOLD", "Anisotropy threshold (0-1)", 0.3, minimum=0.0, maximum=1.0),
                _dbl("ANGLE_TOLERANCE", "Orientation tolerance (degrees)", 22.5, minimum=0.0),
                _int("MAX_STEPS", "Maximum steps", 200, minimum=1)],
    ),
    dict(
        id="threshold_to_nan", name="Threshold to NaN",
        group=GROUP_UTIL, fft=False, compute=_c_nan, operation="Threshold to NaN",
        help="Sets values above, below, or between two thresholds to no data.",
        params=[_opt("CONDITION", "Set to NaN", ["above", "below", "between"]),
                _dbl("ABOVE", "Above (upper threshold)", 1.0),
                _dbl("BELOW", "Below (lower threshold)", -1.0)],
    ),
]
