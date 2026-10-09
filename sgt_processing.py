"""
SGTool filters as QGIS Processing algorithms.

Registered by the plugin as the "SGTool" provider, so every filter is available
from the Processing Toolbox, the graphical modeler, the batch-processing
dialog, and Python scripts (processing.run("sgtool:derivative", {...})).

The algorithms call the same calculation code as the plugin dialog
(GeophysicalProcessor, ConvolutionFilter, SG_Util, SpatialStats, ...), so a
result from the toolbox matches the one from the dialog for the same settings.
Like the dialog, each result is a GeoTIFF with its provenance (operation,
parameters, source file and its history) embedded.

Every algorithm can be cancelled from the Processing dialog. Cancelling takes
effect between steps and inside long loops (Euler deconvolution, B-spline
levels); a single FFT call cannot be interrupted part-way.
"""

import gc
import os

import numpy as np
from osgeo import gdal, osr

from qgis.PyQt.QtGui import QIcon
from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingContext,
    QgsProcessingException,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterEnum,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterField,
    QgsProcessingParameterFileDestination,
    QgsProcessingParameterNumber,
    QgsProcessingParameterRasterDestination,
    QgsProcessingOutputNumber,
    QgsProcessingParameterRasterLayer,
    QgsProcessingProvider,
)

from .calcs.ConvolutionFilter import ConvolutionFilter
from .calcs.GeophysicalProcessor import GeophysicalProcessor
from .calcs.PCAICA import PCAICA
from .calcs.SG_Util import SG_Util
from .calcs.SpatialStats import SpatialStats
from .calcs.euler.euler_python import euler_deconv
from .calcs.saga_mba_gridding import mba_gridding
from .calcs.sgt_cancel import OperationCancelled
from .calcs.sgt_metadata import write_sgt_metadata

PROVIDER_ID = "sgtool"


# ---------------------------------------------------------------- QGIS 3 / 4
def _enum(owner, scoped, name):
    """owner.<scoped>.<name> (QGIS 4) or owner.<name> (QGIS 3)."""
    holder = getattr(owner, scoped, None)
    if holder is not None and hasattr(holder, name):
        return getattr(holder, name)
    return getattr(owner, name)


def _number_type(name):
    return _enum(QgsProcessingParameterNumber, "Type", name)


def _point_source_type():
    return _enum(QgsProcessing, "SourceType", "TypeVectorPoint")


def _numeric_field_type():
    return _enum(QgsProcessingParameterField, "DataType", "Numeric")


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
        raise QgsProcessingException(f"Could not open the raster {path}")
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
        raise QgsProcessingException(f"Could not create {path}")
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
        return self._get("c", lambda: ConvolutionFilter(self.arr))

    @property
    def sg_util(self):
        return self._get("u", lambda: SG_Util(self.arr))

    @property
    def spatial(self):
        return self._get("s", lambda: SpatialStats(self.arr))

    def require_projected(self, what):
        if self.geographic:
            raise QgsProcessingException(
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
    "threshold_to_nan": lambda p: "_Clean",
    "pca": lambda p: "_PCA",
    "ica": lambda p: "_ICA",
}


def auto_output_path(source_path, suffix, extension=".tif"):
    """Where the dialog would save a result: next to the source, named after it
    plus the processing step."""
    base, _ext = os.path.splitext(source_path)
    return base + suffix + extension


def _check_length(g, length):
    """Lengths in degrees for geographic grids (as the plugin dialog checks)."""
    if g.geographic and length > 100:
        raise QgsProcessingException(
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
    return g.processor.tilt_angle(g.arr, buffer_size=g.buffer)


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
        help="Tilt angle of the field.", params=[],
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
        id="threshold_to_nan", name="Threshold to NaN",
        group=GROUP_UTIL, fft=False, compute=_c_nan, operation="Threshold to NaN",
        help="Sets values above, below, or between two thresholds to no data.",
        params=[_opt("CONDITION", "Set to NaN", ["above", "below", "between"]),
                _dbl("ABOVE", "Above (upper threshold)", 1.0),
                _dbl("BELOW", "Below (lower threshold)", -1.0)],
    ),
]


def _add_parameter(alg, spec):
    kind = spec["kind"]
    name, label = spec["name"], spec["label"]
    if kind in ("double", "int"):
        param = QgsProcessingParameterNumber(
            name, label,
            _number_type("Double" if kind == "double" else "Integer"),
            spec["default"],
            spec.get("optional", False),
        )
        if spec["minimum"] is not None:
            param.setMinimum(spec["minimum"])
        if spec["maximum"] is not None:
            param.setMaximum(spec["maximum"])
    elif kind == "enum":
        param = QgsProcessingParameterEnum(
            name, label, spec["options"], False, spec["default"]
        )
    else:
        param = QgsProcessingParameterBoolean(name, label, spec["default"])
    alg.addParameter(param)


def _read_values(alg, specs, parameters, context):
    """Parameter values by name (enums as their option text) and a version of
    the same for the provenance record."""
    values = {}
    for spec in specs:
        kind, name = spec["kind"], spec["name"]
        if kind == "double":
            values[name] = alg.parameterAsDouble(parameters, name, context)
        elif kind == "int":
            values[name] = alg.parameterAsInt(parameters, name, context)
        elif kind == "enum":
            values[name] = spec["options"][alg.parameterAsEnum(parameters, name, context)]
        else:
            values[name] = alg.parameterAsBoolean(parameters, name, context)
    return values


class _SGToolAlgorithm(QgsProcessingAlgorithm):
    """Shared plumbing: names, groups, icon, help."""

    def __init__(self, spec):
        super().__init__()
        self.spec = spec

    def createInstance(self):
        return type(self)(self.spec)

    def name(self):
        return self.spec["id"]

    def displayName(self):
        return self.spec["name"]

    def group(self):
        return self.spec["group"][1]

    def groupId(self):
        return self.spec["group"][0]

    def tags(self):
        return ["sgtool", "geophysics", "potential field", "grid"]

    OUTPUT_HELP = (
        "\n\nOutput: choose a file, or leave it empty to save next to the input "
        "named after it plus the processing step (for example grid_d1z.tif), "
        "and add it to the project, as the SGTool dialog does."
    )

    def helpString(self):
        return self.spec["help"] + self.OUTPUT_HELP

    def shortHelpString(self):
        return self.spec["help"] + self.OUTPUT_HELP

    def _auto_target(self, parameters, context):
        """The automatic output path for these settings, or None. Subclasses
        that save a grid next to their input provide it."""
        return None

    def prepareAlgorithm(self, parameters, context, feedback):
        """Runs on the main thread before the calculation. If the result will
        replace a grid that is open in QGIS, remove that layer from the project
        first (as the SGTool dialog does) so the file is not locked."""
        output_key = getattr(self, "OUTPUT", None)  # (some algorithms have none)
        if output_key and parameters.get(output_key) in (None, ""):
            try:
                target = self._auto_target(parameters, context)
            except Exception:
                target = None
            if target and os.path.exists(target):
                project = context.project()
                if project is not None:
                    wanted = os.path.normcase(os.path.abspath(target))
                    for lyr in list(project.mapLayers().values()):
                        source = lyr.source().split("|")[0]
                        if os.path.normcase(os.path.abspath(source)) == wanted:
                            project.removeMapLayer(lyr.id())
                    gc.collect()  # let go of the file handle
        return True

    def _resolve_output(self, parameters, context, auto_path, layer_name, raster=True):
        """(path to write, chosen automatically?)

        A file picked by the user is used as given. If the output is left
        empty the result goes to auto_path (the dialog's naming) and, for
        grids, is queued to be loaded into the project under layer_name.
        An existing file of that name is replaced, as in the dialog."""
        if raster:
            chosen = self.parameterAsOutputLayer(parameters, self.OUTPUT, context)
        else:
            chosen = self.parameterAsFileOutput(parameters, self.OUTPUT, context)
        if chosen:
            return chosen, False
        if os.path.exists(auto_path):
            try:
                os.remove(auto_path)
            except OSError:
                raise QgsProcessingException(
                    f"{auto_path} already exists and could not be replaced. If it is "
                    "open in QGIS, remove that layer first (or choose another output file)."
                )
            for stale in (auto_path + ".aux.xml", auto_path + ".sgt.xml"):
                if os.path.exists(stale):
                    try:
                        os.remove(stale)
                    except OSError:
                        pass
        project = context.project()
        if raster and project is not None:
            context.addLayerToLoadOnCompletion(
                auto_path,
                QgsProcessingContext.LayerDetails(layer_name, project, self.OUTPUT),
            )
        return auto_path, True

    def icon(self):
        return QIcon(os.path.join(os.path.dirname(os.path.abspath(__file__)), "icon.png"))


class RasterFilterAlgorithm(_SGToolAlgorithm):
    """A filter from FILTER_SPECS: raster in, raster out."""

    INPUT = "INPUT"
    BUFFER = "BUFFER"
    OUTPUT = "OUTPUT"

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer(self.INPUT, "Input grid"))
        for spec in self.spec["params"]:
            _add_parameter(self, spec)
        if self.spec["fft"]:
            buf = QgsProcessingParameterNumber(
                self.BUFFER, "FFT buffer (pixels, 0 = automatic)",
                _number_type("Integer"), 0, False,
            )
            buf.setMinimum(0)
            self.addParameter(buf)
        # optional, not created by default: left empty it is saved next to the
        # input with the dialog's naming (see _resolve_output)
        self.addParameter(
            QgsProcessingParameterRasterDestination(
                self.OUTPUT, "Output grid (optional: default is next to the input)",
                None, True, False,
            )
        )

    def _auto_target(self, parameters, context):
        layer = self.parameterAsRasterLayer(parameters, self.INPUT, context)
        if layer is None:
            return None
        values = _read_values(self, self.spec["params"], parameters, context)
        return auto_output_path(
            layer.source().split("|")[0], SUFFIXES[self.spec["id"]](values)
        )

    def processAlgorithm(self, parameters, context, feedback):
        layer = self.parameterAsRasterLayer(parameters, self.INPUT, context)
        if layer is None:
            raise QgsProcessingException("Select an input grid")
        source_path = layer.source().split("|")[0]
        values = _read_values(self, self.spec["params"], parameters, context)
        suffix = SUFFIXES[self.spec["id"]](values)
        buffer = (
            self.parameterAsInt(parameters, self.BUFFER, context) if self.spec["fft"] else 0
        )

        try:
            feedback.setProgress(5)
            arr, gt, projection = read_grid(source_path)
            if feedback.isCanceled():
                return {}
            grid = GridContext(arr, gt, layer, buffer)
            feedback.setProgress(20)
            feedback.pushInfo(f"{self.spec['name']}: {arr.shape[1]} x {arr.shape[0]} cells")
            result = self.spec["compute"](grid, values)
            if feedback.isCanceled():
                return {}
            if result is None:
                raise QgsProcessingException("Nothing was calculated for these settings")
            feedback.setProgress(85)
            out_path, _auto = self._resolve_output(
                parameters, context,
                auto_output_path(source_path, suffix), layer.name() + suffix,
            )
            write_grid(out_path, np.asarray(result), gt, projection)
        except OperationCancelled:
            return {}
        write_sgt_metadata(
            out_path, source_path, self.spec["operation"], values, _plugin_version()
        )
        feedback.setProgress(100)
        return {self.OUTPUT: out_path}


class EulerAlgorithm(_SGToolAlgorithm):
    """Euler deconvolution solutions for one structural index, as a CSV."""

    INPUT = "INPUT"
    SI = "STRUCTURAL_INDEX"
    WINDOW = "WINDOW"
    KEEP = "KEEP"
    OUTPUT = "OUTPUT"
    SI_OPTIONS = ["0 - prism / contact", "1 - line of poles", "2 - single pole", "3 - dipole"]
    SI_VALUES = [0.001, 1, 2, 3]  # 0.001 stands in for 0, as the plugin does
    HEADER = "y_source, x_source, z_source, base_level, std_dfdz"

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer(self.INPUT, "Input grid"))
        self.addParameter(
            QgsProcessingParameterEnum(self.SI, "Structural index", self.SI_OPTIONS, False, 1)
        )
        window = QgsProcessingParameterNumber(
            self.WINDOW, "Moving window size (pixels)", _number_type("Integer"), 10, False
        )
        window.setMinimum(3)
        self.addParameter(window)
        keep = QgsProcessingParameterNumber(
            self.KEEP, "Fraction of solutions kept (those with the largest std of df/dz)",
            _number_type("Double"), 0.1, False,
        )
        keep.setMinimum(0.0)
        keep.setMaximum(1.0)
        self.addParameter(keep)
        self.addParameter(
            QgsProcessingParameterFileDestination(
                self.OUTPUT, "Euler solutions (optional: default is next to the input)",
                "CSV files (*.csv)", None, True, False,
            )
        )

    def processAlgorithm(self, parameters, context, feedback):
        layer = self.parameterAsRasterLayer(parameters, self.INPUT, context)
        if layer is None:
            raise QgsProcessingException("Select an input grid")
        source_path = layer.source().split("|")[0]
        si_index = self.parameterAsEnum(parameters, self.SI, context)
        SI = self.SI_VALUES[si_index]
        winsize = self.parameterAsInt(parameters, self.WINDOW, context)
        filt = self.parameterAsDouble(parameters, self.KEEP, context)

        arr, gt, _projection = read_grid(source_path)
        grid = GridContext(arr, gt, layer, 0)
        grid.require_projected("Euler deconvolution")
        try:
            data, _mask = grid.processor.fill_nan(arr)
            shape = (data.shape[0], data.shape[1])
            xmin, ymin, xmax, ymax = grid.extent
            area = [ymin, ymax, xmin, xmax]  # south, north, west, east
            rows, cols = data.shape
            x_coords_1d = np.linspace(area[2], area[3], cols)
            y_coords_1d = np.linspace(area[0], area[1], rows)
            yi, xi = np.meshgrid(x_coords_1d, y_coords_1d)
            zi = np.ones(shape)
            xi, yi, zi = xi.flatten(), yi.flatten(), zi.flatten()
            data = np.flipud(data).flatten()

            def callback(fraction=None):
                if feedback.isCanceled():
                    raise OperationCancelled()
                if fraction is not None:
                    feedback.setProgress(5 + 90 * float(fraction))

            result = euler_deconv(
                data, xi, yi, zi, shape, area, SI, winsize, filt, callback=callback
            )
            # same x-coordinate flip the plugin dialog applies to its solutions
            result[1, :] = area[0] - result[1, :] + area[1]
        except OperationCancelled:
            return {}

        # the dialog's name for these solutions: <grid>_estimates_SI_<n>, here as .csv
        out_path, _auto = self._resolve_output(
            parameters, context,
            auto_output_path(source_path, f"_estimates_SI_{si_index}", ".csv"),
            "", raster=False,
        )
        np.savetxt(out_path, result, delimiter=",", header=self.HEADER, comments="")
        write_sgt_metadata(
            out_path, source_path, "Euler deconvolution",
            {
                "structural_index": SI,
                "window_size": winsize,
                "fraction_of_solutions_kept": filt,
                "columns": self.HEADER,
            },
            _plugin_version(),
        )
        feedback.setProgress(100)
        return {self.OUTPUT: out_path}


class ComponentAnalysisAlgorithm(_SGToolAlgorithm):
    """PCA or ICA of a multiband grid; the result has one band per component."""

    INPUT = "INPUT"
    COMPONENTS = "COMPONENTS"
    OUTPUT = "OUTPUT"

    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterRasterLayer(self.INPUT, "Input multiband grid")
        )
        n = QgsProcessingParameterNumber(
            self.COMPONENTS, "Number of components", _number_type("Integer"), 3, False
        )
        n.setMinimum(1)
        self.addParameter(n)
        self.addParameter(
            QgsProcessingParameterRasterDestination(
                self.OUTPUT, "Output grid (optional: default is next to the input)",
                None, True, False,
            )
        )

    def _auto_target(self, parameters, context):
        layer = self.parameterAsRasterLayer(parameters, self.INPUT, context)
        if layer is None:
            return None
        return auto_output_path(
            layer.source().split("|")[0], SUFFIXES[self.spec["id"]]({})
        )

    def processAlgorithm(self, parameters, context, feedback):
        try:
            import sklearn  # noqa: F401
        except ImportError:
            raise QgsProcessingException(
                "scikit-learn is required for PCA/ICA; install it with "
                "'pip3 install scikit-learn' in the QGIS Python environment"
            )
        layer = self.parameterAsRasterLayer(parameters, self.INPUT, context)
        if layer is None:
            raise QgsProcessingException("Select an input grid")
        source_path = layer.source().split("|")[0]
        n = self.parameterAsInt(parameters, self.COMPONENTS, context)
        suffix = SUFFIXES[self.spec["id"]]({})
        # the analysis writes the file itself, so the output is settled first
        out_path, _auto = self._resolve_output(
            parameters, context,
            auto_output_path(source_path, suffix), layer.name() + suffix,
        )
        feedback.setProgress(10)
        analysis = PCAICA([[0.0]])  # the methods work from the file paths
        if self.spec["id"] == "pca":
            ok = analysis.pca_with_nans(source_path, out_path, n)[0] is not None
        else:
            ok = analysis.ica_with_nans(source_path, out_path, n)[0] is not None
        if feedback.isCanceled():
            return {}
        if not ok:
            raise QgsProcessingException("The analysis produced no components")
        write_sgt_metadata(
            out_path, source_path, self.spec["operation"], {"n_components": n},
            _plugin_version(),
        )
        feedback.setProgress(100)
        return {self.OUTPUT: out_path}


class BSplineGriddingAlgorithm(_SGToolAlgorithm):
    """Multilevel B-spline gridding of scattered points."""

    INPUT = "INPUT"
    FIELD = "FIELD"
    CELL_SIZE = "CELL_SIZE"
    EPSILON = "EPSILON"
    LEVELS = "LEVELS"
    IGNORE_BELOW = "IGNORE_BELOW"
    OUTPUT = "OUTPUT"

    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.INPUT, "Input points", [_point_source_type()]
            )
        )
        self.addParameter(
            QgsProcessingParameterField(
                self.FIELD, "Value field", parentLayerParameterName=self.INPUT,
                type=_numeric_field_type(),
            )
        )
        cell = QgsProcessingParameterNumber(
            self.CELL_SIZE, "Cell size (map units)", _number_type("Double"), 100.0, False
        )
        cell.setMinimum(0.0)
        self.addParameter(cell)
        self.addParameter(
            QgsProcessingParameterNumber(
                self.EPSILON, "Threshold error (epsilon)", _number_type("Double"), 0.0001, False
            )
        )
        levels = QgsProcessingParameterNumber(
            self.LEVELS, "Maximum number of levels", _number_type("Integer"), 11, False
        )
        levels.setMinimum(1)
        self.addParameter(levels)
        self.addParameter(
            QgsProcessingParameterNumber(
                self.IGNORE_BELOW, "Ignore points with a value below (leave empty to use all)",
                _number_type("Double"), None, True,
            )
        )
        self.addParameter(
            QgsProcessingParameterRasterDestination(
                self.OUTPUT, "Output grid (optional: default is next to the input points)",
                None, True, False,
            )
        )

    def _points_naming(self, parameters, context):
        """(points layer name, its file or None, automatic output path): the
        dialog's name <points>_<field>_bspline, next to the points file (the
        temp folder if the points are not a file, e.g. a memory layer)."""
        field = self.parameterAsString(parameters, self.FIELD, context)
        points_layer = self.parameterAsVectorLayer(parameters, self.INPUT, context)
        layer_name, file_path = "points", None
        if points_layer is not None:
            layer_name = points_layer.name()
            candidate = points_layer.source().split("|")[0]
            if os.path.isfile(candidate):
                file_path = candidate
        suffix = f"_{field}_bspline"
        if file_path:
            return layer_name, file_path, auto_output_path(file_path, suffix)
        import tempfile

        return layer_name, None, os.path.join(
            tempfile.gettempdir(), layer_name + suffix + ".tif"
        )

    def _auto_target(self, parameters, context):
        return self._points_naming(parameters, context)[2]

    def processAlgorithm(self, parameters, context, feedback):
        source = self.parameterAsSource(parameters, self.INPUT, context)
        if source is None:
            raise QgsProcessingException("Select an input points layer")
        field = self.parameterAsString(parameters, self.FIELD, context)
        cell_size = self.parameterAsDouble(parameters, self.CELL_SIZE, context)
        epsilon = self.parameterAsDouble(parameters, self.EPSILON, context)
        level_max = self.parameterAsInt(parameters, self.LEVELS, context)
        ignore_below = (
            self.parameterAsDouble(parameters, self.IGNORE_BELOW, context)
            if parameters.get(self.IGNORE_BELOW) not in (None, "")
            else None
        )

        xs, ys, zs = [], [], []
        for feat in source.getFeatures():
            geom = feat.geometry()
            if geom is None or geom.isEmpty():
                continue
            val = feat[field]
            if val is None:
                continue
            try:
                zval = float(val)
            except (TypeError, ValueError):
                continue
            points = geom.asMultiPoint() if geom.isMultipart() else [geom.asPoint()]
            for pt in points:
                xs.append(pt.x())
                ys.append(pt.y())
                zs.append(zval)
            if feedback.isCanceled():
                return {}
        xs, ys, zs = (np.asarray(v, dtype=float) for v in (xs, ys, zs))
        if ignore_below is not None:
            keep = zs >= ignore_below
            xs, ys, zs = xs[keep], ys[keep], zs[keep]
        if xs.size < 3:
            raise QgsProcessingException("Need at least 3 valid points with a numeric value")

        # grid geometry as the plugin dialog and SAGA: lower-left node snapped down
        xmin = float(np.floor(xs.min() / cell_size) * cell_size)
        ymin = float(np.floor(ys.min() / cell_size) * cell_size)
        nx = 2 + int((xs.max() - xmin) / cell_size)
        ny = 2 + int((ys.max() - ymin) / cell_size)
        if nx < 2 or ny < 2:
            raise QgsProcessingException("Cell size too large for the point extent")

        def callback(fraction=None):
            if feedback.isCanceled():
                raise OperationCancelled()
            if fraction is not None:
                feedback.setProgress(5 + 85 * float(fraction))

        try:
            grid = mba_gridding(
                xs, ys, zs, cellsize=cell_size, xmin=xmin, ymin=ymin, nx=nx, ny=ny,
                epsilon=epsilon, level_max=level_max, refinement=False, callback=callback,
            )
        except OperationCancelled:
            return {}
        grid = np.flipud(grid)  # row 0 = ymin from the gridder; GeoTIFF is north-up
        gt = (xmin - cell_size / 2.0, cell_size, 0.0,
              ymin + (ny - 0.5) * cell_size, 0.0, -cell_size)

        # source file and name of the points layer, for provenance and naming
        layer_name, source_layer_path, auto_path = self._points_naming(parameters, context)
        out_path, _auto = self._resolve_output(
            parameters, context, auto_path, layer_name + f"_{field}_bspline"
        )
        write_grid(out_path, grid, gt, source.sourceCrs().toWkt())

        params = {
            "data_field": field, "cell_size": cell_size, "epsilon": epsilon,
            "max_levels": level_max, "points_used": int(xs.size),
        }
        if ignore_below is not None:
            params["ignore_values_below"] = ignore_below
        write_sgt_metadata(
            out_path, source_layer_path, "Multilevel B-spline gridding", params,
            _plugin_version(),
        )
        feedback.setProgress(100)
        return {self.OUTPUT: out_path}


COMPONENT_SPECS = [
    dict(id="pca", name="Principal component analysis", group=GROUP_MULTI,
         operation="Principal component analysis",
         help="Principal components of a multiband grid (needs scikit-learn)."),
    dict(id="ica", name="Independent component analysis", group=GROUP_MULTI,
         operation="Independent component analysis",
         help="Independent components of a multiband grid (needs scikit-learn)."),
]

EULER_SPEC = dict(
    id="euler_deconvolution", name="Euler deconvolution", group=GROUP_EULER,
    help=(
        "Euler deconvolution of a magnetic grid for one structural index. Writes a "
        "CSV of solutions (y, x, depth, base level and the standard deviation of "
        "df/dz of the window each came from)."
    ),
)

BSPLINE_SPEC = dict(
    id="bspline_gridding", name="Multilevel B-spline gridding", group=GROUP_GRID,
    help="Grids scattered points with the multilevel B-spline method (after SAGA).",
)


class SGToolProvider(QgsProcessingProvider):
    """The SGTool Processing provider."""

    def id(self):
        return PROVIDER_ID

    def name(self):
        return "SGTool"

    def longName(self):
        return "SGTool: structural geophysics grid tools"

    def icon(self):
        return QIcon(os.path.join(os.path.dirname(os.path.abspath(__file__)), "icon.png"))

    def loadAlgorithms(self):
        for spec in FILTER_SPECS:
            self.addAlgorithm(RasterFilterAlgorithm(spec))
        self.addAlgorithm(EulerAlgorithm(EULER_SPEC))
        for spec in COMPONENT_SPECS:
            self.addAlgorithm(ComponentAnalysisAlgorithm(spec))
        self.addAlgorithm(BSplineGriddingAlgorithm(BSPLINE_SPEC))
        # imported here: the replay module itself builds on this one
        from .sgt_replay_algorithm import ReplayHistoryAlgorithm, REPLAY_SPEC

        self.addAlgorithm(ReplayHistoryAlgorithm(REPLAY_SPEC))
