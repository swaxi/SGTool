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

from .calcs.PCAICA import PCAICA
from .calcs.euler.euler_python import euler_deconv
from .calcs.layer_shim import LayerShim
from .calcs.mba_suggest import describe as describe_mba_suggestion, suggest_mba_parameters
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
# The filters and grid helpers live in sgt_filters.py (no QGIS dependency, so the
# ArcGIS Pro GUI can use them too); a failure there is a FilterError, shown here as
# a Processing error with the same message.
from .sgt_filters import (  # noqa: E402,F401
    FilterError,
    FILTER_SPECS,
    SUFFIXES,
    GridContext,
    auto_output_path,
    _plugin_version,
    _fmt,
    _dbl,
    _int,
    _opt,
    _bool,
    _DIRECTIONS,
    GROUP_GRAV,
    GROUP_FREQ,
    GROUP_GRAD,
    GROUP_CONV,
    GROUP_STATS,
    GROUP_UTIL,
    GROUP_GRID,
    GROUP_EULER,
    GROUP_MULTI,
    read_grid as _read_grid,
    write_grid as _write_grid,
)


def read_grid(path):
    """Band 1 of a raster as a north-up float array (NaN for no data), with its
    geotransform and projection, as the plugin reads grids."""
    try:
        return _read_grid(path)
    except FilterError as e:
        raise QgsProcessingException(str(e))


def write_grid(path, arr, gt, projection):
    try:
        _write_grid(path, arr, gt, projection)
    except FilterError as e:
        raise QgsProcessingException(str(e))


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
                self.OUTPUT, "Output grid (optional: default is standard SGTool file naming convention in same directory as input file)",
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
            try:
                result = self.spec["compute"](grid, values)
            except FilterError as e:
                raise QgsProcessingException(str(e))
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
                self.OUTPUT, "Euler solutions (optional: default is standard SGTool file naming convention in same directory as input file)",
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
                self.OUTPUT, "Output grid (optional: default is standard SGTool file naming convention in same directory as input file)",
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
    AUTO = "AUTO"
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
            QgsProcessingParameterBoolean(
                self.AUTO,
                "Choose the cell size and number of levels from the point spacing "
                "(the two settings above are then ignored)",
                False,
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.IGNORE_BELOW, "Ignore points with a value below (leave empty to use all)",
                _number_type("Double"), None, True,
            )
        )
        self.addParameter(
            QgsProcessingParameterRasterDestination(
                self.OUTPUT, "Output grid (optional: default is standard SGTool file naming convention in same directory as input points)",
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

        auto = self.parameterAsBoolean(parameters, self.AUTO, context)
        if auto:
            try:
                suggestion = suggest_mba_parameters(
                    xs, ys, geographic=source.sourceCrs().isGeographic()
                )
            except ValueError as e:
                raise QgsProcessingException(f"Cannot choose the settings automatically: {e}")
            cell_size, level_max = suggestion["cell_size"], suggestion["max_levels"]
            feedback.pushInfo(
                f"Automatic settings: cell size {cell_size:g}, {level_max} levels. "
                + describe_mba_suggestion(suggestion)
            )

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
        if auto:
            params["settings_chosen_automatically"] = True
        if ignore_below is not None:
            params["ignore_values_below"] = ignore_below
        write_sgt_metadata(
            out_path, source_layer_path, "Multilevel B-spline gridding", params,
            _plugin_version(),
        )
        feedback.setProgress(100)
        return {self.OUTPUT: out_path}


class WormsAlgorithm(_SGToolAlgorithm):
    """BSDWorms: multilevel wavelet edges of a gravity or RTP/RTE + vertical integration grid."""

    INPUT = "INPUT"
    LEVELS = "LEVELS"
    BASE = "BASE"
    INCREMENT = "INCREMENT"
    SHAPEFILE = "SHAPEFILE"
    OUTPUT = "OUTPUT"

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer(self.INPUT, "Input grid"))
        levels = QgsProcessingParameterNumber(
            self.LEVELS, "Number of levels", _number_type("Integer"), 10, False
        )
        levels.setMinimum(1)
        self.addParameter(levels)
        base = QgsProcessingParameterNumber(
            self.BASE, "Base level: upward continuation height of the first level (map units)",
            _number_type("Double"), 1000.0, False,
        )
        base.setMinimum(0.0)
        self.addParameter(base)
        inc = QgsProcessingParameterNumber(
            self.INCREMENT, "Increment: spacing between levels (map units)",
            _number_type("Double"), 1000.0, False,
        )
        inc.setMinimum(0.0)
        self.addParameter(inc)
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.SHAPEFILE, "Also save the worms as a polyline shapefile (needs scikit-learn)", False
            )
        )
        self.addParameter(
            QgsProcessingParameterFileDestination(
                self.OUTPUT, "Worm points (optional: default is standard SGTool file naming convention in same directory as input file)",
                "CSV files (*.csv)", None, True, False,
            )
        )

    def _auto_target(self, parameters, context):
        return None  # nothing here replaces a grid that could be open in QGIS

    def processAlgorithm(self, parameters, context, feedback):
        layer = self.parameterAsRasterLayer(parameters, self.INPUT, context)
        if layer is None:
            raise QgsProcessingException("Select an input grid")
        source_path = layer.source().split("|")[0]
        num_levels = self.parameterAsInt(parameters, self.LEVELS, context)
        bottom_level = self.parameterAsDouble(parameters, self.BASE, context)
        delta_z = self.parameterAsDouble(parameters, self.INCREMENT, context)
        shapefile = self.parameterAsBoolean(parameters, self.SHAPEFILE, context)

        if layer.crs().isGeographic():
            raise QgsProcessingException(
                "Worms need a grid in a projected (metre-based) coordinate system, "
                "not latitude / longitude"
            )
        authid = layer.crs().authid()
        if not authid.upper().startswith("EPSG:"):
            raise QgsProcessingException(
                "The grid needs a coordinate system with an EPSG code (the worm shapefile is "
                "written in it): assign one to the layer first"
            )
        if shapefile:
            try:
                import sklearn  # noqa: F401
            except ImportError:
                raise QgsProcessingException(
                    "scikit-learn is required for the worm shapefile; install it with "
                    "'pip3 install scikit-learn' in the QGIS Python environment, or untick "
                    "the shapefile option"
                )

        arr, gt, _projection = read_grid(source_path)
        grid = GridContext(arr, gt, layer, 0)
        ext = layer.extent()
        shim = LayerShim((ext.xMinimum(), ext.yMinimum(), ext.xMaximum(), ext.yMaximum()),
                         authid, False)

        # the worms are written next to the grid under the dialog's names; a CSV chosen
        # elsewhere is moved there afterwards
        default_csv = auto_output_path(source_path, "_worms", ".csv")
        out_csv, _auto = self._resolve_output(
            parameters, context, default_csv, "", raster=False
        )
        base_, ext_ = os.path.splitext(source_path)
        padded = base_ + "_padded" + ext_
        shp_path = auto_output_path(source_path, "_worms", ".shp")

        def callback(fraction=None):
            if feedback.isCanceled():
                raise OperationCancelled()
            if fraction is not None:
                feedback.setProgress(5 + 90 * float(fraction))

        try:
            grid.processor.bsdwormer(
                arr, shim, source_path, num_levels, bottom_level, delta_z, shapefile,
                int(authid.split(":")[1]), callback=callback,
            )
        except OperationCancelled:
            return {}

        params = {
            "levels": num_levels, "bottom_level": bottom_level,
            "level_increment": delta_z, "shapefile": shapefile,
        }
        if os.path.abspath(out_csv) != os.path.abspath(default_csv) and os.path.exists(default_csv):
            os.replace(default_csv, out_csv)
        written = [
            (out_csv, "Worms (BSDWormer) points"),
            (padded, "Worms padded and tapered input grid"),
        ]
        if shapefile:
            written.append((shp_path, "Worms (BSDWormer) polylines"))
        for path, operation in written:
            if os.path.exists(path):
                write_sgt_metadata(path, source_path, operation, params, _plugin_version())
        feedback.pushInfo(
            f"Worms written to {out_csv}; the padded grid is {padded}"
            + (f" and the polylines {shp_path}" if shapefile else "")
        )
        feedback.setProgress(100)
        return {self.OUTPUT: out_csv}


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


WORMS_SPEC = dict(
    id="worms", name="BSDWorms (multilevel wavelet edges)", group=("wavelets", "Wavelets"),
    help=(
        "Builds wavelet-transform worms (Horowitz's bsdwormer): the grid is upward continued to "
        "a series of heights (levels), and the edges at each level are found from the wavelet "
        "transform of the horizontal gradient. The source depth of an edge is approximately half "
        "the upward continuation height of its level (Jacobsen, 1987), so edges seen at a 2000 m "
        "level come from sources at about 1000 m depth; this is a guide, not a precise depth. "
        "Use a gravity grid, or an RTP/RTE magnetic grid after vertical integration, in a "
        "metre-based coordinate system. Writes <grid>_worms.csv (and optionally <grid>_worms.shp) "
        "and a padded copy of the grid <grid>_padded.tif next to the input."
    ),
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
        self.addAlgorithm(WormsAlgorithm(WORMS_SPEC))
        # imported here: the replay module itself builds on this one
        from .sgt_replay_algorithm import ReplayHistoryAlgorithm, REPLAY_SPEC

        self.addAlgorithm(ReplayHistoryAlgorithm(REPLAY_SPEC))
