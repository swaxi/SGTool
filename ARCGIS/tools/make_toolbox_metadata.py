"""
Write the ArcGIS Pro help for the SGTool Python toolbox.

ArcGIS Pro reads a Python toolbox tool's help (the summary shown when you hover
over the tool, its usage notes, and the description shown for each parameter)
from <toolbox>.<ToolClass>.pyt.xml beside the .pyt. This script writes those
files, and SGTool.pyt.xml for the toolbox itself, from the text below; the
parameter names and labels are read from SGTool.pyt so they always match.

Run it with any Python 3 after changing a tool or its text:
    python ARCGIS/tools/make_toolbox_metadata.py
"""

import html
import os
import sys
import types
from datetime import datetime

HERE = os.path.dirname(os.path.abspath(__file__))
ARCGIS = os.path.dirname(HERE)
PYT = os.path.join(ARCGIS, "SGTool.pyt")

TOOLBOX_SUMMARY = (
    "Structural Geophysics Tool (SGTool) for ArcGIS Pro: potential-field grid processing "
    "(reductions, continuation, filters, derivatives), convolution filters and windowed "
    "statistics, apparent susceptibility and density, Euler deconvolution, worms, gridding "
    "and utilities. The same calculations as the SGTool QGIS plugin. The Launch SGTool GUI "
    "tool opens a multi-tab window with live preview; the other tools run one calculation each."
)

# shared parameter text
IN_RASTER = "Grid to process (a single-band raster; no-data cells are left as no-data in the result)."
OUT_RASTER = "Output raster (GeoTIFF), with the input's extent, cell size and coordinate system."
BUFFER = ("FFT buffer in pixels added around the grid before the Fourier transform, to reduce edge "
          "effects. 0 (the default) chooses it as in QGIS: the shorter side of the grid, up to 5000. "
          "A value greater than 0 is used as that upper limit instead.")
MAG_NOTE = "Best applied to a magnetic grid that has been reduced to the pole (RTP) or the equator (RTE)."

TOOLS = {
    "LaunchGUI": dict(
        summary="Open the SGTool window: the multi-tab interface of the SGTool QGIS plugin, running "
                "beside ArcGIS Pro, with live preview of the filters.",
        usage=["The window runs as a separate program using the active ArcGIS Pro Python environment; "
               "run Set Up Python Environment first if packages are missing.",
               "Results are written beside the input grid; run Add SGTool Outputs to Map to load them."],
        params={}),
    "LoadSGToolOutputs": dict(
        summary="Add the rasters written by the SGTool window since the last run to the first map of "
                "the current project.",
        usage=["Run it after processing in the SGTool window to see the results in ArcGIS Pro."],
        params={}),
    "SetUpPythonEnvironment": dict(
        summary="Create an ArcGIS Pro Python environment with the packages SGTool, Tomofast-x and "
                "GeoSubSampler need, and make it the environment ArcGIS Pro uses.",
        usage=["Creates a NEW environment as a copy of arcgispro-py3; it does not change arcgispro-py3 "
               "or any other environment you already have.",
               "Needs an internet connection and takes about 10 minutes. Restart ArcGIS Pro when it "
               "has finished to use the new environment.",
               "The same tool is in the SGTool, Tomofast-x and GeoSubSampler toolboxes; one run covers "
               "all three. Running it again only adds what is missing.",
               "It cannot change or delete the environment ArcGIS Pro is using at the time: to recreate "
               "it, switch to another environment in Project > Package Manager and restart first."],
        params={
            "env_name": "Name of the environment to create (or to add missing packages to, if it exists).",
            "activate": "Make the new environment the one ArcGIS Pro uses (restart ArcGIS Pro afterwards). "
                        "Untick to activate it later in Project > Package Manager.",
            "recreate": "Delete the environment and make it again from a fresh copy of arcgispro-py3, "
                        "instead of only adding what is missing.",
        }),
    "ReductionToPoleEquator": dict(
        summary="Reduce a total magnetic intensity grid to the pole (RTP) or the equator (RTE): the "
                "anomalies are transformed to look as if measured where the field is vertical (pole) "
                "or horizontal (equator), so they sit more directly above their sources.",
        usage=["RTP suits magnetic latitudes above about 20 degrees; at lower latitudes use "
               "Pole (low latitude) or Equator.",
               "Pole (low latitude) is a stable RTP for low latitudes (Taylor series iterative "
               "method, Hao et al. 2018).",
               "Differential RTP (spatially varying inclination and declination) is in the SGTool window."],
        params={
            "rtp_rte": "Pole: reduction to the pole. Pole (low latitude): stable RTP for low magnetic "
                       "latitudes. Equator: reduction to the equator.",
            "inclination": "Inclination of the geomagnetic field, degrees from horizontal "
                           "(negative in the southern hemisphere).",
            "declination": "Declination of the geomagnetic field, degrees clockwise from north.",
            "iterations": "Low-latitude RTP only: number of iterations. 0 chooses it automatically, as the "
                          "best match to the analytic signal. Fewer iterations give less striping along "
                          "the declination but a less complete reduction to the pole.",
            "stripe_width": "Low-latitude RTP only: suppresses stripes parallel to the declination "
                            "(degrees, 0 = off). Within about this angle of the direction perpendicular "
                            "to the declination the iterations fade to none. 10-20 degrees suits "
                            "inclinations within a few degrees of the equator; use less further from it.",
        }),
    "UpwardContinuation": dict(
        summary="Upward continue a potential-field grid: the field as it would be measured higher up, "
                "which attenuates short wavelengths (shallow sources and noise).",
        usage=["The height is in map units, so the grid should be in a projected (metre-based) "
               "coordinate system."],
        params={"height": "Height to continue upward by, in map units (metres for a projected grid)."}),
    "DownwardContinuation": dict(
        summary="Downward continue a potential-field grid: the field as it would be measured closer to "
                "the sources, which enhances short wavelengths.",
        usage=["Downward continuation amplifies noise strongly: keep the height small, a few cell sizes "
               "at most, and consider smoothing first."],
        params={"height": "Height to continue downward by, in map units (metres for a projected grid)."}),
    "VerticalIntegration": dict(
        summary="Vertical integration: converts magnetic anomalies into gravity-like (pseudo-gravity) "
                "anomalies with the same decay with distance from the source, for comparison or joint "
                "interpretation with gravity. " + MAG_NOTE,
        usage=["Also useful for stitched grids with very different line spacing.",
               "Requires a projected (metre-based) coordinate system."],
        params={
            "max_wn": "Highest wavenumber kept (radians per map unit); higher wavenumbers are clipped to "
                      "it to limit noise. 0 = no limit.",
            "min_wn": "Lowest wavenumber used (radians per map unit), so the integration filter (1/k) "
                      "does not blow up at the longest wavelengths. 0.0001 matches the QGIS plugin.",
        }),
    "BandPassFilter": dict(
        summary="Band-pass filter: keeps wavelengths between the low and high cut-offs, isolating "
                "features of a given size.",
        usage=["Wavelengths are in map units. A transition width of 0 gives a sharp cut-off; increase "
               "it (start with the cut-off value) to reduce ringing."],
        params={
            "low_cut": "Shortest wavelength kept (map units).",
            "high_cut": "Longest wavelength kept (map units).",
            "high_width": "Width of the cosine roll-off at the short-wavelength cut-off (map units, "
                          "0 = sharp).",
            "low_width": "Width of the cosine roll-off at the long-wavelength cut-off (map units, "
                         "0 = sharp).",
        }),
    "HighPassFilter": dict(
        summary="High-pass filter: keeps wavelengths shorter than the cut-off, isolating short-"
                "wavelength (shallow) features.",
        usage=["A transition width of 0 gives a sharp cut-off; increase it to reduce ringing."],
        params={"cutoff": "Cut-off wavelength (map units).",
                "width": "Width of the cosine roll-off (map units, 0 = sharp)."}),
    "LowPassFilter": dict(
        summary="Low-pass filter: keeps wavelengths longer than the cut-off, isolating long-wavelength "
                "(regional or deep) features.",
        usage=["A transition width of 0 gives a sharp cut-off; increase it to reduce ringing."],
        params={"cutoff": "Cut-off wavelength (map units).",
                "width": "Width of the cosine roll-off (map units, 0 = sharp)."}),
    "AutomaticGainControl": dict(
        summary="Automatic gain control (amplitude normalisation): divides the grid by its local "
                "amplitude, highlighting short-wavelength and low-amplitude features.",
        usage=["The local amplitude is a Gaussian-smoothed absolute value of the grid."],
        params={"window_size": "Width of the smoothing window (Gaussian standard deviation, in pixels)."}),
    "RemoveRegional": dict(
        summary="Remove a regional trend: fits a 1st-order (plane) or 2nd-order polynomial surface to "
                "the valid cells by least squares and subtracts it.",
        usage=["Fitted in the space domain, so the FFT buffer is not used."],
        params={"order": "1st Order: a plane. 2nd Order: a quadratic surface.",
                "buffer_size": "Not used by this tool (the fit is in the space domain)."}),
    "DirectionalButterworthFilter": dict(
        summary="Directional Butterworth band-pass filter: keeps a band of wavelengths within a wedge "
                "of directions about an azimuth.",
        usage=["To remove survey line noise, use Remove Line Noise, which estimates the noise this way "
               "and subtracts it."],
        params={
            "low_cut": "Shortest wavelength kept (map units).",
            "high_cut": "Longest wavelength kept (map units).",
            "direction": "Azimuth of the wedge, degrees clockwise from north.",
            "direction_width": "Half-width of the wedge, in degrees.",
            "order": "Butterworth filter order: higher is a sharper cut-off.",
        }),
    "RemoveLineNoise": dict(
        summary="Remove survey line (acquisition) noise. A zero-centred noise estimate is made from "
                "wavelengths between 2 x the smallest and 10 x the largest line spacing, within a "
                "wedge of directions about the line azimuth, scaled, and subtracted from the grid.",
        usage=["The same calculation as Remove line noise in the QGIS plugin."],
        params={
            "azimuth": "Azimuth of the survey lines (degrees clockwise from north): the noise varies "
                       "across them.",
            "line_spacing_min": "Smallest line spacing in the survey (map units). The noise band starts at "
                                "2 x this spacing.",
            "line_spacing_max": "Largest line spacing (map units). Leave blank if the spacing is constant. "
                                "The noise band ends at 10 x this spacing.",
            "scale": "Multiplier applied to the noise estimate before it is subtracted.",
            "wedge": "Half-width of the directional wedge in degrees: noise within this angle of the "
                     "azimuth is kept in the estimate, with a cosine taper. Wider keeps more of the "
                     "noise (and more of the geology).",
        }),
    "ComputeDerivative": dict(
        summary="Derivative of a potential-field grid in x, y or z (vertical), to a chosen power. "
                "Highlights near-surface, short-wavelength features.",
        usage=[MAG_NOTE],
        params={"direction": "x (east), y (north) or z (vertical) derivative.",
                "order": "Power of the derivative: 1 for the first derivative, 2 for the second; "
                         "fractional powers are allowed."}),
    "TotalHorizontalGradient": dict(
        summary="Total horizontal gradient: the magnitude of the horizontal gradient, which peaks over "
                "near-vertical contacts.",
        usage=[MAG_NOTE], params={}),
    "TiltAngle": dict(
        summary="Tilt angle (tilt derivative): arctan of the vertical derivative over the total "
                "horizontal gradient, in degrees (-90 to 90). Enhances edges and shallow sources "
                "independently of their amplitude.",
        usage=[MAG_NOTE, "Tends to over-connect structural features."], params={}),
    "AnalyticSignal": dict(
        summary="Analytic signal amplitude: combines the horizontal and vertical derivatives "
                "(sqrt(dx^2 + dy^2 + dz^2)) to highlight anomaly edges and amplitude, independently "
                "of direction. Always positive.",
        usage=["Largely independent of the magnetisation direction, so it can be used on grids that "
               "have not been reduced to the pole."], params={}),
    "MeanFilter": dict(
        summary="Mean filter: the mean of the values around each cell (a moving average), which "
                "smooths the grid. No-data cells are left out of each mean.",
        usage=[], params={"kernel_size": "Width of the square window in pixels (odd)."}),
    "MedianFilter": dict(
        summary="Median filter: the median of the values around each cell, which removes spikes and "
                "high-frequency noise while keeping edges. No-data cells are left out.",
        usage=[], params={"kernel_size": "Width of the square window in pixels (odd)."}),
    "GaussianFilter": dict(
        summary="Gaussian filter: smooths the grid with a Gaussian kernel. No-data cells are left out.",
        usage=[], params={"sigma": "Standard deviation of the Gaussian, in pixels: larger smooths more."}),
    "DirectionalConvFilter": dict(
        summary="Directional enhancement: highlights short-wavelength features in a chosen compass "
                "direction.",
        usage=["A kernel size of 3 matches the QGIS plugin."],
        params={"direction": "Direction of enhancement.",
                "kernel_size": "Width of the kernel in pixels (odd)."}),
    "SunShading": dict(
        summary="Sun shading (hillshade) of a grid: a shaded-relief image lit from the given azimuth "
                "and altitude, for visualising topography and enhancing features.",
        usage=["GRASS-style uses the GRASS r.relief algorithm for softer shading.",
               "For a grid in degrees (geographic), the shading treats the cell size as degrees; the "
               "QGIS plugin converts it to metres, so prefer a projected grid."],
        params={
            "method": "Simple: standard hillshade. GRASS-style: the GRASS r.relief algorithm (softer).",
            "altitude": "Sun altitude above the horizon in degrees (90 = directly overhead).",
            "azimuth": "Sun azimuth in degrees clockwise from north (315 = north-west).",
            "z_scale": "GRASS-style only: vertical exaggeration applied to the grid values.",
        }),
    "SpatialStatistics": dict(
        summary="Windowed statistics: the variance, standard deviation, skewness, kurtosis, minimum, "
                "maximum, mean or median of the values in a circular window around each cell.",
        usage=["The window is a circle whose diameter is the window size; no-data cells are left out, "
               "and a window needs at least 3 valid values.",
               "Skewness and kurtosis are the biased estimates of scipy.stats (kurtosis is the excess "
               "kurtosis, 0 for a normal distribution); a near-constant window gives no-data."],
        params={"stat_type": "Statistic to calculate in each window.",
                "window_size": "Diameter of the circular window in pixels (odd; an even value is "
                               "rounded up)."}),
    "MRVBFTool": dict(
        summary="Multiresolution Valley Bottom Flatness (MRVBF) and Ridge Top Flatness (MRRTF) after "
                "Gallant & Dowling (2003), from a DEM, with slope and a 3-band composite.",
        usage=["Needs a projected DEM with cell sizes in metres.",
               "Writes the MRVBF grid to the output, and beside it <output>_MRRTF, <output>_slope and a "
               "3-band <output>_RGB composite (R = MRRTF, G = MRVBF, B = slope)."],
        params={
            "in_raster": "Digital elevation model in a projected coordinate system (metres).",
            "out_raster": "Output MRVBF raster; the MRRTF, slope and RGB grids are written beside it.",
            "t_slope": "Initial threshold for slope (percent) at the finest resolution.",
            "t_pctl_v": "Threshold for the lowness percentile (valley bottom flatness).",
            "t_pctl_r": "Threshold for the upness percentile (ridge top flatness).",
            "p_slope": "Shape parameter for the slope transformation.",
            "p_pctl": "Shape parameter for the elevation percentile transformation.",
            "max_res": "Coarsest resolution used, as a percentage of the grid diagonal.",
        }),
    "PCAAnalysis": dict(
        summary="Principal component analysis of a multi-band raster: reduces the bands to the "
                "components that hold most of the variance, for identifying patterns in large datasets.",
        usage=["Only works on multi-band rasters. Needs scikit-learn (installed by Set Up Python "
               "Environment)."],
        params={"in_raster": "Multi-band raster to analyse.",
                "out_raster": "Multi-band GeoTIFF of the principal components.",
                "n_components": "Number of components to keep (0 = all)."}),
    "ICAAnalysis": dict(
        summary="Independent component analysis of a multi-band raster: separates the bands into "
                "additive, independent components, for separating mixed signals.",
        usage=["Only works on multi-band rasters. Needs scikit-learn (installed by Set Up Python "
               "Environment)."],
        params={"in_raster": "Multi-band raster to analyse.",
                "out_raster": "Multi-band GeoTIFF of the independent components.",
                "n_components": "Number of components to keep (0 = all)."}),
    "EulerDeconvolution": dict(
        summary="Euler deconvolution: estimates the position and depth of gravity or magnetic sources "
                "from the grid's gradients in moving windows, keeping the most reliable solutions "
                "(Reliable Euler Deconvolution, Melo & Barbosa 2020).",
        usage=["The same solver as the QGIS plugin.",
               "The output CSV has one row per kept solution: X, Y, Depth, B (base level) and StdDfDz "
               "(the standard deviation of the vertical derivative in the window, used to rank them)."],
        params={
            "out_csv": "CSV file of the kept solutions.",
            "si": "Structural index: 0 contact, 1 dyke or sill edge, 2 pipe or horizontal cylinder, "
                  "3 sphere or dipole.",
            "window_size": "Width of the moving window in pixels (odd).",
            "threshold": "Fraction of the solutions to keep (0-1), those with the largest standard "
                         "deviation of the vertical derivative: 0.1 keeps the best 10%.",
        }),
    "ApparentSusceptibility": dict(
        summary="Apparent susceptibility (SI) from a magnetic anomaly grid (nT, IGRF removed) that has "
                "been reduced to the pole (or approximately, the equator). The ground is modelled as a "
                "layer of vertical prisms magnetised by induction, and the layer's response is divided "
                "out in the Fourier domain.",
        usage=["The grid's absolute level is not used: the result is a variation about the background "
               "susceptibility.",
               "Requires a projected (metre-based) coordinate system."],
        params={
            "in_raster": "Magnetic anomaly grid in nT, reduced to the pole (or the equator).",
            "field": "Strength of the geomagnetic field (nT). Leave blank to use the IGRF at the centre "
                     "of the grid in the survey year.",
            "year": "Survey year (decimal), for the IGRF field strength.",
            "depth": "Depth to the top of the layer (m).",
            "thickness": "Thickness of the layer (m); 0 = infinite thickness (the usual choice).",
            "background": "Background susceptibility (SI) that the result varies about.",
            "high_cut": "Short wavelengths are removed below this wavelength (m), so the inverse filter "
                        "does not amplify noise. 0 = 4 x the cell size, -1 = no cut.",
        }),
    "ApparentDensity": dict(
        summary="Apparent density (g/cc) from a Bouguer gravity grid. The ground is modelled as a layer "
                "of vertical prisms of finite thickness, and the layer's response is divided out in the "
                "Fourier domain.",
        usage=["A Bouguer grid's absolute level depends on its datum and the regional field, so it is "
               "not used: the result's mean is the background density.",
               "Requires a projected (metre-based) coordinate system."],
        params={
            "in_raster": "Bouguer gravity grid in mGal or um/s2.",
            "units": "Units of the gravity grid (1 mGal = 10 um/s2 = 10 gu).",
            "depth": "Depth to the top of the layer (m).",
            "thickness": "Thickness of the layer (m); must be greater than 0.",
            "background": "Background density (g/cc): the mean of the result.",
            "high_cut": "Short wavelengths are removed below this wavelength (m), so the inverse filter "
                        "does not amplify noise. 0 = 4 x the cell size, -1 = no cut.",
        }),
    "IDWGridding": dict(
        summary="Grid point data with inverse distance weighting (ArcGIS Spatial Analyst).",
        usage=["Needs the Spatial Analyst extension. Multilevel B-Spline Gridding does not."],
        params={"in_points": "Point layer to grid.", "z_field": "Field holding the values to grid.",
                "cell_size": "Cell size of the output grid (map units).",
                "power": "IDW power: higher gives nearer points more weight.",
                "out_raster": "Output raster."}),
    "SplineGridding": dict(
        summary="Grid point data with a regularised or tension spline (ArcGIS Spatial Analyst).",
        usage=["Needs the Spatial Analyst extension. Multilevel B-Spline Gridding does not."],
        params={"in_points": "Point layer to grid.", "z_field": "Field holding the values to grid.",
                "cell_size": "Cell size of the output grid (map units).",
                "spline_type": "REGULARIZED: a smooth surface. TENSION: a stiffer surface, closer to "
                               "the data.",
                "out_raster": "Output raster."}),
    "MultilevelBSplineGridding": dict(
        summary="Grid point data with the multilevel B-spline method (a translation of SAGA's "
                "grid_spline tool, as in the SGTool QGIS plugin).",
        usage=["Does not need the Spatial Analyst extension."],
        params={
            "in_points": "Point layer to grid.",
            "z_field": "Field holding the values to grid.",
            "cell_size": "Cell size of the output grid (map units).",
            "epsilon": "Threshold error: refinement stops once every point's residual is below it.",
            "max_levels": "Maximum number of refinement levels.",
            "auto": "Choose the cell size and maximum number of levels from the spacing of the points "
                    "(the two settings above are then ignored).",
            "ignore_below": "Points with a value below this are left out (for example no-data values "
                            "such as -99999). Blank = use all points.",
            "out_raster": "Output raster.",
        }),
    "BSDWorms": dict(
        summary="Worms: multi-scale edges of a potential-field grid. The grid is upward continued to a "
                "series of levels and the edges (maxima of the horizontal gradient) are found at each, "
                "using Frank Horowitz's bsdwormer code.",
        usage=["The source depth of an edge is roughly half the continuation height of the level it "
               "appears at (Jacobsen 1987).",
               "Writes <grid>_worms.csv (points with the level as z) and the padded grid to the output "
               "folder, and optionally <grid>_worms.shp: 3D polylines (PolylineZ) with the level as Z.",
               "The grid needs a coordinate system with an EPSG code."],
        params={
            "out_dir": "Folder for the worms CSV, the padded grid and the optional shapefile.",
            "num_levels": "Number of continuation levels.",
            "base_level": "Height of the lowest level (map units).",
            "delta": "Height between levels (map units).",
            "save_shp": "Also write the worms as a 3D polyline shapefile (can be slow; best with the "
                        "base level at 2000 m or more).",
        }),
    "WTMMAnalysis": dict(
        summary="Wavelet transform modulus maxima (WTMM): multi-scale edge positions of a grid from the "
                "maxima of its wavelet transform.",
        usage=[],
        params={"out_raster": "Output raster of the wavelet transform modulus.",
                "num_scales": "Number of wavelet scales.",
                "threshold_rel": "Maxima below this fraction of the largest are ignored.",
                "min_distance": "Smallest distance between maxima (pixels)."}),
    "ThresholdToNaN": dict(
        summary="Set background values to no-data: cells above, below or between thresholds become "
                "no-data, to mask outliers or noisy borders.",
        usage=[],
        params={"condition": "above: values above the upper threshold. below: values below the lower "
                             "threshold. between: values between the two.",
                "above_value": "Upper threshold value.",
                "below_value": "Lower threshold value."}),
    "CreateBoundaryPolygon": dict(
        summary="Outline the area with data: traces the boundary between data and no-data cells and "
                "saves it as polygons.",
        usage=[], params={"out_shapefile": "Output polygon shapefile."}),
    "ClipAllLayersToPolygon": dict(
        summary="Clip every layer of a map in the open project to a polygon layer (or its selected "
                "polygons), and save a copy of the project that uses the clipped layers.",
        usage=["Clipped vector layers go into one file geodatabase, clipped grids are GeoTIFFs, and "
               "each layer's symbology is saved as a .lyrx file, named <layer>_clip, with a copy of the "
               "project (.aprx) beside them.",
               "Run it from inside ArcGIS Pro with the project open. Web layers, annotation and tables "
               "are left as they are; the tool stops when you press Cancel."],
        params={"in_polygon": "Polygon layer to clip to (only its selected polygons if some are selected).",
                "out_folder": "Folder for the clipped data and the project copy.",
                "map_name": "Map whose layers are clipped (blank = the active map)."}),
    "NormaliseGrids": dict(
        summary="Normalise a folder of GeoTIFFs so they can be merged: removes a 1st- or 2nd-order "
                "regional from each, scales them to the standard deviation of the alphabetically first "
                "grid, and removes the mean.",
        usage=["Assumes the grids have the same processing level, flight heights normalised by "
               "continuation, and the same coordinate system."],
        params={"in_dir": "Folder of GeoTIFFs to normalise.",
                "out_dir": "Folder for the normalised GeoTIFFs.",
                "order": "1st Order: remove a plane. 2nd Order: remove a quadratic surface (slower)."}),
    "ConvertRGBToGrayscale": dict(
        summary="Convert a colour (RGB) image of a grid back to values: each pixel takes the value of "
                "the nearest colour in a colour sequence, spread evenly between a minimum and a maximum.",
        usage=["Do not use on images with shading applied."],
        params={"in_raster": "RGB raster (at least 3 bands).",
                "out_raster": "Output single-band raster of values.",
                "colour_list": "Colours from the lowest to the highest value, one per line, as CSS colour "
                               "names or R,G,B triplets.",
                "min_value": "Value given to the first colour.",
                "max_value": "Value given to the last colour."}),
}


def _load_toolbox():
    """The Toolbox from SGTool.pyt, with arcpy replaced by a stand-in that only records parameters."""
    class Param:
        def __init__(self, **kw):
            self.__dict__.update(kw)
            self.filter = types.SimpleNamespace(type=None, list=[])
            self.value = None
            self.parameterDependencies = []

    arcpy = types.ModuleType("arcpy")
    arcpy.Parameter = Param
    arcpy.env = types.SimpleNamespace(isCancelled=False)
    arcpy.ExecuteError = Exception
    arcpy.__getattr__ = lambda name: (lambda *a, **k: None)
    saved = {n: sys.modules.get(n) for n in ("arcpy", "arcpy.mp", "arcpy.da", "arcpy.sa")}
    sys.modules["arcpy"] = arcpy
    for n in ("arcpy.mp", "arcpy.da", "arcpy.sa"):
        sys.modules[n] = types.ModuleType(n)
    sys.path[:0] = [ARCGIS, os.path.dirname(ARCGIS)]
    try:
        ns = {"__file__": PYT, "__name__": "sgtool_pyt"}
        exec(compile(open(PYT, encoding="utf-8").read(), PYT, "exec"), ns)
        return ns["Toolbox"]()
    finally:
        for n, m in saved.items():
            if m is None:
                sys.modules.pop(n, None)
            else:
                sys.modules[n] = m


def _rich(*paragraphs):
    """Esri metadata rich text: escaped HTML paragraphs."""
    body = "".join(f"<P><SPAN>{html.escape(p, quote=False)}</SPAN></P>" for p in paragraphs if p)
    return html.escape(f"<DIV STYLE=\"text-align:Left;\"><DIV>{body}</DIV></DIV>", quote=False)


def _param_text(tool, p):
    text = TOOLS.get(tool, {}).get("params", {}).get(p.name)
    if text:
        return text
    if p.name in ("in_raster",):
        return IN_RASTER
    if p.name in ("out_raster",):
        return OUT_RASTER
    if p.name == "buffer_size":
        return BUFFER
    return None


def _esri(stamp):
    return (f"<Esri><CreaDate>{stamp:%Y%m%d}</CreaDate><CreaTime>{stamp:%H%M%S}00</CreaTime>"
            "<ArcGISFormat>1.0</ArcGISFormat><SyncOnce>TRUE</SyncOnce>"
            f"<ModDate>{stamp:%Y%m%d}</ModDate><ModTime>{stamp:%H%M%S}00</ModTime></Esri>")


def write(toolbox, out_dir=ARCGIS):
    stamp = datetime.now()
    alias = toolbox.alias
    written, missing = [], []
    for cls in toolbox.tools:
        tool = cls()
        name = cls.__name__
        info = TOOLS.get(name)
        if info is None:
            missing.append(name)
            info = dict(summary=getattr(tool, "description", ""), usage=[], params={})
        params = tool.getParameterInfo() or []
        xml_params = []
        for p in params:
            text = _param_text(name, p) or p.displayName
            ptype = getattr(p, "parameterType", "Required")
            xml_params.append(
                f'<param name="{html.escape(p.name)}" displayname="{html.escape(p.displayName)}" '
                f'type="{ptype}" direction="{getattr(p, "direction", "Input")}" '
                f'datatype="{html.escape(str(getattr(p, "datatype", "")))}" expression="{html.escape(p.name)}">'
                f"<dialogReference>{_rich(text)}</dialogReference>"
                f"<pythonReference>{_rich(text)}</pythonReference></param>")
        usage = info.get("usage") or []
        xml = (
            '<?xml version="1.0" encoding="UTF-8"?>\n'
            f'<metadata xml:lang="en">{_esri(stamp)}'
            f'<tool name="{name}" displayname="{html.escape(tool.label)}" toolboxalias="{alias}" xmlns="">'
            "<arcToolboxHelpPath>c:\\program files\\arcgis\\pro\\Resources\\Help\\gp</arcToolboxHelpPath>"
            f"<summary>{_rich(info['summary'])}</summary>"
            + (f"<usage>{_rich(*usage)}</usage>" if usage else "")
            + f"<parameters>{''.join(xml_params)}</parameters></tool>"
            f"<dataIdInfo><idCitation><resTitle>{html.escape(tool.label)}</resTitle></idCitation>"
            f"<idAbs>{_rich(info['summary'])}</idAbs>"
            "<searchKeys><keyword>SGTool</keyword><keyword>geophysics</keyword></searchKeys>"
            f"<idPurp>{html.escape(info['summary'])}</idPurp></dataIdInfo>"
            "<distInfo><distributor><distorFormat><formatName>ArcToolbox Tool</formatName>"
            "</distorFormat></distributor></distInfo>"
            '<mdHrLv><ScopeCd value="005"/></mdHrLv></metadata>\n')
        path = os.path.join(out_dir, f"SGTool.{name}.pyt.xml")
        with open(path, "w", encoding="utf-8", newline="\n") as f:
            f.write(xml)
        written.append(path)

    tb_xml = (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        f'<metadata xml:lang="en">{_esri(stamp)}'
        f'<toolbox name="SGTool" alias="{alias}">'
        "<arcToolboxHelpPath>c:\\program files\\arcgis\\pro\\Resources\\Help\\gp</arcToolboxHelpPath>"
        "<toolsets/></toolbox>"
        f"<dataIdInfo><idCitation><resTitle>{html.escape(toolbox.label)}</resTitle></idCitation>"
        f"<idAbs>{_rich(TOOLBOX_SUMMARY)}</idAbs><idPurp>{html.escape(TOOLBOX_SUMMARY)}</idPurp>"
        "<searchKeys><keyword>SGTool</keyword><keyword>geophysics</keyword></searchKeys></dataIdInfo>"
        "<distInfo><distributor><distorFormat><formatName>ArcToolbox Toolbox</formatName>"
        "</distorFormat></distributor></distInfo></metadata>\n")
    with open(os.path.join(out_dir, "SGTool.pyt.xml"), "w", encoding="utf-8", newline="\n") as f:
        f.write(tb_xml)
    return written, missing


if __name__ == "__main__":
    tb = _load_toolbox()
    files, missing = write(tb)
    print(f"wrote {len(files)} tool metadata files and SGTool.pyt.xml")
    if missing:
        print("no help text yet for:", ", ".join(missing))
    unused = sorted(set(TOOLS) - {c.__name__ for c in tb.tools})
    if unused:
        print("help text for tools no longer in the toolbox:", ", ".join(unused))
