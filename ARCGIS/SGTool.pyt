# -*- coding: utf-8 -*-
"""
SGTool.pyt  –  ArcGIS Pro Python Toolbox
Full-feature geophysical processing toolbox mirroring the SGTool QGIS plugin.

Toolsets (categories) match the QGIS dockwidget tabs:
  • FFT Filters
  • Gradient Filters
  • Convolution + Statistics
  • Grid + Wavelets
  • Utils
"""

import arcpy
import os
import sys
import numpy as np

# ---------------------------------------------------------------------------
# Bootstrap path so calcs/ package is importable
# ---------------------------------------------------------------------------
_HERE    = os.path.dirname(os.path.abspath(__file__))
_SGTOOL  = os.path.dirname(_HERE)
for _p in (_SGTOOL, os.path.join(_SGTOOL, "calcs"),
           os.path.join(_SGTOOL, "calcs", "worms")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

# Import utility helpers from sibling module
sys.path.insert(0, _HERE)
from arcgis_utils import (
    raster_to_numpy, numpy_to_raster, run_raster_tool, make_processor, NODATA
)


# ============================================================================
# Toolbox class
# ============================================================================
class Toolbox:
    def __init__(self):
        self.label   = "SGTool – Geophysical Processing"
        self.alias   = "sgtool"
        self.tools   = [
            # ---- Launch GUI ----
            LaunchGUI,
            LoadSGToolOutputs,
            # ---- FFT Filters ----
            ReductionToPoleEquator,
            UpwardContinuation,
            DownwardContinuation,
            VerticalIntegration,
            BandPassFilter,
            HighPassFilter,
            LowPassFilter,
            AutomaticGainControl,
            RemoveRegional,
            DirectionalButterworthFilter,
            # ---- Gradient Filters ----
            ComputeDerivative,
            TotalHorizontalGradient,
            TiltAngle,
            AnalyticSignal,
            # ---- Convolution + Statistics ----
            MeanFilter,
            MedianFilter,
            GaussianFilter,
            DirectionalConvFilter,
            SunShading,
            SpatialStatistics,
            DTMCurvatureClassifier,
            PCAAnalysis,
            ICAAnalysis,
            EulerDeconvolution,
            # ---- Grid + Wavelets ----
            IDWGridding,
            SplineGridding,
            BSDWorms,
            WTMMAnalysis,
            # ---- Utils ----
            ThresholdToNaN,
            CreateBoundaryPolygon,
            NormaliseGrids,
            ConvertRGBToGrayscale,
        ]


# ============================================================================
# Base helper
# ============================================================================
class _RasterTool:
    """Base class providing common raster-in / raster-out execute pattern."""

    canRunInBackground = False

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        pass

    def updateMessages(self, parameters):
        pass

    def _read(self, input_path):
        """Return (arr, nodata, ll_x, ll_y, cx, cy, sr)."""
        return raster_to_numpy(str(input_path))

    def _write(self, arr, output_path, ll_x, ll_y, cx, cy, sr, nodata=NODATA):
        numpy_to_raster(arr, str(output_path), ll_x, ll_y, cx, cy, sr, nodata)

    def _proc(self, cx, cy, buf):
        return make_processor(cx, cy, buf)

    @staticmethod
    def _param(name, label, dtype="GPRasterLayer",
                direction="Input", required=True, default=None):
        p = arcpy.Parameter(
            displayName=label,
            name=name,
            datatype=dtype,
            parameterType="Required" if required else "Optional",
            direction=direction,
        )
        if default is not None:
            p.value = default
        return p

    @staticmethod
    def _out(name, label, dtype="DERasterDataset"):
        return arcpy.Parameter(
            displayName=label,
            name=name,
            datatype=dtype,
            parameterType="Required",
            direction="Output",
        )


# ============================================================================
# ---- LAUNCH GUI ---------------------------------------------------------
# ============================================================================

class LaunchGUI:
    """Launch the SGTool full-featured tkinter GUI as a separate process."""
    category           = "Launch"
    label              = "Launch SGTool GUI"
    description        = ("Open the SGTool multi-tab graphical interface. "
                          "Runs as a separate window alongside ArcGIS Pro.")
    canRunInBackground = False

    def getParameterInfo(self):
        return []

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        pass

    def updateMessages(self, parameters):
        pass

    def execute(self, parameters, messages):
        import subprocess, json, tempfile, time

        gui_path = os.path.join(_HERE, "SGTool_GUI.py")
        if not os.path.exists(gui_path):
            arcpy.AddError(f"SGTool_GUI.py not found at: {gui_path}")
            return

        # sys.executable in ArcGIS Pro is ArcGISPro.exe, not python.exe.
        # Resolve the real interpreter from the active conda environment prefix.
        python_exe = os.path.join(sys.exec_prefix, "python.exe")
        if not os.path.exists(python_exe):
            python_exe = os.path.normpath(
                os.path.join(os.path.dirname(os.__file__), "..", "python.exe"))
        if not os.path.exists(python_exe):
            arcpy.AddError(
                f"Cannot locate python.exe. Tried: {python_exe}\n"
                "Run manually from the ArcGIS Pro Python window:\n"
                f"  import subprocess; subprocess.Popen([r'{python_exe}', r'{gui_path}'])"
            )
            return

        # Collect current raster layers and .aprx path (runs on GP thread — always works)
        info = {"aprx": "", "layers": {}}
        try:
            proj = arcpy.mp.ArcGISProject("CURRENT")
            info["aprx"] = proj.filePath or ""
            for m in proj.listMaps():
                for lyr in m.listLayers():
                    if lyr.isRasterLayer:
                        display = f"{lyr.name}  [{m.name}]"
                        info["layers"][display] = lyr.dataSource
        except Exception:
            pass

        info_file    = os.path.join(tempfile.gettempdir(), "sgtool_info.json")
        pending_file = os.path.join(tempfile.gettempdir(), "sgtool_pending.json")
        refresh_req  = os.path.join(tempfile.gettempdir(), "sgtool_refresh_req.json")

        with open(info_file, "w") as f:
            json.dump(info, f)
        with open(pending_file, "w") as f:
            json.dump([], f)

        proc = subprocess.Popen([python_exe, gui_path, info_file],
                                creationflags=subprocess.CREATE_NEW_CONSOLE)
        n = len(info["layers"])
        messages.addMessage(
            f"SGTool GUI launched via {python_exe}  "
            f"({n} raster layer(s) passed to dropdown)\n"
            "Outputs will be loaded into the map automatically as they are created.\n"
            "This tool will appear as 'running' until the SGTool GUI window is closed."
        )

        # Poll on the GP thread — arcpy.mp works here, no COM issues.
        seen  = set()
        cycle = 0

        while proc.poll() is None:
            time.sleep(2)
            cycle += 1

            # --- Layer refresh (on GUI request or every ~60 s) ---
            do_refresh = os.path.exists(refresh_req) or (cycle % 30 == 0)
            if do_refresh:
                try:
                    if os.path.exists(refresh_req):
                        os.remove(refresh_req)
                    proj   = arcpy.mp.ArcGISProject("CURRENT")
                    layers = {}
                    for m_obj in proj.listMaps():
                        for lyr in m_obj.listLayers():
                            if lyr.isRasterLayer:
                                display = f"{lyr.name}  [{m_obj.name}]"
                                layers[display] = lyr.dataSource
                    if layers:
                        try:
                            with open(info_file) as f:
                                cur = json.load(f)
                        except Exception:
                            cur = {}
                        cur["layers"] = layers
                        with open(info_file, "w") as f:
                            json.dump(cur, f)
                except Exception:
                    pass

            # --- Auto-load new output rasters ---
            if os.path.exists(pending_file):
                try:
                    with open(pending_file) as f:
                        paths = json.load(f)
                    new_paths = [p for p in paths
                                 if p not in seen and os.path.exists(p)]
                    if new_paths:
                        proj = arcpy.mp.ArcGISProject("CURRENT")
                        maps = proj.listMaps()
                        if maps:
                            m = maps[0]
                            for p in new_paths:
                                try:
                                    m.addDataFromPath(p)
                                    seen.add(p)
                                    messages.addMessage(
                                        f"Added to map: {os.path.basename(p)}")
                                except Exception:
                                    pass
                            arcpy.RefreshActiveView()
                            arcpy.RefreshTOC()
                except Exception:
                    pass

        messages.addMessage("SGTool GUI closed.")


# ============================================================================
# ---- ADD GUI OUTPUTS TO MAP ---------------------------------------------
# ============================================================================

class LoadSGToolOutputs:
    """Add rasters written by the SGTool GUI to the current ArcGIS Pro map."""
    category           = "Launch"
    label              = "Add SGTool Outputs to Map"
    description        = ("Loads rasters written by the SGTool GUI into the active map. "
                          "Run this after processing in the GUI to see results in ArcGIS Pro.")
    canRunInBackground = False

    def getParameterInfo(self):
        return []

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        pass

    def updateMessages(self, parameters):
        pass

    def execute(self, parameters, messages):
        import json, tempfile
        pending_file = os.path.join(tempfile.gettempdir(), "sgtool_pending.json")
        if not os.path.exists(pending_file):
            messages.addMessage("No pending SGTool outputs found.")
            return
        try:
            with open(pending_file) as f:
                paths = json.load(f)
        except Exception as e:
            arcpy.AddError(f"Could not read pending outputs file: {e}")
            return
        if not paths:
            messages.addMessage("No new outputs to add.")
            return

        proj = arcpy.mp.ArcGISProject("CURRENT")
        maps = proj.listMaps()
        if not maps:
            arcpy.AddError("No maps found in the current project.")
            return
        m = maps[0]

        added = 0
        for path in paths:
            if os.path.exists(path):
                try:
                    m.addDataFromPath(path)
                    messages.addMessage(f"Added: {os.path.basename(path)}")
                    added += 1
                except Exception as e:
                    messages.addWarning(f"Could not add {path}: {e}")
            else:
                messages.addWarning(f"File not found (skipped): {path}")

        # Clear the queue
        with open(pending_file, "w") as f:
            json.dump([], f)

        messages.addMessage(f"Added {added} layer(s) to map '{m.name}'.")


# ============================================================================
# ---- FFT FILTERS --------------------------------------------------------
# ============================================================================

class ReductionToPoleEquator(_RasterTool):
    category    = "FFT Filters"
    label       = "Reduction to Pole / Equator"
    description = ("Transform a magnetic grid to the pole (RTP) or equator (RTE) "
                   "using the specified geomagnetic field parameters.")

    def getParameterInfo(self):
        in_r  = self._param("in_raster",   "Input Raster")
        out_r = self._out  ("out_raster",   "Output Raster")

        rtp   = arcpy.Parameter("rtp_rte", "Transform Type",
                                "GPString", "Required", "Input")
        rtp.filter.type = "ValueList"
        rtp.filter.list = ["Pole", "Equator"]
        rtp.value = "Pole"

        inc   = self._param("inclination",  "Inclination (°)",  "GPDouble", default=-60.0)
        dec   = self._param("declination",  "Declination (°)",  "GPDouble", default=0.0)
        buf   = self._param("buffer_size",  "FFT Buffer (pixels)", "GPLong", required=False, default=10)
        return [in_r, out_r, rtp, inc, dec, buf]

    def execute(self, parameters, messages):
        in_r  = parameters[0].valueAsText
        out_r = parameters[1].valueAsText
        mode  = parameters[2].valueAsText
        inc   = float(parameters[3].value)
        dec   = float(parameters[4].value)
        buf   = int(parameters[5].value) if parameters[5].value else 10

        arr, nodata, ll_x, ll_y, cx, cy, sr = self._read(in_r)
        proc = self._proc(cx, cy, buf)
        if mode == "Pole":
            result = proc.reduction_to_pole(arr, inc, dec)
        else:
            result = proc.reduction_to_equator(arr, inc, dec)
        self._write(result, out_r, ll_x, ll_y, cx, cy, sr, nodata)
        messages.addMessage(f"RTP/RTE ({mode}) complete → {out_r}")


class UpwardContinuation(_RasterTool):
    category    = "FFT Filters"
    label       = "Upward Continuation"
    description = "Continue a potential-field grid upward by the specified height."

    def getParameterInfo(self):
        return [
            self._param("in_raster",  "Input Raster"),
            self._out ("out_raster",  "Output Raster"),
            self._param("height",     "Continuation Height (map units)", "GPDouble", default=1000.0),
            self._param("buffer_size","FFT Buffer (pixels)", "GPLong", required=False, default=10),
        ]

    def execute(self, parameters, messages):
        in_r   = parameters[0].valueAsText
        out_r  = parameters[1].valueAsText
        height = float(parameters[2].value)
        buf    = int(parameters[3].value) if parameters[3].value else 10

        arr, nodata, ll_x, ll_y, cx, cy, sr = self._read(in_r)
        proc   = self._proc(cx, cy, buf)
        result = proc.upward_continuation(arr, height)
        self._write(result, out_r, ll_x, ll_y, cx, cy, sr, nodata)
        messages.addMessage(f"Upward continuation by {height} units complete → {out_r}")


class DownwardContinuation(_RasterTool):
    category    = "FFT Filters"
    label       = "Downward Continuation"
    description = ("Continue a potential-field grid downward. "
                   "Use with caution – amplifies high-frequency noise.")

    def getParameterInfo(self):
        return [
            self._param("in_raster",  "Input Raster"),
            self._out ("out_raster",  "Output Raster"),
            self._param("height",     "Continuation Height (map units)", "GPDouble", default=500.0),
            self._param("buffer_size","FFT Buffer (pixels)", "GPLong", required=False, default=10),
        ]

    def execute(self, parameters, messages):
        in_r   = parameters[0].valueAsText
        out_r  = parameters[1].valueAsText
        height = float(parameters[2].value)
        buf    = int(parameters[3].value) if parameters[3].value else 10

        arr, nodata, ll_x, ll_y, cx, cy, sr = self._read(in_r)
        proc   = self._proc(cx, cy, buf)
        result = proc.downward_continuation(arr, height)
        self._write(result, out_r, ll_x, ll_y, cx, cy, sr, nodata)
        messages.addMessage(f"Downward continuation by {height} units complete → {out_r}")


class VerticalIntegration(_RasterTool):
    category    = "FFT Filters"
    label       = "Vertical Integration (Pseudo-Gravity)"
    description = ("Integrate the vertical field to produce pseudo-gravity. "
                   "Converts magnetic anomaly to equivalent gravity.")

    def getParameterInfo(self):
        return [
            self._param("in_raster",    "Input Raster"),
            self._out ("out_raster",    "Output Raster"),
            self._param("max_wn",       "Max Wavenumber (0 = auto)", "GPDouble",
                        required=False, default=0.0),
            self._param("min_wn",       "Min Wavenumber (stability floor)",
                        "GPDouble", required=False, default=1e-6),
            self._param("buffer_size",  "FFT Buffer (pixels)", "GPLong",
                        required=False, default=10),
        ]

    def execute(self, parameters, messages):
        in_r   = parameters[0].valueAsText
        out_r  = parameters[1].valueAsText
        max_wn = float(parameters[2].value) if parameters[2].value else None
        if max_wn == 0.0:
            max_wn = None
        min_wn = float(parameters[3].value) if parameters[3].value else 1e-6
        buf    = int(parameters[4].value) if parameters[4].value else 10

        arr, nodata, ll_x, ll_y, cx, cy, sr = self._read(in_r)
        proc   = self._proc(cx, cy, buf)
        result = proc.vertical_integration(arr, max_wavenumber=max_wn,
                                           min_wavenumber=min_wn)
        self._write(result, out_r, ll_x, ll_y, cx, cy, sr, nodata)
        messages.addMessage(f"Vertical integration complete → {out_r}")


class BandPassFilter(_RasterTool):
    category    = "FFT Filters"
    label       = "Band-Pass Filter"
    description = ("Pass wavelengths between low_cut and high_cut, "
                   "with cosine roll-off at both ends.")

    def getParameterInfo(self):
        return [
            self._param("in_raster",   "Input Raster"),
            self._out ("out_raster",   "Output Raster"),
            self._param("low_cut",     "Low-cut Wavelength (map units)",  "GPDouble", default=5000.0),
            self._param("high_cut",    "High-cut Wavelength (map units)", "GPDouble", default=50000.0),
            self._param("high_width",  "High-transition Width (0=step)",  "GPDouble", default=0.0),
            self._param("low_width",   "Low-transition Width  (0=step)",  "GPDouble", default=0.0),
            self._param("buffer_size", "FFT Buffer (pixels)", "GPLong",   required=False, default=10),
        ]

    def execute(self, parameters, messages):
        in_r  = parameters[0].valueAsText
        out_r = parameters[1].valueAsText
        lc    = float(parameters[2].value)
        hc    = float(parameters[3].value)
        hw    = float(parameters[4].value)
        lw    = float(parameters[5].value)
        buf   = int(parameters[6].value) if parameters[6].value else 10

        arr, nodata, ll_x, ll_y, cx, cy, sr = self._read(in_r)
        proc   = self._proc(cx, cy, buf)
        result = proc.band_pass_filter(arr, lc, hc, hw, lw)
        self._write(result, out_r, ll_x, ll_y, cx, cy, sr, nodata)
        messages.addMessage(f"Band-pass filter ({lc}–{hc} m) complete → {out_r}")


class HighPassFilter(_RasterTool):
    category    = "FFT Filters"
    label       = "High-Pass Filter"
    description = "Pass wavelengths shorter than the cutoff wavelength."

    def getParameterInfo(self):
        return [
            self._param("in_raster",   "Input Raster"),
            self._out ("out_raster",   "Output Raster"),
            self._param("cutoff",      "Cutoff Wavelength (map units)", "GPDouble", default=10000.0),
            self._param("width",       "Transition Width (0=step)",     "GPDouble", default=0.0),
            self._param("buffer_size", "FFT Buffer (pixels)", "GPLong", required=False, default=10),
        ]

    def execute(self, parameters, messages):
        in_r   = parameters[0].valueAsText
        out_r  = parameters[1].valueAsText
        cutoff = float(parameters[2].value)
        width  = float(parameters[3].value)
        buf    = int(parameters[4].value) if parameters[4].value else 10

        arr, nodata, ll_x, ll_y, cx, cy, sr = self._read(in_r)
        proc   = self._proc(cx, cy, buf)
        result = proc.high_pass_filter(arr, cutoff, width)
        self._write(result, out_r, ll_x, ll_y, cx, cy, sr, nodata)
        messages.addMessage(f"High-pass filter (cutoff={cutoff} m) complete → {out_r}")


class LowPassFilter(_RasterTool):
    category    = "FFT Filters"
    label       = "Low-Pass Filter"
    description = "Pass wavelengths longer than the cutoff wavelength (regional)."

    def getParameterInfo(self):
        return [
            self._param("in_raster",   "Input Raster"),
            self._out ("out_raster",   "Output Raster"),
            self._param("cutoff",      "Cutoff Wavelength (map units)", "GPDouble", default=10000.0),
            self._param("width",       "Transition Width (0=step)",     "GPDouble", default=0.0),
            self._param("buffer_size", "FFT Buffer (pixels)", "GPLong", required=False, default=10),
        ]

    def execute(self, parameters, messages):
        in_r   = parameters[0].valueAsText
        out_r  = parameters[1].valueAsText
        cutoff = float(parameters[2].value)
        width  = float(parameters[3].value)
        buf    = int(parameters[4].value) if parameters[4].value else 10

        arr, nodata, ll_x, ll_y, cx, cy, sr = self._read(in_r)
        proc   = self._proc(cx, cy, buf)
        result = proc.low_pass_filter(arr, cutoff, width)
        self._write(result, out_r, ll_x, ll_y, cx, cy, sr, nodata)
        messages.addMessage(f"Low-pass filter (cutoff={cutoff} m) complete → {out_r}")


class AutomaticGainControl(_RasterTool):
    category    = "FFT Filters"
    label       = "Automatic Gain Control (AGC)"
    description = ("Normalise amplitude using a sliding window. "
                   "Equalises near- and far-field sources.")

    def getParameterInfo(self):
        return [
            self._param("in_raster",   "Input Raster"),
            self._out ("out_raster",   "Output Raster"),
            self._param("window_size", "Window Size (pixels, odd)", "GPLong", default=15),
        ]

    def execute(self, parameters, messages):
        in_r   = parameters[0].valueAsText
        out_r  = parameters[1].valueAsText
        win    = int(parameters[2].value)

        arr, nodata, ll_x, ll_y, cx, cy, sr = self._read(in_r)
        proc   = self._proc(cx, cy, 10)
        result = proc.automatic_gain_control(arr, win)
        self._write(result, out_r, ll_x, ll_y, cx, cy, sr, nodata)
        messages.addMessage(f"AGC (window={win}) complete → {out_r}")


class RemoveRegional(_RasterTool):
    category    = "FFT Filters"
    label       = "Remove Regional Trend"
    description = ("Fit and remove a 1st- or 2nd-order polynomial trend "
                   "(frequency-domain method).")

    def getParameterInfo(self):
        in_r  = self._param("in_raster",  "Input Raster")
        out_r = self._out  ("out_raster",  "Output Raster")

        order = arcpy.Parameter("order", "Polynomial Order",
                                "GPString", "Required", "Input")
        order.filter.type = "ValueList"
        order.filter.list = ["1st Order", "2nd Order"]
        order.value = "1st Order"

        buf = self._param("buffer_size", "FFT Buffer (pixels)", "GPLong",
                          required=False, default=10)
        return [in_r, out_r, order, buf]

    def execute(self, parameters, messages):
        in_r  = parameters[0].valueAsText
        out_r = parameters[1].valueAsText
        order = parameters[2].valueAsText
        buf   = int(parameters[3].value) if parameters[3].value else 10

        arr, nodata, ll_x, ll_y, cx, cy, sr = self._read(in_r)
        proc  = self._proc(cx, cy, buf)
        mask  = ~np.isnan(arr)
        if "2nd" in order:
            result = proc.remove_2o_gradient(arr, mask)
        else:
            result = proc.remove_gradient(arr, mask)
        self._write(result, out_r, ll_x, ll_y, cx, cy, sr, nodata)
        messages.addMessage(f"Regional trend removal ({order}) complete → {out_r}")


class DirectionalButterworthFilter(_RasterTool):
    category    = "FFT Filters"
    label       = "Directional Butterworth Band-Pass Filter"
    description = ("Remove or enhance a specific wavelength band within a "
                   "defined azimuth corridor. Useful for suppressing linear "
                   "acquisition noise.")

    def getParameterInfo(self):
        return [
            self._param("in_raster",      "Input Raster"),
            self._out ("out_raster",      "Output Raster"),
            self._param("low_cut",        "Low-cut Wavelength (map units)",   "GPDouble", default=2000.0),
            self._param("high_cut",       "High-cut Wavelength (map units)",  "GPDouble", default=20000.0),
            self._param("direction",      "Azimuth Direction (°, N=0)",       "GPDouble", default=0.0),
            self._param("direction_width","Direction Width (°)",              "GPDouble", default=30.0),
            self._param("order",          "Butterworth Order",                "GPLong",   default=4),
            self._param("buffer_size",    "FFT Buffer (pixels)", "GPLong",    required=False, default=10),
        ]

    def execute(self, parameters, messages):
        in_r  = parameters[0].valueAsText
        out_r = parameters[1].valueAsText
        lc    = float(parameters[2].value)
        hc    = float(parameters[3].value)
        dirn  = float(parameters[4].value)
        dw    = float(parameters[5].value)
        ordr  = int(parameters[6].value)
        buf   = int(parameters[7].value) if parameters[7].value else 10

        arr, nodata, ll_x, ll_y, cx, cy, sr = self._read(in_r)
        proc   = self._proc(cx, cy, buf)
        result = proc.directional_butterworth_band_pass(
            arr, lc, hc, dirn, dw, order=ordr)
        self._write(result, out_r, ll_x, ll_y, cx, cy, sr, nodata)
        messages.addMessage(f"Directional Butterworth ({dirn}°, {lc}–{hc} m) complete → {out_r}")


# ============================================================================
# ---- GRADIENT FILTERS ---------------------------------------------------
# ============================================================================

class ComputeDerivative(_RasterTool):
    category    = "Gradient Filters"
    label       = "Compute Derivative"
    description = "Compute x, y or z (vertical) derivative of a potential-field grid."

    def getParameterInfo(self):
        in_r  = self._param("in_raster",  "Input Raster")
        out_r = self._out  ("out_raster",  "Output Raster")

        dirn = arcpy.Parameter("direction", "Derivative Direction",
                               "GPString", "Required", "Input")
        dirn.filter.type = "ValueList"
        dirn.filter.list = ["x", "y", "z"]
        dirn.value = "z"

        order = self._param("order",       "Order (power)", "GPDouble", default=1.0)
        buf   = self._param("buffer_size", "FFT Buffer (pixels)", "GPLong",
                            required=False, default=10)
        return [in_r, out_r, dirn, order, buf]

    def execute(self, parameters, messages):
        in_r  = parameters[0].valueAsText
        out_r = parameters[1].valueAsText
        dirn  = parameters[2].valueAsText
        order = float(parameters[3].value)
        buf   = int(parameters[4].value) if parameters[4].value else 10

        arr, nodata, ll_x, ll_y, cx, cy, sr = self._read(in_r)
        proc   = self._proc(cx, cy, buf)
        result = proc.compute_derivative(arr, direction=dirn, order=order)
        self._write(result, out_r, ll_x, ll_y, cx, cy, sr, nodata)
        messages.addMessage(f"Derivative ({dirn}, order={order}) complete → {out_r}")


class TotalHorizontalGradient(_RasterTool):
    category    = "Gradient Filters"
    label       = "Total Horizontal Gradient (THG)"
    description = ("Compute the total horizontal gradient magnitude. "
                   "Peaks over vertical contacts.")

    def getParameterInfo(self):
        return [
            self._param("in_raster",   "Input Raster"),
            self._out ("out_raster",   "Output Raster"),
            self._param("buffer_size", "FFT Buffer (pixels)", "GPLong",
                        required=False, default=10),
        ]

    def execute(self, parameters, messages):
        in_r  = parameters[0].valueAsText
        out_r = parameters[1].valueAsText
        buf   = int(parameters[2].value) if parameters[2].value else 10

        arr, nodata, ll_x, ll_y, cx, cy, sr = self._read(in_r)
        proc   = self._proc(cx, cy, buf)
        result = proc.total_hz_grad(arr)
        self._write(result, out_r, ll_x, ll_y, cx, cy, sr, nodata)
        messages.addMessage(f"Total Horizontal Gradient complete → {out_r}")


class TiltAngle(_RasterTool):
    category    = "Gradient Filters"
    label       = "Tilt Angle"
    description = ("Compute the tilt angle (arctan(dz/THG)). "
                   "Result in degrees. Useful for edge detection.")

    def getParameterInfo(self):
        return [
            self._param("in_raster",   "Input Raster"),
            self._out ("out_raster",   "Output Raster"),
            self._param("buffer_size", "FFT Buffer (pixels)", "GPLong",
                        required=False, default=10),
        ]

    def execute(self, parameters, messages):
        in_r  = parameters[0].valueAsText
        out_r = parameters[1].valueAsText
        buf   = int(parameters[2].value) if parameters[2].value else 10

        arr, nodata, ll_x, ll_y, cx, cy, sr = self._read(in_r)
        proc   = self._proc(cx, cy, buf)
        result_rad = proc.tilt_angle(arr)
        result     = np.degrees(result_rad)
        self._write(result, out_r, ll_x, ll_y, cx, cy, sr, nodata)
        messages.addMessage(f"Tilt angle (degrees) complete → {out_r}")


class AnalyticSignal(_RasterTool):
    category    = "Gradient Filters"
    label       = "Analytic Signal"
    description = ("Compute the analytic signal amplitude: "
                   "sqrt(dx² + dy² + dz²). Always positive.")

    def getParameterInfo(self):
        return [
            self._param("in_raster",   "Input Raster"),
            self._out ("out_raster",   "Output Raster"),
            self._param("buffer_size", "FFT Buffer (pixels)", "GPLong",
                        required=False, default=10),
        ]

    def execute(self, parameters, messages):
        in_r  = parameters[0].valueAsText
        out_r = parameters[1].valueAsText
        buf   = int(parameters[2].value) if parameters[2].value else 10

        arr, nodata, ll_x, ll_y, cx, cy, sr = self._read(in_r)
        proc   = self._proc(cx, cy, buf)
        result = proc.analytic_signal(arr)
        self._write(result, out_r, ll_x, ll_y, cx, cy, sr, nodata)
        messages.addMessage(f"Analytic signal complete → {out_r}")


# ============================================================================
# ---- CONVOLUTION + STATISTICS -------------------------------------------
# ============================================================================

class MeanFilter(_RasterTool):
    category    = "Convolution + Statistics"
    label       = "Mean Filter"
    description = "Apply a moving average (boxcar) filter. NaN-aware."

    def getParameterInfo(self):
        return [
            self._param("in_raster",   "Input Raster"),
            self._out ("out_raster",   "Output Raster"),
            self._param("kernel_size", "Kernel Size (pixels, odd)", "GPLong", default=3),
        ]

    def execute(self, parameters, messages):
        in_r  = parameters[0].valueAsText
        out_r = parameters[1].valueAsText
        n     = int(parameters[2].value)

        arr, nodata, ll_x, ll_y, cx, cy, sr = self._read(in_r)
        from calcs.ConvolutionFilter import ConvolutionFilter
        result = ConvolutionFilter(arr).mean_filter(n)
        self._write(result, out_r, ll_x, ll_y, cx, cy, sr, nodata)
        messages.addMessage(f"Mean filter (kernel={n}) complete → {out_r}")


class MedianFilter(_RasterTool):
    category    = "Convolution + Statistics"
    label       = "Median Filter"
    description = "Apply a moving median filter. Robust to spike noise."

    def getParameterInfo(self):
        return [
            self._param("in_raster",   "Input Raster"),
            self._out ("out_raster",   "Output Raster"),
            self._param("kernel_size", "Kernel Size (pixels, odd)", "GPLong", default=3),
        ]

    def execute(self, parameters, messages):
        in_r  = parameters[0].valueAsText
        out_r = parameters[1].valueAsText
        n     = int(parameters[2].value)

        arr, nodata, ll_x, ll_y, cx, cy, sr = self._read(in_r)
        from calcs.ConvolutionFilter import ConvolutionFilter
        result = ConvolutionFilter(arr).median_filter(n)
        self._write(result, out_r, ll_x, ll_y, cx, cy, sr, nodata)
        messages.addMessage(f"Median filter (kernel={n}) complete → {out_r}")


class GaussianFilter(_RasterTool):
    category    = "Convolution + Statistics"
    label       = "Gaussian Filter"
    description = "Smooth the grid with a Gaussian kernel."

    def getParameterInfo(self):
        return [
            self._param("in_raster", "Input Raster"),
            self._out ("out_raster", "Output Raster"),
            self._param("sigma",     "Sigma (pixels)", "GPDouble", default=1.0),
        ]

    def execute(self, parameters, messages):
        in_r  = parameters[0].valueAsText
        out_r = parameters[1].valueAsText
        sigma = float(parameters[2].value)

        arr, nodata, ll_x, ll_y, cx, cy, sr = self._read(in_r)
        from calcs.ConvolutionFilter import ConvolutionFilter
        result = ConvolutionFilter(arr).gaussian_filter(sigma)
        self._write(result, out_r, ll_x, ll_y, cx, cy, sr, nodata)
        messages.addMessage(f"Gaussian filter (σ={sigma}) complete → {out_r}")


class DirectionalConvFilter(_RasterTool):
    category    = "Convolution + Statistics"
    label       = "Directional Enhancement Filter"
    description = "Enhance features in a chosen compass direction."

    def getParameterInfo(self):
        in_r  = self._param("in_raster",   "Input Raster")
        out_r = self._out  ("out_raster",   "Output Raster")

        dirn = arcpy.Parameter("direction", "Enhancement Direction",
                               "GPString", "Required", "Input")
        dirn.filter.type = "ValueList"
        dirn.filter.list = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"]
        dirn.value = "N"

        n = self._param("kernel_size", "Kernel Size (pixels, odd)",
                        "GPLong", default=3)
        return [in_r, out_r, dirn, n]

    def execute(self, parameters, messages):
        in_r  = parameters[0].valueAsText
        out_r = parameters[1].valueAsText
        dirn  = parameters[2].valueAsText
        n     = int(parameters[3].value)

        arr, nodata, ll_x, ll_y, cx, cy, sr = self._read(in_r)
        from calcs.ConvolutionFilter import ConvolutionFilter
        result = ConvolutionFilter(arr).directional_filter(dirn, n)
        self._write(result, out_r, ll_x, ll_y, cx, cy, sr, nodata)
        messages.addMessage(f"Directional filter ({dirn}, n={n}) complete → {out_r}")


class SunShading(_RasterTool):
    category    = "Convolution + Statistics"
    label       = "Sun Shading (Hillshade)"
    description = "Compute an illuminated relief image from a potential-field or DEM grid."

    def getParameterInfo(self):
        in_r  = self._param("in_raster",  "Input Raster")
        out_r = self._out  ("out_raster",  "Output Raster")

        method = arcpy.Parameter("method", "Method",
                                 "GPString", "Required", "Input")
        method.filter.type = "ValueList"
        method.filter.list = ["Simple", "GRASS-style"]
        method.value = "Simple"

        alt   = self._param("altitude",  "Sun Altitude (°)",  "GPDouble", default=45.0)
        azim  = self._param("azimuth",   "Sun Azimuth (°)",   "GPDouble", default=315.0)
        zscl  = self._param("z_scale",   "Z Scale Factor",    "GPDouble", required=False, default=1.0)
        return [in_r, out_r, method, alt, azim, zscl]

    def execute(self, parameters, messages):
        in_r   = parameters[0].valueAsText
        out_r  = parameters[1].valueAsText
        method = parameters[2].valueAsText
        alt    = float(parameters[3].value)
        azim   = float(parameters[4].value)
        zscl   = float(parameters[5].value) if parameters[5].value else 1.0

        arr, nodata, ll_x, ll_y, cx, cy, sr = self._read(in_r)
        from calcs.ConvolutionFilter import ConvolutionFilter
        cf = ConvolutionFilter(arr)
        if method == "Simple":
            result = cf.sun_shading_filter(alt, azim)
        else:
            result = cf.sun_shading_filter_grass(alt, azim, cy, cx, 1, zscl)
        self._write(result, out_r, ll_x, ll_y, cx, cy, sr, nodata)
        messages.addMessage(f"Sun shading ({method}, alt={alt}°, az={azim}°) complete → {out_r}")


class SpatialStatistics(_RasterTool):
    category    = "Convolution + Statistics"
    label       = "Windowed Spatial Statistics"
    description = ("Compute a windowed statistic at each grid cell. "
                   "Options: variance, std, skewness, kurtosis, min, max, mean, median.")

    def getParameterInfo(self):
        in_r  = self._param("in_raster",    "Input Raster")
        out_r = self._out  ("out_raster",    "Output Raster")

        stat = arcpy.Parameter("stat_type", "Statistic Type",
                               "GPString", "Required", "Input")
        stat.filter.type = "ValueList"
        stat.filter.list = ["variance", "std", "skewness", "kurtosis",
                            "min", "max", "mean", "median"]
        stat.value = "variance"

        win = self._param("window_size", "Window Size (pixels, odd)",
                          "GPLong", default=5)
        return [in_r, out_r, stat, win]

    def execute(self, parameters, messages):
        in_r  = parameters[0].valueAsText
        out_r = parameters[1].valueAsText
        stat  = parameters[2].valueAsText
        win   = int(parameters[3].value)

        arr, nodata, ll_x, ll_y, cx, cy, sr = self._read(in_r)
        from calcs.SpatialStats import SpatialStats
        result = SpatialStats(arr).calculate_windowed_stats(win, stat)
        self._write(result, out_r, ll_x, ll_y, cx, cy, sr, nodata)
        messages.addMessage(f"Spatial {stat} (window={win}) complete → {out_r}")


class DTMCurvatureClassifier(_RasterTool):
    category    = "Convolution + Statistics"
    label       = "DTM Curvature Classifier"
    description = ("Classify terrain into concave (-1), flat (0), convex (1) "
                   "and cliff (2) categories.")

    def getParameterInfo(self):
        return [
            self._param("in_raster",   "Input Raster (DEM)"),
            self._out ("out_raster",   "Output Classification Raster"),
            self._param("curv_thresh", "Curvature Threshold",       "GPDouble", default=0.01),
            self._param("slope_thresh","Slope/Cliff Threshold (°)", "GPDouble", default=45.0),
            self._param("window_size", "Window Size (pixels, odd)", "GPLong",   default=3),
            self._param("sigma",       "Smoothing Sigma",           "GPDouble", default=1.0),
        ]

    def execute(self, parameters, messages):
        in_r   = parameters[0].valueAsText
        out_r  = parameters[1].valueAsText
        ct     = float(parameters[2].value)
        st     = float(parameters[3].value)
        win    = int(parameters[4].value)
        sigma  = float(parameters[5].value)

        arr, nodata, ll_x, ll_y, cx, cy, sr = self._read(in_r)
        from calcs.SpatialStats import SpatialStats
        result = SpatialStats(arr).classify_terrain_with_cell_size(
            cx, cy, ct, st, win, sigma)
        self._write(result.astype(np.float32), out_r, ll_x, ll_y, cx, cy, sr, nodata)
        messages.addMessage(f"DTM classification complete → {out_r}")


class PCAAnalysis(_RasterTool):
    category    = "Convolution + Statistics"
    label       = "Principal Component Analysis (PCA)"
    description = ("Perform PCA on a multi-band raster. "
                   "Requires scikit-learn. Output is a multi-band GeoTIFF.")

    canRunInBackground = False

    def getParameterInfo(self):
        in_r  = arcpy.Parameter("in_raster",  "Input Multi-band Raster",
                                "GPRasterLayer", "Required", "Input")
        out_r = arcpy.Parameter("out_raster",  "Output Multi-band Raster",
                                "DERasterDataset", "Required", "Output")
        n     = self._param("n_components", "Number of Components (0=all)",
                            "GPLong", required=False, default=0)
        return [in_r, out_r, n]

    def execute(self, parameters, messages):
        in_r  = parameters[0].valueAsText
        out_r = parameters[1].valueAsText
        n     = int(parameters[2].value) if parameters[2].value else 0

        from calcs.PCAICA import PCAICA
        dummy = np.zeros((2, 2))
        pcaica = PCAICA(dummy)
        pcaica.pca_with_nans(in_r, out_r, n_components=n if n > 0 else None)
        messages.addMessage(f"PCA complete → {out_r}")


class ICAAnalysis(_RasterTool):
    category    = "Convolution + Statistics"
    label       = "Independent Component Analysis (ICA)"
    description = ("Perform ICA on a multi-band raster. "
                   "Requires scikit-learn. Output is a multi-band GeoTIFF.")

    canRunInBackground = False

    def getParameterInfo(self):
        in_r  = arcpy.Parameter("in_raster",  "Input Multi-band Raster",
                                "GPRasterLayer", "Required", "Input")
        out_r = arcpy.Parameter("out_raster",  "Output Multi-band Raster",
                                "DERasterDataset", "Required", "Output")
        n     = self._param("n_components", "Number of Components (0=all)",
                            "GPLong", required=False, default=0)
        return [in_r, out_r, n]

    def execute(self, parameters, messages):
        in_r  = parameters[0].valueAsText
        out_r = parameters[1].valueAsText
        n     = int(parameters[2].value) if parameters[2].value else 0

        from calcs.PCAICA import PCAICA
        dummy = np.zeros((2, 2))
        pcaica = PCAICA(dummy)
        pcaica.ica_with_nans(in_r, out_r, n_components=n if n > 0 else None)
        messages.addMessage(f"ICA complete → {out_r}")


class EulerDeconvolution(_RasterTool):
    category    = "Convolution + Statistics"
    label       = "Euler Deconvolution"
    description = ("Estimate source depths from a potential-field grid using "
                   "the Euler deconvolution method. Outputs source location CSV.")

    def getParameterInfo(self):
        in_r  = self._param("in_raster",   "Input Raster")
        out_csv = arcpy.Parameter("out_csv", "Output CSV",
                                  "DEFile", "Required", "Output")

        si = arcpy.Parameter("si", "Structural Index (SI)",
                             "GPDouble", "Required", "Input")
        si.filter.type = "ValueList"
        si.filter.list = [0.0, 1.0, 2.0, 3.0]
        si.value = 1.0

        win   = self._param("window_size", "Window Size (pixels, odd)", "GPLong",   default=5)
        filt  = self._param("threshold",   "Percent Uncertainty Filter (0–1)",
                            "GPDouble", default=0.1)
        return [in_r, out_csv, si, win, filt]

    def execute(self, parameters, messages):
        in_r   = parameters[0].valueAsText
        out_c  = parameters[1].valueAsText
        SI     = float(parameters[2].value)
        win    = int(parameters[3].value)
        filt   = float(parameters[4].value)

        arr, nodata, ll_x, ll_y, cx, cy, sr = self._read(in_r)
        rows, cols = arr.shape

        from calcs.euler.euler_python_optimised import euler_deconv_optimized
        from calcs.euler.estimates_statistics import statistics_euler

        x_max = ll_x + cols * cx
        y_max = ll_y + rows * cy
        area  = (ll_y, y_max, ll_x, x_max)
        xs    = np.linspace(ll_x, x_max, cols)
        ys    = np.linspace(ll_y, y_max, rows)
        XI, YI = np.meshgrid(xs, ys)
        ZI     = np.zeros_like(arr)

        filled = np.where(np.isnan(arr), 0.0, arr)
        results = euler_deconv_optimized(
            filled, XI, YI, ZI,
            filled.shape, area, SI, win, filt
        )
        if results is not None:
            est_x, est_y, est_z, est_b, stdz = results
            stat_x, stat_y, stat_z = statistics_euler(
                est_x, est_y, est_z, stdz, filt)

            import csv
            with open(out_c, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["X", "Y", "Depth", "BaseLevel", "StdZ"])
                for i in range(len(stat_x)):
                    writer.writerow([stat_x[i], stat_y[i], stat_z[i], 0, 0])
            messages.addMessage(
                f"Euler deconvolution (SI={SI}) → {len(stat_x)} estimates → {out_c}")
        else:
            messages.addWarning("Euler deconvolution produced no results.")


# ============================================================================
# ---- GRID + WAVELETS ----------------------------------------------------
# ============================================================================

class IDWGridding(_RasterTool):
    category    = "Grid + Wavelets"
    label       = "IDW Gridding (from points)"
    description = ("Interpolate scattered point data to a raster using "
                   "Inverse Distance Weighting (ArcGIS Spatial Analyst).")

    def getParameterInfo(self):
        pts   = arcpy.Parameter("in_points",  "Input Points Feature Class",
                                "GPFeatureLayer", "Required", "Input")
        field = arcpy.Parameter("z_field",    "Z-Value Field",
                                "Field", "Required", "Input")
        field.parameterDependencies = ["in_points"]

        cell  = self._param("cell_size",  "Output Cell Size (map units)",
                            "GPDouble", default=100.0)
        power = self._param("power",      "IDW Power", "GPDouble",
                            required=False, default=2.0)
        out_r = self._out  ("out_raster", "Output Raster")
        return [pts, field, cell, power, out_r]

    def execute(self, parameters, messages):
        pts   = parameters[0].valueAsText
        field = parameters[1].valueAsText
        cell  = float(parameters[2].value)
        power = float(parameters[3].value) if parameters[3].value else 2.0
        out_r = parameters[4].valueAsText

        arcpy.CheckOutExtension("Spatial")
        from arcpy.sa import Idw
        result = Idw(pts, field, cell, power)
        result.save(out_r)
        arcpy.CheckInExtension("Spatial")
        messages.addMessage(f"IDW gridding complete → {out_r}")


class SplineGridding(_RasterTool):
    category    = "Grid + Wavelets"
    label       = "Spline Gridding (from points)"
    description = ("Interpolate scattered point data to a raster using "
                   "tension spline (ArcGIS Spatial Analyst).")

    def getParameterInfo(self):
        pts   = arcpy.Parameter("in_points",  "Input Points Feature Class",
                                "GPFeatureLayer", "Required", "Input")
        field = arcpy.Parameter("z_field",    "Z-Value Field",
                                "Field", "Required", "Input")
        field.parameterDependencies = ["in_points"]

        cell  = self._param("cell_size",  "Output Cell Size (map units)",
                            "GPDouble", default=100.0)
        sptype = arcpy.Parameter("spline_type", "Spline Type",
                                 "GPString", "Optional", "Input")
        sptype.filter.type = "ValueList"
        sptype.filter.list = ["REGULARIZED", "TENSION"]
        sptype.value = "REGULARIZED"

        out_r = self._out("out_raster", "Output Raster")
        return [pts, field, cell, sptype, out_r]

    def execute(self, parameters, messages):
        pts    = parameters[0].valueAsText
        field  = parameters[1].valueAsText
        cell   = float(parameters[2].value)
        sptype = parameters[3].valueAsText or "REGULARIZED"
        out_r  = parameters[4].valueAsText

        arcpy.CheckOutExtension("Spatial")
        from arcpy.sa import Spline
        result = Spline(pts, field, cell, sptype)
        result.save(out_r)
        arcpy.CheckInExtension("Spatial")
        messages.addMessage(f"Spline gridding complete → {out_r}")


class BSDWorms(_RasterTool):
    category    = "Grid + Wavelets"
    label       = "BSD Worms (Multi-scale Edge Extraction)"
    description = ("Extract multi-scale maxima of the analytic signal at multiple "
                   "upward-continuation levels. Produces CSV files and optionally "
                   "shapefiles.")

    def getParameterInfo(self):
        in_r  = self._param("in_raster",   "Input Raster")
        out_d = arcpy.Parameter("out_dir",  "Output Directory",
                                "DEFolder", "Required", "Output")
        levels = self._param("num_levels", "Number of Continuation Levels",
                             "GPLong",   default=10)
        base   = self._param("base_level", "Base Continuation Level (map units)",
                             "GPDouble", default=1000.0)
        delta  = self._param("delta",      "Level Increment (map units)",
                             "GPDouble", default=1000.0)
        shp    = arcpy.Parameter("save_shp", "Also Save Shapefiles",
                                 "GPBoolean", "Optional", "Input")
        shp.value = False
        return [in_r, out_d, levels, base, delta, shp]

    def execute(self, parameters, messages):
        in_r   = parameters[0].valueAsText
        out_d  = parameters[1].valueAsText
        levels = int(parameters[2].value)
        base   = float(parameters[3].value)
        delta  = float(parameters[4].value)
        shp    = bool(parameters[5].value)

        arr, nodata, ll_x, ll_y, cx, cy, sr = self._read(in_r)
        proc = self._proc(cx, cy, 10)

        # Get CRS string
        try:
            crs = sr.exportToString() if hasattr(sr, 'exportToString') else str(sr)
        except Exception:
            crs = ""

        os.makedirs(out_d, exist_ok=True)
        layer_name = os.path.splitext(os.path.basename(in_r))[0]
        proc.bsdwormer(arr, layer_name, out_d,
                       levels, base, delta, shp, crs)
        messages.addMessage(f"BSD Worms ({levels} levels) complete → {out_d}")


class WTMMAnalysis(_RasterTool):
    category    = "Grid + Wavelets"
    label       = "Wavelet Transform Modulus Maxima (WTMM)"
    description = ("Detect multi-scale edge positions using the WTMM method "
                   "on a 2-D grid.")

    def getParameterInfo(self):
        in_r  = self._param("in_raster",     "Input Raster")
        out_r = self._out  ("out_raster",     "Output Modulus Raster")
        scales= self._param("num_scales",     "Number of Scales",    "GPLong",   default=5)
        thresh= self._param("threshold_rel",  "Relative Threshold",  "GPDouble", default=0.2)
        mindist=self._param("min_distance",   "Min Distance (pixels)","GPLong",  default=5)
        return [in_r, out_r, scales, thresh, mindist]

    def execute(self, parameters, messages):
        in_r   = parameters[0].valueAsText
        out_r  = parameters[1].valueAsText
        scales = int(parameters[2].value)
        thresh = float(parameters[3].value)
        mdist  = int(parameters[4].value)

        arr, nodata, ll_x, ll_y, cx, cy, sr = self._read(in_r)
        from calcs.WTMM import WTMM
        wtmm   = WTMM()
        res    = wtmm.wtmm_2d(
            arr, num_scales=scales,
            threshold_rel=thresh, min_distance=mdist)
        # wtmm_2d returns a dict; extract the modulus maxima grid
        if isinstance(res, dict):
            modulus = res.get("modulus", arr)
        else:
            modulus = arr if res is None else res
        self._write(modulus, out_r, ll_x, ll_y, cx, cy, sr, nodata)
        messages.addMessage(f"WTMM analysis complete → {out_r}")


# ============================================================================
# ---- UTILS --------------------------------------------------------------
# ============================================================================

class ThresholdToNaN(_RasterTool):
    category    = "Utils"
    label       = "Threshold to NaN"
    description = ("Replace values above, below, or between thresholds with NaN "
                   "(output nodata). Useful to mask outliers or noisy borders.")

    def getParameterInfo(self):
        in_r  = self._param("in_raster",  "Input Raster")
        out_r = self._out  ("out_raster",  "Output Raster")

        cond = arcpy.Parameter("condition", "Condition",
                               "GPString", "Required", "Input")
        cond.filter.type = "ValueList"
        cond.filter.list = ["above", "below", "between"]
        cond.value = "above"

        above = self._param("above_value", "Above Threshold Value", "GPDouble", default=1e10)
        below = self._param("below_value", "Below Threshold Value", "GPDouble", default=-1e10)
        return [in_r, out_r, cond, above, below]

    def execute(self, parameters, messages):
        in_r   = parameters[0].valueAsText
        out_r  = parameters[1].valueAsText
        cond   = parameters[2].valueAsText
        above  = float(parameters[3].value)
        below  = float(parameters[4].value)

        arr, nodata, ll_x, ll_y, cx, cy, sr = self._read(in_r)
        from calcs.SG_Util import SG_Util
        result = SG_Util(arr).Threshold2Nan(arr, cond, above, below)
        self._write(result, out_r, ll_x, ll_y, cx, cy, sr, nodata)
        messages.addMessage(f"Threshold to NaN ({cond}) complete → {out_r}")


class CreateBoundaryPolygon(_RasterTool):
    category    = "Utils"
    label       = "Create Data Boundary Polygon"
    description = ("Trace the boundary between data and NoData regions "
                   "and save as a polygon shapefile.")

    def getParameterInfo(self):
        in_r  = self._param("in_raster",  "Input Raster")
        out_shp = arcpy.Parameter("out_shapefile", "Output Shapefile",
                                  "DEShapefile", "Required", "Output")
        return [in_r, out_shp]

    def execute(self, parameters, messages):
        in_r    = parameters[0].valueAsText
        out_shp = parameters[1].valueAsText

        from calcs.SG_Util import SG_Util
        dummy  = np.zeros((2, 2))
        sg     = SG_Util(dummy)
        result = sg.create_data_boundary_lines(in_r, out_shp)
        messages.addMessage(f"Boundary polygon created → {result or out_shp}")


class NormaliseGrids(_RasterTool):
    category    = "Utils"
    label       = "Normalise Grids (Batch)"
    description = ("Normalise all GeoTIFFs in an input folder by fitting a 1st- "
                   "or 2nd-order polynomial trend and dividing by it.")

    def getParameterInfo(self):
        in_d  = arcpy.Parameter("in_dir",   "Input Directory",
                                "DEFolder", "Required", "Input")
        out_d = arcpy.Parameter("out_dir",  "Output Directory",
                                "DEFolder", "Required", "Output")

        order = arcpy.Parameter("order", "Polynomial Order",
                                "GPString", "Required", "Input")
        order.filter.type = "ValueList"
        order.filter.list = ["1st Order", "2nd Order"]
        order.value = "1st Order"

        return [in_d, out_d, order]

    def execute(self, parameters, messages):
        in_d   = parameters[0].valueAsText
        out_d  = parameters[1].valueAsText
        order  = 1 if "1st" in parameters[2].valueAsText else 2

        os.makedirs(out_d, exist_ok=True)
        proc = make_processor(1.0, 1.0, 10)
        proc.normalise_geotiffs(in_d, out_d, order)
        messages.addMessage(f"Grid normalisation complete → {out_d}")


class ConvertRGBToGrayscale(_RasterTool):
    category    = "Utils"
    label       = "Convert RGB Grid to Grayscale"
    description = ("Convert an RGB colour-coded raster to grayscale values using "
                   "a CSS colour list or user-defined RGB triplets.")

    def getParameterInfo(self):
        in_r   = arcpy.Parameter("in_raster",  "Input RGB Raster",
                                 "GPRasterLayer", "Required", "Input")
        out_r  = self._out("out_raster", "Output Grayscale Raster")
        colours = arcpy.Parameter("colour_list",
                                  "Colour List (CSS names or R,G,B triplets, one per line)",
                                  "GPString", "Required", "Input")
        colours.multiValue = False
        min_v  = self._param("min_value", "Output Minimum Value", "GPDouble", default=0.0)
        max_v  = self._param("max_value", "Output Maximum Value", "GPDouble", default=1000.0)
        return [in_r, out_r, colours, min_v, max_v]

    def execute(self, parameters, messages):
        in_r     = parameters[0].valueAsText
        out_r    = parameters[1].valueAsText
        col_str  = parameters[2].valueAsText
        min_v    = float(parameters[3].value)
        max_v    = float(parameters[4].value)

        from osgeo import gdal
        ds     = gdal.Open(in_r, gdal.GA_ReadOnly)
        if ds is None:
            arcpy.AddError(f"Cannot open: {in_r}")
            return

        bands  = [ds.GetRasterBand(i + 1).ReadAsArray().astype(np.float32)
                  for i in range(min(ds.RasterCount, 3))]
        gt     = ds.GetGeoTransform()
        proj   = ds.GetProjection()
        ds     = None

        if len(bands) < 3:
            arcpy.AddError("Input raster must have at least 3 bands (R, G, B).")
            return

        R, G, B = bands[0], bands[1], bands[2]

        # Parse colour list
        lut = []
        for line in col_str.strip().splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                import matplotlib.colors as mcolors
                rgba = mcolors.to_rgba(line)
                lut.append((int(rgba[0] * 255),
                             int(rgba[1] * 255),
                             int(rgba[2] * 255)))
            except (ValueError, KeyError):
                parts = [p.strip() for p in line.split(",")]
                if len(parts) == 3:
                    try:
                        lut.append(tuple(int(p) for p in parts))
                    except ValueError:
                        pass

        if not lut:
            arcpy.AddError("No valid colours parsed from colour list.")
            return

        # Map each pixel to nearest LUT colour and assign normalised value
        n      = len(lut)
        lut_a  = np.array(lut, dtype=np.float32)
        rows, cols = R.shape
        result = np.full((rows, cols), np.nan, dtype=np.float32)

        for idx, (lr, lg, lb) in enumerate(lut):
            dist = np.sqrt((R - lr) ** 2 + (G - lg) ** 2 + (B - lb) ** 2)
            if idx == 0:
                best_dist = dist
                best_idx  = np.zeros_like(R, dtype=np.int32)
            else:
                better = dist < best_dist
                best_dist[better] = dist[better]
                best_idx[better]  = idx

        result = min_v + (best_idx.astype(np.float32) / max(n - 1, 1)) * (max_v - min_v)

        # Write output
        driver = gdal.GetDriverByName("GTiff")
        out_ds = driver.Create(out_r, cols, rows, 1, gdal.GDT_Float32)
        out_ds.SetGeoTransform(gt)
        if proj:
            out_ds.SetProjection(proj)
        band = out_ds.GetRasterBand(1)
        band.WriteArray(result)
        band.SetNoDataValue(NODATA)
        out_ds.FlushCache()
        out_ds = None
        messages.addMessage(f"RGB→Grayscale conversion complete → {out_r}")
