#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SGTool_GUI.py  –  Standalone multi-tab GUI for SGTool (ArcGIS Pro / standalone).

Changes from v1:
  • Light grey Windows-style theme
  • Input selector shows layers already loaded in the ArcGIS Pro map (+Refresh)
  • No output-path field: output is auto-named  <input>_<suffix>.tif beside input
  • FFT buffer is now passed to every FFT method call (fixes vignetting)
"""

import logging
import os
import re
import sys
import webbrowser
import json
import tempfile
import threading
import traceback
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

# ---------------------------------------------------------------------------
# Bootstrap path
# ---------------------------------------------------------------------------
_HERE   = os.path.dirname(os.path.abspath(__file__))
_SGTOOL = os.path.dirname(_HERE)
for _p in (_HERE, _SGTOOL,
           os.path.join(_SGTOOL, "calcs"),
           os.path.join(_SGTOOL, "calcs", "worms")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from arcgis_utils import (suggest_bspline_settings, raster_to_numpy, numpy_to_raster, make_processor, NODATA,
                          write_provenance, step_record, raster_header_info,
                          grid_info_text, plugin_version, step_suffix,
                          write_mrvbf_outputs, MRVBF_DEFAULTS,
                          list_point_fields, read_points_ogr, grid_points_to_file,
                          map_view_to_cells, read_attributes_ogr, read_polygon_categories)
from calcs.sgt_cancel import OperationCancelled
from sgt_preview_window import PreviewWindow, block_mean
from calcs.sgt_metadata import (read_sgt_metadata, format_sgt_metadata_text,
                                sgt_metadata_xml_text, save_sgt_metadata_xml)
from calcs.workflow import history_steps
from sgt_workflow import plan_steps, run_replay, uniform_field, check_igrf_possible
from calcs import Petrophysics
import numpy as np

# Petrophys plots are built in a worker thread and shown in this GUI's own Tk windows
Petrophysics.USE_PYPLOT = False

# ---------------------------------------------------------------------------
# Light grey palette
# ---------------------------------------------------------------------------
BG        = "#f0f0f0"
FG        = "#1c1c1c"
ACCENT    = "#0078d4"
FRAME_BG  = "#e0e0e0"
ENTRY_BG  = "#ffffff"
BTN_BG    = "#d0d0d0"
BTN_ABG   = "#b8d4f0"
SECTION_FG= "#0060b0"
FONT      = ("Segoe UI", 9)
FONT_BOLD = ("Segoe UI", 9, "bold")
FONT_H    = ("Segoe UI", 10, "bold")
FONT_MONO = ("Consolas",  8)


_wheel_installed = set()


def _install_wheel_scrolling(frame):
    """Scroll the scrollable frame under the pointer with the mouse wheel (one handler for the
    whole window; a Text box, such as the Help tab, scrolls itself and is left alone)."""
    root = frame.winfo_toplevel()
    if root in _wheel_installed:
        return
    _wheel_installed.add(root)

    def on_wheel(e):
        w = root.winfo_containing(e.x_root, e.y_root)
        while w is not None:
            if isinstance(w, tk.Text):
                return
            canvas = getattr(w, "_canvas", None)
            if canvas is not None:
                canvas.yview_scroll(int(-1 * (e.delta / 120)), "units")
                return
            w = w.master

    root.bind_all("<MouseWheel>", on_wheel)


# ============================================================================
# Scrollable frame
# ============================================================================
class ScrollFrame(tk.Frame):
    def __init__(self, parent, **kwargs):
        super().__init__(parent, bg=BG)
        canvas = tk.Canvas(self, bg=BG, highlightthickness=0, **kwargs)
        vsb    = ttk.Scrollbar(self, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=vsb.set)
        vsb.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)
        self.inner = tk.Frame(canvas, bg=BG)
        win_id = canvas.create_window((0, 0), window=self.inner, anchor="nw")

        def _cfg(e):
            canvas.configure(scrollregion=canvas.bbox("all"))
            canvas.itemconfig(win_id, width=canvas.winfo_width())

        self.inner.bind("<Configure>", _cfg)
        canvas.bind("<Configure>", _cfg)
        self._canvas = canvas
        _install_wheel_scrolling(self)


# ============================================================================
# Styled widget helpers
# ============================================================================
def lbl(p, text, bold=False, **kw):
    return tk.Label(p, text=text, bg=BG, fg=FG,
                    font=FONT_BOLD if bold else FONT, **kw)

def entry(p, var=None, w=10, chk_var=None, **kw):
    e = tk.Entry(p, textvariable=var, width=w,
                 bg=ENTRY_BG, fg=FG, insertbackground=FG,
                 relief="sunken", bd=1, font=FONT, **kw)
    if chk_var is not None:
        e.bind("<FocusIn>", lambda _: chk_var.set(True), add=True)
        e.bind("<Key>",     lambda _: chk_var.set(True), add=True)
    return e

def combo(p, values, var=None, w=12, chk_var=None, **kw):
    c = ttk.Combobox(p, values=values, textvariable=var,
                     width=w, state="readonly", font=FONT, **kw)
    if chk_var is not None:
        c.bind("<<ComboboxSelected>>", lambda _: chk_var.set(True), add=True)
    return c

def check(p, text, var, **kw):
    return tk.Checkbutton(p, text=text, variable=var,
                          bg=BG, fg=FG, selectcolor=ENTRY_BG,
                          activebackground=BG, activeforeground=FG,
                          font=FONT, **kw)

def btn(p, text, cmd, w=18, **kw):
    return tk.Button(p, text=text, command=cmd,
                     bg=BTN_BG, fg=FG, activebackground=BTN_ABG,
                     activeforeground=FG, relief="raised", bd=1,
                     font=FONT, width=w, **kw)

def section(p, text):
    return tk.LabelFrame(p, text=text, bg=BG, fg=SECTION_FG,
                         font=FONT_BOLD, relief="groove", bd=1)

def _browse_file(var, filetypes=None, save=False):
    ft = filetypes or [("Raster files", "*.tif *.tiff *.img *.grd *.ers *.mag *.grv"),
                       ("All files", "*.*")]
    path = (filedialog.asksaveasfilename(filetypes=ft, defaultextension=".tif")
            if save else filedialog.askopenfilename(filetypes=ft))
    if path:
        var.set(path)

def _browse_dir(var):
    d = filedialog.askdirectory()
    if d:
        var.set(d)


# ============================================================================
# Tooltip helper
# ============================================================================
class ToolTip:
    """Simple tooltip for any tkinter widget."""
    def __init__(self, widget, text):
        self.widget = widget
        self.text   = text
        self._tw    = None
        widget.bind("<Enter>", self._show, add=True)
        widget.bind("<Leave>", self._hide, add=True)

    def _show(self, event=None):
        if not self.text:
            return
        x = self.widget.winfo_rootx() + 20
        y = self.widget.winfo_rooty() + self.widget.winfo_height() + 4
        self._tw = tk.Toplevel(self.widget)
        self._tw.wm_overrideredirect(True)
        self._tw.wm_geometry(f"+{x}+{y}")
        tk.Label(self._tw, text=self.text, bg="#ffffe0", fg="#1c1c1c",
                 font=("Segoe UI", 8), relief="solid", bd=1,
                 justify="left", padx=4, pady=2).pack()

    def _hide(self, event=None):
        if self._tw:
            self._tw.destroy()
            self._tw = None


def tip(widget, text):
    """Attach a tooltip and return the widget (for chaining)."""
    ToolTip(widget, text)
    return widget


# ---- the Help tab: the QGIS plugin's Help text; **bold** and [label](url) are marked up in it
HELP_TEXT = [
    "**Simple Potential Field Calcs to assist WAXI/Agate Structural Geophysics Course**",
    "[https://waxi4.org](https://waxi4.org)   [https://agate-project.org](https://agate-project.org)",
    "**Example magnetic grid from OMRG Mauritania:**",
    "[grid](http://tectonique.net/sgtools_data/ogrm_usgs_mag_tmi.tif)",
    "**Noddy Imports:**",
    "Existing Noddy mag and grav files can be found at the Atlas of Structural Geophysics: "
    "[https://tectonique.net/asg/](https://tectonique.net/asg/)",
    "New Noddy models can be calculated using the Windows version at "
    "[https://tectonique.net/noddy/OpenNoddy_installer.exe](https://tectonique.net/noddy/OpenNoddy_installer.exe) "
    "or a python wrapper at [https://github.com/cgre-aachen/pynoddy](https://github.com/cgre-aachen/pynoddy)",
    "**Code Repository**",
    "This plugin is available directly within QGIS using the Plugin Manager, however the latest version "
    "(including this ArcGIS Pro version) will always be at this GitHub site:",
    "[https://github.com/swaxi/SGTool](https://github.com/swaxi/SGTool)",
    "**Documentation**",
    "[Help File](https://tectonique.net/sgtools_data/Structural%20Geophysics%20Tools.pdf)",
    "Help File reflects latest changes on GitHub. If a feature is not available in the version you are "
    "using, go to Code Repository for latest version",
    "**Code development**",
    "- Calcs Claude, ChatGPT and Mark Jessell",
    "- Plugin construction - Mark Jessell using QGIS Plugin Builder Plugin "
    "[https://g-sherman.github.io/Qgis-Plugin-Builder/](https://g-sherman.github.io/Qgis-Plugin-Builder/)",
    "- IGRF calculation - using pyIGRF [https://github.com/ciaranbe/pyIGRF](https://github.com/ciaranbe/pyIGRF)",
    "- GRD Loader & Radially averaged power spectrum Fatiando a Terra "
    "[https://www.fatiando.org/](https://www.fatiando.org/)",
    "- Example geophysics data above courtesy of Mauritania Govt. "
    "[https://anarpam.mr/en/](https://anarpam.mr/en/)",
    "- Worming of grids uses Frank Horowitz's bsdwormer "
    "[https://bitbucket.org/fghorow/bsdwormer/](https://bitbucket.org/fghorow/bsdwormer/)",
    "- Multilevel B-Spline (MBA) Gridding — native Python translation of SAGA's grid_spline algorithm "
    "(no external plugin required)",
    "- Differential (Variable) RTP — G.R.J. Cooper & D.R. Cowan (2005) Taylor-series method",
    "- MRVBF/MRRTF (Multiresolution Valley Bottom/Ridge Top Flatness) — J.C. Gallant & T.I. Dowling (2003)",
    "- Euler Deconvolution uses Felipe F. Melo and Valéria C.F. Barbosa's Reliable Euler method "
    "[https://github.com/ffigura/Euler-deconvolution-python](https://github.com/ffigura/Euler-deconvolution-python)",
]

_HELP_MARKUP = re.compile(r"(\*\*[^*]+\*\*|\[[^\]]+\]\([^)]+\))")
_HELP_LINK = re.compile(r"\[([^\]]+)\]\(([^)]+)\)")


def _help_pieces(paragraph):
    """(kind, text, url) pieces of a marked-up help paragraph; kind is "text", "bold" or "link"."""
    pieces = []
    for part in _HELP_MARKUP.split(paragraph):
        if not part:
            continue
        if part.startswith("**"):
            pieces.append(("bold", part[2:-2], None))
        elif part.startswith("["):
            label, url = _HELP_LINK.match(part).groups()
            pieces.append(("link", label, url))
        else:
            pieces.append(("text", part, None))
    return pieces


# ============================================================================
# Main application
# ============================================================================
class SGToolApp(tk.Tk):

    def __init__(self):
        super().__init__()
        self.title("SGTool – Geophysical Processing (ArcGIS)")
        self.configure(bg=BG)
        self.geometry("760x900")
        self.minsize(650, 820)
        self.resizable(True, True)

        self._layer_map      = {}   # display_name → file_path
        self._info_file      = sys.argv[1] if len(sys.argv) > 1 else None
        self._pending_file   = os.path.join(tempfile.gettempdir(), "sgtool_pending.json")
        self._refresh_req    = os.path.join(tempfile.gettempdir(), "sgtool_refresh_req.json")
        self._tls            = threading.local()   # the cancel token of the job a thread runs
        # external preview (see 'External preview' below)
        self.v_preview       = tk.BooleanVar(value=False)
        self._preview_win    = None     # the preview window, while previewing
        self._preview_tab    = None     # tab (0 FFT, 1 Conv) the preview belongs to
        self._preview_gen    = 0        # newest request: older results are dropped
        self._preview_timer  = None
        self._preview_src    = None     # (path, mtime, (array, cell x, cell y))
        self._preview_path   = None
        self._preview_enforcing = False
        self._preview_ending = False    # Apply is ending the preview: keep the ticks
        self._apply_btns     = {}       # tab -> its Apply button (renamed while previewing)
        self.v_preview.trace_add("write", self._preview_var_changed)
        self._active_token   = None                # cancel token of the running job
        self._apply_theme()
        self._build_top_bar()
        self.v_input.trace_add("write", lambda *_: self._update_grid_info())
        self._update_grid_info()
        self._build_notebook()
        self._hook_preview_vars(0)
        for seq in ("<KeyRelease>", "<ButtonRelease-1>", "<<ComboboxSelected>>"):
            self.bind_all(seq, self._preview_event, add="+")
        self._build_status_bar()
        self.after(200, lambda: threading.Thread(
            target=self._refresh_layers, daemon=True).start())

    # ------------------------------------------------------------------
    # Theme
    # ------------------------------------------------------------------
    def _apply_theme(self):
        s = ttk.Style(self)
        s.theme_use("clam")
        s.configure("TNotebook",      background=BG, borderwidth=0)
        s.configure("TNotebook.Tab",  background=FRAME_BG, foreground=FG,
                    padding=[8, 4], font=FONT_BOLD)
        s.map("TNotebook.Tab",
              background=[("selected", ACCENT)],
              foreground=[("selected", "#ffffff")])
        s.configure("TCombobox",
                    fieldbackground=ENTRY_BG, background=BTN_BG,
                    foreground=FG, selectbackground=ACCENT, arrowcolor=FG)
        s.configure("Vertical.TScrollbar",
                    background=BTN_BG, troughcolor=BG, arrowcolor=FG)

    # ------------------------------------------------------------------
    # Top bar: layer selector (ArcGIS Pro map layers) + file browse
    # ------------------------------------------------------------------
    def _build_top_bar(self):
        top = tk.Frame(self, bg=FRAME_BG, pady=5)
        top.pack(fill="x", padx=4, pady=(4, 0))

        self.v_input = tk.StringVar()

        # ---- Row 0: layer dropdown ----
        r0 = tk.Frame(top, bg=FRAME_BG)
        r0.pack(fill="x", padx=6, pady=(2, 1))

        lbl(r0, "Input Layer/File:").pack(side="left")
        self._layer_combo = ttk.Combobox(r0, textvariable=self.v_input,
                                          width=54, font=FONT)
        self._layer_combo.pack(side="left", padx=(4, 2))
        self._layer_combo.bind("<<ComboboxSelected>>", self._on_layer_selected)
        ToolTip(self._layer_combo, "Layer selected for processing\n(raster layers currently loaded in the ArcGIS Pro map)")

        tip(btn(r0, "⟳ Refresh",
            lambda: threading.Thread(target=self._refresh_layers, daemon=True).start(),
            w=10), "Refresh the list of raster layers from the active ArcGIS Pro map").pack(side="left", padx=2)
        tip(btn(r0, "Browse File…", self._browse_input, w=12),
            "Load new file for processing\n*.tif *.tiff *.img *.grd *.ers *.mag *.grv\n"
            ".grd files are automatically converted to GeoTIFF").pack(side="left", padx=2)

        # ---- Row 1: output preview ----
        r1 = tk.Frame(top, bg=FRAME_BG)
        r1.pack(fill="x", padx=6, pady=(1, 3))
        lbl(r1, "Output (auto):").pack(side="left")
        self.v_output_preview = tk.StringVar(value="(select input and enable operations)")
        tk.Label(r1, textvariable=self.v_output_preview,
                 bg=FRAME_BG, fg="#555555", font=FONT_MONO,
                 anchor="w").pack(side="left", padx=4, fill="x", expand=True)

        # ---- Row 2: version and the selected grid's units, CRS, pixels, pixel size ----
        r2 = tk.Frame(top, bg=FRAME_BG)
        r2.pack(fill="x", padx=6, pady=(0, 2))
        self.v_grid_info = tk.StringVar(value="")
        tk.Label(r2, textvariable=self.v_grid_info, bg=FRAME_BG, fg="#333333",
                 font=FONT, anchor="w").pack(side="left")

        top.columnconfigure(0, weight=1)

    def _refresh_layers(self):
        """Populate the layer dropdown — safe to call from any thread."""
        import time
        layer_map = {}
        names     = []

        # Strategy 1: signal the in-process watcher to refresh the layer list,
        # then wait up to 5 s for it to update the info file.
        if self._info_file and os.path.exists(self._info_file):
            try:
                mtime_before = os.path.getmtime(self._info_file)
                with open(self._refresh_req, "w") as f:
                    json.dump({"requested": True}, f)
                for _ in range(10):                    # up to 5 s (10 × 0.5 s)
                    time.sleep(0.5)
                    if os.path.getmtime(self._info_file) > mtime_before:
                        break
            except Exception as e:
                print(f"SGTool: layer refresh signal failed: {e}")

        # Strategy 2: arcpy CURRENT (works only inside the ArcGIS Pro process)
        def _harvest(proj):
            for m in proj.listMaps():
                for lyr in m.listLayers():
                    if lyr.isRasterLayer:
                        display = f"{lyr.name}  [{m.name}]"
                        layer_map[display] = lyr.dataSource
                        names.append(display)

        try:
            import arcpy
            _harvest(arcpy.mp.ArcGISProject("CURRENT"))
        except Exception as e:
            # Expected to fail when not running inside ArcGIS Pro
            print(f"SGTool: arcpy CURRENT project harvest unavailable: {e}")

        # Strategy 3: read the info file (kept fresh by the watcher every ~30 s)
        if not names and self._info_file and os.path.exists(self._info_file):
            try:
                for k, v in json.load(open(self._info_file)).get("layers", {}).items():
                    layer_map[k] = v
                    names.append(k)
            except Exception as e:
                print(f"SGTool: failed to read layer info file: {e}")

        # All tkinter updates must happen on the main thread
        self.after(0, lambda: self._apply_layer_list(layer_map, names))

    def _apply_layer_list(self, layer_map, names):
        self._layer_map = layer_map
        self._layer_combo["values"] = names
        self._petro_refresh_rasters()
        if names and not self.v_input.get():
            self._layer_combo.set(names[0])
            self.v_input.set(self._layer_map[names[0]])
        self._status(f"Layer list refreshed – {len(names)} raster(s) found." if names
                     else "No raster layers found. Use Browse File or re-launch from ArcGIS Pro.")

    def _update_grid_info(self):
        """Version, units, CRS, pixel counts and pixel size of the selected grid."""
        path = self._layer_map.get(self.v_input.get(), self.v_input.get())
        text = ""
        if path and os.path.exists(path):
            text = grid_info_text(raster_header_info(path))
        version = plugin_version()
        self.v_grid_info.set((f"Version {version}   " if version else "") + text)

    def _on_layer_selected(self, event=None):
        sel = self._layer_combo.get()
        if sel in self._layer_map:
            self.v_input.set(self._layer_map[sel])

    def _browse_input(self):
        path = filedialog.askopenfilename(
            filetypes=[("Raster files", "*.tif *.tiff *.img *.grd *.ers *.mag *.grv"),
                       ("All files", "*.*")])
        if not path:
            return
        suffix = os.path.splitext(path)[1].lower()
        if suffix == ".grd":
            converted = self._convert_grd_to_tif(path)
            if converted is None:
                return
            path = converted
        name = os.path.basename(path)
        self._layer_map[name] = path
        current = list(self._layer_combo["values"])
        if name not in current:
            current.insert(0, name)
            self._layer_combo["values"] = current
        self._layer_combo.set(name)
        self.v_input.set(path)
        # Add the browsed file to the ArcGIS Pro map via the pending queue
        self._record_pending(path)

    def _convert_grd_to_tif(self, grd_path):
        """Convert an Oasis Montaj .grd file to GeoTIFF; return tif path or None on error."""
        try:
            from calcs.geosoft_grid_parser import (
                load_oasis_montaj_grid_optimized, extract_proj_str)
            from osgeo import gdal, osr

            epsg = 4326
            xml_path = grd_path + ".xml"
            if os.path.exists(xml_path):
                found = extract_proj_str(xml_path)
                if found:
                    epsg = int(found)

            grid, header, Gdata_type = load_oasis_montaj_grid_optimized(grd_path)
            if Gdata_type == -1:
                messagebox.showwarning("GRD Error",
                    "Cannot read SHORT or INT data types from this .grd file.")
                return None

            base     = os.path.splitext(grd_path)[0]
            tif_path = base + ".tif"

            ordering  = header.get("ordering", 1)
            if ordering == 1:
                nrows, ncols = header["shape_v"], header["shape_e"]
            else:
                nrows, ncols = header["shape_e"], header["shape_v"]
            dx = header["spacing_e"]
            dy = header["spacing_v"]
            x0 = header["x_origin"] - dx / 2
            y0 = header["y_origin"] - dy / 2
            geotransform = [x0, dx, 0, y0, 0, dy]

            driver = gdal.GetDriverByName("GTiff")
            if os.path.exists(tif_path):
                try:
                    os.remove(tif_path)
                except OSError:
                    pass
            ds = driver.Create(tif_path, ncols, nrows, 1, Gdata_type)
            ds.SetGeoTransform(geotransform)
            srs = osr.SpatialReference()
            srs.ImportFromEPSG(epsg)
            ds.SetProjection(srs.ExportToWkt())
            ds.GetRasterBand(1).WriteArray(grid)
            ds.GetRasterBand(1).SetNoDataValue(NODATA)
            ds.FlushCache()
            ds = None
            self._status(f"GRD converted → {os.path.basename(tif_path)}")
            return tif_path
        except Exception as e:
            messagebox.showerror("GRD Conversion Error", str(e))
            return None

    def _browse_import_file(self):
        """Browse for a point/line data file and auto-populate X/Y column dropdowns."""
        path = filedialog.askopenfilename(
            filetypes=[("Data files", "*.csv *.dat *.xyz *.txt"), ("All", "*.*")])
        if not path:
            return
        self.v_import_file.set(path)
        try:
            import csv as csvmod
            with open(path, newline="", encoding="utf-8-sig") as f:
                sample = f.read(4096)
                f.seek(0)
                try:
                    dialect = csvmod.Sniffer().sniff(sample, delimiters=",\t ")
                except csvmod.Error:
                    dialect = csvmod.excel
                reader  = csvmod.DictReader(f, dialect=dialect)
                headers = list(reader.fieldnames or [])
            self._xcol_combo["values"] = headers
            self._ycol_combo["values"] = headers
            # Auto-select likely X/Y columns
            for h in headers:
                if h.lower() in ("x", "lon", "longitude", "easting", "long_x", "long", "e"):
                    self.v_x_col.set(h)
                    break
            for h in headers:
                if h.lower() in ("y", "lat", "latitude", "northing", "lat_y", "n"):
                    self.v_y_col.set(h)
                    break
        except Exception as e:
            print(f"SGTool: failed to parse import file headers: {e}")

    # ------------------------------------------------------------------
    # Pending-layer queue (watcher in SGTool.pyt picks this up in-process)
    # ------------------------------------------------------------------
    def _record_pending(self, path):
        """Append a written raster path to the pending queue."""
        try:
            if os.path.exists(self._pending_file):
                with open(self._pending_file) as f:
                    paths = json.load(f)
            else:
                paths = []
            if path not in paths:
                paths.append(path)
            with open(self._pending_file, "w") as f:
                json.dump(paths, f)
        except Exception as e:
            print(f"SGTool: failed to record pending layer path: {e}")

    def _on_output_written(self, path):
        """Update status bar and add newly written raster to the layer dropdown."""
        name = os.path.basename(path)
        self._status(f"Written → {name}")
        self._layer_map[name] = path
        vals = list(self._layer_combo["values"])
        if name not in vals:
            vals.insert(0, name)
            self._layer_combo["values"] = vals
            self._petro_refresh_rasters()

    # ------------------------------------------------------------------
    # Auto output-path derivation
    # ------------------------------------------------------------------
    # ------------------------------------------------------------------
    # External preview
    #
    # ArcGIS Pro gives a separate Python window no map hooks, so the preview has its own
    # window (sgt_preview_window.PreviewWindow). One filter at a time is ticked; the filter
    # is run on what is in view (or on the whole grid, reduced) by the very same code that
    # Apply uses (_fft_steps / _conv_steps), with a _save that keeps the first result
    # instead of writing it. Preview follows you between the FFT Filters and Conv + Stats
    # tabs; Apply Preview saves the ticked filter at full resolution.
    # ------------------------------------------------------------------
    PREVIEW_FILTERS = {
        0: ("v_rtp", "v_vint", "v_cont", "v_dirclean", "v_regrem", "v_bp", "v_hlp", "v_agc",
            "v_deriv", "v_thg", "v_ta", "v_as_"),
        1: ("v_mean", "v_med", "v_gauss", "v_dconv", "v_sun", "v_ss_min", "v_ss_max", "v_ss_var",
            "v_ss_std", "v_ss_skew", "v_ss_kurt", "v_ss_aniso", "v_ss_chain", "v_ss_stream",
            "v_mrvbf", "v_pca", "v_ica", "v_euler"),
        3: ("v_app_sus", "v_app_dens"),  # Petrophys tab: at most one ticked
        4: ("v_thresh_en",),  # Utils tab: Threshold to NaN
    }
    PREVIEW_NOT_POSSIBLE = ("v_ss_aniso", "v_ss_chain", "v_ss_stream", "v_mrvbf",
                            "v_pca", "v_ica", "v_euler")
    PREVIEW_VIEW_SIDE = 1000    # most cells (longest side) calculated for the view
    PREVIEW_WHOLE_SIDE = 600    # the whole grid is reduced to about this
    APPLY_TEXT = {0: "Apply FFT Processing", 1: "Apply Conv + Stats", 3: "Apply Petrophysics",
                  4: "Apply Threshold \u2192 NaN"}

    def _filter_var(self, name):
        return getattr(self, name, None)

    @staticmethod
    def _untick(var):
        """Untick a filter's tick box."""
        if isinstance(var, tk.BooleanVar):
            var.set(False)

    def _ticked(self, tab):
        return [n for n in self.PREVIEW_FILTERS[tab]
                if self._filter_var(n) is not None and self._filter_var(n).get()]

    def _hook_preview_vars(self, tab):
        """Watch a tab's filter tick boxes (called once the tab is built)."""
        for name in self.PREVIEW_FILTERS.get(tab, ()):  # (only two tabs can preview)
            var = self._filter_var(name)
            if var is not None:
                var.trace_add("write", lambda *_a, t=tab, n=name: self._preview_filter_changed(t, n))
        self._update_apply_text()

    def _update_apply_text(self):
        on = self.v_preview.get()
        for tab, button in self._apply_btns.items():
            button.config(text="Apply Preview" if on else self.APPLY_TEXT[tab])

    def _current_tab(self):
        try:
            return self.nb.index(self.nb.select())
        except (AttributeError, tk.TclError):
            return 0

    def _preview_var_changed(self, *_):
        self._update_apply_text()
        if self.v_preview.get():
            self._preview_start()
        else:
            self._preview_stop(untick=not self._preview_ending)

    def _preview_start(self):
        tab = self._current_tab()
        self._preview_tab = tab if tab in self.PREVIEW_FILTERS else None
        if self._preview_tab is not None:
            self._enforce_single(self._preview_tab, None)
        if self._preview_win is None or not self._preview_win.winfo_exists():
            self._preview_win = PreviewWindow(self, self._preview_view_changed, self._preview_closed,
                                              self._preview_match_map)
        self._preview_path = None
        self._preview_schedule(50)

    def _preview_stop(self, untick):
        """Close the preview. With untick, also untick the filters (Apply keeps them)."""
        if self._preview_timer is not None:
            self.after_cancel(self._preview_timer)
            self._preview_timer = None
        self._preview_gen += 1
        win, self._preview_win = self._preview_win, None
        if win is not None:
            try:
                win.destroy()
            except tk.TclError:
                pass
        if untick:
            self._preview_enforcing = True
            try:
                for tab in self.PREVIEW_FILTERS:
                    for name in self.PREVIEW_FILTERS[tab]:
                        var = self._filter_var(name)
                        if var is not None:
                            self._untick(var)
            finally:
                self._preview_enforcing = False
        self._preview_tab = None

    def _preview_closed(self):
        """The user closed the preview window."""
        self.v_preview.set(False)

    def _end_preview_for_apply(self):
        """Apply while previewing: close the preview but keep the ticked filter, which
        Apply is about to run (and unticks when it finishes)."""
        if self.v_preview.get():
            self._preview_ending = True
            try:
                self.v_preview.set(False)
            finally:
                self._preview_ending = False

    def _enforce_single(self, tab, keep):
        """Only one filter stays ticked while previewing: `keep`, else the first ticked
        of the tab; the filters of every other tab are unticked too (e.g. a Petrophys
        calculation ticked before Preview was turned on)."""
        ticked = self._ticked(tab)
        if keep is None and ticked:
            keep = ticked[0]
        self._preview_enforcing = True
        try:
            for name in ticked:
                if name != keep:
                    self._filter_var(name).set(False)
            for other_tab, names in self.PREVIEW_FILTERS.items():
                if other_tab != tab:
                    for name in names:
                        var = self._filter_var(name)
                        if var is not None:
                            self._untick(var)
        finally:
            self._preview_enforcing = False

    def _preview_filter_changed(self, tab, name):
        if (not self.v_preview.get() or self._preview_tab != tab or self._preview_enforcing):
            return
        if self._filter_var(name).get():
            self._enforce_single(tab, name)
        self._preview_schedule()

    def _preview_tab_changed(self, idx):
        """Leaving a previewed tab unticks its filters; Preview stays on and is shown for
        the tab arrived at (tabs without a preview just wait)."""
        self._update_apply_text()
        if not self.v_preview.get():
            return
        old = self._preview_tab
        if old is not None and old != idx:
            self._preview_enforcing = True
            try:
                for name in self.PREVIEW_FILTERS[old]:
                    var = self._filter_var(name)
                    if var is not None:
                        self._untick(var)
            finally:
                self._preview_enforcing = False
        self._preview_tab = idx if idx in self.PREVIEW_FILTERS else None
        if self._preview_tab is not None:
            self._enforce_single(self._preview_tab, None)
        self._preview_schedule(50)

    def _preview_event(self, event):
        """A key or click anywhere in the app: settings may have changed."""
        if not self.v_preview.get() or self._preview_tab is None:
            return
        try:
            if self._preview_win is not None and event.widget.winfo_toplevel() is self._preview_win:
                return
        except (AttributeError, tk.TclError, KeyError):
            pass
        self._preview_schedule()

    def _map_view(self, path):
        """(the ArcGIS Pro map view as a box of cells of the grid at `path`, the map canvas's
        width / height), or None. The toolbox tool that launched the GUI writes the view to a
        file every couple of seconds."""
        try:
            with open(os.path.join(tempfile.gettempdir(), "sgtool_view.json")) as f:
                view = json.load(f)
            if map_view_to_cells(view, path) is None:   # not looking at this grid
                return None
            xmin, ymin, xmax, ymax = view["extent"]
            return (map_view_to_cells(view, path, clip=False),
                    (xmax - xmin) / (ymax - ymin))
        except Exception:
            return None

    def _preview_show_map_view(self, path):
        """Show in the preview window what the ArcGIS Pro map shows of this grid, in a window
        of the same shape as the map canvas. False if there is no map view over the grid."""
        win = self._preview_win
        found = self._map_view(path) if win is not None else None
        if found is None:
            return False
        box, aspect = found
        win.match_canvas_aspect(aspect)
        win.set_view(box)
        return True

    def _preview_match_map(self):
        """The Match map view button."""
        win = self._preview_win
        if win is None or not win.winfo_exists():
            return
        path = self._resolve_input()
        if not (path and os.path.exists(path) and self._preview_show_map_view(path)):
            win.set_status("No ArcGIS Pro map view over this grid (the GUI must be launched "
                           "from the SGTool toolbox, with a map open over the grid)")

    def _preview_view_changed(self):
        self._preview_schedule(250)

    def _preview_schedule(self, delay=400):
        if self._preview_timer is not None:
            self.after_cancel(self._preview_timer)
        self._preview_timer = self.after(delay, self._preview_refresh)

    def _preview_refresh(self):
        self._preview_timer = None
        win = self._preview_win
        if win is None or not win.winfo_exists():
            return
        tab = self._preview_tab
        if tab is None:
            win.clear("The preview is on the FFT Filters, Conv + Stats, Petrophys and Utils "
                      "(threshold) tabs")
            win.set_status("")
            return
        path = self._resolve_input()
        if not path or not os.path.exists(path):
            win.clear("Select an input raster")
            return
        if path != self._preview_path:  # a new grid: show all of it
            header = raster_header_info(path)
            if header is None:
                win.clear("Cannot read " + os.path.basename(path))
                return
            self._preview_path = path
            win.set_grid(header["height"], header["width"], force=True,
                         ratio=header["cell_y"] / header["cell_x"])
            self._preview_show_map_view(path)  # start where the ArcGIS Pro map is looking
        ticked = self._ticked(tab)
        if not ticked:
            win.clear("Tick a filter to preview")
            win.set_status("")
            return
        if ticked[0] in self.PREVIEW_NOT_POSSIBLE or (
                ticked[0] == "v_rtp" and self.v_rtp_type.get() == "Diff. RTP"):
            win.clear("This filter cannot be previewed (Apply Preview still runs it)")
            win.set_status("")
            return
        self._preview_gen += 1
        gen = self._preview_gen
        win.set_status("Calculating…")
        threading.Thread(target=self._preview_worker,
                         args=(gen, tab, path, win.v_mode.get(), win.view_box()),
                         daemon=True).start()

    def _preview_worker(self, gen, tab, path, mode, view):
        try:
            result, error = self._compute_preview(tab, path, mode, view), None
        except Exception as e:
            result, error = None, str(e) or type(e).__name__
        try:
            self.after(0, self._preview_done, gen, result, error)
        except (RuntimeError, tk.TclError):
            pass  # the app is closing

    def _preview_done(self, gen, result, error):
        win = self._preview_win
        if gen != self._preview_gen or win is None or not win.winfo_exists():
            return  # a newer request is under way, or the preview was closed
        if error:
            win.clear(error)
            win.set_status("")
            return
        win.show_result(result["filtered"], result["original"], result["box"], result["std_clip"])
        win.set_status(result["text"])

    def _preview_source(self, path):
        """(array, cell x, cell y) of the input grid, read once per file."""
        try:
            stamp = os.path.getmtime(path)
        except OSError:
            stamp = None
        if self._preview_src is not None and self._preview_src[:2] == (path, stamp):
            return self._preview_src[2]
        arr, _nodata, _llx, _lly, cx, cy, _sr = raster_to_numpy(path)
        self._preview_src = (path, stamp, (arr, cx, cy))
        return self._preview_src[2]

    def _compute_preview(self, tab, path, mode, view):
        """Run the ticked filter of `tab` on the part of the grid the preview shows.

        mode "view" calculates what is in `view` = (x0, y0, x1, y1) cells at full detail
        (reduced if it is large), "whole" the whole grid reduced. Returns a dict with the
        filtered grid, the input it came from, where it lies (box, in cells of the full
        grid), whether to stretch by standard deviation, and a line of text."""
        import math
        import time

        t0 = time.time()
        arr, cx, cy = self._preview_source(path)
        rows, cols = arr.shape
        if mode == "whole":
            x0, y0, x1, y1 = 0, 0, cols, rows
            side = self.PREVIEW_WHOLE_SIDE
        else:
            vx0, vy0, vx1, vy1 = view
            x0, y0 = int(max(0, math.floor(vx0))), int(max(0, math.floor(vy0)))
            x1, y1 = int(min(cols, math.ceil(vx1))), int(min(rows, math.ceil(vy1)))
            if x1 - x0 < 16:
                x0, x1 = max(0, min(x0, cols - 16)), min(cols, max(x1, 16))
            if y1 - y0 < 16:
                y0, y1 = max(0, min(y0, rows - 16)), min(rows, max(y1, 16))
            side = self.PREVIEW_VIEW_SIDE
        crop = arr[y0:y1, x0:x1]
        k = max(1, int(math.ceil(max(crop.shape) / side)))
        sub = block_mean(crop, k)
        ecx, ecy = cx * k, cy * k
        buf = max(1, min(int(self.v_buf.get() or 10), *sub.shape))
        proc = make_processor(ecx, ecy, buf)

        class _Done(Exception):
            """The first result is all the preview needs."""

        got = {}

        def keep(result, suffix, step=None):
            got["result"] = np.asarray(result, dtype=float)
            got["step"] = step
            raise _Done()

        count = [0]
        try:
            if tab == 0:
                self._fft_steps(path, sub, ecx, ecy, buf, proc, count, keep)
            elif tab == 3:
                self._petro_steps(sub, ecx, ecy, buf, proc, keep)
            elif tab == 4:
                self._threshold_step(sub, keep)
            else:
                self._conv_steps(path, sub, NODATA, 0.0, 0.0, ecx, ecy, None,
                                 os.path.splitext(path)[0], ".tif", count, keep)
        except _Done:
            pass
        if "result" not in got:
            raise RuntimeError("Nothing to preview for these settings")
        operation = got["step"][0] if got["step"] else "Result"
        text = f"{operation}:  {sub.shape[1]} x {sub.shape[0]} cells"
        if k > 1:
            text += f", block-averaged {k} x ({ecx:g} map units per cell)"
        text += f"   ({time.time() - t0:.1f} s)"
        try:
            if self.v_dirclean.get() and ecx > float(self.v_dc_wl.get()):
                text += "   - cell size exceeds the line spacing: zoom in"
        except (AttributeError, ValueError):
            pass
        return {"filtered": got["result"], "original": sub, "box": (x0, y0, x1, y1),
                "std_clip": operation != "Sun shading", "text": text}

    def _planned_suffixes(self):
        """File-name suffixes (QGIS naming) of the results the ticked operations will
        write, so the 'Output (auto)' line shows the real names."""
        def val(name, default=None):
            var = getattr(self, name, None)
            return var.get() if var is not None else default

        def num(name, default):
            try:
                return float(val(name, default))
            except (TypeError, ValueError):
                return default

        out = []

        def add(spec_id, **values):
            out.append(step_suffix(step_record(spec_id, **values)))

        try:
            if val("v_rtp"):
                rt = val("v_rtp_type", "Pole")
                if rt == "Diff. RTP":
                    out += ["_DRTP", "_DRTP_inc", "_DRTP_dec"]
                else:
                    add("reduction_to_equator" if rt == "Equator" else "reduction_to_pole")
            if val("v_vint"):
                add("vertical_integration")
            if val("v_cont"):
                add("continuation", DIRECTION=val("v_cont_dir", "up"), HEIGHT=num("v_cont_h", 500.0))
            if val("v_dirclean"):
                add("directional_line_noise")
            if val("v_regrem"):
                add("remove_regional",
                    ORDER="1st order" if val("v_regorder", "1st") == "1st" else "2nd order")
            if val("v_bp"):
                add("band_pass", LOW_CUT=num("v_bp_lo", 50000.0), HIGH_CUT=num("v_bp_hi", 5000.0))
            if val("v_hlp"):
                add("high_low_pass", TYPE=val("v_hlp_type", "High"), CUTOFF=num("v_hlp_cut", 5000.0))
            if val("v_agc"):
                add("automatic_gain_control")
            if val("v_deriv"):
                add("derivative", DIRECTION=val("v_deriv_dir", "z"), POWER=num("v_deriv_pow", 1.0))
            if val("v_thg"):
                add("total_horizontal_gradient")
            if val("v_ta"):
                add("tilt_angle")
            if val("v_as_"):
                add("analytic_signal")
            if val("v_mean"):
                add("mean_filter")
            if val("v_med"):
                add("median_filter")
            if val("v_gauss"):
                add("gaussian_filter")
            if val("v_dconv"):
                add("directional_filter")
            if val("v_sun"):
                add("sun_shading")
            for attr, stat in (("v_ss_min", "min"), ("v_ss_max", "max"), ("v_ss_var", "variance"),
                               ("v_ss_std", "std"), ("v_ss_skew", "skewness"),
                               ("v_ss_kurt", "kurtosis")):
                if val(attr):
                    add("windowed_statistic", STATISTIC=stat)
            if val("v_ss_aniso"):
                out += ["_SS_AnisoMag", "_SS_AnisoOrient"]
            if val("v_ss_chain"):
                out.append("_SS_ChainLen")
            if val("v_ss_stream"):
                out.append("_SS_StreamLen")
            if val("v_mrvbf"):
                out += ["_MRVBF", "_MRRTF", "_slope", "_MRVBF_RGB"]
            if val("v_pca"):
                out.append("_PCA")
            if val("v_ica"):
                out.append("_ICA")
            if val("v_euler"):
                out.append(f"_estimates_SI_{int(num('v_eu_si', 1.0))}")
            if self._current_tab() == 3:  # Petrophys: at most one ticked
                if val("v_app_sus"):
                    out.append("_AppSus")
                elif val("v_app_dens"):
                    out.append("_AppDens")
        except Exception:
            logging.getLogger(__name__).debug("Ignored: a half-typed entry: show what could be worked out", exc_info=True)
        return [p for p in out if p]

    def _compute_output_path(self):
        """Build output path from input path + suffix describing active operations."""
        in_r = self.v_input.get()
        # Resolve layer name → path if needed
        if in_r in self._layer_map:
            in_r = self._layer_map[in_r]
        if not in_r or not os.path.exists(in_r):
            return None

        parts = self._planned_suffixes()

        base, ext = os.path.splitext(in_r)
        if not ext:
            ext = ".tif"
        bn = os.path.basename(base)
        if parts:
            if len(parts) <= 3:
                preview = ",  ".join(f"{bn}{p}{ext}" for p in parts)
            else:
                preview = (f"{len(parts)} outputs: "
                           + ", ".join(parts[:3]) + f"  +{len(parts)-3} more")
        else:
            preview = "(no operations selected)"
        self.v_output_preview.set(preview)
        return in_r

    # ------------------------------------------------------------------
    # Notebook
    # ------------------------------------------------------------------
    def _build_notebook(self):
        self.nb = ttk.Notebook(self)
        self.nb.pack(fill="both", expand=True, padx=4, pady=4)

        self._tab_specs = [
            ("FFT Filters",     self._tab_fft),
            ("Conv + Stats",    self._tab_conv_stats),
            ("Grid + Wavelets", self._tab_grid_wavelets),
            ("Petrophys",       self._tab_petro),
            ("Utils",           self._tab_utils),
            ("Help",            self._tab_help),
        ]

        # Build first tab now so the window appears with content
        title, builder = self._tab_specs[0]
        sf = ScrollFrame(self.nb)
        builder(sf.inner)
        self.nb.add(sf, text=title)

        # Add placeholder frames for remaining tabs; build them on first switch
        self._tab_frames   = {0: sf}
        self._tab_built    = {0}
        for i, (ttl, _) in enumerate(self._tab_specs[1:], start=1):
            if ttl == "Help":  # a text box with its own scroll bar that fills the tab
                ph = tk.Frame(self.nb, bg=BG)
                ph.inner = ph
            else:
                ph = ScrollFrame(self.nb)
            self.nb.add(ph, text=ttl)
            self._tab_frames[i] = ph

        self.nb.bind("<<NotebookTabChanged>>", self._on_tab_changed)

    def _on_tab_changed(self, event=None):
        idx = self.nb.index(self.nb.select())
        if idx not in self._tab_built:
            self._tab_built.add(idx)
            _, builder = self._tab_specs[idx]
            builder(self._tab_frames[idx].inner)
            self._hook_preview_vars(idx)
        self._preview_tab_changed(idx)

    # ------------------------------------------------------------------
    # Status bar
    # ------------------------------------------------------------------
    def _build_status_bar(self):
        bar = tk.Frame(self, bg=FRAME_BG, height=22, relief="sunken", bd=1)
        bar.pack(fill="x", side="bottom")
        self.v_status = tk.StringVar(value="Ready.")
        tk.Label(bar, textvariable=self.v_status,
                 bg=FRAME_BG, fg=FG, font=FONT, anchor="w").pack(side="left", padx=6)
        self._cancel_btn = tip(btn(bar, "Cancel", self._cancel_job, w=8),
            "Stop the running calculation now. A step that is part-way through an FFT\n"
            "finishes in the background, but its result is discarded and nothing is saved.")
        self._cancel_btn.config(state="disabled")
        self._cancel_btn.pack(side="right", padx=4, pady=1)

    def _status(self, msg):
        self.v_status.set(msg)
        self.update_idletasks()

    # ======================================================================
    # TAB 1 – FFT FILTERS
    # ======================================================================
    def _tab_fft(self, p):
        # ---- Grav/Mag ----
        gm = section(p, "Grav/Mag Filters")
        gm.pack(fill="x", padx=8, pady=4)

        r0 = tk.Frame(gm, bg=BG); r0.pack(fill="x", padx=4, pady=2)
        self.v_rtp      = tk.BooleanVar()
        self.v_rtp_type = tk.StringVar(value="Pole")
        tip(check(r0, "RTP/RTE", self.v_rtp),
            "Reduction to pole or equator\n"
            "Transforms magnetic data so anomalies appear directly above sources\n"
            "Corrects asymmetry caused by the Earth's inclined field").pack(side="left")
        tip(combo(r0, ["Pole", "Equator", "Diff. RTP"], self.v_rtp_type, w=10, chk_var=self.v_rtp),
            "Choose Pole (high latitudes >20°), Equator (low latitudes <20°),\n"
            "or Diff. RTP (variable RTP accounting for spatial inc/dec variation\n"
            "across the survey area — Cooper & Cowan 2005)\n"
            "Diff. RTP uses the IGRF Year field; Inc/Dec are computed automatically from corners").pack(side="left", padx=4)

        r1 = tk.Frame(gm, bg=BG); r1.pack(fill="x", padx=4, pady=2)
        self.v_inc = tk.StringVar(value="-60.0")
        self.v_dec = tk.StringVar(value="0.0")
        self.v_intensity = tk.StringVar(value="50000")
        lbl(r1, "Inc:").pack(side="left")
        tip(entry(r1, self.v_inc, 7),
            "Magnetic inclination [degrees from horizontal]").pack(side="left", padx=2)
        lbl(r1, "Dec:").pack(side="left", padx=(8,0))
        tip(entry(r1, self.v_dec, 7),
            "Magnetic declination [degrees clockwise from North]").pack(side="left", padx=2)
        lbl(r1, "Int (nT):").pack(side="left", padx=(8,0))
        tip(entry(r1, self.v_intensity, 8),
            "Survey intensity in nT").pack(side="left", padx=2)

        r2 = tk.Frame(gm, bg=BG); r2.pack(fill="x", padx=4, pady=2)
        self.v_igrf_lat  = tk.StringVar(value="-25.0")
        self.v_igrf_lon  = tk.StringVar(value="130.0")
        self.v_igrf_year = tk.StringVar(value="2025.0")
        lbl(r2, "Lat:").pack(side="left")
        entry(r2, self.v_igrf_lat, 7).pack(side="left", padx=2)
        lbl(r2, "Lon:").pack(side="left", padx=(8,0))
        entry(r2, self.v_igrf_lon, 7).pack(side="left", padx=2)
        lbl(r2, "Year:").pack(side="left", padx=(8,0))
        entry(r2, self.v_igrf_year, 7).pack(side="left", padx=2)
        tip(btn(r2, "Calc IGRF", self._calc_igrf, w=10),
            "Calculate IGRF Inclination & Declination\nbased on centroid of selected grid and specified date").pack(side="left", padx=8)

        r3 = tk.Frame(gm, bg=BG); r3.pack(fill="x", padx=4, pady=2)
        self.v_vint = tk.BooleanVar()
        tip(check(r3, "Vertical Integration (Pseudo-gravity if applied to RTP/RTE data)", self.v_vint),
            "Vertical Integration:\nWhen applied to RTE/P result converts magnetic anomalies into gravity-like anomalies\n"
            "Also good for stitched grids with very different line spacing\n"
            "Requires a metre-based projection").pack(side="left")

        r4 = tk.Frame(gm, bg=BG); r4.pack(fill="x", padx=4, pady=2)
        self.v_cont     = tk.BooleanVar()
        self.v_cont_dir = tk.StringVar(value="up")
        self.v_cont_h   = tk.StringVar(value="1000.0")
        tip(check(r4, "Continuation", self.v_cont),
            "Upward or downward continuation\n"
            "Upward attenuates high-frequency noise and shallow features\n"
            "Downward enhances shallow or high-frequency anomalies").pack(side="left")
        tip(combo(r4, ["up", "down"], self.v_cont_dir, w=6, chk_var=self.v_cont),
            "Select direction of continuation").pack(side="left", padx=4)
        lbl(r4, "Height:").pack(side="left")
        tip(entry(r4, self.v_cont_h, 9, chk_var=self.v_cont),
            "Amount of continuation [m only]").pack(side="left", padx=2)

        # ---- Frequency ----
        ff = section(p, "Frequency Filters")
        ff.pack(fill="x", padx=8, pady=4)

        r5 = tk.Frame(ff, bg=BG); r5.pack(fill="x", padx=4, pady=2)
        tip(btn(r5, "Radial Power Spectrum", self._show_power_spectrum, w=22),
            "Displays Radial Averaged Power Spectrum\n(wavenumber vs ln power)").pack(side="left")

        r6 = tk.Frame(ff, bg=BG); r6.pack(fill="x", padx=4, pady=2)
        self.v_dirclean  = tk.BooleanVar()
        self.v_dc_az     = tk.StringVar(value="0.0")
        self.v_dc_wl     = tk.StringVar(value="2000.0")
        self.v_dc_wl_max = tk.StringVar(value="")
        self.v_dc_scale  = tk.StringVar(value="1.0")
        self.v_dc_wedge  = tk.StringVar(value="45")
        tip(check(r6, "Directional Butterworth", self.v_dirclean),
            "Filter (DirCos + Butterworth) to remove a specific direction and wavelength\n"
            "Useful for filtering flight line noise").pack(side="left")
        lbl(r6, "Az:").pack(side="left", padx=(8,0))
        tip(entry(r6, self.v_dc_az, 5, chk_var=self.v_dirclean),
            "Azimuth of high-frequency noise to be filtered\n"
            "(degrees clockwise from North)").pack(side="left", padx=2)
        lbl(r6, "Line spacing min/max:").pack(side="left", padx=(6,0))
        tip(entry(r6, self.v_dc_wl, 7, chk_var=self.v_dirclean),
            "Minimum line spacing in the survey (map units)\n"
            "Noise band is 2× min spacing to 10× max spacing").pack(side="left", padx=2)
        tip(entry(r6, self.v_dc_wl_max, 7, chk_var=self.v_dirclean),
            "Maximum line spacing in the survey (map units)\n"
            "Leave blank if spacing is constant (uses the min value)").pack(side="left", padx=2)
        lbl(r6, "Scale:").pack(side="left")
        tip(entry(r6, self.v_dc_scale, 5, chk_var=self.v_dirclean),
            "Multiplier applied to the noise estimate before it is\n"
            "subtracted from the original grid").pack(side="left", padx=2)
        lbl(r6, "Wedge:").pack(side="left", padx=(6,0))
        tip(entry(r6, self.v_dc_wedge, 4, chk_var=self.v_dirclean),
            "Half-width of the directional wedge in degrees: noise whose direction is\n"
            "within this angle of the azimuth is kept in the noise estimate, with a\n"
            "cosine taper. Wider keeps more of the noise (and more of the geology).").pack(side="left", padx=2)

        r7 = tk.Frame(ff, bg=BG); r7.pack(fill="x", padx=4, pady=2)
        self.v_regrem   = tk.BooleanVar()
        self.v_regorder = tk.StringVar(value="1st")
        tip(check(r7, "Remove Regional", self.v_regrem),
            "Remove regional based on 1st or 2nd order polynomial").pack(side="left")
        combo(r7, ["1st", "2nd"], self.v_regorder, w=5, chk_var=self.v_regrem).pack(side="left", padx=4)
        lbl(r7, "order").pack(side="left")

        r8 = tk.Frame(ff, bg=BG); r8.pack(fill="x", padx=4, pady=2)
        self.v_bp    = tk.BooleanVar()
        self.v_bp_lo = tk.StringVar(value="5000.0")
        self.v_bp_hi = tk.StringVar(value="50000.0")
        self.v_bp_w  = tk.StringVar(value="5000.0")
        tip(check(r8, "Band Pass", self.v_bp),
            "Band pass filter (BP)\nIsolates specific wavelength features.").pack(side="left")
        lbl(r8, "Low:").pack(side="left", padx=(8,0))
        tip(entry(r8, self.v_bp_lo, 9, chk_var=self.v_bp),
            "Low wavelength cutoff [m or other length unit]").pack(side="left", padx=2)
        lbl(r8, "High:").pack(side="left")
        tip(entry(r8, self.v_bp_hi, 9, chk_var=self.v_bp),
            "High wavelength cutoff [m or other length unit]").pack(side="left", padx=2)
        lbl(r8, "Width:").pack(side="left")
        tip(entry(r8, self.v_bp_w, 5, chk_var=self.v_bp),
            "Width of cosine rolloff [m or other length unit]\n"
            "Start with cutoff value and increase to reduce ringing\n"
            "0 = step cutoff").pack(side="left", padx=2)

        r9 = tk.Frame(ff, bg=BG); r9.pack(fill="x", padx=4, pady=2)
        self.v_hlp      = tk.BooleanVar()
        self.v_hlp_type = tk.StringVar(value="High")
        self.v_hlp_cut  = tk.StringVar(value="10000.0")
        self.v_hlp_w    = tk.StringVar(value="5000.0")
        tip(check(r9, "High/Low Pass", self.v_hlp),
            "High or Low pass filter\n"
            "High pass (HP) isolates short wavelength features\n"
            "Low pass (LP) isolates long wavelength features").pack(side="left")
        tip(combo(r9, ["High", "Low"], self.v_hlp_type, w=5, chk_var=self.v_hlp),
            "Cut off type").pack(side="left", padx=4)
        lbl(r9, "Cutoff:").pack(side="left")
        tip(entry(r9, self.v_hlp_cut, 9, chk_var=self.v_hlp),
            "Cutoff wavelength [m or other length unit]").pack(side="left", padx=2)
        lbl(r9, "Width:").pack(side="left")
        tip(entry(r9, self.v_hlp_w, 5, chk_var=self.v_hlp),
            "Width of cosine rolloff [m or other length unit]\n"
            "Start with cutoff value and increase to reduce ringing\n"
            "0 = step cutoff").pack(side="left", padx=2)

        r10 = tk.Frame(ff, bg=BG); r10.pack(fill="x", padx=4, pady=2)
        self.v_agc     = tk.BooleanVar()
        self.v_agc_win = tk.StringVar(value="15")
        self._setup_odd_trace(self.v_agc_win)
        tip(check(r10, "AGC", self.v_agc),
            "Automatic Gain Control (AGC) or Amplitude Normalisation\n"
            "Highlights short wavelength/low amplitude features").pack(side="left")
        lbl(r10, "Window (px, odd):").pack(side="left", padx=(8,0))
        tip(entry(r10, self.v_agc_win, 5, chk_var=self.v_agc),
            "Window size for normalisation").pack(side="left", padx=2)

        # ---- Gradient ----
        gf = section(p, "Gradient Filters")
        gf.pack(fill="x", padx=8, pady=4)

        r11 = tk.Frame(gf, bg=BG); r11.pack(fill="x", padx=4, pady=2)
        self.v_deriv     = tk.BooleanVar()
        self.v_deriv_dir = tk.StringVar(value="z")
        self.v_deriv_pow = tk.StringVar(value="1.0")
        tip(check(r11, "Derivative", self.v_deriv),
            "Calculate derivative (d + power + direction) parallel to x, y or z\n"
            "Highlights near-surface/short-wavelength features").pack(side="left")
        tip(combo(r11, ["x","y","z"], self.v_deriv_dir, w=4, chk_var=self.v_deriv),
            "Select derivative direction").pack(side="left", padx=4)
        lbl(r11, "Power:").pack(side="left")
        tip(entry(r11, self.v_deriv_pow, 5, chk_var=self.v_deriv),
            "Power of derivative").pack(side="left", padx=2)

        self.v_thg = tk.BooleanVar()
        self.v_ta  = tk.BooleanVar()
        self.v_as_ = tk.BooleanVar()
        r12a = tk.Frame(gf, bg=BG); r12a.pack(fill="x", padx=4, pady=2)
        tip(check(r12a, "Total Horizontal Gradient", self.v_thg),
            "Total Horizontal Gradient Calculation (THG)").pack(side="left")
        r12b = tk.Frame(gf, bg=BG); r12b.pack(fill="x", padx=4, pady=2)
        tip(check(r12b, "Tilt Angle", self.v_ta),
            "Tilt Derivative (TD)\nEnhances edges and detects shallow sources\n"
            "Tends to overconnect structural features").pack(side="left")
        r12c = tk.Frame(gf, bg=BG); r12c.pack(fill="x", padx=4, pady=2)
        tip(check(r12c, "Analytic Signal", self.v_as_),
            "Analytic Signal (AS)\nCombines horizontal and vertical derivatives to highlight\n"
            "anomaly edges and amplitude variations, independent of direction").pack(side="left")

        # ---- Buffer + Apply ----
        ba = tk.Frame(p, bg=BG)
        ba.pack(fill="x", padx=8, pady=8)
        lbl(ba, "FFT Buffer (px):").pack(side="left")
        self.v_buf = tk.StringVar(value="5000")
        tip(entry(ba, self.v_buf, 5),
            "Maximum buffer applied to grid to reduce edge effects (vignetting)\n"
            "Increase for noisy or clipped-edge grids").pack(side="left", padx=4)
        b_fft = btn(ba, "Apply FFT Processing", self._apply_fft, w=22)
        self._apply_btns[0] = b_fft
        tip(b_fft, "Apply selected processing steps to the selected grid").pack(side="right", padx=4)
        tip(check(ba, "Preview", self.v_preview),
            "Show the selected filter in a separate preview window that updates as you change\n"
            "settings (zoom and pan to see detail). Only one filter can be ticked while\n"
            "previewing. Apply Preview then saves the previewed filter at full resolution").pack(side="right", padx=10)

    # ======================================================================
    # TAB 2 – CONV + STATS
    # ======================================================================
    def _tab_conv_stats(self, p):
        cf = section(p, "Convolution Filters")
        cf.pack(fill="x", padx=8, pady=4)

        def crow(par, lbl_txt, v_chk, v_size, size_lbl="Size (px, odd):"):
            r = tk.Frame(par, bg=BG); r.pack(fill="x", padx=4, pady=2)
            check(r, lbl_txt, v_chk).pack(side="left")
            lbl(r, size_lbl).pack(side="left", padx=(8,0))
            entry(r, v_size, 5, chk_var=v_chk).pack(side="left", padx=2)

        self.v_mean,  self.v_mean_n  = tk.BooleanVar(), tk.StringVar(value="3")
        self.v_med,   self.v_med_n   = tk.BooleanVar(), tk.StringVar(value="3")
        self.v_gauss, self.v_gauss_s = tk.BooleanVar(), tk.StringVar(value="1.0")
        self._setup_odd_trace(self.v_mean_n)
        self._setup_odd_trace(self.v_med_n)
        # v_gauss_s is a continuous sigma value — no odd constraint
        crow(cf, "Mean",     self.v_mean,  self.v_mean_n)
        crow(cf, "Median",   self.v_med,   self.v_med_n)
        crow(cf, "Gaussian", self.v_gauss, self.v_gauss_s, "Sigma:")

        rd = tk.Frame(cf, bg=BG); rd.pack(fill="x", padx=4, pady=2)
        self.v_dconv     = tk.BooleanVar()
        self.v_dconv_dir = tk.StringVar(value="N")
        self.v_dconv_n   = tk.StringVar(value="3")
        self._setup_odd_trace(self.v_dconv_n)
        tip(check(rd, "Directional", self.v_dconv),
            "Directional enhancement\nHighlights high-frequency data in a particular direction").pack(side="left")
        combo(rd, ["N","NE","E","SE","S","SW","W","NW"],
              self.v_dconv_dir, w=5, chk_var=self.v_dconv).pack(side="left", padx=4)
        lbl(rd, "Size (odd):").pack(side="left")
        entry(rd, self.v_dconv_n, 4, chk_var=self.v_dconv).pack(side="left", padx=2)

        rs = tk.Frame(cf, bg=BG); rs.pack(fill="x", padx=4, pady=2)
        self.v_sun       = tk.BooleanVar()
        self.v_sun_alt   = tk.StringVar(value="45.0")
        self.v_sun_az    = tk.StringVar(value="315.0")
        self.v_sun_relief= tk.BooleanVar()
        tip(check(rs, "Sun Shading", self.v_sun),
            "Sun shading of grid\nUses azimuth and altitude angles to create a shaded relief effect").pack(side="left")
        lbl(rs, "Alt:").pack(side="left", padx=(8,0))
        tip(entry(rs, self.v_sun_alt, 5, chk_var=self.v_sun),
            "Altitude angle for sun shading\n90° = directly overhead, 0° = horizon").pack(side="left", padx=2)
        lbl(rs, "Az:").pack(side="left")
        tip(entry(rs, self.v_sun_az, 5, chk_var=self.v_sun),
            "Azimuth angle for sun shading\n0 = North, 90 = East, 180 = South, 270 = West").pack(side="left", padx=2)
        tip(check(rs, "relief", self.v_sun_relief),
            "Uses Grass-like shading algorithm for softer shading").pack(side="left", padx=(8,0))

        ss = section(p, "Spatial Statistics")
        ss.pack(fill="x", padx=8, pady=4)

        rw = tk.Frame(ss, bg=BG); rw.pack(fill="x", padx=4, pady=2)
        lbl(rw, "Window (px, odd):").pack(side="left")
        self.v_ss_win = tk.StringVar(value="5")
        self._setup_odd_trace(self.v_ss_win)
        entry(rw, self.v_ss_win, 5).pack(side="left", padx=4)

        self.v_ss_min  = tk.BooleanVar()
        self.v_ss_max  = tk.BooleanVar()
        self.v_ss_var  = tk.BooleanVar()
        self.v_ss_std  = tk.BooleanVar()
        self.v_ss_skew = tk.BooleanVar()
        self.v_ss_kurt = tk.BooleanVar()
        ra = tk.Frame(ss, bg=BG); ra.pack(fill="x", padx=4, pady=2)
        rb = tk.Frame(ss, bg=BG); rb.pack(fill="x", padx=4, pady=2)
        _ss_tips = {
            "Min":      "Calculate minimum of values around central pixel",
            "Max":      "Calculate maximum of values around central pixel",
            "Variance": "Calculate variance of values around central pixel\nWill be slow for larger grids/windows",
            "StdDev":   "Calculate standard deviation of values around central pixel\nWill be slow for larger grids/windows",
            "Skewness": "Calculate skewness of values around central pixel\nWill be VERY slow for larger grids/windows",
            "Kurtosis": "Calculate kurtosis of values around central pixel\nWill be VERY slow for larger grids/windows",
        }
        for i, (txt, v) in enumerate([("Min", self.v_ss_min), ("Max", self.v_ss_max),
                                       ("Variance", self.v_ss_var), ("StdDev", self.v_ss_std),
                                       ("Skewness", self.v_ss_skew), ("Kurtosis", self.v_ss_kurt)]):
            tip(check(ra if i < 3 else rb, txt, v), _ss_tips[txt]).pack(side="left", padx=(0,10))

        # ---- MRVBF / MRRTF / slope (needs a projected DEM, cell size in metres) ----
        rmv1 = tk.Frame(ss, bg=BG); rmv1.pack(fill="x", padx=4, pady=2)
        rmv2 = tk.Frame(ss, bg=BG); rmv2.pack(fill="x", padx=4, pady=(0, 2))
        self.v_mrvbf = tk.BooleanVar()
        d_ = MRVBF_DEFAULTS
        self.v_mrvbf_ts = tk.StringVar(value=f"{d_['t_slope']:g}")
        self.v_mrvbf_pv = tk.StringVar(value=f"{d_['t_pctl_v']:.2f}")
        self.v_mrvbf_pr = tk.StringVar(value=f"{d_['t_pctl_r']:.2f}")
        self.v_mrvbf_ps = tk.StringVar(value=f"{d_['p_slope']:g}")
        self.v_mrvbf_pp = tk.StringVar(value=f"{d_['p_pctl']:g}")
        self.v_mrvbf_mr = tk.StringVar(value=f"{d_['max_res']:g}")
        tip(check(rmv1, "MRVBF / MRRTF / Slope", self.v_mrvbf),
            "Multiresolution Valley Bottom Flatness (MRVBF) and Ridge Top Flatness (MRRTF)\n"
            "after Gallant & Dowling 2003. Writes _MRVBF, _MRRTF and _slope grids plus a\n"
            "3-band _MRVBF_RGB composite (R=MRRTF, G=MRVBF, B=slope).\n"
            "Needs a projected DEM with cell size in metres").pack(side="left")
        for frame, label, var, text in (
                (rmv1, "Slope threshold:", self.v_mrvbf_ts, "Initial threshold for slope (percent) at the finest resolution"),
                (rmv1, "Pctl lowness:", self.v_mrvbf_pv, "Threshold for lowness percentile (valley bottom flatness)"),
                (rmv1, "Pctl upness:", self.v_mrvbf_pr, "Threshold for upness percentile (ridge top flatness)"),
                (rmv2, "Slope shape:", self.v_mrvbf_ps, "Shape parameter for the slope transformation"),
                (rmv2, "Pctl shape:", self.v_mrvbf_pp, "Shape parameter for the elevation percentile transformation"),
                (rmv2, "Max res %:", self.v_mrvbf_mr, "Maximum resolution as a percentage of the grid diagonal")):
            lbl(frame, label).pack(side="left", padx=(8, 0))
            tip(entry(frame, var, 6, chk_var=self.v_mrvbf), text).pack(side="left", padx=2)

        # ---- Local Anisotropy / Chain Length / Streamline Length ----
        ra1 = tk.Frame(ss, bg=BG); ra1.pack(fill="x", padx=4, pady=2)
        self.v_ss_aniso          = tk.BooleanVar()
        self.v_ss_aniso_win_type = tk.StringVar(value="Gaussian")
        tip(check(ra1, "Local Anisotropy", self.v_ss_aniso),
            "Compute anisotropy magnitude (0–1) and dominant orientation (0–180°)\n"
            "via the structure tensor of Sobel gradients averaged over the window.\n"
            "Returns two layers: _SS_AnisoMag and _SS_AnisoOrient").pack(side="left")
        tip(combo(ra1, ["Gaussian", "Box"], self.v_ss_aniso_win_type, w=9, chk_var=self.v_ss_aniso),
            "Window type for smoothing the structure tensor\n"
            "Gaussian: weighted by distance (recommended)\n"
            "Box: uniform mean").pack(side="left", padx=6)

        ra2 = tk.Frame(ss, bg=BG); ra2.pack(fill="x", padx=4, pady=2)
        self.v_ss_chain  = tk.BooleanVar()
        self.v_ss_stream = tk.BooleanVar()
        tip(check(ra2, "Chain Length", self.v_ss_chain),
            "Score each active pixel by the size of the connected anisotropy component\n"
            "it belongs to.  Two active pixels link when they are within Search r pixels\n"
            "and their orientations differ by at most Angle tol°.\n"
            "Returns: _SS_ChainLen").pack(side="left", padx=(0, 10))
        tip(check(ra2, "Streamline Length", self.v_ss_stream),
            "From each active pixel trace forward and backward along the orientation field\n"
            "using sub-pixel bilinear interpolation.  Total path length is the score.\n"
            "Returns: _SS_StreamLen").pack(side="left")

        ra3 = tk.Frame(ss, bg=BG); ra3.pack(fill="x", padx=4, pady=2)
        self.v_ss_aniso_thresh = tk.StringVar(value="0.3")
        self.v_ss_angle_tol    = tk.StringVar(value="22.5")
        self.v_ss_search_r     = tk.StringVar(value="2")
        self.v_ss_max_steps    = tk.StringVar(value="200")
        lbl(ra3, "Aniso thresh:").pack(side="left")
        tip(entry(ra3, self.v_ss_aniso_thresh, 5),
            "Minimum anisotropy (0–1) for a pixel to be treated as part of a feature\n"
            "0.3 is a reasonable starting point").pack(side="left", padx=2)
        lbl(ra3, "Angle tol°:").pack(side="left", padx=(8, 0))
        tip(entry(ra3, self.v_ss_angle_tol, 5),
            "Maximum orientation difference (degrees) allowed between linked pixels\n"
            "22.5° = half of the 45° diagonal sector").pack(side="left", padx=2)
        lbl(ra3, "Search r:").pack(side="left", padx=(8, 0))
        tip(entry(ra3, self.v_ss_search_r, 4),
            "Disk search radius in pixels for Chain Length connectivity\n"
            "Larger values connect across wider features and small gaps").pack(side="left", padx=2)
        lbl(ra3, "Max steps:").pack(side="left", padx=(8, 0))
        tip(entry(ra3, self.v_ss_max_steps, 5),
            "Maximum trace steps for Streamline Length\n"
            "Each step is 0.5 px; 200 steps ≈ 100 px maximum path length").pack(side="left", padx=2)

        mv = section(p, "Multivariate Analysis")
        mv.pack(fill="x", padx=8, pady=4)

        rp = tk.Frame(mv, bg=BG); rp.pack(fill="x", padx=4, pady=2)
        self.v_pca   = tk.BooleanVar()
        self.v_pca_n = tk.StringVar(value="0")
        tip(check(rp, "PCA", self.v_pca),
            "Principal Component Analysis (PCA)\nOnly works on multiband grids\n"
            "Reduces dimensionality while preserving variance").pack(side="left")
        lbl(rp, "Components (0=all):").pack(side="left", padx=(8,0))
        tip(entry(rp, self.v_pca_n, 4, chk_var=self.v_pca),
            "Number of components to keep after PCA\nSet to 0 to keep all components").pack(side="left", padx=2)

        ri = tk.Frame(mv, bg=BG); ri.pack(fill="x", padx=4, pady=2)
        self.v_ica   = tk.BooleanVar()
        self.v_ica_n = tk.StringVar(value="0")
        tip(check(ri, "ICA", self.v_ica),
            "Independent Component Analysis (ICA)\nOnly works on multiband grids\n"
            "Separates a multivariate signal into additive, independent components").pack(side="left")
        lbl(ri, "Components (0=all):").pack(side="left", padx=(8,0))
        tip(entry(ri, self.v_ica_n, 4, chk_var=self.v_ica),
            "Number of components to keep after ICA\nSet to 0 to keep all components").pack(side="left", padx=2)

        eu = section(p, "Euler Deconvolution")
        eu.pack(fill="x", padx=8, pady=4)

        re1 = tk.Frame(eu, bg=BG); re1.pack(fill="x", padx=4, pady=2)
        self.v_euler   = tk.BooleanVar()
        self.v_eu_si   = tk.StringVar(value="1.0")
        self.v_eu_win  = tk.StringVar(value="5")
        self.v_eu_filt = tk.StringVar(value="0.1")
        tip(check(re1, "Euler Deconvolution", self.v_euler),
            "Euler deconvolution depth estimation\nEstimates source depths from potential field data").pack(side="left")
        lbl(re1, "SI:").pack(side="left", padx=(8,0))
        tip(combo(re1, ["0.0","1.0","2.0","3.0"], self.v_eu_si, w=5, chk_var=self.v_euler),
            "Structural Index (0=sphere, 1=pipe, 2=dyke, 3=fault)").pack(side="left", padx=2)
        lbl(re1, "Win:").pack(side="left", padx=(6,0))
        tip(entry(re1, self.v_eu_win, 4, chk_var=self.v_euler),
            "Window size in pixels for Euler deconvolution").pack(side="left", padx=2)
        lbl(re1, "Thresh:").pack(side="left")
        tip(entry(re1, self.v_eu_filt, 5, chk_var=self.v_euler),
            "Percentage uncertainty threshold for depth estimates (0–1)").pack(side="left", padx=2)

        ba = tk.Frame(p, bg=BG); ba.pack(fill="x", padx=8, pady=8)
        b_conv = btn(ba, "Apply Conv + Stats", self._apply_conv_stats, w=22)
        self._apply_btns[1] = b_conv
        tip(b_conv, "Apply selected processing steps to the selected grid").pack(side="right", padx=4)
        tip(check(ba, "Preview", self.v_preview),
            "Show the selected filter in a separate preview window that updates as you change\n"
            "settings (zoom and pan to see detail). Only one filter can be ticked while\n"
            "previewing. Apply Preview then saves the previewed filter at full resolution").pack(side="right", padx=10)
        self._update_apply_text()

    # ======================================================================
    # TAB 3 – GRID + WAVELETS
    # ======================================================================
    def _tab_grid_wavelets(self, p):
        im = section(p, "Import Point / Line Data")
        im.pack(fill="x", padx=8, pady=4)

        ri1 = tk.Frame(im, bg=BG); ri1.pack(fill="x", padx=4, pady=2)
        lbl(ri1, "Input file (CSV/DAT/XYZ):").pack(side="left")
        self.v_import_file = tk.StringVar()
        tip(entry(ri1, self.v_import_file, 34),
            "Select CSV, DAT or XYZ format points file").pack(side="left", padx=2)
        tip(btn(ri1, "Browse…", self._browse_import_file, w=8),
            "Select CSV, DAT (ASEG-GDF2 or legacy DAT) or XYZ format points file").pack(side="left")

        ri2 = tk.Frame(im, bg=BG); ri2.pack(fill="x", padx=4, pady=2)
        lbl(ri2, "X col:").pack(side="left")
        self.v_x_col = tk.StringVar(value="")
        self._xcol_combo = ttk.Combobox(ri2, textvariable=self.v_x_col,
                                        width=12, font=FONT)
        self._xcol_combo.pack(side="left", padx=2)
        ToolTip(self._xcol_combo, "Define X coordinate column (for csv & dat files)")
        lbl(ri2, "Y col:").pack(side="left", padx=(8,0))
        self.v_y_col = tk.StringVar(value="")
        self._ycol_combo = ttk.Combobox(ri2, textvariable=self.v_y_col,
                                        width=12, font=FONT)
        self._ycol_combo.pack(side="left", padx=2)
        ToolTip(self._ycol_combo, "Define Y coordinate column (for csv & dat files)")
        lbl(ri2, "EPSG:").pack(side="left", padx=(8,0))
        self.v_import_epsg = tk.StringVar(value="4326")
        _common_epsg = ["4326","4283","7844",
                        "28348","28349","28350","28351","28352","28353","28354","28355","28356",
                        "7848","7849","7850","7851","7852","7853","7854","7855","7856",
                        "32754","32755","32756","32757","32758","3857"]
        _epsg_cb = ttk.Combobox(ri2, values=_common_epsg,
                                textvariable=self.v_import_epsg, width=8, font=FONT)
        _epsg_cb.pack(side="left", padx=2)
        ToolTip(_epsg_cb, "Define Coordinate System of point data as EPSG code\n"
                          "(4326=WGS84, 4283=GDA94, 7844=GDA2020, 28348-56=MGA94 zones)")

        ri3 = tk.Frame(im, bg=BG); ri3.pack(fill="x", padx=4, pady=2)
        tip(btn(ri3, "Import Data", self._import_data, w=12),
            "Load points file and convert to shapefile layer").pack(side="left")

        bs = section(p, "Multilevel B-Spline Gridding")
        bs.pack(fill="x", padx=8, pady=4)
        rb1 = tk.Frame(bs, bg=BG); rb1.pack(fill="x", padx=4, pady=2)
        lbl(rb1, "Points file:").pack(side="left")
        self.v_bs_file = tk.StringVar()
        tip(entry(rb1, self.v_bs_file, 40),
            "Point file to grid (shapefile or GeoPackage).\n"
            "Filled in after Import Data above").pack(side="left", padx=2)
        tip(btn(rb1, "Browse…", self._browse_bspline_points, w=8),
            "Select a point shapefile / GeoPackage").pack(side="left")
        rb2 = tk.Frame(bs, bg=BG); rb2.pack(fill="x", padx=4, pady=2)
        lbl(rb2, "Value field:").pack(side="left")
        self.v_bs_field = tk.StringVar()
        self._bs_field_combo = ttk.Combobox(rb2, textvariable=self.v_bs_field, width=16, font=FONT)
        self._bs_field_combo.pack(side="left", padx=2)
        ToolTip(self._bs_field_combo, "Data field to grid (numbers, or text that holds numbers)")
        lbl(rb2, "Cell size:").pack(side="left", padx=(8, 0))
        self.v_bs_cell = tk.StringVar(value="100.0")
        tip(entry(rb2, self.v_bs_cell, 8), "Output cell size (map units)").pack(side="left", padx=2)
        rb3 = tk.Frame(bs, bg=BG); rb3.pack(fill="x", padx=4, pady=2)
        lbl(rb3, "Threshold error:").pack(side="left")
        self.v_bs_eps = tk.StringVar(value="0.0001")
        tip(entry(rb3, self.v_bs_eps, 8),
            "MBA threshold error (EPSILON): stop refining once every point residual is below this value").pack(side="left", padx=2)
        lbl(rb3, "Max levels:").pack(side="left", padx=(8, 0))
        self.v_bs_levels = tk.StringVar(value="11")
        tip(entry(rb3, self.v_bs_levels, 4),
            "MBA maximum number of refinement levels (LEVEL_MAX)").pack(side="left", padx=2)
        self.v_bs_auto = tk.BooleanVar(value=False)
        self._bs_auto_busy = False
        tip(check(rb2, "Auto", self.v_bs_auto, command=self._bs_auto_toggled),
            "Choose the cell size and Max levels from the spacing of the points.\n"
            "Typing either value yourself switches Auto off").pack(side="left", padx=(12, 0))
        self.v_bs_cell.trace_add("write", self._bs_user_edit)
        self.v_bs_levels.trace_add("write", self._bs_user_edit)
        self.v_bs_auto_info = tk.StringVar(value="")
        tk.Label(bs, textvariable=self.v_bs_auto_info, bg=BG, fg="#555555", font=FONT,
                 wraplength=620, justify="left", anchor="w").pack(fill="x", padx=8)
        self.v_bs_ignore = tk.BooleanVar(value=False)
        self.v_bs_ignore_val = tk.StringVar(value="-99999")
        tip(check(rb3, "Ignore values <", self.v_bs_ignore),
            "Drop points whose data value is below the threshold before gridding\n"
            "(useful for removing nodata sentinels such as -99999 or -999999.9)").pack(side="left", padx=(12, 0))
        entry(rb3, self.v_bs_ignore_val, 8, chk_var=self.v_bs_ignore).pack(side="left", padx=2)
        rb4 = tk.Frame(bs, bg=BG); rb4.pack(fill="x", padx=4, pady=2)
        tip(btn(rb4, "Multilevel B-Spline Gridding", self._grid_bspline, w=28),
            "Grid the points with the multilevel B-spline method (native translation of SAGA's\n"
            "grid_spline tool). Writes <points>_<field>_bspline.tif beside the points file").pack(side="left")

        wr = section(p, "BSD Worms")
        wr.pack(fill="x", padx=8, pady=4)

        rw1 = tk.Frame(wr, bg=BG); rw1.pack(fill="x", padx=4, pady=2)
        lbl(rw1, "# Levels:").pack(side="left")
        self.v_worm_levels = tk.StringVar(value="10")
        entry(rw1, self.v_worm_levels, 5).pack(side="left", padx=2)
        lbl(rw1, "Base level:").pack(side="left", padx=(8,0))
        self.v_worm_base = tk.StringVar(value="1000.0")
        entry(rw1, self.v_worm_base, 9).pack(side="left", padx=2)
        lbl(rw1, "Increment:").pack(side="left")
        self.v_worm_inc = tk.StringVar(value="1000.0")
        entry(rw1, self.v_worm_inc, 9).pack(side="left", padx=2)

        rw2 = tk.Frame(wr, bg=BG); rw2.pack(fill="x", padx=4, pady=2)
        self.v_worm_shp = tk.BooleanVar(value=False)
        tip(check(rw2, "Also save shapefiles", self.v_worm_shp),
            "Convert worms to polyline shapefile\n(Can be slow; best to start worms at ≥2000 m)").pack(side="left")
        tip(btn(rw2, "Calculate Worms", self._calc_worms, w=16),
            "Create CSV file of worms using BSDWormer code").pack(side="right", padx=4)


    # ======================================================================
    # TAB 4 – PETROPHYS
    # ======================================================================
    APPARENT_DENSITY = ("Apparent density (g/cc)", False)    # axis label, log axis
    APPARENT_SUS = ("Apparent susceptibility (SI)", True)
    VECTOR_FILES = [("Vector files", "*.shp *.gpkg *.geojson *.json *.gml *.kml *.tab"),
                    ("All files", "*.*")]

    def _tab_petro(self, p):
        def row(parent, indent=0):
            r = tk.Frame(parent, bg=BG)
            r.pack(fill="x", padx=(4 + indent, 4), pady=2)
            return r

        def picks(kind, *entries):
            """Typing in an entry ticks its calculation."""
            for e in entries:
                e.bind("<FocusIn>", lambda _e: self._app_tick(kind), add=True)
                e.bind("<Key>", lambda _e: self._app_tick(kind), add=True)
            return entries

        cut_tip = ("High-cut wavelength (m). Dividing out the layer response amplifies short\n"
                   "wavelengths, so they are removed with a cosine roll-off: wavelengths longer\n"
                   "than twice this value are kept, those at this value or shorter removed.\n"
                   "Leave blank for 4 x the cell size; 0 for no high-cut.")

        # ---- Apparent Petrophysics ----
        ap = section(p, "Apparent Petrophysics")
        ap.pack(fill="x", padx=8, pady=4)
        # at most one calculation ticked: ticking one unticks the other (both may be
        # unticked, so arriving on this tab with Preview on calculates nothing)
        self.v_app_sus = tk.BooleanVar(value=False)
        self.v_app_dens = tk.BooleanVar(value=False)
        self.v_app_sus.trace_add("write", lambda *_: self.v_app_sus.get() and self.v_app_dens.set(False))
        self.v_app_dens.trace_add("write", lambda *_: self.v_app_dens.get() and self.v_app_sus.set(False))

        r = row(ap)
        tip(check(r, "Apparent Susceptibility", self.v_app_sus),
            "Apparent susceptibility (SI) from a magnetic anomaly grid (nT, IGRF removed) that\n"
            "has already been reduced to the pole (RTP), or to the equator (RTE: the result is\n"
            "then only approximate). The ground is modelled as a layer of vertical prisms\n"
            "magnetised by induction in the ambient field, and the layer's response is divided\n"
            "out in the Fourier domain. The grid's absolute level is not used: the result is a\n"
            "variation about the background susceptibility. Output: <grid>_AppSus.").pack(side="left")
        r = row(ap, 22)
        self.v_app_field = tk.StringVar(value="50000")
        self.v_app_year = tk.StringVar(value=self.v_igrf_year.get())
        lbl(r, "Field (nT):").pack(side="left")
        e_field = tip(entry(r, self.v_app_field, 8),
                      "Ambient field strength (nT): the full IGRF total intensity, not anything\n"
                      "from the anomaly grid. Filled in from the IGRF at the centre of the input\n"
                      "grid in the year shown whenever the input changes or IGRF is pressed (or\n"
                      "from a Noddy grid's own metadata). It can also be typed in.")
        e_field.pack(side="left", padx=2)
        tip(btn(r, "IGRF", self._petro_igrf, w=6),
            "Set the field strength from the IGRF at the centre of the input grid").pack(side="left", padx=4)
        lbl(r, "Year:").pack(side="left", padx=(6, 0))
        e_year = tip(entry(r, self.v_app_year, 7), "Survey year (decimal) for the IGRF")
        e_year.pack(side="left", padx=2)
        r = row(ap, 22)
        self.v_app_sus_depth = tk.StringVar(value="0")
        self.v_app_sus_finite = tk.BooleanVar(value=False)
        self.v_app_sus_thick = tk.StringVar(value="1000")
        lbl(r, "Depth to top (m):").pack(side="left")
        e_sd = tip(entry(r, self.v_app_sus_depth, 7),
                   "Depth from the sensor to the top of the magnetic layer (m).\n"
                   "0 puts the top of the sources at the sensor. For airborne data set it to\n"
                   "the terrain clearance (plus any cover) so the sources sit at the ground.")
        e_sd.pack(side="left", padx=2)
        tip(check(r, "Thickness (m):", self.v_app_sus_finite,
                  command=lambda: self._app_tick("sus")),
            "Untick for a layer of infinite thickness (the usual choice); tick to give\n"
            "the layer a finite thickness").pack(side="left", padx=(8, 0))
        e_st = tip(entry(r, self.v_app_sus_thick, 7), "Layer thickness (m), used when Thickness is ticked")
        e_st.pack(side="left", padx=2)
        r = row(ap, 22)
        self.v_app_sus_bg = tk.StringVar(value="0.001")
        self.v_app_sus_cut = tk.StringVar(value="")
        lbl(r, "Background (SI):").pack(side="left")
        e_sb = tip(entry(r, self.v_app_sus_bg, 7),
                   "Background susceptibility (SI) added to the result. A magnetic anomaly grid\n"
                   "has no absolute level, so the calculated susceptibility is relative to an\n"
                   "unknown background and is partly negative. A typical background (e.g. 0.001\n"
                   "SI) keeps most values positive; Henkel plots and histograms leave out values\n"
                   "<= 0, which a log axis can't show.")
        e_sb.pack(side="left", padx=2)
        lbl(r, "High-cut (m):").pack(side="left", padx=(8, 0))
        e_sc = tip(entry(r, self.v_app_sus_cut, 8), cut_tip)
        e_sc.pack(side="left", padx=2)
        tk.Label(r, text="(blank = 4 x cell)", bg=BG, fg="#666666", font=FONT).pack(side="left")
        picks("sus", e_field, e_year, e_sd, e_st, e_sb, e_sc)

        r = row(ap)
        tip(check(r, "Apparent Density", self.v_app_dens),
            "Apparent density (g/cc) from a Bouguer gravity grid. The ground is modelled as a\n"
            "layer of vertical prisms of finite thickness, and the layer's response is divided\n"
            "out in the Fourier domain. A Bouguer grid's absolute level depends on its datum\n"
            "and the regional field, so it is not used: the result's mean is the background\n"
            "density. Output: <grid>_AppDens.").pack(side="left")
        self.v_app_units = tk.StringVar(value="mGal")
        lbl(r, "Grid units:").pack(side="left", padx=(12, 0))
        c_units = tip(combo(r, ["mGal", "um/s2 (gu)"], self.v_app_units, w=10),
                      "Units of the gravity grid (1 mGal = 10 um/s2 = 10 gu). Geoscience Australia's\n"
                      "national grids are in um/s2.")
        c_units.pack(side="left", padx=4)
        c_units.bind("<<ComboboxSelected>>", lambda _e: self._app_tick("dens"), add=True)
        r = row(ap, 22)
        self.v_app_dens_depth = tk.StringVar(value="0")
        self.v_app_dens_thick = tk.StringVar(value="5000")
        lbl(r, "Depth to top (m):").pack(side="left")
        e_dd = tip(entry(r, self.v_app_dens_depth, 7),
                   "Depth from the sensor to the top of the layer (m).\n"
                   "0 for ground gravity; the terrain clearance for airborne gravity.")
        e_dd.pack(side="left", padx=2)
        lbl(r, "Thickness (m):").pack(side="left", padx=(8, 0))
        e_dt = tip(entry(r, self.v_app_dens_thick, 7),
                   "Layer thickness (m); must be greater than 0. Thin layers exaggerate the\n"
                   "density contrasts needed to explain long-wavelength anomalies.")
        e_dt.pack(side="left", padx=2)
        r = row(ap, 22)
        self.v_app_dens_bg = tk.StringVar(value="2.67")
        self.v_app_dens_cut = tk.StringVar(value="")
        lbl(r, "Background (g/cc):").pack(side="left")
        e_db = tip(entry(r, self.v_app_dens_bg, 7), "Background density (g/cc) added to the density contrast")
        e_db.pack(side="left", padx=2)
        lbl(r, "High-cut (m):").pack(side="left", padx=(8, 0))
        e_dc = tip(entry(r, self.v_app_dens_cut, 8), cut_tip)
        e_dc.pack(side="left", padx=2)
        tk.Label(r, text="(blank = 4 x cell)", bg=BG, fg="#666666", font=FONT).pack(side="left")
        picks("dens", e_dd, e_dt, e_db, e_dc)

        r = row(ap)
        b_app = btn(r, "Apply Petrophysics", self._apply_petro, w=22)
        self._apply_btns[3] = b_app
        tip(b_app, "Calculate the chosen grid (apparent susceptibility or density) from the\n"
                   "input grid, using the FFT buffer set on the FFT Filters tab").pack(side="right", padx=4)
        tip(check(r, "Preview", self.v_preview),
            "Show the chosen calculation in a separate preview window that updates as you\n"
            "change settings (zoom and pan to see detail). Apply Preview then saves it at\n"
            "full resolution").pack(side="right", padx=10)
        self._update_apply_text()

        # ---- Apparent Henkel Plot / Histogram ----
        hk = section(p, "Apparent Henkel Plot / Histogram")
        hk.pack(fill="x", padx=8, pady=4)
        self.v_hk_sus = tk.StringVar()
        self.v_hk_dens = tk.StringVar()
        self._petro_raster_combos = []
        for text, var, tiptext in (
                ("App. Susceptibility:", self.v_hk_sus,
                 "Apparent susceptibility grid (SI): a layer of the map, or a file. Leave\n"
                 "blank, with a density grid chosen, for density histograms."),
                ("App. Density:", self.v_hk_dens,
                 "Apparent density grid (g/cc): a layer of the map, or a file. Leave blank,\n"
                 "with a susceptibility grid chosen, for susceptibility histograms.")):
            r = row(hk)
            lbl(r, text, width=17, anchor="e").pack(side="left")
            c = ttk.Combobox(r, textvariable=var, width=44, font=FONT)
            c.pack(side="left", padx=2)
            ToolTip(c, tiptext)
            self._petro_raster_combos.append(c)
            btn(r, "Browse…", lambda v=var: self._petro_browse_raster(v), w=8).pack(side="left", padx=2)
        self._petro_refresh_rasters()
        r = row(hk)
        tip(btn(r, "Plot Apparent", self._plot_apparent, w=16),
            "Heat map of apparent density (linear x axis) against apparent susceptibility\n"
            "(log y axis), counting grid cells per bin on a log colour scale. Each cell of\n"
            "the coarser grid (usually gravity) is paired with the nearest cell of the finer\n"
            "grid, reprojecting if the CRSs differ; cells where either grid has no data are\n"
            "left out, as are susceptibilities <= 0 (counted in the title). Axes span the\n"
            "0.1-99.9 percentiles, to keep FFT edge effects from squashing the plot.\n"
            "With only one grid chosen: a histogram of it instead.").pack(side="left")
        ttk.Separator(hk, orient="horizontal").pack(fill="x", padx=4, pady=4)
        self.v_hk_poly = tk.StringVar()
        self.v_hk_poly_field = tk.StringVar()
        r = row(hk)
        lbl(r, "Polygons:", width=17, anchor="e").pack(side="left")
        tip(entry(r, self.v_hk_poly, 46),
            "Polygon file, e.g. a geology map (shapefile, GeoPackage ...)").pack(side="left", padx=2)
        btn(r, "Browse…", lambda: self._petro_browse_vector(
            self.v_hk_poly, [(self._hk_poly_field_combo, False)]), w=8).pack(side="left", padx=2)
        r = row(hk)
        lbl(r, "Category field:", width=17, anchor="e").pack(side="left")
        self._hk_poly_field_combo = tip(combo(r, [], self.v_hk_poly_field, w=20),
                                        "Field holding the category (e.g. lithology) of each polygon")
        self._hk_poly_field_combo.pack(side="left", padx=2)
        r = row(hk)
        self.v_hk_position = tk.BooleanVar(value=False)
        tip(btn(r, "Plot per Category", self._plot_categories, w=16),
            "A grid of Henkel heat maps (or histograms, with only one grid chosen), one per\n"
            "category, using the grids above clipped to the polygons of that category.\n"
            "Categories with no cells where the grids have data are left out. All panels\n"
            "share their axes and colour scale; each shows the fraction of its own category's\n"
            "cells per bin, so small and large categories compare fairly.").pack(side="left")
        tip(check(r, "Colour by position", self.v_hk_position),
            "For Henkel plots per category: colour each bin by where it is on the plot, red\n"
            "to green along the density axis and red to blue up the susceptibility axis (low/\n"
            "low red, dense green, magnetic blue, both cyan), with the fraction of the\n"
            "category's cells per bin shown by opacity instead of a colour scale.").pack(side="left", padx=10)

        # ---- Measured Petrophysics ----
        ms = section(p, "Measured Petrophysics")
        ms.pack(fill="x", padx=8, pady=4)
        self.v_ms_points = tk.StringVar()
        points_row = row(ms)
        lbl(points_row, "Points:", width=17, anchor="e").pack(side="left")
        tip(entry(points_row, self.v_ms_points, 46),
            "Point file of samples with measured properties, e.g. density and susceptibility\n"
            "(shapefile, GeoPackage ...)").pack(side="left", padx=2)
        self._ms_field_combos = []
        field_tip = ("Field plotted on the {axis} axis, and named on it: numeric, or text holding\n"
                     "numbers (values that aren't numbers are skipped). Leave blank, with the\n"
                     "{other} field chosen, for histograms of that field.")
        scale_tip = ("Multiplier applied to the field's values, e.g. 0.001 for density in kg/m3\n"
                     "to g/cc, or 0.00001 for susceptibility in x10^-5 SI to SI")
        log_tip = ("Log scale for the {axis} axis (and log bins in histograms). Values <= 0\n"
                   "can't be shown on a log axis, so they are left out (their number is given\n"
                   "in the title).")
        for axis, other, log_default in (("x", "Y", False), ("y", "X", True)):
            field, scale, log = tk.StringVar(), tk.StringVar(value="1"), tk.BooleanVar(value=log_default)
            setattr(self, f"v_ms_{axis}", field)
            setattr(self, f"v_ms_{axis}_scale", scale)
            setattr(self, f"v_ms_{axis}_log", log)
            r = row(ms)
            lbl(r, f"{axis.upper()} field:", width=17, anchor="e").pack(side="left")
            c = tip(combo(r, [""], field, w=20), field_tip.format(axis=axis, other=other))
            c.pack(side="left", padx=2)
            self._ms_field_combos.append((c, True))
            lbl(r, "x").pack(side="left", padx=(6, 0))
            tip(entry(r, scale, 8), scale_tip).pack(side="left", padx=2)
            tip(check(r, "Log", log), log_tip.format(axis=axis)).pack(side="left", padx=6)
        self.v_ms_cat = tk.StringVar()
        r = row(ms)
        lbl(r, "Category field:", width=17, anchor="e").pack(side="left")
        c = tip(combo(r, [""], self.v_ms_cat, w=20),
                "Optional: field holding a category (e.g. lithology). With Points, each\n"
                "category gets its own symbol and a legend; with only one of the X and Y\n"
                "fields, a histogram per category.")
        c.pack(side="left", padx=2)
        self._ms_field_combos.append((c, True))
        btn(points_row, "Browse…",
            lambda: self._petro_browse_vector(self.v_ms_points, self._ms_field_combos),
            w=8).pack(side="left", padx=2)
        r = row(ms)
        self.v_ms_style = tk.StringVar(value="points")
        tip(btn(r, "Plot Samples", self._plot_samples, w=16),
            "The Y field against the X field for every sample (a Henkel plot with density\n"
            "as X and susceptibility as Y, log). The axes are named after the fields.\n"
            "With only one field chosen: a histogram of it instead, one per category if a\n"
            "category field is chosen.").pack(side="left")
        for text, value, tiptext in (
                ("Points", "points",
                 "Plot every sample, with a colour and marker per category if one is chosen:\n"
                 "8 colours x 8 markers give up to 64 categories, assigned to the categories\n"
                 "with the most samples (beyond 64, the rest are grouped as Other)."),
                ("Heat map", "heat",
                 "Samples per bin on a log colour scale, for large datasets (ignores category)")):
            tip(tk.Radiobutton(r, text=text, variable=self.v_ms_style, value=value, bg=BG, fg=FG,
                               selectcolor=ENTRY_BG, activebackground=BG, activeforeground=FG,
                               font=FONT), tiptext).pack(side="left", padx=(10, 0))

        # ---- Save Plots ----
        sp = section(p, "Save Plots")
        sp.pack(fill="x", padx=8, pady=4)
        r = row(sp)
        self.v_plot_png = tk.BooleanVar(value=False)
        self.v_plot_pdf = tk.BooleanVar(value=False)
        self.v_plot_dir = tk.StringVar()
        tip(check(r, "PNG", self.v_plot_png), "Also save each plot in this format").pack(side="left")
        tip(check(r, "PDF", self.v_plot_pdf), "Also save each plot in this format").pack(side="left")
        lbl(r, "Folder:").pack(side="left", padx=(8, 0))
        tip(entry(r, self.v_plot_dir, 36),
            "Folder for saved plots. Leave blank to save them beside the input grid or file.\n"
            "Plots can also be saved from the plot window's toolbar.").pack(side="left", padx=2)
        btn(r, "Browse…", lambda: _browse_dir(self.v_plot_dir), w=8).pack(side="left", padx=2)

        # the field strength follows the input grid
        self.v_input.trace_add("write", lambda *_: self.after(300, self._petro_igrf, True))
        self._petro_igrf(True)

    # ---- Petrophys: inputs ----
    def _petro_refresh_rasters(self):
        """List the map's raster layers (and a blank) in the apparent grid dropdowns."""
        names = [""] + list(self._layer_combo["values"])
        for c in getattr(self, "_petro_raster_combos", []):
            c["values"] = names

    def _petro_browse_raster(self, var):
        path = filedialog.askopenfilename(
            filetypes=[("Raster files", "*.tif *.tiff *.img *.ers"), ("All files", "*.*")])
        if path:
            var.set(path)

    def _petro_browse_vector(self, var, combos):
        """Choose a vector file and list its fields in combos: (combobox, blank allowed)."""
        path = filedialog.askopenfilename(filetypes=self.VECTOR_FILES)
        if not path:
            return
        var.set(path)
        fields = list_point_fields(path)
        for c, blank in combos:
            c["values"] = ([""] if blank else []) + fields

    def _petro_raster(self, var):
        """File path of an apparent grid dropdown, or None if it is blank."""
        name = var.get().strip()
        if not name:
            return None
        path = self._layer_map.get(name, name)
        if not os.path.exists(path):
            raise ValueError(f"Cannot find the grid {name}")
        return path

    def _petro_igrf(self, quiet=False):
        """Field strength for apparent susceptibility: from a Noddy grid's metadata,
        otherwise the IGRF at the centre of the input grid in the year set."""
        try:
            from arcgis_utils import calc_igrf, raster_center_latlon, get_noddy_mag_params

            in_r = self._resolve_input()
            if not in_r or not os.path.exists(in_r):
                if not quiet:
                    messagebox.showwarning("IGRF", "Select an input raster first.")
                return
            noddy = get_noddy_mag_params(in_r)
            if noddy is not None:
                field = noddy[2]
            else:
                centre = raster_center_latlon(in_r)
                res = calc_igrf(centre[0], centre[1], 0.0, float(self.v_app_year.get())) if centre else None
                if res is None:
                    if not quiet:
                        self._status("IGRF failed – check the grid's coordinate system and the year.")
                    return
                field = res[2]
            self.v_app_field.set(f"{field:.0f}")
            if not quiet:
                self._status(f"IGRF field strength {field:.0f} nT")
        except Exception as e:
            if not quiet:
                messagebox.showerror("IGRF Error", str(e))

    # ---- Petrophys: apparent susceptibility / density ----
    def _petro_values(self):
        """(spec id, values) of the chosen calculation, named as in the shared filter
        table (sgt_filters): thickness 0 = infinite, high cut 0 = 4 x cell, -1 = none."""
        def num(var, what, default=None):
            text = var.get().strip().replace(",", ".")
            if not text:
                if default is None:
                    raise ValueError(f"Set the {what}")
                return default
            try:
                return float(text)
            except ValueError:
                raise ValueError(f"The {what} must be a number")

        def cut(var):
            value = num(var, "high-cut wavelength", default=0.0) if var.get().strip() else None
            return 0.0 if value is None else (-1.0 if value == 0 else value)

        if not (self.v_app_sus.get() or self.v_app_dens.get()):
            raise ValueError("Tick Apparent Susceptibility or Apparent Density")
        if self.v_app_sus.get():
            values = dict(
                FIELD=num(self.v_app_field, "field strength (nT)"),
                DEPTH=num(self.v_app_sus_depth, "depth to top", 0.0),
                THICKNESS=(num(self.v_app_sus_thick, "thickness")
                           if self.v_app_sus_finite.get() else 0.0),
                BACKGROUND=num(self.v_app_sus_bg, "background susceptibility", 0.0),
                HIGH_CUT=cut(self.v_app_sus_cut),
            )
            if values["FIELD"] <= 0:
                raise ValueError("Set the field strength (nT) for apparent susceptibility")
            return "apparent_susceptibility", values
        values = dict(
            UNITS=self.v_app_units.get(),
            DEPTH=num(self.v_app_dens_depth, "depth to top", 0.0),
            THICKNESS=num(self.v_app_dens_thick, "thickness"),
            BACKGROUND=num(self.v_app_dens_bg, "background density", 0.0),
            HIGH_CUT=cut(self.v_app_dens_cut),
        )
        if values["THICKNESS"] <= 0:
            raise ValueError("The apparent density layer thickness must be greater than 0")
        return "apparent_density", values

    def _petro_steps(self, arr, cx, cy, buf, proc, _save):
        """The chosen apparent calculation; the result goes to _save(result, suffix, step)."""
        spec_id, v = self._petro_values()
        cutoff = 4 * max(cx, cy) if v["HIGH_CUT"] == 0 else (None if v["HIGH_CUT"] < 0 else v["HIGH_CUT"])
        if spec_id == "apparent_susceptibility":
            self.after(0, lambda: self._status("Apparent susceptibility…"))
            result = proc.apparent_susceptibility(
                arr, v["FIELD"], depth=v["DEPTH"], thickness=v["THICKNESS"] or None,
                background=v["BACKGROUND"], cutoff_wavelength=cutoff, buffer_size=buf)
            suffix = "_AppSus"
        else:
            self.after(0, lambda: self._status("Apparent density…"))
            to_mgal = 0.1 if v["UNITS"].startswith("um") else 1.0  # 1 mGal = 10 um/s2
            result = proc.apparent_density(
                arr * to_mgal, depth=v["DEPTH"], thickness=v["THICKNESS"],
                background=v["BACKGROUND"], cutoff_wavelength=cutoff, buffer_size=buf)
            suffix = "_AppDens"
        _save(result, suffix, step_record(spec_id, **v))

    def _app_tick(self, kind):
        """Tick apparent susceptibility ("sus") or density ("dens"), unticking the other."""
        (self.v_app_sus if kind == "sus" else self.v_app_dens).set(True)

    @staticmethod
    def _is_geographic(path):
        from osgeo import gdal, osr

        ds = gdal.Open(str(path))
        wkt = ds.GetProjection() if ds is not None else ""
        ds = None
        if not wkt:
            return False
        srs = osr.SpatialReference()
        srs.ImportFromWkt(wkt)
        return bool(srs.IsGeographic())

    def _apply_petro(self):
        in_r = self._compute_output_path()
        if not in_r:
            messagebox.showwarning("No Input", "Select an input raster.")
            return
        try:
            self._petro_values()
            if self._is_geographic(in_r):
                raise ValueError("Apparent susceptibility and density need a grid in a "
                                 "projected (metre-based) coordinate system")
        except ValueError as e:
            messagebox.showwarning("Petrophysics", str(e))
            return
        self._end_preview_for_apply()
        self._status("Starting processing…")
        self._run_in_thread(self._do_apply_petro, in_r)

    def _do_apply_petro(self, in_r):
        arr, nodata, ll_x, ll_y, cx, cy, sr = raster_to_numpy(in_r)
        buf = min(int(self.v_buf.get() or 10), *arr.shape)
        proc = make_processor(cx, cy, buf)
        base, ext = os.path.splitext(in_r)
        ext = ext or ".tif"
        written = []

        def _save(result, suffix, step=None):
            self._check_cancel()
            out_r = f"{base}{step_suffix(step) or suffix}{ext}"
            numpy_to_raster(result, out_r, ll_x, ll_y, cx, cy, sr, nodata)
            write_provenance(out_r, in_r, *(step or (suffix, {})))
            self._record_pending(out_r)
            written.append(out_r)
            self.after(0, self._on_output_written, out_r)

        self._petro_steps(arr, cx, cy, buf, proc, _save)
        self.after(0, lambda: self._status(f"Done — {len(written)} output(s) written."))
        self.after(0, lambda: (self.v_app_sus.set(False), self.v_app_dens.set(False)))

    # ---- Petrophys: plots ----
    def _petro_plot(self, gather, build):
        """gather() reads the settings (main thread); build(settings) makes the figure in
        a worker thread and returns (figure, file name, default folder). The figure is
        shown in its own window and saved as PNG / PDF if ticked."""
        try:
            settings = gather()
        except ValueError as e:
            messagebox.showwarning("Petrophysics", str(e))
            return
        settings["formats"] = [f for f, v in (("png", self.v_plot_png), ("pdf", self.v_plot_pdf)) if v.get()]
        settings["folder"] = self.v_plot_dir.get().strip()
        self._status("Plotting…")

        def work():
            try:
                fig, name, folder = build(settings)
            except ValueError as e:
                self.after(0, lambda m=str(e): (messagebox.showwarning("Petrophysics", m),
                                                self._status(m)))
                return
            saved = []
            if settings["formats"]:
                folder = settings["folder"] or folder
                if folder and os.path.isdir(folder):
                    base = os.path.join(folder, re.sub(r"[^\w\-]+", "_", name))
                    saved = Petrophysics.save_figure(fig, base, settings["formats"])
            self.after(0, self._show_figure, fig, name, saved)

        self._run_in_thread(work)

    def _show_figure(self, fig, name, saved):
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

        win = tk.Toplevel(self)
        win.title(name)
        w, h = fig.get_size_inches()
        win.geometry(f"{int(w * fig.dpi)}x{int(h * fig.dpi) + 40}")
        canvas = FigureCanvasTkAgg(fig, master=win)
        canvas.draw()
        toolbar = NavigationToolbar2Tk(canvas, win)
        toolbar.update()
        canvas.get_tk_widget().pack(fill="both", expand=True)
        self._status("Saved " + ", ".join(saved) if saved else "Plot displayed.")

    @staticmethod
    def _stem(path):
        return os.path.splitext(os.path.basename(path))[0]

    def _gather_apparent(self):
        dens, sus = self._petro_raster(self.v_hk_dens), self._petro_raster(self.v_hk_sus)
        if dens is None and sus is None:
            raise ValueError("Choose an apparent susceptibility and/or density grid")
        if dens is not None and sus is not None and os.path.abspath(dens) == os.path.abspath(sus):
            raise ValueError("The apparent susceptibility and density grids are the same file")
        return {"dens": dens, "sus": sus}

    def _single_grid(self, s, density, sus):
        """For a histogram of whichever grid was chosen: (values, label, log, path)."""
        if density is not None:
            return (density,) + self.APPARENT_DENSITY + (s["dens"],)
        return (sus,) + self.APPARENT_SUS + (s["sus"],)

    def _plot_apparent(self):
        def build(s):
            density, sus, _gt, _wkt = Petrophysics.paired_grids(s["dens"], s["sus"])
            if density is not None and sus is not None:
                names = f"{self._stem(s['dens'])} vs {self._stem(s['sus'])}"
                fig = Petrophysics.heatmap(
                    density.ravel(), sus.ravel(), f"Apparent Henkel plot: {names}",
                    self.APPARENT_DENSITY[0], self.APPARENT_SUS[0], what="cells")
                return (fig, f"{self._stem(s['dens'])}_{self._stem(s['sus'])}_henkel",
                        os.path.dirname(s["dens"]))
            values, label, log, path = self._single_grid(s, density, sus)
            fig = Petrophysics.histogram(values.ravel(), f"{label.split(' (')[0]}: {self._stem(path)}",
                                         label, log=log, what="cells")
            return fig, f"{self._stem(path)}_hist", os.path.dirname(path)

        self._petro_plot(self._gather_apparent, build)

    def _plot_categories(self):
        def gather():
            s = self._gather_apparent()
            s["poly"], s["field"] = self.v_hk_poly.get().strip(), self.v_hk_poly_field.get()
            s["position"] = self.v_hk_position.get()
            if not s["poly"] or not s["field"]:
                raise ValueError("Choose a polygon file and its category field")
            return s

        def build(s):
            density, sus, gt, wkt = Petrophysics.paired_grids(s["dens"], s["sus"])
            shapes = read_polygon_categories(s["poly"], s["field"], wkt)
            if not shapes:
                raise ValueError("The polygon file has no polygons with a category")
            names = sorted({name for name, _ in shapes}, key=Petrophysics._natural_key)
            index = {name: i for i, name in enumerate(names)}
            shape = (density if density is not None else sus).shape
            cats = Petrophysics.rasterize_categories(
                [(geom, index[name]) for name, geom in shapes], wkt, gt, shape).ravel()
            field, polys = s["field"], self._stem(s["poly"])
            if density is not None and sus is not None:
                fig = Petrophysics.heatmap_by_category(
                    density.ravel(), sus.ravel(), cats, names,
                    f"Apparent Henkel plots by {field} ({polys})",
                    self.APPARENT_DENSITY[0], self.APPARENT_SUS[0],
                    position_colours=s["position"])
                return (fig, f"{self._stem(s['dens'])}_{self._stem(s['sus'])}_henkel_by_{field}",
                        os.path.dirname(s["dens"]))
            values, label, log, path = self._single_grid(s, density, sus)
            fig = Petrophysics.histogram_by_category(
                values.ravel(), cats, names, f"{label.split(' (')[0]} by {field} ({polys})",
                label, log=log)
            return fig, f"{self._stem(path)}_hist_by_{field}", os.path.dirname(path)

        self._petro_plot(gather, build)

    def _plot_samples(self):
        def gather():
            s = {"points": self.v_ms_points.get().strip(), "x": self.v_ms_x.get(),
                 "y": self.v_ms_y.get(), "cat": self.v_ms_cat.get(),
                 "x_log": self.v_ms_x_log.get(), "y_log": self.v_ms_y_log.get(),
                 "heat": self.v_ms_style.get() == "heat"}
            if not s["points"] or not (s["x"] or s["y"]):
                raise ValueError("Choose a point file and its X and/or Y field")
            try:
                s["x_scale"] = float(self.v_ms_x_scale.get().replace(",", "."))
                s["y_scale"] = float(self.v_ms_y_scale.get().replace(",", "."))
            except ValueError:
                raise ValueError("The X and Y multipliers must be numbers")
            return s

        def build(s):
            x_field, y_field, cat_field = s["x"], s["y"], s["cat"]
            raw = read_attributes_ogr(s["points"], [x_field, y_field, cat_field])
            xs = (np.array([Petrophysics.to_number(v) for v in raw[x_field]]) * s["x_scale"]
                  if x_field else None)
            ys = (np.array([Petrophysics.to_number(v) for v in raw[y_field]]) * s["y_scale"]
                  if y_field else None)
            cats = [Petrophysics.category_label(v) for v in raw[cat_field]] if cat_field else None
            layer, folder = self._stem(s["points"]), os.path.dirname(s["points"])

            if not (x_field and y_field):
                values, field, log = (xs, x_field, s["x_log"]) if x_field else (ys, y_field, s["y_log"])
                if cat_field:
                    index, names = Petrophysics.index_categories(cats)
                    fig = Petrophysics.histogram_by_category(
                        values, index, names, f"{field} by {cat_field}: {layer}", field,
                        log=log, robust=False, what="samples")
                    return fig, f"{layer}_{field}_hist_by_{cat_field}", folder
                fig = Petrophysics.histogram(values, f"{field}: {layer}", field,
                                             log=log, robust=False, what="samples")
                return fig, f"{layer}_{field}_hist", folder

            title = f"{y_field} vs {x_field}: {layer}"
            if s["heat"]:
                n = int(np.count_nonzero(np.isfinite(xs) & np.isfinite(ys)))
                fig = Petrophysics.heatmap(xs, ys, title, x_field, y_field,
                                           x_log=s["x_log"], y_log=s["y_log"], robust=False,
                                           bins=Petrophysics.heat_bins(n), what="samples")
            else:
                fig = Petrophysics.scatter(xs, ys, title, x_field, y_field,
                                           x_log=s["x_log"], y_log=s["y_log"],
                                           categories=cats, legend_title=cat_field or None)
            return fig, f"{layer}_{y_field}_vs_{x_field}", folder

        self._petro_plot(gather, build)

    # ======================================================================
    # TAB 5 – UTILS
    # ======================================================================
    def _tab_utils(self, p):
        md = section(p, "SGTool Metadata")
        md.pack(fill="x", padx=8, pady=4)
        rm1 = tk.Frame(md, bg=BG); rm1.pack(fill="x", padx=4, pady=2)
        tip(btn(rm1, "Read metadata of input grid", self._show_metadata, w=26),
            "Show how the input grid was made: when, from which source file(s) and with\n"
            "which operation and parameters, back through every earlier SGTool step.\n"
            "This is stored inside GeoTIFFs saved by SGTool.").pack(side="left")
        tip(btn(rm1, "Save as XML", self._save_metadata_xml, w=12),
            "Save the input grid's SGTool metadata as <grid file>.sgt.xml\n"
            "in the same folder as the grid").pack(side="left", padx=4)
        tip(btn(rm1, "Apply same steps to another grid…", self._open_replay, w=30),
            "Repeat the processing steps in the input grid's history on another grid,\n"
            "each on the result of the one before. Reduction to the pole / equator and\n"
            "differential RTP get a new inclination and declination from IGRF").pack(side="left", padx=4)

        th = section(p, "Threshold to NaN")
        th.pack(fill="x", padx=8, pady=4)
        rt1 = tk.Frame(th, bg=BG); rt1.pack(fill="x", padx=4, pady=2)
        self.v_thresh_en    = tk.BooleanVar()
        self.v_thresh_cond  = tk.StringVar(value="above")
        self.v_thresh_above = tk.StringVar(value="9.99e9")
        self.v_thresh_below = tk.StringVar(value="-9.99e9")
        tip(check(rt1, "Enable", self.v_thresh_en),
            "Threshold background values to NaN").pack(side="left")
        tip(combo(rt1, ["above","below","between"], self.v_thresh_cond, w=8, chk_var=self.v_thresh_en),
            "Threshold background values above, below or between set values to NaN").pack(side="left", padx=4)
        lbl(rt1, "Above:").pack(side="left")
        tip(entry(rt1, self.v_thresh_above, 10, chk_var=self.v_thresh_en),
            "Upper threshold value").pack(side="left", padx=2)
        lbl(rt1, "Below:").pack(side="left")
        tip(entry(rt1, self.v_thresh_below, 10, chk_var=self.v_thresh_en),
            "Lower threshold value").pack(side="left", padx=2)
        b_thr = btn(rt1, "Apply Threshold → NaN", self._apply_threshold, w=22)
        self._apply_btns[4] = b_thr
        tip(b_thr, "Apply threshold and set matching values to NaN (no-data)").pack(side="left", padx=8)
        tip(check(rt1, "Preview", self.v_preview),
            "Show the threshold in a separate preview window that updates as you change\n"
            "settings (zoom and pan to see detail). Apply Preview then saves it at full resolution").pack(side="left", padx=4)
        self._update_apply_text()

        bp = section(p, "Create Data Boundary Polygon")
        bp.pack(fill="x", padx=8, pady=4)
        rb1 = tk.Frame(bp, bg=BG); rb1.pack(fill="x", padx=4, pady=2)
        tip(btn(rb1, "Create Polygon", self._create_boundary, w=14),
            "Create grid outline polygon(s) layer\nResult is loaded into the active ArcGIS Pro map").pack(side="left")

        cw = section(p, "Clip to Window")
        cw.pack(fill="x", padx=8, pady=4)
        rc1 = tk.Frame(cw, bg=BG); rc1.pack(fill="x", padx=4, pady=2)
        tip(btn(rc1, "Clip to Window", self._clip_to_window, w=16),
            "Clip the input grid to the area the ArcGIS Pro map is showing.\n"
            "Saved beside the original as <name>_clip.tif (or _clip_01, _clip_02 ...\n"
            "if that name is already used) and added to the layer list").pack(side="left")

        ng = section(p, "Normalise Grids (Batch)")
        ng.pack(fill="x", padx=8, pady=4)
        rn1 = tk.Frame(ng, bg=BG); rn1.pack(fill="x", padx=4, pady=2)
        lbl(rn1, "Input dir:").pack(side="left")
        self.v_norm_in = tk.StringVar()
        entry(rn1, self.v_norm_in, 34).pack(side="left", padx=2)
        btn(rn1, "Browse…", lambda: _browse_dir(self.v_norm_in), w=8).pack(side="left")
        rn3 = tk.Frame(ng, bg=BG); rn3.pack(fill="x", padx=4, pady=2)
        self.v_norm_order = tk.StringVar(value="1st")
        lbl(rn3, "Polynomial order:").pack(side="left")
        combo(rn3, ["1st","2nd"], self.v_norm_order, w=5).pack(side="left", padx=4)
        btn(rn3, "Normalise Grids", self._normalise_grids, w=16).pack(side="left", padx=10)

        rg = section(p, "Convert RGB Grid to Grayscale")
        rg.pack(fill="x", padx=8, pady=4)
        rrg1 = tk.Frame(rg, bg=BG); rrg1.pack(fill="x", padx=4, pady=2)
        lbl(rrg1, "Input RGB raster:").pack(side="left")
        self.v_rgb_in = tk.StringVar()
        entry(rrg1, self.v_rgb_in, 32).pack(side="left", padx=2)
        btn(rrg1, "Browse…", self._browse_rgb_in, w=8).pack(side="left")
        lbl(rg, "Colour list (CSS names or R,G,B – comma separated on one line):").pack(anchor="w", padx=4)
        self.v_rgb_colours = tk.Text(rg, height=5, bg=ENTRY_BG, fg=FG,
                                     font=FONT_MONO, relief="sunken", bd=1)
        self.v_rgb_colours.pack(fill="x", padx=4, pady=2)
        self.v_rgb_colours.insert("1.0", "")
        rrg2 = tk.Frame(rg, bg=BG); rrg2.pack(fill="x", padx=4, pady=2)
        tip(btn(rrg2, "Pick RGB from image", self._pick_rgb_from_image, w=20),
            "Open a preview of the RGB raster and click anywhere to sample\n"
            "the R,G,B values at that pixel — appended to the colour list").pack(side="left")
        rrg3 = tk.Frame(rg, bg=BG); rrg3.pack(fill="x", padx=4, pady=2)
        lbl(rrg3, "Min val:").pack(side="left")
        self.v_rgb_min = tk.StringVar(value="0.0")
        entry(rrg3, self.v_rgb_min, 7).pack(side="left", padx=2)
        lbl(rrg3, "Max val:").pack(side="left", padx=(8,0))
        self.v_rgb_max = tk.StringVar(value="1000.0")
        entry(rrg3, self.v_rgb_max, 7).pack(side="left", padx=2)
        btn(rrg3, "Convert Grid", self._convert_rgb, w=14).pack(side="left", padx=12)

    def _selected_metadata(self):
        """(path, provenance) of the input grid, or None after saying why not."""
        path = self._resolve_input()
        if not path or not os.path.exists(path):
            messagebox.showwarning("No Input", "Select an input raster.")
            return None
        root = read_sgt_metadata(path)
        if root is None:
            messagebox.showinfo("SGTool metadata",
                                f"No SGTool metadata found for {os.path.basename(path)}")
            return None
        return path, root

    def _show_metadata(self):
        found = self._selected_metadata()
        if not found:
            return
        path, root = found
        win = tk.Toplevel(self)
        win.title(f"SGTool metadata: {os.path.basename(path)}")
        win.configure(bg=BG)
        win.geometry("760x600")
        nb = ttk.Notebook(win)
        nb.pack(fill="both", expand=True, padx=4, pady=4)
        for title, text in (("History", format_sgt_metadata_text(root)),
                            ("XML", sgt_metadata_xml_text(root))):
            frame = tk.Frame(nb, bg=BG)
            box = tk.Text(frame, wrap="none", bg=ENTRY_BG, fg=FG, font=FONT_MONO,
                          relief="sunken", bd=1)
            sb_y = ttk.Scrollbar(frame, orient="vertical", command=box.yview)
            sb_x = ttk.Scrollbar(frame, orient="horizontal", command=box.xview)
            box.configure(yscrollcommand=sb_y.set, xscrollcommand=sb_x.set)
            sb_y.pack(side="right", fill="y")
            sb_x.pack(side="bottom", fill="x")
            box.pack(side="left", fill="both", expand=True)
            box.insert("1.0", text)
            box.configure(state="disabled")
            nb.add(frame, text=title)

    def _save_metadata_xml(self):
        found = self._selected_metadata()
        if not found:
            return
        path, root = found
        try:
            out = save_sgt_metadata_xml(root, path)
        except OSError as e:
            messagebox.showerror("SGTool metadata", f"Could not save metadata: {e}")
            return
        self._status(f"Saved metadata → {os.path.basename(out)}")

    def _open_replay(self):
        found = self._selected_metadata()
        if not found:
            return
        path, root = found
        _origin, steps = history_steps(root)
        plan = plan_steps(steps)
        if not any(p["replayable"] for p in plan):
            messagebox.showinfo("SGTool metadata",
                                f"The history of {os.path.basename(path)} has no steps that can be repeated.")
            return
        ReplayWindow(self, os.path.basename(path), plan)

    # ======================================================================
    # TAB 5 – HELP
    # ======================================================================
    def _tab_help(self, p):
        """The same text as the QGIS plugin's Help tab; the links open in the web browser."""
        txt = tk.Text(p, bg=BG, fg=FG, font=FONT, wrap="word",
                      relief="flat", padx=12, pady=8)
        vsb = ttk.Scrollbar(p, orient="vertical", command=txt.yview)
        txt.configure(yscrollcommand=vsb.set)
        txt.tag_configure("bold", font=FONT_BOLD)
        txt.tag_configure("link", foreground="#0563c1", underline=True)
        n_links = 0
        for paragraph in HELP_TEXT:
            for kind, text, url in _help_pieces(paragraph):
                if kind == "link":
                    n_links += 1
                    tag = f"link{n_links}"
                    txt.insert("end", text, ("link", tag))
                    txt.tag_bind(tag, "<Button-1>", lambda _e, u=url: webbrowser.open(u))
                    txt.tag_bind(tag, "<Enter>", lambda _e: txt.configure(cursor="hand2"))
                    txt.tag_bind(tag, "<Leave>", lambda _e: txt.configure(cursor=""))
                else:
                    txt.insert("end", text, ("bold",) if kind == "bold" else ())
            txt.insert("end", "\n\n")
        txt.configure(state="disabled")
        vsb.pack(side="right", fill="y", pady=4)
        txt.pack(side="left", fill="both", expand=True, padx=4, pady=4)
        self._help_text = txt

    # ======================================================================
    # PROCESSING HELPERS
    # ======================================================================
    @staticmethod
    def _ensure_odd(n):
        """Return n if odd, n+1 if even (minimum 1)."""
        n = max(1, int(n))
        return n if n % 2 == 1 else n + 1

    @staticmethod
    def _setup_odd_trace(var):
        """Attach a write-trace to var that immediately corrects even integers to the next odd value."""
        _lock = [False]
        def _enforce(*_):
            if _lock[0]:
                return
            try:
                n = int(var.get())
                if n > 0 and n % 2 == 0:
                    _lock[0] = True
                    var.set(str(n + 1))
                    _lock[0] = False
            except ValueError:
                pass
        var.trace_add("write", _enforce)

    def _run_in_thread(self, func, *args):
        token = threading.Event()          # set to cancel this job
        self._active_token = token
        self._cancel_btn.config(state="normal")

        def _wrap():
            self._tls.token = token
            try:
                func(*args)
            except OperationCancelled:
                pass                       # cancelled: the status bar already says so
            except Exception as e:
                msg = str(e)
                tb  = traceback.format_exc()
                self.after(0, lambda m=msg, t=tb: messagebox.showerror("Error", f"{m}\n\n{t}"))
                self.after(0, lambda m=msg: self._status(f"Error: {m}"))
            finally:
                self.after(0, self._job_finished, token)
        threading.Thread(target=_wrap, daemon=True).start()

    def _job_finished(self, token):
        if self._active_token is token:
            self._active_token = None
            self._cancel_btn.config(state="disabled")

    def _cancel_job(self):
        """End the running job now. The calculation it is part-way through finishes
        in the background but is thrown away: it stops at its next checkpoint without
        saving anything."""
        if self._active_token is not None:
            self._active_token.set()
            self._active_token = None
            self._cancel_btn.config(state="disabled")
            self._status("Cancelled.")

    def _check_cancel(self):
        """Called by a job between steps and before saving: stop if it was cancelled."""
        token = getattr(self._tls, "token", None)
        if token is not None and token.is_set():
            raise OperationCancelled()

    def _resolve_input(self):
        """Return actual file path for whatever is in the input field."""
        val = self.v_input.get()
        return self._layer_map.get(val, val)

    # ======================================================================
    # FFT APPLY  (buffer passed to every method)
    # ======================================================================
    def _apply_fft(self):
        in_r = self._compute_output_path()
        if not in_r:
            messagebox.showwarning("No Input", "Select an input raster.")
            return
        self._end_preview_for_apply()
        self._status("Starting processing…")
        self._run_in_thread(self._do_apply_fft, in_r)

    def _do_apply_fft(self, in_r):
        arr, nodata, ll_x, ll_y, cx, cy, sr = raster_to_numpy(in_r)
        max_buf = int(self.v_buf.get() or 10)
        rows, cols = arr.shape
        buf  = min(max_buf, rows, cols)
        proc = make_processor(cx, cy, buf)
        base, ext = os.path.splitext(in_r)
        if not ext:
            ext = ".tif"
        count = [0]

        def _save(result, suffix, step=None):
            self._check_cancel()
            suffix = step_suffix(step) or suffix   # QGIS naming: grid_RTP.tif, grid_d1z.tif ...
            out_r = f"{base}{suffix}{ext}"
            numpy_to_raster(result, out_r, ll_x, ll_y, cx, cy, sr, nodata)
            # how it was made, inside the GeoTIFF: (operation, parameters)
            write_provenance(out_r, in_r, *(step or (suffix, {})))
            self._record_pending(out_r)
            count[0] += 1
            self.after(0, self._on_output_written, out_r)

        self._fft_steps(in_r, arr, cx, cy, buf, proc, count, _save)

        self.after(0, lambda: self._status(f"Done — {count[0]} output(s) written."))
        self.after(0, self._uncheck_fft)

    def _fft_steps(self, in_r, arr, cx, cy, buf, proc, count, _save):
        """The ticked FFT filters; each result goes to _save(result, suffix, step)
        (which writes it, or for the preview just keeps it)."""
        if self.v_rtp.get():
            self._check_cancel()
            rt  = self.v_rtp_type.get()
            inc = float(self.v_inc.get())
            dec = float(self.v_dec.get())
            if rt == "Diff. RTP":
                self.after(0, lambda: self._status("Differential (variable) RTP…"))
                year = float(self.v_igrf_year.get())
                inc_c, dec_c, inc_ctr, dec_ctr = self._get_igrf_corners(in_r, year)
                if inc_c is None:
                    self.after(0, lambda: messagebox.showwarning(
                        "Diff. RTP",
                        "Could not compute IGRF at grid corners.\n"
                        "Check that the raster has a valid CRS and the IGRF Year is correct."))
                else:
                    from calcs.Cooper_Cowan_rtpvariable import rtpvariable as cooper_rtpvariable
                    nr, nc = arr.shape
                    data = np.where(np.isnan(arr), 0.0, arr)
                    _, varrtp, incv_deg, decv_deg = cooper_rtpvariable(
                        np.flipud(data), nr, nc, inc_c, dec_c,
                        inc_center=inc_ctr, dec_center=dec_ctr,
                    )
                    when = self._year_to_date(year)
                    _save(np.flipud(varrtp),    "_DRTP",
                          ("Differential reduction to pole", {"date": when}))
                    _save(np.flipud(incv_deg),  "_DRTP_inc",
                          ("Differential RTP: interpolated inclination field", {"date": when}))
                    _save(np.flipud(decv_deg),  "_DRTP_dec",
                          ("Differential RTP: interpolated declination field", {"date": when}))
            elif rt == "Pole":
                self.after(0, lambda: self._status("RTP…"))
                _save(proc.reduction_to_pole(arr, inc, dec, buffer_size=buf), "rtp",
                      step_record("reduction_to_pole", INCLINATION=inc, DECLINATION=dec))
            else:
                self.after(0, lambda: self._status("RTE…"))
                _save(proc.reduction_to_equator(arr, inc, dec, buffer_size=buf), "rte",
                      step_record("reduction_to_equator", INCLINATION=inc, DECLINATION=dec))

        if self.v_vint.get():
            self._check_cancel()
            self.after(0, lambda: self._status("Vertical integration…"))
            _save(proc.vertical_integration(arr, buffer_size=buf), "vi",
                  step_record("vertical_integration"))

        if self.v_cont.get():
            self._check_cancel()
            h = float(self.v_cont_h.get())
            self.after(0, lambda: self._status(f"Continuation {h}…"))
            if self.v_cont_dir.get() == "up":
                _save(proc.upward_continuation(arr, h, buffer_size=buf), "cont_up",
                      step_record("continuation", DIRECTION="up", HEIGHT=h))
            else:
                _save(proc.downward_continuation(arr, h, buffer_size=buf), "cont_dn",
                      step_record("continuation", DIRECTION="down", HEIGHT=h))

        if self.v_dirclean.get():
            self._check_cancel()
            self.after(0, lambda: self._status("Directional Butterworth…"))
            scale  = float(self.v_dc_scale.get() or "1.0")
            wedge  = float(self.v_dc_wedge.get() or "45")
            min_sp = float(self.v_dc_wl.get())
            try:
                max_sp = float(self.v_dc_wl_max.get())
            except ValueError:
                max_sp = min_sp  # blank/invalid: constant line spacing
            max_sp = max(max_sp, min_sp)
            # zero-centred line-noise estimate (2x min to 10x max spacing,
            # within a 45 degree wedge about the azimuth), scaled and
            # subtracted from the original grid
            noise = proc.directional_butterworth_band_pass(
                arr,
                2 * min_sp, 10 * max_sp,
                direction_angle=float(self.v_dc_az.get()),
                direction_width=wedge,
                order=4,
                buffer_size=buf,
                preserve_dc=False)
            _save(arr - noise * scale, "dir",
                  step_record("directional_line_noise",
                              AZIMUTH=float(self.v_dc_az.get()),
                              LINE_SPACING_MIN=min_sp, LINE_SPACING_MAX=max_sp,
                              SCALE=scale, WEDGE=wedge, OUTPUT_TYPE="Corrected grid"))

        if self.v_regrem.get():
            self._check_cancel()
            self.after(0, lambda: self._status("Remove regional…"))
            mask = ~np.isnan(arr)
            if self.v_regorder.get() == "1st":
                _save(proc.remove_gradient(arr, mask), "reg1",
                      step_record("remove_regional", ORDER="1st order"))
            else:
                _save(proc.remove_2o_gradient(arr, mask), "reg2",
                      step_record("remove_regional", ORDER="2nd order"))

        if self.v_bp.get():
            self._check_cancel()
            self.after(0, lambda: self._status("Band-pass filter…"))
            w = float(self.v_bp_w.get()) or None
            _save(proc.band_pass_filter(
                arr, float(self.v_bp_lo.get()), float(self.v_bp_hi.get()),
                w, w, buffer_size=buf), "bp",
                step_record("band_pass", LOW_CUT=float(self.v_bp_lo.get()),
                            HIGH_CUT=float(self.v_bp_hi.get()), WIDTH=w or 0.0))

        if self.v_hlp.get():
            self._check_cancel()
            self.after(0, lambda: self._status("High/Low-pass filter…"))
            cut = float(self.v_hlp_cut.get())
            wid = float(self.v_hlp_w.get())
            if self.v_hlp_type.get() == "High":
                _save(proc.high_pass_filter(arr, cut, wid, buffer_size=buf), "hp",
                      step_record("high_low_pass", TYPE="High", CUTOFF=cut, WIDTH=wid))
            else:
                _save(proc.low_pass_filter(arr, cut, wid, buffer_size=buf), "lp",
                      step_record("high_low_pass", TYPE="Low", CUTOFF=cut, WIDTH=wid))

        if self.v_agc.get():
            self._check_cancel()
            self.after(0, lambda: self._status("AGC…"))
            agc_n = self._ensure_odd(self.v_agc_win.get())
            _save(proc.automatic_gain_control(arr, agc_n), "agc",
                  step_record("automatic_gain_control", WINDOW=agc_n))

        if self.v_deriv.get():
            self._check_cancel()
            self.after(0, lambda: self._status("Derivative…"))
            d = self.v_deriv_dir.get()
            _save(proc.compute_derivative(
                arr, direction=d,
                order=float(self.v_deriv_pow.get()), buffer_size=buf), f"d{d}",
                step_record("derivative", DIRECTION=d, POWER=float(self.v_deriv_pow.get())))

        if self.v_thg.get():
            self._check_cancel()
            self.after(0, lambda: self._status("Total Horizontal Gradient…"))
            _save(proc.total_hz_grad(arr, buffer_size=buf), "thg",
                  step_record("total_horizontal_gradient"))

        if self.v_ta.get():
            self._check_cancel()
            self.after(0, lambda: self._status("Tilt angle…"))
            _save(np.degrees(proc.tilt_angle(arr, buffer_size=buf)), "ta",
                  step_record("tilt_angle"))

        if self.v_as_.get():
            self._check_cancel()
            self.after(0, lambda: self._status("Analytic signal…"))
            _save(proc.analytic_signal(arr, buffer_size=buf), "as",
                  step_record("analytic_signal"))

    def _uncheck_fft(self):
        for v in (self.v_rtp, self.v_vint, self.v_cont, self.v_dirclean,
                  self.v_regrem, self.v_bp, self.v_hlp, self.v_agc,
                  self.v_deriv, self.v_thg, self.v_ta, self.v_as_):
            v.set(False)

    # ======================================================================
    # CONV + STATS APPLY
    # ======================================================================
    def _apply_conv_stats(self):
        in_r = self._compute_output_path()
        if not in_r:
            messagebox.showwarning("No Input", "Select an input raster.")
            return
        self._end_preview_for_apply()
        self._run_in_thread(self._do_apply_conv_stats, in_r)

    def _do_apply_conv_stats(self, in_r):
        arr, nodata, ll_x, ll_y, cx, cy, sr = raster_to_numpy(in_r)
        from calcs.ConvolutionFilter import ConvolutionFilter
        from calcs.SpatialStats      import SpatialStats
        base, ext = os.path.splitext(in_r)
        if not ext:
            ext = ".tif"
        count = [0]

        def _save(result, suffix, step=None):
            self._check_cancel()
            suffix = step_suffix(step) or suffix   # QGIS naming: grid_RTP.tif, grid_d1z.tif ...
            out_r = f"{base}{suffix}{ext}"
            numpy_to_raster(result, out_r, ll_x, ll_y, cx, cy, sr, nodata)
            # how it was made, inside the GeoTIFF: (operation, parameters)
            write_provenance(out_r, in_r, *(step or (suffix, {})))
            self._record_pending(out_r)
            count[0] += 1
            self.after(0, self._on_output_written, out_r)

        self._conv_steps(in_r, arr, nodata, ll_x, ll_y, cx, cy, sr, base, ext, count, _save)

        self.after(0, lambda: self._status(f"Done — {count[0]} output(s) written."))
        self.after(0, self._uncheck_conv)

    def _conv_steps(self, in_r, arr, nodata, ll_x, ll_y, cx, cy, sr, base, ext, count, _save):
        """The ticked Conv + Stats operations; each result goes to _save(result, suffix, step)."""
        from calcs.ConvolutionFilter import ConvolutionFilter
        from calcs.SpatialStats      import SpatialStats
        if self.v_mean.get():
            self._check_cancel()
            self.after(0, lambda: self._status("Mean filter…"))
            n = self._ensure_odd(self.v_mean_n.get())
            _save(ConvolutionFilter(arr).mean_filter(n), f"mean{n}",
                  step_record("mean_filter", SIZE=n))

        if self.v_med.get():
            self._check_cancel()
            self.after(0, lambda: self._status("Median filter…"))
            n = self._ensure_odd(self.v_med_n.get())
            _save(ConvolutionFilter(arr).median_filter(n), f"med{n}",
                  step_record("median_filter", SIZE=n))

        if self.v_gauss.get():
            self._check_cancel()
            self.after(0, lambda: self._status("Gaussian filter…"))
            sigma = float(self.v_gauss_s.get())
            _save(ConvolutionFilter(arr).gaussian_filter(sigma), "gauss",
                  step_record("gaussian_filter", SIGMA=sigma))

        if self.v_dconv.get():
            self._check_cancel()
            self.after(0, lambda: self._status("Directional filter…"))
            d = self.v_dconv_dir.get().lower()
            _save(ConvolutionFilter(arr).directional_filter(
                self.v_dconv_dir.get(), self._ensure_odd(self.v_dconv_n.get())), f"dir{d}",
                step_record("directional_filter", DIRECTION=self.v_dconv_dir.get()))

        if self.v_sun.get():
            self._check_cancel()
            self.after(0, lambda: self._status("Sun shading…"))
            alt = float(self.v_sun_alt.get())
            az  = float(self.v_sun_az.get())
            cf  = ConvolutionFilter(arr)
            result = (cf.sun_shading_filter_grass(arr, cy, cx, altitude=alt, azimuth=az)
                      if self.v_sun_relief.get()
                      else cf.sun_shading_filter(arr, alt, az))
            _save(result, "shade",
                  step_record("sun_shading", AZIMUTH=az, ZENITH=alt,
                              RELIEF=bool(self.v_sun_relief.get())))

        win = self._ensure_odd(self.v_ss_win.get() or 5)
        for v_chk, stat, sfx in [
            (self.v_ss_min,  "min",      "min"),
            (self.v_ss_max,  "max",      "max"),
            (self.v_ss_var,  "variance", "var"),
            (self.v_ss_std,  "std",      "std"),
            (self.v_ss_skew, "skewness", "skew"),
            (self.v_ss_kurt, "kurtosis", "kurt"),
        ]:
            if v_chk.get():
                self._check_cancel()
                self.after(0, lambda s=stat: self._status(f"Spatial {s}…"))
                _save(SpatialStats(arr).calculate_windowed_stats(win, stat), sfx,
                      step_record("windowed_statistic", STATISTIC=stat, WINDOW=win))

        if self.v_mrvbf.get():
            self._check_cancel()
            hdr = raster_header_info(in_r)
            if hdr and hdr["geographic"]:
                self.after(0, lambda: messagebox.showwarning(
                    "MRVBF", "MRVBF needs a projected raster with cell sizes in metres.\n"
                             "Reproject the grid to a projected CRS and try again."))
            else:
                self.after(0, lambda: self._status("MRVBF / MRRTF / slope…"))
                for path in write_mrvbf_outputs(
                        in_r, base, ext, arr, ll_x, ll_y, cx, cy, sr, nodata,
                        t_slope=float(self.v_mrvbf_ts.get()), t_pctl_v=float(self.v_mrvbf_pv.get()),
                        t_pctl_r=float(self.v_mrvbf_pr.get()), p_slope=float(self.v_mrvbf_ps.get()),
                        p_pctl=float(self.v_mrvbf_pp.get()), max_res=float(self.v_mrvbf_mr.get()),
                        callback=lambda fraction=None: self._check_cancel()):
                    self._check_cancel()
                    self._record_pending(path)
                    self.after(0, self._on_output_written, path)
                    count[0] += 1

        if self.v_ss_aniso.get() or self.v_ss_chain.get() or self.v_ss_stream.get():
            self._check_cancel()
            from calcs.SpatialStats import SpatialStats as _SS
            win_type  = 'gaussian' if self.v_ss_aniso_win_type.get() == 'Gaussian' else 'box'
            a_thresh  = float(self.v_ss_aniso_thresh.get() or "0.3")
            ang_tol   = float(self.v_ss_angle_tol.get()    or "22.5")
            search_r  = int(float(self.v_ss_search_r.get() or "2"))
            max_steps = int(float(self.v_ss_max_steps.get() or "200"))
            self.after(0, lambda: self._status("Computing local anisotropy…"))
            ss_inst = _SS(arr)
            aniso_map, orient_map = ss_inst.calculate_local_anisotropy(
                arr, window_size=win, window_type=win_type)
            if self.v_ss_aniso.get():
                self._check_cancel()
                aniso_step = ("Local anisotropy", {"window_size": win,
                                                   "window_type": win_type})
                _save(aniso_map,  "_SS_AnisoMag", aniso_step)
                _save(orient_map, "_SS_AnisoOrient", aniso_step)
            if self.v_ss_chain.get():
                self._check_cancel()
                self.after(0, lambda: self._status("Chain length…"))
                _save(ss_inst.calculate_anisotropy_chain_length(
                    aniso_map, orient_map,
                    aniso_threshold=a_thresh,
                    angle_tolerance=ang_tol,
                    search_radius=search_r,
                ), "_SS_ChainLen",
                ("Anisotropy chain length", {
                    "window_size": win, "anisotropy_threshold": a_thresh,
                    "angle_tolerance": ang_tol, "search_radius": search_r}))
            if self.v_ss_stream.get():
                self._check_cancel()
                self.after(0, lambda: self._status("Streamline length…"))
                _save(ss_inst.calculate_anisotropy_streamline_length(
                    aniso_map, orient_map,
                    aniso_threshold=a_thresh,
                    angle_tolerance=ang_tol,
                    max_steps=max_steps,
                ), "_SS_StreamLen",
                ("Anisotropy streamline length", {
                    "window_size": win, "anisotropy_threshold": a_thresh,
                    "angle_tolerance": ang_tol, "max_steps": max_steps}))

        if self.v_pca.get():
            self._check_cancel()
            self.after(0, lambda: self._status("PCA…"))
            from calcs.PCAICA import PCAICA
            n = int(self.v_pca_n.get() or 0)
            PCAICA(np.zeros((2, 2))).pca_with_nans(
                in_r, f"{base}_PCA.tif", n_components=n if n > 0 else None)
            write_provenance(f"{base}_PCA.tif", in_r, "Principal component analysis",
                             {"n_components": n})
            self._record_pending(f"{base}_PCA.tif")
            self.after(0, self._on_output_written, f"{base}_PCA.tif")
            count[0] += 1

        if self.v_ica.get():
            self._check_cancel()
            self.after(0, lambda: self._status("ICA…"))
            from calcs.PCAICA import PCAICA
            n = int(self.v_ica_n.get() or 0)
            PCAICA(np.zeros((2, 2))).ica_with_nans(
                in_r, f"{base}_ICA.tif", n_components=n if n > 0 else None)
            write_provenance(f"{base}_ICA.tif", in_r, "Independent component analysis",
                             {"n_components": n})
            self._record_pending(f"{base}_ICA.tif")
            self.after(0, self._on_output_written, f"{base}_ICA.tif")
            count[0] += 1

        if self.v_euler.get():
            self._check_cancel()
            self.after(0, lambda: self._status("Euler deconvolution…"))
            self._do_euler(arr, nodata, ll_x, ll_y, cx, cy, in_r)
            count[0] += 1

    def _uncheck_conv(self):
        for v in (self.v_mean, self.v_med, self.v_gauss, self.v_dconv,
                  self.v_sun, self.v_sun_relief,
                  self.v_ss_min, self.v_ss_max, self.v_ss_var, self.v_ss_std,
                  self.v_ss_skew, self.v_ss_kurt,
                  self.v_ss_aniso, self.v_ss_chain, self.v_ss_stream,
                  self.v_mrvbf, self.v_pca, self.v_ica, self.v_euler):
            v.set(False)

    def _do_euler(self, arr, nodata, ll_x, ll_y, cx, cy, in_r):
        from calcs.euler.euler_python_optimised import euler_deconv_optimized
        import csv

        rows, cols = arr.shape
        x_max = ll_x + cols * cx
        y_max = ll_y + rows * cy
        area  = (ll_y, y_max, ll_x, x_max)
        xs = np.linspace(ll_x, x_max, cols)
        ys = np.linspace(ll_y, y_max, rows)
        XI, YI = np.meshgrid(xs, ys)
        ZI     = np.zeros_like(arr)
        filled = np.where(np.isnan(arr), 0.0, arr)

        # Returns array shape (n_keep, 5): columns [X, Y, Depth, B, StdDfDz]
        res = euler_deconv_optimized(
            filled, XI, YI, ZI, filled.shape, area,
            float(self.v_eu_si.get()),
            int(self.v_eu_win.get()),
            float(self.v_eu_filt.get()))
        if res is None or len(res) == 0:
            self.after(0, lambda: self._status("Euler: no estimates produced."))
            return

        # QGIS naming: <grid>_estimates_SI_<n>.txt, n being the structural index
        out_csv = os.path.splitext(in_r)[0] + f"_estimates_SI_{int(float(self.v_eu_si.get()))}.txt"
        with open(out_csv, "w", newline="") as f:
            wr = csv.writer(f)
            wr.writerow(["X", "Y", "Depth", "B", "StdDfDz"])
            for row in res:
                wr.writerow([row[0], row[1], row[2], row[3], row[4]])
        write_provenance(out_csv, in_r, "Euler deconvolution", {
            "structural_index": self.v_eu_si.get(), "window_size": self.v_eu_win.get(),
            "fraction_kept": self.v_eu_filt.get()})
        n = len(res)
        self.after(0, lambda: self._status(f"Euler: {n} estimates → {os.path.basename(out_csv)}"))

    # ======================================================================
    # GRID + WAVELETS
    # ======================================================================
    def _import_data(self):
        in_f = self.v_import_file.get()
        if not in_f:
            messagebox.showwarning("No File", "Select an input data file.")
            return
        out_shp = os.path.splitext(in_f)[0] + ".shp"
        self._run_in_thread(self._do_import_data, in_f, out_shp)

    def _do_import_data(self, in_f, out_shp):
        self.after(0, lambda: self._status("Importing point data…"))
        import csv as csvmod
        rows_data = []
        with open(in_f, newline="", encoding="utf-8-sig") as f:
            sample = f.read(4096); f.seek(0)
            dialect = csvmod.Sniffer().sniff(sample, delimiters=",\t ")
            for row in csvmod.DictReader(f, dialect=dialect):
                rows_data.append(row)
        if not rows_data:
            self.after(0, lambda: self._status("No data found."))
            return

        headers = list(rows_data[0].keys())
        x_col   = self.v_x_col.get() or "long_x"
        y_col   = self.v_y_col.get() or "lat_y"
        for c in headers:
            if x_col not in headers and c.lower() in ("x","lon","longitude","easting"):
                x_col = c
            if y_col not in headers and c.lower() in ("y","lat","latitude","northing"):
                y_col = c

        if out_shp:
            from osgeo import ogr, osr
            driver = ogr.GetDriverByName("ESRI Shapefile")
            if os.path.exists(out_shp):
                driver.DeleteDataSource(out_shp)
            srs = osr.SpatialReference()
            srs.ImportFromEPSG(int(self.v_import_epsg.get() or 4326))
            ds    = driver.CreateDataSource(out_shp)
            layer = ds.CreateLayer("points", srs, ogr.wkbPoint)
            for h in headers:
                if h not in (x_col, y_col):
                    layer.CreateField(ogr.FieldDefn(h[:10], ogr.OFTString))
            ldef = layer.GetLayerDefn()
            for row in rows_data:
                try:
                    x, y = float(row[x_col]), float(row[y_col])
                except (ValueError, KeyError):
                    continue
                feat = ogr.Feature(ldef)
                pt   = ogr.Geometry(ogr.wkbPoint); pt.AddPoint(x, y)
                feat.SetGeometry(pt)
                for h in headers:
                    if h not in (x_col, y_col):
                        feat.SetField(h[:10], str(row.get(h, "")))
                layer.CreateFeature(feat)
            ds = None
        self.after(0, lambda: self._status(
            f"Imported {len(rows_data)} points" +
            (f" → {out_shp}" if out_shp else "")))
        if out_shp and os.path.exists(out_shp):
            self.after(0, self._set_bspline_points, out_shp)


    # ------------------------------------------------------------------
    # Multilevel B-spline gridding
    # ------------------------------------------------------------------
    def _browse_bspline_points(self):
        path = filedialog.askopenfilename(
            filetypes=[("Point files", "*.shp *.gpkg *.geojson"), ("All files", "*.*")])
        if path:
            self._set_bspline_points(path)

    def _set_bspline_points(self, path):
        """Use `path` as the points to grid and list its fields."""
        self.v_bs_file.set(path)
        fields = list_point_fields(path)
        self._bs_field_combo["values"] = fields
        if fields and self.v_bs_field.get() not in fields:
            self.v_bs_field.set(fields[-1])
        if self.v_bs_auto.get():
            self._bs_auto_apply()

    # ---- Auto cell size and Max levels
    def _bs_auto_toggled(self):
        if self.v_bs_auto.get():
            self._bs_auto_apply()
        else:
            self.v_bs_auto_info.set("")

    def _bs_user_edit(self, *_args):
        """Typing a cell size or a number of levels takes over from Auto."""
        if not self._bs_auto_busy and self.v_bs_auto.get():
            self.v_bs_auto.set(False)
            self.v_bs_auto_info.set("")

    def _bs_auto_apply(self):
        points = self.v_bs_file.get()
        if not points or not os.path.exists(points):
            self.v_bs_auto_info.set("Choose the points file first.")
            return
        self.v_bs_auto_info.set("Looking at the point spacing…")

        def work():
            try:
                result = suggest_bspline_settings(points)
            except Exception as e:  # unreadable file, too few points ...
                result = str(e)
            self.after(0, self._bs_auto_done, points, result)

        threading.Thread(target=work, daemon=True).start()

    def _bs_auto_done(self, points, result):
        if points != self.v_bs_file.get() or not self.v_bs_auto.get():
            return  # the file was changed or Auto switched off while this was running
        if isinstance(result, str):
            self.v_bs_auto_info.set("Cannot choose automatically: " + result)
            return
        self._bs_auto_busy = True
        try:
            self.v_bs_cell.set(f"{result['cell_size']:g}")
            self.v_bs_levels.set(str(result["max_levels"]))
        finally:
            self._bs_auto_busy = False
        self.v_bs_auto_info.set(result["description"])

    def _grid_bspline(self):
        points = self.v_bs_file.get()
        field = self.v_bs_field.get()
        if not points or not os.path.exists(points) or not field:
            messagebox.showwarning("Gridding", "Choose a points file and the value field to grid.")
            return
        try:
            cell = float(self.v_bs_cell.get())
            eps = float(self.v_bs_eps.get())
            levels = int(float(self.v_bs_levels.get()))
            ignore = float(self.v_bs_ignore_val.get()) if self.v_bs_ignore.get() else None
        except ValueError:
            messagebox.showwarning("Gridding", "Cell size, threshold error and levels must be numbers.")
            return
        if cell <= 0:
            messagebox.showwarning("Gridding", "The cell size must be greater than zero.")
            return
        self._status("Multilevel B-spline gridding…")
        self._run_in_thread(self._do_bspline, points, field, cell, eps, levels, ignore)

    def _do_bspline(self, points, field, cell, eps, levels, ignore):
        xs, ys, zs, wkt = read_points_ogr(points, field)
        out_r = os.path.splitext(points)[0] + f"_{field}_bspline.tif"  # as in QGIS
        self._check_cancel()
        n = grid_points_to_file(
            xs, ys, zs, out_r, points, field, cell, epsilon=eps, level_max=levels,
            ignore_below=ignore, spatial_ref=wkt or None,
            callback=lambda fraction=None: self._check_cancel())
        self._record_pending(out_r)
        self.after(0, self._on_output_written, out_r)
        self.after(0, lambda: self._status(
            f"Gridded {n} points → {os.path.basename(out_r)}"))

    def _calc_worms(self):
        in_r = self._resolve_input()
        if not in_r:
            messagebox.showwarning("Missing", "Select an input raster.")
            return
        stem  = os.path.splitext(os.path.basename(in_r))[0]
        out_d = os.path.join(os.path.dirname(in_r), f"{stem}_worms")
        self._run_in_thread(self._do_calc_worms, in_r, out_d)

    def _do_calc_worms(self, in_r, out_d):
        self.after(0, lambda: self._status("Calculating BSD Worms…"))
        arr, nodata, ll_x, ll_y, cx, cy, sr = raster_to_numpy(in_r)
        proc = make_processor(cx, cy, 10)
        try:
            crs = sr.exportToString() if hasattr(sr, "exportToString") else str(sr)
        except Exception:
            crs = ""
        os.makedirs(out_d, exist_ok=True)
        proc.bsdwormer(arr,
                       os.path.splitext(os.path.basename(in_r))[0], in_r,
                       int(self.v_worm_levels.get() or 10),
                       float(self.v_worm_base.get() or 1000.0),
                       float(self.v_worm_inc.get()  or 1000.0),
                       bool(self.v_worm_shp.get()), crs)
        self.after(0, lambda: self._status(f"BSD Worms done → {out_d}"))

    # ======================================================================
    # UTILS
    # ======================================================================
    def _apply_threshold(self):
        in_r = self._resolve_input()
        if not in_r or not self.v_thresh_en.get():
            messagebox.showwarning("Missing", "Enable threshold and select input raster.")
            return
        out_r = os.path.splitext(in_r)[0] + "_Clean.tif"
        was_previewing = self.v_preview.get()
        self._end_preview_for_apply()
        self._run_in_thread(self._do_threshold, in_r, out_r, was_previewing)

    def _do_threshold(self, in_r, out_r, untick=False):
        self.after(0, lambda: self._status("Threshold to NaN…"))
        arr, nodata, ll_x, ll_y, cx, cy, sr = raster_to_numpy(in_r)

        def _save(result, suffix, step=None):
            numpy_to_raster(result, out_r, ll_x, ll_y, cx, cy, sr, nodata)
            write_provenance(out_r, in_r, *step)
            self._record_pending(out_r)
            self.after(0, self._on_output_written, out_r)

        self._threshold_step(arr, _save)
        if untick:  # an Apply Preview: leave nothing ticked, as the other tabs do
            self.after(0, lambda: self.v_thresh_en.set(False))

    def _threshold_step(self, arr, _save):
        """Threshold to NaN with the settings on the Utils tab; the result goes to
        _save(result, suffix, step) (which writes it, or for the preview keeps it)."""
        from calcs.SG_Util import SG_Util
        cond = self.v_thresh_cond.get()
        above, below = float(self.v_thresh_above.get()), float(self.v_thresh_below.get())
        _save(SG_Util(arr).Threshold2Nan(arr, cond, above, below), "_Clean",
              step_record("threshold_to_nan", CONDITION=cond, ABOVE=above, BELOW=below))

    def _clip_to_window(self):
        in_r = self._resolve_input()
        if not in_r:
            messagebox.showwarning("Missing", "Select an input raster.")
            return
        try:
            with open(os.path.join(tempfile.gettempdir(), "sgtool_view.json")) as f:
                view = json.load(f)
            cells = map_view_to_cells(view, in_r)
        except Exception:
            cells = None
        if cells is None:
            messagebox.showwarning(
                "Clip to Window",
                "The ArcGIS Pro map is not showing this grid (or the map view is not "
                "available: launch SGTool from the toolbox and keep the map open).")
            return
        x0, y0 = int(np.floor(cells[0])), int(np.floor(cells[1]))
        x1, y1 = int(np.ceil(cells[2])), int(np.ceil(cells[3]))
        # never overwrite an earlier clip: _clip, then _clip_01, _clip_02 ...
        base, n = os.path.splitext(in_r)[0], 0
        while True:
            out_r = base + ("_clip" if n == 0 else f"_clip_{n:02d}") + ".tif"
            if not os.path.exists(out_r):
                break
            n += 1
        self._run_in_thread(self._do_clip_to_window, in_r, out_r, x0, y0, x1, y1)

    def _do_clip_to_window(self, in_r, out_r, x0, y0, x1, y1):
        from osgeo import gdal
        self.after(0, lambda: self._status("Clipping to window…"))
        ds = gdal.Open(in_r)
        out = gdal.Translate(out_r, ds, srcWin=[x0, y0, x1 - x0, y1 - y0], format="GTiff")
        if out is None:
            raise RuntimeError("GDAL could not write the clipped grid")
        out = ds = None
        write_provenance(out_r, in_r, "Clip to map window",
                         {"window": f"cols {x0}-{x1}, rows {y0}-{y1}"})
        self._record_pending(out_r)
        self.after(0, self._on_output_written, out_r)

    def _create_boundary(self):
        in_r = self._resolve_input()
        if not in_r:
            messagebox.showwarning("Missing", "Select an input raster.")
            return
        out_shp = os.path.splitext(in_r)[0] + "_boundary.shp"
        self._run_in_thread(self._do_create_boundary, in_r, out_shp)

    def _do_create_boundary(self, in_r, out_shp):
        self.after(0, lambda: self._status("Creating boundary polygon…"))
        from calcs.SG_Util import SG_Util
        result = SG_Util(np.zeros((2, 2))).create_data_boundary_lines(in_r, out_shp)
        shp_path = result or out_shp
        self.after(0, lambda: self._status(f"Boundary → {shp_path}"))
        self.after(0, lambda: self._load_shapefile_to_arcpro(shp_path))

    def _load_shapefile_to_arcpro(self, shp_path):
        """Add a shapefile to the active ArcGIS Pro map."""
        if not shp_path or not os.path.exists(shp_path):
            return
        try:
            import arcpy
            proj = arcpy.mp.ArcGISProject("CURRENT")
            m    = proj.activeMap or (proj.listMaps() or [None])[0]
            if m:
                m.addDataFromPath(shp_path)
                self._status(f"Boundary polygon added to map: {os.path.basename(shp_path)}")
            else:
                self._status(f"Boundary created (no active map found): {os.path.basename(shp_path)}")
        except Exception:
            self._status(f"Boundary created (open in ArcPro manually): {os.path.basename(shp_path)}")

    def _normalise_grids(self):
        in_d = self.v_norm_in.get()
        if not in_d:
            messagebox.showwarning("Missing", "Set input directory.")
            return
        out_d = in_d.rstrip("/\\") + "_norm"
        self._run_in_thread(self._do_normalise, in_d, out_d,
                            1 if "1st" in self.v_norm_order.get() else 2)

    def _do_normalise(self, in_d, out_d, order):
        self.after(0, lambda: self._status("Normalising grids…"))
        os.makedirs(out_d, exist_ok=True)
        make_processor(1.0, 1.0, 10).normalise_geotiffs(in_d, out_d, order)
        self.after(0, lambda: self._status(f"Normalisation done → {out_d}"))

    def _convert_rgb(self):
        in_r = self.v_rgb_in.get()
        if not in_r:
            messagebox.showwarning("Missing", "Set RGB input path.")
            return
        out_r   = os.path.splitext(in_r)[0] + "_grey.tif"
        col_str = self.v_rgb_colours.get("1.0", "end").strip()
        self._run_in_thread(self._do_convert_rgb, in_r, out_r, col_str,
                            float(self.v_rgb_min.get()), float(self.v_rgb_max.get()))

    def _do_convert_rgb(self, in_r, out_r, col_str, min_v, max_v):
        # Direct port of QGIS convert_RGB_to_grey — algorithm unchanged.
        self.after(0, lambda: self._status("RGB→Grayscale…"))
        from osgeo import gdal
        from scipy.spatial import cKDTree

        dataset = gdal.Open(in_r, gdal.GA_ReadOnly)
        if not dataset:
            self.after(0, lambda: self._status(f"Cannot open {in_r}"))
            return
        if dataset.RasterCount < 3:
            self.after(0, lambda: self._status("Data file must have at least 3 layers."))
            return

        red   = dataset.GetRasterBand(1).ReadAsArray().astype(float)
        green = dataset.GetRasterBand(2).ReadAsArray().astype(float)
        blue  = dataset.GetRasterBand(3).ReadAsArray().astype(float)
        transform  = dataset.GetGeoTransform()
        projection = dataset.GetProjection()
        dataset = None

        rgb_raster = np.dstack((red, green, blue))

        LUT = col_str.replace(" ", "")
        lut = self._parse_lut_string(LUT, num_entries=1024)
        if not lut:
            self.after(0, lambda: self._status("Couldn't generate LUT — check colour list."))
            return

        scalar_values, lut_colors = zip(*lut)
        lut_colors = np.array(lut_colors) / 255.0  # normalise LUT to [0, 1]

        white_mask = (rgb_raster == [255, 255, 255]).all(axis=2)
        black_mask = (rgb_raster == [0, 0, 0]).all(axis=2)

        normalized_rgb = rgb_raster / 255.0
        reshaped_rgb   = normalized_rgb.reshape(-1, 3)

        lut_tree = cKDTree(lut_colors)
        distances, indices = lut_tree.query(reshaped_rgb)

        scalar_grid = np.array(scalar_values)[indices].reshape(rgb_raster.shape[:2])

        scalar_grid[white_mask] = np.nan
        scalar_grid[black_mask] = np.nan

        scalar_grid = (scalar_grid * (max_v - min_v)) + min_v

        rows, cols = red.shape
        driver = gdal.GetDriverByName("GTiff")
        out_ds = driver.Create(out_r, cols, rows, 1, gdal.GDT_Float32)
        out_ds.SetGeoTransform(transform)
        out_ds.SetProjection(projection)
        out_ds.GetRasterBand(1).WriteArray(scalar_grid)
        out_ds.GetRasterBand(1).SetNoDataValue(np.nan)
        out_ds.FlushCache()
        out_ds = None
        self._record_pending(out_r)
        self.after(0, lambda: self._status(f"RGB→Grayscale done → {out_r}"))

    def _browse_rgb_in(self):
        path = filedialog.askopenfilename(
            filetypes=[("Raster files", "*.tif *.tiff *.img"), ("All files", "*.*")])
        if not path:
            return
        self.v_rgb_in.set(path)
        self._record_pending(path)
        name = os.path.basename(path)
        self._layer_map[name] = path
        vals = list(self._layer_combo["values"])
        if name not in vals:
            vals.insert(0, name)
            self._layer_combo["values"] = vals
        self._status(f"RGB raster loaded: {name}")

    def _pick_rgb_from_image(self):
        in_r = self.v_rgb_in.get()
        if not in_r or not os.path.exists(in_r):
            messagebox.showwarning("No input", "Set the RGB input raster path first.")
            return

        try:
            from osgeo import gdal
            ds = gdal.Open(in_r, gdal.GA_ReadOnly)
            if ds is None or ds.RasterCount < 3:
                messagebox.showwarning("Error", "Cannot open raster or fewer than 3 bands.")
                return
            r_arr = ds.GetRasterBand(1).ReadAsArray()
            g_arr = ds.GetRasterBand(2).ReadAsArray()
            b_arr = ds.GetRasterBand(3).ReadAsArray()
            orig_h, orig_w = r_arr.shape
            ds = None
        except Exception as e:
            messagebox.showerror("Error opening raster", str(e))
            return

        try:
            from PIL import Image, ImageTk
        except ImportError:
            messagebox.showerror(
                "Missing library",
                "Pillow is required for image preview.\n"
                "Install it with:  pip install Pillow")
            return

        r_u8 = np.clip(r_arr, 0, 255).astype(np.uint8)
        g_u8 = np.clip(g_arr, 0, 255).astype(np.uint8)
        b_u8 = np.clip(b_arr, 0, 255).astype(np.uint8)
        rgb_img = Image.fromarray(np.dstack([r_u8, g_u8, b_u8]))

        max_size = 600
        scale = min(max_size / orig_w, max_size / orig_h, 1.0)
        thumb_w = max(1, int(orig_w * scale))
        thumb_h = max(1, int(orig_h * scale))
        thumb = rgb_img.resize((thumb_w, thumb_h), Image.LANCZOS)

        win = tk.Toplevel(self)
        win.title("Pick RGB – click image to sample colour")
        win.configure(bg=BG)
        win.resizable(False, False)

        tk.Label(win,
                 text="Click anywhere on the image to append that pixel's R,G,B to the colour list.\n"
                      "Click multiple times to pick several colours.",
                 bg=BG, fg=FG, font=FONT, justify="left").pack(padx=10, pady=(8, 2))

        canvas = tk.Canvas(win, width=thumb_w, height=thumb_h, cursor="crosshair",
                           bd=0, highlightthickness=0)
        canvas.pack(padx=10, pady=4)
        tk_img = ImageTk.PhotoImage(thumb)
        canvas.create_image(0, 0, anchor="nw", image=tk_img)
        canvas.image = tk_img  # keep reference

        info_var = tk.StringVar(value="")
        tk.Label(win, textvariable=info_var, bg=BG, fg=SECTION_FG, font=FONT).pack(pady=2)

        def _on_click(event):
            px = max(0, min(orig_w - 1, int(event.x / scale)))
            py = max(0, min(orig_h - 1, int(event.y / scale)))
            rv = int(r_arr[py, px])
            gv = int(g_arr[py, px])
            bv = int(b_arr[py, px])
            current = self.v_rgb_colours.get("1.0", "end").strip()
            new_triplet = f"{rv}, {gv}, {bv}"
            new_text = (current + ", " + new_triplet) if current else new_triplet
            self.v_rgb_colours.delete("1.0", "end")
            self.v_rgb_colours.insert("1.0", new_text)
            info_var.set(f"Pixel ({px}, {py}) → R={rv}  G={gv}  B={bv}  — added to colour list")

        canvas.bind("<Button-1>", _on_click)
        btn(win, "Done", win.destroy, w=10).pack(pady=(2, 10))

    def _parse_lut_string(self, LUT, num_entries=1024):
        """Port of QGIS parse_lut_string."""
        import re
        LUT = LUT.replace(" ", "")

        if "(" in LUT:
            rgb_pattern = r"\(([^)]+)\)"
            matches = re.findall(rgb_pattern, LUT)
            if matches:
                rgb_list = []
                for match in matches:
                    values = match.split(",")
                    try:
                        if len(values) >= 3:
                            rgb = tuple(float(v) for v in values[:3])
                            rgb_list.append(rgb)
                        else:
                            raise ValueError
                    except ValueError:
                        break
                else:
                    rgb_list.reverse()
                    return self._generate_rgb_lut_from_rgb(rgb_list, num_entries)

        elements = LUT.split(",")
        try:
            numbers = [float(elem) for elem in elements]
            if len(numbers) >= 3 and len(numbers) % 3 == 0:
                rgb_list = []
                for i in range(0, len(numbers), 3):
                    rgb_list.append(tuple(numbers[i:i+3]))
                rgb_list.reverse()
                return self._generate_rgb_lut_from_rgb(rgb_list, num_entries)
        except ValueError:
            pass

        css_color_list = elements
        css_color_list.reverse()
        return self._generate_rgb_lut(css_color_list, num_entries)

    def _generate_rgb_lut_from_rgb(self, rgb_list, num_entries=1024):
        """Port of QGIS generate_rgb_lut_from_rgb."""
        import matplotlib.colors as mcolors
        normalized_rgb_list = []
        for rgb in rgb_list:
            r, g, b = rgb[0], rgb[1], rgb[2]
            if any(val > 1 for val in [r, g, b]):
                normalized_rgb_list.append((r / 255.0, g / 255.0, b / 255.0))
            else:
                normalized_rgb_list.append((r, g, b))
        decimal_indices = np.linspace(0, 1, num_entries)
        try:
            cmap = mcolors.LinearSegmentedColormap.from_list("custom_cmap", normalized_rgb_list)
        except Exception:
            return False
        rgb_colors = [cmap(i)[:3] for i in decimal_indices]
        rgb_colors_255 = [(int(r * 255), int(g * 255), int(b * 255)) for r, g, b in rgb_colors]
        return [[round(di, 6), rgb] for di, rgb in zip(decimal_indices, rgb_colors_255)]

    def _generate_rgb_lut(self, css_color_list, num_entries=1024):
        """Port of QGIS generate_rgb_lut."""
        import matplotlib.colors as mcolors
        decimal_indices = np.linspace(0, 1, num_entries)
        try:
            cmap = mcolors.LinearSegmentedColormap.from_list("custom_cmap", css_color_list)
        except Exception:
            return False
        rgb_colors = [cmap(i)[:3] for i in decimal_indices]
        rgb_colors_255 = [(int(r * 255), int(g * 255), int(b * 255)) for r, g, b in rgb_colors]
        return [[round(di, 6), rgb] for di, rgb in zip(decimal_indices, rgb_colors_255)]

    # ======================================================================
    # IGRF  +  Power Spectrum
    # ======================================================================
    @staticmethod
    def _year_to_date(year):
        """A decimal year (2024.5) as an ISO date, for the provenance record."""
        import datetime
        y = int(year)
        days = (datetime.date(y + 1, 1, 1) - datetime.date(y, 1, 1)).days
        return (datetime.date(y, 1, 1)
                + datetime.timedelta(days=round((year - y) * days))).isoformat()

    def _calc_igrf(self):
        try:
            from arcgis_utils import (calc_igrf, raster_center_latlon,
                                      get_noddy_mag_params)

            in_r = self._resolve_input()

            # Noddy GeoTIFF: use the field stored in its metadata
            noddy = get_noddy_mag_params(in_r) if in_r else None
            if noddy is not None:
                inc, dec, F = noddy
                self.v_inc.set(f"{inc:.2f}")
                self.v_dec.set(f"{dec:.2f}")
                self.v_intensity.set(f"{F:.1f}")
                self._status(f"Noddy metadata: inc={inc:.2f}°  dec={dec:.2f}°  "
                             f"F={F:.1f} nT")
                return

            # Auto-populate Lat/Lon from the selected raster centre
            if in_r:
                centre = raster_center_latlon(in_r)
                if centre is not None:
                    lat, lon = centre
                    self.v_igrf_lat.set(f"{lat:.4f}")
                    self.v_igrf_lon.set(f"{lon:.4f}")

            res = calc_igrf(float(self.v_igrf_lat.get()),
                            float(self.v_igrf_lon.get()),
                            0.0, float(self.v_igrf_year.get()))
            if res:
                inc, dec, F = res
                self.v_inc.set(f"{inc:.2f}")
                self.v_dec.set(f"{dec:.2f}")
                self.v_intensity.set(f"{F:.1f}")
                self._status(f"IGRF: inc={inc:.2f}°  dec={dec:.2f}°  F={F:.1f} nT")
            else:
                self._status("IGRF failed – check parameters.")
        except Exception as e:
            messagebox.showerror("IGRF Error", str(e))

    def _get_igrf_corners(self, in_r, year):
        """Return (inc_corners, dec_corners, inc_center, dec_center) for Cooper/Cowan Diff. RTP.

        Corner order: [NW, NE, SW, SE] — the axis-xy convention expected by
        Cooper_Cowan_rtpvariable (row nr = geographic north, row 1 = south).
        Returns (None, None, None, None) on any failure.
        """
        from arcgis_utils import calc_igrf, get_noddy_mag_params
        noddy = get_noddy_mag_params(in_r)
        if noddy is not None:
            # Noddy GeoTIFF: uniform field, so centre equals corners
            inc, dec, _ = noddy
            return [inc] * 4, [dec] * 4, inc, dec
        try:
            from osgeo import gdal, osr
            ds = gdal.Open(str(in_r), gdal.GA_ReadOnly)
            if ds is None:
                return None, None, None, None
            gt   = ds.GetGeoTransform()
            cols = ds.RasterXSize
            rows = ds.RasterYSize
            wkt  = ds.GetProjection()
            ds   = None

            # Five points in raster CRS: NW, NE, SW, SE, centre
            pts_xy = [
                (gt[0],                                    gt[3]),                                    # NW
                (gt[0] + gt[1]*cols,                       gt[3] + gt[4]*cols),                       # NE
                (gt[0] + gt[2]*rows,                       gt[3] + gt[5]*rows),                       # SW
                (gt[0] + gt[1]*cols + gt[2]*rows,          gt[3] + gt[4]*cols + gt[5]*rows),           # SE
                (gt[0] + gt[1]*cols/2 + gt[2]*rows/2,      gt[3] + gt[4]*cols/2 + gt[5]*rows/2),      # centre
            ]

            src_srs = osr.SpatialReference()
            if wkt:
                src_srs.ImportFromWkt(wkt)
            else:
                src_srs.ImportFromEPSG(4326)

            if src_srs.IsGeographic():
                def to_latlon(x, y):
                    return float(y), float(x)
            else:
                geo_srs = osr.SpatialReference()
                geo_srs.ImportFromEPSG(4326)
                geo_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
                ct = osr.CoordinateTransformation(src_srs, geo_srs)
                def to_latlon(x, y):
                    lon, lat, _ = ct.TransformPoint(x, y)
                    return float(lat), float(lon)

            inc_corners, dec_corners = [], []
            for x, y in pts_xy[:4]:
                lat, lon = to_latlon(x, y)
                res = calc_igrf(lat, lon, 0.0, year)
                if res is None:
                    return None, None, None, None
                inc, dec, _ = res
                inc_corners.append(inc)
                dec_corners.append(dec)

            ctr_x, ctr_y = pts_xy[4]
            lat_c, lon_c = to_latlon(ctr_x, ctr_y)
            res_c = calc_igrf(lat_c, lon_c, 0.0, year)
            if res_c is None:
                return None, None, None, None
            inc_center, dec_center, _ = res_c

            return inc_corners, dec_corners, inc_center, dec_center
        except Exception:
            return None, None, None, None

    def _show_power_spectrum(self):
        in_r = self._resolve_input()
        if not in_r:
            messagebox.showwarning("No Input", "Select an input raster first.")
            return
        self._run_in_thread(self._do_power_spectrum, in_r)

    def _do_power_spectrum(self, in_r):
        self.after(0, lambda: self._status("Computing radial power spectrum…"))
        try:
            arr, _, ll_x, ll_y, cx, cy, sr = raster_to_numpy(in_r)
            filled, _ = make_processor(cx, cy, 10).fill_nan(arr)
            fft2   = np.fft.fft2(filled)
            power  = np.abs(np.fft.fftshift(fft2))**2
            nr, nc = power.shape
            yr, xr = np.indices(power.shape)
            r      = np.sqrt((xr - nc//2)**2 + (yr - nr//2)**2).astype(int)
            r_max  = min(nc//2, nr//2)
            bins   = np.arange(1, r_max)
            rad_p  = np.array([power[r==ri].mean() for ri in bins])
            cell   = min(cx, cy)
            dim    = min(nr, nc)
            wn     = bins / (dim * cell)   # wavenumber (cycles per map unit)

            def _plot():
                from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
                from matplotlib.figure import Figure

                win = tk.Toplevel(self)
                win.title(f"Radial Power Spectrum – {os.path.basename(in_r)}")
                win.configure(bg=BG)
                win.geometry("720x440")

                fig = Figure(figsize=(7, 4), tight_layout=True)
                ax  = fig.add_subplot(111)
                ax.plot(wn, np.log(rad_p + 1e-10), color="#0078d4")
                ax.set_xlabel("Wavenumber (cycles / map unit)")
                ax.set_ylabel("ln(Power)")
                ax.set_title(f"Radial Power Spectrum – {os.path.basename(in_r)}")
                ax.grid(True, alpha=0.4)

                canvas = FigureCanvasTkAgg(fig, master=win)
                canvas.draw()
                toolbar = NavigationToolbar2Tk(canvas, win)
                toolbar.update()
                canvas.get_tk_widget().pack(fill="both", expand=True)

            self.after(0, _plot)
            self.after(0, lambda: self._status("Power spectrum displayed."))
        except Exception as e:
            self.after(0, lambda msg=str(e): self._status(f"Power spectrum error: {msg}"))


# ============================================================================
class ReplayWindow(tk.Toplevel):
    """Apply the steps in a grid's history to another grid.

    Each step has an Apply box and a Save box (keep that step's result as a file; an
    unticked step is deleted once the next has been made from it, and the last result
    is always kept). Reduction to the pole / equator and differential RTP get their
    inclination and declination from IGRF for the new grid and the survey date."""

    TICK, UNTICK = "\u2611", "\u2610"

    def __init__(self, app, history_name, plan):
        super().__init__(app)
        self.app = app
        self.plan = plan
        self.apply = {i: p["replayable"] for i, p in enumerate(plan)}
        self.save = {i: p["replayable"] for i, p in enumerate(plan)}
        self.title("Apply the same steps to another grid")
        self.configure(bg=BG)
        self.geometry("860x470")

        usable = sum(1 for p in plan if p["replayable"])
        lbl(self, f"History of {history_name}: {len(plan)} steps, {usable} can be repeated. "
                  "Steps run oldest first, each on the result of the one before."
            ).pack(anchor="w", padx=8, pady=(8, 2))

        row = tk.Frame(self, bg=BG); row.pack(fill="x", padx=8, pady=3)
        lbl(row, "Apply to grid:").pack(side="left")
        self.v_target = tk.StringVar()
        ttk.Combobox(row, textvariable=self.v_target, values=list(app._layer_map),
                     width=58, font=FONT).pack(side="left", padx=4)
        btn(row, "Browse…", self._browse, w=9).pack(side="left")

        cols = ("apply", "save", "step", "settings", "note")
        self.tree = ttk.Treeview(self, columns=cols, show="headings", height=11,
                                 selectmode="none")
        for col, text, width, anchor in (
                ("apply", "Apply", 50, "center"), ("save", "Save", 50, "center"),
                ("step", "Step", 190, "w"), ("settings", "Settings recorded", 290, "w"),
                ("note", "Note", 240, "w")):
            self.tree.heading(col, text=text)
            self.tree.column(col, width=width, anchor=anchor, stretch=col in ("settings", "note"))
        for i, p in enumerate(plan):
            self.tree.insert("", "end", iid=str(i), values=(
                self.TICK if p["replayable"] else "", self.TICK if p["replayable"] else "",
                p["operation"] or "Unnamed step", p["settings"], p["note"]))
        self.tree.pack(fill="both", expand=True, padx=8, pady=4)
        self.tree.bind("<Button-1>", self._on_click)

        self.needs_igrf = any(p["igrf"] for p in plan)
        recorded = next((p["date"] for p in plan if p["date"] is not None), None)
        self.v_date = tk.StringVar(value=(recorded or __import__("datetime").date.today()).isoformat())
        if self.needs_igrf:
            rd = tk.Frame(self, bg=BG); rd.pack(fill="x", padx=8, pady=2)
            lbl(rd, "Inclination and declination for the reductions to the pole / equator are "
                    "recalculated from IGRF for the new grid.   Survey date (yyyy-mm-dd):"
                ).pack(side="left")
            entry(rd, self.v_date, 11).pack(side="left", padx=4)

        rb = tk.Frame(self, bg=BG); rb.pack(fill="x", padx=8, pady=(4, 8))
        tip(btn(rb, "Apply steps", self._apply, w=14),
            "Run the ticked steps on the chosen grid. Every ticked Save step is kept as\n"
            "<grid>_<suffix>.tif beside it, with its own history").pack(side="left")
        btn(rb, "Close", self.destroy, w=10).pack(side="left", padx=6)
        current = app._resolve_input()
        names = [n for n, pth in app._layer_map.items() if pth != current]
        if names:
            self.v_target.set(names[0])

    def _browse(self):
        path = filedialog.askopenfilename(
            filetypes=[("Raster files", "*.tif *.tiff *.img *.ers"), ("All files", "*.*")])
        if path:
            self.v_target.set(path)

    def _on_click(self, event):
        if self.tree.identify_region(event.x, event.y) != "cell":
            return
        row = self.tree.identify_row(event.y)
        col = self.tree.identify_column(event.x)
        if not row or col not in ("#1", "#2"):
            return
        i = int(row)
        if not self.plan[i]["replayable"]:
            return
        state = self.apply if col == "#1" else self.save
        state[i] = not state[i]
        self.tree.set(row, "apply" if col == "#1" else "save",
                      self.TICK if state[i] else self.UNTICK)

    def _apply(self):
        items = [dict(p, save=self.save[i]) for i, p in enumerate(self.plan)
                 if p["replayable"] and self.apply[i]]
        if not items:
            messagebox.showinfo("Apply steps", "No steps were chosen.", parent=self)
            return
        target = self.app._layer_map.get(self.v_target.get(), self.v_target.get())
        if not target or not os.path.exists(target):
            messagebox.showwarning("Apply steps", "Choose the grid to apply the steps to.",
                                   parent=self)
            return
        info = raster_header_info(target)
        if info is None:
            messagebox.showerror("Apply steps", f"Could not read {target}", parent=self)
            return
        try:
            when = __import__("datetime").date.fromisoformat(self.v_date.get().strip())
        except ValueError:
            messagebox.showwarning("Apply steps", "The survey date must be written yyyy-mm-dd.",
                                   parent=self)
            return
        if any(i["igrf"] for i in items) and uniform_field(target) is None:
            why = check_igrf_possible(info)
            if why:
                messagebox.showwarning("Apply steps", why, parent=self)
                return
        self.app._status(f"Applying {len(items)} steps to {os.path.basename(target)}…")
        self.app._run_in_thread(self._work, items, target, info, when)
        self.destroy()

    def _work(self, items, target, info, when):
        app = self.app

        class _Task:  # lets run_replay report progress and stop when Cancel is pressed
            @staticmethod
            def check(percent=None):
                app._check_cancel()
                if percent is not None:
                    app.after(0, app._status, f"Applying steps… {percent:.0f}%")

        result = run_replay(target, items, info, when=when, task=_Task())
        for path in result["outputs"]:
            app._record_pending(path)
            app.after(0, app._on_output_written, path)
        kept = len(result["outputs"])
        app.after(0, app._status,
                  f"Applied {result['n_steps']} steps → {os.path.basename(result['final'])} "
                  f"({kept} result{'s' if kept != 1 else ''} saved)")


def main():
    app = SGToolApp()
    app.mainloop()

if __name__ == "__main__":
    main()
