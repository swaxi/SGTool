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

import os
import sys
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

from arcgis_utils import raster_to_numpy, numpy_to_raster, make_processor, NODATA
import numpy as np

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
        canvas.bind_all("<MouseWheel>",
                        lambda e: canvas.yview_scroll(int(-1*(e.delta/120)), "units"))


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


# ============================================================================
# Main application
# ============================================================================
class SGToolApp(tk.Tk):

    def __init__(self):
        super().__init__()
        self.title("SGTool – Geophysical Processing (ArcGIS)")
        self.configure(bg=BG)
        self.minsize(700, 640)
        self.resizable(True, True)

        self._layer_map      = {}   # display_name → file_path
        self._info_file      = sys.argv[1] if len(sys.argv) > 1 else None
        self._pending_file   = os.path.join(tempfile.gettempdir(), "sgtool_pending.json")
        self._refresh_req    = os.path.join(tempfile.gettempdir(), "sgtool_refresh_req.json")
        self._apply_theme()
        self._build_top_bar()
        self._build_notebook()
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
            except Exception:
                pass

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
        except Exception:
            pass

        # Strategy 3: read the info file (kept fresh by the watcher every ~30 s)
        if not names and self._info_file and os.path.exists(self._info_file):
            try:
                for k, v in json.load(open(self._info_file)).get("layers", {}).items():
                    layer_map[k] = v
                    names.append(k)
            except Exception:
                pass

        # All tkinter updates must happen on the main thread
        self.after(0, lambda: self._apply_layer_list(layer_map, names))

    def _apply_layer_list(self, layer_map, names):
        self._layer_map = layer_map
        self._layer_combo["values"] = names
        if names and not self.v_input.get():
            self._layer_combo.set(names[0])
            self.v_input.set(self._layer_map[names[0]])
        self._status(f"Layer list refreshed – {len(names)} raster(s) found." if names
                     else "No raster layers found. Use Browse File or re-launch from ArcGIS Pro.")

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
        except Exception:
            pass

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
        except Exception:
            pass

    def _on_output_written(self, path):
        """Update status bar and add newly written raster to the layer dropdown."""
        name = os.path.basename(path)
        self._status(f"Written → {name}")
        self._layer_map[name] = path
        vals = list(self._layer_combo["values"])
        if name not in vals:
            vals.insert(0, name)
            self._layer_combo["values"] = vals

    # ------------------------------------------------------------------
    # Auto output-path derivation
    # ------------------------------------------------------------------
    def _compute_output_path(self):
        """Build output path from input path + suffix describing active operations."""
        in_r = self.v_input.get()
        # Resolve layer name → path if needed
        if in_r in self._layer_map:
            in_r = self._layer_map[in_r]
        if not in_r or not os.path.exists(in_r):
            return None

        parts = []
        # -- FFT --
        if getattr(self, "v_rtp",      tk.BooleanVar()).get():
            _rt = getattr(self, "v_rtp_type", tk.StringVar(value="Pole")).get()
            parts.append("rte" if _rt == "Equator" else ("DRTP" if _rt == "Diff. RTP" else "rtp"))
        if getattr(self, "v_vint",     tk.BooleanVar()).get(): parts.append("vi")
        if getattr(self, "v_cont",     tk.BooleanVar()).get():
            parts.append(f"cont_{getattr(self,'v_cont_dir',tk.StringVar(value='up')).get()}")
        if getattr(self, "v_dirclean", tk.BooleanVar()).get(): parts.append("dir")
        if getattr(self, "v_regrem",   tk.BooleanVar()).get():
            parts.append(f"reg{getattr(self,'v_regorder',tk.StringVar(value='1st')).get()[0]}")
        if getattr(self, "v_bp",       tk.BooleanVar()).get(): parts.append("bp")
        if getattr(self, "v_hlp",      tk.BooleanVar()).get():
            parts.append("hp" if getattr(self, "v_hlp_type",
                          tk.StringVar(value="High")).get() == "High" else "lp")
        if getattr(self, "v_agc",      tk.BooleanVar()).get(): parts.append("agc")
        if getattr(self, "v_deriv",    tk.BooleanVar()).get():
            parts.append(f"d{getattr(self,'v_deriv_dir',tk.StringVar(value='z')).get()}")
        if getattr(self, "v_thg",      tk.BooleanVar()).get(): parts.append("thg")
        if getattr(self, "v_ta",       tk.BooleanVar()).get(): parts.append("ta")
        if getattr(self, "v_as_",      tk.BooleanVar()).get(): parts.append("as")
        # -- Conv --
        if getattr(self, "v_mean",     tk.BooleanVar()).get():
            parts.append(f"mean{getattr(self,'v_mean_n',tk.StringVar(value='3')).get()}")
        if getattr(self, "v_med",      tk.BooleanVar()).get():
            parts.append(f"med{getattr(self,'v_med_n',tk.StringVar(value='3')).get()}")
        if getattr(self, "v_gauss",    tk.BooleanVar()).get(): parts.append("gauss")
        if getattr(self, "v_dconv",    tk.BooleanVar()).get():
            parts.append(f"dir{getattr(self,'v_dconv_dir',tk.StringVar(value='N')).get().lower()}")
        if getattr(self, "v_sun",      tk.BooleanVar()).get(): parts.append("shade")
        for attr, tag in [("v_ss_min","min"),("v_ss_max","max"),("v_ss_var","var"),
                          ("v_ss_std","std"),("v_ss_skew","skew"),("v_ss_kurt","kurt"),
                          ("v_ss_aniso","aniso"),("v_ss_chain","chain"),("v_ss_stream","stream")]:
            if getattr(self, attr, tk.BooleanVar()).get(): parts.append(tag)
        if getattr(self, "v_dtm",      tk.BooleanVar()).get(): parts.append("dtm")
        if getattr(self, "v_pca",      tk.BooleanVar()).get(): parts.append("pca")
        if getattr(self, "v_ica",      tk.BooleanVar()).get(): parts.append("ica")
        if getattr(self, "v_euler",    tk.BooleanVar()).get(): parts.append("euler")

        base, ext = os.path.splitext(in_r)
        if not ext:
            ext = ".tif"
        bn = os.path.basename(base)
        if parts:
            if len(parts) <= 3:
                preview = ",  ".join(f"{bn}_{p}{ext}" for p in parts)
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

    # ------------------------------------------------------------------
    # Status bar
    # ------------------------------------------------------------------
    def _build_status_bar(self):
        bar = tk.Frame(self, bg=FRAME_BG, height=22, relief="sunken", bd=1)
        bar.pack(fill="x", side="bottom")
        self.v_status = tk.StringVar(value="Ready.")
        tk.Label(bar, textvariable=self.v_status,
                 bg=FRAME_BG, fg=FG, font=FONT, anchor="w").pack(side="left", padx=6)

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
        self.v_dc_scale  = tk.StringVar(value="1.0")
        tip(check(r6, "Directional Butterworth", self.v_dirclean),
            "Filter (DirCos + Butterworth) to remove a specific direction and wavelength\n"
            "Useful for filtering flight line noise").pack(side="left")
        lbl(r6, "Az:").pack(side="left", padx=(8,0))
        tip(entry(r6, self.v_dc_az, 5, chk_var=self.v_dirclean),
            "Azimuth of high-frequency noise to be filtered\n"
            "(degrees clockwise from North)").pack(side="left", padx=2)
        lbl(r6, "Wavelength:").pack(side="left", padx=(6,0))
        tip(entry(r6, self.v_dc_wl, 9, chk_var=self.v_dirclean),
            "Wavelength of high-frequency noise to be filtered\n"
            "Set to 4× line spacing").pack(side="left", padx=2)
        lbl(r6, "Scale:").pack(side="left")
        tip(entry(r6, self.v_dc_scale, 5, chk_var=self.v_dirclean),
            "Multiplier applied to the filtered result").pack(side="left", padx=2)

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
        tip(btn(ba, "Apply FFT Processing", self._apply_fft, w=22),
            "Apply selected processing steps to the selected grid").pack(side="right", padx=4)

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

        rdm = tk.Frame(ss, bg=BG); rdm.pack(fill="x", padx=4, pady=2)
        self.v_dtm       = tk.BooleanVar()
        self.v_dtm_curv  = tk.StringVar(value="0.01")
        self.v_dtm_slope = tk.StringVar(value="45.0")
        self.v_dtm_sigma = tk.StringVar(value="1.0")
        tip(check(rdm, "DTM Curvature Classify", self.v_dtm),
            "Calculate DTM classification based on curvature and slope\n"
            "-1 = concave up  0 = flat  1 = convex up  2 = steep slope").pack(side="left")
        lbl(rdm, "Curv:").pack(side="left", padx=(8,0))
        tip(entry(rdm, self.v_dtm_curv, 5, chk_var=self.v_dtm),
            "Curvature threshold for DTM classification\nPositive = hill, negative = valley").pack(side="left", padx=2)
        lbl(rdm, "Slope°:").pack(side="left")
        tip(entry(rdm, self.v_dtm_slope, 5, chk_var=self.v_dtm),
            "Slope threshold for Steep Slope DTM classification").pack(side="left", padx=2)
        lbl(rdm, "σ:").pack(side="left")
        tip(entry(rdm, self.v_dtm_sigma, 5, chk_var=self.v_dtm),
            "Smoothing parameter for DTM classification\nHigher values smooth the data more").pack(side="left", padx=2)

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
        tip(btn(ba, "Apply Conv + Stats", self._apply_conv_stats, w=22),
            "Apply selected processing steps to the selected grid").pack(side="right", padx=4)

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
    # TAB 4 – UTILS
    # ======================================================================
    def _tab_utils(self, p):
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
        tip(btn(rt1, "Apply Threshold → NaN", self._apply_threshold, w=22),
            "Apply threshold and set matching values to NaN (no-data)").pack(side="left", padx=8)

        bp = section(p, "Create Data Boundary Polygon")
        bp.pack(fill="x", padx=8, pady=4)
        rb1 = tk.Frame(bp, bg=BG); rb1.pack(fill="x", padx=4, pady=2)
        tip(btn(rb1, "Create Polygon", self._create_boundary, w=14),
            "Create grid outline polygon(s) layer\nResult is loaded into the active ArcGIS Pro map").pack(side="left")

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
        btn(rrg1, "Browse…", lambda: _browse_file(self.v_rgb_in), w=8).pack(side="left")
        lbl(rg, "Colour list (CSS names or R,G,B – comma separated on one line):").pack(anchor="w", padx=4)
        self.v_rgb_colours = tk.Text(rg, height=5, bg=ENTRY_BG, fg=FG,
                                     font=FONT_MONO, relief="sunken", bd=1)
        self.v_rgb_colours.pack(fill="x", padx=4, pady=2)
        self.v_rgb_colours.insert("1.0", "45,133,186,251,248,183,215,26,29")
        rrg3 = tk.Frame(rg, bg=BG); rrg3.pack(fill="x", padx=4, pady=2)
        lbl(rrg3, "Min val:").pack(side="left")
        self.v_rgb_min = tk.StringVar(value="0.0")
        entry(rrg3, self.v_rgb_min, 7).pack(side="left", padx=2)
        lbl(rrg3, "Max val:").pack(side="left", padx=(8,0))
        self.v_rgb_max = tk.StringVar(value="1000.0")
        entry(rrg3, self.v_rgb_max, 7).pack(side="left", padx=2)
        btn(rrg3, "Convert Grid", self._convert_rgb, w=14).pack(side="left", padx=12)

    # ======================================================================
    # TAB 5 – HELP
    # ======================================================================
    def _tab_help(self, p):
        txt = tk.Text(p, bg=BG, fg=FG, font=FONT, wrap="word",
                      relief="flat", padx=12, pady=8)
        txt.insert("1.0", (
            "SGTool – Geophysical Processing (ArcGIS Pro)\n"
            "=============================================\n\n"
            "INPUT SELECTION\n"
            "  The dropdown at the top lists all raster layers currently loaded\n"
            "  in the active ArcGIS Pro map.  Click ⟳ Refresh to update the list\n"
            "  after adding layers, or Browse File… to use any file directly.\n\n"
            "OUTPUT NAMING\n"
            "  Output files are written to the same folder as the input, with a\n"
            "  suffix describing the operations applied, e.g.:\n"
            "    mag.tif  →  mag_rtp.tif  or  mag_rtp_agc.tif\n\n"
            "FFT BUFFER\n"
            "  The FFT Buffer (px) value is applied as a mirror-padded border\n"
            "  before every FFT operation to suppress edge effects (vignetting).\n"
            "  Increase it for noisy or clipped-edge grids.\n\n"
            "TAB SUMMARY\n"
            "  FFT Filters    – RTP/RTE/Diff.RTP, continuation, band/high/low pass,\n"
            "                   AGC, regional removal, directional Butterworth,\n"
            "                   derivative, THG, tilt angle, analytic signal\n"
            "  Conv + Stats   – mean/median/Gaussian/directional/sun-shading,\n"
            "                   windowed statistics, DTM curvature,\n"
            "                   local anisotropy + chain/streamline length, PCA, ICA,\n"
            "                   Euler deconvolution depth estimation\n"
            "\n"
            "DIFF. RTP (Variable RTP)\n"
            "  Accounts for spatial variation of inclination/declination across the\n"
            "  survey area (Cooper & Cowan 2005 Taylor-series method).  The IGRF is\n"
            "  computed at the four grid corners and the centre using the IGRF Year\n"
            "  field; Inc/Dec need not be entered.  Outputs: _DRTP, _DRTP_inc, _DRTP_dec.\n"
            "\n"
            "LOCAL ANISOTROPY\n"
            "  Structure-tensor analysis of Sobel gradients averaged over the window.\n"
            "  AnisoMag (0–1): saliency = λ1−λ2 normalised to 99th percentile.\n"
            "  AnisoOrient (0–180°): dominant strike direction.\n"
            "\n"
            "CHAIN LENGTH / STREAMLINE LENGTH\n"
            "  Both require active pixels (AnisoMag ≥ Aniso threshold).\n"
            "  Chain Length: disk-search Union-Find connected components — all pixels\n"
            "    in a connected group receive the same score (component size).\n"
            "  Streamline Length: sub-pixel bilinear trace along the orientation field;\n"
            "    total forward + backward path length in pixels.\n"
            "  Grid + Wavelets– point import (CSV/DAT/XYZ → shapefile),\n"
            "                   BSD Worms\n"
            "  Utils          – threshold→NaN, boundary polygon, normalise,\n"
            "                   RGB→grayscale LUT conversion\n\n"
            "REFERENCES\n"
            "  BSDWorms: Holden et al. (2016)\n"
            "  Euler deconvolution: Melo & Barbosa (2019)\n"
            "  IGRF: pyIGRF / Beggan (BGS)\n"
            "  Fatiando a Terra / NumPy / SciPy\n\n"
            "Author: Mark Jessell  |  GUI: Claude (Anthropic)\n"
        ))
        txt.configure(state="disabled")
        txt.pack(fill="both", expand=True, padx=4, pady=4)

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
        def _wrap():
            try:
                func(*args)
            except Exception as e:
                msg = str(e)
                tb  = traceback.format_exc()
                self.after(0, lambda m=msg, t=tb: messagebox.showerror("Error", f"{m}\n\n{t}"))
                self.after(0, lambda m=msg: self._status(f"Error: {m}"))
        threading.Thread(target=_wrap, daemon=True).start()

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
        self._status("Starting processing…")
        self._run_in_thread(self._do_apply_fft, in_r)

    def _do_apply_fft(self, in_r):
        arr, nodata, ll_x, ll_y, cx, cy, sr = raster_to_numpy(in_r)
        buf  = int(self.v_buf.get() or 10)
        proc = make_processor(cx, cy, buf)
        base, ext = os.path.splitext(in_r)
        if not ext:
            ext = ".tif"
        count = [0]

        def _save(result, suffix):
            out_r = f"{base}_{suffix}{ext}"
            numpy_to_raster(result, out_r, ll_x, ll_y, cx, cy, sr, nodata)
            self._record_pending(out_r)
            count[0] += 1
            self.after(0, self._on_output_written, out_r)

        if self.v_rtp.get():
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
                    _save(np.flipud(varrtp),    "DRTP")
                    _save(np.flipud(incv_deg),  "DRTP_inc")
                    _save(np.flipud(decv_deg),  "DRTP_dec")
            elif rt == "Pole":
                self.after(0, lambda: self._status("RTP…"))
                _save(proc.reduction_to_pole(arr, inc, dec, buffer_size=buf), "rtp")
            else:
                self.after(0, lambda: self._status("RTE…"))
                _save(proc.reduction_to_equator(arr, inc, dec, buffer_size=buf), "rte")

        if self.v_vint.get():
            self.after(0, lambda: self._status("Vertical integration…"))
            _save(proc.vertical_integration(arr, buffer_size=buf), "vi")

        if self.v_cont.get():
            h = float(self.v_cont_h.get())
            self.after(0, lambda: self._status(f"Continuation {h}…"))
            if self.v_cont_dir.get() == "up":
                _save(proc.upward_continuation(arr, h, buffer_size=buf), "cont_up")
            else:
                _save(proc.downward_continuation(arr, h, buffer_size=buf), "cont_dn")

        if self.v_dirclean.get():
            self.after(0, lambda: self._status("Directional Butterworth…"))
            scale  = float(self.v_dc_scale.get() or "1.0")
            result = proc.directional_butterworth_band_pass(
                arr,
                1e-8, float(self.v_dc_wl.get()),
                direction_angle=float(self.v_dc_az.get()),
                direction_width=20,
                buffer_size=buf)
            _save(result * scale, "dir")

        if self.v_regrem.get():
            self.after(0, lambda: self._status("Remove regional…"))
            mask = ~np.isnan(arr)
            if self.v_regorder.get() == "1st":
                _save(proc.remove_gradient(arr, mask), "reg1")
            else:
                _save(proc.remove_2o_gradient(arr, mask), "reg2")

        if self.v_bp.get():
            self.after(0, lambda: self._status("Band-pass filter…"))
            w = float(self.v_bp_w.get()) or None
            _save(proc.band_pass_filter(
                arr, float(self.v_bp_lo.get()), float(self.v_bp_hi.get()),
                w, w, buffer_size=buf), "bp")

        if self.v_hlp.get():
            self.after(0, lambda: self._status("High/Low-pass filter…"))
            cut = float(self.v_hlp_cut.get())
            wid = float(self.v_hlp_w.get())
            if self.v_hlp_type.get() == "High":
                _save(proc.high_pass_filter(arr, cut, wid, buffer_size=buf), "hp")
            else:
                _save(proc.low_pass_filter(arr, cut, wid, buffer_size=buf), "lp")

        if self.v_agc.get():
            self.after(0, lambda: self._status("AGC…"))
            _save(proc.automatic_gain_control(arr, self._ensure_odd(self.v_agc_win.get())), "agc")

        if self.v_deriv.get():
            self.after(0, lambda: self._status("Derivative…"))
            d = self.v_deriv_dir.get()
            _save(proc.compute_derivative(
                arr, direction=d,
                order=float(self.v_deriv_pow.get()), buffer_size=buf), f"d{d}")

        if self.v_thg.get():
            self.after(0, lambda: self._status("Total Horizontal Gradient…"))
            _save(proc.total_hz_grad(arr, buffer_size=buf), "thg")

        if self.v_ta.get():
            self.after(0, lambda: self._status("Tilt angle…"))
            _save(np.degrees(proc.tilt_angle(arr, buffer_size=buf)), "ta")

        if self.v_as_.get():
            self.after(0, lambda: self._status("Analytic signal…"))
            _save(proc.analytic_signal(arr, buffer_size=buf), "as")

        self.after(0, lambda: self._status(f"Done — {count[0]} output(s) written."))
        self.after(0, self._uncheck_fft)

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
        self._run_in_thread(self._do_apply_conv_stats, in_r)

    def _do_apply_conv_stats(self, in_r):
        arr, nodata, ll_x, ll_y, cx, cy, sr = raster_to_numpy(in_r)
        from calcs.ConvolutionFilter import ConvolutionFilter
        from calcs.SpatialStats      import SpatialStats
        base, ext = os.path.splitext(in_r)
        if not ext:
            ext = ".tif"
        count = [0]

        def _save(result, suffix):
            out_r = f"{base}_{suffix}{ext}"
            numpy_to_raster(result, out_r, ll_x, ll_y, cx, cy, sr, nodata)
            self._record_pending(out_r)
            count[0] += 1
            self.after(0, self._on_output_written, out_r)

        if self.v_mean.get():
            self.after(0, lambda: self._status("Mean filter…"))
            n = self._ensure_odd(self.v_mean_n.get())
            _save(ConvolutionFilter(arr).mean_filter(n), f"mean{n}")

        if self.v_med.get():
            self.after(0, lambda: self._status("Median filter…"))
            n = self._ensure_odd(self.v_med_n.get())
            _save(ConvolutionFilter(arr).median_filter(n), f"med{n}")

        if self.v_gauss.get():
            self.after(0, lambda: self._status("Gaussian filter…"))
            _save(ConvolutionFilter(arr).gaussian_filter(float(self.v_gauss_s.get())), "gauss")

        if self.v_dconv.get():
            self.after(0, lambda: self._status("Directional filter…"))
            d = self.v_dconv_dir.get().lower()
            _save(ConvolutionFilter(arr).directional_filter(
                self.v_dconv_dir.get(), self._ensure_odd(self.v_dconv_n.get())), f"dir{d}")

        if self.v_sun.get():
            self.after(0, lambda: self._status("Sun shading…"))
            alt = float(self.v_sun_alt.get())
            az  = float(self.v_sun_az.get())
            cf  = ConvolutionFilter(arr)
            result = (cf.sun_shading_filter_grass(arr, cy, cx, altitude=alt, azimuth=az)
                      if self.v_sun_relief.get()
                      else cf.sun_shading_filter(arr, alt, az))
            _save(result, "shade")

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
                self.after(0, lambda s=stat: self._status(f"Spatial {s}…"))
                _save(SpatialStats(arr).calculate_windowed_stats(win, stat), sfx)

        if self.v_dtm.get():
            self.after(0, lambda: self._status("DTM curvature classification…"))
            _save(SpatialStats(arr).classify_terrain_with_cell_size(
                cx, cy,
                float(self.v_dtm_curv.get()), float(self.v_dtm_slope.get()),
                self._ensure_odd(self.v_ss_win.get() or 3), float(self.v_dtm_sigma.get()),
            ).astype(np.float64), "dtm")

        if self.v_ss_aniso.get() or self.v_ss_chain.get() or self.v_ss_stream.get():
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
                _save(aniso_map,  "SS_AnisoMag")
                _save(orient_map, "SS_AnisoOrient")
            if self.v_ss_chain.get():
                self.after(0, lambda: self._status("Chain length…"))
                _save(ss_inst.calculate_anisotropy_chain_length(
                    aniso_map, orient_map,
                    aniso_threshold=a_thresh,
                    angle_tolerance=ang_tol,
                    search_radius=search_r,
                ), "SS_ChainLen")
            if self.v_ss_stream.get():
                self.after(0, lambda: self._status("Streamline length…"))
                _save(ss_inst.calculate_anisotropy_streamline_length(
                    aniso_map, orient_map,
                    aniso_threshold=a_thresh,
                    angle_tolerance=ang_tol,
                    max_steps=max_steps,
                ), "SS_StreamLen")

        if self.v_pca.get():
            self.after(0, lambda: self._status("PCA…"))
            from calcs.PCAICA import PCAICA
            n = int(self.v_pca_n.get() or 0)
            PCAICA(np.zeros((2, 2))).pca_with_nans(
                in_r, f"{base}_pca.tif", n_components=n if n > 0 else None)
            self._record_pending(f"{base}_pca.tif")
            self.after(0, self._on_output_written, f"{base}_pca.tif")
            count[0] += 1

        if self.v_ica.get():
            self.after(0, lambda: self._status("ICA…"))
            from calcs.PCAICA import PCAICA
            n = int(self.v_ica_n.get() or 0)
            PCAICA(np.zeros((2, 2))).ica_with_nans(
                in_r, f"{base}_ica.tif", n_components=n if n > 0 else None)
            self._record_pending(f"{base}_ica.tif")
            self.after(0, self._on_output_written, f"{base}_ica.tif")
            count[0] += 1

        if self.v_euler.get():
            self.after(0, lambda: self._status("Euler deconvolution…"))
            self._do_euler(arr, nodata, ll_x, ll_y, cx, cy, in_r)
            count[0] += 1

        self.after(0, lambda: self._status(f"Done — {count[0]} output(s) written."))
        self.after(0, self._uncheck_conv)

    def _uncheck_conv(self):
        for v in (self.v_mean, self.v_med, self.v_gauss, self.v_dconv,
                  self.v_sun, self.v_sun_relief,
                  self.v_ss_min, self.v_ss_max, self.v_ss_var, self.v_ss_std,
                  self.v_ss_skew, self.v_ss_kurt,
                  self.v_ss_aniso, self.v_ss_chain, self.v_ss_stream,
                  self.v_dtm, self.v_pca, self.v_ica, self.v_euler):
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

        # Returns array shape (n_keep, 4): columns [X, Y, Depth, B]
        res = euler_deconv_optimized(
            filled, XI, YI, ZI, filled.shape, area,
            float(self.v_eu_si.get()),
            int(self.v_eu_win.get()),
            float(self.v_eu_filt.get()))
        if res is None or len(res) == 0:
            self.after(0, lambda: self._status("Euler: no estimates produced."))
            return

        out_csv = os.path.splitext(in_r)[0] + "_euler.csv"
        with open(out_csv, "w", newline="") as f:
            wr = csv.writer(f)
            wr.writerow(["X", "Y", "Depth", "B"])
            for row in res:
                wr.writerow([row[0], row[1], row[2], row[3]])
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
        out_r = os.path.splitext(in_r)[0] + "_nan.tif"
        self._run_in_thread(self._do_threshold, in_r, out_r)

    def _do_threshold(self, in_r, out_r):
        self.after(0, lambda: self._status("Threshold to NaN…"))
        arr, nodata, ll_x, ll_y, cx, cy, sr = raster_to_numpy(in_r)
        from calcs.SG_Util import SG_Util
        result = SG_Util(arr).Threshold2Nan(
            arr, self.v_thresh_cond.get(),
            float(self.v_thresh_above.get()), float(self.v_thresh_below.get()))
        numpy_to_raster(result, out_r, ll_x, ll_y, cx, cy, sr, nodata)
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
    def _calc_igrf(self):
        try:
            from arcgis_utils import calc_igrf, raster_center_latlon

            # Auto-populate Lat/Lon from the selected raster centre
            in_r = self._resolve_input()
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
        from arcgis_utils import calc_igrf
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
            self.after(0, lambda: self._status(f"Power spectrum error: {e}"))


# ============================================================================
def main():
    app = SGToolApp()
    app.mainloop()

if __name__ == "__main__":
    main()
