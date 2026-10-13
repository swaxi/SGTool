"""
External preview window for the ArcGIS Pro GUI.

ArcGIS Pro offers no hook into its map, so the preview is shown in its own window
(like the colour picker of the RGB tool): the filtered grid, drawn with the usual
mean +/- 2 standard deviation grey stretch, that you can zoom with the mouse wheel and
pan by dragging. The GUI recalculates the filter for what is in view (or for the whole
grid, reduced) whenever a setting changes.

This module only draws and handles the mouse; SGToolApp (SGTool_GUI.py) does the
calculation and tells the window what to show.
"""

import warnings
import tkinter as tk

import numpy as np

BG = "#f0f0f0"
FG = "#1c1c1c"
FONT = ("Segoe UI", 9)
NODATA_GREY = 217


def block_mean(arr, k):
    """arr reduced k-fold by averaging k x k blocks (ignoring NaN), as the QGIS plugin
    does for its subsampled preview. k <= 1 returns arr unchanged."""
    if k <= 1:
        return arr
    rows, cols = (arr.shape[0] // k) * k, (arr.shape[1] // k) * k
    if rows == 0 or cols == 0:
        return arr
    blocks = arr[:rows, :cols].reshape(rows // k, k, cols // k, k)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN blocks stay NaN
        return np.nanmean(blocks, axis=(1, 3))


def stretch_to_grey(a, std_clip=True):
    """8-bit RGB array for a grid: grey from mean +/- 2 std (or min to max), no data light grey."""
    finite = np.isfinite(a)
    if finite.any():
        values = a[finite]
        if std_clip:
            lo = float(values.mean() - 2.0 * values.std())
            hi = float(values.mean() + 2.0 * values.std())
        else:
            lo, hi = float(values.min()), float(values.max())
    else:
        lo, hi = 0.0, 1.0
    if hi <= lo:
        hi = lo + 1.0
    grey = np.zeros(a.shape, dtype=np.uint8)
    grey[finite] = np.clip((a[finite] - lo) / (hi - lo) * 255.0, 0, 255).astype(np.uint8)
    rgb = np.dstack([grey, grey, grey])
    rgb[~finite] = NODATA_GREY
    return rgb


class PreviewWindow(tk.Toplevel):
    """The window. Callbacks: on_view_changed() after the view moved or the mode was
    switched, on_close() when the user closes the window."""

    ZOOM_STEP = 1.25
    MAX_PIXELS_PER_CELL = 100.0   # deepest zoom: 100 canvas pixels per cell

    def __init__(self, app, on_view_changed, on_close, on_match_map=None):
        super().__init__(app)
        self.title("SGTool preview")
        self.configure(bg=BG)
        self.geometry("780x680")
        self.minsize(360, 300)
        self._on_view_changed = on_view_changed
        self._on_match_map = on_match_map
        self.protocol("WM_DELETE_WINDOW", on_close)

        bar = tk.Frame(self, bg=BG)
        bar.pack(fill="x", padx=6, pady=(6, 2))
        self.v_mode = tk.StringVar(value="view")
        self.v_original = tk.BooleanVar(value=False)
        for text, value, tip in (
                ("View extent", "view",
                 "Recalculate at full detail for the part of the grid in view\n"
                 "(zoom and pan to update; edges show filter edge effects)"),
                ("Whole grid (reduced)", "whole",
                 "Preview the whole grid, block-averaged to a coarser cell size.\n"
                 "Fast, but filters set in pixels act on the coarser cells, and line\n"
                 "noise is lost if the cell size exceeds the line spacing")):
            rb = tk.Radiobutton(bar, text=text, value=value, variable=self.v_mode,
                                command=self._mode_changed, bg=BG, fg=FG,
                                activebackground=BG, selectcolor="#ffffff", font=FONT)
            rb.pack(side="left", padx=4)
        tk.Checkbutton(bar, text="Show original", variable=self.v_original,
                       command=self._redraw, bg=BG, fg=FG, activebackground=BG,
                       selectcolor="#ffffff", font=FONT).pack(side="left", padx=12)
        tk.Label(bar, text="wheel: zoom   drag: pan   double-click: whole grid",
                 bg=BG, fg="#666666", font=FONT).pack(side="right", padx=4)
        if on_match_map is not None:
            tk.Button(bar, text="Match map view", command=on_match_map, font=FONT, bg="#d0d0d0",
                      fg=FG, relief="raised", bd=1).pack(side="left", padx=(6, 0))

        self.canvas = tk.Canvas(self, bg="#d9d9d9", highlightthickness=0, cursor="fleur")
        self.canvas.pack(fill="both", expand=True, padx=6, pady=2)
        self.v_status = tk.StringVar(value="")
        tk.Label(self, textvariable=self.v_status, bg=BG, fg="#0060b0", font=FONT,
                 anchor="w").pack(fill="x", padx=8, pady=(0, 6))

        self.rows = self.cols = 0
        self.ratio = 1.0              # height of a cell relative to its width (map units)
        self.center = (0.0, 0.0)      # view centre, in cells (x = column, y = row)
        self.cpp = 1.0                # cells per canvas pixel
        self._result = None           # what is shown
        self._fit_box = None          # a box kept in view while the window is resized
        self._photo = None
        self._drag = None
        self._message = "Tick a filter to preview"

        c = self.canvas
        c.bind("<Configure>", lambda e: self._resized())
        c.bind("<MouseWheel>", self._wheel)
        c.bind("<ButtonPress-1>", self._press)
        c.bind("<B1-Motion>", self._motion)
        c.bind("<ButtonRelease-1>", self._release)
        c.bind("<Double-Button-1>", lambda e: self.reset_view())

    # ------------------------------------------------------------------ geometry
    def _canvas_size(self):
        return max(self.canvas.winfo_width(), 50), max(self.canvas.winfo_height(), 50)

    def set_grid(self, rows, cols, force=False, ratio=1.0):
        """Tell the window how big the grid is (cells, and cell height / cell width so that
        cells that are not square are drawn true); the view is reset if it changed (or force,
        for another grid of the same size)."""
        if force or (rows, cols, ratio) != (self.rows, self.cols, self.ratio):
            self.rows, self.cols, self.ratio = rows, cols, ratio
            self._result = None
            self.reset_view()

    def reset_view(self):
        if not self.rows:
            return
        self._fit_box = None
        cw, ch = self._canvas_size()
        self.cpp = max(self.cols / cw, self.rows * self.ratio / ch)
        self.center = (self.cols / 2.0, self.rows / 2.0)
        self._redraw()
        self._on_view_changed()

    def match_canvas_aspect(self, aspect):
        """Resize the window so the drawing area has the shape (width / height) of the
        ArcGIS Pro map canvas, keeping about the same area (and fitting on the screen)."""
        self.update_idletasks()
        cw, ch = self.canvas.winfo_width(), self.canvas.winfo_height()
        if cw < 50 or ch < 50 or not aspect or aspect <= 0:
            return
        extra_w, extra_h = self.winfo_width() - cw, self.winfo_height() - ch
        width = (cw * ch * aspect) ** 0.5
        height = width / aspect
        shrink = min(1.0, (self.winfo_screenwidth() * 0.85 - extra_w) / width,
                     (self.winfo_screenheight() * 0.85 - extra_h) / height)
        width, height = max(width * shrink, 200), max(height * shrink, 150)
        self.geometry(f"{int(width + extra_w)}x{int(height + extra_h)}")
        self.update_idletasks()

    def set_view(self, box):
        """Show the part of the grid in box = (x0, y0, x1, y1), in cells (all of it fits)."""
        if not self.rows:
            return
        self._fit_box = box
        self._fit()
        self._redraw()
        self._on_view_changed()

    def _fit(self):
        """Fit the kept box into the canvas (it stays fitted if the window is resized, until
        the user zooms or pans)."""
        if self._fit_box is None or not self.rows:
            return
        x0, y0, x1, y1 = self._fit_box
        cw, ch = self._canvas_size()
        self.cpp = max((x1 - x0) / cw, (y1 - y0) * self.ratio / ch)
        self.center = ((x0 + x1) / 2.0, (y0 + y1) / 2.0)
        self._clamp()

    def _visible(self):
        cw, ch = self._canvas_size()
        half_w, half_h = cw * self.cpp / 2.0, ch * self.cpp / (2.0 * self.ratio)
        return (self.center[0] - half_w, self.center[1] - half_h,
                self.center[0] + half_w, self.center[1] + half_h)

    def view_box(self):
        """The part of the grid in view as (x0, y0, x1, y1) in cells, clipped to the grid."""
        x0, y0, x1, y1 = self._visible()
        return (max(0.0, x0), max(0.0, y0), min(float(self.cols), x1), min(float(self.rows), y1))

    def _clamp(self):
        cw, ch = self._canvas_size()
        fit = max(self.cols / cw, self.rows * self.ratio / ch)
        self.cpp = min(max(self.cpp, 1.0 / self.MAX_PIXELS_PER_CELL), fit * 1.5)
        self.center = (min(max(self.center[0], 0.0), float(self.cols)),
                       min(max(self.center[1], 0.0), float(self.rows)))

    # ------------------------------------------------------------------ mouse
    def zoom(self, factor, px, py):
        """Zoom by factor about canvas point (px, py), keeping the cell under it fixed."""
        if not self.rows:
            return
        self._fit_box = None
        cw, ch = self._canvas_size()
        cell_x = self.center[0] + (px - cw / 2.0) * self.cpp
        cell_y = self.center[1] + (py - ch / 2.0) * self.cpp / self.ratio
        self.cpp /= factor
        self._clamp()
        self.center = (cell_x - (px - cw / 2.0) * self.cpp,
                       cell_y - (py - ch / 2.0) * self.cpp / self.ratio)
        self._clamp()
        self._redraw()
        self._on_view_changed()

    def _wheel(self, event):
        self.zoom(self.ZOOM_STEP if event.delta > 0 else 1.0 / self.ZOOM_STEP, event.x, event.y)
        return "break"  # the main window's own wheel binding must not scroll as well

    def _press(self, event):
        self._fit_box = None
        self._drag = (event.x, event.y, self.center)

    def _motion(self, event):
        if self._drag is None:
            return
        x0, y0, (cx, cy) = self._drag
        self.center = (cx - (event.x - x0) * self.cpp, cy - (event.y - y0) * self.cpp / self.ratio)
        self._clamp()
        self._redraw()

    def _release(self, event):
        if self._drag is not None:
            moved = abs(event.x - self._drag[0]) + abs(event.y - self._drag[1]) > 2
            self._drag = None
            if moved:
                self._on_view_changed()

    def _resized(self):
        self._fit()
        self._redraw()
        self._on_view_changed()

    def _mode_changed(self):
        self._on_view_changed()

    # ------------------------------------------------------------------ content
    def set_status(self, text):
        self.v_status.set(text)

    def clear(self, message):
        """Show a message instead of an image."""
        self._result = None
        self._message = message
        self._redraw()

    def show_result(self, filtered, original, box, std_clip=True):
        """Show a filtered grid (and the input it came from) covering box = (x0, y0, x1, y1)
        in cells of the full grid."""
        from PIL import Image

        self._result = {
            "box": box,
            "filtered": Image.fromarray(stretch_to_grey(filtered, std_clip)),
            "original": Image.fromarray(stretch_to_grey(original, True)) if original is not None else None,
        }
        self._redraw()

    def _redraw(self):
        c = self.canvas
        c.delete("all")
        cw, ch = self._canvas_size()
        res = self._result
        if res is None:
            c.create_text(cw / 2, ch / 2, text=self._message, fill="#555555", font=("Segoe UI", 11))
            return
        from PIL import Image, ImageTk

        img = res["original"] if (self.v_original.get() and res["original"] is not None) else res["filtered"]
        vx0, vy0, _vx1, _vy1 = self._visible()
        bx0, by0, bx1, by1 = res["box"]
        # where the image sits on the canvas (it may be bigger or smaller than the canvas)
        left, top = (bx0 - vx0) / self.cpp, (by0 - vy0) * self.ratio / self.cpp
        width, height = (bx1 - bx0) / self.cpp, (by1 - by0) * self.ratio / self.cpp
        if width < 1 or height < 1:
            return
        ix0, iy0 = max(0.0, left), max(0.0, top)
        ix1, iy1 = min(float(cw), left + width), min(float(ch), top + height)
        if ix1 <= ix0 or iy1 <= iy0:
            return
        # only the visible part of the image is cut out and resized
        sx0, sx1 = (ix0 - left) / width * img.width, (ix1 - left) / width * img.width
        sy0, sy1 = (iy0 - top) / height * img.height, (iy1 - top) / height * img.height
        piece = img.crop((int(sx0), int(sy0), max(int(sx0) + 1, int(round(sx1))),
                          max(int(sy0) + 1, int(round(sy1)))))
        piece = piece.resize((max(1, int(round(ix1 - ix0))), max(1, int(round(iy1 - iy0)))),
                             Image.BILINEAR)
        self._photo = ImageTk.PhotoImage(piece)  # keep a reference
        c.create_image(ix0, iy0, anchor="nw", image=self._photo)
