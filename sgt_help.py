# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2024-2026 Mark Jessell
"""
Ctrl+click (Cmd+click on macOS too) on a filter's tickbox or button opens its
section of the SGTool reference (docs/reference.md) on GitHub.

The docs are not in the plugin zip (export-ignore), so the guide is read
online, at the heading's anchor (sgt_reference.py, shared with the ArcGIS Pro GUI).
"""

from qgis.PyQt.QtCore import QEvent, QObject, Qt, QTimer, QUrl
from qgis.PyQt.QtGui import QCursor, QDesktopServices
from qgis.PyQt.QtWidgets import QApplication

_CTRL_MODS = Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.MetaModifier

from .sgt_reference import HINT, REFERENCE_URL, anchor  # noqa: F401 (re-exported)
from .sgt_reference import section_url as _section_url

# dialog attribute -> heading in docs/reference.md (the exact heading text, or
# a function of the dialog giving it)
SECTIONS = {
    # FFT Filters+
    "checkBox_4_RTE_P": "RTE/P (reduction to the pole or equator)",
    "checkBox_RemIGRF": "Remove IGRF",
    "checkBox_4_PGrav": "Vertical Integration ●",
    "checkBox_9_continuation": lambda d: (
        "Variable continuation ●"
        if d.comboBox_2_continuationDirection.currentText() == "variable"
        else "Continuation ●"
    ),
    "widget_VarCont": "Variable continuation ●",
    "checkBox_3_DirClean": "Direction Cos/Butterw. (line-noise removal)",
    "pushButton_rad_power_spectrum": "Radial Power Spectrum",
    "checkBox_5_regional": "Remove Regional",
    "checkBox_10_bandPass": "Band Pass ●",
    "checkBox_10_freqCut": "High/Low Pass ●",
    "checkBox_11_1vd_agc": "Automatic Gain Control ●",
    "checkBox_6_derivative": "Derivative ●",
    "checkBox_11_tot_hz_grad": "Total Horizontal Gradient ●",
    "checkBox_7_tiltDerivative": "Tilt Angle ●",
    "checkBox_8_analyticSignal": "Analytic Signal",
    "checkBox_NSTD": "Normalised Standard Deviation (NSTD) ●",
    "checkBox_NAVEdge": "NAV-Edge ●",
    "checkBox_ED": "Euler Deconvolution",
    "checkBox_ED_Stats": "Euler Deconvolution",
    # Conv + Stats
    "checkBox_Mean": "Mean",
    "checkBox_Median": "Median",
    "checkBox_Gaussian": "Gaussian",
    "checkBox_Directional": "Directional",
    "checkBox_SunShading": "Sun Shading",
    "checkBox_SNN": "Symmetric Nearest Neighbour (SNN)",
    "checkBox_Kuwahara": "Kuwahara",
    "checkBox_SS_Min": "Spatial Statistics",
    "checkBox_SS_Max": "Spatial Statistics",
    "checkBox_SS_Variance": "Spatial Statistics",
    "checkBox_SS_StdDev": "Spatial Statistics",
    "checkBox_SS_Skewness": "Spatial Statistics",
    "checkBox_SS_Kurtosis": "Spatial Statistics",
    "checkBox_SS_Anisotropy": "Local Anisotropy",
    "checkBox_SS_ChainLength": "Local Anisotropy",
    "checkBox_SS_Streamline": "Local Anisotropy",
    "checkBox_MRVBF": "MRVBF / MRRTF / Slope",
    # Multiband
    "checkBox_Efimov": "Efimov F-parameter",
    "checkBox_Kd": "Kd and Ud",
    "checkBox_Ud": "Kd and Ud",
    "checkBox_TernGrad": "Ternary Gradient",
    "checkBox_LogRatio": "Multivariate Statistical Analysis",
    "checkBox_PCA": "Multivariate Statistical Analysis",
    "checkBox_ICA": "Multivariate Statistical Analysis",
    # Grid + Wavelets
    "pushButton_load_point_data": "Import point or line data",
    "pushButton_idw_2": "Gridding",
    "pushButton_bspline_3": "Gridding",
    "checkBox_bspline_auto": "Choosing the B-spline settings (Auto)",
    "pushButton_worms": "BSDWorms",
    # Petrophys
    "checkBox_AppSus": "Apparent susceptibility",
    "checkBox_AppDens": "Apparent density",
    # Utils
    "pushButton_save_metadata_xml": "SGTool Metadata (provenance)",
    "checkBox_NaN": "Threshold to NaN",
    "checkBox_polygons": "Create Clipping Polygon",
    "pushButton_clip_to_window": "Clip to Window",
    "pushButton_clip_all_layers": "Clip All Layers",
    "pushButton_normalise": "Normalise Grids",
    "pushButton_3_applyProcessing_Conv_2": "Convert LUT to Grayscale",
}

# controls without a line in SECTIONS, more specific than their group's section
EXTRA = {
    "checkBox_relief": "Sun Shading",
    "checkBox_Kuwahara_generalized": "Kuwahara",
    "radioButton_RR_1st": "Remove Regional",
    "radioButton_RR_2nd": "Remove Regional",
    "pushButton_4_calcIGRF": "RTE/P (reduction to the pole or equator)",
    "checkBox_AppSus_thickness": "Apparent susceptibility",
    "pushButton_AppSus_IGRF": "Apparent susceptibility",
    "checkBox_Henkel_group": "Grouping similar categories",
    "checkBox_Petro_group": "Grouping similar categories",
    "pushButton_replay_history": "Apply the same steps to another grid",
    "checkBox_replace_layers": "Output names",
}

# controls found on every tab (outside the groups), by attribute name prefix
COMMON = (
    ("checkBox_preview", "Live preview"),
    ("pushButton_apply_window", "Live preview"),
    ("pushButton_cancel", "Background calculations"),
    ("pushButton_3_applyProcessing", "Controls common to all tabs"),
    ("pushButton_2_selectGrid", "Controls common to all tabs"),
)

# group titles whose section heading differs (None: no section, not watched)
GROUP_HEADINGS = {
    "SGTool Metadata": "SGTool Metadata (provenance)",
    "Wavelet Transform Modulus Maxima": None,
}


def section_url(heading):
    return QUrl(_section_url(heading))


class _CtrlClickHelp(QObject):
    """Event filter: a Ctrl+click on a watched widget opens its section, and
    is swallowed so a tickbox does not toggle (or a button run)."""

    def __init__(self, dialog, parent=None):
        super().__init__(parent)
        self._dialog = dialog
        self._heading = {}  # widget -> heading or function of the dialog
        self._groups = set()  # group boxes: their title, and controls inside them
        self._bars = set()  # other parents of controls (e.g. the Apply / Cancel bar)
        # a question-mark cursor while Ctrl is held over the dialog. Ctrl's key
        # events go to whichever widget has the focus (often the map), but the
        # application is NOT watched for them: a Python event filter on the
        # QApplication sees every event of every object in QGIS, and wrapping an
        # object that is being deleted crashes QGIS (in SIP, before any Python
        # runs). Instead, while the mouse is over the dialog (its Enter / Leave),
        # a timer polls the keyboard modifiers.
        self._cursor_on = False
        self._poll = QTimer(self)
        self._poll.setInterval(100)
        self._poll.timeout.connect(self._poll_ctrl)
        dialog.installEventFilter(self)  # its Enter / Leave

    def uninstall(self):
        """Stop the Ctrl polling and restore the cursor (e.g. when the plugin
        is unloaded)."""
        self._poll.stop()
        self._set_cursor(False)

    def _over_dialog(self):
        d = self._dialog
        return d.isVisible() and d.rect().contains(d.mapFromGlobal(QCursor.pos()))

    def _set_cursor(self, on):
        if on and not self._cursor_on:
            QApplication.setOverrideCursor(QCursor(Qt.CursorShape.WhatsThisCursor))
            self._cursor_on = True
        elif not on and self._cursor_on:
            QApplication.restoreOverrideCursor()
            self._cursor_on = False

    def _poll_ctrl(self):
        """Question mark while Ctrl is down, the mouse is over the dialog and
        QGIS is the active application; polling stops once the mouse leaves."""
        try:
            over = self._over_dialog()
            held = bool(QApplication.queryKeyboardModifiers() & _CTRL_MODS)
            active = QApplication.activeWindow() is not None
        except RuntimeError:  # the dialog has been deleted
            self._poll.stop()
            self._set_cursor(False)
            return
        self._set_cursor(over and held and active)
        if not over:
            self._poll.stop()

    def watch(self, widget, heading, group=False):
        self._heading[widget] = heading
        if group:
            self._groups.add(widget)
        widget.installEventFilter(self)

    def watch_parent(self, widget):
        """A disabled control gets no mouse events: they go to its parent. So
        the parent of a control outside the groups is watched too, and a Ctrl+
        click there is resolved to the control under the mouse."""
        parent = widget.parentWidget()
        if parent is not None and parent not in self._heading and parent not in self._bars:
            self._bars.add(parent)
            parent.installEventFilter(self)

    def _target(self, obj, event):
        """The watched widget a Ctrl+click on `obj` is for (None: none). A
        click arriving at a group box or a bar is resolved to the watched
        control under the mouse (a disabled one sends its clicks there)."""
        if obj not in self._groups and obj not in self._bars:
            return obj
        pos = event.position().toPoint() if hasattr(event, "position") else event.pos()
        child = obj.childAt(pos)
        while child is not None and child is not obj:
            if child in self._heading:
                return child
            child = child.parentWidget()
        # the title or an empty part of a group: its section; of a bar: nothing
        return obj if obj in self._groups else None

    def eventFilter(self, obj, event):
        kind = event.type()
        if obj is self._dialog and kind in (QEvent.Type.Enter, QEvent.Type.Leave):
            if kind == QEvent.Type.Enter:
                self._poll.start()
                self._poll_ctrl()
            else:
                self._poll.stop()
                self._set_cursor(False)
            return False  # only watched: Enter / Leave go on as usual
        if obj not in self._heading and obj not in self._bars:
            return False
        kind = event.type()
        if kind in (QEvent.Type.MouseButtonPress, QEvent.Type.MouseButtonRelease,
                    QEvent.Type.MouseButtonDblClick):
            # Ctrl is Cmd on macOS; the Ctrl key itself arrives as Meta there,
            # and Qt may turn a Ctrl+click into a right click
            mods = event.modifiers()
            if event.button() in (
                Qt.MouseButton.LeftButton, Qt.MouseButton.RightButton
            ) and mods & (
                Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.MetaModifier
            ):
                target = self._target(obj, event)
                if target is None:
                    return False  # an empty part of a bar: not ours
                if kind == QEvent.Type.MouseButtonRelease:
                    heading = self._heading[target]
                    text = heading(self._dialog) if callable(heading) else heading
                    QDesktopServices.openUrl(section_url(text))
                return True
        return False


def _group_heading(group):
    title = group.title().strip()
    return GROUP_HEADINGS.get(title, title) if title else None


def install(dialog):
    """Watch every group box (Ctrl+click on its title) and every tickbox and
    button of the dialog, each with its section: from SECTIONS, else EXTRA,
    else its group's, else the common controls' (COMMON). Adds the hint to the
    controls' tooltips. Returns the filter (keep a reference to it)."""
    from qgis.PyQt.QtWidgets import QAbstractButton, QGroupBox

    helper = _CtrlClickHelp(dialog, dialog)
    names = {}  # widget -> attribute name of the dialog
    for name, value in vars(dialog).items():
        if isinstance(value, (QAbstractButton, QGroupBox)):
            names.setdefault(value, name)

    for group in dialog.findChildren(QGroupBox):
        heading = _group_heading(group)
        if heading:
            helper.watch(group, heading, group=True)

    # everything SECTIONS names, whatever its type (e.g. the variable
    # continuation panel)
    for name, heading in SECTIONS.items():
        widget = getattr(dialog, name, None)
        if widget is not None and widget not in helper._groups:
            helper.watch(widget, heading)
            tip = widget.toolTip()
            if HINT not in tip:
                widget.setToolTip(tip + ("\n\n" if tip else "") + HINT)

    def enclosing_group(widget):
        parent = widget.parentWidget()
        while parent is not None and parent is not dialog:
            if parent in helper._groups:
                return helper._heading[parent]
            parent = parent.parentWidget()
        return None

    for widget in dialog.findChildren(QAbstractButton):
        if (widget in helper._heading
                or type(widget).__name__ == "QgsGroupBoxCollapseButton"):
            continue  # done above, or a group's fold arrow
        name = names.get(widget, widget.objectName())
        heading = SECTIONS.get(name) or EXTRA.get(name) or enclosing_group(widget)
        if heading is None:
            heading = next((h for prefix, h in COMMON if name.startswith(prefix)), None)
        if heading is None:
            continue  # e.g. a widget the dialog builds without a name
        helper.watch(widget, heading)
        if enclosing_group(widget) is None:
            helper.watch_parent(widget)  # for its clicks while disabled
        tip = widget.toolTip()
        if HINT not in tip:
            widget.setToolTip(tip + ("\n\n" if tip else "") + HINT)
    return helper
