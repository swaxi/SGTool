"""
Ctrl+click (Cmd+click on macOS too) on a filter's tickbox or button opens its
section of the SGTool reference (docs/reference.md) on GitHub.

The docs are not in the plugin zip (export-ignore), so the guide is read
online, at the heading's anchor (sgt_reference.py, shared with the ArcGIS Pro GUI).
"""

from qgis.PyQt.QtCore import QEvent, QObject, Qt, QUrl
from qgis.PyQt.QtGui import QDesktopServices

from .sgt_reference import HINT, REFERENCE_URL, anchor  # noqa: F401 (re-exported)
from .sgt_reference import section_url as _section_url

# dialog attribute -> heading in docs/reference.md (the exact heading text, or
# a function of the dialog giving it)
SECTIONS = {
    # FFT Filters+
    "checkBox_4_RTE_P": "RTE/P (reduction to the pole or equator)",
    "pushButton_removeIGRF": "Remove IGRF",
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

def section_url(heading):
    return QUrl(_section_url(heading))


class _CtrlClickHelp(QObject):
    """Event filter: a Ctrl+click on a watched widget opens its section, and
    is swallowed so a tickbox does not toggle (or a button run)."""

    def __init__(self, dialog, parent=None):
        super().__init__(parent)
        self._dialog = dialog
        self._heading = {}  # widget -> heading or function of the dialog

    def watch(self, widget, heading):
        self._heading[widget] = heading
        widget.installEventFilter(self)

    def eventFilter(self, obj, event):
        heading = self._heading.get(obj)
        if heading is None:
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
                if kind == QEvent.Type.MouseButtonRelease:
                    text = heading(self._dialog) if callable(heading) else heading
                    QDesktopServices.openUrl(section_url(text))
                return True
        return False


def install(dialog):
    """Watch every widget in SECTIONS that the dialog has, and add the hint
    to its tooltip. Returns the filter (keep a reference to it)."""
    helper = _CtrlClickHelp(dialog, dialog)
    for name, heading in SECTIONS.items():
        widget = getattr(dialog, name, None)
        if widget is None:
            continue
        helper.watch(widget, heading)
        tip = widget.toolTip()
        if HINT not in tip:
            widget.setToolTip(tip + ("\n\n" if tip else "") + HINT)
    return helper
