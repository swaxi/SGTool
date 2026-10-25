# -*- coding: utf-8 -*-
from qgis.PyQt import QtWidgets
from qgis.PyQt.QtCore import pyqtSignal, Qt, QDate, QCoreApplication
from qgis.PyQt.QtWidgets import (
    QDockWidget, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QScrollArea, QTabWidget, QPushButton, QLabel, QLineEdit,
    QCheckBox, QRadioButton, QComboBox, QGroupBox, QDateEdit,
    QTextEdit, QTextBrowser, QSpinBox, QDoubleSpinBox, QFrame,
    QButtonGroup,
)
from qgis.gui import (
    QgsMapLayerComboBox, QgsProjectionSelectionWidget,
    QgsFieldComboBox, QgsSpinBox, QgsDoubleSpinBox,
)

# PyQt6 removed flat enum aliases; restore them so call sites need no changes.
if not hasattr(Qt, 'AlignRight'):
    for _n in ('AlignLeft', 'AlignRight', 'AlignHCenter', 'AlignVCenter',
               'AlignCenter', 'AlignTop', 'AlignBottom'):
        setattr(Qt, _n, getattr(Qt.AlignmentFlag, _n))
if not hasattr(Qt, 'ScrollBarAlwaysOff'):
    for _n in ('ScrollBarAlwaysOff', 'ScrollBarAsNeeded', 'ScrollBarAlwaysOn'):
        setattr(Qt, _n, getattr(Qt.ScrollBarPolicy, _n))
if not hasattr(QFrame, 'NoFrame'):
    for _n in ('NoFrame', 'Box', 'Panel', 'StyledPanel', 'HLine', 'VLine', 'WinPanel'):
        setattr(QFrame, _n, getattr(QFrame.Shape, _n))

TAB_STYLE = """
QTabWidget::tab-bar { alignment: left; }
QTabBar::tab {
    min-width: 116px; min-height: 23px;
    background-color: rgb(220,220,220);
    border-radius: 2px;
    border: 1px solid #000000;
    font-family: Arial; font-size: 13px;
}
QTabBar::tab:selected { background-color: rgb(190,190,190); }
QTabBar::tab:hover    { background-color: rgb(229,241,251); }
"""

TAB_BG = "background-color: rgb(240,240,240);"


def _tr(text):
    """Translate text using the SGTool Qt translation context."""
    return QCoreApplication.translate("SGTool", text)


class SGToolDockWidget(QDockWidget):

    closingPlugin = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(_tr("SG Tool"))
        self._build_ui()

    def closeEvent(self, event):
        self.closingPlugin.emit()
        event.accept()

    # ------------------------------------------------------------------
    # Top-level UI assembly
    # ------------------------------------------------------------------
    def _build_ui(self):
        root = QWidget()
        root_layout = QVBoxLayout(root)
        root_layout.setContentsMargins(0, 0, 0, 0)
        root_layout.setSpacing(0)

        self.tabWidget = QTabWidget()
        self.tabWidget.setTabShape(QTabWidget.TabShape.Triangular)
        self.tabWidget.setStyleSheet(TAB_STYLE)
        self.tabWidget.setUsesScrollButtons(False)

        self.tabWidget.addTab(self._tab_fft_filters(),   _tr("FFT Filters"))
        self.tabWidget.addTab(self._tab_conv_stats(),    _tr("Conv + Stats"))
        self.tabWidget.addTab(self._tab_grid_wavelets(), _tr("Grid + Wavelets"))
        self.tabWidget.addTab(self._tab_petrophysics(),  _tr("Petrophys"))
        self.tabWidget.addTab(self._tab_utils(),         _tr("Utils"))
        self.tabWidget.addTab(self._tab_help(),          _tr("Help"))

        root_layout.addWidget(self.tabWidget)

        # Hidden legacy buttons that exist in SGTool.py references
        self.autoinc_pushButton = QPushButton(_tr("Set behaviour"))
        self.autoinc_pushButton.setVisible(False)
        self.pushButton_19 = QPushButton(_tr("RESET THE WINDOW"))
        self.pushButton_19.setVisible(False)
        self.pushButton_22 = QPushButton(_tr("RESET THE WINDOW"))
        self.pushButton_22.setVisible(False)
        self.label_31 = QLabel("© 2024 West African Exploration Initiative.")
        self.label_31.setVisible(False)

        # The whole plugin sits in a scroll area, so the dock can be made much
        # narrower than its contents (a horizontal scroll bar then reaches what
        # no longer fits), leaving more room for the map canvas on small screens.
        shell = QScrollArea()
        shell.setWidgetResizable(True)
        shell.setFrameShape(QFrame.Shape.NoFrame)
        shell.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        shell.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        shell.setWidget(root)
        self.setWidget(shell)

    # ------------------------------------------------------------------
    # Helper: wrap a widget list in a QScrollArea
    # ------------------------------------------------------------------
    @staticmethod
    def _scroll_wrap(content_widget):
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setWidget(content_widget)
        return scroll

    # ------------------------------------------------------------------
    # Helper: Apply button + live-preview controls (fixed bar above the
    # scroll area). `suffix` keeps the widget names unique per tab.
    # ------------------------------------------------------------------
    def _apply_preview_bar(self, apply_button, suffix=""):
        chk = QCheckBox(_tr("Preview"))
        chk.setToolTip(
            _tr(
                "Live preview of one filter in a single temporary layer that "
                "updates as you change parameters.\n"
                "Only one filter can be ticked while previewing, and filters "
                "that can't be previewed are disabled.\n"
                "Unchecking discards the preview; press Apply Processing to "
                "save the full-resolution result as a permanent layer."
            )
        )
        rad_extent = QRadioButton(_tr("Map extent"))
        rad_extent.setToolTip(
            _tr("Preview at full resolution for the current map canvas extent\n"
                "(pan/zoom to update; edges show filter edge effects)")
        )
        rad_extent.setChecked(True)
        rad_sub = QRadioButton(_tr("Subsampled grid"))
        rad_sub.setToolTip(
            _tr("Preview the whole grid, block-averaged to a coarser cell size.\n"
                "Fast, but filters set in pixels (convolution, stats, AGC) act on\n"
                "the coarser cells, and line noise is lost if the cell size\n"
                "exceeds the line spacing.")
        )
        # explicit group so these don't share exclusivity with other radios
        group = QButtonGroup(self)
        group.addButton(rad_extent)
        group.addButton(rad_sub)
        cancel = QPushButton(_tr("Cancel"))
        cancel.setToolTip(
            _tr("Stop the running calculation. It stops at its next checkpoint "
                "(between steps, and inside long loops such as Euler "
                "deconvolution), and nothing from it is saved.")
        )
        cancel.setEnabled(False)  # only while a calculation is running

        setattr(self, "pushButton_cancel" + suffix, cancel)
        setattr(self, "checkBox_preview" + suffix, chk)
        setattr(self, "radioButton_preview_extent" + suffix, rad_extent)
        setattr(self, "radioButton_preview_sub" + suffix, rad_sub)
        setattr(self, "buttonGroup_preview" + suffix, group)

        bar = QWidget()
        bl = QHBoxLayout(bar)
        bl.setContentsMargins(0, 0, 0, 0)
        bl.setSpacing(6)
        bl.addWidget(apply_button)
        bl.addWidget(cancel)
        bl.addSpacing(10)
        bl.addWidget(chk)
        bl.addWidget(rad_extent)
        bl.addWidget(rad_sub)
        bl.addStretch()
        return bar

    # ------------------------------------------------------------------
    # TAB 1 – FFT Filters
    # ------------------------------------------------------------------
    def _tab_fft_filters(self):
        tab = QWidget()
        tab.setObjectName("Import_data")
        tab.setStyleSheet("#Import_data { " + TAB_BG + " }")

        outer = QVBoxLayout(tab)
        outer.setContentsMargins(6, 6, 6, 6)
        outer.setSpacing(4)

        # ── Fixed top: Select Grid + Load Grid ──────────────────────
        top = QWidget()
        tg = QGridLayout(top)
        tg.setContentsMargins(0, 0, 0, 0)
        tg.setSpacing(4)

        self.label_11 = QLabel(_tr("Select Grid"))
        self.label_11.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.mMapLayerComboBox_selectGrid = QgsMapLayerComboBox()

        self.label_38 = QLabel(_tr("Load Grid"))
        self.label_38.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.lineEdit_2_loadGridPath = QLineEdit()
        self.pushButton_2_selectGrid = QPushButton("...")
        self.pushButton_2_selectGrid.setFixedWidth(30)

        tg.addWidget(self.label_11,                    0, 0)
        tg.addWidget(self.mMapLayerComboBox_selectGrid, 0, 1, 1, 2)
        tg.addWidget(self.label_38,                    1, 0)
        tg.addWidget(self.lineEdit_2_loadGridPath,     1, 1)
        tg.addWidget(self.pushButton_2_selectGrid,     1, 2)
        tg.setColumnStretch(1, 1)

        outer.addWidget(top)

        # ── Fixed: Apply Processing + live preview controls ─────────
        self.pushButton_3_applyProcessing = QPushButton(_tr("Apply Processing"))
        self.pushButton_3_applyProcessing.setStyleSheet(
            "font-weight: bold;")
        outer.addWidget(
            self._apply_preview_bar(self.pushButton_3_applyProcessing)
        )

        # ── Scrollable middle ────────────────────────────────────────
        scroll_content = QWidget()
        sl = QVBoxLayout(scroll_content)
        sl.setContentsMargins(2, 2, 2, 2)
        sl.setSpacing(6)
        sl.addWidget(self._group_grav_mag())
        sl.addWidget(self._group_freq_filters())
        sl.addWidget(self._group_gradient())
        sl.addStretch()

        outer.addWidget(self._scroll_wrap(scroll_content), 1)

        # ── Fixed bottom: buffer / version row ──────────────────────
        bot = QWidget()
        bl = QHBoxLayout(bot)
        bl.setContentsMargins(0, 2, 0, 0)
        bl.setSpacing(6)

        self.label_37 = QLabel(_tr("Max FFT Buffer (pixels)"))
        self.label_37.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.lineEdit_13_max_buffer = QLineEdit("5000")
        self.lineEdit_13_max_buffer.setFixedWidth(55)

        self.label_27 = QLabel(_tr("Version"))
        self.label_27.setAlignment(Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignRight)
        self.version_label = QLabel("0.1")
        self.version_label.setAlignment(Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft)

        self.label_41_units = QLabel(_tr("Units"))

        bl.addWidget(self.label_37)
        bl.addWidget(self.lineEdit_13_max_buffer)
        bl.addWidget(self.label_27)
        bl.addWidget(self.version_label)
        bl.addWidget(self.label_41_units)
        bl.addStretch()  # keep the row at the left of the tab
        self._info_labels = [self.label_41_units]
        self._version_labels = [self.version_label]

        outer.addWidget(bot)
        return tab

    def _grid_info_bar(self):
        """The version and grid information row (units, CRS, pixels, pixel size)
        for the bottom of a tab. Each tab needs its own widgets, which the plugin
        updates together with those of the first tab."""
        bar = QWidget()
        bl = QHBoxLayout(bar)
        bl.setContentsMargins(0, 2, 0, 0)
        bl.setSpacing(6)
        version_text = QLabel(_tr("Version"))
        version_text.setAlignment(Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignRight)
        version = QLabel(self._version_labels[0].text())
        version.setAlignment(Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft)
        info = QLabel(self._info_labels[0].text())
        bl.addWidget(version_text)
        bl.addWidget(version)
        bl.addWidget(info)
        bl.addStretch()
        self._version_labels.append(version)
        self._info_labels.append(info)
        return bar

    @staticmethod
    def _bullet_check(bullet, checkbox):
        """The dot-point marker and its checkbox side by side at the left edge, so
        the checkbox lines up with those of filters that have no marker."""
        box = QHBoxLayout()
        box.setContentsMargins(0, 0, 0, 0)
        box.setSpacing(2)
        bullet.setFixedWidth(10)
        box.addWidget(bullet)
        box.addWidget(checkbox)
        box.addStretch()
        return box

    # ------------------------------------------------------------------
    # Group: Grav/Mag Filters
    # ------------------------------------------------------------------
    def _group_grav_mag(self):
        gb = QGroupBox(_tr("Grav/Mag Filters"))
        g = QGridLayout(gb)
        g.setSpacing(4)

        # Row 0: RTE/P + combo + IGRF
        self.checkBox_4_RTE_P = QCheckBox("RTE/P")
        self.comboBox_3_rte_p_list = QComboBox()
        self.pushButton_4_calcIGRF = QPushButton("IGRF")

        g.addWidget(self.checkBox_4_RTE_P,        0, 0)
        g.addWidget(self.comboBox_3_rte_p_list,   0, 1, 1, 2)
        g.addWidget(self.pushButton_4_calcIGRF,   0, 3)

        # Row 1: Inc / Dec / Int / date
        self.label_15 = QLabel(_tr("Inc"))
        self.label_15.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.lineEdit_6_inc = QLineEdit("0")
        self.lineEdit_6_inc.setFixedWidth(45)

        self.label_14 = QLabel(_tr("Dec"))
        self.label_14.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.lineEdit_5_dec = QLineEdit("0")
        self.lineEdit_5_dec.setFixedWidth(45)

        self.label_33 = QLabel(_tr("Int"))
        self.label_33.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.lineEdit_6_int = QLineEdit("50000")
        self.lineEdit_6_int.setFixedWidth(60)

        self.dateEdit = QDateEdit()
        self.dateEdit.setMinimumDate(QDate(1900, 1, 1))
        self.dateEdit.setMaximumDate(QDate(2030, 12, 31))
        self.dateEdit.setCalendarPopup(True)

        row1 = QHBoxLayout()
        row1.addWidget(self.label_15)
        row1.addWidget(self.lineEdit_6_inc)
        row1.addWidget(self.label_14)
        row1.addWidget(self.lineEdit_5_dec)
        row1.addWidget(self.label_33)
        row1.addWidget(self.lineEdit_6_int)
        row1.addWidget(self.dateEdit)
        row1.addStretch()
        g.addLayout(row1, 1, 0, 1, 5)

        # Row 2: ● PGrav
        self.label_26 = QLabel("●")
        self.label_26.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.checkBox_4_PGrav = QCheckBox(
            _tr("Vertical Integration (apply to RTE/P result to get Pseudo Gravity)"))
        g.addLayout(self._bullet_check(self.label_26, self.checkBox_4_PGrav), 2, 0, 1, 5)

        # Row 3: ● Continuation / Direction / Height
        self.label_28 = QLabel("●")
        self.label_28.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.checkBox_9_continuation = QCheckBox(_tr("Continuation"))

        self.label_20 = QLabel(_tr("Direction"))
        self.label_20.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.comboBox_2_continuationDirection = QComboBox()
        self.comboBox_2_continuationDirection.setFixedWidth(70)

        self.label_21 = QLabel(_tr("Height (m)"))
        self.label_21.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.lineEdit_10_continuationHeight = QLineEdit("500")
        self.lineEdit_10_continuationHeight.setFixedWidth(70)

        row3 = QHBoxLayout()
        row3.addWidget(self.label_28)
        row3.addWidget(self.checkBox_9_continuation)
        row3.addWidget(self.label_20)
        row3.addWidget(self.comboBox_2_continuationDirection)
        row3.addWidget(self.label_21)
        row3.addWidget(self.lineEdit_10_continuationHeight)
        row3.addStretch()
        g.addLayout(row3, 3, 0, 1, 5)

        g.setColumnStretch(4, 1)
        return gb

    # ------------------------------------------------------------------
    # Group: Frequency Filters
    # ------------------------------------------------------------------
    def _group_freq_filters(self):
        """One filter per block, its settings on indented rows beneath it, so the
        group stays narrow (the tab has plenty of vertical room)."""
        gb = QGroupBox(_tr("Frequency Filters"))
        v = QVBoxLayout(gb)
        v.setSpacing(4)

        def row(*items, indent=True):
            """A row of widgets (or layouts) in a line, indented under its filter."""
            box = QHBoxLayout()
            box.setContentsMargins(0, 0, 0, 0)
            box.setSpacing(4)
            if indent:
                box.addSpacing(22)
            for item in items:
                if isinstance(item, int):
                    box.addSpacing(item)
                elif hasattr(item, "addWidget"):
                    box.addLayout(item)
                else:
                    box.addWidget(item)
            box.addStretch()
            v.addLayout(box)

        self.pushButton_rad_power_spectrum = QPushButton(_tr("Radial Power Spectrum"))
        row(self.pushButton_rad_power_spectrum, indent=False)

        # Direction filter
        self.checkBox_3_DirClean = QCheckBox(_tr("Direction Cos/Butterw."))
        self.label_12 = QLabel(_tr("Azimuth"))
        self.label_12.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.lineEdit_3_azimuth = QLineEdit("0")
        self.lineEdit_3_azimuth.setFixedWidth(35)
        self.label_13 = QLabel(_tr("Line Spacing min/max (proj units)"))
        self.label_13.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.lineEdit_3_DC_wavelength = QLineEdit("1000")
        self.lineEdit_3_DC_wavelength.setFixedWidth(55)
        self.lineEdit_3_DC_maxspacing = QLineEdit("")
        self.lineEdit_3_DC_maxspacing.setPlaceholderText(_tr("max"))
        self.lineEdit_3_DC_maxspacing.setFixedWidth(55)
        self.label_16 = QLabel(_tr("Scale"))
        self.label_16.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        # (1 is about right since the filter keeps both halves of the spectrum;
        # it used to be 3 to make up for losing half the amplitude)
        self.lineEdit_3_DC_scale = QLineEdit("1")
        self.lineEdit_3_DC_scale.setFixedWidth(35)
        self.label_wedge = QLabel(_tr("Wedge"))
        self.label_wedge.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.lineEdit_3_DC_width = QLineEdit("45")
        self.lineEdit_3_DC_width.setFixedWidth(35)
        self.lineEdit_3_DC_width.setToolTip(
            _tr("Half-width of the directional wedge in degrees: noise whose "
                "direction is within this angle of the azimuth is kept in the "
                "noise estimate, with a cosine taper. Wider keeps more of the "
                "noise (and more of the geology).")
        )

        row(self.checkBox_3_DirClean, indent=False)
        row(self.label_12, self.lineEdit_3_azimuth, 6, self.label_16, self.lineEdit_3_DC_scale,
            6, self.label_wedge, self.lineEdit_3_DC_width)
        row(self.label_13, self.lineEdit_3_DC_wavelength, self.lineEdit_3_DC_maxspacing)

        # Remove Regional
        self.checkBox_5_regional = QCheckBox(_tr("Remove Regional"))
        self.radioButton_RR_1st = QRadioButton(_tr("1st order fit"))
        self.radioButton_RR_1st.setChecked(True)
        self.radioButton_RR_2nd = QRadioButton(_tr("2nd order fit"))
        row(self.checkBox_5_regional, 8, self.radioButton_RR_1st, self.radioButton_RR_2nd,
            indent=False)

        # Band Pass
        self.label_29 = QLabel("●")
        self.label_29.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.checkBox_10_bandPass = QCheckBox(_tr("Band Pass"))
        self.label_22 = QLabel(_tr("Low"))
        self.label_22.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.lineEdit_12_bandPassLow = QLineEdit("50000")
        self.lineEdit_12_bandPassLow.setFixedWidth(60)
        self.label_23 = QLabel(_tr("High (proj units)"))
        self.label_23.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.lineEdit_11_bandPassHigh = QLineEdit("5000")
        self.lineEdit_11_bandPassHigh.setFixedWidth(55)
        self.label_39 = QLabel(_tr("Width"))
        self.label_39.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.lineEdit_3_BP_width = QLineEdit("5000")
        self.lineEdit_3_BP_width.setFixedWidth(55)

        v.addLayout(self._bullet_check(self.label_29, self.checkBox_10_bandPass))
        row(self.label_22, self.lineEdit_12_bandPassLow, 6, self.label_23, self.lineEdit_11_bandPassHigh)
        row(self.label_39, self.lineEdit_3_BP_width)

        # High/Low Pass
        self.label_32 = QLabel("●")
        self.label_32.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.checkBox_10_freqCut = QCheckBox(_tr("High/Low Pass"))
        self.label_35 = QLabel(_tr("Accept"))
        self.label_35.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.comboBox_2_FreqCutType = QComboBox()
        self.comboBox_2_FreqCutType.setFixedWidth(80)
        self.label_34 = QLabel(_tr("Cutoff (proj units)"))
        self.label_34.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.lineEdit_12_FreqPass = QLineEdit("5000")
        self.lineEdit_12_FreqPass.setFixedWidth(55)
        self.label_69 = QLabel(_tr("Width"))
        self.label_69.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.lineEdit_3_HLP_width = QLineEdit("5000")
        self.lineEdit_3_HLP_width.setFixedWidth(55)

        v.addLayout(self._bullet_check(self.label_32, self.checkBox_10_freqCut))
        row(self.label_35, self.comboBox_2_FreqCutType, 6, self.label_34, self.lineEdit_12_FreqPass)
        row(self.label_69, self.lineEdit_3_HLP_width)

        # AGC
        self.label_30 = QLabel("●")
        self.label_30.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.checkBox_11_1vd_agc = QCheckBox(_tr("Automatic Gain Control"))
        self.label_24 = QLabel(_tr("Window (pixels)"))
        self.label_24.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.lineEdit_13_agc_window = QLineEdit("10")
        self.lineEdit_13_agc_window.setFixedWidth(45)

        v.addLayout(self._bullet_check(self.label_30, self.checkBox_11_1vd_agc))
        row(self.label_24, self.lineEdit_13_agc_window)

        return gb

    # ------------------------------------------------------------------
    # Group: Gradient Filters
    # ------------------------------------------------------------------
    def _group_gradient(self):
        gb = QGroupBox(_tr("Gradient Filters"))
        g = QGridLayout(gb)
        g.setSpacing(4)

        # Row 0: Derivative
        self.label_25 = QLabel("●")
        self.label_25.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.checkBox_6_derivative = QCheckBox(_tr("Derivative"))
        self.label_18 = QLabel(_tr("Direction"))
        self.label_18.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.comboBox_derivDirection = QComboBox()
        self.comboBox_derivDirection.setFixedWidth(60)
        self.label_19 = QLabel(_tr("Power"))
        self.label_19.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.lineEdit_9_derivePower = QLineEdit("1")
        self.lineEdit_9_derivePower.setFixedWidth(55)

        g.addWidget(self.label_25,              0, 0)
        g.addWidget(self.checkBox_6_derivative, 0, 1)
        g.addWidget(self.label_18,              0, 2)
        g.addWidget(self.comboBox_derivDirection, 0, 3)
        g.addWidget(self.label_19,              0, 4)
        g.addWidget(self.lineEdit_9_derivePower, 0, 5)

        # Row 1: Total Horizontal Gradient
        self.label_40 = QLabel("●")
        self.label_40.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.checkBox_11_tot_hz_grad = QCheckBox(_tr("Total Horizontal Gradient"))
        g.addWidget(self.label_40,               1, 0)
        g.addWidget(self.checkBox_11_tot_hz_grad, 1, 1, 1, 3)

        # Row 2: Tilt Angle
        self.label_41 = QLabel("●")
        self.label_41.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.checkBox_7_tiltDerivative = QCheckBox(_tr("Tilt Angle"))
        g.addWidget(self.label_41,               2, 0)
        g.addWidget(self.checkBox_7_tiltDerivative, 2, 1, 1, 3)

        # Row 3: Analytic Signal
        self.checkBox_8_analyticSignal = QCheckBox(_tr("Analytic Signal"))
        g.addWidget(self.checkBox_8_analyticSignal, 3, 0, 1, 3)

        g.setColumnStretch(6, 1)
        return gb

    # ------------------------------------------------------------------
    # TAB 2 – Conv + Stats
    # ------------------------------------------------------------------
    def _tab_conv_stats(self):
        tab = QWidget()
        tab.setObjectName("Fieldwork_preparation")
        tab.setStyleSheet("#Fieldwork_preparation { " + TAB_BG + " }")

        outer = QVBoxLayout(tab)
        outer.setContentsMargins(6, 6, 6, 6)
        outer.setSpacing(4)

        # ── Fixed top: Select Grid ───────────────────────────────────
        top = QWidget()
        tl = QHBoxLayout(top)
        tl.setContentsMargins(0, 0, 0, 0)
        self.label_47 = QLabel(_tr("Select Grid"))
        self.label_47.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.mMapLayerComboBox_selectGrid_Conv = QgsMapLayerComboBox()
        tl.addWidget(self.label_47)
        tl.addWidget(self.mMapLayerComboBox_selectGrid_Conv, 1)
        outer.addWidget(top)

        # ── Fixed: Apply Processing + live preview controls ─────────
        self.pushButton_3_applyProcessing_Conv = QPushButton(_tr("Apply Processing"))
        self.pushButton_3_applyProcessing_Conv.setStyleSheet(
            "font-weight: bold;")
        outer.addWidget(
            self._apply_preview_bar(
                self.pushButton_3_applyProcessing_Conv, suffix="_conv"
            )
        )

        # ── Scrollable middle ────────────────────────────────────────
        sc = QWidget()
        sl = QVBoxLayout(sc)
        sl.setContentsMargins(2, 2, 2, 2)
        sl.setSpacing(6)
        sl.addWidget(self._group_conv_filters())
        sl.addWidget(self._group_spatial_stats())
        sl.addWidget(self._group_multivar_stats())
        sl.addWidget(self._group_euler())
        sl.addStretch()
        outer.addWidget(self._scroll_wrap(sc), 1)
        outer.addWidget(self._grid_info_bar())
        return tab

    # ------------------------------------------------------------------
    # Group: Convolution Filters
    # ------------------------------------------------------------------
    def _group_conv_filters(self):
        gb = QGroupBox(_tr("Convolution Filters"))
        g = QGridLayout(gb)
        g.setSpacing(4)

        hdr_lbl = QLabel(_tr("Filter size (odd numbers)"))
        hdr_lbl.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)

        # Mean
        self.checkBox_Mean = QCheckBox(_tr("Mean"))
        self.lineEdit_Mean_size = QLineEdit("3")
        self.lineEdit_Mean_size.setFixedWidth(50)
        g.addWidget(self.checkBox_Mean,       0, 0)
        g.addWidget(hdr_lbl,                  0, 1)
        g.addWidget(self.lineEdit_Mean_size,  0, 2)

        # Median
        self.checkBox_Median = QCheckBox(_tr("Median"))
        lbl_med = QLabel(_tr("Filter size (odd numbers)"))
        lbl_med.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.lineEdit_Median_size = QLineEdit("3")
        self.lineEdit_Median_size.setFixedWidth(50)
        g.addWidget(self.checkBox_Median,     1, 0)
        g.addWidget(lbl_med,                  1, 1)
        g.addWidget(self.lineEdit_Median_size, 1, 2)

        # Gaussian
        self.checkBox_Gaussian = QCheckBox(_tr("Gaussian"))
        self.label_65 = QLabel(_tr("Sigma"))
        self.label_65.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.lineEdit_Gaussian_Sigma = QLineEdit("1")
        self.lineEdit_Gaussian_Sigma.setFixedWidth(50)
        g.addWidget(self.checkBox_Gaussian,   2, 0)
        g.addWidget(self.label_65,            2, 1)
        g.addWidget(self.lineEdit_Gaussian_Sigma, 2, 2)

        # Directional
        self.checkBox_Directional = QCheckBox(_tr("Directional"))
        self.label_36 = QLabel(_tr("Direction"))
        self.label_36.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.comboBox_Dir_dir = QComboBox()
        self.comboBox_Dir_dir.setFixedWidth(70)
        g.addWidget(self.checkBox_Directional, 3, 0)
        g.addWidget(self.label_36,             3, 1)
        g.addWidget(self.comboBox_Dir_dir,     3, 2)

        # Sun Shading
        self.checkBox_SunShading = QCheckBox(_tr("Sun Shading"))
        self.label_64 = QLabel(_tr("Azimuth"))
        self.label_64.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.lineEdit_SunSh_Az = QLineEdit("-45")
        self.lineEdit_SunSh_Az.setFixedWidth(50)
        self.label_68 = QLabel(_tr("Zenith"))
        self.label_68.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.lineEdit_SunSh_Zn = QLineEdit("45")
        self.lineEdit_SunSh_Zn.setFixedWidth(50)
        self.checkBox_relief = QCheckBox(_tr("relief"))
        self.checkBox_relief.setChecked(True)

        g.addWidget(self.checkBox_SunShading, 4, 0)
        g.addWidget(self.label_64,            4, 1)
        g.addWidget(self.lineEdit_SunSh_Az,   4, 2)
        g.addWidget(self.label_68,            4, 3)
        g.addWidget(self.lineEdit_SunSh_Zn,   4, 4)
        g.addWidget(self.checkBox_relief,     4, 5)

        g.setColumnStretch(6, 1)
        return gb

    # ------------------------------------------------------------------
    # Group: Spatial Statistics
    # ------------------------------------------------------------------
    def _group_spatial_stats(self):
        gb = QGroupBox(_tr("Spatial Statistics"))
        g = QGridLayout(gb)
        g.setSpacing(4)

        # Window size row
        lbl_ws = QLabel(_tr("Window size for all filters (odd numbers)"))
        lbl_ws.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.label_78 = lbl_ws
        self.lineEdit_SS_Window = QLineEdit("3")
        self.lineEdit_SS_Window.setFixedWidth(50)
        g.addWidget(lbl_ws,                   0, 0, 1, 3)
        g.addWidget(self.lineEdit_SS_Window,  0, 3)

        # Stats checkboxes
        self.checkBox_SS_Min = QCheckBox(_tr("Min"))
        self.checkBox_SS_Max = QCheckBox(_tr("Max"))
        self.checkBox_SS_Variance = QCheckBox(_tr("Variance"))
        self.checkBox_SS_StdDev = QCheckBox(_tr("Standard Deviation"))
        self.checkBox_SS_Skewness = QCheckBox(_tr("Skewness"))
        self.checkBox_SS_Kurtosis = QCheckBox(_tr("Kurtosis"))
        self.checkBox_SS_Anisotropy = QCheckBox(_tr("Local Anisotropy"))
        self.comboBox_SS_Anisotropy_WinType = QComboBox()
        self.comboBox_SS_Anisotropy_WinType.addItems([_tr("Gaussian"), _tr("Square")])
        self.comboBox_SS_Anisotropy_WinType.setFixedWidth(90)

        g.addWidget(self.checkBox_SS_Min,     1, 0)
        g.addWidget(self.checkBox_SS_Max,     1, 2)
        g.addWidget(self.checkBox_SS_Variance, 2, 0)
        g.addWidget(self.checkBox_SS_StdDev,  2, 2)
        g.addWidget(self.checkBox_SS_Skewness, 3, 0)
        g.addWidget(self.checkBox_SS_Kurtosis, 3, 2)
        g.addWidget(self.checkBox_SS_Anisotropy,         4, 0, 1, 2)
        g.addWidget(self.comboBox_SS_Anisotropy_WinType, 4, 2, 1, 2)

        # Anisotropy connectivity analysis
        self.checkBox_SS_ChainLength  = QCheckBox(_tr("Chain Length"))
        self.checkBox_SS_Streamline   = QCheckBox(_tr("Streamline Length"))

        lbl_aniso_thr = QLabel(_tr("Aniso thresh:"))
        lbl_aniso_thr.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.doubleSpinBox_SS_AnisoThresh = QDoubleSpinBox()
        self.doubleSpinBox_SS_AnisoThresh.setRange(0.0, 1.0)
        self.doubleSpinBox_SS_AnisoThresh.setSingleStep(0.05)
        self.doubleSpinBox_SS_AnisoThresh.setValue(0.3)
        self.doubleSpinBox_SS_AnisoThresh.setFixedWidth(55)

        lbl_ang_tol = QLabel(_tr("Angle tol°:"))
        lbl_ang_tol.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.doubleSpinBox_SS_AngleTol = QDoubleSpinBox()
        self.doubleSpinBox_SS_AngleTol.setRange(1.0, 90.0)
        self.doubleSpinBox_SS_AngleTol.setSingleStep(5.0)
        self.doubleSpinBox_SS_AngleTol.setValue(22.5)
        self.doubleSpinBox_SS_AngleTol.setFixedWidth(55)

        lbl_steps = QLabel(_tr("Max steps:"))
        lbl_steps.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.spinBox_SS_MaxSteps = QSpinBox()
        self.spinBox_SS_MaxSteps.setRange(10, 2000)
        self.spinBox_SS_MaxSteps.setSingleStep(50)
        self.spinBox_SS_MaxSteps.setValue(200)
        self.spinBox_SS_MaxSteps.setFixedWidth(55)

        lbl_search_r = QLabel(_tr("Search r:"))
        lbl_search_r.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.spinBox_SS_SearchRadius = QSpinBox()
        self.spinBox_SS_SearchRadius.setRange(1, 10)
        self.spinBox_SS_SearchRadius.setSingleStep(1)
        self.spinBox_SS_SearchRadius.setValue(2)
        self.spinBox_SS_SearchRadius.setFixedWidth(45)

        g.addWidget(self.checkBox_SS_ChainLength,          5, 0, 1, 2)
        g.addWidget(self.checkBox_SS_Streamline,           5, 2, 1, 2)
        g.addWidget(lbl_aniso_thr,                         6, 0)
        g.addWidget(self.doubleSpinBox_SS_AnisoThresh,     6, 1)
        g.addWidget(lbl_ang_tol,                           6, 2)
        g.addWidget(self.doubleSpinBox_SS_AngleTol,        6, 3)
        g.addWidget(lbl_steps,                             6, 4)
        g.addWidget(self.spinBox_SS_MaxSteps,              6, 5)
        g.addWidget(lbl_search_r,                          6, 6)
        g.addWidget(self.spinBox_SS_SearchRadius,          6, 7)

        # MRVBF / MRRTF (Multiresolution Valley Bottom / Ridge Top Flatness)
        self.checkBox_MRVBF = QCheckBox(_tr("MRVBF / MRRTF / Slope"))

        self.label_mrvbf_tslope = QLabel(_tr("Slope Threshold"))
        self.label_mrvbf_tslope.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.lineEdit_MRVBF_tslope = QLineEdit("16")
        self.lineEdit_MRVBF_tslope.setFixedWidth(55)

        self.label_mrvbf_tpctlv = QLabel(_tr("Pctl Lowness"))
        self.label_mrvbf_tpctlv.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.lineEdit_MRVBF_tpctlv = QLineEdit("0.40")
        self.lineEdit_MRVBF_tpctlv.setFixedWidth(55)

        self.label_mrvbf_tpctlr = QLabel(_tr("Pctl Upness"))
        self.label_mrvbf_tpctlr.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.lineEdit_MRVBF_tpctlr = QLineEdit("0.35")
        self.lineEdit_MRVBF_tpctlr.setFixedWidth(55)

        self.label_mrvbf_pslope = QLabel(_tr("Slope Shape"))
        self.label_mrvbf_pslope.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.lineEdit_MRVBF_pslope = QLineEdit("4")
        self.lineEdit_MRVBF_pslope.setFixedWidth(55)

        self.label_mrvbf_ppctl = QLabel(_tr("Pctl Shape"))
        self.label_mrvbf_ppctl.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.lineEdit_MRVBF_ppctl = QLineEdit("3")
        self.lineEdit_MRVBF_ppctl.setFixedWidth(55)

        self.label_mrvbf_maxres = QLabel(_tr("Max Res %"))
        self.label_mrvbf_maxres.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.lineEdit_MRVBF_maxres = QLineEdit("100")
        self.lineEdit_MRVBF_maxres.setFixedWidth(55)

        g.addWidget(self.checkBox_MRVBF,        7, 0, 1, 2)
        g.addWidget(self.label_mrvbf_tslope,    7, 2)
        g.addWidget(self.lineEdit_MRVBF_tslope, 7, 3)
        g.addWidget(self.label_mrvbf_tpctlv,    7, 4)
        g.addWidget(self.lineEdit_MRVBF_tpctlv, 7, 5)
        g.addWidget(self.label_mrvbf_tpctlr,    8, 2)
        g.addWidget(self.lineEdit_MRVBF_tpctlr, 8, 3)
        g.addWidget(self.label_mrvbf_pslope,    8, 4)
        g.addWidget(self.lineEdit_MRVBF_pslope, 8, 5)
        g.addWidget(self.label_mrvbf_ppctl,     9, 2)
        g.addWidget(self.lineEdit_MRVBF_ppctl,  9, 3)
        g.addWidget(self.label_mrvbf_maxres,    9, 4)
        g.addWidget(self.lineEdit_MRVBF_maxres, 9, 5)

        g.setColumnStretch(8, 1)
        return gb

    # ------------------------------------------------------------------
    # Group: Multivariate Statistical Analysis
    # ------------------------------------------------------------------
    def _group_multivar_stats(self):
        gb = QGroupBox(_tr("Multivariate Statistical Analysis"))
        g = QGridLayout(gb)
        g.setSpacing(4)

        lbl_comp = QLabel(_tr("Number of Components (0=all)"))
        lbl_comp.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        lbl_comp2 = QLabel(_tr("Number of Components (0=all)"))
        lbl_comp2.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)

        self.checkBox_PCA = QCheckBox(_tr("Principal Component Analysis"))
        self.mQgsSpinBox_PCA = QgsSpinBox()
        self.label = lbl_comp

        self.checkBox_ICA = QCheckBox(_tr("Independent Component Analysis"))
        self.mQgsSpinBox_ICA = QgsSpinBox()
        self.label_2 = lbl_comp2

        g.addWidget(self.checkBox_PCA,     0, 0)
        g.addWidget(lbl_comp,              0, 1)
        g.addWidget(self.mQgsSpinBox_PCA,  0, 2)
        g.addWidget(self.checkBox_ICA,     1, 0)
        g.addWidget(lbl_comp2,             1, 1)
        g.addWidget(self.mQgsSpinBox_ICA,  1, 2)

        g.setColumnStretch(3, 1)
        return gb

    # ------------------------------------------------------------------
    # Group: Euler Deconvolution
    # ------------------------------------------------------------------
    def _group_euler(self):
        gb = QGroupBox(_tr("Euler Deconvolution"))
        g = QGridLayout(gb)
        g.setSpacing(4)

        self.checkBox_ED = QCheckBox(_tr("Euler Deconvolution"))
        self.checkBox_ED_Stats = QCheckBox(_tr("Stats"))
        self.label_79 = QLabel(_tr("Threshold"))
        self.label_79.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.doubleSpinBox_ED_Threshold = QDoubleSpinBox()
        self.doubleSpinBox_ED_Threshold.setDecimals(5)
        self.doubleSpinBox_ED_Threshold.setMinimum(0.00001)
        self.doubleSpinBox_ED_Threshold.setMaximum(1.0)
        self.doubleSpinBox_ED_Threshold.setSingleStep(0.1)
        self.doubleSpinBox_ED_Threshold.setValue(0.1)
        self.label_80 = QLabel(_tr("Window (odd)"))
        self.label_80.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.lineEdit_ED_Window = QLineEdit("5")
        self.lineEdit_ED_Window.setFixedWidth(45)

        g.addWidget(self.checkBox_ED,               0, 0)
        g.addWidget(self.checkBox_ED_Stats,          0, 1)
        g.addWidget(self.label_79,                   0, 2)
        g.addWidget(self.doubleSpinBox_ED_Threshold, 0, 3)
        g.addWidget(self.label_80,                   0, 4)
        g.addWidget(self.lineEdit_ED_Window,         0, 5)

        g.setColumnStretch(6, 1)
        return gb

    # ------------------------------------------------------------------
    # TAB 3 – Grid + Wavelets
    # ------------------------------------------------------------------
    def _tab_grid_wavelets(self):
        tab = QWidget()
        tab.setObjectName("Data_management")
        tab.setStyleSheet("#Data_management { " + TAB_BG + " }")

        outer = QVBoxLayout(tab)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        # Fixed bar: gridding and worms run as background tasks
        bar = QWidget()
        bar_layout = QHBoxLayout(bar)
        bar_layout.setContentsMargins(6, 4, 6, 0)
        self.pushButton_cancel_grid = QPushButton(_tr("Cancel running calculation"))
        self.pushButton_cancel_grid.setToolTip(
            _tr("Stop the running calculation. It stops at its next checkpoint "
                "and nothing from it is saved.")
        )
        self.pushButton_cancel_grid.setEnabled(False)  # only while one is running
        bar_layout.addWidget(self.pushButton_cancel_grid)
        bar_layout.addStretch()
        outer.addWidget(bar)

        sc = QWidget()
        sl = QVBoxLayout(sc)
        sl.setContentsMargins(6, 6, 6, 6)
        sl.setSpacing(6)
        sl.addWidget(self._group_import_points())
        sl.addWidget(self._group_gridding())
        sl.addWidget(self._group_worm_header())
        sl.addWidget(self._group_worms())
        # WTMM (Wavelet Transform Modulus Maxima) is retained in code but its
        # UI is hidden — build the widgets so SGTool.py references stay valid,
        # then keep the group off-screen.
        self._wtmm_group = self._group_wtmm()
        self._wtmm_group.setVisible(False)
        sl.addWidget(self._wtmm_group)
        sl.addStretch()

        outer.addWidget(self._scroll_wrap(sc), 1)
        outer.addWidget(self._grid_info_bar())
        return tab

    # ------------------------------------------------------------------
    # Group: Import point or line data
    # ------------------------------------------------------------------
    def _group_import_points(self):
        gb = QGroupBox(_tr("Import point or line data"))
        g = QGridLayout(gb)
        g.setSpacing(4)

        # Row 0: point/line data file
        self.label_53 = QLabel(_tr("Point/Line data"))
        self.label_53.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.lineEdit_loadPointsPath = QLineEdit()
        self.pushButton_selectPoints = QPushButton("...")
        self.pushButton_selectPoints.setFixedWidth(30)

        g.addWidget(self.label_53,              0, 0)
        g.addWidget(self.lineEdit_loadPointsPath, 0, 1)
        g.addWidget(self.pushButton_selectPoints, 0, 2)

        # Row 1: CRS
        self.mQgsProjectionSelectionWidget = QgsProjectionSelectionWidget()
        self.mQgsProjectionSelectionWidget.setEnabled(False)
        g.addWidget(self.mQgsProjectionSelectionWidget, 1, 0, 1, 3)

        # Row 2: X column
        self.label_44 = QLabel("long_x")
        self.label_44.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.comboBox_grid_x = QComboBox()
        self.comboBox_grid_x.setEnabled(False)
        g.addWidget(self.label_44,       2, 0)
        g.addWidget(self.comboBox_grid_x, 2, 1, 1, 2)

        # Row 3: Y column
        self.label_45 = QLabel("lat_y")
        self.label_45.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.comboBox_grid_y = QComboBox()
        self.comboBox_grid_y.setEnabled(False)

        self.checkBox_load_tie_lines = QCheckBox(_tr("Load Tie Lines from XYZ files"))
        self.pushButton_load_point_data = QPushButton(_tr("Import Data"))
        self.pushButton_load_point_data.setStyleSheet("font-weight: bold;")
        self.pushButton_load_point_data.setEnabled(False)

        g.addWidget(self.label_45,       3, 0)
        g.addWidget(self.comboBox_grid_y, 3, 1)
        g.addWidget(self.checkBox_load_tie_lines, 3, 2)

        g.addWidget(self.pushButton_load_point_data, 4, 1, 1, 2)

        g.setColumnStretch(1, 1)
        return gb

    # ------------------------------------------------------------------
    # Group: Gridding
    # ------------------------------------------------------------------
    def _group_gridding(self):
        gb = QGroupBox(_tr("Gridding"))
        g = QGridLayout(gb)
        g.setSpacing(4)

        # Select points layer + data field
        self.label_54 = QLabel(_tr("Select points layer"))
        self.label_54.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.mMapLayerComboBox_selectGrid_3 = QgsMapLayerComboBox()
        self.label_55 = QLabel(_tr("Select Data Field"))
        self.label_55.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.comboBox_select_grid_data_field = QComboBox()

        g.addWidget(self.label_54,                       0, 0)
        g.addWidget(self.mMapLayerComboBox_selectGrid_3, 0, 1, 1, 3)
        g.addWidget(self.label_55,                       1, 0)
        g.addWidget(self.comboBox_select_grid_data_field, 1, 1, 1, 3)

        # Cell size + # cells
        self.label_77 = QLabel(_tr("Cell Size\n(proj Units)"))
        self.label_77.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.doubleSpinBox_cellsize = QDoubleSpinBox()
        self.doubleSpinBox_cellsize.setDecimals(6)
        self.doubleSpinBox_cellsize.setMinimum(0.000001)
        self.doubleSpinBox_cellsize.setMaximum(9999999.99)
        self.doubleSpinBox_cellsize.setSingleStep(10.0)
        self.doubleSpinBox_cellsize.setValue(100.0)

        self.label_51 = QLabel(_tr("# Cells"))
        self.label_51.setStyleSheet("font-weight: bold;")
        self.label_49 = QLabel("X")
        self.label_49.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.nx_label = QLabel("0")
        self.label_48 = QLabel("Y")
        self.label_48.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.ny_label = QLabel("0")

        g.addWidget(self.label_77,            2, 0)
        g.addWidget(self.doubleSpinBox_cellsize, 2, 1)
        g.addWidget(self.label_51,            2, 2)
        g.addWidget(self.label_49,            2, 3)
        g.addWidget(self.nx_label,            2, 4)
        g.addWidget(self.label_48,            3, 3)
        g.addWidget(self.ny_label,            3, 4)

        # ── IDW gridding (its own row) ───────────────────────────────
        self.label_idw_hdr = QLabel(_tr("Inverse Distance Weighting"))
        self.label_idw_hdr.setStyleSheet("font-weight: bold;")
        g.addWidget(self.label_idw_hdr, 4, 0, 1, 4)
        self.pushButton_idw_2 = QPushButton(_tr("IDW Gridding"))
        self.pushButton_idw_2.setStyleSheet("font-weight: bold;")
        g.addWidget(self.pushButton_idw_2, 5, 1, 1, 3)

        # ── Multilevel B-Spline (MBA) block ──────────────────────────
        self.label_bspline_hdr = QLabel(_tr("Multilevel B-Spline (MBA)"))
        self.label_bspline_hdr.setStyleSheet("font-weight: bold;")
        g.addWidget(self.label_bspline_hdr, 6, 0, 1, 4)

        self.label_bspline_eps = QLabel(_tr("Threshold Error"))
        self.label_bspline_eps.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.lineEdit_bspline_epsilon = QLineEdit("0.0001")
        self.lineEdit_bspline_epsilon.setFixedWidth(70)

        self.label_bspline_lvl = QLabel(_tr("Max Levels"))
        self.label_bspline_lvl.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.spinBox_bspline_levels = QSpinBox()
        self.spinBox_bspline_levels.setMinimum(1)
        self.spinBox_bspline_levels.setMaximum(20)
        self.spinBox_bspline_levels.setValue(11)
        self.spinBox_bspline_levels.setFixedWidth(55)

        self.checkBox_bspline_ignore = QCheckBox(_tr("Ignore values <"))
        self.lineEdit_bspline_ignore = QLineEdit("-99999")
        self.lineEdit_bspline_ignore.setFixedWidth(70)

        g.addWidget(self.label_bspline_eps,        7, 0)
        g.addWidget(self.lineEdit_bspline_epsilon, 7, 1)
        g.addWidget(self.label_bspline_lvl,        7, 2)
        g.addWidget(self.spinBox_bspline_levels,   7, 3)
        g.addWidget(self.checkBox_bspline_ignore,  8, 0, 1, 2)
        g.addWidget(self.lineEdit_bspline_ignore,  8, 2)

        self.checkBox_bspline_auto = QCheckBox(_tr("Auto"))
        self.checkBox_bspline_auto.setToolTip(_tr(
            "Choose the cell size and Max Levels from the spacing of the points.\n"
            "Editing either value yourself switches Auto off."))
        g.addWidget(self.checkBox_bspline_auto, 6, 3, 1, 2)

        self.label_bspline_auto_info = QLabel("")
        self.label_bspline_auto_info.setWordWrap(True)
        self.label_bspline_auto_info.setStyleSheet("color: #555555;")
        self.label_bspline_auto_info.setVisible(False)
        g.addWidget(self.label_bspline_auto_info, 10, 0, 1, 5)

        self.pushButton_bspline_3 = QPushButton(_tr("BSpline Gridding"))
        self.pushButton_bspline_3.setStyleSheet("font-weight: bold;")
        g.addWidget(self.pushButton_bspline_3, 9, 1, 1, 3)

        g.setColumnStretch(1, 1)
        return gb

    # ------------------------------------------------------------------
    # Worm header (select grid)
    # ------------------------------------------------------------------
    def _group_worm_header(self):
        w = QWidget()
        h = QHBoxLayout(w)
        h.setContentsMargins(4, 0, 4, 0)
        self.label_52 = QLabel(_tr("Select Grid"))
        self.label_52.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.mMapLayerComboBox_selectGrid_worms = QgsMapLayerComboBox()
        h.addWidget(self.label_52)
        h.addWidget(self.mMapLayerComboBox_selectGrid_worms, 1)
        return w

    # ------------------------------------------------------------------
    # Group: BSDWorms
    # ------------------------------------------------------------------
    def _group_worms(self):
        gb = QGroupBox("BSDWorms")
        self.groupBox_8 = gb
        g = QGridLayout(gb)
        g.setSpacing(4)

        self.label_58 = QLabel(_tr("# Levels"))
        self.label_58.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.spinBox_levels = QSpinBox()
        self.spinBox_levels.setMinimum(1)
        self.spinBox_levels.setValue(10)
        self.checkBox_worms_shp = QCheckBox(_tr("Also save to shapefile"))

        self.label_60 = QLabel(_tr("Base Level"))
        self.label_60.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.doubleSpinBox_base = QDoubleSpinBox()
        self.doubleSpinBox_base.setDecimals(0)
        self.doubleSpinBox_base.setMinimum(0)
        self.doubleSpinBox_base.setMaximum(100000)
        self.doubleSpinBox_base.setSingleStep(100)
        self.doubleSpinBox_base.setValue(1000)

        self.label_62 = QLabel(_tr("Increment"))
        self.label_62.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.doubleSpinBox_inc = QDoubleSpinBox()
        self.doubleSpinBox_inc.setDecimals(0)
        self.doubleSpinBox_inc.setMinimum(0)
        self.doubleSpinBox_inc.setMaximum(100000)
        self.doubleSpinBox_inc.setSingleStep(100)
        self.doubleSpinBox_inc.setValue(1000)

        self.pushButton_worms = QPushButton(_tr("Calculate Worms"))
        self.pushButton_worms.setStyleSheet("font-weight: bold;")

        g.addWidget(self.label_58,         0, 0)
        g.addWidget(self.spinBox_levels,   0, 1)
        g.addWidget(self.checkBox_worms_shp, 0, 2, 1, 2)
        g.addWidget(self.label_60,         1, 0)
        g.addWidget(self.doubleSpinBox_base, 1, 1)
        g.addWidget(self.label_62,         2, 0)
        g.addWidget(self.doubleSpinBox_inc, 2, 1)
        g.addWidget(self.pushButton_worms, 2, 2, 1, 2)

        g.setColumnStretch(3, 1)
        return gb

    # ------------------------------------------------------------------
    # Group: Wavelet Transform Modulus Maxima
    # ------------------------------------------------------------------
    def _group_wtmm(self):
        gb = QGroupBox(_tr("Wavelet Transform Modulus Maxima"))
        g = QGridLayout(gb)
        g.setSpacing(4)

        # Select points/line layer
        self.label_57 = QLabel(_tr("Select points/line"))
        self.label_57.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.mMapLayerComboBox_selectVectors = QgsMapLayerComboBox()

        self.label_59 = QLabel(_tr("Select feature"))
        self.label_59.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.mFieldComboBox_feature = QComboBox()

        self.label_76 = QLabel(_tr("Select data field"))
        self.label_76.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.mFieldComboBox_data = QgsFieldComboBox()

        self.label_61 = QLabel(_tr("Spacing\n(Proj Units)"))
        self.label_61.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.doubleSpinBox_wtmm_spacing = QDoubleSpinBox()
        self.doubleSpinBox_wtmm_spacing.setDecimals(8)
        self.doubleSpinBox_wtmm_spacing.setMinimum(0)
        self.doubleSpinBox_wtmm_spacing.setMaximum(100000)
        self.doubleSpinBox_wtmm_spacing.setSingleStep(100)
        self.doubleSpinBox_wtmm_spacing.setValue(100)

        self.pushButton_wtmm = QPushButton(_tr("Calculate WTMM"))
        self.pushButton_wtmm.setStyleSheet("font-weight: bold;")

        g.addWidget(self.label_57,                    0, 0)
        g.addWidget(self.mMapLayerComboBox_selectVectors, 0, 1, 1, 3)
        g.addWidget(self.label_59,                    1, 0)
        g.addWidget(self.mFieldComboBox_feature,      1, 1, 1, 3)
        g.addWidget(self.label_76,                    2, 0)
        g.addWidget(self.mFieldComboBox_data,         2, 1, 1, 3)
        g.addWidget(self.label_61,                    3, 0)
        g.addWidget(self.doubleSpinBox_wtmm_spacing,  3, 1)
        g.addWidget(self.pushButton_wtmm,             3, 2, 1, 2)

        g.setColumnStretch(3, 1)
        return gb

    # ------------------------------------------------------------------
    # TAB 4 – Petrophysics
    # ------------------------------------------------------------------
    def _tab_petrophysics(self):
        tab = QWidget()
        tab.setObjectName("Petrophysics")
        tab.setStyleSheet("#Petrophysics { " + TAB_BG + " }")

        outer = QVBoxLayout(tab)
        outer.setContentsMargins(6, 6, 6, 6)
        outer.setSpacing(4)

        # ── Fixed top: Select Grid ───────────────────────────────────
        top = QWidget()
        tl = QHBoxLayout(top)
        tl.setContentsMargins(0, 0, 0, 0)
        lbl = QLabel(_tr("Select Grid"))
        lbl.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.mMapLayerComboBox_selectGrid_Petro = QgsMapLayerComboBox()
        self.mMapLayerComboBox_selectGrid_Petro.setToolTip(
            _tr("Grid the Apparent Petrophysics calculations are applied to")
        )
        tl.addWidget(lbl)
        tl.addWidget(self.mMapLayerComboBox_selectGrid_Petro, 1)
        outer.addWidget(top)

        # ── Fixed: Apply Processing + live preview controls ─────────
        self.pushButton_3_applyProcessing_Petro = QPushButton(_tr("Apply Processing"))
        self.pushButton_3_applyProcessing_Petro.setStyleSheet("font-weight: bold;")
        self.pushButton_3_applyProcessing_Petro.setToolTip(
            _tr("Run the ticked Apparent Petrophysics calculations")
        )
        outer.addWidget(
            self._apply_preview_bar(
                self.pushButton_3_applyProcessing_Petro, suffix="_petro"
            )
        )

        # ── Scrollable middle ────────────────────────────────────────
        sc = QWidget()
        sl = QVBoxLayout(sc)
        sl.setContentsMargins(2, 2, 2, 2)
        sl.setSpacing(6)
        sl.addWidget(self._group_apparent_petro())
        sl.addWidget(self._group_apparent_henkel())
        sl.addWidget(self._group_measured_petro())
        sl.addWidget(self._group_save_plots())
        sl.addStretch()
        outer.addWidget(self._scroll_wrap(sc), 1)
        outer.addWidget(self._grid_info_bar())
        return tab

    @staticmethod
    def _petro_row(layout, *items, indent=True):
        """A row of widgets in a line, indented under its checkbox."""
        box = QHBoxLayout()
        box.setContentsMargins(0, 0, 0, 0)
        box.setSpacing(4)
        if indent:
            box.addSpacing(22)
        for item in items:
            if isinstance(item, str):
                label = QLabel(item)
                label.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
                box.addWidget(label)
            else:
                box.addWidget(item)
        box.addStretch()
        layout.addLayout(box)

    @staticmethod
    def _plot_button(text):
        """A button that draws a plot as soon as it is pressed."""
        button = QPushButton(text)
        button.setStyleSheet("font-weight: bold;")
        return button

    @staticmethod
    def _short_edit(text, width=60, tip=None):
        edit = QLineEdit(text)
        edit.setFixedWidth(width)
        if tip:
            edit.setToolTip(tip)
        return edit

    # ------------------------------------------------------------------
    # Group: Apparent Petrophysics
    # ------------------------------------------------------------------
    def _group_apparent_petro(self):
        gb = QGroupBox(_tr("Apparent Petrophysics"))
        v = QVBoxLayout(gb)
        v.setSpacing(4)
        row = lambda *items, indent=True: self._petro_row(v, *items, indent=indent)
        cut_tip = _tr(
            "High-cut wavelength (m). Dividing out the layer response amplifies short\n"
            "wavelengths, so they are removed with a cosine roll-off: wavelengths longer\n"
            "than twice this value are kept, those at this value or shorter removed.\n"
            "Leave blank for 4 x the cell size; 0 for no high-cut."
        )

        # Apparent susceptibility
        self.radioButton_AppSus = QRadioButton(_tr("Apparent Susceptibility"))
        self.radioButton_AppSus.setChecked(True)
        self.radioButton_AppSus.setToolTip(
            _tr("Apparent susceptibility (SI) from a magnetic grid (nT) that has already been\n"
                "reduced to the pole (RTP), or to the equator (RTE: the result is then only\n"
                "approximate).\n"
                "The ground is modelled as a layer of vertical prisms magnetised by induction\n"
                "in the ambient field, and the layer's response is divided out in the\n"
                "Fourier domain. With a finite thickness the mean level cannot be\n"
                "recovered, so values are relative to the background susceptibility set.\n"
                "Output: <grid>_AppSus.")
        )
        row(self.radioButton_AppSus, indent=False)

        self.lineEdit_AppSus_field = self._short_edit(
            "50000", 60,
            _tr("Ambient field strength (nT).\n"
                "Filled in from the IGRF at the centre of the selected grid on the date\n"
                "shown whenever a grid is selected or IGRF is pressed (or from the grid's\n"
                "own metadata if it holds the field intensity). You can also type a value.")
        )
        self.pushButton_AppSus_IGRF = QPushButton("IGRF")
        self.pushButton_AppSus_IGRF.setToolTip(
            _tr("Set the field strength from the IGRF model at the centre of the selected\n"
                "grid on the survey date")
        )
        self.dateEdit_AppSus = QDateEdit(QDate.currentDate())
        self.dateEdit_AppSus.setMinimumDate(QDate(1900, 1, 1))
        self.dateEdit_AppSus.setMaximumDate(QDate(2030, 12, 31))
        self.dateEdit_AppSus.setCalendarPopup(True)
        self.dateEdit_AppSus.setToolTip(_tr("Survey date used for the IGRF field strength"))
        row(_tr("Field (nT)"), self.lineEdit_AppSus_field,
            self.pushButton_AppSus_IGRF, self.dateEdit_AppSus)

        self.lineEdit_AppSus_depth = self._short_edit(
            "0", 60,
            _tr("Depth from the sensor to the top of the magnetic layer (m).\n"
                "0 puts the top of the sources at the sensor. For airborne data set it to\n"
                "the terrain clearance (plus any cover) so the sources sit at the ground.")
        )
        self.checkBox_AppSus_thickness = QCheckBox(_tr("Thickness (m)"))
        self.checkBox_AppSus_thickness.setToolTip(
            _tr("Untick for a layer of infinite thickness (the usual choice); tick to give\n"
                "the layer a finite thickness")
        )
        self.lineEdit_AppSus_thickness = self._short_edit(
            "1000", 60, _tr("Layer thickness (m), used when Thickness is ticked")
        )
        self.lineEdit_AppSus_thickness.setEnabled(False)
        self.checkBox_AppSus_thickness.toggled.connect(
            self.lineEdit_AppSus_thickness.setEnabled
        )
        row(_tr("Depth to top (m)"), self.lineEdit_AppSus_depth,
            self.checkBox_AppSus_thickness, self.lineEdit_AppSus_thickness)

        self.lineEdit_AppSus_background = self._short_edit(
            "0.001", 50,
            _tr("Background susceptibility (SI) added to the result.\n"
                "A magnetic anomaly grid (IGRF removed) has no absolute level, so the\n"
                "calculated susceptibility is relative to an unknown background and is\n"
                "partly negative. A typical background (e.g. 0.001 SI) keeps most values\n"
                "positive; Henkel plots and histograms leave out values <= 0, which a log axis can't show.")
        )
        self.lineEdit_AppSus_cut = self._short_edit("", 70, cut_tip)
        self.lineEdit_AppSus_cut.setPlaceholderText(_tr("4 x cell"))
        row(_tr("Background (SI)"), self.lineEdit_AppSus_background,
            _tr("High-cut (m)"), self.lineEdit_AppSus_cut)

        # Apparent density
        self.radioButton_AppDens = QRadioButton(_tr("Apparent Density"))
        self.radioButton_AppDens.setToolTip(
            _tr("Apparent density (g/cc) from a Bouguer gravity grid.\n"
                "The ground is modelled as a layer of vertical prisms of finite thickness,\n"
                "and the layer's response is divided out in the Fourier domain. A Bouguer\n"
                "grid's absolute level depends on its datum and the regional field, so it\n"
                "is not used: the result's mean is the background density.\n"
                "Output: <grid>_AppDens.")
        )
        self.comboBox_AppDens_units = QComboBox()
        self.comboBox_AppDens_units.addItems(["mGal", "µm/s² (gu)"])
        self.comboBox_AppDens_units.setToolTip(
            _tr("Units of the gravity grid (1 mGal = 10 µm/s² = 10 gu)")
        )
        row(self.radioButton_AppDens, _tr("Grid units"), self.comboBox_AppDens_units,
            indent=False)

        self.lineEdit_AppDens_depth = self._short_edit(
            "0", 60,
            _tr("Depth from the sensor to the top of the layer (m).\n"
                "0 for ground gravity; the terrain clearance for airborne gravity.")
        )
        self.lineEdit_AppDens_thickness = self._short_edit(
            "5000", 60,
            _tr("Layer thickness (m); must be greater than 0. Thin layers exaggerate the\n"
                "density contrasts needed to explain long-wavelength anomalies.")
        )
        row(_tr("Depth to top (m)"), self.lineEdit_AppDens_depth,
            _tr("Thickness (m)"), self.lineEdit_AppDens_thickness)

        self.lineEdit_AppDens_background = self._short_edit(
            "2.67", 50, _tr("Background density (g/cc) added to the density contrast")
        )
        self.lineEdit_AppDens_cut = self._short_edit("", 70, cut_tip)
        self.lineEdit_AppDens_cut.setPlaceholderText(_tr("4 x cell"))
        row(_tr("Background (g/cc)"), self.lineEdit_AppDens_background,
            _tr("High-cut (m)"), self.lineEdit_AppDens_cut)

        # one calculation at a time; Apply Processing on this tab runs it
        group = QButtonGroup(self)
        group.addButton(self.radioButton_AppSus)
        group.addButton(self.radioButton_AppDens)
        self.buttonGroup_AppPetro = group
        return gb

    # ------------------------------------------------------------------
    # Group: Apparent Plots (Henkel plot, or histogram of one grid)
    # ------------------------------------------------------------------
    def _group_apparent_henkel(self):
        gb = QGroupBox(_tr("Apparent Henkel Plot / Histogram"))
        g = QGridLayout(gb)
        g.setSpacing(4)

        lbl_sus = QLabel(_tr("App. Susceptibility"))
        lbl_sus.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.mMapLayerComboBox_Henkel_sus = QgsMapLayerComboBox()
        self.mMapLayerComboBox_Henkel_sus.setToolTip(
            _tr("Apparent susceptibility grid (SI). Leave blank, with a density grid\n"
                "chosen, for density histograms.")
        )
        lbl_dens = QLabel(_tr("App. Density"))
        lbl_dens.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.mMapLayerComboBox_Henkel_dens = QgsMapLayerComboBox()
        self.mMapLayerComboBox_Henkel_dens.setToolTip(
            _tr("Apparent density grid (g/cc). Leave blank, with a susceptibility grid\n"
                "chosen, for susceptibility histograms.")
        )

        self.pushButton_Henkel = self._plot_button(_tr("Plot Apparent"))
        self.pushButton_Henkel.setToolTip(
            _tr("Heat map of apparent density (linear x axis) against apparent susceptibility\n"
                "(log y axis), counting grid cells per bin on a log colour scale.\n"
                "Each cell of the coarser grid (usually gravity) is paired with the nearest\n"
                "cell of the finer grid, reprojecting if the CRSs differ; cells where either\n"
                "grid has no data are left out, as are susceptibilities <= 0, which a log\n"
                "axis can't show (their number is given in the title). Axes span the\n"
                "0.1-99.9 percentiles, to keep FFT edge effects from squashing the plot.\n"
                "With only one grid chosen: a histogram of it instead (linear bins for\n"
                "density, log bins for susceptibility).")
        )
        # Per-category plots of the same two grids, clipped to polygons
        lbl_poly = QLabel(_tr("Polygons"))
        lbl_poly.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.mMapLayerComboBox_Petro_polygons = QgsMapLayerComboBox()
        self.mMapLayerComboBox_Petro_polygons.setToolTip(
            _tr("Polygon layer, e.g. a geology map. If some polygons are selected, only\n"
                "those are used.")
        )
        lbl_cat = QLabel(_tr("Category field"))
        lbl_cat.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.mFieldComboBox_Petro_polyCat = QgsFieldComboBox()
        self.mFieldComboBox_Petro_polyCat.setToolTip(
            _tr("Field holding the category (e.g. lithology) of each polygon")
        )
        self.pushButton_Henkel_categories = self._plot_button(_tr("Plot per Category"))
        self.pushButton_Henkel_categories.setToolTip(
            _tr("A grid of Henkel heat maps (or histograms, with only one grid chosen),\n"
                "one per category, using the grids above\n"
                "clipped to the polygons of that category. Categories\n"
                "with no cells where both grids have data are left out. All panels\n"
                "share their axes and colour scale; each shows the fraction of its own\n"
                "category's cells per bin, so small and large categories compare fairly.\n"
                "With only one grid chosen: a histogram of it per category instead.")
        )
        divider = QFrame()
        divider.setFrameShape(QFrame.Shape.HLine)
        divider.setFrameShadow(QFrame.Shadow.Sunken)

        g.addWidget(lbl_sus,                               0, 0)
        g.addWidget(self.mMapLayerComboBox_Henkel_sus,     0, 1)
        g.addWidget(lbl_dens,                              1, 0)
        g.addWidget(self.mMapLayerComboBox_Henkel_dens,    1, 1)
        g.addWidget(self.pushButton_Henkel,                2, 0, 1, 2, Qt.AlignmentFlag.AlignLeft)
        g.addWidget(divider,                               3, 0, 1, 2)
        g.addWidget(lbl_poly,                              4, 0)
        g.addWidget(self.mMapLayerComboBox_Petro_polygons, 4, 1)
        g.addWidget(lbl_cat,                               5, 0)
        g.addWidget(self.mFieldComboBox_Petro_polyCat,     5, 1)
        self.checkBox_Henkel_position = QCheckBox(_tr("Colour by position"))
        self.checkBox_Henkel_position.setToolTip(
            _tr("For Henkel plots per category: colour each bin by where it is on the plot,\n"
                "red to green along the density axis and red to blue up the susceptibility\n"
                "axis (low/low red, dense green, magnetic blue, both cyan), with the\n"
                "fraction of the category's cells per bin shown by opacity instead of a\n"
                "colour scale.")
        )
        category_row = QHBoxLayout()
        category_row.setSpacing(10)
        category_row.addWidget(self.pushButton_Henkel_categories)
        category_row.addWidget(self.checkBox_Henkel_position)
        category_row.addStretch()
        g.addLayout(category_row,                          6, 0, 1, 2)
        g.setColumnStretch(1, 1)
        return gb

    # ------------------------------------------------------------------
    # Group: Measured Petrophysics
    # ------------------------------------------------------------------
    def _group_measured_petro(self):
        gb = QGroupBox(_tr("Measured Petrophysics"))
        g = QGridLayout(gb)
        g.setSpacing(4)
        field_tip = _tr(
            "Field plotted on the {axis} axis, and named on it: numeric, or text holding\n"
            "numbers (values that aren't numbers are skipped). Leave blank, with the\n"
            "{other} field chosen, for histograms of that field."
        )
        scale_tip = _tr(
            "Multiplier applied to the field's values, e.g. 0.001 for density in kg/m3\n"
            "to g/cc, or 0.00001 for susceptibility in x10^-5 SI to SI"
        )
        log_tip = _tr(
            "Log scale for the {axis} axis (and log bins in histograms). Values <= 0\n"
            "can't be shown on a log axis, so they are left out (their number is given\n"
            "in the title)."
        )

        def label(text):
            lbl = QLabel(text)
            lbl.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            return lbl

        self.mMapLayerComboBox_Petro_points = QgsMapLayerComboBox()
        self.mMapLayerComboBox_Petro_points.setToolTip(
            _tr("Point layer of samples with measured properties, e.g. density and\n"
                "susceptibility. If some points are selected, only those are used.")
        )
        self.mFieldComboBox_Petro_x = QgsFieldComboBox()
        self.mFieldComboBox_Petro_x.setToolTip(field_tip.format(axis="x", other="Y"))
        self.lineEdit_Petro_x_scale = self._short_edit("1", 55, scale_tip)
        self.checkBox_Petro_x_log = QCheckBox(_tr("Log"))
        self.checkBox_Petro_x_log.setToolTip(log_tip.format(axis="x"))
        self.mFieldComboBox_Petro_y = QgsFieldComboBox()
        self.mFieldComboBox_Petro_y.setToolTip(field_tip.format(axis="y", other="X"))
        self.lineEdit_Petro_y_scale = self._short_edit("1", 55, scale_tip)
        self.checkBox_Petro_y_log = QCheckBox(_tr("Log"))
        self.checkBox_Petro_y_log.setChecked(True)  # Henkel plot: log susceptibility
        self.checkBox_Petro_y_log.setToolTip(log_tip.format(axis="y"))
        self.mFieldComboBox_Petro_cat = QgsFieldComboBox()
        self.mFieldComboBox_Petro_cat.setToolTip(
            _tr("Optional: field holding a category (e.g. lithology). With Points, each\n"
                "category gets its own symbol and a legend; with only one of the X and Y\n"
                "fields, a histogram per category.")
        )

        self.pushButton_Henkel_points = self._plot_button(_tr("Plot Samples"))
        self.pushButton_Henkel_points.setToolTip(
            _tr("The Y field against the X field for every sample (a Henkel plot with\n"
                "density as X and susceptibility as Y, log). The axes are named after\n"
                "the fields.\n"
                "With only one field chosen: a histogram of it instead, one per category\n"
                "if a category field is chosen.")
        )
        self.radioButton_Petro_points = QRadioButton(_tr("Points"))
        self.radioButton_Petro_points.setChecked(True)
        self.radioButton_Petro_points.setToolTip(
            _tr("Plot every sample, with a colour and marker per category if one is\n"
                "chosen: 8 colours x 8 markers give up to 64 categories, assigned to the\n"
                "categories with the most samples (beyond 64, the rest are grouped as\n"
                "Other).")
        )
        self.radioButton_Petro_heat = QRadioButton(_tr("Heat map"))
        self.radioButton_Petro_heat.setToolTip(
            _tr("Samples per bin on a log colour scale, for large datasets (ignores category)")
        )
        group = QButtonGroup(self)
        group.addButton(self.radioButton_Petro_points)
        group.addButton(self.radioButton_Petro_heat)
        self.buttonGroup_Petro_style = group

        g.addWidget(label(_tr("Points")),                0, 0)
        g.addWidget(self.mMapLayerComboBox_Petro_points, 0, 1, 1, 4)
        g.addWidget(label(_tr("X field")),               1, 0)
        g.addWidget(self.mFieldComboBox_Petro_x,         1, 1)
        g.addWidget(label(_tr("x")),                     1, 2)
        g.addWidget(self.lineEdit_Petro_x_scale,         1, 3)
        g.addWidget(self.checkBox_Petro_x_log,           1, 4)
        g.addWidget(label(_tr("Y field")),               2, 0)
        g.addWidget(self.mFieldComboBox_Petro_y,         2, 1)
        g.addWidget(label(_tr("x")),                     2, 2)
        g.addWidget(self.lineEdit_Petro_y_scale,         2, 3)
        g.addWidget(self.checkBox_Petro_y_log,           2, 4)
        g.addWidget(label(_tr("Category field")),        3, 0)
        g.addWidget(self.mFieldComboBox_Petro_cat,       3, 1)

        style_row = QHBoxLayout()
        style_row.setSpacing(6)
        style_row.addWidget(self.pushButton_Henkel_points)
        style_row.addSpacing(10)
        style_row.addWidget(self.radioButton_Petro_points)
        style_row.addWidget(self.radioButton_Petro_heat)
        style_row.addStretch()
        g.addLayout(style_row, 4, 0, 1, 5)
        g.setColumnStretch(1, 1)
        return gb

    # ------------------------------------------------------------------
    # Group: Save Plots
    # ------------------------------------------------------------------
    def _group_save_plots(self):
        gb = QGroupBox(_tr("Save Plots"))
        h = QHBoxLayout(gb)
        h.setSpacing(4)
        self.checkBox_plot_png = QCheckBox("PNG")
        self.checkBox_plot_pdf = QCheckBox("PDF")
        for box in (self.checkBox_plot_png, self.checkBox_plot_pdf):
            box.setToolTip(_tr("Also save each plot in this format"))
        lbl = QLabel(_tr("Folder"))
        lbl.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.lineEdit_plot_folder = QLineEdit()
        self.lineEdit_plot_folder.setPlaceholderText(_tr("beside the input layer"))
        self.lineEdit_plot_folder.setToolTip(
            _tr("Folder for saved plots. Leave blank to save them beside the input layer.\n"
                "Plots can also be saved from the plot window's toolbar.")
        )
        self.pushButton_plot_folder = QPushButton("...")
        self.pushButton_plot_folder.setFixedWidth(30)
        h.addWidget(self.checkBox_plot_png)
        h.addWidget(self.checkBox_plot_pdf)
        h.addWidget(lbl)
        h.addWidget(self.lineEdit_plot_folder, 1)
        h.addWidget(self.pushButton_plot_folder)
        return gb

    # ------------------------------------------------------------------
    # TAB 5 – Utils
    # ------------------------------------------------------------------
    def _tab_utils(self):
        tab = QWidget()
        tab.setObjectName("Utils")
        tab.setStyleSheet("#Utils { " + TAB_BG + " }")

        outer = QVBoxLayout(tab)
        outer.setContentsMargins(6, 6, 6, 6)
        outer.setSpacing(4)

        # ── Fixed top: Select Grid ───────────────────────────────────
        top_row = QWidget()
        tl = QHBoxLayout(top_row)
        tl.setContentsMargins(0, 0, 0, 0)
        self.label_56 = QLabel(_tr("Select Grid"))
        self.label_56.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.mMapLayerComboBox_selectGrid_Conv_2 = QgsMapLayerComboBox()
        tl.addWidget(self.label_56)
        tl.addWidget(self.mMapLayerComboBox_selectGrid_Conv_2, 1)
        outer.addWidget(top_row)

        # ── Fixed: Apply Processing + live preview controls ─────────
        self.pushButton_3_applyProcessing_Conv_3 = QPushButton(_tr("Apply Processing"))
        self.pushButton_3_applyProcessing_Conv_3.setStyleSheet("font-weight: bold;")
        outer.addWidget(
            self._apply_preview_bar(
                self.pushButton_3_applyProcessing_Conv_3, suffix="_utils"
            )
        )

        # ── Scrollable middle ────────────────────────────────────────
        sc = QWidget()
        sl = QVBoxLayout(sc)
        sl.setContentsMargins(2, 2, 2, 2)
        sl.setSpacing(6)

        sl.addWidget(self._group_metadata())
        sl.addWidget(self._group_nan_threshold())
        sl.addWidget(self._group_clip_polygon())
        clip_row = QHBoxLayout()
        clip_row.setSpacing(6)
        clip_row.addWidget(self._group_clip_to_window(), 1)
        clip_row.addWidget(self._group_clip_all_layers(), 1)
        sl.addLayout(clip_row)
        sl.addWidget(self._group_normalise())
        sl.addWidget(self._group_lut_convert())
        sl.addStretch()

        outer.addWidget(self._scroll_wrap(sc), 1)
        outer.addWidget(self._grid_info_bar())
        return tab

    # ------------------------------------------------------------------
    # Group: SGTool provenance metadata reader
    # ------------------------------------------------------------------
    def _group_metadata(self):
        gb = QGroupBox(_tr("SGTool Metadata"))
        g = QGridLayout(gb)
        g.setSpacing(4)

        self.pushButton_read_metadata = QPushButton(
            _tr("Read metadata of selected grid")
        )
        self.pushButton_read_metadata.setToolTip(
            _tr(
                "Show how the selected grid was made: when, from which source "
                "file(s), and with which operation and parameters, back through "
                "every earlier SGTool step.\n"
                "This provenance is stored inside GeoTIFFs saved by SGTool."
            )
        )
        self.pushButton_save_metadata_xml = QPushButton(_tr("Save as XML"))
        self.pushButton_save_metadata_xml.setToolTip(
            _tr(
                "Save the selected grid's SGTool metadata as an XML file "
                "called <grid file>.sgt.xml in the same folder as the grid "
                "(replacing any existing file of that name)"
            )
        )
        self.pushButton_replay_history = QPushButton(
            _tr("Apply the same steps to another grid...")
        )
        self.pushButton_replay_history.setToolTip(
            _tr(
                "Repeat the processing steps in the selected grid's history on "
                "another grid, each on the result of the one before.\n"
                "Magnetic reductions (RTP, RTE, differential RTP) get a new "
                "inclination and declination from the IGRF model for the new grid."
            )
        )
        g.addWidget(self.pushButton_read_metadata, 0, 0)
        g.addWidget(self.pushButton_save_metadata_xml, 0, 1)
        g.addWidget(self.pushButton_replay_history, 0, 2)  # right of Save as XML
        g.setColumnStretch(3, 1)
        return gb

    # ------------------------------------------------------------------
    # Group: Threshold to NaN
    # ------------------------------------------------------------------
    def _group_nan_threshold(self):
        gb = QGroupBox(_tr("Threshold to NaN"))
        g = QGridLayout(gb)
        g.setSpacing(4)

        self.checkBox_NaN = QCheckBox()
        self.radioButton_NaN_Above = QRadioButton(_tr("Above"))
        self.radioButton_NaN_Above.setChecked(True)
        self.radioButton_NaN_Below = QRadioButton(_tr("Below"))
        self.radioButton_NaN_Both = QRadioButton(_tr("Between"))

        self.label_70 = QLabel(_tr("Above"))
        self.label_70.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.doubleSpinBox_NaN_Above = QDoubleSpinBox()
        self.doubleSpinBox_NaN_Above.setDecimals(8)
        self.doubleSpinBox_NaN_Above.setMinimum(-1e30)
        self.doubleSpinBox_NaN_Above.setMaximum(1e30)
        self.doubleSpinBox_NaN_Above.setValue(1.0)

        self.label_71 = QLabel(_tr("Below"))
        self.label_71.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.doubleSpinBox_NaN_Below = QDoubleSpinBox()
        self.doubleSpinBox_NaN_Below.setDecimals(8)
        self.doubleSpinBox_NaN_Below.setMinimum(-1e30)
        self.doubleSpinBox_NaN_Below.setMaximum(1e30)
        self.doubleSpinBox_NaN_Below.setValue(-1.0)

        g.addWidget(self.checkBox_NaN,         0, 0)
        g.addWidget(self.radioButton_NaN_Above, 0, 1)
        g.addWidget(self.radioButton_NaN_Below, 0, 2)
        g.addWidget(self.radioButton_NaN_Both,  0, 3)
        g.addWidget(self.label_70,              1, 0)
        g.addWidget(self.doubleSpinBox_NaN_Above, 1, 1, 1, 2)
        g.addWidget(self.label_71,              2, 0)
        g.addWidget(self.doubleSpinBox_NaN_Below, 2, 1, 1, 2)

        g.setColumnStretch(4, 1)
        return gb

    # ------------------------------------------------------------------
    # Group: Create Clipping Polygon
    # ------------------------------------------------------------------
    def _group_clip_polygon(self):
        gb = QGroupBox(_tr("Create Clipping Polygon"))
        g = QGridLayout(gb)
        g.setSpacing(4)
        self.checkBox_polygons = QCheckBox(_tr("Create polygon based on grid"))
        g.addWidget(self.checkBox_polygons, 0, 0)
        g.setColumnStretch(1, 1)
        return gb

    # ------------------------------------------------------------------
    # Group: Clip to Window
    # ------------------------------------------------------------------
    def _group_clip_to_window(self):
        gb = QGroupBox(_tr("Clip to Window"))
        g = QGridLayout(gb)
        g.setSpacing(4)
        self.pushButton_clip_to_window = QPushButton(_tr("Clip to Window"))
        self.pushButton_clip_to_window.setStyleSheet("font-weight: bold;")
        self.pushButton_clip_to_window.setToolTip(
            _tr("Clip the selected grid to the area currently shown in the map window.\n"
                "The result is saved beside the original as <name>_clip.tif (or _clip_01, "
                "_clip_02 ... if that name is already used) and loaded into the project.")
        )
        g.addWidget(self.pushButton_clip_to_window, 0, 0)
        g.setRowStretch(1, 1)
        return gb

    # ------------------------------------------------------------------
    # Group: Clip All Layers
    # ------------------------------------------------------------------
    def _group_clip_all_layers(self):
        gb = QGroupBox(_tr("Clip All Layers"))
        g = QGridLayout(gb)
        g.setSpacing(4)
        self.label_clip_polygon = QLabel(_tr("Clip polygon"))
        self.label_clip_polygon.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.mMapLayerComboBox_clip_polygon = QgsMapLayerComboBox()
        self.mMapLayerComboBox_clip_polygon.setToolTip(
            _tr("Polygon layer to clip every layer of the project to. If some of its "
                "polygons are selected, only those are used.")
        )
        self.pushButton_clip_all_layers = QPushButton(_tr("Clip all Layers"))
        self.pushButton_clip_all_layers.setStyleSheet("font-weight: bold;")
        self.pushButton_clip_all_layers.setToolTip(
            _tr("Choose a folder, then clip every layer of the project to the polygon.\n"
                "Each clipped layer is saved there as <name>_clip (GeoPackage or GeoTIFF) "
                "with a .qml file holding its style, and a copy of the project "
                "(.qgz) with the same groups and order, using the clipped layers, is saved "
                "beside them.")
        )
        g.addWidget(self.label_clip_polygon, 0, 0)
        g.addWidget(self.mMapLayerComboBox_clip_polygon, 0, 1)
        g.addWidget(self.pushButton_clip_all_layers, 1, 1)
        g.setColumnStretch(1, 1)
        return gb

    # ------------------------------------------------------------------
    # Group: Normalise Grids
    # ------------------------------------------------------------------
    def _group_normalise(self):
        gb = QGroupBox(_tr("Normalise Grids"))
        self.groupBox_10 = gb
        g = QGridLayout(gb)
        g.setSpacing(4)

        self.label_66 = QLabel(_tr("Input Directory"))
        self.label_66.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.lineEdit_loadPointsPath_normalise_in = QLineEdit()
        self.pushButton_select_normalise_in = QPushButton("...")
        self.pushButton_select_normalise_in.setFixedWidth(30)
        self.radioButton_normalise_1st = QRadioButton(_tr("1st order fit"))
        self.radioButton_normalise_1st.setChecked(True)

        self.label_72 = QLabel(_tr("Output Directory"))
        self.label_72.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.lineEdit_loadPointsPath_normalise_out = QLineEdit()
        self.pushButton_select_normalise_out = QPushButton("...")
        self.pushButton_select_normalise_out.setFixedWidth(30)
        self.radioButton_normalise_2nd = QRadioButton(_tr("2nd order fit"))

        self.pushButton_normalise = QPushButton(_tr("Normalise Grids"))
        self.pushButton_normalise.setStyleSheet("font-weight: bold;")

        g.addWidget(self.label_66,                           0, 0)
        g.addWidget(self.lineEdit_loadPointsPath_normalise_in, 0, 1)
        g.addWidget(self.pushButton_select_normalise_in,     0, 2)
        g.addWidget(self.radioButton_normalise_1st,          0, 3)

        g.addWidget(self.label_72,                            1, 0)
        g.addWidget(self.lineEdit_loadPointsPath_normalise_out, 1, 1)
        g.addWidget(self.pushButton_select_normalise_out,    1, 2)
        g.addWidget(self.radioButton_normalise_2nd,          1, 3)

        g.addWidget(self.pushButton_normalise,               2, 1, 1, 2)

        g.setColumnStretch(1, 1)
        return gb

    # ------------------------------------------------------------------
    # Group: Convert LUT to Grayscale
    # ------------------------------------------------------------------
    def _group_lut_convert(self):
        gb = QGroupBox(_tr("Convert LUT to Grayscale"))
        gb.setObjectName("groupBox_7")
        self.groupBox_7 = gb
        g = QGridLayout(gb)
        g.setSpacing(4)

        # RGB grid selector
        self.label_42 = QLabel(_tr("Select RGB Grid"))
        self.label_42.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.lineEdit_2_loadGridPath_2 = QLineEdit()
        self.pushButton_2_selectGrid_RGB = QPushButton("...")
        self.pushButton_2_selectGrid_RGB.setFixedWidth(30)

        g.addWidget(self.label_42,                0, 0)
        g.addWidget(self.lineEdit_2_loadGridPath_2, 0, 1)
        g.addWidget(self.pushButton_2_selectGrid_RGB, 0, 2)

        # CSS colour list
        self.label_43 = QLabel(_tr("CSS Colour List\nor\nRGB triplets"))
        self.label_43.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.textEdit_2_colour_list = QTextEdit()
        self.textEdit_2_colour_list.setPlaceholderText(
            _tr("comma separated list of CSS colour names or comma separated RGB triplets"))
        self.textEdit_2_colour_list.setFixedHeight(70)

        self.pushButton_pick_rgb = QPushButton("Pick RGB")
        self.pushButton_pick_rgb.setToolTip(
            _tr("Pick RGB from map: click button then click a point on the loaded RGB raster\nClick points equally spaced from low to high values"))

        g.addWidget(self.label_43,               1, 0)
        g.addWidget(self.textEdit_2_colour_list, 1, 1)
        g.addWidget(self.pushButton_pick_rgb,    1, 2)

        # Example line
        self.lineEdit = QLineEdit(
            "e.g.  teal, lemonchiffon,red  or 45,133,186,251,248,183,215,26,29")
        g.addWidget(self.lineEdit,               2, 1, 1, 2)

        # CSS colours button
        self.pushButton_CSSS_Colours = QPushButton(
            _tr("Full list of CSS Colour Names"))
        g.addWidget(self.pushButton_CSSS_Colours, 3, 1, 1, 2)

        # Min / Max
        self.label_46 = QLabel(_tr("Min"))
        self.label_46.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.mQgsDoubleSpinBox_LUT_min = QgsDoubleSpinBox()
        self.mQgsDoubleSpinBox_LUT_min.setMinimum(-1e20)
        self.mQgsDoubleSpinBox_LUT_min.setMaximum(1e20)
        self.mQgsDoubleSpinBox_LUT_min.setSingleStep(100)

        self.label_50 = QLabel(_tr("Max"))
        self.label_50.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.mQgsDoubleSpinBox_LUT_max = QgsDoubleSpinBox()
        self.mQgsDoubleSpinBox_LUT_max.setMinimum(-1e20)
        self.mQgsDoubleSpinBox_LUT_max.setMaximum(1e20)
        self.mQgsDoubleSpinBox_LUT_max.setSingleStep(100)
        self.mQgsDoubleSpinBox_LUT_max.setValue(1000)

        g.addWidget(self.label_46,                   4, 0)
        g.addWidget(self.mQgsDoubleSpinBox_LUT_min,  4, 1)
        g.addWidget(self.label_50,                   5, 0)
        g.addWidget(self.mQgsDoubleSpinBox_LUT_max,  5, 1)

        # Convert button
        self.pushButton_3_applyProcessing_Conv_2 = QPushButton(_tr("Convert Grid"))
        self.pushButton_3_applyProcessing_Conv_2.setStyleSheet("font-weight: bold;")
        g.addWidget(self.pushButton_3_applyProcessing_Conv_2, 6, 1, 1, 2)

        g.setColumnStretch(1, 1)
        return gb

    # ------------------------------------------------------------------
    # TAB 6 – Help
    # ------------------------------------------------------------------
    def _tab_help(self):
        tab = QWidget()
        tab.setObjectName("Help")
        tab.setStyleSheet("#Help { " + TAB_BG + " }")

        outer = QVBoxLayout(tab)
        outer.setContentsMargins(6, 6, 6, 6)

        self.textBrowser = QTextBrowser()
        self.textBrowser.setOpenExternalLinks(True)
        self.textBrowser.setHtml(self._help_html())
        outer.addWidget(self.textBrowser)
        return tab

    @staticmethod
    def _help_html():
        return """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0//EN">
<html><body style="font-family:'MS Shell Dlg 2'; font-size:8pt;">
<p><b>Simple Potential Field Calcs to assist WAXI/Agate Structural Geophysics Course</b></p>
<p><a href="https://waxi4.org">https://waxi4.org</a> &nbsp;
   <a href="https://agate-project.org">https://agate-project.org</a></p>
<p><b>Example magnetic grid from OMRG Mauritania:</b></p>
<p><a href="http://tectonique.net/sgtools_data/ogrm_usgs_mag_tmi.tif">grid</a></p>
<p><b>Noddy Imports:</b></p>
<p>Existing Noddy mag and grav files can be found at the Atlas of Structural Geophysics:
<a href="https://tectonique.net/asg/">https://tectonique.net/asg/</a></p>
<p>New Noddy models can be calculated using the Windows version at
<a href="https://tectonique.net/noddy/OpenNoddy_installer.exe">https://tectonique.net/noddy/OpenNoddy_installer.exe</a>
or a python wrapper at <a href="https://github.com/cgre-aachen/pynoddy">https://github.com/cgre-aachen/pynoddy</a></p>
<p><b>Code Repository</b></p>
<p>This plugin is available directly within QGIS using the Plugin Manager,</p>
<p>however the latest version will always be at this GitHub site:</p>
<p><a href="https://github.com/swaxi/SGTool"><b>https://github.com/swaxi/SGTool</b></a></p>
<p><b>Documentation</b></p>
<p><a href="https://tectonique.net/sgtools_data/Structural%20Geophysics%20Tools.pdf"><b>Help File</b></a></p>
<p>Help File reflects latest changes on GitHub. If a feature is not available in the version
you are using, go to Code Repository for latest version</p>
<p><b>Code development</b></p>
<p>- Calcs Claude, ChatGPT and Mark Jessell</p>
<p>- Plugin construction - Mark Jessell using QGIS Plugin Builder Plugin
<a href="https://g-sherman.github.io/Qgis-Plugin-Builder/">https://g-sherman.github.io/Qgis-Plugin-Builder/</a></p>
<p>- IGRF calculation - using pyIGRF
<a href="https://github.com/ciaranbe/pyIGRF">https://github.com/ciaranbe/pyIGRF</a></p>
<p>- GRD Loader &amp; Radially averaged power spectrum Fatiando a Terra
<a href="https://www.fatiando.org/">https://www.fatiando.org/</a></p>
<p>- Example geophysics data above courtesy of Mauritania Govt.
<a href="https://anarpam.mr/en/">https://anarpam.mr/en/</a></p>
<p>- Worming of grids uses Frank Horowitz's bsdwormer
<a href="https://bitbucket.org/fghorow/bsdwormer/">https://bitbucket.org/fghorow/bsdwormer/</a></p>
<p>- Multilevel B-Spline (MBA) Gridding — native Python translation of SAGA's grid_spline algorithm (no external plugin required)</p>
<p>- Differential (Variable) RTP — G.R.J. Cooper & D.R. Cowan (2005) Taylor-series method</p>
<p>- MRVBF/MRRTF (Multiresolution Valley Bottom/Ridge Top Flatness) — J.C. Gallant & T.I. Dowling (2003)</p>
<p>- Euler Deconvolution uses Felipe F. Melo and Valéria C.F. Barbosa's Reliable Euler method https://github.com/ffigura/Euler-deconvolution-python</p>
</body></html>"""
