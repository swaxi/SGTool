# -*- coding: utf-8 -*-
from qgis.PyQt import QtWidgets
from qgis.PyQt.QtCore import pyqtSignal, Qt, QDate, QCoreApplication
from qgis.PyQt.QtWidgets import (
    QDockWidget, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QScrollArea, QTabWidget, QPushButton, QLabel, QLineEdit,
    QCheckBox, QRadioButton, QComboBox, QGroupBox, QDateEdit,
    QTextEdit, QTextBrowser, QSpinBox, QDoubleSpinBox, QFrame,
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
        try:
            self.tabWidget.setTabShape(QTabWidget.TabShape.Triangular)
        except AttributeError:
            self.tabWidget.setTabShape(QTabWidget.Triangular)
        self.tabWidget.setStyleSheet(TAB_STYLE)
        self.tabWidget.setUsesScrollButtons(False)

        self.tabWidget.addTab(self._tab_fft_filters(),   _tr("FFT Filters"))
        self.tabWidget.addTab(self._tab_conv_stats(),    _tr("Conv + Stats"))
        self.tabWidget.addTab(self._tab_grid_wavelets(), _tr("Grid + Wavelets"))
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

        self.setWidget(root)

    # ------------------------------------------------------------------
    # Helper: wrap a widget list in a QScrollArea (vertical only)
    # ------------------------------------------------------------------
    @staticmethod
    def _scroll_wrap(content_widget):
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        scroll.setFrameShape(QFrame.NoFrame)
        scroll.setWidget(content_widget)
        return scroll

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
        self.label_11.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.mMapLayerComboBox_selectGrid = QgsMapLayerComboBox()

        self.label_38 = QLabel(_tr("Load Grid"))
        self.label_38.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
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

        # ── Fixed bottom: Apply Processing row ──────────────────────
        bot = QWidget()
        bl = QHBoxLayout(bot)
        bl.setContentsMargins(0, 2, 0, 0)
        bl.setSpacing(6)

        self.pushButton_3_applyProcessing = QPushButton(_tr("Apply Processing"))
        self.pushButton_3_applyProcessing.setStyleSheet(
            "font-weight: bold;")

        self.label_37 = QLabel(_tr("Max FFT Buffer (pixels)"))
        self.label_37.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.lineEdit_13_max_buffer = QLineEdit("5000")
        self.lineEdit_13_max_buffer.setFixedWidth(55)

        self.label_27 = QLabel(_tr("Version"))
        self.label_27.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
        self.version_label = QLabel("0.1")
        self.version_label.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)

        self.label_41_units = QLabel(_tr("Units"))

        bl.addWidget(self.pushButton_3_applyProcessing)
        bl.addStretch()
        bl.addWidget(self.label_37)
        bl.addWidget(self.lineEdit_13_max_buffer)
        bl.addWidget(self.label_27)
        bl.addWidget(self.version_label)
        bl.addWidget(self.label_41_units)

        outer.addWidget(bot)
        return tab

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
        self.label_15.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.lineEdit_6_inc = QLineEdit("0")
        self.lineEdit_6_inc.setFixedWidth(45)

        self.label_14 = QLabel(_tr("Dec"))
        self.label_14.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.lineEdit_5_dec = QLineEdit("0")
        self.lineEdit_5_dec.setFixedWidth(45)

        self.label_33 = QLabel(_tr("Int"))
        self.label_33.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
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
        self.label_26.setAlignment(Qt.AlignCenter)
        self.checkBox_4_PGrav = QCheckBox(
            _tr("Vertical Integration (apply to RTE/P result to get Pseudo Gravity)"))
        g.addWidget(self.label_26,        2, 0)
        g.addWidget(self.checkBox_4_PGrav, 2, 1, 1, 4)

        # Row 3: ● Continuation / Direction / Height
        self.label_28 = QLabel("●")
        self.label_28.setAlignment(Qt.AlignCenter)
        self.checkBox_9_continuation = QCheckBox(_tr("Continuation"))

        self.label_20 = QLabel(_tr("Direction"))
        self.label_20.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.comboBox_2_continuationDirection = QComboBox()
        self.comboBox_2_continuationDirection.setFixedWidth(70)

        self.label_21 = QLabel(_tr("Height (m)"))
        self.label_21.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
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
        gb = QGroupBox(_tr("Frequency Filters"))
        g = QGridLayout(gb)
        g.setSpacing(4)

        # Row 0: Radial Power Spectrum button (right-aligned)
        self.pushButton_rad_power_spectrum = QPushButton(_tr("Radial Power Spectrum"))
        row0 = QHBoxLayout()
        row0.addStretch()
        row0.addWidget(self.pushButton_rad_power_spectrum)
        g.addLayout(row0, 0, 0, 1, 9)

        # Row 1: Direction filter
        self.checkBox_3_DirClean = QCheckBox(_tr("Direction Cos/Butterw."))
        self.label_12 = QLabel(_tr("Azimuth"))
        self.label_12.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.lineEdit_3_azimuth = QLineEdit("0")
        self.lineEdit_3_azimuth.setFixedWidth(35)
        self.label_13 = QLabel(_tr("Wavelength"))
        self.label_13.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.lineEdit_3_DC_wavelength = QLineEdit("1000")
        self.lineEdit_3_DC_wavelength.setFixedWidth(55)
        self.label_16 = QLabel(_tr("Scale"))
        self.label_16.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.lineEdit_3_DC_scale = QLineEdit("3")
        self.lineEdit_3_DC_scale.setFixedWidth(35)

        g.addWidget(self.checkBox_3_DirClean,   1, 0, 1, 2)
        g.addWidget(self.label_12,              1, 2)
        g.addWidget(self.lineEdit_3_azimuth,    1, 3)
        g.addWidget(self.label_13,              1, 4)
        g.addWidget(self.lineEdit_3_DC_wavelength, 1, 5)
        g.addWidget(self.label_16,              1, 6)
        g.addWidget(self.lineEdit_3_DC_scale,   1, 7)

        # Row 2: Remove Regional
        self.checkBox_5_regional = QCheckBox(_tr("Remove Regional"))
        self.radioButton_RR_1st = QRadioButton(_tr("1st order fit"))
        self.radioButton_RR_1st.setChecked(True)
        self.radioButton_RR_2nd = QRadioButton(_tr("2nd order fit"))
        g.addWidget(self.checkBox_5_regional,  2, 0, 1, 2)
        g.addWidget(self.radioButton_RR_1st,   2, 3, 1, 2)
        g.addWidget(self.radioButton_RR_2nd,   2, 5, 1, 2)

        # Row 3: Band Pass
        self.label_29 = QLabel("●")
        self.label_29.setAlignment(Qt.AlignCenter)
        self.checkBox_10_bandPass = QCheckBox(_tr("Band Pass"))
        self.label_22 = QLabel(_tr("Low"))
        self.label_22.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.lineEdit_12_bandPassLow = QLineEdit("50000")
        self.lineEdit_12_bandPassLow.setFixedWidth(60)
        self.label_23 = QLabel(_tr("High (proj units)"))
        self.label_23.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.lineEdit_11_bandPassHigh = QLineEdit("5000")
        self.lineEdit_11_bandPassHigh.setFixedWidth(55)
        self.label_39 = QLabel(_tr("Width"))
        self.label_39.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.lineEdit_3_BP_width = QLineEdit("5000")
        self.lineEdit_3_BP_width.setFixedWidth(55)

        g.addWidget(self.label_29,              3, 0)
        g.addWidget(self.checkBox_10_bandPass,  3, 1)
        g.addWidget(self.label_22,              3, 2)
        g.addWidget(self.lineEdit_12_bandPassLow, 3, 3)
        g.addWidget(self.label_23,              3, 4)
        g.addWidget(self.lineEdit_11_bandPassHigh, 3, 5)
        g.addWidget(self.label_39,              3, 6)
        g.addWidget(self.lineEdit_3_BP_width,   3, 7)

        # Row 4: High/Low Pass
        self.label_32 = QLabel("●")
        self.label_32.setAlignment(Qt.AlignCenter)
        self.checkBox_10_freqCut = QCheckBox(_tr("High/Low Pass"))
        self.label_35 = QLabel(_tr("Accept"))
        self.label_35.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.comboBox_2_FreqCutType = QComboBox()
        self.comboBox_2_FreqCutType.setFixedWidth(80)
        self.label_34 = QLabel(_tr("Cutoff (proj units)"))
        self.label_34.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.lineEdit_12_FreqPass = QLineEdit("5000")
        self.lineEdit_12_FreqPass.setFixedWidth(55)
        self.label_69 = QLabel(_tr("Width"))
        self.label_69.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.lineEdit_3_HLP_width = QLineEdit("5000")
        self.lineEdit_3_HLP_width.setFixedWidth(55)

        g.addWidget(self.label_32,              4, 0)
        g.addWidget(self.checkBox_10_freqCut,   4, 1)
        g.addWidget(self.label_35,              4, 2)
        g.addWidget(self.comboBox_2_FreqCutType, 4, 3)
        g.addWidget(self.label_34,              4, 4)
        g.addWidget(self.lineEdit_12_FreqPass,  4, 5)
        g.addWidget(self.label_69,              4, 6)
        g.addWidget(self.lineEdit_3_HLP_width,  4, 7)

        # Row 5: AGC
        self.label_30 = QLabel("●")
        self.label_30.setAlignment(Qt.AlignCenter)
        self.checkBox_11_1vd_agc = QCheckBox(_tr("Automatic Gain Control"))
        self.label_24 = QLabel(_tr("Window (pixels)"))
        self.label_24.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.lineEdit_13_agc_window = QLineEdit("10")
        self.lineEdit_13_agc_window.setFixedWidth(45)

        g.addWidget(self.label_30,              5, 0)
        g.addWidget(self.checkBox_11_1vd_agc,   5, 1, 1, 2)
        g.addWidget(self.label_24,              5, 4)
        g.addWidget(self.lineEdit_13_agc_window, 5, 5)

        g.setColumnStretch(8, 1)
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
        self.label_25.setAlignment(Qt.AlignCenter)
        self.checkBox_6_derivative = QCheckBox(_tr("Derivative"))
        self.label_18 = QLabel(_tr("Direction"))
        self.label_18.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.comboBox_derivDirection = QComboBox()
        self.comboBox_derivDirection.setFixedWidth(60)
        self.label_19 = QLabel(_tr("Power"))
        self.label_19.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
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
        self.label_40.setAlignment(Qt.AlignCenter)
        self.checkBox_11_tot_hz_grad = QCheckBox(_tr("Total Horizontal Gradient"))
        g.addWidget(self.label_40,               1, 0)
        g.addWidget(self.checkBox_11_tot_hz_grad, 1, 1, 1, 3)

        # Row 2: Tilt Angle
        self.label_41 = QLabel("●")
        self.label_41.setAlignment(Qt.AlignCenter)
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
        self.label_47.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.mMapLayerComboBox_selectGrid_Conv = QgsMapLayerComboBox()
        tl.addWidget(self.label_47)
        tl.addWidget(self.mMapLayerComboBox_selectGrid_Conv, 1)
        outer.addWidget(top)

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

        # ── Fixed bottom: Apply Processing ───────────────────────────
        bot = QWidget()
        bl = QHBoxLayout(bot)
        bl.setContentsMargins(0, 2, 0, 0)
        self.pushButton_3_applyProcessing_Conv = QPushButton(_tr("Apply Processing"))
        self.pushButton_3_applyProcessing_Conv.setStyleSheet(
            "font-weight: bold;")
        bl.addWidget(self.pushButton_3_applyProcessing_Conv)
        bl.addStretch()
        outer.addWidget(bot)
        return tab

    # ------------------------------------------------------------------
    # Group: Convolution Filters
    # ------------------------------------------------------------------
    def _group_conv_filters(self):
        gb = QGroupBox(_tr("Convolution Filters"))
        g = QGridLayout(gb)
        g.setSpacing(4)

        hdr_lbl = QLabel(_tr("Filter size (odd numbers)"))
        hdr_lbl.setAlignment(Qt.AlignRight | Qt.AlignVCenter)

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
        lbl_med.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.lineEdit_Median_size = QLineEdit("3")
        self.lineEdit_Median_size.setFixedWidth(50)
        g.addWidget(self.checkBox_Median,     1, 0)
        g.addWidget(lbl_med,                  1, 1)
        g.addWidget(self.lineEdit_Median_size, 1, 2)

        # Gaussian
        self.checkBox_Gaussian = QCheckBox(_tr("Gaussian"))
        self.label_65 = QLabel(_tr("Sigma"))
        self.label_65.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.lineEdit_Gaussian_Sigma = QLineEdit("1")
        self.lineEdit_Gaussian_Sigma.setFixedWidth(50)
        g.addWidget(self.checkBox_Gaussian,   2, 0)
        g.addWidget(self.label_65,            2, 1)
        g.addWidget(self.lineEdit_Gaussian_Sigma, 2, 2)

        # Directional
        self.checkBox_Directional = QCheckBox(_tr("Directional"))
        self.label_36 = QLabel(_tr("Direction"))
        self.label_36.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.comboBox_Dir_dir = QComboBox()
        self.comboBox_Dir_dir.setFixedWidth(70)
        g.addWidget(self.checkBox_Directional, 3, 0)
        g.addWidget(self.label_36,             3, 1)
        g.addWidget(self.comboBox_Dir_dir,     3, 2)

        # Sun Shading
        self.checkBox_SunShading = QCheckBox(_tr("Sun Shading"))
        self.label_64 = QLabel(_tr("Azimuth"))
        self.label_64.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.lineEdit_SunSh_Az = QLineEdit("-45")
        self.lineEdit_SunSh_Az.setFixedWidth(50)
        self.label_68 = QLabel(_tr("Zenith"))
        self.label_68.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
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
        lbl_ws.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
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
        lbl_aniso_thr.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.doubleSpinBox_SS_AnisoThresh = QDoubleSpinBox()
        self.doubleSpinBox_SS_AnisoThresh.setRange(0.0, 1.0)
        self.doubleSpinBox_SS_AnisoThresh.setSingleStep(0.05)
        self.doubleSpinBox_SS_AnisoThresh.setValue(0.3)
        self.doubleSpinBox_SS_AnisoThresh.setFixedWidth(55)

        lbl_ang_tol = QLabel(_tr("Angle tol°:"))
        lbl_ang_tol.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.doubleSpinBox_SS_AngleTol = QDoubleSpinBox()
        self.doubleSpinBox_SS_AngleTol.setRange(1.0, 90.0)
        self.doubleSpinBox_SS_AngleTol.setSingleStep(5.0)
        self.doubleSpinBox_SS_AngleTol.setValue(22.5)
        self.doubleSpinBox_SS_AngleTol.setFixedWidth(55)

        lbl_steps = QLabel(_tr("Max steps:"))
        lbl_steps.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.spinBox_SS_MaxSteps = QSpinBox()
        self.spinBox_SS_MaxSteps.setRange(10, 2000)
        self.spinBox_SS_MaxSteps.setSingleStep(50)
        self.spinBox_SS_MaxSteps.setValue(200)
        self.spinBox_SS_MaxSteps.setFixedWidth(55)

        lbl_search_r = QLabel(_tr("Search r:"))
        lbl_search_r.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
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

        # DTM Curvature Classifier
        self.checkBox_DTM_Class = QCheckBox(_tr("DTM Curvature Classifier"))
        self.label_74 = QLabel(_tr("Curve Threshold"))
        self.label_74.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.lineEdit_DTM_Curve = QLineEdit(".0001")
        self.lineEdit_DTM_Curve.setFixedWidth(55)
        self.label_73 = QLabel(_tr("Cliff Angle"))
        self.label_73.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.lineEdit_DTM_Cliff = QLineEdit("5")
        self.lineEdit_DTM_Cliff.setFixedWidth(40)
        self.label_75 = QLabel(_tr("Sigma"))
        self.label_75.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.lineEdit_DTM_Sigma = QLineEdit("0")
        self.lineEdit_DTM_Sigma.setFixedWidth(40)

        g.addWidget(self.checkBox_DTM_Class,  7, 0, 1, 2)
        g.addWidget(self.label_74,            7, 2)
        g.addWidget(self.lineEdit_DTM_Curve,  7, 3)
        g.addWidget(self.label_73,            7, 4)
        g.addWidget(self.lineEdit_DTM_Cliff,  7, 5)
        g.addWidget(self.label_75,            8, 2)
        g.addWidget(self.lineEdit_DTM_Sigma,  8, 3)

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
        lbl_comp.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        lbl_comp2 = QLabel(_tr("Number of Components (0=all)"))
        lbl_comp2.setAlignment(Qt.AlignRight | Qt.AlignVCenter)

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
        self.label_79.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.doubleSpinBox_ED_Threshold = QDoubleSpinBox()
        self.doubleSpinBox_ED_Threshold.setDecimals(5)
        self.doubleSpinBox_ED_Threshold.setMinimum(0.00001)
        self.doubleSpinBox_ED_Threshold.setMaximum(1.0)
        self.doubleSpinBox_ED_Threshold.setSingleStep(0.1)
        self.doubleSpinBox_ED_Threshold.setValue(0.1)
        self.label_80 = QLabel(_tr("Window (odd)"))
        self.label_80.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
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

        sc = QWidget()
        sl = QVBoxLayout(sc)
        sl.setContentsMargins(6, 6, 6, 6)
        sl.setSpacing(6)
        sl.addWidget(self._group_import_points())
        sl.addWidget(self._group_gridding())
        sl.addWidget(self._group_worm_header())
        sl.addWidget(self._group_worms())
        sl.addWidget(self._group_wtmm())
        sl.addStretch()

        outer.addWidget(self._scroll_wrap(sc))
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
        self.label_53.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
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
        self.label_44.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.comboBox_grid_x = QComboBox()
        self.comboBox_grid_x.setEnabled(False)
        g.addWidget(self.label_44,       2, 0)
        g.addWidget(self.comboBox_grid_x, 2, 1, 1, 2)

        # Row 3: Y column
        self.label_45 = QLabel("lat_y")
        self.label_45.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
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
        self.label_54.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.mMapLayerComboBox_selectGrid_3 = QgsMapLayerComboBox()
        self.label_55 = QLabel(_tr("Select Data Field"))
        self.label_55.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.comboBox_select_grid_data_field = QComboBox()

        g.addWidget(self.label_54,                       0, 0)
        g.addWidget(self.mMapLayerComboBox_selectGrid_3, 0, 1, 1, 3)
        g.addWidget(self.label_55,                       1, 0)
        g.addWidget(self.comboBox_select_grid_data_field, 1, 1, 1, 3)

        # Cell size + # cells
        self.label_77 = QLabel(_tr("Cell Size\n(proj Units)"))
        self.label_77.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.doubleSpinBox_cellsize = QDoubleSpinBox()
        self.doubleSpinBox_cellsize.setDecimals(6)
        self.doubleSpinBox_cellsize.setMinimum(0.000001)
        self.doubleSpinBox_cellsize.setMaximum(9999999.99)
        self.doubleSpinBox_cellsize.setSingleStep(10.0)
        self.doubleSpinBox_cellsize.setValue(100.0)

        self.label_51 = QLabel(_tr("# Cells"))
        self.label_51.setStyleSheet("font-weight: bold;")
        self.label_49 = QLabel("X")
        self.label_49.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.nx_label = QLabel("0")
        self.label_48 = QLabel("Y")
        self.label_48.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.ny_label = QLabel("0")

        g.addWidget(self.label_77,            2, 0)
        g.addWidget(self.doubleSpinBox_cellsize, 2, 1)
        g.addWidget(self.label_51,            2, 2)
        g.addWidget(self.label_49,            2, 3)
        g.addWidget(self.nx_label,            2, 4)
        g.addWidget(self.label_48,            3, 3)
        g.addWidget(self.ny_label,            3, 4)

        # Gridding buttons
        self.pushButton_idw_2 = QPushButton(_tr("IDW Gridding"))
        self.pushButton_idw_2.setStyleSheet("font-weight: bold;")
        self.pushButton_bspline_3 = QPushButton(_tr("BSpline Gridding"))
        self.pushButton_bspline_3.setStyleSheet("font-weight: bold;")

        g.addWidget(self.pushButton_idw_2,    4, 1)
        g.addWidget(self.pushButton_bspline_3, 4, 2, 1, 2)

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
        self.label_52.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
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
        self.label_58.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.spinBox_levels = QSpinBox()
        self.spinBox_levels.setMinimum(1)
        self.spinBox_levels.setValue(10)
        self.checkBox_worms_shp = QCheckBox(_tr("Also save to shapefile"))

        self.label_60 = QLabel(_tr("Base Level"))
        self.label_60.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.doubleSpinBox_base = QDoubleSpinBox()
        self.doubleSpinBox_base.setDecimals(0)
        self.doubleSpinBox_base.setMinimum(0)
        self.doubleSpinBox_base.setMaximum(100000)
        self.doubleSpinBox_base.setSingleStep(100)
        self.doubleSpinBox_base.setValue(1000)

        self.label_62 = QLabel(_tr("Increment"))
        self.label_62.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
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
        self.label_57.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.mMapLayerComboBox_selectVectors = QgsMapLayerComboBox()

        self.label_59 = QLabel(_tr("Select feature"))
        self.label_59.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.mFieldComboBox_feature = QComboBox()

        self.label_76 = QLabel(_tr("Select data field"))
        self.label_76.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.mFieldComboBox_data = QgsFieldComboBox()

        self.label_61 = QLabel(_tr("Spacing\n(Proj Units)"))
        self.label_61.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
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
    # TAB 4 – Utils
    # ------------------------------------------------------------------
    def _tab_utils(self):
        tab = QWidget()
        tab.setObjectName("Utils")
        tab.setStyleSheet("#Utils { " + TAB_BG + " }")

        outer = QVBoxLayout(tab)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        sc = QWidget()
        sl = QVBoxLayout(sc)
        sl.setContentsMargins(6, 6, 6, 6)
        sl.setSpacing(6)

        # Select Grid (top of Utils)
        top_row = QWidget()
        tl = QHBoxLayout(top_row)
        tl.setContentsMargins(0, 0, 0, 0)
        self.label_56 = QLabel(_tr("Select Grid"))
        self.label_56.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.mMapLayerComboBox_selectGrid_Conv_2 = QgsMapLayerComboBox()
        tl.addWidget(self.label_56)
        tl.addWidget(self.mMapLayerComboBox_selectGrid_Conv_2, 1)
        sl.addWidget(top_row)

        sl.addWidget(self._group_nan_threshold())
        sl.addWidget(self._group_clip_polygon())

        # Apply Processing
        apply_row = QWidget()
        al = QHBoxLayout(apply_row)
        al.setContentsMargins(0, 0, 0, 0)
        self.pushButton_3_applyProcessing_Conv_3 = QPushButton(_tr("Apply Processing"))
        self.pushButton_3_applyProcessing_Conv_3.setStyleSheet("font-weight: bold;")
        al.addWidget(self.pushButton_3_applyProcessing_Conv_3)
        al.addStretch()
        sl.addWidget(apply_row)

        sl.addWidget(self._group_normalise())
        sl.addWidget(self._group_lut_convert())
        sl.addStretch()

        outer.addWidget(self._scroll_wrap(sc))
        return tab

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
        self.label_70.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.doubleSpinBox_NaN_Above = QDoubleSpinBox()
        self.doubleSpinBox_NaN_Above.setDecimals(8)
        self.doubleSpinBox_NaN_Above.setMinimum(-1e30)
        self.doubleSpinBox_NaN_Above.setMaximum(1e30)
        self.doubleSpinBox_NaN_Above.setValue(1.0)

        self.label_71 = QLabel(_tr("Below"))
        self.label_71.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
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
    # Group: Normalise Grids
    # ------------------------------------------------------------------
    def _group_normalise(self):
        gb = QGroupBox(_tr("Normalise Grids"))
        self.groupBox_10 = gb
        g = QGridLayout(gb)
        g.setSpacing(4)

        self.label_66 = QLabel(_tr("Input Directory"))
        self.label_66.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.lineEdit_loadPointsPath_normalise_in = QLineEdit()
        self.pushButton_select_normalise_in = QPushButton("...")
        self.pushButton_select_normalise_in.setFixedWidth(30)
        self.radioButton_normalise_1st = QRadioButton(_tr("1st order fit"))
        self.radioButton_normalise_1st.setChecked(True)

        self.label_72 = QLabel(_tr("Output Directory"))
        self.label_72.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
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
        self.label_42.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.lineEdit_2_loadGridPath_2 = QLineEdit()
        self.pushButton_2_selectGrid_RGB = QPushButton("...")
        self.pushButton_2_selectGrid_RGB.setFixedWidth(30)

        g.addWidget(self.label_42,                0, 0)
        g.addWidget(self.lineEdit_2_loadGridPath_2, 0, 1)
        g.addWidget(self.pushButton_2_selectGrid_RGB, 0, 2)

        # CSS colour list
        self.label_43 = QLabel(_tr("CSS Colour List\nor\nRGB triplets"))
        self.label_43.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.textEdit_2_colour_list = QTextEdit()
        self.textEdit_2_colour_list.setPlaceholderText(
            _tr("comma separated list of CSS colour names or comma separated RGB triplets"))
        self.textEdit_2_colour_list.setFixedHeight(70)

        self.pushButton_pick_rgb = QPushButton("Pick RGB")
        self.pushButton_pick_rgb.setToolTip(
            _tr("Pick RGB from map: click button then click a point on the loaded RGB raster"))

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
        self.label_46.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.mQgsDoubleSpinBox_LUT_min = QgsDoubleSpinBox()
        self.mQgsDoubleSpinBox_LUT_min.setMinimum(-1e20)
        self.mQgsDoubleSpinBox_LUT_min.setMaximum(1e20)
        self.mQgsDoubleSpinBox_LUT_min.setSingleStep(100)

        self.label_50 = QLabel(_tr("Max"))
        self.label_50.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
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
    # TAB 5 – Help
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
<p><b>Example magnetic grid from Mauritania:</b></p>
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
<p>- Calcs ChatGPT and Mark Jessell</p>
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
</body></html>"""
