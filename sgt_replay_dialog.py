"""
The "Apply the same steps to another grid" dialog: shows the steps found in a
grid's history, lets the user choose which to repeat and which grid to apply
them to, and exposes the settings the replay needs (the survey date for the
IGRF) and which intermediate results to keep).
"""

from qgis.PyQt.QtCore import Qt, QDate
from qgis.PyQt.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QTableWidget, QTableWidgetItem,
    QCheckBox, QDateEdit, QDialogButtonBox, QHeaderView,
)
from qgis.core import QgsMapLayerProxyModel
from qgis.gui import QgsMapLayerComboBox

from .sgt_workflow import step_label

def _enum(owner, scope, name):
    """An enum value under either the Qt6 (scoped) or the Qt5 (flat) spelling."""
    try:
        return getattr(getattr(owner, scope), name)
    except AttributeError:
        return getattr(owner, name)


_CHECKABLE = _enum(Qt, "ItemFlag", "ItemIsUserCheckable")
_ENABLED = _enum(Qt, "ItemFlag", "ItemIsEnabled")
_CHECKED = _enum(Qt, "CheckState", "Checked")
_UNCHECKED = _enum(Qt, "CheckState", "Unchecked")
_RESIZE_TO_CONTENTS = _enum(QHeaderView, "ResizeMode", "ResizeToContents")
_STRETCH = _enum(QHeaderView, "ResizeMode", "Stretch")
_OK = _enum(QDialogButtonBox, "StandardButton", "Ok")
_CANCEL = _enum(QDialogButtonBox, "StandardButton", "Cancel")


class ReplayDialog(QDialog):
    """Choose the steps, the target grid and the replay settings."""

    COLUMNS = ["Apply", "Save", "Step", "Settings recorded", "Note"]

    def __init__(self, history_name, plan, default_date, parent=None, tr=lambda s: s):
        super().__init__(parent)
        self.plan = plan
        self.setWindowTitle(tr("Apply the same steps to another grid"))
        self.resize(820, 460)
        layout = QVBoxLayout(self)

        usable = sum(1 for p in plan if p["replayable"])
        layout.addWidget(QLabel(
            tr("History of %s: %d steps, %d can be repeated. Steps run oldest first, "
               "each on the result of the one before.") % (history_name, len(plan), usable)
        ))

        row = QHBoxLayout()
        row.addWidget(QLabel(tr("Apply to grid:")))
        self.target = QgsMapLayerComboBox()
        self.target.setFilters(QgsMapLayerProxyModel.RasterLayer)
        row.addWidget(self.target, 1)
        layout.addLayout(row)

        self.table = QTableWidget(len(plan), len(self.COLUMNS))
        self.table.setHorizontalHeaderLabels([tr(c) for c in self.COLUMNS])
        self.table.verticalHeader().setVisible(False)
        self.table.setWordWrap(True)
        for r, item in enumerate(plan):
            check = QTableWidgetItem()
            if item["replayable"]:
                check.setFlags(_ENABLED | _CHECKABLE)
                check.setCheckState(_CHECKED)
            else:
                check.setFlags(_ENABLED & ~_ENABLED)  # no flags: cannot be ticked
                check.setCheckState(_UNCHECKED)
            self.table.setItem(r, 0, check)
            keep = QTableWidgetItem()
            if item["replayable"]:
                keep.setFlags(_ENABLED | _CHECKABLE)
                keep.setCheckState(_CHECKED)  # keep every intermediate result by default
            else:
                keep.setFlags(_ENABLED & ~_ENABLED)
                keep.setCheckState(_UNCHECKED)
            self.table.setItem(r, 1, keep)
            for c, text in enumerate(
                (step_label(item), item["settings"], item["note"]), start=2
            ):
                cell = QTableWidgetItem(text)
                cell.setFlags(_ENABLED)
                if not item["replayable"]:
                    cell.setForeground(self.palette().mid())  # greyed: cannot be repeated
                self.table.setItem(r, c, cell)
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(0, _RESIZE_TO_CONTENTS)
        header.setSectionResizeMode(1, _RESIZE_TO_CONTENTS)
        header.setSectionResizeMode(2, _RESIZE_TO_CONTENTS)
        header.setSectionResizeMode(3, _STRETCH)
        header.setSectionResizeMode(4, _STRETCH)
        self.table.horizontalHeaderItem(1).setToolTip(tr(
            "Keep this step's result as a file. Untick to delete it once the next step "
            "has been made from it (the last result is always kept)."
        ))
        self.table.resizeRowsToContents()
        layout.addWidget(self.table, 1)

        # IGRF: only shown if a magnetic reduction is among the steps
        self.needs_igrf = any(p["igrf"] for p in plan)
        self.date = QDateEdit()
        self.date.setCalendarPopup(True)
        self.date.setDisplayFormat("yyyy-MM-dd")
        self.date.setDate(default_date or QDate.currentDate())
        self.date.setToolTip(tr("Date of the survey, for the IGRF model"))
        igrf_row = QHBoxLayout()
        igrf_row.addWidget(QLabel(tr(
            "Inclination and declination for the reductions to the pole / equator "
            "are recalculated from IGRF for the new grid. Survey date:")))
        igrf_row.addStretch()
        igrf_row.addWidget(self.date)
        if self.needs_igrf:
            layout.addLayout(igrf_row)

        buttons = QDialogButtonBox(_OK | _CANCEL)
        buttons.button(_OK).setText(tr("Apply steps"))
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    # ------------------------------------------------------------ results
    def selected_items(self):
        """Recipe items ticked in the table, in order, each with "save" set from
        its Save box (whether to keep that step's result as a file)."""
        chosen = []
        for r, item in enumerate(self.plan):
            cell = self.table.item(r, 0)
            if item["replayable"] and cell.checkState() == _CHECKED:
                chosen.append(dict(item, save=self.table.item(r, 1).checkState() == _CHECKED))
        return chosen

    def target_layer(self):
        return self.target.currentLayer()

    def settings(self):
        """The survey date (datetime.date) for the IGRF."""
        return self.date.date().toPyDate()
