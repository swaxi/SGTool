"""
Processing algorithm: apply the steps in a grid's processing history to another
grid (the same as the dialog's "Apply the same steps to another grid" button).
"""

import os
from datetime import date as _date

from qgis.core import (
    QgsProcessingContext,
    QgsProcessingException,
    QgsProcessingOutputRasterLayer,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterFile,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterString,
)

from .calcs.sgt_metadata import read_sgt_metadata
from .calcs.workflow import history_steps
from .sgt_processing import _SGToolAlgorithm, GROUP_UTIL
from .sgt_workflow import (
    ReplayError, plan_steps, run_replay, output_paths, final_output_path,
    uniform_field, check_igrf_possible,
)
from .calcs.sgt_cancel import OperationCancelled

REPLAY_SPEC = dict(
    id="replay_processing_history", name="Apply the steps of a processing history to a grid",
    group=GROUP_UTIL,
    help=(
        "Repeats the processing steps recorded in one grid's SGTool history on another "
        "grid, each step on the result of the one before. Reduction to the pole or "
        "equator and differential RTP get a new inclination and declination from the "
        "IGRF model for the input grid's location and the survey date. Every step is "
        "saved next to the input with the usual names (grid_RTP.tif, grid_RTP_d1z.tif) "
        "and the last one is added to the project. Steps that cannot be repeated "
        "(gridding, Euler, component analysis) are skipped."
    ),
)


def _read_history(path):
    """Provenance of a grid, or of an exported .sgt.xml file."""
    if path.lower().endswith(".xml"):
        import xml.etree.ElementTree as ET

        try:
            return ET.parse(path).getroot()
        except (OSError, ET.ParseError) as e:
            raise QgsProcessingException(f"Could not read {path}: {e}")
    return read_sgt_metadata(path)


class ReplayHistoryAlgorithm(_SGToolAlgorithm):
    INPUT = "INPUT"
    HISTORY = "HISTORY"
    DATE = "DATE"
    OUTPUT = "OUTPUT"

    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterRasterLayer(self.INPUT, "Grid to apply the steps to")
        )
        self.addParameter(
            QgsProcessingParameterFile(
                self.HISTORY,
                "Grid (or exported .sgt.xml file) whose processing history to repeat",
            )
        )
        self.addParameter(
            QgsProcessingParameterString(
                self.DATE,
                "Survey date for the IGRF, yyyy-mm-dd (empty = the date recorded in the "
                "history, else today)",
                "", False, True,
            )
        )
        self.addOutput(QgsProcessingOutputRasterLayer(self.OUTPUT, "Final grid"))

    def shortHelpString(self):
        return self.spec["help"]

    def helpString(self):
        return self.spec["help"]

    def prepareAlgorithm(self, parameters, context, feedback):
        """Main thread: work out the recipe and release layers it will overwrite."""
        layer = self.parameterAsRasterLayer(parameters, self.INPUT, context)
        if layer is None:
            raise QgsProcessingException("Select the grid to apply the steps to")
        history_path = self.parameterAsFile(parameters, self.HISTORY, context)
        root = _read_history(history_path)
        if root is None:
            raise QgsProcessingException(f"No SGTool processing history found in {history_path}")
        _origin, steps = history_steps(root)
        self._plan = [p for p in plan_steps(steps) if p["replayable"]]
        skipped = [p for p in plan_steps(steps) if not p["replayable"]]
        for p in skipped:
            feedback.pushWarning(f"Skipping {p['operation']}: {p['note']}")
        if not self._plan:
            raise QgsProcessingException("The history has no steps that can be repeated")

        text = (self.parameterAsString(parameters, self.DATE, context) or "").strip()
        if text:
            try:
                self._when = _date.fromisoformat(text)
            except ValueError:
                raise QgsProcessingException("The survey date must be written yyyy-mm-dd")
        else:
            recorded = next((p["date"] for p in self._plan if p["date"] is not None), None)
            self._when = recorded or _date.today()

        self._source = layer.source().split("|")[0]
        ext = layer.extent()
        self._grid_info = {
            "geographic": layer.crs().isGeographic(),
            "authid": layer.crs().authid(),
            "extent": (ext.xMinimum(), ext.yMinimum(), ext.xMaximum(), ext.yMaximum()),
        }
        needs_field = any(
            p["igrf"] for p in self._plan
        )
        if needs_field and uniform_field(self._source) is None:
            why = check_igrf_possible(self._grid_info)
            if why:
                raise QgsProcessingException(why)

        paths = output_paths(self._source, self._plan)
        project = context.project()
        if project is not None:
            wanted = {os.path.normcase(os.path.abspath(p)) for p in paths}
            for lyr in list(project.mapLayers().values()):
                if os.path.normcase(os.path.abspath(lyr.source().split("|")[0])) in wanted:
                    project.removeMapLayer(lyr.id())
        import gc

        gc.collect()
        return True

    def processAlgorithm(self, parameters, context, feedback):
        class _Task:  # lets run_replay report progress and honour Cancel
            def check(_self, percent=None):
                if percent is not None:
                    feedback.setProgress(percent)
                if feedback.isCanceled():
                    raise OperationCancelled()

        try:
            result = run_replay(
                self._source, self._plan, self._grid_info, when=self._when,
                task=_Task(),
            )
        except OperationCancelled:
            return {}
        except ReplayError as e:
            raise QgsProcessingException(str(e))
        for line in result["log"]:
            feedback.pushInfo(line)
        final = result["final"]
        project = context.project()
        if project is not None:
            for path in result["outputs"]:  # every result kept, the last on top
                context.addLayerToLoadOnCompletion(
                    path,
                    QgsProcessingContext.LayerDetails(
                        os.path.splitext(os.path.basename(path))[0], project, self.OUTPUT
                    ),
                )
        feedback.pushInfo(f"{result['n_steps']} steps applied; result {final}")
        return {self.OUTPUT: final}
