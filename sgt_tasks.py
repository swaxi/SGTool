"""
Background tasks for SGTool.

Calculations run in a worker thread through QGIS's task manager so the QGIS
interface stays responsive, progress shows in the task manager / status bar,
and the user can cancel (from SGTool's Cancel button or the task manager).

Rule: the work function runs off the main thread, so it must only compute
(numpy / scipy / GDAL file work). Anything that touches QGIS layers, the
project, widgets or the message bar happens beforehand (preparing the job) or
afterwards in the finished callback, both on the main thread.

Cancelling is cooperative: it takes effect the next time the work calls
task.check() (between steps, and inside long loops such as Euler
deconvolution). A single FFT or filter call cannot be interrupted part-way.
"""

import traceback

from qgis.core import QgsTask

from .calcs.sgt_cancel import OperationCancelled


def _can_cancel_flag():
    try:
        return QgsTask.Flag.CanCancel  # QGIS 4 / Qt6 style
    except AttributeError:
        return QgsTask.CanCancel


class SGToolTask(QgsTask):
    """Run work(task) in the background, then on_done(task, ok) on the main thread.

    After the task ends:
      task.result      what work() returned (None if it did not finish)
      task.cancelled   True if the user cancelled it
      task.error       formatted traceback if work() raised, else None
    """

    def __init__(self, description, work, on_done):
        super().__init__(description, _can_cancel_flag())
        self._work = work
        self._on_done = on_done
        self.result = None
        self.error = None
        self.cancelled = False

    # -- worker thread -------------------------------------------------
    def run(self):
        try:
            self.result = self._work(self)
        except OperationCancelled:
            self.cancelled = True
            return False
        except Exception:
            self.error = traceback.format_exc()
            return False
        if self.isCanceled():
            self.cancelled = True
            return False
        return True

    def check(self, percent=None):
        """Call from the work function: report progress and stop if cancelled."""
        if percent is not None:
            self.setProgress(max(0.0, min(100.0, float(percent))))
        if self.isCanceled():
            raise OperationCancelled()

    def callback(self, start=0.0, end=100.0):
        """A progress callback for calculation code: callback(fraction 0..1)
        maps onto start..end percent and raises OperationCancelled if the task
        has been cancelled."""

        def _cb(fraction=None):
            pct = None
            if fraction is not None:
                pct = start + (end - start) * float(fraction)
            self.check(pct)

        return _cb

    # -- main thread ---------------------------------------------------
    def finished(self, ok):
        if self.isCanceled():
            self.cancelled = True
        self._on_done(self, ok)
