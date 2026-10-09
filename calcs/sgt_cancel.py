"""
Cooperative cancellation for long SGTool calculations.

Calculation code (which must not import QGIS) accepts an optional callback and
calls it as it works, e.g. ``callback(fraction_done)``. The callback may raise
OperationCancelled to stop the calculation; whoever started the calculation
(a background task, a Processing algorithm) catches it and discards the
partial result.
"""


class OperationCancelled(Exception):
    """Raised from a progress callback to abandon the running calculation."""
