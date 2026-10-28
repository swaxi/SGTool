"""
Whether a QGIS coordinate reference system is in degrees.

QgsCoordinateReferenceSystem.isGeographic() is False for some lat/long grids
QGIS reads without an EPSG code (for example an ER Mapper .ers grid with
Projection "GEODETIC" and CoordinateType LATLONG), although their map units
are degrees. Treated as projected, their cell sizes in degrees would be used as
metres (relief shading, continuation heights, wavelength checks ...).
"""

from qgis.core import Qgis


def _degrees_unit():
    try:
        return Qgis.DistanceUnit.Degrees  # QGIS 3.30+
    except AttributeError:
        from qgis.core import QgsUnitTypes

        return QgsUnitTypes.DistanceDegrees


def is_geographic(crs):
    """True if crs is geographic or its map units are degrees."""
    if crs.isGeographic():
        return True
    try:
        return crs.isValid() and crs.mapUnits() == _degrees_unit()
    except Exception:
        return False
