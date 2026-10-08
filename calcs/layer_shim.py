"""
A plain-Python stand-in for the few facts calculation code reads from a QGIS
raster layer (its extent and CRS code).

QGIS layers must only be touched on the main thread, but some calculation code
(for example the worms processing, which writes a padded grid) takes a layer to
georeference its output. Background tasks pass one of these instead, built on
the main thread from the layer, so the worker never touches QGIS.
"""


class LayerShim:
    """Quacks like the parts of QgsRasterLayer the calculation code uses:
    layer.dataProvider().extent() -> width()/height()/xMinimum()/yMaximum()...
    and layer.crs().authid()."""

    def __init__(self, extent, authid):
        """extent = (xmin, ymin, xmax, ymax); authid e.g. 'EPSG:28350'."""
        self._xmin, self._ymin, self._xmax, self._ymax = extent
        self._authid = authid

    # QgsRasterLayer / QgsRasterDataProvider
    def dataProvider(self):
        return self

    def extent(self):
        return self

    def crs(self):
        return self

    # QgsRectangle
    def xMinimum(self):
        return self._xmin

    def xMaximum(self):
        return self._xmax

    def yMinimum(self):
        return self._ymin

    def yMaximum(self):
        return self._ymax

    def width(self):
        return self._xmax - self._xmin

    def height(self):
        return self._ymax - self._ymin

    # QgsCoordinateReferenceSystem
    def authid(self):
        return self._authid
