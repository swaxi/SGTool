"""
Stand-ins for the osgeo modules where GDAL is not installed (pysgtool from pip).

The calculations work on arrays and do not need GDAL; only the functions that read
or write files do. With these stand-ins a module imports without GDAL, and only a
call that really needs it fails, with a message saying so.
"""


class _Missing:
    def __init__(self, name):
        self._name = name

    def __getattr__(self, attr):
        if attr.startswith("__"):
            raise AttributeError(attr)
        raise ImportError(
            f"osgeo.{self._name} (GDAL) is needed for this function and is not installed. "
            "The pysgtool functions read and write files without it: "
            "pip install pysgtool[io]"
        )

    def __bool__(self):
        return False


gdal = _Missing("gdal")
ogr = _Missing("ogr")
osr = _Missing("osr")
gdalnumeric = _Missing("gdalnumeric")
gdalconst = _Missing("gdalconst")
