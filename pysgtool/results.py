# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2024-2026 Mark Jessell
"""What the calculations that do not give a grid return: points, lines and polygons."""

import os

import numpy as np

from .calcs import raster_io, sgt_metadata
from .calcs.filter_table import _plugin_version, auto_output_path
from .calcs.text_io import savetxt, write_rows


class _Result:
    def __init__(self, grid, operation, parameters):
        self.crs = grid.crs
        self._operation, self._parameters = operation, parameters
        self._source, self._known = grid._source()
        self._beside = grid.path

    def _path(self, path, suffix, ext):
        if path is None:
            if self._beside is None:
                raise ValueError("Give a path: the data did not come from a file")
            path = auto_output_path(self._beside, suffix, ext)
        return os.path.abspath(str(path))

    def _record(self, path):
        sgt_metadata.write_sgt_metadata(
            path, self._source, self._operation, self._parameters, _plugin_version(),
            source_provenance=self._known)
        return path


class Features(_Result):
    """Lines or polygons: .geometries (shapely), .fields and save()."""

    def __init__(self, geometries, geometry_type, fields, grid, operation, parameters,
                 suffix="_boundary"):
        super().__init__(grid, operation, parameters)
        self.geometries, self.geometry_type, self.fields = list(geometries), geometry_type, fields
        self._suffix = suffix

    def __len__(self):
        return len(self.geometries)

    def save(self, path=None):
        """Save as .shp, .gpkg or .geojson; returns the path."""
        path = self._path(path, self._suffix, ".shp")
        raster_io.write_vector(path, [g.wkb for g in self.geometries], self.geometry_type,
                               self.fields, self.crs)
        return self._record(path)


class EulerSolutions(_Result):
    """Euler solutions: x, y, depth (positive down), elevation (-depth), base_level and
    std_dfdz arrays, largest std_dfdz first."""

    HEADER = "y_source, x_source, z_source, base_level, std_dfdz"

    def __init__(self, grid, table, structural_index, window, keep):
        super().__init__(grid, "Euler deconvolution", {
            "structural_index": structural_index, "window_size": window,
            "fraction_of_solutions_kept": keep, "columns": self.HEADER,
        })
        self.table = table  # x, y, depth, base level, std of df/dz
        self.structural_index = structural_index

    x = property(lambda self: self.table[:, 0])
    y = property(lambda self: self.table[:, 1])
    depth = property(lambda self: self.table[:, 2])
    elevation = property(lambda self: -self.table[:, 2])
    base_level = property(lambda self: self.table[:, 3])
    std_dfdz = property(lambda self: self.table[:, 4])

    def __len__(self):
        return len(self.table)

    def save(self, path=None):
        """Save as the plugin's CSV (z as elevation, so negative below the surface);
        returns the path."""
        si = 0 if self.structural_index == 0.001 else f"{self.structural_index:g}"
        path = self._path(path, f"_estimates_SI_{si}", ".csv")
        savetxt(path, np.column_stack([self.y, self.x, self.elevation, self.base_level,
                                       self.std_dfdz]), header=self.HEADER)
        return self._record(path)


class Worms(_Result):
    """Worm points: x, y, z (height of the level) and value arrays."""

    def __init__(self, grid, table, parameters):
        super().__init__(grid, "Worms (BSDWormer) points", dict(parameters, shapefile=False))
        self.table = table  # y, x, z, value
        self._grid = grid

    x = property(lambda self: self.table[:, 1])
    y = property(lambda self: self.table[:, 0])
    z = property(lambda self: self.table[:, 2])
    value = property(lambda self: self.table[:, 3])

    def __len__(self):
        return len(self.table)

    def save(self, path=None):
        """Save the points as the plugin's CSV (y, x, z, val); returns the path."""
        path = self._path(path, "_worms", ".csv")
        with open(path, "w") as f:
            f.write("y,x,z,val\n")
            write_rows(f, self.table, "%s")
        return self._record(path)

    def polylines(self, callback=None):
        """The worms joined into 3D lines, as Features with the level in field Z (needs
        scikit-learn and shapely). Levels below 1 are left out, as in the plugin."""
        import shapely

        from .calcs.GeophysicalProcessor import GeophysicalProcessor

        grid = self._grid
        lines = GeophysicalProcessor(grid.dx, grid.dy, 0).points_to_polylines(
            self.table[:, [1, 0, 2]], grid.dx * 1.3, callback)
        geometries = [shapely.LineString([(x, y, float(z)) for x, y in line.coords])
                      for line, z in lines]
        return Features(geometries, "LineString Z", {"Z": [int(z) for _line, z in lines]}, grid,
                        "Worms (BSDWormer) polylines", dict(self._parameters, shapefile=True),
                        suffix="_worms")
