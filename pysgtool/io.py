# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2024-2026 Mark Jessell
"""Reading grids and points."""

import os

import numpy as np

from .calcs import raster_io
from .grid import Grid, Stack


def read(path, band=1, crs=None):
    """Read one band of a raster as a Grid: GeoTIFF, ER Mapper, and whatever else GDAL
    reads, Geosoft .grd and Noddy .mag / .grv. crs: replaces the file's coordinate
    system ('EPSG:28350')."""
    path = os.path.abspath(str(path))
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    if path.lower().endswith((".mag", ".grv")):
        return _read_noddy(path, crs)
    if path.lower().endswith(".grd"):
        grid = _read_geosoft(path)
        if grid is not None:
            if crs is not None:
                grid.crs, grid._geographic = raster_io.crs_to_wkt(crs), None
            return grid
    r = raster_io.read(path, [band])
    grid = Grid(r["stack"][0], gt=r["gt"], crs=crs if crs is not None else r["wkt"],
                name=_stem(path))
    grid.metadata = r["metadata"]
    grid.path, grid._on_disk = path, True
    return grid


def read_stack(path, bands=None):
    """Read all (or some) bands of a raster as a Stack. Cells an alpha band or mask marks
    as invalid are no data."""
    path = os.path.abspath(str(path))
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    r = raster_io.read(path, bands, masks=True)
    numbers = bands or range(1, len(r["names"]) + 1)
    stack = Stack(r["stack"], gt=r["gt"], crs=r["wkt"], name=_stem(path),
                  names=[n or f"Band {b}" for n, b in zip(r["names"], numbers)])
    stack.path, stack._on_disk = path, True
    return stack


def _stem(path):
    return os.path.splitext(os.path.basename(path))[0]


def _read_noddy(path, crs):
    """A Noddy ASCII grid, as the plugin imports it: in UTM 29N unless crs is given,
    with the header (the model's field ...) in .metadata."""
    from .calcs import sgt_metadata
    from .calcs.filter_table import _plugin_version
    from .calcs.line_data import NODDY_EPSG, noddy_output_name, read_noddy

    noddy = read_noddy(path)
    out = noddy_output_name(path, noddy["data_type"])
    grid = Grid(noddy["values"], gt=noddy["gt"],
                crs=crs if crs is not None else f"EPSG:{NODDY_EPSG}", name=_stem(out))
    grid.metadata = dict(noddy["metadata"])
    grid._file_items = dict(noddy["metadata"])
    grid.path = out
    grid._provenance = sgt_metadata._build(
        out, path, "Import Noddy grid", dict(noddy["metadata"], data_type=noddy["data_type"]),
        _plugin_version(), file_stats=False)
    return grid


def _read_geosoft(path):
    """A Geosoft grid, or None if the .grd file is of another kind."""
    from .calcs.geosoft_grid_parser import extract_proj_str, load_oasis_montaj_grid_optimized

    with open(path, "rb") as f:
        head = f.read(5 * 4)
    if len(head) < 20:
        return None
    size, _sign, _ne, _nv, order = np.frombuffer(head, dtype="<i4")
    if size not in (1, 2, 4, 8, 1025, 1026, 1028, 1032) or order not in (-1, 1):
        return None  # e.g. Surfer or netCDF, which GDAL reads
    values, header, _type = load_oasis_montaj_grid_optimized(path)
    if not isinstance(values, np.ndarray):
        raise ValueError("Geosoft grids of byte or short integers cannot be read")
    if header["rotation"]:
        raise ValueError("Rotated Geosoft grids are not supported")
    # KX = 1: each vector is a row; KX = -1: each vector is a column
    dx, dy = header["spacing_e"], header["spacing_v"]
    if header["ordering"] == -1:
        dx, dy = dy, dx
    values = np.flipud(values)  # row 0 of the file is the south
    gt = (header["x_origin"] - dx / 2, dx, 0.0,
          header["y_origin"] - dy / 2 + values.shape[0] * dy, 0.0, -dy)
    code = extract_proj_str(path + ".xml") if os.path.exists(path + ".xml") else None
    grid = Grid(values, gt=gt, crs=f"EPSG:{code}" if code else None, name=_stem(path))
    grid.path, grid._on_disk = path, True
    return grid


class Points:
    """Points with fields: x, y, fields {name: array}, crs. grid(field) interpolates them."""

    def __init__(self, x, y, fields, crs=None, name="points"):
        self.x, self.y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
        self.fields = {n: np.asarray(v) for n, v in fields.items()}
        self.crs = raster_io.crs_to_wkt(crs)
        self.name = name
        self.path = None
        self._import = None  # (operation, parameters) of the plugin's import, for save()

    def __len__(self):
        return len(self.x)

    def __repr__(self):
        return f"<Points {self.name!r}: {len(self)} points, fields {', '.join(self.fields)}>"

    def grid(self, field, cell_size=None, **kwargs):
        """Multilevel B-spline grid of a field: see pysgtool.grid_points."""
        from .gridding import grid_points
        from .calcs.Petrophysics import to_number

        if field not in self.fields:
            raise ValueError(f"There is no field {field}: {', '.join(self.fields)}")
        z = np.array([to_number(v) for v in self.fields[field]])
        out = grid_points(self.x, self.y, z, cell_size, crs=self.crs, field=field,
                          name=f"{self.name}_{field}_bspline", source=self.path, **kwargs)
        return out

    def profile(self, field, line=None, spacing="auto"):
        """One line's points resampled to an even spacing, as the plugin's WTMM does: put
        in order along x or y (whichever spans more) and interpolated linearly. line:
        the LINE_ID (or LINE) of the points to take, None for all. spacing: a distance
        or "auto" (the median distance between the points).
        Returns (coords (n, 2), values, spacing used)."""
        from .calcs.line_data import regularize_profile

        keep = np.ones(len(self), dtype=bool)
        if line is not None:
            key = next((k for k in ("LINE_ID", "LINE") if k in self.fields), None)
            if key is None:
                raise ValueError("The points have no LINE_ID or LINE field")
            keep = self.fields[key] == line
            if not keep.any():
                raise ValueError(f"No points have {key} = {line}")
        if field not in self.fields:
            raise ValueError(f"There is no field {field}: {', '.join(self.fields)}")
        coords = np.column_stack([self.x[keep], self.y[keep]])
        return regularize_profile(coords, self.fields[field][keep].astype(float), spacing)

    def save(self, path=None):
        """Save as a point file (.shp, .gpkg or .geojson); by default <file>.shp beside the
        file read, as the plugin does. Returns the path."""
        import struct

        from .calcs import sgt_metadata
        from .calcs.filter_table import _plugin_version

        if path is None:
            if self.path is None:
                raise ValueError("Give a path: the points did not come from a file")
            path = self.path + ".shp"
        path = os.path.abspath(str(path))
        wkb = [struct.pack("<BIdd", 1, 1, px, py) for px, py in zip(self.x, self.y)]
        raster_io.write_vector(path, wkb, "Point", self.fields, self.crs)
        operation, parameters = self._import or ("Points", {})
        sgt_metadata.write_sgt_metadata(path, self.path, operation, parameters,
                                        _plugin_version())
        return path


def read_points(path, x=None, y=None, crs=None, load_ties=True):
    """Read points as Points from a vector file (.shp, .gpkg, .geojson ...), an
    ASEG-GDF2 .dat file (its .dfn beside it), delimited text (.csv, .txt) or an XYZ line
    file (.xyz). x, y: the coordinate fields of a .dat, .csv or .txt file. The points of
    an XYZ file have fields LINE_ID, data_0, data_1 ... (the values after x and y);
    load_ties: keep its tie lines."""
    path = os.path.abspath(str(path))
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    ext = os.path.splitext(path)[1].lower()
    if ext == ".dat":
        points = _read_aseg(path, x, y, crs)
    elif ext == ".xyz":
        points = _read_xyz(path, crs, load_ties)
    elif ext in (".csv", ".txt"):
        points = _read_delimited(path, x, y, crs)
    else:
        import struct

        v = raster_io.read_vector(path)
        xs, ys, rows = [], [], []
        for i, wkb in enumerate(v["wkb"]):
            for px, py in _wkb_points(wkb, struct):
                xs.append(px)
                ys.append(py)
                rows.append(i)
        if not xs:
            raise ValueError(f"{path} has no points")
        fields = {n: _column([values[i] for i in rows]) for n, values in v["fields"].items()}
        points = Points(xs, ys, fields, crs if crs is not None else v["wkt"], _stem(path))
    points.path = path
    return points


def _column(values):
    """Field values as numbers (NaN for none) if they all are, else as they are."""
    if all(v is None or isinstance(v, (int, float, np.integer, np.floating)) for v in values):
        return np.array([np.nan if v is None else v for v in values])
    return np.array(values, dtype=object)


def _wkb_points(wkb, struct):
    """The (x, y) of a point or multipoint geometry ([] for anything else or empty)."""
    if wkb is None:
        return []
    wkb = bytes(wkb)
    order = "<" if wkb[0] == 1 else ">"
    kind = struct.unpack(order + "I", wkb[1:5])[0]
    flags, kind = kind & 0xE0000000, kind & 0x1FFFFFFF
    base = kind % 1000
    dims = 2 + bool(flags & 0x80000000) + bool(flags & 0x40000000)
    if kind >= 1000:
        dims += 2 if kind >= 3000 else 1
    if base == 1:
        px, py = struct.unpack(order + "dd", wkb[5:21])
        return [(px, py)] if np.isfinite(px) and np.isfinite(py) else []
    if base == 4:
        n, step, out = struct.unpack(order + "I", wkb[5:9])[0], 5 + 8 * dims, []
        for i in range(n):
            out += _wkb_points(wkb[9 + i * step: 9 + (i + 1) * step], struct)
        return out
    return []


def _read_xyz(path, crs, load_ties):
    from .calcs.line_data import read_xyz

    rows = read_xyz(path, load_ties)
    width = len(rows[0])
    # as the plugin: the number of values is that of the first row; missing ones are NaN
    table = np.full((len(rows), width), np.nan)
    for i, row in enumerate(rows):
        x, y, *values, line_id = row
        values = values[:width - 3]
        table[i, :2] = x, y
        table[i, 2:2 + len(values)] = values
        table[i, -1] = line_id
    fields = {"LINE_ID": table[:, -1].astype(np.int64)}
    fields.update({f"data_{i}": table[:, 2 + i] for i in range(width - 3)})
    points = Points(table[:, 0], table[:, 1], fields, crs, _stem(path))
    points._import = ("Import line (XYZ) data", {
        "crs": _crs_label(crs), "layer_name": _stem(path), "load_tie_lines": load_ties})
    return points


def _read_delimited(path, x, y, crs):
    from .calcs.line_data import read_delimited

    names, columns = read_delimited(path)
    if x not in columns or y not in columns:
        raise ValueError(f"Give the x and y coordinate fields, two of: {', '.join(names)}")
    xs, ys = columns[x], columns[y]
    if xs.dtype == object or ys.dtype == object:
        raise ValueError(f"The fields {x} and {y} are not all numbers")
    keep = np.isfinite(xs) & np.isfinite(ys)
    fields = {}
    for n in names:
        v = columns[n][keep]
        # whole numbers as integers, as QGIS reads them
        if v.dtype.kind == "f" and v.size and np.isfinite(v).all() and (v == np.round(v)).all():
            v = v.astype(np.int64)
        fields[n] = v
    points = Points(xs[keep], ys[keep], fields, crs, os.path.basename(path))
    points._import = ("Import points from delimited text",
                      {"x_field": x, "y_field": y, "crs": _crs_label(crs)})
    return points


def _crs_label(crs):
    """'EPSG:nnnn' for a coordinate system given as a number or code, else as given."""
    if crs is None:
        return None
    if isinstance(crs, (int, np.integer)):
        return f"EPSG:{int(crs)}"
    if len(str(crs)) < 30:
        return str(crs)
    code = raster_io.epsg(raster_io.crs_to_wkt(crs))
    return f"EPSG:{code}" if code else None


def _read_aseg(path, x, y, crs):
    from .calcs.aseggdf2parser import AsegGdf2Parser

    base = os.path.splitext(path)[0]
    dfn = next((base + e for e in (".dfn", ".DFN", ".Dfn") if os.path.exists(base + e)), None)
    if dfn is None:
        raise FileNotFoundError(f"{base}.dfn (the definition of {os.path.basename(path)})")
    parser = AsegGdf2Parser()
    header, defs, code = parser.parse_dfn_file(dfn)
    rows = parser.parse_dat_file(path, header, defs)
    if not rows:
        raise ValueError(f"No data could be read from {path}")
    names = list(rows[0])
    if x not in names or y not in names:
        raise ValueError(f"Give the x and y coordinate fields, two of: {', '.join(names)}")
    fields = {n: _column([r.get(n) for r in rows]) for n in names}
    xs, ys = (np.array([np.nan if not isinstance(v, (int, float)) else v for v in
                        (r.get(n) for r in rows)], dtype=float) for n in (x, y))
    keep = np.isfinite(xs) & np.isfinite(ys)
    points = Points(xs[keep], ys[keep], {n: v[keep] for n, v in fields.items()},
                    crs if crs is not None else (f"EPSG:{code}" if code else None), _stem(path))
    points._import = ("Import ASEG-GDF2 (.dat) points",
                      {"x_field": x, "y_field": y, "epsg": code or _crs_label(crs)})
    return points
