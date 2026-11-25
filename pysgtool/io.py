"""Reading grids and points."""

import os

import numpy as np

from .calcs import raster_io
from .grid import Grid, Stack


def read(path, band=1, crs=None):
    """Read one band of a raster as a Grid: GeoTIFF, ER Mapper, and whatever else GDAL
    reads, and Geosoft .grd. crs: replaces the file's coordinate system ('EPSG:28350')."""
    path = os.path.abspath(str(path))
    if not os.path.exists(path):
        raise FileNotFoundError(path)
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


def read_points(path, x=None, y=None, crs=None):
    """Read points from a vector file (.shp, .gpkg, .geojson ...) or an ASEG-GDF2 .dat file
    (its .dfn beside it) as Points. x, y: the coordinate fields of a .dat file."""
    path = os.path.abspath(str(path))
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    if path.lower().endswith(".dat"):
        points = _read_aseg(path, x, y, crs)
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
    return Points(xs[keep], ys[keep], {n: v[keep] for n, v in fields.items()},
                  crs if crs is not None else (f"EPSG:{code}" if code else None), _stem(path))
