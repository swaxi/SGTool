"""
Reading line (XYZ) data, delimited text points and Noddy grids, and profiles
along a line: shared by the QGIS plugin and pysgtool (numpy and scipy only).
"""

import csv
import logging
import os
import re

import numpy as np

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------- XYZ lines
def read_xyz(path, load_ties=True):
    """
    The rows of an XYZ line file: a list of [x, y, value..., line number].

    A row belongs to the line of the last "LINE: 10" / "Line 10" marker above it,
    or of the last "TIE: 900" / "Tie 900" marker if load_ties (otherwise the rows
    of tie lines are left out). Blank rows and rows starting with "/" are
    comments; words in a row (dates, times) are left out, and rows with fewer
    than two numbers are skipped.

    Raises ValueError if no row has coordinates.
    """
    rows = []
    current = None
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("/"):
                continue
            if re.match(r"(?i)^line\b", line):
                current = int(re.search(r"\d+", line).group())
            elif re.match(r"(?i)^tie\b", line):
                current = int(re.search(r"\d+", line).group()) if load_ties else None
            elif current is not None:
                parts = []
                for token in line.split():
                    try:
                        parts.append(float(token))
                    except ValueError:
                        continue
                if len(parts) >= 2:
                    rows.append(parts + [current])
            elif load_ties:
                logger.debug("Invalid line: %s", line)
    if not rows:
        raise ValueError(
            f"No usable coordinate data found in {path}. "
            "Expected 'LINE:'/'Line' (and optionally 'TIE:'/'Tie') markers "
            "followed by rows of numeric X Y [values...]."
        )
    return rows


def xyz_value_count(path):
    """The number of values after x and y in the first data row of an XYZ file,
    or None if it has none."""
    try:
        return len(read_xyz(path)[0]) - 3
    except ValueError:
        return None


def xyz_lines(rows):
    """{line number: [(x, y, [values]), ...]} of the rows of read_xyz, in the
    order the lines first appear."""
    lines = {}
    for row in rows:
        x, y, *values, line_id = row
        lines.setdefault(line_id, []).append((x, y, values))
    return lines


# ---------------------------------------------------------------- delimited text
def read_delimited(path, delimiter=","):
    """(field names, {name: column}) of a delimited text file with a header row.
    A column is numbers (float, NaN where empty) if every value is one, else text."""
    with open(path, "r", newline="") as f:
        reader = csv.reader(f, delimiter=delimiter)
        names = next(reader)
        rows = [r for r in reader if any(v.strip() for v in r)]
    columns = {}
    for i, name in enumerate(names):
        values = [r[i].strip() if i < len(r) else "" for r in rows]
        try:
            columns[name] = np.array([float(v) if v else np.nan for v in values])
        except ValueError:
            columns[name] = np.array(values, dtype=object)
    return names, columns


def read_header(path, delimiter=","):
    """The field names in the first row of a delimited text file (empty ones left out)."""
    with open(path, "r") as f:
        return list(filter(None, f.readline().strip().split(delimiter)))


# ---------------------------------------------------------------- Noddy grids
def read_noddy(path):
    """
    A Noddy ASCII grid (.mag or .grv).

    Returns a dict: values (2-D array, row 0 = north), gt (GDAL geotransform),
    data_type ("mag" for file code 333, else "grv") and metadata (the header
    values as strings, as the plugin writes them into the GeoTIFF).
    """
    with open(path, "r") as f:
        lines = f.readlines()

    # the file code is the first line that is not blank or a comment
    header_start = 0
    for i, line in enumerate(lines):
        if not line.strip().startswith("#") and line.strip():
            header_start = i
            break
    file_code = int(lines[header_start].strip())
    data_type = "mag" if file_code == 333 else "grv"

    at = header_start + 1
    while at < len(lines) and lines[at].strip().startswith("#"):
        at += 1
    calc_range, cols, rows, layers = map(int, lines[at].strip().split())
    inclination, declination, intensity = map(float, lines[at + 1].strip().split())
    sw_x, sw_y, sw_z = map(float, lines[at + 2].strip().split())
    ne_x, ne_y, ne_z = map(float, lines[at + 3].strip().split())
    cell_size, survey_height = map(float, lines[at + 4].strip().split())

    values = np.loadtxt(path, skiprows=8)
    n_rows, n_cols = values.shape
    gt = [sw_x, (ne_x - sw_x) / n_cols, 0, ne_y, 0, -(ne_y - sw_y) / n_rows]
    metadata = {
        "file_code": str(file_code),
        "calculation_range": str(calc_range),
        "original_rows": str(rows),
        "original_cols": str(cols),
        "layers": str(layers),
        "intensity": str(intensity),
        "inclination": str(inclination),
        "declination": str(declination),
        "sw_corner_x": str(sw_x),
        "sw_corner_y": str(sw_y),
        "sw_corner_z": str(sw_z),
        "ne_corner_x": str(ne_x),
        "ne_corner_y": str(ne_y),
        "ne_corner_z": str(ne_z),
        "cube_size": str(cell_size),
        "survey_height": str(survey_height),
    }
    return {"values": values, "gt": gt, "data_type": data_type, "metadata": metadata}


NODDY_EPSG = 32629  # the plugin gives Noddy grids UTM 29N


def noddy_output_name(path, data_type):
    """<file>_mag.tif or <file>_grv.tif beside a Noddy grid."""
    return f"{os.path.splitext(path)[0]}_{data_type}.tif"


# ---------------------------------------------------------------- profiles
def regularize_profile(coords, values, spacing="auto", num_points=None):
    """
    Points of a line resampled to an even spacing.

    coords: (n, 2) x, y; values: n values. The points are put in order along x or
    y, whichever spans more, and the values are interpolated linearly along the
    path through them. spacing: a distance, or "auto" for the median distance
    between the points; if neither, num_points (or the number of points).

    Returns (new coords (m, 2), new values (m,), spacing used).
    """
    from scipy.interpolate import interp1d

    coords = np.asarray(coords, dtype=float)
    values = np.asarray(values)
    x_range = np.max(coords[:, 0]) - np.min(coords[:, 0])
    y_range = np.max(coords[:, 1]) - np.min(coords[:, 1])
    if x_range >= y_range:
        logger.debug("Sorting points along X axis (dominant)")
        order = np.argsort(coords[:, 0])
    else:
        logger.debug("Sorting points along Y axis (dominant)")
        order = np.argsort(coords[:, 1])
    coords_sorted = coords[order]
    values_sorted = values[order]

    distances = [0]
    steps = []
    for i in range(1, len(coords_sorted)):
        prev, curr = coords_sorted[i - 1], coords_sorted[i]
        d = np.sqrt((curr[0] - prev[0]) ** 2 + (curr[1] - prev[1]) ** 2)
        steps.append(d)
        distances.append(distances[-1] + d)
    total_length = distances[-1]

    if spacing == "auto":
        used_spacing = np.median(steps)
        logger.debug("Using auto-calculated median spacing: %s", used_spacing)
        count = int(total_length / used_spacing) + 1
    elif isinstance(spacing, (int, float)):
        count = int(total_length / spacing) + 1
        used_spacing = spacing
    else:
        count = num_points if num_points is not None else len(coords_sorted)
        used_spacing = total_length / (count - 1) if count > 1 else 0

    new_distances = np.linspace(0, total_length, count)
    x_interp = interp1d(distances, coords_sorted[:, 0], kind="linear")
    y_interp = interp1d(distances, coords_sorted[:, 1], kind="linear")
    value_interp = interp1d(distances, values_sorted, kind="linear")
    new_coords = np.column_stack((x_interp(new_distances), y_interp(new_distances)))
    return new_coords, value_interp(new_distances), used_spacing


def points_along_line(vertices, spacing):
    """
    Points at 0, spacing, 2 x spacing ... along a polyline, up to its length (the
    end point only if the length is a whole number of spacings), as QGIS's
    Points Along Geometry does.

    Returns (distances (m,), points (m, 2)).
    """
    v = np.asarray(vertices, dtype=float)[:, :2]
    if spacing <= 0:
        raise ValueError("spacing must be greater than 0")
    seg = np.hypot(*np.diff(v, axis=0).T)
    along = np.concatenate([[0.0], np.cumsum(seg)])
    total = along[-1]
    distances = np.arange(int(np.floor(total / spacing * (1 + 1e-12))) + 1) * spacing
    x = np.interp(distances, along, v[:, 0])
    y = np.interp(distances, along, v[:, 1])
    return distances, np.column_stack([x, y])


def sample_cells(values, gt, points):
    """The value of the cell each point falls in (NaN outside the grid), as QGIS's
    raster sampling does. values: 2-D array, row 0 at gt's origin."""
    points = np.asarray(points, dtype=float)
    col = np.floor((points[:, 0] - gt[0]) / gt[1]).astype(int)
    row = np.floor((points[:, 1] - gt[3]) / gt[5]).astype(int)
    ny, nx = values.shape
    inside = (col >= 0) & (col < nx) & (row >= 0) & (row < ny)
    out = np.full(len(points), np.nan)
    out[inside] = values[row[inside], col[inside]]
    return out
