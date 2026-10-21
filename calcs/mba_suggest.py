"""
Suggested Multilevel B-Spline gridding settings from the points themselves.

Pure NumPy / SciPy (no QGIS), so the QGIS dialog, the Processing algorithm and the ArcGIS Pro
GUI all use the same rule.

Cell size
---------
What limits the detail a set of points can support is not the along-line sampling but the gaps
between the points. The median distance from a random location inside the survey footprint to
the nearest data point measures that gap directly:

* for flight lines it is a quarter of the line spacing (a location is uniformly 0 to half a line
  spacing from the nearest line), which is the usual rule of thumb for a gridding cell size;
* for scattered points it is about half the mean point spacing.

The suggested cell size is that gap distance for line-like data (the points are much closer along
the lines than between them) and twice it for scattered data. It is rounded to two significant
figures, which is well inside the precision of the rule.

Max Levels
----------
Level L of the multilevel B-spline has a control lattice cell of range / 2**(L - 1), range being
the longer side of the grid. Levels with a lattice finer than half the output cell cannot add
anything the grid can show, and on noisy data they start to fit the noise (a held-out test on 824k
aeromagnetic points gave the lowest error at lattice = cell / 2 and a slightly worse one beyond
it), so the suggestion is lattice = cell / 2:  L = round(2 + log2(range / cell)).
(The Threshold Error does not help choose the level: the loop only stops early when no point is
above it, which noisy data never allows.)
"""

from __future__ import annotations

import math

import numpy as np

METRES_PER_DEGREE = 111320.0


def round_sig(value: float, digits: int = 2) -> float:
    """value rounded to a number of significant figures (0.001134 -> 0.0011)."""
    if not np.isfinite(value) or value == 0:
        return float(value)
    return float(f"{value:.{digits}g}")


def levels_for_cell(extent_range: float, cell_size: float) -> int:
    """Max Levels for which the finest lattice is half the output cell (3 to 20)."""
    if cell_size <= 0 or extent_range <= 0:
        return 3
    return int(min(20, max(3, round(2 + math.log2(extent_range / cell_size)))))


def _hull_mask(points: np.ndarray, locations: np.ndarray) -> np.ndarray:
    """Which of the locations lie inside the convex hull of the points."""
    try:
        from scipy.spatial import ConvexHull

        hull = ConvexHull(points)
        normals, offsets = hull.equations[:, :2], hull.equations[:, 2]
        return np.all(locations @ normals.T + offsets <= 1e-12, axis=1)
    except Exception:  # fewer than three non-collinear points: use the bounding box
        return np.ones(len(locations), dtype=bool)


def suggest_mba_parameters(x, y, geographic: bool = False, cell_size: float | None = None,
                           n_locations: int = 20000, seed: int = 0,
                           max_points: int = 3_000_000) -> dict:
    """Settings for gridding the points (x, y).

    geographic: the coordinates are degrees. East-west distances are then shortened by
        cos(latitude) so the gap is measured on the ground; the cell size stays in degrees.
    cell_size: if given it is used as it is and only Max Levels (and the warnings) follow from it.

    Returns a dict: n_points, x_range, y_range, extent_range, nn_spacing (median distance to
    the nearest other point), gap_distance (median distance from a random place in the footprint
    to the nearest point), line_like, cell_size_raw, cell_size, max_levels, n_cols, n_rows,
    units ("degrees" or "map units"), cell_ground_m (geographic only) and warnings (list of text).
    """
    from scipy.spatial import cKDTree

    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    keep = np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    if x.size < 10:
        raise ValueError("need at least 10 points to suggest settings")
    rng = np.random.default_rng(seed)
    if x.size > max_points:
        pick = rng.choice(x.size, max_points, replace=False)
        x, y = x[pick], y[pick]

    scale = math.cos(math.radians(float(np.mean(y)))) if geographic else 1.0
    pts = np.column_stack([x * scale, y])
    tree = cKDTree(pts)

    sample = rng.choice(len(pts), min(len(pts), 50000), replace=False)
    nn = tree.query(pts[sample], k=2)[0][:, 1]
    nn = nn[nn > 0]  # duplicated points have no spacing
    nn_spacing = float(np.median(nn)) if nn.size else 0.0

    lo, hi = pts.min(axis=0), pts.max(axis=0)
    locations = np.empty((0, 2))
    for _ in range(8):
        cand = rng.uniform(lo, hi, size=(4 * n_locations, 2))
        cand = cand[_hull_mask(pts[rng.choice(len(pts), min(len(pts), 20000), replace=False)], cand)]
        locations = np.vstack([locations, cand])
        if len(locations) >= n_locations:
            break
    locations = locations[:n_locations]
    gap = tree.query(locations)[0]
    gap_distance = float(np.median(gap))
    trimmed = gap[gap < 5.0 * gap_distance] if gap_distance > 0 else gap
    if trimmed.size:  # a survey outline with big voids in its hull must not inflate the gap
        gap_distance = float(np.median(trimmed))

    line_like = nn_spacing > 0 and gap_distance / nn_spacing >= 2.0
    cell_raw = gap_distance * (1.0 if line_like else 2.0)
    warnings = []
    if cell_raw <= 0:
        raise ValueError("the points do not have a usable spacing")

    x_range, y_range = float(x.max() - x.min()), float(y.max() - y.min())
    extent_range = max(x_range, y_range)
    suggested = round_sig(cell_raw, 2)
    used = float(cell_size) if cell_size else suggested
    if cell_size:
        ratio = used / suggested
        if ratio > 2.0:
            warnings.append(f"the cell size is {ratio:.1f} times the suggested one: it is coarser "
                            "than the point spacing supports and will smooth detail away")
        elif ratio < 0.4:
            warnings.append(f"the cell size is {1 / ratio:.1f} times smaller than the suggested "
                            "one: most cells have no data nearby and the result is interpolated")
    n_cols = 2 + int(x_range / used)
    n_rows = 2 + int(y_range / used)
    if n_cols * n_rows > 5e7:
        warnings.append(f"the grid would have {n_cols * n_rows / 1e6:.0f} million cells")
    if geographic and abs(float(np.mean(y))) > 60:
        warnings.append("at this latitude cells are much narrower east-west than north-south "
                        "on the ground")

    return {
        "n_points": int(len(pts)),
        "x_range": x_range, "y_range": y_range, "extent_range": extent_range,
        "nn_spacing": nn_spacing, "gap_distance": gap_distance, "line_like": bool(line_like),
        "cell_size_raw": cell_raw, "cell_size": suggested,
        "max_levels": levels_for_cell(extent_range, used),
        "n_cols": n_cols, "n_rows": n_rows,
        "units": "degrees" if geographic else "map units",
        "cell_ground_m": used * METRES_PER_DEGREE if geographic else None,
        "warnings": warnings,
    }


def describe(suggestion: dict) -> str:
    """One short paragraph for a dialog label."""
    unit = "°" if suggestion["units"] == "degrees" else ""
    kind = "line data" if suggestion["line_like"] else "scattered points"
    text = (f"{suggestion['n_points']:,} points ({kind}); distance from a place in the survey to "
            f"the nearest point {suggestion['gap_distance']:.3g}{unit}.")
    if suggestion["cell_ground_m"]:
        text += f" Cell is about {suggestion['cell_ground_m']:.0f} m north-south."
    for w in suggestion["warnings"]:
        text += " Note: " + w + "."
    return text
