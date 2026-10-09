"""
Replay a processing history on another grid.

The provenance SGTool embeds in every grid it saves lists the operation and
parameters of each step. plan_steps() turns that history into a recipe, and
run_replay() applies the recipe to a different grid, step by step, saving every
intermediate grid next to the input with the usual file names
(grid.tif -> grid_RTP.tif -> grid_RTP_d1z.tif) and its own provenance, so the
new grid can itself be replayed or inspected.

Magnetic reductions (reduction to the pole or equator, differential RTP) depend
on the field direction at the survey, so their inclination and declination are
recalculated from the IGRF model for the new grid's location and date rather
than copied from the original.

run_replay() only computes and writes files (no QGIS layers, widgets or
messages), so it can run in a background task or a Processing algorithm.
"""

import os
from datetime import date as _date

import numpy as np

from .calcs.layer_shim import LayerShim
from .calcs.sgt_metadata import write_sgt_metadata
from .calcs.Cooper_Cowan_rtpvariable import rtpvariable as cooper_rtpvariable
from .calcs.igrf.igrf_field import grid_field, decimal_year
from .sgt_processing import (
    FILTER_SPECS, SUFFIXES, GridContext, auto_output_path, read_grid, write_grid,
    _plugin_version,
)

DRTP_OPERATION = "Differential reduction to pole"
LINE_NOISE_OPERATION = "Directional Cosine/Butterworth line-noise removal"

_SPEC_BY_ID = {s["id"]: s for s in FILTER_SPECS}

# the dialog's wording, and the algorithms' (which is the spec's own operation)
_OPERATION_TO_SPEC = {s["operation"]: s["id"] for s in FILTER_SPECS}
_OPERATION_TO_SPEC.update({
    "Upward continuation": "continuation",
    "Downward continuation": "continuation",
    "Low pass": "high_low_pass",
    "High pass": "high_low_pass",
})
# settings implied by the operation's name
_FIXED = {
    "Upward continuation": {"DIRECTION": "up"},
    "Downward continuation": {"DIRECTION": "down"},
    "Low pass": {"TYPE": "Low"},
    "High pass": {"TYPE": "High"},
}
# the dialog records some settings under other names than the algorithms
_RENAMES = {
    "direction_width_deg": "WEDGE",
    "transition_width": "WIDTH",
    "cutoff_wavelength": "CUTOFF",
    "window_pixels": "WINDOW",
    "window_size": "WINDOW",
    "filter_size": "SIZE",
    "polynomial_order": "ORDER",
    "relief_shading": "RELIEF",
}
_IGRF_SPECS = ("reduction_to_pole", "reduction_to_equator")


class ReplayError(Exception):
    """A problem that stops the replay, with a message for the user."""


# ---------------------------------------------------------------- planning
def _convert(spec_param, text):
    """The recorded text as the parameter's value; ValueError if it won't fit."""
    kind = spec_param["kind"]
    text = str(text).strip()
    if kind == "double":
        return float(text)
    if kind == "int":
        return int(round(float(text)))
    if kind == "bool":
        if text.lower() in ("true", "1", "yes"):
            return True
        if text.lower() in ("false", "0", "no"):
            return False
        raise ValueError(text)
    options = spec_param["options"]  # enum: by its text, or by its start ("1" -> "1st order")
    for option in options:
        if option.lower() == text.lower():
            return option
    for option in options:
        if option.lower().startswith(text.lower()) and text:
            return option
    raise ValueError(text)


def _settings_text(parameters):
    return ", ".join(f"{k} {v}" for k, v in parameters.items() if v != "")


def _plan_one(index, step):
    operation = step["operation"]
    item = {
        "index": index,
        "operation": operation,
        "kind": "unsupported",
        "spec_id": None,
        "values": {},
        "igrf": False,
        "line_noise": False,
        "settings": _settings_text(step["parameters"]),
        "note": "",
        "replayable": False,
        "date": None,
    }
    if step["n_sources"] > 1:
        item["note"] = "made from several grids, so it cannot be replayed on one"
        return item

    if operation == DRTP_OPERATION:
        item.update(kind="drtp", igrf=True, replayable=True)
        recorded = step["parameters"].get("date", "")
        try:
            y, m, d = (int(p) for p in recorded.split("-"))
            item["date"] = _date(y, m, d)
        except ValueError:
            pass
        item["note"] = "inclination and declination recalculated from IGRF"
        return item

    spec_id = _OPERATION_TO_SPEC.get(operation)
    if spec_id is None:
        if step["parameters"] or operation:
            item["note"] = (
                "not a step that can be replayed (for example gridding, "
                "Euler deconvolution or component analysis)"
            )
        return item

    spec = _SPEC_BY_ID[spec_id]
    recorded = {}
    for name, text in step["parameters"].items():
        key = _RENAMES.get(name.lower(), name.upper())
        recorded[key] = text
    recorded.update(_FIXED.get(operation, {}))

    values, defaults = {}, []
    for p in spec["params"]:
        name = p["name"]
        if name in recorded:
            try:
                values[name] = _convert(p, recorded[name])
                continue
            except ValueError:
                pass
        values[name] = p["default"] if p["kind"] != "enum" else p["options"][p["default"]]
        defaults.append(name)
    notes = []
    if defaults:
        notes.append("not recorded, default used for " + ", ".join(defaults))
    if spec_id in _IGRF_SPECS:
        notes.append("inclination and declination recalculated from IGRF")
    item.update(
        kind="filter", spec_id=spec_id, values=values, replayable=True,
        igrf=spec_id in _IGRF_SPECS,
        line_noise=spec_id == "directional_line_noise",
        note="; ".join(notes),
    )
    return item


def plan_steps(steps):
    """A recipe item for every step of a history (see history_steps)."""
    return [_plan_one(i, s) for i, s in enumerate(steps)]


def step_label(item):
    return item["operation"] or "Unnamed step"


# ---------------------------------------------------------------- naming
def _suffix(item, values):
    if item["kind"] == "drtp":
        return "_DRTP"
    return SUFFIXES[item["spec_id"]](values)


def output_paths(source_path, items):
    """Every file the replay of these recipe items will write, in order."""
    paths, current = [], source_path
    for item in items:
        out = auto_output_path(current, _suffix(item, item["values"]))
        paths.append(out)
        if item["kind"] == "drtp":
            paths.append(auto_output_path(current, "_DRTP_inc"))
            paths.append(auto_output_path(current, "_DRTP_dec"))
        current = out
    return paths


def final_output_path(source_path, items):
    current = source_path
    for item in items:
        current = auto_output_path(current, _suffix(item, item["values"]))
    return current


def clear_outputs(paths):
    """Delete earlier results (and their side files) about to be replaced."""
    for path in paths:
        for stale in (path, path + ".aux.xml", path + ".sgt.xml"):
            if os.path.exists(stale):
                try:
                    os.remove(stale)
                except OSError:
                    raise ReplayError(
                        f"{stale} already exists and could not be replaced. If it is "
                        "open in QGIS, remove that layer first."
                    )


# ---------------------------------------------------------------- IGRF
def uniform_field(path):
    """(inclination, declination) stored in a GeoTIFF's metadata (as the Noddy
    import does for a synthetic model), or None."""
    try:
        from osgeo import gdal

        ds = gdal.Open(str(path))
        md = ds.GetMetadata() if ds is not None else {}
        ds = None
        if all(k in md for k in ("inclination", "declination", "intensity")):
            return float(md["inclination"]), float(md["declination"])
    except Exception:
        pass
    return None


def check_igrf_possible(grid_info):
    """None if the field direction can be worked out for this grid, else why not."""
    authid = grid_info.get("authid") or ""
    if not authid.upper().startswith("EPSG:"):
        return (
            "The new grid needs a coordinate system with an EPSG code so that the "
            "IGRF field direction can be calculated for its location"
        )
    return None


def _field_for_grid(source_path, grid_info, when):
    """(inc, dec) at the grid centre, and for differential RTP the corner values."""
    fixed = uniform_field(source_path)
    if fixed is not None:
        return fixed, ((fixed[0],) * 4, (fixed[1],) * 4, fixed[0], fixed[1])
    year = decimal_year(when.year, when.month, when.day)
    inc_c, dec_c = grid_field(grid_info["extent"], grid_info["authid"], year)
    corners = grid_field(grid_info["extent"], grid_info["authid"], year, corners=True)
    return (inc_c, dec_c), corners


# ---------------------------------------------------------------- replay
def run_replay(
    source_path, items, grid_info, when=None,
    task=None,
):
    """Apply recipe items (from plan_steps) to the grid at source_path.

    grid_info: facts about the grid read on the main thread (extent, authid,
    geographic). when: date of the survey (datetime.date) for the IGRF.
    task: optional object with check(percent) to report progress / cancel.

    An item with "save" set to False is an intermediate result that is deleted
    once the next step has been made from it (the last result is always kept;
    the provenance of a deleted step stays inside the files made from it).

    Returns {"outputs": [paths of the grids kept, in order],
             "final": path of the last one, "n_steps": steps applied,
             "log": [messages]}.
    """
    when = when or _date.today()
    log = []
    shim = LayerShim(grid_info["extent"], grid_info["authid"], grid_info.get("geographic", False))
    current = source_path
    outputs = []
    discard = []  # unsaved intermediate waiting for the next step to read it
    min_max = []  # results to display stretched min to max (sun shading)
    n = max(len(items), 1)

    def check(done):
        if task is not None:
            task.check(100.0 * done / n)

    for k, item in enumerate(items):
        check(k)
        arr, gt, projection = read_grid(current)
        values = dict(item["values"])
        record = None
        grid = GridContext(arr, gt, shim, 0)

        if item["kind"] == "drtp":
            record = _replay_drtp(item, current, source_path, grid_info, when,
                                  arr, gt, projection, log)
            result = record.pop("grid")
            suffix = "_DRTP"
        else:
            spec = _SPEC_BY_ID[item["spec_id"]]
            if item["igrf"]:
                (inc, dec), _corners = _field_for_grid(source_path, grid_info, when)
                inc, dec = float(inc), float(dec)
                values["INCLINATION"], values["DECLINATION"] = inc, dec
                record = dict(values, igrf_date=when.isoformat())
                log.append(
                    f"{item['operation']}: inclination {inc:.1f}, declination {dec:.1f} "
                    f"from IGRF for {when.isoformat()}"
                )
            result = spec["compute"](grid, values)
            if result is None:
                raise ReplayError(f"{item['operation']}: nothing was calculated")
            suffix = SUFFIXES[item["spec_id"]](values)
            record = record or values

        if task is not None:
            task.check()
        out = auto_output_path(current, suffix)
        clear_outputs([out])
        write_grid(out, np.asarray(result), gt, projection)
        write_sgt_metadata(out, current, item["operation"], record, _plugin_version())
        # the step before is no longer needed once this one has been made from it
        _delete_files(discard)
        discard = []
        if item.get("spec_id") == "sun_shading":
            min_max.append(out)  # shown min to max rather than mean +/- 2 std
        last = k == len(items) - 1
        extras = (
            [auto_output_path(current, "_DRTP_inc"), auto_output_path(current, "_DRTP_dec")]
            if item["kind"] == "drtp" else []
        )
        if item.get("save", True) or last:
            outputs.append(out)
        else:
            discard = [out] + extras
        current = out
    check(n)
    return {"outputs": outputs, "final": current, "n_steps": len(items), "log": log,
            "min_max": min_max}


def _delete_files(paths):
    """Delete grids (and their side files) that were only needed on the way."""
    for path in paths:
        for f in (path, path + ".aux.xml", path + ".sgt.xml"):
            try:
                if os.path.exists(f):
                    os.remove(f)
            except OSError:
                pass  # leave it: a stray intermediate is harmless


def _replay_drtp(item, current, source_path, grid_info, when,
                 arr, gt, projection, log):
    """Differential RTP for the grid, plus the inclination / declination grids
    it used (written, as the dialog does, but not loaded)."""
    _centre, corners = _field_for_grid(source_path, grid_info, when)
    inc_corners, dec_corners, inc_c, dec_c = corners
    nr, nc = arr.shape
    _, varrtp, incv, decv = cooper_rtpvariable(
        np.flipud(arr), nr, nc, inc_corners, dec_corners,
        inc_center=inc_c, dec_center=dec_c,
    )
    log.append(
        f"{item['operation']}: field from IGRF for {when.isoformat()} "
        f"(centre inclination {inc_c:.1f}, declination {dec_c:.1f})"
    )
    for field, tag, what in ((incv, "_DRTP_inc", "inclination"), (decv, "_DRTP_dec", "declination")):
        path = auto_output_path(current, tag)
        clear_outputs([path])
        write_grid(path, np.flipud(field), gt, projection)
        write_sgt_metadata(
            path, current, f"Differential RTP: interpolated {what} field",
            {"date": when.isoformat()}, _plugin_version(),
        )
    return {"grid": np.flipud(varrtp), "date": when.isoformat()}
