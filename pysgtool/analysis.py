# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2024-2026 Mark Jessell
"""
Calculations on several grids, on profiles, and the processing history: radiometric
ratios, IGRF values, the 1-D wavelet transform (WTMM) and replay.
"""

import os

import numpy as np

from .calcs import Radiometrics as _R
from .calcs import sgt_metadata
from .calcs import replay as _replay
from .calcs.workflow import history_steps
from .grid import Grid, _decimal, _Raster, as_date


# ---------------------------------------------------------------- radiometrics
def _radiometric(kind, grids, record, calculate):
    first = next(iter(grids.values()))
    for name, g in grids.items():
        if not isinstance(g, Grid):
            raise TypeError(f"{name} must be a Grid")
    stack = np.stack([g.values if g is first else g.resample_to(first).values
                      for g in grids.values()])
    operation, suffix = _R.RADIOMETRIC_OPERATIONS[kind]
    record = dict({name: g.path or g.name for name, g in grids.items()}, **record)
    return first._derive(calculate(stack, first), operation, record, suffix)


def efimov(k, u, th):
    """Efimov F-parameter K * eU / eTh, on the cells of the K grid."""
    return _radiometric("efimov", {"K": k, "eU": u, "eTh": th}, {},
                        lambda s, g: _R.efimov_f(s[0], s[1], s[2]))


def kd(k, th, background="ratio"):
    """Thorium-normalised deviation of K. background: predicted from eTh by the 'ratio' of
    the means or by 'regression'."""
    return _deviation("kd", "K", k, th, background)


def ud(u, th, background="ratio"):
    """Thorium-normalised deviation of eU (see kd)."""
    return _deviation("ud", "eU", u, th, background)


def _deviation(kind, name, x, th, background):
    if background not in ("ratio", "regression"):
        raise ValueError("background is 'ratio' or 'regression'")
    return _radiometric(kind, {name: x, "eTh": th}, {"background": background},
                        lambda s, g: _R.thorium_normalised(s[0], s[1], background))


def ternary_gradient(k, u, th, scaling="zscore"):
    """Multichannel (Di Zenzo) Sobel gradient of K, eU and eTh. scaling: 'zscore' weighs
    the three equally, 'clr' follows changes of their ratios only."""
    if scaling not in ("zscore", "clr"):
        raise ValueError("scaling is 'zscore' or 'clr'")
    return _radiometric("ternary_gradient", {"K": k, "eU": u, "eTh": th}, {"scaling": scaling},
                        lambda s, g: _R.ternary_gradient(s, scaling, dx=g.dx, dy=g.dy))


# ---------------------------------------------------------------- IGRF
def igrf(lat, lon, date, altitude_km=0.0):
    """(inclination, declination, intensity in nT) of the IGRF at geodetic positions
    (numbers or arrays) on a date ('yyyy-mm-dd', a date or a decimal year)."""
    from .calcs.igrf.igrf_field import calc_igrf

    return calc_igrf(_decimal(date), altitude_km, lat, lon)


# ---------------------------------------------------------------- WTMM
def wtmm(profile, wavelet="mexh", scales=None, num_scales=15, threshold_rel=0.05,
         min_distance=3, chaining_threshold=2.0):
    """1-D wavelet transform modulus maxima of an evenly spaced profile (needs PyWavelets),
    with the plugin's settings. Returns a dict: scales, cwt, modulus, maxima, maxima_lines."""
    from .calcs.WTMM import WTMM

    profile = np.asarray(profile, dtype=float).ravel()
    if not np.isfinite(profile).all():
        raise ValueError("The profile has gaps (NaN): fill them first")
    if np.ptp(profile) == 0:
        raise ValueError("The profile is constant")
    return WTMM().wtmm_1d(profile, wavelet, scales, num_scales, threshold_rel, min_distance,
                          chaining_threshold)


def plot_wtmm(profile, results, name="profile", line_number=None, save_path=None):
    """Figure of the profile, its scalogram, the maxima lines and the cladogram."""
    from .calcs.WTMM import WTMM

    return WTMM().visualize_wtmm_1d(np.asarray(profile, dtype=float).ravel(), name, results,
                                    line_number=line_number, save_path=save_path)


def multifractal_spectrum(results, q_values=None):
    """(h, D(h)) of the results of wtmm()."""
    from .calcs.WTMM import WTMM

    return WTMM().simplified_multifractal_spectrum(results, q_values)


# ---------------------------------------------------------------- history
def _provenance(source):
    if isinstance(source, _Raster):
        root, what = source.provenance, source.name
    else:
        path = str(source)
        if not os.path.exists(path):
            raise FileNotFoundError(path)
        root = (sgt_metadata.parse_xml_file(path) if path.lower().endswith(".xml")
                else sgt_metadata.read_sgt_metadata(path))
        what = path
    if root is None:
        raise ValueError(f"No SGTool processing history found in {what}")
    return root


def history(source):
    """The processing steps of a Grid, a file SGTool saved, or an exported .sgt.xml,
    oldest first: dicts with operation, parameters, created, version ..."""
    return history_steps(_provenance(source))[1]


def replay(source, target, date=None, save_intermediate=True):
    """Apply the steps of a history to another grid file, saving the results beside it
    with SGTool's usual names, as the plugin's "Apply the same steps to another grid".

    source: a Grid, a file SGTool saved or an exported .sgt.xml. target: the path of a
    grid. date: survey date for the IGRF (by default the one in the history, else today).
    Returns the dict of calcs.replay.run_replay: outputs, final, log ..., with
    "skipped", the steps that cannot be repeated."""
    from .io import read

    plan = _replay.plan_steps(history(source))
    items = [p for p in plan if p["replayable"]]
    if not items:
        raise ValueError("The history has no steps that can be repeated")
    for item in items:
        item["save"] = save_intermediate
    if date is not None:
        when = as_date(date)
    else:
        when = next((p["date"] for p in items if p["date"] is not None), None)
    grid = read(target)
    code = grid.epsg
    info = {"extent": grid.extent, "authid": f"EPSG:{code}" if code else "",
            "geographic": grid.geographic}
    if any(p["igrf"] for p in items) and grid._own_field() is None:
        why = _replay.check_igrf_possible(info)
        if why:
            raise ValueError(why)
    result = _replay.run_replay(grid.path, items, info, when=when)
    result["skipped"] = [f"{p['operation']}: {p['note']}" for p in plan if not p["replayable"]]
    return result
