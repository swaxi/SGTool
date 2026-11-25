"""
The filters of calcs/filter_table.py as methods of Grid, made from that table: a filter
added to it is here too, with the same name, settings and result as in QGIS.
"""

import inspect
import math

import numpy as np

from .calcs.filter_table import FILTER_SPECS, SUFFIXES
from .calcs.replay import _IGRF_SPECS, _convert
from .grid import Grid, as_date

# settings whose default in the dialog is only a placeholder: here they are given, or
# worked out from the IGRF for a date
_FROM_FIELD = {spec_id: ("INCLINATION", "DECLINATION") for spec_id in _IGRF_SPECS}
_FROM_FIELD["apparent_susceptibility"] = ("FIELD",)


def _default(p):
    return p["options"][p["default"]] if p["kind"] == "enum" else p["default"]


def _value(p, value):
    """The value of a setting, checked against the table."""
    name = p["name"].lower()
    kind = p["kind"]
    if kind == "bool":
        if not isinstance(value, (bool, np.bool_)):
            raise ValueError(f"{name} is True or False")
        return bool(value)
    if kind == "enum":
        try:
            return _convert(p, value)
        except ValueError:
            raise ValueError(f"{name} is one of {', '.join(p['options'])}, not {value!r}") from None
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be a number, not {value!r}") from None
    if not math.isfinite(number) or (kind == "int" and number != int(number)):
        raise ValueError(f"{name} must be a {'whole' if kind == 'int' else 'finite'} number")
    for limit, bad, text in (("minimum", lambda m: number < m, "at least"),
                             ("maximum", lambda m: number > m, "at most")):
        if p.get(limit) is not None and bad(p[limit]):
            raise ValueError(f"{name} must be {text} {p[limit]:g}")
    return int(number) if kind == "int" else number


def _fill_from_field(grid, spec, values, date):
    names = _FROM_FIELD[spec["id"]]
    missing = [n for n in names if values[n] is None]
    if not missing:
        return {}
    if date is None:
        raise ValueError(
            f"{spec['id']} needs {' and '.join(n.lower() for n in missing)}, or the survey "
            "date to work it out from the IGRF")
    inc, dec, intensity = grid.field_direction(date)
    values.update({n: {"INCLINATION": inc, "DECLINATION": dec, "FIELD": intensity}[n]
                   for n in missing})
    return {"igrf_date": as_date(date).isoformat()}


def _method(spec):
    from_field = _FROM_FIELD.get(spec["id"], ())
    P = inspect.Parameter
    parameters = [
        P(p["name"].lower(), P.POSITIONAL_OR_KEYWORD,
          default=None if p["name"] in from_field else _default(p))
        for p in spec["params"]
    ]
    if from_field:
        parameters.append(P("date", P.KEYWORD_ONLY, default=None))
    if spec["fft"]:
        parameters.append(P("buffer", P.KEYWORD_ONLY, default=0))
    signature = inspect.Signature([P("self", P.POSITIONAL_OR_KEYWORD)] + parameters)

    def method(self, *args, **kwargs):
        bound = signature.bind(self, *args, **kwargs)
        bound.apply_defaults()
        given = bound.arguments
        values = {
            p["name"]: None if p["name"] in from_field and given[p["name"].lower()] is None
            else _value(p, given[p["name"].lower()])
            for p in spec["params"]
        }
        extra = _fill_from_field(self, spec, values, given["date"]) if from_field else {}
        result = spec["compute"](self._context(given.get("buffer", 0)), values)
        if result is None:
            raise ValueError("Nothing was calculated for these settings")
        return self._derive(np.asarray(result), spec["operation"], dict(values, **extra),
                            SUFFIXES[spec["id"]](values))

    lines = [spec["help"], ""]
    for p in spec["params"]:
        text = f"{p['name'].lower()}: {p['label']}"
        if p["kind"] == "enum":
            text += ": " + ", ".join(repr(o) for o in p["options"])
        lines.append(text)
    if from_field:
        lines.append("date: survey date, to take what is not given above from the IGRF")
    if spec["fft"]:
        lines.append("buffer: padding in cells for the Fourier transform (0 = automatic)")
    method.__name__ = spec["id"]
    method.__qualname__ = f"Grid.{spec['id']}"
    method.__doc__ = "\n".join(lines)
    method.__signature__ = signature
    return method


for _spec in FILTER_SPECS:
    if hasattr(Grid, _spec["id"]):
        raise RuntimeError(f"The filter {_spec['id']} would hide a method of Grid")
    setattr(Grid, _spec["id"], _method(_spec))


def filters():
    """{name: {"group", "help", "settings": {name: default}}} of every filter method."""
    out = {}
    for spec in FILTER_SPECS:
        from_field = _FROM_FIELD.get(spec["id"], ())
        out[spec["id"]] = {
            "group": spec["group"][1],
            "help": spec["help"],
            "settings": {p["name"].lower(): None if p["name"] in from_field else _default(p)
                         for p in spec["params"]},
        }
    return out
