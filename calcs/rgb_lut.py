"""
RGB image to a grey (scalar) grid through a colour look-up table.

Shared by the QGIS plugin, the ArcGIS Pro GUI and pysgtool. Needs matplotlib
(for the colour names and the colour ramp) and scipy.
"""

import re

import numpy as np


def parse_lut_string(LUT, num_entries=1024):
    """
    The look-up table of a colour list, from the lowest value to the highest.

    LUT (str): CSS colour names or RGB triplets, from the highest value to the
        lowest, e.g. "red,green,blue", "(255,0,0),(0,255,0),(0,0,255)",
        "255,0,0,0,255,0,0,0,255" or "(1.0,0,0),(0,1.0,0),(0,0,1.0)".

    Returns a list of [value 0..1, (R, G, B) 0..255], or False if the colours
    cannot be read. Raises ImportError without matplotlib.
    """
    LUT = LUT.replace(" ", "")

    # parentheses: (R,G,B),(R,G,B)
    if "(" in LUT:
        matches = re.findall(r"\(([^)]+)\)", LUT)
        if matches:
            rgb_list = []
            for match in matches:
                values = match.split(",")
                try:
                    if len(values) >= 3:
                        rgb_list.append(tuple(float(v) for v in values[:3]))
                    else:
                        raise ValueError("Not enough values for RGB")
                except ValueError:
                    break  # not numbers: read them as colour names below
            else:
                rgb_list.reverse()
                return lut_from_rgb(rgb_list, num_entries)

    elements = LUT.split(",")

    # a flat list of numbers, three to a colour
    try:
        numbers = [float(elem) for elem in elements]
        if len(numbers) >= 3 and len(numbers) % 3 == 0:
            rgb_list = [tuple(numbers[i:i + 3]) for i in range(0, len(numbers), 3)]
            rgb_list.reverse()
            return lut_from_rgb(rgb_list, num_entries)
    except ValueError:
        pass

    # colour names
    css_color_list = elements
    css_color_list.reverse()
    return lut_from_names(css_color_list, num_entries)


def lut_from_rgb(rgb_list, num_entries=1024):
    """Look-up table of RGB triplets (0..255, or 0..1 if no value is above 1):
    a list of [value 0..1, (R, G, B) 0..255], or False."""
    normalized = []
    for rgb in rgb_list:
        r, g, b = rgb[0], rgb[1], rgb[2]
        if any(val > 1 for val in (r, g, b)):
            normalized.append((r / 255.0, g / 255.0, b / 255.0))
        else:
            normalized.append((r, g, b))
    return _ramp(normalized, num_entries)


def lut_from_names(css_color_list, num_entries=1024):
    """Look-up table of colour names matplotlib knows: a list of
    [value 0..1, (R, G, B) 0..255], or False."""
    return _ramp(css_color_list, num_entries)


def _ramp(colours, num_entries):
    import matplotlib.colors as mcolors

    decimal_indices = np.linspace(0, 1, num_entries)
    try:
        cmap = mcolors.LinearSegmentedColormap.from_list("custom_cmap", colours)
    except Exception:
        return False
    rgb_colors = [cmap(i)[:3] for i in decimal_indices]
    rgb_colors_255 = [(int(r * 255), int(g * 255), int(b * 255)) for r, g, b in rgb_colors]
    return [[round(decimal_index, 6), rgb]
            for decimal_index, rgb in zip(decimal_indices, rgb_colors_255)]


def rgb_to_scalar(red, green, blue, lut, lut_min, lut_max):
    """
    The value of each cell of an RGB image: that of the nearest colour of the
    look-up table (parse_lut_string), scaled from 0..1 to lut_min..lut_max.
    White (255, 255, 255) and black (0, 0, 0) cells are NaN.
    """
    from scipy.spatial import cKDTree

    rgb_raster = np.dstack((np.asarray(red, dtype=float), np.asarray(green, dtype=float),
                            np.asarray(blue, dtype=float)))
    scalar_values, lut_colors = zip(*lut)
    lut_colors = np.array(lut_colors) / 255.0

    white_mask = (rgb_raster == [255, 255, 255]).all(axis=2)
    black_mask = (rgb_raster == [0, 0, 0]).all(axis=2)

    _distances, indices = cKDTree(lut_colors).query((rgb_raster / 255.0).reshape(-1, 3))
    scalar_grid = np.array(scalar_values)[indices].reshape(rgb_raster.shape[:2])
    scalar_grid[white_mask] = np.nan
    scalar_grid[black_mask] = np.nan
    return (scalar_grid * (lut_max - lut_min)) + lut_min
