# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 Mark Jessell
"""
Padding a grid for Fourier filtering, as Wormer.buildPaddedRaster does it: the
grid is centred in a larger array of zeros, and a taper (apodisation mask) takes
it smoothly to zero at its edges.

Written for SGTool to replace the libtim functions of the same names, which were
under a licence (CC BY-SA 3.0) that cannot be combined with the GPL. Only what
BSDWorms uses is provided.
"""

import numpy as np

# taper shapes: each maps 0 (outer edge of the taper) .. 1 (inner edge) to a weight
TAPERS = {
    "hann": lambda t: 0.5 - 0.5 * np.cos(np.pi * t),
    "hamm": lambda t: 0.54 - 0.46 * np.cos(np.pi * t),
    "cos": lambda t: np.sin(0.5 * np.pi * t),
    "sin": lambda t: np.sin(0.5 * np.pi * t),
    "lanc": lambda t: np.sinc(t - 1.0),
}


def _taper_function(apod_f):
    if callable(apod_f):
        return apod_f
    if isinstance(apod_f, str):
        name = apod_f.lower()
        for key, f in TAPERS.items():
            if name.startswith(key):
                return f
    raise ValueError(f"unknown taper {apod_f!r}: use hann, hamming, cosine, lanczos or a function")


def _per_axis(value, ndim, name):
    if np.ndim(value) == 0:
        return (value,) * ndim
    if len(value) != ndim:
        raise ValueError(f"{name} needs one value, or one per dimension")
    return tuple(value)


def mk_apod_mask(masksz, apodsz=None, wsize=0, apod_f="hann"):
    """
    A rectangular taper of shape masksz (cells), centred in it: 1 in the middle,
    falling to 0 over wsize cells at the edges of a box apodsz cells across, and
    0 outside that box. apodsz and wsize are in cells (one value, or one per
    axis); apodsz defaults to masksz. apod_f is the shape of the fall: "hann"
    (the default), "hamming", "cosine", "lanczos" or a function of t in [0, 1].
    The taper in 2-D is the product of the taper along each axis.
    """
    masksz = tuple(int(n) for n in masksz)
    if min(masksz) < 1:
        raise ValueError("every dimension of masksz must be at least 1")
    ndim = len(masksz)
    boxes = _per_axis(masksz if apodsz is None else apodsz, ndim, "apodsz")
    widths = _per_axis(wsize, ndim, "wsize")
    taper = _taper_function(apod_f)

    mask = np.ones(masksz)
    for axis, (n, box, width) in enumerate(zip(masksz, boxes, widths)):
        # distance inside the box's edge, in units of half the taper width
        # (a width of 0 gives a near-step, over a thousandth of a cell)
        from_centre = np.abs(np.arange(n) - (n - 1) / 2.0)
        t = (box / 2.0 - from_centre) / (width / 2.0 if width else 0.001)
        weight = np.where(t >= 1, 1.0, np.where(t <= 0, 0.0, taper(np.clip(t, 0, 1))))
        shape = [1] * ndim
        shape[axis] = n
        mask = mask * weight.reshape(shape)
    return mask


def embed_data(indata, big_shape, pad_size=None):
    """
    indata, less its mean, centred in an array of zeros of shape big_shape. With
    pad_size, the pad_size cells around it hold the nearest edge value of the
    data (less the mean). Returns the array, and the row and column slices that
    hold the data.
    """
    if big_shape is None or None in tuple(big_shape):
        raise ValueError("big_shape is needed")
    ny, nx = indata.shape
    big_y, big_x = big_shape
    if big_y < ny or big_x < nx:
        raise ValueError("big_shape must be at least the shape of indata")

    y0, x0 = (big_y - ny) // 2, (big_x - nx) // 2
    rows, cols = slice(y0, y0 + ny, 1), slice(x0, x0 + nx, 1)
    centred = indata - np.mean(indata)
    out = np.zeros(big_shape, dtype=indata.dtype)
    if pad_size is None:
        out[rows, cols] = centred
    else:
        p = pad_size
        out[y0 - p:y0 + ny + p, x0 - p:x0 + nx + p] = np.pad(centred, p, mode="edge")
    return out, rows, cols
