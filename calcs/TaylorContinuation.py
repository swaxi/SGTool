"""
Variable-height continuation (Cordell & Grauch, 1985).

Data g observed on an uneven surface z_obs(x, y) are continued to another
surface z_out(x, y), either of which may be flat, in three steps, each skipped
when it is not needed:

1. flatten: data on z_obs -> the level z1 = min(z_obs), by the Cordell &
   Grauch iteration

       f1  <-  f1 + (g - forward(f1)),     starting with f1 = g,

   where forward(f1) is the field of the level at each observation point.
   As every point is at or above the level, forward is an upward
   continuation, calculated exactly: the level is continued by FFT to a stack
   of levels about a cell apart and interpolated to each point's height
   (cubic Lagrange). At a wavenumber k and height dz above the level the
   iteration multiplies the error by 1 - exp(-k dz), which lies between 0
   and 1, so it converges for every wavelength without any damping (short
   wavelengths high above the level only slowly, but those are the ones that
   matter least once the field is continued upwards).
2. continue the level from z1 to z2 = mean(z_out)        (FFT, up or down)
3. drape: level z2 -> surface z_out by Taylor series about the level

       f(x, y, z2 + dz) = sum_n  dz^n / n!  *  d^n f / dz^n (x, y, z2)

   with the vertical derivatives by FFT (d/dz = -|k| upwards).

Downward continuation (step 2 when z2 < z1, and step 3 where the output is
below z2) amplifies short wavelengths, so there the derivatives are
low-passed and the continuation Tikhonov-stabilised above a smoothing
wavelength. A draped survey continued upwards to a constant height therefore
needs no damping at all; a barometric survey draped onto the topography does.

No QGIS dependency: numpy and scipy only.
"""

import math

import numpy as np

try:  # as part of the plugin package
    from .detrend import robust_plane
except ImportError:  # calcs on the path as top-level modules
    from detrend import robust_plane


MAX_GAIN = 20.0  # largest gain of an automatically stabilised downward continuation
MAX_LEVELS = 64  # most levels in the stack the flattening interpolates between


class TaylorContinuationError(Exception):
    """The continuation cannot be run with these inputs (the message says why)."""


def _fill_nearest(a):
    """NaN cells filled with the nearest valid value."""
    from scipy.ndimage import distance_transform_edt

    bad = ~np.isfinite(a)
    if not bad.any():
        return a
    idx = distance_transform_edt(bad, return_distances=False, return_indices=True)
    return a[tuple(idx)]


def _as_grid(z, shape, what):
    if np.isscalar(z):
        return np.full(shape, float(z))
    z = np.asarray(z, dtype=float)
    if z.shape != shape:
        raise TaylorContinuationError(f"The {what} grid is not the same size as the data")
    return z


class _Spectral:
    """FFT filtering of grids of one shape: mirror padding tapered to zero within
    the padding, so the edges neither wrap round nor ring."""

    def __init__(self, shape, dx, dy, pad, check=None):
        from scipy.fft import next_fast_len, rfftfreq, fftfreq

        ny, nx = shape
        self.shape = shape
        self.pad = pad
        self.check = check
        self.full = (
            next_fast_len(ny + 2 * pad, real=True),
            next_fast_len(nx + 2 * pad, real=True),
        )
        ky = fftfreq(self.full[0], dy) * 2 * np.pi
        kx = rfftfreq(self.full[1], dx) * 2 * np.pi
        self.k = np.sqrt(kx[None, :] ** 2 + ky[:, None] ** 2)
        # cosine taper from 1 at the grid edge to 0 at the end of the padding
        def taper(n, total):
            t = np.ones(total)
            if pad:
                ramp = 0.5 * (1 + np.cos(np.pi * (np.arange(1, pad + 1) / (pad + 1))))
                t[pad - 1 :: -1][:pad] = ramp
                t[pad + n : pad + n + pad] = ramp
                t[pad + n + pad :] = 0.0
            return t

        self.window = np.outer(taper(ny, self.full[0]), taper(nx, self.full[1]))

    def transform(self, f):
        from scipy.fft import rfft2

        if self.check is not None:
            self.check()
        p = self.pad
        g = np.pad(f, p, mode="reflect") if p else f
        buf = np.zeros(self.full)
        buf[: g.shape[0], : g.shape[1]] = g
        return rfft2(buf * self.window, workers=-1)

    def inverse(self, F):
        from scipy.fft import irfft2

        if self.check is not None:
            self.check()
        p = self.pad
        ny, nx = self.shape
        return irfft2(F, s=self.full, workers=-1)[p : p + ny, p : p + nx]


def _lowpass(k, wavelength):
    """Butterworth-like low pass, 0.5 at the given wavelength (none if 0)."""
    if not wavelength or wavelength <= 0:
        return np.ones_like(k)
    kc = 2 * np.pi / wavelength
    return 1.0 / (1.0 + (k / kc) ** 8)


def _continuation(k, h, wavelength):
    """Continuation by h (m, up > 0). Downwards it is stabilised as Tikhonov's
    exp(k|h|) / (1 + exp(2|h|(k - kc))): the gain rises to about
    exp(kc |h|) / 2 at the smoothing wavelength and falls after it."""
    if h >= 0:
        return np.exp(-k * h)
    h = -h
    kc = 2 * np.pi / wavelength
    with np.errstate(over="ignore"):
        return 1.0 / (np.exp(-k * h) + np.exp(h * (k - 2 * kc)))


class _Series:
    """The Taylor series about a level, for a fixed dz grid."""

    def __init__(self, spec, dz, order, wavelength):
        self.spec = spec
        self.dz = dz
        self.order = int(order)
        lp = _lowpass(spec.k, wavelength)
        # (-|k|)^n / n! * low pass, for n = 1..order
        self.ops = [
            ((-spec.k) ** n / math.factorial(n)) * lp for n in range(1, self.order + 1)
        ]

    def forward(self, f0):
        """Field on the surface level + dz, from the field f0 on the level."""
        F = self.spec.transform(f0)
        out = f0.copy()
        dzn = np.ones_like(self.dz)
        for op in self.ops:
            dzn = dzn * self.dz
            out += dzn * self.spec.inverse(F * op)
        return out


def _lagrange4(u):
    """Cubic Lagrange weights at nodes 0, 1, 2, 3 for positions u (array)."""
    return np.stack(
        [
            -(u - 1) * (u - 2) * (u - 3) / 6.0,
            u * (u - 2) * (u - 3) / 2.0,
            -u * (u - 1) * (u - 3) / 2.0,
            u * (u - 1) * (u - 2) / 6.0,
        ]
    )


class _Upward:
    """Exact upward continuation from a level to points at heights dz >= 0
    above it: the level continued by FFT to a stack of levels (about a cell
    apart, from 0 up to past the highest point) and interpolated to each
    point's height between the four nearest."""

    def __init__(self, spec, dz, step):
        self.spec = spec
        top = max(float(dz.max()), 0.0)
        n = int(np.ceil(top / step)) + 3
        if n > MAX_LEVELS:  # a lot of relief: levels further apart
            step = top / (MAX_LEVELS - 3)
            n = MAX_LEVELS
        self.heights = step * np.arange(n)
        t = np.clip(dz, 0.0, None) / step
        self.i0 = np.clip(np.floor(t).astype(int) - 1, 0, n - 4)
        self.w = _lagrange4(t - self.i0)  # (4, rows, cols)

    def forward(self, f0):
        F = self.spec.transform(f0)
        out = np.zeros_like(f0)
        for level, h in enumerate(self.heights):
            j = level - self.i0
            use = (j >= 0) & (j < 4)
            if not use.any():
                continue
            fl = f0 if h == 0 else self.spec.inverse(F * np.exp(-self.spec.k * h))
            w = np.take_along_axis(self.w, np.clip(j, 0, 3)[None], 0)[0]
            out[use] += w[use] * fl[use]
        return out


def _flatten(g, valid, series, max_iterations, tolerance, progress):
    """Cordell & Grauch iteration: the field on the level whose forward
    continuation reproduces g. Returns (best field, iterations, misfit,
    diverged)."""
    f = g.copy()
    best, best_misfit, worse, it = f, np.inf, 0, 0
    for it in range(1, int(max_iterations) + 1):
        r = g - series.forward(f)
        misfit = float(np.sqrt(np.mean(r[valid] ** 2)))
        if progress is not None:
            progress(it, max_iterations, misfit)
        if not np.isfinite(misfit):
            break
        if misfit < best_misfit:
            best, best_misfit, worse = f, misfit, 0
        else:
            worse += 1
            if worse >= 2:
                break
        if misfit < tolerance:
            break
        # gap cells follow the correction of the nearest data, rather than
        # being fitted to the filled values as if they were data
        f = f + (r if valid.all() else _fill_nearest(np.where(valid, r, np.nan)))
    # stopped getting better well short of the tolerance: diverging
    diverged = best_misfit > max(10 * tolerance, 0.01)
    return best, it, best_misfit, diverged


def variable_continuation(
    data,
    z_obs,
    z_out,
    dx,
    dy,
    order=6,
    smoothing=0.0,
    max_iterations=50,
    tolerance=1e-3,
    pad=None,
    cancel_check=None,
    progress=None,
):
    """
    Continue a grid from one (possibly uneven) surface to another
    (Cordell & Grauch, 1985; see the module notes for the steps).

    Parameters
    ----------
    data : 2-D array
        The field, north-up, NaN where there are no data.
    z_obs, z_out : number or 2-D array
        Height of the observations and of the output surface at each node
        (m, same datum, e.g. above sea level). NaN cells are NaN in the result.
    dx, dy : float
        Cell sizes in metres.
    order : int
        Number of terms (derivatives) in the Taylor series of the draping step.
    smoothing : float
        Wavelength (m) below which downward continuation and the draping
        series are damped, as they blow up short wavelengths. 0 chooses it:
        pi x the largest height difference of the output surface from its
        mean, and enough to keep the gain of a downward continuation below
        about 20; at least four cells. The flattening is never damped.
    max_iterations : int
        Most iterations of the flattening step.
    tolerance : float
        Stop iterating when the misfit to the data falls below this fraction of
        their standard deviation.
    pad : int, optional
        Cells of mirror padding (default a quarter of the larger side, 16 to 512).
    cancel_check : callable, optional
        Called often; raises to abandon the calculation.
    progress : callable, optional
        progress(iteration, max_iterations, relative_misfit) after each iteration.

    Returns
    -------
    (2-D array, dict)
        The field on the output surface, and facts about the run: the levels
        (level_obs, level_out), iterations, relative_misfit (of the flattened
        level continued back onto the observations, as a fraction of the data's
        standard deviation), diverged (the iteration stopped getting better
        well short of the tolerance: the result is the best one found),
        smoothing (the wavelength used for downward continuation and draping,
        or 0 when neither is needed).
    """
    data = np.asarray(data, dtype=float)
    if data.ndim != 2:
        raise TaylorContinuationError("The data must be a grid")
    shape = data.shape
    z_obs = _as_grid(z_obs, shape, "observation height")
    z_out = _as_grid(z_out, shape, "output height")
    dx, dy = float(abs(dx)), float(abs(dy))
    if not (dx > 0 and dy > 0):
        raise TaylorContinuationError("The cell size must be positive")

    valid = np.isfinite(data) & np.isfinite(z_obs) & np.isfinite(z_out)
    if valid.sum() < 16:
        raise TaylorContinuationError(
            "Fewer than 16 cells have data, an observation height and an output height"
        )
    cell = math.sqrt(dx * dy)
    # heights everywhere on the grid, also in the gaps in the data (real
    # heights where the height grid or DEM has them, else the nearest)
    z_obs_all = _fill_nearest(z_obs)
    z1 = float(z_obs_all.min())  # every cell at or above it: nothing clipped
    z2 = float(np.mean(z_out[valid]))
    dz_obs = z_obs_all - z1
    dz_out = _fill_nearest(z_out) - z2
    flat_obs = np.ptp(dz_obs) < 1e-6 * cell
    flat_out = np.ptp(dz_out[valid]) < 1e-6 * cell

    if not smoothing or smoothing <= 0:
        # downward continuation by d: gain at most about MAX_GAIN
        dc = max(z1 - z2, 0.0) * 2 * np.pi / np.log(2 * MAX_GAIN)
        sm_out = max(4 * cell, np.pi * np.abs(dz_out[valid]).max(), dc)
    else:
        sm_out = float(smoothing)
    damped = z2 < z1 - 1e-9 * cell or not flat_out

    # fill the gaps, then remove the regional plane (harmonic: continued
    # unchanged), robust and over the whole filled rectangle (see detrend)
    g = _fill_nearest(np.where(valid, data, np.nan))
    trend, _coef = robust_plane(g)
    g = g - trend
    scale = float(np.std(g[valid])) or 1.0
    g = g / scale

    if pad is None:
        pad = int(min(max(16, max(shape) // 4), 512))
    pad = max(0, int(pad))
    spec = _Spectral(shape, dx, dy, pad, cancel_check)

    # 1. flatten: data on z_obs -> level z1 (all upward: no damping)
    iterations, misfit, diverged = 0, 0.0, False
    if flat_obs:
        f = g
    else:
        f, iterations, misfit, diverged = _flatten(
            g, valid, _Upward(spec, dz_obs, cell), max_iterations, tolerance, progress,
        )

    # 2. level z1 -> level z2 (downward: stabilised at the smoothing wavelength)
    h = z2 - z1
    if abs(h) > 1e-9 * cell:
        f = spec.inverse(spec.transform(f) * _continuation(spec.k, h, sm_out))

    # 3. drape: level z2 -> surface z_out
    if not flat_out:
        f = _Series(spec, dz_out, order, sm_out).forward(f)

    out = f * scale + trend
    out[~valid] = np.nan
    info = {
        "level_obs": z1,
        "level_out": z2,
        "iterations": iterations,
        "relative_misfit": misfit,
        "diverged": diverged,
        "smoothing": float(sm_out) if damped else 0.0,
    }
    return out, info
