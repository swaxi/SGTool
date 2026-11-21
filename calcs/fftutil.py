"""
Fast Fourier filtering of real grids, shared by the FFT filters.

- Real transforms (rfft2 / irfft2): a real grid's spectrum is Hermitian, so only
  half of it is calculated, which halves the time and memory of the transforms.
- All cores: every transform is given workers=-1 (SciPy is single-threaded by
  default).
- Smooth sizes: fast_len() pads to a size whose only prime factors are 2, 3 and 5
  (scipy.fft.next_fast_len), where the transforms are quickest.

A filter is applied as the code using full complex transforms applied it: the
real part of ifft2(fft2(grid) * H). That real part equals the inverse transform
of the grid's spectrum times the Hermitian part of the filter,
(H(k) + conj(H(-k))) / 2, which is what is applied to the half spectrum here.
So the results are those of the full transforms (to rounding), also for filters
that are not Hermitian themselves (e.g. a fractional x or y derivative).
"""

import numpy as np


def fast_len(n):
    """The smallest size >= n with no prime factors but 2, 3 and 5."""
    from scipy.fft import next_fast_len

    return next_fast_len(int(n), real=True)


class RealSpectrum:
    """The half spectrum of a real 2-D grid, to apply filters to.

    The wavenumbers are given to the filters as in the full-grid code:
    k0 along axis 0 (rows), k1 along axis 1 (columns), in radians per unit
    length, k0 negated if flip0 (so that it is positive north on a north-up
    grid), on the grid of the first nx // 2 + 1 columns of the full transform
    (the Nyquist column, for an even width, at its full-grid negative value).
    """

    def __init__(self, grid, d0, d1, flip0=False, workers=-1):
        from scipy.fft import fftfreq, rfft2

        grid = np.asarray(grid, dtype=float)
        self.shape = grid.shape
        self.workers = workers
        n0, n1 = self.shape
        self.F = rfft2(grid, workers=workers)
        h = self.F.shape[1]  # n1 // 2 + 1
        k0 = fftfreq(n0, d0) * 2 * np.pi
        k1 = (fftfreq(n1, d1) * 2 * np.pi)[:h]
        # wavenumbers of the mirror bins (-k): the Nyquist bins mirror onto themselves
        k0_m = k0[(-np.arange(n0)) % n0]
        cols = np.arange(h)
        self_mirror = (cols == 0) | (2 * cols == n1)
        k1_m = np.where(self_mirror, k1, -k1)
        if flip0:
            k0, k0_m = -k0, -k0_m
        # the filters get (axis-1 wavenumber, axis-0 wavenumber) grids
        self.k1, self.k0 = np.meshgrid(k1, k0)
        self.k1_m, self.k0_m = np.meshgrid(k1_m, k0_m)

    def response(self, filter_function, hermitian=False):
        """The Hermitian part of filter_function(k1, k0) on the half spectrum.

        hermitian: the filter is Hermitian as a function, H(-k) = conj(H(k))
        (true of any filter of a real kernel: continuation, reductions to the
        pole, low / high / band pass, integer-order derivatives...). Then it is
        its own Hermitian part except on the Nyquist row and column, whose bins
        mirror onto themselves, so only those are evaluated twice."""
        # as numbers: a filter may be a mask (k <= cutoff), and True + True is True
        H = np.broadcast_to(np.asarray(filter_function(self.k1, self.k0)) + 0.0,
                            self.k1.shape)
        if not hermitian:
            H_m = np.asarray(filter_function(self.k1_m, self.k0_m)) + 0.0
            return 0.5 * (H + np.conj(H_m)) if np.iscomplexobj(H_m) else 0.5 * (H + H_m)
        out = H.copy()
        n0, n1 = self.shape
        lines = []
        if n0 % 2 == 0:
            lines.append((n0 // 2, slice(None)))  # the Nyquist row
        if n1 % 2 == 0:
            lines.append((slice(None), H.shape[1] - 1))  # the Nyquist column
        for idx in lines:
            H_m = np.asarray(filter_function(self.k1_m[idx], self.k0_m[idx])) + 0.0
            out[idx] = 0.5 * (H[idx] + np.conj(H_m))
        return out

    def apply(self, filter_function, hermitian=False):
        """The real grid filtered by filter_function(k1, k0), full size."""
        from scipy.fft import irfft2

        return irfft2(self.F * self.response(filter_function, hermitian), s=self.shape,
                      workers=self.workers)
