from unittest import result
import numpy as np
from numpy.polynomial.polynomial import polyval2d
from math import ceil, fabs
import numpy as np
import matplotlib.pyplot as plt

from osgeo import ogr, osr

import csv

import glob
import os
from osgeo import gdal


from .worms.wormer import Wormer
from .worms.Utility import (
    GetExtent,
    numpy_array_to_raster,
    insert_text_before_extension,
    fill_nan,
)
from .text_io import write_rows
from .sgt_cancel import OperationCancelled


import os


class GeophysicalProcessor:
    def __init__(self, dx, dy, buffer):
        """
        Initialize the processor with grid spacing.

        Parameters:
            dx (float): Grid spacing in the x-direction.
            dy (float): Grid spacing in the y-direction.
        """

        self.dx = dx
        self.dy = dy
        self.buffer = buffer

    # --- Utility Methods ---
    def create_wavenumber_grids(self, nx, ny):
        """
        Create wavenumber grids for Fourier-based operations.
        """
        kx = np.fft.fftfreq(nx, self.dx) * 2 * np.pi
        ky = np.fft.fftfreq(ny, self.dy) * 2 * np.pi
        ky = -ky
        kx_grid, ky_grid = np.meshgrid(kx, ky)
        return kx_grid, ky_grid
        # return np.meshgrid(kx, ky)

    def fill_nan(self, data):
        """
        Replace NaN values using 2D interpolation WITHOUT directional bias.
        Uses griddata for true 2D interpolation instead of sequential 1D passes.
        """
        from scipy.ndimage import distance_transform_edt

        # Create a copy of the data
        filled_data = np.copy(data)
        nan_mask = np.isnan(data)

        # If there are no NaNs or all NaNs, handle accordingly
        if not np.any(nan_mask) or np.all(nan_mask):
            return filled_data, nan_mask

        # Nearest-valid-cell fill. The distance transform returns, for every
        # cell, the index of the nearest non-NaN cell, in linear time and
        # without building point lists or a KD-tree (much faster and lighter
        # than griddata on large grids).
        idx = distance_transform_edt(
            nan_mask, return_distances=False, return_indices=True
        )
        filled_data[nan_mask] = data[idx[0][nan_mask], idx[1][nan_mask]]

        return filled_data, nan_mask

    def restore_nan(self, data, nan_mask):
        """
        Restore NaN values to their original positions in the grid.
        """
        data_with_nan = np.copy(data)
        data_with_nan[nan_mask] = np.nan
        return data_with_nan

    def find_optimal_fft_size(self, dimension, min_padding_each_side):
        """
        Find the optimal size for FFT computation that is a product of powers of 2, 3, and 5,
        ensuring sufficient padding on each side of the grid.

        Args:
            dimension (int): Original dimension size
            min_padding_each_side (int): Minimum padding required on each side

        Returns:
            tuple: (optimal_size, left_padding, right_padding)
        """
        # Minimum size needed with padding on both sides
        min_size = dimension + 2 * min_padding_each_side

        # Find the next size that is a product of powers of 2, 3, and 5
        optimal_size = min_size

        while True:
            # Test if the current size is a product of powers of 2, 3, and 5 only
            n = optimal_size
            for p in [2, 3, 5]:
                while n % p == 0:
                    n //= p

            # If n becomes 1, it means size was a product of powers of 2, 3, 5 only
            if n == 1:
                break

            # Try the next size
            optimal_size += 1

        # Calculate the padding on each side
        total_padding = optimal_size - dimension

        # Handle asymmetric padding (for odd total_padding)
        left_padding = total_padding // 2
        right_padding = total_padding - left_padding

        return optimal_size, left_padding, right_padding

    def optimize_grid_size(self, x_dim, y_dim, min_padding):
        """
        Optimize a 2D grid size for FFT computation with padding on all sides.

        Args:
            x_dim (int): Original X dimension
            y_dim (int): Original Y dimension
            min_padding (int): Minimum required padding on each side

        Returns:
            dict: Contains optimal dimensions and padding information
        """
        # Find optimal sizes and padding for each dimension
        optimal_x, left_x_pad, right_x_pad = self.find_optimal_fft_size(
            x_dim, min_padding
        )
        optimal_y, top_y_pad, bottom_y_pad = self.find_optimal_fft_size(
            y_dim, min_padding
        )

        return {
            "original_dims": (x_dim, y_dim),
            "optimal_dims": (optimal_x, optimal_y),
            "padding": {
                "left": left_x_pad,
                "right": right_x_pad,
                "top": top_y_pad,
                "bottom": bottom_y_pad,
            },
            "total_padding": {
                "x_axis": left_x_pad + right_x_pad,
                "y_axis": top_y_pad + bottom_y_pad,
            },
        }

    def add_buffer(self, data, buffer_size, method="mirror"):
        """
        Add a buffer around the edges to reduce edge effects.
        """
        self.data_mean = np.nanmean(data)
        # data_zero = data - self.data_mean

        padding_info = self.optimize_grid_size(
            data.shape[1], data.shape[0], buffer_size
        )

        # Pad the image with zeros

        padding = padding_info["padding"]
        padded_image = np.pad(
            data,
            pad_width=(
                (padding["top"], padding["bottom"]),
                (padding["left"], padding["right"]),
            ),
            mode="linear_ramp",
        )

        # Create a 2D Hanning window
        window_y = np.hanning(padded_image.shape[0])
        window_x = np.hanning(padded_image.shape[1])
        hanning_window = np.outer(window_y, window_x)

        # Lower alpha means less windowing effect
        alpha = 0.5  # Adjust this value to control strength
        modified_window = 1.0 - alpha * (1.0 - hanning_window)

        # Apply the modified window
        # padded_image = padded_image * modified_window

        alpha = 0.50  # Only taper 50% from each edge
        window_y = self.tukey_window(padded_image.shape[0], alpha)
        window_x = self.tukey_window(padded_image.shape[1], alpha)
        tukey_window = np.outer(window_y, window_x)

        # Apply the Tukey window
        padded_image = padded_image * tukey_window  # * modified_window
        return padded_image, padding_info

    def tukey_window(self, M, alpha=0.25):
        """Create a Tukey window with length M and shape parameter alpha."""
        if alpha <= 0:
            return np.ones(M)  # Rectangular window
        elif alpha >= 1:
            return np.hanning(M)  # Hanning window

        # Create the window
        x = np.linspace(0, 1, M)
        w = np.ones(M)

        # Left taper
        left_taper = x < alpha / 2
        w[left_taper] = 0.5 * (
            1 + np.cos(2 * np.pi / alpha * (x[left_taper] - alpha / 2))
        )

        # Right taper
        right_taper = x >= (1 - alpha / 2)
        w[right_taper] = 0.5 * (
            1 + np.cos(2 * np.pi / alpha * (x[right_taper] - 1 + alpha / 2))
        )

        return w

    def remove_buffer(self, data, buffer_size, padding_info):
        """
        Remove the buffer from the edges of the grid.
        """
        data = data  # + self.data_mean

        # Remove the buffer
        return data[
            padding_info["padding"]["top"] : padding_info["padding"]["top"]
            + padding_info["original_dims"][1],
            padding_info["padding"]["left"] : padding_info["padding"]["left"]
            + padding_info["original_dims"][0],
        ]

    # --- Trend Removal ---
    def remove_polynomial_trend(self, x, y, z, grid_x, grid_y, degree=2):
        """
        Remove a regional trend using polynomial fitting.
        """
        X = np.column_stack(
            [x**i * y**j for i in range(degree + 1) for j in range(degree + 1 - i)]
        )
        coeffs, _, _, _ = np.linalg.lstsq(X, z, rcond=None)

        # Evaluate polynomial on grid
        trend = polyval2d(grid_x, grid_y, coeffs.reshape(degree + 1, -1))
        residual = z - polyval2d(x, y, coeffs.reshape(degree + 1, -1))

        return residual, trend

    def remove_regional_trend_fourier(
        self, data, cutoff_wavelength, buffer_size=10, buffer_method="mirror"
    ):
        """
        Remove the regional trend using a low-pass filter in the Fourier domain.
        """
        return self.low_pass_filter(data, cutoff_wavelength, buffer_size, buffer_method)

    def vertical_integration(
        self,
        data,
        max_wavenumber=None,
        min_wavenumber=1e-6,
        buffer_size=10,
        buffer_method="mirror",
    ):
        """
        Perform vertical integration of a field in the frequency domain with stabilization.

        Parameters:
            data (numpy.ndarray): 2D array of the field data.
            max_wavenumber (float): Optional maximum wavenumber for filtering (to stabilize high-frequency noise).
            min_wavenumber (float): Minimum wavenumber to prevent amplification of low-frequency components.
            buffer_size (int): Buffer size for edge handling.
            buffer_method (str): Buffering method ('mirror' or 'zero').

        Returns:
            numpy.ndarray: Vertically integrated data.
        """

        def filter_function(kx, ky):
            # Compute wavenumber magnitude
            k = np.sqrt(kx**2 + ky**2)

            # Stabilize low and high wavenumbers
            k = np.maximum(k, min_wavenumber)  # Apply low-wavenumber stabilization
            if max_wavenumber is not None:
                k = np.minimum(k, max_wavenumber)  # Apply high-wavenumber cutoff

            # Vertical integration filter
            vertical_integration_filter = 1 / k
            return vertical_integration_filter

        # Apply the Fourier filter
        return self._apply_fourier_filter(
            data, filter_function, buffer_size, buffer_method
        )

    # --- Apparent petrophysics ---
    @staticmethod
    def _high_cut_taper(k, cutoff_wavelength):
        """1 for wavelengths longer than twice the cutoff, cosine roll-off to 0 at
        the cutoff wavelength itself (stops the inverse filters blowing up)."""
        if not cutoff_wavelength or cutoff_wavelength <= 0:
            return np.ones_like(k)
        k_cut = 2 * np.pi / cutoff_wavelength
        k_pass = k_cut / 2
        taper = np.ones_like(k)
        taper[k >= k_cut] = 0.0
        roll = (k > k_pass) & (k < k_cut)
        taper[roll] = 0.5 * (1 + np.cos(np.pi * (k[roll] - k_pass) / (k_cut - k_pass)))
        return taper

    def _apply_layer_inverse(self, data, inverse, buffer_size, buffer_method):
        """Run a layer-inversion filter. inverse(k) gives the filter for radial
        wavenumber k. The k=0 term is dropped: the absolute level of a gravity
        or magnetic anomaly grid (datum, regional, IGRF removal) is arbitrary, so
        the result has a zero mean and the caller adds a background value."""

        def filter_function(kx, ky):
            k = np.sqrt(kx**2 + ky**2)
            return np.where(k > 0, inverse(k), 0.0)

        return self._apply_fourier_filter(
            data, filter_function, buffer_size, buffer_method, preserve_dc=False
        )

    def apparent_susceptibility(
        self,
        data,
        field_nT,
        depth=0.0,
        thickness=None,
        background=0.0,
        cutoff_wavelength=None,
        buffer_size=10,
        buffer_method="mirror",
    ):
        """
        Apparent susceptibility (SI) from a reduced-to-pole TMI grid (nT), for a
        layer of vertical prisms magnetised by induction only.

        The pole response of a layer with top at `depth` below the sensor and
        thickness t is  T(k) = chi(k) F/2 exp(-k depth) (1 - exp(-k t))
        (t infinite: no last term), which is divided out here.

        Parameters:
            field_nT (float): Ambient (IGRF) field strength in nT.
            depth (float): Depth from the sensor to the top of the layer (m).
            thickness (float or None): Layer thickness (m), None for infinite.
            background (float): Susceptibility (SI) added to the result, since
                an anomaly grid has no absolute level.
            cutoff_wavelength (float or None): High-cut wavelength (m).

        Returns:
            numpy.ndarray: Apparent susceptibility (SI), with the background
            as its mean.
        """

        def inverse(k):
            taper = self._high_cut_taper(k, cutoff_wavelength)
            # upward continuation term, inverted; capped so exp never overflows
            response = 2.0 / field_nT * np.exp(np.minimum(k * depth, 700.0))
            if thickness:
                with np.errstate(divide="ignore"):
                    response = np.where(k > 0, response / (1 - np.exp(-k * thickness)), 0.0)
            return np.where(taper > 0, response * taper, 0.0)

        return background + self._apply_layer_inverse(
            data, inverse, buffer_size, buffer_method
        )

    def apparent_density(
        self,
        data,
        depth=0.0,
        thickness=1000.0,
        background=2.67,
        cutoff_wavelength=None,
        buffer_size=10,
        buffer_method="mirror",
    ):
        """
        Apparent density (g/cc) from a Bouguer gravity grid (mGal), for a layer of
        vertical prisms.

        The response of a layer with top at `depth` below the sensor and
        thickness t is  g(k) = 2 pi G drho(k) exp(-k depth) (1 - exp(-k t)) / k,
        which is divided out here. The mean level is not used (a Bouguer grid's
        level depends on its datum and the regional field).

        Parameters:
            depth (float): Depth from the sensor to the top of the layer (m).
            thickness (float): Layer thickness (m), must be > 0.
            background (float): Density (g/cc) added to the density contrast.
            cutoff_wavelength (float or None): High-cut wavelength (m).

        Returns:
            numpy.ndarray: Apparent density (g/cc), with the background as its mean.
        """
        G = 6.674e-11
        # mGal -> m/s^2, and kg/m^3 -> g/cc
        scale = 1e-5 / 1000.0

        def inverse(k):
            taper = self._high_cut_taper(k, cutoff_wavelength)
            with np.errstate(divide="ignore", invalid="ignore"):
                slab = np.where(
                    k > 0, (1 - np.exp(-k * thickness)) / np.where(k > 0, k, 1.0), thickness
                )
                response = scale * np.exp(np.minimum(k * depth, 700.0)) / (2 * np.pi * G * slab)
            return np.where(taper > 0, response * taper, 0.0)

        return background + self._apply_layer_inverse(
            data, inverse, buffer_size, buffer_method
        )

    # --- Continuation ---
    def upward_continuation(self, data, height, buffer_size=10, buffer_method="mirror"):
        """
        Perform upward continuation to attenuate high-frequency signals.
        """

        def filter_function(kx, ky):
            k = np.sqrt(kx**2 + ky**2)
            return np.exp(-k * height)

        return self._apply_fourier_filter(
            data, filter_function, buffer_size, buffer_method
        )

    def downward_continuation(
        self, data, height, buffer_size=10, buffer_method="mirror"
    ):
        """
        Perform downward continuation to enhance high-frequency signals.
        """

        def filter_function(kx, ky):
            k = np.sqrt(kx**2 + ky**2)
            return np.exp(k * height)

        return self._apply_fourier_filter(
            data, filter_function, buffer_size, buffer_method
        )

    def total_hz_grad(self, data, buffer_size=10, buffer_method="mirror"):
        """
        Compute the total horizontal gradient (THG) of a 2D grid using Fourier filtering.

        Parameters:
            data (numpy.ndarray): Input 2D grid data.
            buffer_size (int): Buffer size for edge handling.
            buffer_method (str): Buffering method ('mirror' or 'zero').

        Returns:
            numpy.ndarray: Total horizontal gradient of the input data.
        """

        def filter_function_dx(kx, ky):
            """
            Fourier filter for the x-derivative.
            """
            return 1j * kx

        def filter_function_dy(kx, ky):
            """
            Fourier filter for the y-derivative.
            """
            return 1j * ky

        # Compute derivatives in x and y directions using the Fourier filter
        dx, dy = self._apply_fourier_filters(
            data,
            [filter_function_dx, filter_function_dy],
            buffer_size,
            buffer_method,
            preserve_dc=False,
        )

        print("dx range:", np.nanmin(dx), np.nanmax(dx))
        print("dy range:", np.nanmin(dy), np.nanmax(dy))

        # Compute the total horizontal gradient
        thg = np.sqrt(dx**2 + dy**2)
        print("thg range:", np.nanmin(thg), np.nanmax(thg))

        return thg

    def reduction_to_pole(
        self, data, inclination, declination, buffer_size=10, buffer_method="mirror"
    ):
        """
        Apply reduction to the pole filter in the frequency domain.

        Parameters
        ----------
        data : array-like
            Input magnetic data to be transformed
        inclination : float
            The inclination of the inducing Geomagnetic field in degrees
        declination : float
            The declination of the inducing Geomagnetic field in degrees
        buffer_size : int
            Size of the buffer zone for reducing edge effects
        buffer_method : str
            Method for handling the buffer zone ("mirror", "mean", etc.)

        Returns
        -------
        array-like
            The reduced to pole magnetic data
        """

        def filter_function(kx, ky):
            # Transform to radians
            inc_rad = np.radians(
                -inclination
            )  # Negative sign to correct for convention
            dec_rad = np.radians(declination)
            # Calculate the 3 components
            m_e = np.cos(inc_rad) * np.sin(dec_rad)
            m_n = np.cos(inc_rad) * np.cos(dec_rad)
            m_z = -np.sin(inc_rad)
            f_e, f_n, f_z = m_e, m_n, m_z
            k_squared = ky**2 + kx**2
            k = np.sqrt(k_squared) + 1e-10

            # RTP filter
            # with np.errstate(divide="ignore", invalid="ignore"):
            rtp_filter = (
                k_squared
                * (f_z * k + 1j * (f_e * kx + f_n * ky)) ** (-1)
                * (m_z * k + 1j * (m_e * kx + m_n * ky)) ** (-1)
            )

            rtp_filter[np.abs(rtp_filter) > 1e6] = 0  # Stabilize extreme values
            return rtp_filter

        return self._apply_fourier_filter(
            data, filter_function, buffer_size, buffer_method
        )

    def reduction_to_equator(
        self, data, inclination, declination, buffer_size=10, buffer_method="mirror"
    ):
        inc, dec = np.radians(inclination), np.radians(declination)

        def filter_function(kx, ky):
            return self._rte_operator(kx, ky, inc, dec)

        # mirror padding: the operator changes sharply near the direction
        # perpendicular to the declination and would smear the streaks that
        # ramp padding makes along the declination into the grid
        return self._apply_fourier_filter(
            data, filter_function, buffer_size, buffer_method, mirror_pad=True
        )

    @staticmethod
    def _rte_operator(kx, ky, inc, dec):
        """
        Reduction to the equator operator (Kis 1990; Hao et al. 2018, Eq. 5),
        inc and dec in radians:

            C = -cos^2(D-theta) / [sin I + i cos I cos(D-theta)]^2

        evaluated as -cos^2 * conj(f)^2 / |f|^4 so the denominator is real.
        At inclination 0 it is 1 (no change).
        """
        theta = np.arctan2(kx, ky)  # azimuth of the wavenumber, from north
        cos_term = np.cos(dec - theta)
        field = np.sin(inc) + 1j * np.cos(inc) * cos_term
        power = np.abs(field) ** 2
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(
                power > 1e-12,
                -(cos_term**2) * np.conj(field) ** 2 / power**2,
                1.0,  # inclination 0, perpendicular to the declination: the limit
            )

    # iteration counts tried when choosing m automatically
    TSI_RTP_ITERATIONS = (0, 1, 2, 3, 4, 6, 8, 10, 13, 16, 20, 25, 32, 40, 50,
                          63, 80, 100, 125, 160, 200)

    def reduction_to_pole_low_latitude(
        self, data, inclination, declination, iterations=0, stripe_width=0.0,
        buffer_size=10, buffer_method="mirror",
    ):
        """
        Reduction to the pole at low latitudes by the Taylor series iterative
        (TSI) method of Hao et al. (2018, J. Appl. Geophys. 159, 127-134).

        The data are first reduced to the equator (Eq. 5), then to the pole with
        m iterations of a truncated Taylor-series RTP operator, in closed form
        (Eq. 13):  Tp = [1 - (1 - Cpe*Hea)^(m+1)] * He * Te,
        with He = -sec^2(D-theta), Cpe = -cos^2(D-theta) and
        Hea = -(1 + (D-theta)^2/2)^2. More iterations approach the exact
        (unstable) RTP; the gain perpendicular to the declination is at most
        about 5(m+1).

        Parameters
        ----------
        iterations : int
            Number of iterations m; 0 or less chooses m automatically as the
            one whose result correlates best with the analytic signal
            amplitude of the data (Eq. 14).
        stripe_width : float
            Suppression of the stripes parallel to the declination (degrees,
            0 = none, as in the paper). Within about this angle of the direction
            perpendicular to the declination, where the operator has to amplify
            most and the data hold almost no signal at low latitude, the number
            of iterations fades to 0: m_eff = m * (1 - exp(-(phi / width)^2)),
            phi the angle of the wavenumber from that direction.

        Returns
        -------
        (array-like, int)
            The reduced to pole data, and the number of iterations used.
        """
        from scipy.fft import fft2, ifft2

        inc, dec = np.radians(inclination), np.radians(declination)

        def filter_function(kx, ky, m, width=0.0):
            rte = self._rte_operator(kx, ky, inc, dec)  # Eq. 5
            theta = np.arctan2(kx, ky)
            cos_term = np.cos(dec - theta)
            # angle from the declination, folded to [-90, 90] degrees (sec^2 has period pi)
            x = np.angle(np.exp(2j * (dec - theta))) / 2
            h_ea = -((1 + x**2 / 2) ** 2)  # Eq. 4
            q = cos_term**2 * (1 + x**2 / 2) ** 2  # Cpe * Hea, in [0, 1]
            if width and width > 0:
                phi = np.pi / 2 - np.abs(x)  # angle from the perpendicular to D
                m = m * -np.expm1(-((phi / np.radians(width)) ** 2))
            # [1 - (1-q)^(m+1)] * He = Hea * [1 - (1-q)^(m+1)] / q, -> Hea (m+1) as q -> 0
            with np.errstate(divide="ignore", invalid="ignore"):
                gain = np.where(
                    q > 1e-12,
                    -np.expm1((m + 1) * np.log1p(-np.minimum(q, 1 - 1e-15))) / q,
                    m + 1.0,
                )
            return rte * h_ea * gain

        # Mirror padding, tapered to zero within the padding only: the usual
        # ramp to zero extends each column (or row) straight out from the edge,
        # and where those streaks run along the declination this operator
        # amplifies them into stripes across the whole grid
        median = np.nanmedian(data)
        filled_data, nan_mask = self.fill_nan(data)
        buffered_data, padding_info = self._mirror_buffer(filled_data - median, buffer_size)
        data_fft = fft2(buffered_data, workers=-1)
        ny, nx = buffered_data.shape
        kx, ky = self.create_wavenumber_grids(nx, ny)

        if iterations and iterations > 0:
            m = int(iterations)
            result = np.real(ifft2(
                data_fft * filter_function(kx, ky, m, stripe_width), workers=-1
            ))
            result = self.remove_buffer(result, buffer_size, padding_info)
            return self.restore_nan(result + median, nan_mask), m

        # an inverse transform per trial iteration count. m is chosen without
        # the stripe suppression, which is then applied: the correlation
        # would otherwise just raise m to make up for it
        asa = self.analytic_signal(data, buffer_size, buffer_method)
        valid = ~nan_mask & np.isfinite(asa)
        a = asa[valid] - np.mean(asa[valid])

        best, best_m, best_corr = None, 0, -np.inf
        for m in self.TSI_RTP_ITERATIONS:
            result = np.real(ifft2(data_fft * filter_function(kx, ky, m), workers=-1))
            result = self.remove_buffer(result, buffer_size, padding_info)
            t = result[valid] - np.mean(result[valid])
            denom = np.sqrt(np.sum(t**2) * np.sum(a**2))
            corr = np.sum(t * a) / denom if denom > 0 else -np.inf
            if corr > best_corr:
                best, best_m, best_corr = result, m, corr

        if stripe_width and stripe_width > 0:
            best = np.real(ifft2(
                data_fft * filter_function(kx, ky, best_m, stripe_width), workers=-1
            ))
            best = self.remove_buffer(best, buffer_size, padding_info)

        return self.restore_nan(best + median, nan_mask), best_m

    def _mirror_buffer(self, data, buffer_size):
        """Pad to an FFT-friendly size by mirroring the edges, with a cosine
        taper to zero across the padding (the data themselves are untouched).
        Returns the padded grid and the padding info for remove_buffer."""
        padding_info = self.optimize_grid_size(data.shape[1], data.shape[0], buffer_size)
        pad = padding_info["padding"]

        def taper(n, before, after):
            w = np.ones(n + before + after)
            if before:
                w[:before] = 0.5 * (1 - np.cos(np.pi * np.arange(before) / before))
            if after:
                w[-after:] = 0.5 * (1 + np.cos(np.pi * (np.arange(after) + 1) / after))
            return w

        padded = np.pad(
            data,
            ((pad["top"], pad["bottom"]), (pad["left"], pad["right"])),
            mode="symmetric",
        )
        window = np.outer(
            taper(data.shape[0], pad["top"], pad["bottom"]),
            taper(data.shape[1], pad["left"], pad["right"]),
        )
        return padded * window, padding_info

    # --- Derivatives ---

    def compute_derivative(
        self, data, direction, order=1, buffer_size=10, buffer_method="mirror"
    ):
        """
        Compute derivatives in x, y, or z directions.
        """
        return self.compute_derivatives(
            data, [direction], order, buffer_size, buffer_method
        )[0]

    def compute_derivatives(
        self, data, directions, order=1, buffer_size=10, buffer_method="mirror"
    ):
        """
        Derivatives of the same grid in several directions (each as
        compute_derivative gives it), sharing one forward transform.
        """

        def derivative_filter(direction):
            def filter_function(kx, ky):
                if direction == "x":
                    return (1j * kx) ** order
                elif direction == "y":
                    return (1j * ky) ** order
                elif direction == "z":
                    k = np.sqrt(kx**2 + ky**2)
                    return k**order
                else:
                    raise ValueError("Invalid direction. Choose 'x', 'y', or 'z'.")

            return filter_function

        return self._apply_fourier_filters(
            data,
            [derivative_filter(d) for d in directions],
            buffer_size,
            buffer_method,
            preserve_dc=False,
        )

    def tilt_angle(self, data, buffer_size=10, buffer_method="mirror"):
        """
        Compute the tilt derivative of the data.
        """
        dfdx, dfdy, dfdz = self.compute_derivatives(
            data, "xyz", buffer_size=buffer_size, buffer_method=buffer_method
        )
        horizontal_gradient = np.sqrt(dfdx**2 + dfdy**2)
        return np.arctan2(dfdz, horizontal_gradient)

    def analytic_signal(self, data, buffer_size=10, buffer_method="mirror"):
        """
        Compute the analytic signal amplitude.
        """
        dfdx, dfdy, dfdz = self.compute_derivatives(
            data, "xyz", buffer_size=buffer_size, buffer_method=buffer_method
        )
        return np.sqrt(dfdx**2 + dfdy**2 + dfdz**2)

    def automatic_gain_control(self, grid, window_size):
        """
        Apply automatic gain control (AGC) to a 2D NumPy grid with edge handling.

        Parameters:
            grid (numpy.ndarray): 2D NumPy array representing the grid.
            window_size (int or float): Size of the smoothing window (standard deviation for Gaussian filter).

        Returns:
            numpy.ndarray: Grid after applying AGC.
        """
        from scipy.ndimage import gaussian_filter

        filled_data, nan_mask = self.fill_nan(grid)
        # Take the absolute value of the grid to normalize amplitudes
        abs_grid = np.abs(filled_data)

        # Pad the grid to handle edges using reflection
        padded_abs_grid = np.pad(
            abs_grid, pad_width=int(3 * window_size), mode="reflect"
        )

        # Apply Gaussian smoothing to the padded grid
        smoothed_padded_grid = gaussian_filter(padded_abs_grid, sigma=window_size)

        # Remove padding
        smoothed_grid = smoothed_padded_grid[
            int(3 * window_size) : -int(3 * window_size),
            int(3 * window_size) : -int(3 * window_size),
        ]

        # Prevent division by zero by adding a small constant
        smoothed_grid[smoothed_grid == 0] = np.finfo(float).eps

        # Apply AGC: normalize the original grid by the smoothed grid
        agc_grid = grid / smoothed_grid

        return self.restore_nan(agc_grid, nan_mask)

    def low_pass_filter(
        self,
        data,
        cutoff_wavelength,
        transition_width=None,
        buffer_size=10,
        buffer_method="mirror",
    ):
        """
        Apply a low-pass filter to remove high-frequency noise.

        Parameters
        ----------
        data : array-like
            Input data to be filtered
        cutoff_wavelength : float
            Cutoff wavelength - features with shorter wavelengths will be attenuated
        transition_width : float or None
            Width of transition zone in the same units as cutoff_wavelength.
            If None, no transition is applied (sharp cutoff).
            Only needed if experiencing ringing artifacts.
        buffer_size : int
            Size of the buffer zone for reducing edge effects
        buffer_method : str
            Method for handling the buffer zone ("mirror", "mean", etc.)

        Returns
        -------
        array-like
            The low-pass filtered data
        """

        def filter_function(kx, ky):
            k = np.sqrt(kx**2 + ky**2)
            cutoff_k = 2 * np.pi / (cutoff_wavelength + 1e-10)

            if transition_width is None or transition_width <= 0:
                # Traditional binary filter (no transition)
                return k <= cutoff_k

            # Calculate transition boundaries
            inner_k = 2 * np.pi / (cutoff_wavelength + transition_width / 2 + 1e-10)
            outer_k = 2 * np.pi / (cutoff_wavelength - transition_width / 2 + 1e-10)
            if cutoff_wavelength - transition_width / 2 <= 0:
                outer_k = cutoff_k * 1.5  # Fallback

            # Create filter with smooth transition
            filter_values = np.zeros_like(k)

            # Frequencies below inner_k are fully passed
            mask_pass = k <= inner_k
            filter_values[mask_pass] = 1.0

            # Frequencies above outer_k are fully blocked
            mask_block = k >= outer_k
            filter_values[mask_block] = 0.0

            # Frequencies in transition zone get smooth cosine taper
            mask_transition = np.logical_and(k > inner_k, k < outer_k)
            if np.any(mask_transition):
                # Normalize position within transition zone to [0, 1]
                pos = (k[mask_transition] - inner_k) / (outer_k - inner_k)
                # Apply cosine taper (decreasing from 1 to 0)
                filter_values[mask_transition] = 0.5 * (1 + np.cos(np.pi * pos))

            return filter_values

        return self._apply_fourier_filter(
            data, filter_function, buffer_size, buffer_method
        )

    def low_pass_filterx(
        self, data, cutoff_wavelength, buffer_size=10, buffer_method="mirror"
    ):
        """
        Apply a low-pass filter to suppress high-frequency noise.
        """

        def filter_function(kx, ky):
            k = np.sqrt(kx**2 + ky**2)
            cutoff_k = 2 * np.pi / (cutoff_wavelength + 1e-10)
            return k < cutoff_k  # Low-pass filter mask

        return self._apply_fourier_filter(
            data, filter_function, buffer_size, buffer_method
        )

    def high_pass_filter(
        self,
        data,
        cutoff_wavelength,
        transition_width=20000,
        buffer_size=10,
        buffer_method="mirror",
    ):
        """
        Apply a high-pass filter to remove regional trends with smoothed transition to reduce ringing.

        Parameters:
        -----------
        data : numpy.ndarray
            Input raster data
        cutoff_wavelength : float
            Wavelength at which to cut off low frequencies (in real-world units)
        transition_width : float or None
            Width of transition zone in the same units as cutoff_wavelength.
            If None, defaults to 20% of the cutoff_wavelength.
        buffer_size : int
            Size of the buffer to add around the edges
        buffer_method : str
            Method to use for buffering ('mirror', 'wrap', etc.)
        """

        # Default transition width if not specified
        if transition_width is None:
            transition_width = 0.2 * cutoff_wavelength

        def filter_function(kx, ky):
            k = np.sqrt(kx**2 + ky**2)

            # Convert wavelengths to wave numbers
            cutoff_k = 2 * np.pi / (cutoff_wavelength + 1e-10)
            # Define transition zone boundaries in wave numbers
            lower_k = 2 * np.pi / (cutoff_wavelength + transition_width / 2 + 1e-10)
            upper_k = 2 * np.pi / (cutoff_wavelength - transition_width / 2 + 1e-10)

            # Handle case where transition_width is too large
            if cutoff_wavelength - transition_width / 2 <= 0:
                upper_k = cutoff_k * 1.5  # Fallback to avoid division by zero

            # Create smooth transition using cosine taper
            filter_values = np.ones_like(k)

            # Values below lower_k are 0 (low frequencies)
            mask_low = k <= lower_k
            filter_values[mask_low] = 0

            # Values above upper_k are 1 (high frequencies)
            mask_high = k >= upper_k
            filter_values[mask_high] = 1

            # Values in between get a smooth transition
            mask_transition = np.logical_and(k > lower_k, k < upper_k)
            if np.any(mask_transition):
                # Normalize position within transition zone to [0, 1]
                pos = (k[mask_transition] - lower_k) / (upper_k - lower_k)
                # Apply cosine taper (smoother than linear)
                filter_values[mask_transition] = 0.5 * (1 - np.cos(np.pi * pos))

            return filter_values

        return self._apply_fourier_filter(
            data, filter_function, buffer_size, buffer_method
        )

    def high_pass_filterx(
        self, data, cutoff_wavelength, buffer_size=10, buffer_method="mirror"
    ):
        """
        Apply a high-pass filter to remove regional trends.
        """

        def filter_function(kx, ky):
            k = np.sqrt(kx**2 + ky**2)
            cutoff_k = 2 * np.pi / (cutoff_wavelength + 1e-10)
            return k > cutoff_k  # High-pass filter mask

        return self._apply_fourier_filter(
            data, filter_function, buffer_size, buffer_method
        )

    def band_pass_filterx(
        self, data, low_cut, high_cut, buffer_size=10, buffer_method="mirror"
    ):
        """
        Apply a band-pass filter to isolate anomalies within a wavelength range.
        """

        def filter_function(kx, ky):
            k = np.sqrt(kx**2 + ky**2)
            high_cut_k = 2 * np.pi / low_cut
            low_cut_k = 2 * np.pi / high_cut

            # Create a band-pass mask
            filter_mask = (k <= low_cut_k) & (k >= high_cut_k)

            return filter_mask  # Band-pass filter mask

        return self._apply_fourier_filter(
            data, filter_function, buffer_size, buffer_method
        )

    def band_pass_filter(
        self,
        data,
        low_cut,
        high_cut,
        high_transition_width=None,
        low_transition_width=None,
        buffer_size=10,
        buffer_method="mirror",
    ):
        """
        Apply a band-pass filter to isolate anomalies within a wavelength range,
        with smooth transitions to reduce ringing artifacts.

        Parameters
        ----------
        data : array-like
            Input data to be filtered
        low_cut : float
            Low cut-off wavelength (features with longer wavelengths will be attenuated)
        high_cut : float
            High cut-off wavelength (features with shorter wavelengths will be attenuated)
        high_transition_width : float or None
            Width of transition zone for the high-pass component in the same units as high_cut.
            If None, defaults to 20% of the high_cut wavelength.
        low_transition_width : float or None
            Width of transition zone for the low-pass component in the same units as low_cut.
            If None, defaults to 20% of the low_cut wavelength.
        buffer_size : int
            Size of the buffer zone for reducing edge effects
        buffer_method : str
            Method for handling the buffer zone ("mirror", "mean", etc.)

        Returns
        -------
        array-like
            The band-pass filtered data
        """

        # Default transition widths if not specified
        if high_transition_width is None:
            high_transition_width = 0.2 * high_cut
        if low_transition_width is None:
            low_transition_width = 0.2 * low_cut

        if low_cut - high_transition_width < 0.001:
            low_cut = low_cut + 0.001
        if fabs(low_cut + high_transition_width) < 0.001:
            low_cut = low_cut + 0.001

        if high_cut - low_transition_width < 0.001:
            high_cut = high_cut + 0.001
        if fabs(high_cut + low_transition_width) < 0.001:
            high_cut = high_cut + 0.001

        if low_cut == 0.0:
            low_cut = 0.001
        if high_cut == 0.0:
            high_cut = 0.001

        def filter_function(kx, ky):
            k = np.sqrt(kx**2 + ky**2)

            # Convert wavelength cutoffs to wavenumber cutoffs
            high_cut_k = (
                2 * np.pi / low_cut
            )  # High-pass cutoff (removes long wavelengths)
            low_cut_k = (
                2 * np.pi / high_cut
            )  # Low-pass cutoff (removes short wavelengths)

            # Create transition zone boundaries for high-pass component
            high_cut_k_lower = 2 * np.pi / (low_cut + high_transition_width / 2)
            high_cut_k_upper = 2 * np.pi / (low_cut - high_transition_width / 2)
            if low_cut - high_transition_width / 2 <= 0:
                high_cut_k_upper = high_cut_k * 1.5  # Fallback

            # Create transition zone boundaries for low-pass component
            low_cut_k_lower = 2 * np.pi / (high_cut + low_transition_width / 2)
            low_cut_k_upper = 2 * np.pi / (high_cut - low_transition_width / 2)
            if high_cut - low_transition_width / 2 <= 0:
                low_cut_k_upper = low_cut_k * 1.5  # Fallback

            # Initialize filter with zeros
            filter_values = np.zeros_like(k)

            # Inside the passband (fully passed frequencies)
            mask_passband = np.logical_and(k >= high_cut_k_upper, k <= low_cut_k_lower)
            filter_values[mask_passband] = 1.0

            # High-pass transition zone
            mask_high_transition = np.logical_and(
                k >= high_cut_k_lower, k < high_cut_k_upper
            )
            if np.any(mask_high_transition):
                # Normalize position within transition zone to [0, 1]
                pos = (k[mask_high_transition] - high_cut_k_lower) / (
                    high_cut_k_upper - high_cut_k_lower
                )
                # Apply cosine taper
                filter_values[mask_high_transition] = 0.5 * (1 - np.cos(np.pi * pos))

            # Low-pass transition zone
            mask_low_transition = np.logical_and(
                k > low_cut_k_lower, k <= low_cut_k_upper
            )
            if np.any(mask_low_transition):
                # Normalize position within transition zone to [0, 1]
                pos = (k[mask_low_transition] - low_cut_k_lower) / (
                    low_cut_k_upper - low_cut_k_lower
                )
                # Apply cosine taper (reversed direction)
                filter_values[mask_low_transition] = 0.5 * (1 + np.cos(np.pi * pos))

            return filter_values

        return self._apply_fourier_filter(
            data, filter_function, buffer_size, buffer_method
        )

    def butterworth_band_pass(
        self, data, low_cut, high_cut, order=4, buffer_size=10, buffer_method="mirror"
    ):
        """
        Apply a Butterworth band-pass filter to isolate anomalies within a wavelength range.

        Parameters
        ----------
        data : array-like
            Input data to be filtered
        low_cut : float
            Low cut-off wavelength (features with longer wavelengths will be attenuated)
        high_cut : float
            High cut-off wavelength (features with shorter wavelengths will be attenuated)
        order : int, optional
            Order of the Butterworth filter. Higher orders create sharper transitions
            but may introduce more ringing artifacts. Default is 4.
        buffer_size : int
            Size of the buffer zone for reducing edge effects
        buffer_method : str
            Method for handling the buffer zone ("mirror", "mean", etc.)

        Returns
        -------
        array-like
            The band-pass filtered data
        """

        def filter_function(kx, ky):
            k = np.sqrt(kx**2 + ky**2)

            # Convert wavelength cutoffs to wavenumber cutoffs
            k_high = 2 * np.pi / high_cut  # High-pass cutoff (removes long wavelengths)
            k_low = 2 * np.pi / low_cut  # Low-pass cutoff (removes short wavelengths)

            # Avoid division by zero
            k = np.maximum(k, 1e-10)

            # Butterworth high-pass component (attenuates low frequencies)
            high_pass = 1.0 / (1.0 + (k_high / k) ** (2 * order))

            # Butterworth low-pass component (attenuates high frequencies)
            low_pass = 1.0 / (1.0 + (k / k_low) ** (2 * order))

            # Combine to create band-pass filter
            band_pass = high_pass * low_pass

            return band_pass

        return self._apply_fourier_filter(
            data, filter_function, buffer_size, buffer_method
        )

    def directional_butterworth_band_pass(
        self,
        data,
        low_cut,
        high_cut,
        direction_angle=0,
        direction_width=45,
        order=4,
        buffer_size=10,
        buffer_method="mirror",
        preserve_dc=True,
    ):
        """
        Apply a combined Butterworth band-pass filter with directional filtering.

        preserve_dc : bool, optional
            If True (default) the grid median is added back to the result.
            Use False to get a zero-centred result (e.g. a noise estimate that
            will be scaled and subtracted from the original grid).

        Parameters
        ----------
        data : array-like
            Input data to be filtered
        low_cut : float
            Short-wavelength cut-off (wavelengths shorter than this are attenuated)
        high_cut : float
            Long-wavelength cut-off (wavelengths longer than this are attenuated)
        direction_angle : float, optional
            The primary direction to emphasize, in degrees clockwise from north (0-360)
        direction_width : float, optional
            Angular width parameter controlling directional sensitivity (degrees)
        order : int, optional
            Order of the Butterworth filter. Higher orders create sharper transitions
            but may introduce more ringing artifacts. Default is 4.
        buffer_size : int
            Size of the buffer zone for reducing edge effects
        buffer_method : str
            Method for handling the buffer zone ("mirror", "mean", etc.)

        Returns
        -------
        array-like
            The filtered data with both band-pass and directional filtering applied
        """

        def filter_function(kx, ky):
            # Compute wavenumber magnitude for band-pass filter
            k = np.sqrt(kx**2 + ky**2)

            # Convert wavelength cutoffs to wavenumber cutoffs
            k_high = 2 * np.pi / high_cut  # High-pass cutoff (removes long wavelengths)
            k_low = 2 * np.pi / low_cut  # Low-pass cutoff (removes short wavelengths)

            # Avoid division by zero
            k = np.maximum(k, 1e-10)

            # Butterworth high-pass component (attenuates low frequencies)
            high_pass = 1.0 / (1.0 + (k_high / k) ** (2 * order))

            # Butterworth low-pass component (attenuates high frequencies)
            low_pass = 1.0 / (1.0 + (k / k_low) ** (2 * order))

            # Combine to create band-pass filter
            band_pass = high_pass * low_pass

            # Directional filter component
            # Convert direction angle to radians (0 is north, increases clockwise)
            angle_rad = np.radians(
                -direction_angle
            )  # Convert from N=0 to standard math orientation

            # Calculate wavenumber direction
            # Use arctan2 to get angle in all quadrants
            k_angle = np.arctan2(ky, kx)

            # Find smallest angular difference. A grid's spectrum is symmetric
            # (every wavevector has an opposite one), so the direction is an
            # axis, not an arrow: measure the angle modulo 180 degrees, so that
            # k and -k are kept or removed together. (Measured modulo 360 the
            # opposite lobe was removed and the real part of the result lost
            # half its amplitude.)
            angle_diff = np.abs(np.mod(k_angle - angle_rad + np.pi / 2, np.pi) - np.pi / 2)

            # Convert width to radians
            width_rad = np.radians(direction_width)

            # Create directional cosine filter
            # Use cosine-squared function for smooth transitions
            directional_filter = np.cos(angle_diff * np.pi / (2 * width_rad)) ** 2.0

            # Apply cosine taper for angles beyond the specified width
            directional_filter[angle_diff > width_rad] = 0

            # Combine band-pass and directional filters
            combined_filter = band_pass * directional_filter

            return combined_filter

        return self._apply_fourier_filter(
            data, filter_function, buffer_size, buffer_method, preserve_dc=preserve_dc
        )

    def line_noise_removal(
        self,
        data,
        azimuth,
        min_spacing,
        max_spacing=None,
        scale=1.0,
        buffer_size=10,
        direction_width=45,
        order=4,
        return_noise=False,
    ):
        """
        Remove line-parallel (acquisition) noise from a grid.

        Line noise varies across the lines with wavelengths of at least twice
        the line spacing. A zero-centred noise estimate is made by passing
        roughly 2 x min_spacing to 10 x max_spacing wavelengths through a
        directional wedge about `azimuth` (the long-wavelength cut matters:
        regional field also varies across the lines, so it lies inside the
        wedge and would otherwise pass as offset/trend). The estimate,
        multiplied by `scale`, is subtracted from the data.

        Parameters
        ----------
        data : 2D array (north-up, NaN for no data)
        azimuth : direction of the wedge, degrees clockwise from north
        min_spacing, max_spacing : smallest / largest line spacing, in map
            units (max_spacing defaults to min_spacing for constant spacing)
        scale : multiplier applied to the noise estimate before subtraction
        direction_width, order : wedge half-width (degrees) and Butterworth order
        return_noise : return the scaled noise estimate instead of the
            corrected grid

        Returns
        -------
        2D array with NaN where the input was NaN
        """
        min_spacing = float(min_spacing)
        max_spacing = max(float(max_spacing or min_spacing), min_spacing)
        noise = self.directional_butterworth_band_pass(
            data,
            2 * min_spacing,  # low_cut: suppress wavelengths shorter than this
            10 * max_spacing,  # high_cut: suppress wavelengths longer than this
            direction_angle=float(azimuth),
            direction_width=direction_width,
            order=order,
            buffer_size=buffer_size,
            buffer_method="mirror",
            preserve_dc=False,  # zero-centred, so scale doesn't scale the mean
        )
        nan_mask = np.isnan(noise)
        noise[nan_mask] = 0.0
        result = noise * scale if return_noise else data - noise * scale
        result[nan_mask] = np.nan
        return result

    def total_horizontal_gradient(self, data, buffer_size=10, buffer_method="mirror"):
        """
        Compute the total horizontal gradient (THG) of a 2D grid using Fourier filtering.

        Parameters:
            data (numpy.ndarray): Input 2D grid data.
            buffer_size (int): Buffer size for edge handling.
            buffer_method (str): Buffering method ('mirror' or 'zero').

        Returns:
            numpy.ndarray: Total horizontal gradient of the input data.
        """
        dfdx, dfdy = self.compute_derivatives(
            data, "xy", buffer_size=buffer_size, buffer_method=buffer_method
        )
        horizontal_gradient = np.sqrt(dfdx**2 + dfdy**2)
        return horizontal_gradient

    def directional_cosine_filter(self, kx, ky, center_direction, degree):
        """
        Directional Cosine Filter in the Fourier domain.
        """
        theta = np.arctan2(ky, kx)  # Direction of each frequency component
        center_radians = np.deg2rad(
            center_direction
        )  # Convert center direction to radians
        directional_filter = np.abs(np.cos(theta - center_radians)) ** degree
        return directional_filter

    # --- Internal Fourier Filter Application ---

    def _apply_fourier_filter(
        self,
        data,
        filter_function,
        buffer_size=10,
        buffer_method="mirror",
        preserve_dc=True,
        mirror_pad=False,
    ):
        """
        Apply a Fourier-domain filter with buffering and NaN handling.

        Parameters:
            preserve_dc: If True, add median back after filtering (for RTP, continuation, etc.)
                        If False, return zero-centered result (for derivatives)
            mirror_pad: If True, pad by mirroring the edges (tapered within the
                        padding, see _mirror_buffer) rather than ramping to zero
        """
        return self._apply_fourier_filters(
            data, [filter_function], buffer_size, buffer_method, preserve_dc,
            mirror_pad,
        )[0]

    def _apply_fourier_filters(
        self,
        data,
        filter_functions,
        buffer_size=10,
        buffer_method="mirror",
        preserve_dc=True,
        mirror_pad=False,
    ):
        """
        Apply several Fourier-domain filters to the same grid, sharing the NaN
        fill, buffering and forward transform; returns one grid per filter,
        each as _apply_fourier_filter would give it.
        """
        from scipy.fft import fft2, ifft2

        median = np.nanmedian(data)

        # inpaint NaN values
        filled_data, nan_mask = self.fill_nan(data)

        # remove median to reduce ringing artifacts from windowing
        centered_data = filled_data - median

        # Add buffer
        if mirror_pad:
            buffered_data, padding_info = self._mirror_buffer(centered_data, buffer_size)
        else:
            buffered_data, padding_info = self.add_buffer(
                centered_data, buffer_size, buffer_method
            )

        # Fourier transform
        data_fft = fft2(buffered_data, workers=-1)

        ny, nx = buffered_data.shape
        kx, ky = self.create_wavenumber_grids(nx, ny)

        results = []
        for filter_function in filter_functions:
            # Compute filter
            filter_array = filter_function(kx, ky)

            # Apply filter
            filtered_fft = data_fft * filter_array

            # Inverse transform
            filtered_data = np.real(ifft2(filtered_fft, workers=-1))

            # Remove buffer
            unbuffered_data = self.remove_buffer(
                filtered_data, buffer_size, padding_info
            )

            # Restore DC level only for non-derivative operations
            if preserve_dc:
                unbuffered_data = unbuffered_data + median

            results.append(self.restore_nan(unbuffered_data, nan_mask))
        return results

    def bsdwormer(
        self, image, layer, gridPath, num_levels, bottom_level, delta_z, shps, crs,
        callback=None, out_dir=None,
    ):
        # code borrows heavilly from bsdwormer example ipynb template example
        # adds on the fly calc of padded grid
        # callback(fraction_done) is called before each level and may raise
        # calcs.sgt_cancel.OperationCancelled to stop the calculation.
        # layer: a QgsRasterLayer or LayerShim (extent and EPSG code), crs: the
        # EPSG number; outputs go beside gridPath, or into out_dir if given.

        print(gridPath, num_levels, bottom_level, delta_z)
        # load grid and convert to numpy

        # add padding and save to file

        # Remove the mean
        image_zero_mean = image - np.nanmean(image)
        image_zero_mea_nonan = fill_nan(image_zero_mean)
        pad_y = int(0.2 * image.shape[0])
        pad_x = int(0.2 * image.shape[1])

        pad_x = ceil(pad_x / 2) * 2
        pad_y = ceil(pad_y / 2) * 2
        padded_image = np.pad(
            image_zero_mea_nonan,
            pad_width=((pad_y, pad_y), (pad_x, pad_x)),
            mode="linear_ramp",
            # constant_values=0,
        )
        print("padded_size=", padded_image.shape)

        # Create a 2D Hanning window
        window_y = np.hanning(padded_image.shape[0])
        window_x = np.hanning(padded_image.shape[1])
        hanning_window = np.outer(window_y, window_x)

        # Apply the Hanning window to the padded image
        final_image = padded_image * hanning_window
        # save padded image to file as

        job = Wormer()
        job.importGdalRaster(gridPath)
        extent = job.geomat
        print("job.geomat", job.geomat)
        # extentPadded = GetExtent(job.geomat, -pad_x, -pad_y)
        xmin_padded, ymax_padded = GetExtent(job.geomat, pad_x, pad_y)
        xmin_padded = job.geomat[0] - (pad_x * job.geomat[1])
        ymax_padded = job.geomat[3] - (pad_y * job.geomat[5])
        extentPadded = [
            [None, None],
            [None, None],
            [xmin_padded, ymax_padded],
            [None, None],
        ]
        print("geomat[5]:", job.geomat[5])

        gridPathPadded = insert_text_before_extension(gridPath, "_padded")

        # Separate the file path into directory, base name, and extension
        dir_name, base_name = os.path.split(gridPath)
        if out_dir:
            dir_name = out_dir
            gridPathPadded = os.path.join(out_dir, os.path.basename(gridPathPadded))
        file_name, file_ext = os.path.splitext(base_name)
        wormsPath = dir_name + "/" + file_name + "_worms.csv"

        numpy_array_to_raster(
            final_image,
            gridPathPadded,
            pad_x,
            pad_y,
            dx=job.geomat[1],
            xmin=extentPadded[2][0],
            ymax=extentPadded[2][1],
            reference_layer=layer,
            no_data_value=np.nan,
        )

        job.importExternallyPaddedRaster(gridPathPadded, pad_x, pad_y)

        # Initialize a flag to track if the header has been written
        header_written = False
        print("num_levels", num_levels)

        # the levels take most of the time, unless polylines are made from them too
        level_share = 0.7 if shps else 1.0
        for dz in range(0, num_levels):
            if callback is not None:
                callback(level_share * dz / float(num_levels))
            dzm = (dz * delta_z) + bottom_level
            if dzm == 0:
                dzm = 0.01
            job.wormLevelAsPoints(dz=dzm)

            z = [dzm] * len(job.worm_ys)
            x = (job.worm_xs * job.geomat[1]) + extentPadded[2][0]
            y = (job.worm_ys * job.geomat[5]) + extentPadded[2][1]
            vals = job.worm_vals

            # Combine into a single array
            points = np.column_stack(
                (
                    y[: len(job.worm_ys)],
                    x[: len(job.worm_ys)],
                    z[: len(job.worm_ys)],
                    vals[: len(job.worm_ys)],
                )
            )
            extent2 = [
                min(extent[0], extent[0] + (image.shape[1] * job.geomat[1])),
                max(extent[0], extent[0] + (image.shape[1] * job.geomat[1])),
                min(extent[3], extent[3] + (image.shape[0] * job.geomat[5])),
                max(extent[3], extent[3] + (image.shape[0] * job.geomat[5])),
            ]

            # Apply the filter conditions
            mask = (
                (points[:, 1] >= extent2[0])
                & (points[:, 1] <= extent2[1])  # extent[1]>= y >= extent[0]
                & (points[:, 0] >= extent2[2])
                & (points[:, 0] <= extent2[3])  # extent[3]>= x >= extent[4]
            )

            filtered_points = points[mask]

            if os.path.exists(wormsPath) and not header_written:
                os.remove(wormsPath)

            # Write to file
            with open(wormsPath, "a") as f:
                if not header_written:
                    # Write the header if it hasn't been written yet
                    f.write("y,x,z,val\n")
                    header_written = True
                # Append the filtered points
                write_rows(f, filtered_points, "%s")  # same text as savetxt fmt="%s"

        if shps:
            max_distance = job.geomat[1] * 1.3
            in_path = wormsPath
            out_path = wormsPath.replace(".csv", ".shp")
            polyline_callback = None
            if callback is not None:
                def polyline_callback(fraction=None):
                    callback(None if fraction is None
                             else level_share + (1 - level_share) * fraction)
            self.xyz_to_polylines(in_path, out_path, max_distance, crs,
                                  callback=polyline_callback)

    def display_grid(self, grid):
        import matplotlib.pyplot as plt

        # Plot results
        plt.imshow(grid, cmap="gray")
        plt.colorbar()
        plt.title("Grid")
        plt.show()

    def process_cluster(self, points, z, max_distance, callback=None):
        """Process a cluster of points to create split LineStrings with better connectivity.
        callback(None), if given, is called every 2000 points and may raise
        OperationCancelled to stop."""
        from scipy.spatial import cKDTree
        from shapely.geometry import LineString, Point
        import numpy as np

        if len(points) <= 1:
            return []

        # Create KD-tree for efficient nearest neighbor search
        tree = cKDTree(points)

        # Track visited points
        visited = np.zeros(len(points), dtype=bool)

        # Start with first point
        sorted_points = [points[0]]
        visited[0] = True
        remaining_count = len(points) - 1

        # Find path through points
        while remaining_count > 0:
            if callback is not None and remaining_count % 2000 == 0:
                callback(None)  # a big cluster can take seconds: let Cancel stop it
            last_point = sorted_points[-1]

            # Query more points than needed to have alternatives if the closest are already visited
            # Use a higher k value to have more candidates
            k_value = min(len(points), 10)  # Try up to 10 nearest neighbors
            distances, indices = tree.query(last_point, k=k_value)

            # Handle scalar case
            if np.isscalar(distances):
                distances = [distances]
                indices = [indices]

            # Find the closest unvisited point within max_distance
            found_next = False
            for idx, dist in zip(indices, distances):
                if not visited[idx] and dist <= max_distance:
                    sorted_points.append(points[idx])
                    visited[idx] = True
                    remaining_count -= 1
                    found_next = True
                    break

            # If no point within max_distance is found, try to start a new segment
            if not found_next:
                # Find the closest unvisited point to any point in our current
                # path: the first path point with the smallest distance to an
                # unvisited point (one query per path point against a tree of
                # the unvisited points only), then its nearest unvisited point
                # in the order the full tree lists them
                best_pair = None
                unvisited = np.flatnonzero(~visited)
                if unvisited.size:
                    nearest, _ = cKDTree(points[unvisited]).query(sorted_points, k=1)
                    i = int(np.argmin(nearest))
                    distances, indices = tree.query(sorted_points[i], k=len(points))
                    for idx in np.atleast_1d(indices):
                        if not visited[idx]:
                            best_pair = (i, idx)
                            break

                # If we found an unvisited point, start a new segment
                if best_pair:
                    sorted_points.append(points[best_pair[1]])
                    visited[best_pair[1]] = True
                    remaining_count -= 1
                else:
                    # No more points can be connected
                    break

        # Create and split the polyline
        if len(sorted_points) > 1:
            polyline = LineString(sorted_points)
            split_segments = self.split_long_segments(polyline, max_distance)
            return [(seg, z) for seg in split_segments]

        return []

    def split_long_segments(self, line, max_distance):
        from shapely.geometry import Point, LineString

        """Split a LineString into valid segments and discard only the segments longer than max_distance."""
        coords = list(line.coords)
        new_lines = []
        segment = [coords[0]]

        for i in range(1, len(coords)):
            segment.append(coords[i])
            segment_length = Point(coords[i - 1]).distance(Point(coords[i]))

            if segment_length > max_distance:
                if len(segment) > 2:
                    new_lines.append(LineString(segment[:-1]))
                segment = [coords[i]]

        if len(segment) > 1:
            new_lines.append(LineString(segment))

        return new_lines

    def xyz_to_polylines(self, csv_file, output_shapefile, max_distance, crs, callback=None):
        """
        Converts XYZ points from a CSV file into polylines and saves them as a shapefile.
        Parameters:
        csv_file (str): Path to the input CSV file containing XYZ points.
        output_shapefile (str): Path to the output shapefile where polylines will be saved.
        max_distance (float): Maximum distance between points to be considered part of the same cluster.
        crs (int): Coordinate Reference System (CRS) EPSG code for the output shapefile.
        callback (callable, optional): callback(fraction done), called for each
            cluster (and callback(None) within big ones); it may raise
            OperationCancelled to stop.
        Raises:
        ValueError: If no valid LineStrings are created.
        Notes:
        - The CSV file should have columns named "x", "y", and "z".
        - Points with a Z value less than 1.0 are ignored.
        - Uses DBSCAN clustering to group points into clusters based on the max_distance.
        - Each cluster is processed into a polyline if it contains more than one point.
        - The resulting polylines are saved in a shapefile with an attribute field "Z" indicating the Z value.
        Example:
        >>> processor = GeophysicalProcessor()
        >>> processor.xyz_to_polylines("input.csv", "output.shp", 10.0, 4326)
        """
        from sklearn.cluster import DBSCAN
        from shapely.geometry import Point, LineString

        with open(csv_file, "r") as f:
            reader = csv.DictReader(f)
            data = [
                (float(row["x"]), float(row["y"]), float(row["z"])) for row in reader
            ]

        data = np.array(data)
        unique_z = np.unique(data[:, 2])
        results = []

        for zi, z in enumerate(unique_z):

            mask = data[:, 2] == z
            points = data[mask][:, :2]
            if len(points) < 2 or z < 1.0:
                continue

            print(f"Processing Z-bin: {z}, Number of points: {len(points)}")

            try:
                clustering = DBSCAN(eps=max_distance, min_samples=1).fit(points)
            except Exception as e:
                print(f"Error in DBSCAN at Z-bin {z}: {e}")
                continue

            labels = clustering.labels_
            unique_labels = np.unique(labels)
            unique_labels = unique_labels[unique_labels != -1]

            for ci, label in enumerate(unique_labels):
                if callback is not None:
                    callback((zi + ci / len(unique_labels)) / len(unique_z))
                cluster_points = points[labels == label]
                if len(cluster_points) > 1:
                    try:
                        result = self.process_cluster(
                            cluster_points, z, max_distance, callback=callback)
                        if result:
                            results.extend(result)
                    except OperationCancelled:
                        raise
                    except Exception as e:
                        print(f"Error processing cluster at Z-bin {z}: {e}")

        # Collect all polylines and corresponding Z values
        polylines = [item[0] for item in results if isinstance(item[0], LineString)]
        z_values = [item[1] for item in results if isinstance(item[0], LineString)]

        if not polylines:
            raise ValueError(
                "No valid LineStrings were created. Check your input data or parameters."
            )

        print("Processing complete. Saving to shapefile...")

        # Save all results in a single shapefile
        driver = ogr.GetDriverByName("ESRI Shapefile")
        ds = driver.CreateDataSource(output_shapefile)
        spatial_ref = osr.SpatialReference()
        spatial_ref.ImportFromEPSG(crs)
        # 3D (PolylineZ): every vertex carries the worm level as its Z
        layer = ds.CreateLayer("polylines", spatial_ref, ogr.wkbLineString25D)

        field_defn = ogr.FieldDefn("Z", ogr.OFTInteger)
        layer.CreateField(field_defn)

        for polyline, z in zip(polylines, z_values):
            feature = ogr.Feature(layer.GetLayerDefn())
            feature.SetField("Z", int(z))
            geom = ogr.Geometry(ogr.wkbLineString25D)
            for x, y in polyline.coords:
                geom.AddPoint(x, y, float(z))
            feature.SetGeometry(geom)
            layer.CreateFeature(feature)
            feature = None

        ds = None  # Close and save the shapefile
        print(f"Shapefile saved to {output_shapefile}")

    def remove_2o_gradient(self, data, mask):
        """Fit a second-order polynomial surface to the raster data and remove it."""
        from scipy.optimize import leastsq

        rows, cols = data.shape
        X, Y = np.meshgrid(np.arange(cols), np.arange(rows))

        x, y, z = X[~mask], Y[~mask], data[~mask]
        if len(z) == 0:
            print("Warning: No valid data for gradient removal.")
            return data

        def poly2d(params, x, y):
            a, b, c, d, e, f = params
            return a * x**2 + b * y**2 + c * x * y + d * x + e * y + f

        def residuals(params, x, y, z):
            return poly2d(params, x, y) - z

        params_init = [0, 0, 0, 0, 0, np.mean(z)]
        params_opt, _ = leastsq(residuals, params_init, args=(x, y, z))

        fitted_surface = poly2d(params_opt, X, Y)
        return data - fitted_surface

    def remove_gradient(self, data, mask):
        """Fit a plane to the raster data and remove the first-order gradient."""
        from scipy.optimize import leastsq

        rows, cols = data.shape
        X, Y = np.meshgrid(np.arange(cols), np.arange(rows))

        x, y, z = X[~mask], Y[~mask], data[~mask]
        if len(z) == 0:
            print("Warning: No valid data for gradient removal.")
            return data

        def plane(params, x, y):
            a, b, c = params
            return a * x + b * y + c

        def residuals(params, x, y, z):
            return plane(params, x, y) - z

        params_init = [0, 0, np.mean(z)]
        params_opt, _ = leastsq(residuals, params_init, args=(x, y, z))

        fitted_plane = plane(params_opt, X, Y)
        return data - fitted_plane

    def fix_extreme_values(self, data):
        """Replace extreme NoData values with a safe float32-compatible value."""
        safe_nodata_value = -999999.0
        data[data < -1e38] = safe_nodata_value
        return data, safe_nodata_value

    def normalise_geotiffs(self, input_folder, output_folder, order, callback=None):
        """
        Process multiple GeoTIFF files by normalizing them using parameters from the first GeoTIFF.
        Args:
            input_folder (str): Path to the folder containing input GeoTIFF files.
            output_folder (str): Path to the folder where normalized GeoTIFF files will be saved.
            order (bool): If True, remove first-order gradient; if False, remove second-order gradient.
        Returns:
            None
        The function performs the following steps:
            1. Reads all GeoTIFF files from the input folder.
            2. For each GeoTIFF file:
                a. Opens the file and reads the raster data.
                b. Fixes extreme values in the data.
                c. Removes the gradient (first-order or second-order) from the data.
                d. Normalizes the data using the standard deviation of the valid data from the first GeoTIFF.
                e. Saves the normalized data to a new GeoTIFF file in the output folder.
        """
        tiff_files = sorted(glob.glob(os.path.join(input_folder, "*.tif")))
        if not tiff_files:
            print("No GeoTIFF files found.")
            return

        first = True
        for file_index, tiff_file in enumerate(tiff_files):
            if callback is not None:
                # may raise calcs.sgt_cancel.OperationCancelled to stop
                callback(file_index / float(len(tiff_files)))
            # print("tiff_file", tiff_file)
            ds = gdal.Open(tiff_file, gdal.GA_ReadOnly)
            band = ds.GetRasterBand(1)
            data = band.ReadAsArray()
            geotransform = ds.GetGeoTransform()
            projection = ds.GetProjection()

            data, nodata_value = self.fix_extreme_values(data)
            data = data.astype(np.float32)
            mask = (data == nodata_value) | np.isnan(data)

            if order:
                detrended_data = self.remove_gradient(data, mask)
            else:
                detrended_data = self.remove_2o_gradient(data, mask)

            valid_data = detrended_data[~mask]

            # Check for empty valid_data
            if valid_data.size == 0:
                print("Warning: No valid data found in", tiff_file)
                continue

            # Check for zero standard deviation
            valid_std = np.nanstd(valid_data)
            if valid_std == 0:
                print(
                    "Warning: Standard deviation is zero, skipping normalization for",
                    tiff_file,
                )
                continue

            if first:
                first = False
                reference_std = valid_std  # Use the non-zero std

            normalized_data = (
                reference_std * (detrended_data - np.nanmean(valid_data)) / valid_std
            )

            normalized_data[mask] = nodata_value

            driver = gdal.GetDriverByName("GTiff")
            output_path = os.path.join(output_folder, os.path.basename(tiff_file))
            out_ds = driver.Create(
                output_path, ds.RasterXSize, ds.RasterYSize, 1, gdal.GDT_Float32
            )
            out_ds.SetGeoTransform(geotransform)
            out_ds.SetProjection(projection)
            out_band = out_ds.GetRasterBand(1)
            out_band.WriteArray(normalized_data)
            out_band.SetNoDataValue(nodata_value)
            out_band.FlushCache()

            print(f"Saved: {output_path}")


if __name__ == "__main__":
    # import matplotlib.pyplot as plt

    # Create synthetic data
    nx, ny = 100, 100
    dx, dy = 1, 1
    x = np.linspace(0, 99, nx)
    y = np.linspace(0, 99, ny)
    X, Y = np.meshgrid(x, y)
    data = (
        0.01 * X**2 + 0.02 * Y + np.sin(2 * np.pi * X / 10) + np.cos(2 * np.pi * Y / 15)
    )

    # Add NaNs
    data[40:60, 40:60] = np.nan

    # Initialize processor
    processor = GeophysicalProcessor(dx, dy)

    # Apply upward continuation
    result = processor.upward_continuation(data, height=10)

    # Plot results
    plt.imshow(result, cmap="viridis")
    plt.colorbar()
    plt.title("Upward Continuation")
    plt.show()
    result = processor.tilt_derivative(data)

    plt.imshow(result, cmap="viridis")
    plt.colorbar()
    plt.title("tilt_angle")
    plt.show()

    # Define scattered data points
    x = np.array([0, 10, 20, 30, 40])
    y = np.array([0, 10, 20, 30, 40])
    z = np.array([1, 3, 7, 13, 21])  # Values at data points

    # Define grid
    grid_x, grid_y = np.meshgrid(np.linspace(0, 40, 50), np.linspace(0, 40, 50))

    # Remove polynomial trend
    residual, trend = processor.remove_polynomial_trend(
        x, y, z, grid_x, grid_y, degree=2
    )

    # Plot results
    plt.imshow(trend, cmap="viridis", origin="lower")
    plt.colorbar()
    plt.title("Polynomial Trend")
    plt.show()

    # Remove regional trend using Fourier (low-pass filter)
    regional_removed = processor.remove_regional_trend_fourier(
        data, cutoff_wavelength=30
    )

    # Plot results
    plt.imshow(regional_removed, cmap="viridis", origin="lower")
    plt.colorbar()
    plt.title("Regional Trend Removed (Fourier)")
    plt.show()

    # Apply upward continuation
    upward_result = processor.upward_continuation(data, height=10)

    # Plot results
    plt.imshow(upward_result, cmap="viridis", origin="lower")
    plt.colorbar()
    plt.title("Upward Continuation")
    plt.show()

    # Apply downward continuation
    downward_result = processor.downward_continuation(data, height=10)

    # Plot results
    plt.imshow(downward_result, cmap="viridis", origin="lower")
    plt.colorbar()
    plt.title("Downward Continuation")
    plt.show()

    # Perform Reduction to Pole
    rtp_result = processor.reduction_to_pole(data, inclination=60, declination=30)

    # Plot results
    plt.imshow(rtp_result, cmap="viridis", origin="lower")
    plt.colorbar()
    plt.title("Reduction to Pole")
    plt.show()

    # Perform Reduction to Equator
    rte_result = processor.reduction_to_equator(data, inclination=10, declination=30)

    # Plot results
    plt.imshow(rte_result, cmap="viridis", origin="lower")
    plt.colorbar()
    plt.title("Reduction to Equator")
    plt.show()

    # Compute first horizontal derivative in x-direction
    dx_result = processor.compute_derivative(data, direction="x", order=1)

    # Plot results
    plt.imshow(dx_result, cmap="seismic", origin="lower")
    plt.colorbar()
    plt.title("Horizontal Derivative (X-Direction)")
    plt.show()

    # Compute first vertical derivative
    dz_result = processor.compute_derivative(data, direction="z", order=1)

    # Plot results
    plt.imshow(dz_result, cmap="seismic", origin="lower")
    plt.colorbar()
    plt.title("Vertical Derivative")
    plt.show()

    # Compute tilt derivative
    tilt_result = processor.tilt_angle(data)

    # Plot results
    plt.imshow(
        tilt_result, cmap="seismic", origin="lower", vmin=-np.pi / 2, vmax=np.pi / 2
    )
    plt.colorbar()
    plt.title("Tilt Angle")
    plt.show()

    # Compute analytic signal
    analytic_signal_result = processor.analytic_signal(data)

    # Plot results
    plt.imshow(analytic_signal_result, cmap="viridis", origin="lower")
    plt.colorbar()
    plt.title("Analytic Signal")
    plt.show()

    # Apply Automatic Gain Control
    agc_result = processor.automatic_gain_control(data, window_size=5)

    # Plot results
    plt.imshow(agc_result, cmap="viridis", origin="lower")
    plt.colorbar()
    plt.title("Automatic Gain Control (AGC)")
    plt.show()

    # Apply low-pass filter
    low_pass_result = processor.low_pass_filter(data, cutoff_wavelength=20)

    # Plot results
    plt.imshow(low_pass_result, cmap="viridis", origin="lower")
    plt.colorbar()
    plt.title("Low-Pass Filter")
    plt.show()

    # Apply high-pass filter
    high_pass_result = processor.high_pass_filter(data, cutoff_wavelength=20)

    # Plot results
    plt.imshow(high_pass_result, cmap="viridis", origin="lower")
    plt.colorbar()
    plt.title("High-Pass Filter")
    plt.show()

    # Apply band-pass filter
    band_pass_result = processor.band_pass_filter(data, low_cut=10, high_cut=50)

    # Plot results
    plt.imshow(band_pass_result, cmap="viridis", origin="lower")
    plt.colorbar()
    plt.title("Band-Pass Filter")
    plt.show()
