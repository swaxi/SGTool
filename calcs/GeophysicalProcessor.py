import numpy as np
from numpy.polynomial.polynomial import polyval2d
#from scipy.spatial import cKDTree
#from scipy.fftpack import fft2, ifft2, fftfreq
#from scipy.ndimage import gaussian_filter
#from shapely.geometry import LineString, Point
from math import ceil
import numpy as np


from osgeo import ogr, osr

# from joblib import Parallel, delayed
import csv

import glob
import os
from osgeo import gdal
#from scipy.optimize import leastsq

from ..worms.wormer import Wormer
from ..worms.Utility import (
    GetExtent,
    loadGrid,
    numpy_array_to_raster,
    insert_text_before_extension,
    fill_nan,
)
import os


class GeophysicalProcessor:
    def __init__(self, dx, dy, buffer):
        """
        Initialize the processor with grid spacing.

        Parameters:
            dx (float): Grid spacing in the x-direction.
            dy (float): Grid spacing in the y-direction.
        """
        from scipy.spatial import cKDTree
        from scipy.fftpack import fft2, ifft2, fftfreq
        from scipy.ndimage import gaussian_filter
        from shapely.geometry import LineString, Point
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
        return np.meshgrid(kx, ky)

    def fill_nanyy(self, data):
        """
        Replace NaN values with interpolated values from non-NaN neighbors.
        """
        # Make a copy of the data
        filled_data = np.copy(data)
        nan_mask = np.isnan(data)

        # Get indices of all non-NaN points
        non_nan_indices = np.where(~nan_mask)
        non_nan_values = data[~nan_mask]

        # Get indices of all NaN points
        nan_indices = np.where(nan_mask)

        # Only interpolate if we have some non-NaN values and some NaN values
        if len(non_nan_values) > 0 and np.any(nan_mask):
            # Use griddata for interpolation
            from scipy.interpolate import griddata

            # Interpolate using the non-NaN values
            filled_values = griddata(
                np.array(non_nan_indices).T,  # Points we know
                non_nan_values,  # Values we know
                np.array(nan_indices).T,  # Points to interpolate
                method="linear",  # Linear interpolation
                fill_value=np.nanmedian(
                    data
                ),  # Use median for points outside the convex hull
            )

            # Place the interpolated values back into the array
            filled_data[nan_mask] = filled_values

        return filled_data, nan_mask

    def fill_nan_slow(self, data):
        """
        Replace NaN values with interpolated values from non-NaN neighbors.
        """
        # Make a copy of the data
        filled_data = np.copy(data)
        nan_mask = np.isnan(data)

        # Get indices of all non-NaN points
        non_nan_indices = np.where(~nan_mask)
        non_nan_values = data[~nan_mask]

        # Get indices of all NaN points
        nan_indices = np.where(nan_mask)

        # Only interpolate if we have some non-NaN values and some NaN values
        if len(non_nan_values) > 0 and np.any(nan_mask):
            # Use griddata for interpolation
            from scipy.interpolate import griddata

            # Interpolate using the non-NaN values
            filled_values = griddata(
                np.array(non_nan_indices).T,  # Points we know
                non_nan_values,  # Values we know
                np.array(nan_indices).T,  # Points to interpolate
                method="linear",  # Linear interpolation
                fill_value=np.nanmedian(
                    data
                ),  # Use median for points outside the convex hull
            )

            # Place the interpolated values back into the array
            filled_data[nan_mask] = filled_values

        return filled_data, nan_mask

    def fill_nan(self, data):
        """
        Replace NaN values with interpolated values from non-NaN neighbors
        using a faster approach that works in both x and y directions.
        """
        from scipy import ndimage

        # Create a copy of the data
        filled_data = np.copy(data)
        nan_mask = np.isnan(data)

        # If there are no NaNs or all NaNs, handle accordingly
        if not np.any(nan_mask) or np.all(nan_mask):
            return filled_data, nan_mask

        # First, try a more direct interpolation along both axes
        # Iterate through each row (x-direction)
        for i in range(data.shape[0]):
            row = filled_data[i, :]
            mask = ~np.isnan(row)
            indices = np.arange(len(row))
            if np.any(mask) and not np.all(
                mask
            ):  # If row has some valid values but not all
                # Linear interpolation along x-axis
                filled_data[i, ~mask] = np.interp(
                    indices[~mask], indices[mask], row[mask]
                )

        # Iterate through each column (y-direction)
        for j in range(data.shape[1]):
            col = filled_data[:, j]
            mask = ~np.isnan(col)
            indices = np.arange(len(col))
            if np.any(mask) and not np.all(
                mask
            ):  # If column has some valid values but not all
                # Linear interpolation along y-axis
                filled_data[~mask, j] = np.interp(
                    indices[~mask], indices[mask], col[mask]
                )

        # Some points might still be NaN if entire rows or columns were NaN
        # For remaining NaNs, use the smoother approach
        remaining_nan_mask = np.isnan(filled_data)
        if np.any(remaining_nan_mask):
            # Replace remaining NaNs with a placeholder value
            filled_data[remaining_nan_mask] = 0

            # Create a mask of valid values (non-NaN)
            valid_mask = ~remaining_nan_mask

            # Use a distance-weighted interpolation approach
            # First, generate a distance transform from valid points
            dist = ndimage.distance_transform_edt(~valid_mask)

            # Use a Gaussian filter for interpolation
            sigma = 2.0  # Increased for better coverage

            # Apply weight mask to valid data points
            weights = np.exp(-(dist**2) / (2 * sigma**2))
            weights[~valid_mask] = 0

            # Normalize weights to sum to 1
            weight_sum = ndimage.gaussian_filter(weights, sigma)
            # Avoid division by zero
            weight_sum[weight_sum == 0] = 1

            # Calculate the weighted average
            weighted_data = filled_data * weights
            weighted_avg = ndimage.gaussian_filter(weighted_data, sigma) / weight_sum

            # Only replace remaining NaN values with the interpolated values
            filled_data[remaining_nan_mask] = weighted_avg[remaining_nan_mask]

        return filled_data, nan_mask

    def restore_nan(self, data, nan_mask):
        """
        Restore NaN values to their original positions in the grid.
        """
        data_with_nan = np.copy(data)
        data_with_nan[nan_mask] = np.nan
        return data_with_nan

    def add_buffer(self, data, buffer_size, method="mirror"):
        """
        Add a buffer around the edges to reduce edge effects.
        """
        self.data_mean = np.nanmean(data)
        # data_zero = data - self.data_mean

        # Determine padding sizes (20% of original dimensions)
        pad_y = buffer_size
        pad_x = buffer_size

        # Pad the image with zeros
        padded_image = np.pad(
            data,
            pad_width=((pad_y, pad_y), (pad_x, pad_x)),
            mode="linear_ramp",
            # constant_values=0
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
        return padded_image

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

    def remove_buffer(self, data, buffer_size):
        """
        Remove the buffer from the edges of the grid.
        """
        data = data  # + self.data_mean
        return data[buffer_size:-buffer_size, buffer_size:-buffer_size]

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
        dx = self._apply_fourier_filter(
            data, filter_function_dx, buffer_size, buffer_method
        )
        dy = self._apply_fourier_filter(
            data, filter_function_dy, buffer_size, buffer_method
        )

        # Compute the total horizontal gradient
        thg = np.sqrt(dx**2 + dy**2)

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
        # Convert angles from degrees to radians
        inc, dec = np.radians(inclination), np.radians(declination)

        def filter_function(kx, ky):
            # Transform to radians
            inc_rad = np.radians(inclination)
            dec_rad = np.radians(-declination)
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

    def rte_filter_function(
        kx,
        ky,
        inclination,
        declination,
        magnetization_inclination=None,
        magnetization_declination=None,
    ):
        """
        Generate a Reduction to the Equator (RTE) filter in the frequency domain.

        Parameters
        ----------
        kx : array
            Wavenumber component in x-direction.
        ky : array
            Wavenumber component in y-direction.
        inclination : float
            Inclination of the inducing geomagnetic field in degrees.
        declination : float
            Declination of the inducing geomagnetic field in degrees.
        magnetization_inclination : float, optional
            Inclination of the total magnetization of the source in degrees.
        magnetization_declination : float, optional
            Declination of the total magnetization of the source in degrees.

        Returns
        -------
        numpy.ndarray
            Fourier-domain filter for reduction to the equator.
        """
        # Convert angles to radians
        inc_rad = np.radians(inclination)
        dec_rad = np.radians(declination)

        if magnetization_inclination is None:
            magnetization_inclination = inclination
        if magnetization_declination is None:
            magnetization_declination = declination

        mag_inc_rad = np.radians(magnetization_inclination)
        mag_dec_rad = np.radians(magnetization_declination)

        # Compute wavevector magnitude and direction
        k_squared = kx**2 + ky**2
        k = np.sqrt(k_squared)
        theta_k = np.arctan2(ky, kx)  # Angle of wavevector

        # Avoid division by zero
        epsilon = 1e-10
        k[k == 0] = epsilon

        # Compute field components
        F_x = np.cos(inc_rad) * np.cos(dec_rad)
        F_y = np.cos(inc_rad) * np.sin(dec_rad)
        F_z = np.sin(inc_rad)

        M_x = np.cos(mag_inc_rad) * np.cos(mag_dec_rad)
        M_y = np.cos(mag_inc_rad) * np.sin(mag_dec_rad)
        M_z = np.sin(mag_inc_rad)

        # Compute denominator (instability occurs near zero inclination)
        denominator = F_z * k + 1j * (F_x * kx + F_y * ky)
        numerator = M_z * k + 1j * (M_x * kx + M_y * ky)

        # Apply a taper for numerical stability at low inclinations
        incl_factor = np.abs(np.cos(inc_rad))  # Stability factor
        taper = (
            np.exp(-((1 - incl_factor) * (k / np.max(k)) ** 2))
            if incl_factor < 0.2
            else 1.0
        )

        # Compute RTE filter
        rte_filter = np.zeros_like(k, dtype=complex)
        mask = np.abs(denominator) > epsilon
        rte_filter[mask] = (
            k_squared[mask] / (denominator[mask] * numerator[mask])
        ) * taper

        # Set zero frequency to zero to avoid division instability
        rte_filter[k == 0] = 0
        rte_filter[~np.isfinite(rte_filter)] = 0  # Remove NaNs and infinities
        rte_filter[np.abs(rte_filter) > 1e6] = 0  # Cap large values for stability

        return rte_filter

        return self._apply_fourier_filter(
            data, filter_function, buffer_size, buffer_method
        )

    def reduction_to_equator(
        self, data, inclination, declination, buffer_size=10, buffer_method="mirror"
    ):
        from scipy.fftpack import fft2, ifft2, fftfreq
        # Convert angles from degrees to radians
        inc, dec = np.radians(inclination), np.radians(declination)

        # Handle NaN values
        filled_data, nan_mask = self.fill_nan(data)

        # Add buffer
        buffered_data = self.add_buffer(filled_data, buffer_size, buffer_method)

        # Calculate frequency coordinates
        ny, nx = buffered_data.shape
        kx = fftfreq(nx, self.dx) * 2 * np.pi
        ky = fftfreq(ny, self.dy) * 2 * np.pi
        kx, ky = np.meshgrid(kx, ky, indexing="ij")
        k = np.sqrt(kx**2 + ky**2) + 1e-10  # Avoid division by zero

        # Directional cosines
        cos_inc = np.cos(inc)
        sin_inc = np.sin(inc)
        cos_dec = np.cos(dec)
        sin_dec = np.sin(dec)

        # RTP filter
        # rte_filter = theta_f
        with np.errstate(divide="ignore", invalid="ignore"):
            rte_filter = (
                k * cos_inc * cos_dec + 1j * ky * cos_inc * sin_dec + kx * sin_inc
            ) / (k * cos_inc * cos_dec - 1j * ky * cos_inc * sin_dec + kx * sin_inc)

        rte_filter[np.abs(rte_filter) > 1e6] = 0  # Stabilize extreme values

        # Apply max_wavenumber filtering based on cell size
        # rte_filter[k > 1/(float(self.dx)*2)] = 0

        # Apply RTP filter in the Fourier domain
        data_fft = fft2(buffered_data)
        data_rte_fft = data_fft * rte_filter.T

        data_rte = np.real(ifft2(data_rte_fft))

        # Remove buffer and restore NaNs
        filtered_data = self.remove_buffer(data_rte, buffer_size)
        return self.restore_nan(filtered_data, nan_mask)

    # --- Derivatives ---
    def compute_derivative(
        self, data, direction, order=1, buffer_size=10, buffer_method="mirror"
    ):
        """
        Compute derivatives in x, y, or z directions.
        """

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

        return self._apply_fourier_filter(
            data, filter_function, buffer_size, buffer_method
        )

    def tilt_angle(self, data, buffer_size=10, buffer_method="mirror"):
        """
        Compute the tilt derivative of the data.
        """
        dfdx = self.compute_derivative(
            data, "x", buffer_size=buffer_size, buffer_method=buffer_method
        )
        dfdy = self.compute_derivative(
            data, "y", buffer_size=buffer_size, buffer_method=buffer_method
        )
        horizontal_gradient = np.sqrt(dfdx**2 + dfdy**2)
        dfdz = self.compute_derivative(
            data, "z", buffer_size=buffer_size, buffer_method=buffer_method
        )
        return np.arctan2(dfdz, horizontal_gradient)

    def analytic_signal(self, data, buffer_size=10, buffer_method="mirror"):
        """
        Compute the analytic signal amplitude.
        """
        dfdx = self.compute_derivative(
            data, "x", buffer_size=buffer_size, buffer_method=buffer_method
        )
        dfdy = self.compute_derivative(
            data, "y", buffer_size=buffer_size, buffer_method=buffer_method
        )
        dfdz = self.compute_derivative(
            data, "z", buffer_size=buffer_size, buffer_method=buffer_method
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

    def band_pass_filter(
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
        dfdx = self.compute_derivative(
            data, "x", buffer_size=buffer_size, buffer_method=buffer_method
        )
        dfdy = self.compute_derivative(
            data, "y", buffer_size=buffer_size, buffer_method=buffer_method
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

    def combined_BHP_DirCos_filter(
        self, data, cutoff_wavelength, center_direction, degree, buffer_size
    ):
        """
        Apply the Directional Cosine filter to the data.
        """

        def filter_function_dc(kx, ky):
            # Create the Direction Cosine high-pass filter
            directional = self.directional_cosine_filter(
                kx, ky, center_direction, degree
            )

            return directional

        return self._apply_fourier_filter(
            data, filter_function_dc, buffer_size=buffer_size, buffer_method="mirror"
        )

    # --- Internal Fourier Filter Application ---
    def _apply_fourier_filter(
        self, data, filter_function, buffer_size=10, buffer_method="mirror"
    ):
        """
        Apply a Fourier-domain filter with buffering and NaN handling.
        """
        from scipy.fftpack import fft2, ifft2, fftfreq
        # inpaint NaN values
        # Handle NaN values
        filled_data, nan_mask = self.fill_nan(data)
        # Add buffer
        buffered_data = self.add_buffer(filled_data, buffer_size, buffer_method)
        # self.display_grid(buffered_data)

        # buffered_data=self.apply_window(buffered_data)
        # Fourier transform
        data_fft = fft2(buffered_data)

        # Compute filter
        ny, nx = buffered_data.shape
        kx, ky = self.create_wavenumber_grids(nx, ny)
        filter_array = filter_function(kx, ky)

        # Apply filter
        filtered_fft = data_fft * filter_array

        # Inverse transform
        filtered_data = np.real(ifft2(filtered_fft))

        # Remove buffer and restore NaNs
        filtered_data = self.remove_buffer(filtered_data, buffer_size)
        return self.restore_nan(filtered_data, nan_mask)

    def bsdwormerttt(self, gridPath, num_levels, bottom_level, delta_z, shps, crs):
        from qgis.core import QgsTask, QgsApplication
        import numpy as np
        from math import ceil
        import os

        print(gridPath, num_levels, bottom_level, delta_z)
        # load grid and convert to numpy
        image, layer = loadGrid(gridPath)

        # Process initial grid setup
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
        )
        print("padded_size=", padded_image.shape)

        # Create and apply Hanning window
        window_y = np.hanning(padded_image.shape[0])
        window_x = np.hanning(padded_image.shape[1])
        hanning_window = np.outer(window_y, window_x)
        final_image = padded_image * hanning_window

        # Setup paths
        gridPathPadded = insert_text_before_extension(gridPath, "_padded")
        dir_name, base_name = os.path.split(gridPath)
        file_name, file_ext = os.path.splitext(base_name)
        wormsPath = dir_name + "/" + file_name + "_worms.csv"

        # Initialize primary Wormer instance for setup
        initial_job = Wormer()
        initial_job.importGdalRaster(gridPath)
        extent = initial_job.geomat
        extentPadded = GetExtent(initial_job.geomat, -pad_x, -pad_y)

        # Save padded image
        numpy_array_to_raster(
            final_image,
            gridPathPadded,
            pad_x,
            pad_y,
            dx=initial_job.geomat[1],
            xmin=extentPadded[2][0],
            ymax=extentPadded[2][1],
            reference_layer=layer,
            no_data_value=np.nan,
        )

        # Clear existing worms file
        if os.path.exists(wormsPath):
            os.remove(wormsPath)

        # Write header
        with open(wormsPath, "w") as f:
            f.write("y,x,z,val\n")

        class WormLevelTask(QgsTask):
            def __init__(
                self,
                dz,
                gridPathPadded,
                pad_x,
                pad_y,
                extent,
                extentPadded,
                delta_z,
                bottom_level,
                wormsPath,
            ):
                super().__init__(f"Process Worm Level {dz}", QgsTask.CanCancel)
                self.dz = dz
                self.gridPathPadded = gridPathPadded
                self.pad_x = pad_x
                self.pad_y = pad_y
                self.extent = extent
                self.extentPadded = extentPadded
                self.delta_z = delta_z
                self.bottom_level = bottom_level
                self.wormsPath = wormsPath

            def run(self):
                try:
                    # Create a new Wormer instance for this task
                    job = Wormer()
                    job.importGdalRaster(self.gridPathPadded)
                    job.importExternallyPaddedRaster(
                        self.gridPathPadded, self.pad_x, self.pad_y
                    )

                    dzm = (self.dz * self.delta_z) + self.bottom_level
                    if dzm == 0:
                        dzm = 0.01

                    job.wormLevelAsPoints(dz=dzm)
                    z = [dzm] * len(job.worm_ys)
                    x = (job.worm_xs * job.geomat[1]) + self.extentPadded[2][0]
                    y = (job.worm_ys * job.geomat[5]) + self.extentPadded[2][1]
                    vals = job.worm_vals

                    if len(job.worm_ys) == 0:
                        print(f"No points generated for level {self.dz}")
                        return True

                    points = np.column_stack(
                        (
                            y[: len(job.worm_ys)],
                            x[: len(job.worm_xs)],
                            z[: len(job.worm_ys)],
                            vals[: len(job.worm_ys)],
                        )
                    )

                    # Apply filter conditions
                    mask = (
                        (points[:, 1] >= self.extent[0])
                        & (
                            points[:, 1]
                            <= self.extent[0] + (job.raster_shape[1] * job.geomat[1])
                        )
                        & (
                            points[:, 0]
                            >= self.extent[3] + (job.raster_shape[0] * job.geomat[5])
                        )
                        & (points[:, 0] <= self.extent[3])
                    )

                    filtered_points = points[mask]

                    if len(filtered_points) > 0:
                        # Write to file with lock
                        with open(self.wormsPath, "a") as f:
                            np.savetxt(f, filtered_points, delimiter=",", fmt="%s")

                    return True
                except Exception as e:
                    print(f"Error processing level {self.dz}: {str(e)}")
                    return False

        # Create and queue tasks for each level
        task_manager = QgsApplication.taskManager()
        tasks = []

        for dz in range(0, num_levels):
            task = WormLevelTask(
                dz=dz,
                gridPathPadded=gridPathPadded,
                pad_x=pad_x,
                pad_y=pad_y,
                extent=extent,
                extentPadded=extentPadded,
                delta_z=delta_z,
                bottom_level=bottom_level,
                wormsPath=wormsPath,
            )
            task_manager.addTask(task)
            tasks.append(task)

        # Wait for all tasks to complete
        for task in tasks:
            task.waitForFinished()

        # Process shapefile if requested
        if shps:
            max_distance = initial_job.geomat[1] * 4
            in_path = wormsPath
            out_path = wormsPath.replace(".csv", ".shp")
            self.xyz_to_polylines(in_path, out_path, max_distance, crs)

    def bsdwormertt(self, gridPath, num_levels, bottom_level, delta_z, shps, crs):
        from qgis.core import QgsTask, QgsApplication
        import numpy as np
        from math import ceil
        import os

        print(gridPath, num_levels, bottom_level, delta_z)
        # load grid and convert to numpy
        image, layer = loadGrid(gridPath)

        # Process initial grid setup
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
        )
        print("padded_size=", padded_image.shape)

        # Create and apply Hanning window
        window_y = np.hanning(padded_image.shape[0])
        window_x = np.hanning(padded_image.shape[1])
        hanning_window = np.outer(window_y, window_x)
        final_image = padded_image * hanning_window

        # Setup paths
        gridPathPadded = insert_text_before_extension(gridPath, "_padded")
        dir_name, base_name = os.path.split(gridPath)
        file_name, file_ext = os.path.splitext(base_name)
        wormsPath = dir_name + "/" + file_name + "_worms.csv"

        # Initialize primary Wormer instance for setup
        initial_job = Wormer()
        initial_job.importGdalRaster(gridPath)
        extent = initial_job.geomat
        extentPadded = GetExtent(initial_job.geomat, -pad_x, -pad_y)

        # Save padded image
        numpy_array_to_raster(
            final_image,
            gridPathPadded,
            pad_x,
            pad_y,
            dx=initial_job.geomat[1],
            xmin=extentPadded[2][0],
            ymax=extentPadded[2][1],
            reference_layer=layer,
            no_data_value=np.nan,
        )

        # Clear existing worms file
        if os.path.exists(wormsPath):
            os.remove(wormsPath)

        # Write header
        with open(wormsPath, "w") as f:
            f.write("y,x,z,val\n")

        class WormLevelTask(QgsTask):
            def __init__(
                self,
                dz,
                gridPathPadded,
                pad_x,
                pad_y,
                extent,
                extentPadded,
                delta_z,
                bottom_level,
                wormsPath,
            ):
                super().__init__(f"Process Worm Level {dz}", QgsTask.CanCancel)
                self.dz = dz
                self.gridPathPadded = gridPathPadded
                self.pad_x = pad_x
                self.pad_y = pad_y
                self.extent = extent
                self.extentPadded = extentPadded
                self.delta_z = delta_z
                self.bottom_level = bottom_level
                self.wormsPath = wormsPath

            def run(self):
                try:
                    # Create a new Wormer instance for this task
                    job = Wormer()
                    job.importGdalRaster(self.gridPathPadded)
                    job.importExternallyPaddedRaster(
                        self.gridPathPadded, self.pad_x, self.pad_y
                    )

                    dzm = (self.dz * self.delta_z) + self.bottom_level
                    if dzm == 0:
                        dzm = 0.01

                    job.wormLevelAsPoints(dz=dzm)
                    z = [dzm] * len(job.worm_ys)
                    x = (job.worm_xs * job.geomat[1]) + self.extentPadded[2][0]
                    y = (job.worm_ys * job.geomat[5]) + self.extentPadded[2][1]
                    vals = job.worm_vals

                    if len(job.worm_ys) == 0:
                        print(f"No points generated for level {self.dz}")
                        return True

                    points = np.column_stack(
                        (
                            y[: len(job.worm_ys)],
                            x[: len(job.worm_xs)],
                            z[: len(job.worm_ys)],
                            vals[: len(job.worm_ys)],
                        )
                    )

                    # Apply filter conditions
                    mask = (
                        (points[:, 1] >= self.extent[0])
                        & (
                            points[:, 1]
                            <= self.extent[0] + (job.raster_shape[1] * job.geomat[1])
                        )
                        & (
                            points[:, 0]
                            >= self.extent[3] + (job.raster_shape[0] * job.geomat[5])
                        )
                        & (points[:, 0] <= self.extent[3])
                    )

                    filtered_points = points[mask]

                    if len(filtered_points) > 0:
                        # Write to file with lock
                        with open(self.wormsPath, "a") as f:
                            np.savetxt(f, filtered_points, delimiter=",", fmt="%s")

                    return True
                except Exception as e:
                    print(f"Error processing level {self.dz}: {str(e)}")
                    return False

        # Create and queue tasks for each level
        task_manager = QgsApplication.taskManager()
        tasks = []

        for dz in range(0, num_levels):
            task = WormLevelTask(
                dz=dz,
                gridPathPadded=gridPathPadded,
                pad_x=pad_x,
                pad_y=pad_y,
                extent=extent,
                extentPadded=extentPadded,
                delta_z=delta_z,
                bottom_level=bottom_level,
                wormsPath=wormsPath,
            )
            task_manager.addTask(task)
            tasks.append(task)

        # Wait for all tasks to complete
        for task in tasks:
            task.waitForFinished()

        # Process shapefile if requested
        if shps:
            max_distance = initial_job.geomat[1] * 1.3
            in_path = wormsPath
            out_path = wormsPath.replace(".csv", ".shp")
            self.xyz_to_polylines(in_path, out_path, max_distance, crs)

    def bsdwormerss(self, gridPath, num_levels, bottom_level, delta_z, shps, crs):
        from qgis.core import QgsTask, QgsApplication
        import numpy as np
        from math import ceil
        import os

        print(gridPath, num_levels, bottom_level, delta_z)
        # load grid and convert to numpy
        image, layer = loadGrid(gridPath)

        # Process initial grid setup (same as before)
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
        )
        print("padded_size=", padded_image.shape)

        # Create and apply Hanning window
        window_y = np.hanning(padded_image.shape[0])
        window_x = np.hanning(padded_image.shape[1])
        hanning_window = np.outer(window_y, window_x)
        final_image = padded_image * hanning_window

        # Setup job and paths
        job = Wormer()
        job.importGdalRaster(gridPath)
        extent = job.geomat
        extentPadded = GetExtent(job.geomat, -pad_x, -pad_y)
        gridPathPadded = insert_text_before_extension(gridPath, "_padded")

        dir_name, base_name = os.path.split(gridPath)
        file_name, file_ext = os.path.splitext(base_name)
        wormsPath = dir_name + "/" + file_name + "_worms.csv"

        # Save padded image
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

        # Clear existing worms file
        if os.path.exists(wormsPath):
            os.remove(wormsPath)

        # Write header
        with open(wormsPath, "w") as f:
            f.write("y,x,z,val\n")

        # Create a task for processing each dz level
        class WormLevelTask(QgsTask):
            def __init__(
                self, dz, job, extent, image, wormsPath, delta_z, bottom_level
            ):
                super().__init__(f"Process Worm Level {dz}", QgsTask.CanCancel)
                self.dz = dz
                self.job = job
                self.extent = extent
                self.image = image
                self.wormsPath = wormsPath
                self.delta_z = delta_z
                self.bottom_level = bottom_level

            def run(self):
                try:
                    dzm = (self.dz * self.delta_z) + self.bottom_level
                    if dzm == 0:
                        dzm = 0.01

                    self.job.wormLevelAsPoints(dz=dzm)
                    z = [dzm] * len(self.job.worm_ys)
                    x = (self.job.worm_xs * self.job.geomat[1]) + extentPadded[2][0]
                    y = (self.job.worm_ys * self.job.geomat[5]) + extentPadded[2][1]
                    vals = self.job.worm_vals

                    points = np.column_stack(
                        (
                            y[: len(self.job.worm_ys)],
                            x[: len(self.job.worm_xs)],
                            z[: len(self.job.worm_ys)],
                            vals[: len(self.job.worm_ys)],
                        )
                    )

                    # Apply filter conditions
                    mask = (
                        (points[:, 1] >= self.extent[0])
                        & (
                            points[:, 1]
                            <= self.extent[0]
                            + (self.image.shape[1] * self.job.geomat[1])
                        )
                        & (
                            points[:, 0]
                            >= self.extent[3]
                            + (self.image.shape[0] * self.job.geomat[5])
                        )
                        & (points[:, 0] <= self.extent[3])
                    )

                    filtered_points = points[mask]

                    # Write to file with lock to prevent concurrent writes
                    with open(self.wormsPath, "a") as f:
                        np.savetxt(f, filtered_points, delimiter=",", fmt="%s")

                    return True
                except Exception as e:
                    print(f"Error processing level {self.dz}: {str(e)}")
                    return False

        # Create and queue tasks for each level
        task_manager = QgsApplication.taskManager()
        tasks = []

        for dz in range(0, num_levels):
            task = WormLevelTask(
                dz, job, extent, image, wormsPath, delta_z, bottom_level
            )
            task_manager.addTask(task)
            tasks.append(task)

        # Wait for all tasks to complete
        for task in tasks:
            task.waitForFinished()

        # Process shapefile if requested
        if shps:
            max_distance = job.geomat[1] * 1.3
            in_path = wormsPath
            out_path = wormsPath.replace(".csv", ".shp")
            self.xyz_to_polylines(in_path, out_path, max_distance, crs)

    def bsdwormer(self, gridPath, num_levels, bottom_level, delta_z, shps, crs):
        # code borrows heavilly from bsdwormer example ipynb template example
        # adds on the fly calc of padded grid

        print(gridPath, num_levels, bottom_level, delta_z)
        # load grid and convert to numpy
        image, layer = loadGrid(gridPath)

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
        extentPadded = GetExtent(job.geomat, -pad_x, -pad_y)
        gridPathPadded = insert_text_before_extension(gridPath, "_padded")

        # Separate the file path into directory, base name, and extension
        dir_name, base_name = os.path.split(gridPath)
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

        for dz in range(0, num_levels):
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

            # Apply the filter conditions
            mask = (
                (points[:, 1] >= extent[0])  # x >= extent[0]
                & (
                    points[:, 1] <= extent[0] + (image.shape[1] * job.geomat[1])
                )  # x <= extent[0] + width
                & (
                    points[:, 0] >= extent[3] + (image.shape[0] * job.geomat[5])
                )  # y >= extent[3] + height
                & (points[:, 0] <= extent[3])  # y <= extent[3]
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
                np.savetxt(f, filtered_points, delimiter=",", fmt="%s")

        if shps:
            max_distance = job.geomat[1] * 1.3
            in_path = wormsPath
            out_path = wormsPath.replace(".csv", ".shp")
            self.xyz_to_polylines(in_path, out_path, max_distance, crs)

    def display_grid(self, grid):
        import matplotlib.pyplot as plt

        # Plot results
        plt.imshow(grid, cmap="gray")
        plt.colorbar()
        plt.title("Grid")
        plt.show()

    def split_long_segments(self, line, max_distance):
        """Split a LineString into segments, interpolating points where needed."""
        from shapely.geometry import Point, LineString
        import numpy as np

        coords = list(line.coords)
        new_lines = []

        for i in range(len(coords) - 1):
            p1 = coords[i]
            p2 = coords[i + 1]

            # Calculate distance between consecutive points
            segment_length = Point(p1).distance(Point(p2))

            if segment_length <= max_distance:
                # If segment is within threshold, add it directly
                new_lines.append(LineString([p1, p2]))
            else:
                # If segment exceeds threshold, interpolate points
                num_segments = int(np.ceil(segment_length / max_distance))

                for j in range(num_segments):
                    # Calculate interpolated points
                    t1 = j / num_segments
                    t2 = (j + 1) / num_segments

                    interpolated_p1 = (
                        p1[0] + t1 * (p2[0] - p1[0]),
                        p1[1] + t1 * (p2[1] - p1[1]),
                    )
                    interpolated_p2 = (
                        p1[0] + t2 * (p2[0] - p1[0]),
                        p1[1] + t2 * (p2[1] - p1[1]),
                    )

                    new_lines.append(LineString([interpolated_p1, interpolated_p2]))

        return new_lines

    def process_cluster(self, points, z, max_distance):
        """Process a cluster of points to create split LineStrings with better connectivity."""
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
                # Find the closest unvisited point to any point in our current path
                min_dist = float("inf")
                best_pair = None

                for i, point in enumerate(sorted_points):
                    distances, indices = tree.query(point, k=len(points))
                    for j, (idx, dist) in enumerate(zip(indices, distances)):
                        if not visited[idx] and dist < min_dist:
                            min_dist = dist
                            best_pair = (i, idx)

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

    def xyz_to_polylines(self, csv_file, output_shapefile, max_distance, crs):
        """
        Converts XYZ points from a CSV file into polylines and saves them as a shapefile.
        Parameters:
        csv_file (str): Path to the input CSV file containing XYZ points.
        output_shapefile (str): Path to the output shapefile where polylines will be saved.
        max_distance (float): Maximum distance between points to be considered part of the same cluster.
        crs (int): Coordinate Reference System (CRS) EPSG code for the output shapefile.
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

        for z in unique_z:

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

            for label in unique_labels:
                cluster_points = points[labels == label]
                if len(cluster_points) > 1:
                    try:
                        result = self.process_cluster(cluster_points, z, max_distance)
                        if result:
                            results.extend(result)
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
        layer = ds.CreateLayer("polylines", spatial_ref, ogr.wkbLineString)

        field_defn = ogr.FieldDefn("Z", ogr.OFTInteger)
        layer.CreateField(field_defn)

        for polyline, z in zip(polylines, z_values):
            feature = ogr.Feature(layer.GetLayerDefn())
            feature.SetField("Z", int(z))
            geom = ogr.CreateGeometryFromWkt(polyline.wkt)
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

    def compute_reference_stats(self, reference_tiff):
        """Compute mean and standard deviation from the first (reference) GeoTIFF."""
        ds = gdal.Open(reference_tiff)
        band = ds.GetRasterBand(1)
        data = band.ReadAsArray().astype(np.float32)
        data, nodata_value = self.fix_extreme_values(data)
        mask = (data == nodata_value) | np.isnan(data)

        detrended_data = self.remove_gradient(data, mask)
        valid_data = detrended_data[~mask]
        mean, std = (
            (np.mean(valid_data), np.std(valid_data)) if valid_data.size > 0 else (0, 1)
        )

        return mean, std

    def normalise_geotiffs(self, input_folder, output_folder, order):
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
        for tiff_file in tiff_files:
            # print("tiff_file", tiff_file)
            ds = gdal.Open(tiff_file, gdal.GA_ReadOnly)
            band = ds.GetRasterBand(1)
            data = band.ReadAsArray()
            geotransform = ds.GetGeoTransform()
            projection = ds.GetProjection()

            data, nodata_value = self.fix_extreme_values(data)
            data = data.astype(np.float32)
            mask = (data == nodata_value) | np.isnan(data)

            # print("nodata_value:", nodata_value)
            # print("Unique values in data:", np.unique(data))

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

            # print("nonnan", np.count_nonzero(~np.isnan(valid_data)))
            # print("max min", np.max(valid_data), np.min(valid_data))

            if first:
                first = False
                reference_std = valid_std  # Use the non-zero std

            normalized_data = (
                reference_std * (detrended_data - np.nanmean(valid_data)) / valid_std
            )

            """print(
                "reference_std, np.nanmean(valid_data), np.nanstd(valid_data)",
                reference_std,
                np.nanmean(valid_data),
                valid_std,
            )"""
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
    #import matplotlib.pyplot as plt

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
