import numpy as np
from scipy.fftpack import fft2, ifft2, fftfreq
from scipy.ndimage import gaussian_filter
from numpy.polynomial.polynomial import polyval2d
from scipy.spatial import cKDTree
from math import ceil
from scipy.interpolate import interpn


from osgeo import ogr, osr
from shapely.geometry import LineString, Point

# from joblib import Parallel, delayed
import csv

import glob
import os
from osgeo import gdal
from scipy.optimize import leastsq
import platform

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

    def fill_nan(self, data):
        """
        Replace NaN values with the mean of the non-NaN values.
        """
        nan_mask = np.isnan(data)
        # filled_data = np.copy(data)
        # filled_data[nan_mask] = np.nanmedian(data)

        # Create a grid of coordinates
        x = np.arange(data.shape[0])
        y = np.arange(data.shape[1])
        grid_x, grid_y = np.meshgrid(x, y, indexing="ij")

        # Replace NaN values temporarily with zeros to enable interpolation
        # (valid points will define the interpolation behavior)
        data_temp = np.copy(data)
        data_temp[np.isnan(data)] = (
            0  # Use a placeholder (will be replaced by interpolated values)
        )

        # Interpolation grid for all points
        all_coords = np.array((grid_x.ravel(), grid_y.ravel())).T

        # Perform interpolation
        filled_values = interpn(
            (x, y),
            data_temp,
            all_coords,
            method="linear",
            bounds_error=False,
            fill_value=None,
        )

        # Reshape back to the original array shape
        filled_array = filled_values.reshape(data.shape)

        # self.display_grid(filled_array)

        return filled_array, nan_mask

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

        # Apply the Hanning window to the padded image
        final_image = padded_image * hanning_window
        return final_image

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
        self,
        data,
        inclination,
        declination,
        magnetization_inclination=None,
        magnetization_declination=None,
        buffer_size=10,
        buffer_method="mirror",
    ):
        """
        Apply reduction to the pole filter in the frequency domain.
        (Based on Fatiando a Terra https://github.com/fatiando/harmonica/blob/main/harmonica/filters/_filters.py)

        Parameters
        ----------
        data : array-like
            Input magnetic data to be transformed
        inclination : float
            The inclination of the inducing Geomagnetic field in degrees
        declination : float
            The declination of the inducing Geomagnetic field in degrees
        magnetization_inclination : float or None
            The inclination of the total magnetization of the anomaly source in degrees
        magnetization_declination : float or None
            The declination of the total magnetization of the anomaly source in degrees
        buffer_size : int
            Size of the buffer zone for reducing edge effects
        buffer_method : str
            Method for handling the buffer zone ("mirror", "mean", etc.)

        Returns
        -------
        array-like
            The reduced to pole magnetic data
        """
        if inclination < 0:
            inclination = -inclination
            declination = -declination

        if magnetization_declination is None and magnetization_inclination is None:
            magnetization_inclination = inclination
            magnetization_declination = declination

        # Transform to radians
        inc_rad = np.radians(inclination)
        dec_rad = -np.radians(declination)

        # Calculate the 3 components
        m_e = np.cos(inc_rad) * np.sin(dec_rad)
        m_n = np.cos(inc_rad) * np.cos(dec_rad)
        m_z = np.sin(inc_rad)

        f_e, f_n, f_z = m_e, m_n, m_z

        def filter_function(kx, ky):
            """
            Generate the reduction to pole filter in frequency domain with improved stability.
            """
            k_squared = kx**2 + ky**2
            k = np.sqrt(k_squared)

            # Calculate denominators separately
            denom1 = f_z * k + 1j * (f_e * kx + f_n * ky)
            denom2 = m_z * k + 1j * (m_e * kx + m_n * ky)

            # Small value to prevent division by zero
            epsilon = 1e-10

            # Create masks for small values
            mask1 = np.abs(denom1) > epsilon
            mask2 = np.abs(denom2) > epsilon
            mask = mask1 & mask2

            # Initialize filter array with zeros
            rtp_filter = np.zeros_like(k, dtype=complex)

            # Calculate filter only where denominators are large enough
            rtp_filter[mask] = k_squared[mask] / (denom1[mask] * denom2[mask])

            # Calculate theta (angle between wavenumber vector and north)
            theta = np.abs(np.arctan2(ky, kx) - dec_rad)

            # Apply tapers and stabilization
            alpha = 0.5
            theta_cap = np.radians(65)  # Convert to radians

            # Wavenumber taper
            k_max = np.max(k)
            k_taper = np.exp(-alpha * (k / k_max) ** 2 / 2)

            # Directional taper
            dir_taper = np.ones_like(theta)
            critical_angles = theta < theta_cap
            dir_taper[critical_angles] = 0.2 + 0.8 * np.sin(theta[critical_angles])

            # Combined taper
            taper = k_taper * dir_taper

            # Apply taper
            rtp_filter *= taper

            # Additional stability constraints for low inclinations
            incl_factor = np.abs(np.sin(inc_rad))
            if incl_factor < 0.2:  # Near equator
                rtp_filter *= np.exp(-((1 - incl_factor) * (k / k_max) ** 2))

            # Set specific values to zero
            rtp_filter[k == 0] = 0  # DC component
            rtp_filter[~np.isfinite(rtp_filter)] = 0  # Non-finite values
            rtp_filter[np.abs(rtp_filter) > 1e6] = 0  # Extremely large values

            return rtp_filter

        return self._apply_fourier_filter(
            data, filter_function, buffer_size, buffer_method
        )

    def reduction_to_equator(
        self,
        data,
        inclination,
        declination,
        buffer_size=20,  # Increased buffer size
        buffer_method="mirror",
        epsilon=1e-06,
    ):
        """
        Apply reduction to the equator filter with enhanced anti-striping measures.
        """

        # If inclination/declination of the source is unknown, apply a small perturbation
        inclination2 = (
            inclination  # + 0.1  # Small perturbation to avoid zero numerator
        )
        declination2 = declination  # + 0.1

        # Convert degrees to radians
        inc_rad = np.radians(-inclination)
        dec_rad = np.radians(-declination)
        inclination_rad = np.radians(-inclination2)
        declination_rad = np.radians(-declination2)

        # Get data shape
        ny, nx = data.shape

        # Create properly scaled wavenumber grids
        kx = np.fft.fftfreq(nx, d=self.dx) * (
            2 * np.pi
        )  # Convert to spatial frequency (radians per unit)
        ky = np.fft.fftfreq(ny, d=self.dy) * (2 * np.pi)
        kx, ky = np.meshgrid(kx, ky)

        def filter_function(kx, ky):
            # Compute wavenumber magnitude (ensuring nonzero values)
            k = np.sqrt(kx**2 + ky**2)
            k = np.where(k < epsilon, epsilon, k)  # Avoid zero-division

            # Compute Fourier-space angle theta
            theta = np.arctan2(ky, kx)

            # Compute magnetic field components with correction
            P1 = np.cos(inclination_rad) * np.sin(inc_rad) - np.sin(
                inclination_rad
            ) * np.cos(inc_rad) * np.cos(dec_rad - declination_rad)
            P2 = (
                np.cos(inclination_rad)
                * np.cos(inc_rad)
                * np.sin(dec_rad - declination_rad)
            )

            # Ensure P1 and P2 are not zero by adding a minimal perturbation if needed
            P1 = np.where(np.abs(P1) < epsilon, epsilon, P1)
            P2 = np.where(np.abs(P2) < epsilon, epsilon, P2)

            # Compute denominator and ensure non-zero values
            denom = np.cos(inc_rad) + 1j * np.sin(inc_rad) * np.sin(theta - dec_rad)
            denom = np.where(np.abs(denom) < epsilon, epsilon, denom)

            # Compute the improved RTE filter
            F_rte = (P1 + 1j * P2) / denom
            # Set specific values to zero
            F_rte[k == 0] = 0  # DC component
            F_rte[~np.isfinite(F_rte)] = 0  # Non-finite values
            F_rte[np.abs(F_rte) > 1e6] = 0  # Extremely large values

            return F_rte

        return self._apply_fourier_filter(
            data, filter_function, buffer_size, buffer_method
        )

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

    def tilt_derivative(self, data, buffer_size=10, buffer_method="mirror"):
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
        # inpaint NaN values
        # Handle NaN values
        filled_data, nan_mask = self.fill_nan(data)
        # Add buffer
        buffered_data = self.add_buffer(filled_data, buffer_size, buffer_method)

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
            max_distance = initial_job.geomat[1] * 1.3
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

    def process_cluster(self, points, z, max_distance):
        """Process a cluster of points to create split LineStrings."""
        if len(points) > 1:
            tree = cKDTree(points)
            sorted_points = [points[0]]
            remaining_points = set(range(1, len(points)))

            while remaining_points:
                last_point = sorted_points[-1]
                distances, indices = tree.query(
                    last_point, k=min(len(remaining_points), len(points))
                )
                if np.isscalar(indices):
                    indices = [indices]
                valid_indices = [idx for idx in indices if idx in remaining_points]

                if valid_indices:
                    nearest_idx = valid_indices[0]
                    sorted_points.append(points[nearest_idx])
                    remaining_points.remove(nearest_idx)
                else:
                    break

            if len(sorted_points) > 1:
                polyline = LineString(sorted_points)
                split_segments = self.split_long_segments(polyline, max_distance)
                return [(seg, z) for seg in split_segments]
        return []

    def split_long_segments(self, line, max_distance):
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

    def process_cluster(self, points, z, max_distance):
        """Process a cluster of points to create split LineStrings."""
        if len(points) > 1:
            tree = cKDTree(points)
            sorted_points = [points[0]]
            remaining_points = set(range(1, len(points)))

            while remaining_points:
                last_point = sorted_points[-1]
                distances, indices = tree.query(
                    last_point, k=min(len(remaining_points), len(points))
                )
                if np.isscalar(indices):
                    indices = [indices]
                valid_indices = [idx for idx in indices if idx in remaining_points]

                if valid_indices:
                    nearest_idx = valid_indices[0]
                    sorted_points.append(points[nearest_idx])
                    remaining_points.remove(nearest_idx)
                else:
                    break

            if len(sorted_points) > 1:
                polyline = LineString(sorted_points)
                split_segments = self.split_long_segments(polyline, max_distance)
                return [(seg, z) for seg in split_segments]
        return []

    def xyz_to_polylinestt(self, csv_file, output_shapefile, max_distance, crs):
        """
        Converts XYZ points from a CSV file into polylines and saves them as a shapefile.
        Processes different z-levels in parallel using PyQGIS task management.
        """
        from sklearn.cluster import DBSCAN
        from qgis.core import QgsTask, QgsApplication
        import csv
        from queue import Queue
        import threading

        # Class to handle processing of individual z-levels
        class ZLevelProcessor(QgsTask):
            def __init__(self, z_value, points, max_distance, processor):
                super().__init__(f"Process Z Level {z_value}", QgsTask.CanCancel)
                self.z_value = z_value
                self.points = points
                self.max_distance = max_distance
                self.processor = processor
                self.results = []

            def run(self):
                try:
                    if len(self.points) < 2 or self.z_value < 1.0:
                        return True

                    print(
                        f"Processing Z-bin: {self.z_value}, Number of points: {len(self.points)}"
                    )

                    clustering = DBSCAN(eps=self.max_distance, min_samples=1).fit(
                        self.points
                    )
                    labels = clustering.labels_
                    unique_labels = np.unique(labels)
                    unique_labels = unique_labels[unique_labels != -1]

                    for label in unique_labels:
                        cluster_points = self.points[labels == label]
                        if len(cluster_points) > 1:
                            try:
                                result = self.processor.process_cluster(
                                    cluster_points, self.z_value, self.max_distance
                                )
                                if result:
                                    self.results.extend(result)
                            except Exception as e:
                                print(
                                    f"Error processing cluster at Z-bin {self.z_value}: {e}"
                                )
                                continue

                    return True
                except Exception as e:
                    print(f"Error processing Z level {self.z_value}: {e}")
                    return False

        # Read and validate CSV data
        valid_data = []
        skipped_rows = 0

        with open(csv_file, "r") as f:
            reader = csv.DictReader(f)
            for row_num, row in enumerate(reader, start=2):
                try:
                    if not all(key in row and row[key] for key in ["x", "y", "z"]):
                        skipped_rows += 1
                        continue

                    x = float(row["x"].strip())
                    y = float(row["y"].strip())
                    z = float(row["z"].strip())

                    if not all(map(np.isfinite, [x, y, z])):
                        skipped_rows += 1
                        continue

                    valid_data.append((x, y, z))

                except (ValueError, AttributeError) as e:
                    skipped_rows += 1
                    print(f"Warning: Skipping invalid row {row_num}: {str(e)}")
                    continue

        if skipped_rows > 0:
            print(f"Skipped {skipped_rows} rows due to missing or invalid data")

        if not valid_data:
            raise ValueError("No valid data points found in the CSV file")

        # Convert to numpy array and get unique z values
        data = np.array(valid_data)
        unique_z = np.unique(data[:, 2])

        # Create tasks for each z level
        task_manager = QgsApplication.taskManager()
        tasks = []
        all_results = []

        # Process each z level in parallel
        for z in unique_z:
            mask = data[:, 2] == z
            points = data[mask][:, :2]

            task = ZLevelProcessor(z, points, max_distance, self)
            tasks.append(task)
            task_manager.addTask(task)

        # Wait for all tasks to complete and collect results
        for task in tasks:
            task.waitForFinished()
            if task.results:
                all_results.extend(task.results)

        # Collect all polylines and corresponding Z values
        polylines = [item[0] for item in all_results if isinstance(item[0], LineString)]
        z_values = [item[1] for item in all_results if isinstance(item[0], LineString)]

        if not polylines:
            raise ValueError(
                "No valid LineStrings were created. Check your input data or parameters."
            )

        print("Processing complete. Saving to shapefile...")

        # Save results to shapefile
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

        ds = None
        print(f"Shapefile saved to {output_shapefile}")

    def xyz_to_polylineszz(self, csv_file, output_shapefile, max_distance, crs):
        """
        Converts XYZ points from a CSV file into polylines and saves them as a shapefile.
        Handles missing or invalid data in the CSV file.
        """
        from sklearn.cluster import DBSCAN
        import csv

        valid_data = []
        skipped_rows = 0

        # Read CSV with robust error handling
        with open(csv_file, "r") as f:
            reader = csv.DictReader(f)
            for row_num, row in enumerate(
                reader, start=2
            ):  # start=2 because row 1 is header
                try:
                    # Check if all required fields exist and are not None/empty
                    if not all(key in row and row[key] for key in ["x", "y", "z"]):
                        skipped_rows += 1
                        continue

                    # Try to convert values to float
                    x = float(row["x"].strip())
                    y = float(row["y"].strip())
                    z = float(row["z"].strip())

                    # Skip rows with NaN values
                    if not all(map(np.isfinite, [x, y, z])):
                        skipped_rows += 1
                        continue

                    valid_data.append((x, y, z))

                except (ValueError, AttributeError) as e:
                    skipped_rows += 1
                    print(f"Warning: Skipping invalid row {row_num}: {str(e)}")
                    continue

        if skipped_rows > 0:
            print(f"Skipped {skipped_rows} rows due to missing or invalid data")

        if not valid_data:
            raise ValueError("No valid data points found in the CSV file")

        data = np.array(valid_data)
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
    import matplotlib.pyplot as plt

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
    plt.title("tilt_derivative")
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
    tilt_result = processor.tilt_derivative(data)

    # Plot results
    plt.imshow(
        tilt_result, cmap="seismic", origin="lower", vmin=-np.pi / 2, vmax=np.pi / 2
    )
    plt.colorbar()
    plt.title("Tilt Derivative")
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
