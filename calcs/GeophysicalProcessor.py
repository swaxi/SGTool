import numpy as np
from scipy.fftpack import fft2, ifft2, fftfreq
from scipy.ndimage import gaussian_filter
from numpy.polynomial.polynomial import polyval2d
from scipy.spatial import cKDTree
from math import ceil
from scipy.interpolate import interpn


from osgeo import ogr, osr
from shapely.geometry import LineString, Point
from joblib import Parallel, delayed
import csv

import glob
import os
from osgeo import gdal
from scipy.optimize import leastsq
import platform

if platform.system() == "Windows":
    from sgtool.worms.wormer import Wormer
    from sgtool.worms.Utility import (
        GetExtent,
        loadGrid,
        numpy_array_to_raster,
        insert_text_before_extension,
        fill_nan,
    )
else:
    from sgtool.worms.wormer import Wormer
    from sgtool.worms.Utility import (
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
        self, data, inclination, declination, buffer_size=10, buffer_method="mirror"
    ):
        """
        Perform Reduction to the Pole (RTP) transformation on geophysical data.
        This method transforms magnetic data to what it would look like if the
        magnetic field were vertical. This is useful for interpreting magnetic
        anomalies.
        Parameters:
        -----------
        data : np.ndarray
            2D array of magnetic data to be transformed.
        inclination : float
            Inclination angle of the magnetic field in degrees.
        declination : float
            Declination angle of the magnetic field in degrees.
        buffer_size : int, optional
            Size of the buffer to be added around the data to minimize edge effects,
            by default 10.
        buffer_method : str, optional
            Method to use for buffering the data, by default "mirror". Other options
            might include "constant", "nearest", etc.
        Returns:
        --------
        np.ndarray
            2D array of the RTP transformed data.
        """
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
        with np.errstate(divide="ignore", invalid="ignore"):
            rtp_filter = (
                k * cos_inc * cos_dec + 1j * ky * cos_inc * sin_dec + kx * sin_inc
            ) / k

        rtp_filter[np.abs(rtp_filter) > 1e6] = 0  # Stabilize extreme values

        # Apply max_wavenumber filtering based on cell size
        # rtp_filter[k > 1/(float(self.dx)*2)] = 0

        # Apply RTP filter in the Fourier domain
        data_fft = fft2(buffered_data)
        data_rtp_fft = data_fft * rtp_filter.T

        data_rtp = np.real(ifft2(data_rtp_fft))

        # Remove buffer and restore NaNs
        filtered_data = self.remove_buffer(data_rtp, buffer_size)
        return self.restore_nan(filtered_data, nan_mask)

    def reduction_to_equator(
        self, data, inclination, declination, buffer_size=10, buffer_method="mirror"
    ):
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

            if first:
                first = False
                reference_std = np.std(valid_data)

            normalized_data = (
                reference_std
                * (detrended_data - np.mean(valid_data))
                / np.std(valid_data)
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
