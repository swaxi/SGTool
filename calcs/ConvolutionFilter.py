import numpy as np
from scipy.ndimage import convolve, median_filter, gaussian_filter
from PyQt5.QtGui import QValidator


class OddPositiveIntegerValidator(QValidator):
    def validate(self, input_text, pos):
        """
        Validates the text input to allow only odd positive integers.
        """
        if not input_text:  # Allow empty input (to clear the field)
            return QValidator.Intermediate, input_text, pos

        if not input_text.isdigit():  # Only digits are allowed
            return QValidator.Invalid, input_text, pos

        value = int(input_text)
        if value > 0 and value % 2 == 1:  # Check for positive odd numbers
            return QValidator.Acceptable, input_text, pos
        else:
            return QValidator.Intermediate, input_text, pos

    def fixup(self, input_text):
        """
        Corrects invalid input to the nearest odd positive integer.
        """
        try:
            value = int(input_text)
            if value <= 0:  # Make it a positive number
                return "1"
            elif value % 2 == 0:  # Make it odd
                return str(value + 1)
            else:
                return input_text
        except ValueError:
            return "1"  # Default to 1 if the input cannot be converted


class ConvolutionFilter:
    def __init__(self, grid):
        """
        Initialize the ConvolutionFilter with a grid.

        :param grid: 2D numpy array representing the input grid
        """
        self.grid = np.array(grid, dtype=float)
        self.padded_grid = None

    def apply_padding(self, pad_width):
        """
        Apply reflective padding to the grid.

        :param pad_width: Width of the padding
        :return: Padded grid
        """
        self.padded_grid = np.pad(self.grid, pad_width, mode="reflect")
        return self.padded_grid

    def nan_convolution(self, kernel, mode="reflect"):
        """
        Perform convolution while handling NaN values.

        :param kernel: Convolution kernel
        :param mode: Padding mode (default is 'reflect')
        :return: Convolved grid with NaN handling
        """
        # Create a mask for non-NaN values
        valid_mask = ~np.isnan(self.grid)

        # Replace NaNs with 0 for convolution
        grid_filled = np.nan_to_num(self.grid, nan=0.0)

        # Convolve the filled grid and the valid mask
        convolved_values = convolve(grid_filled, kernel, mode=mode)
        valid_counts = convolve(valid_mask.astype(float), kernel, mode=mode)

        # Avoid division by zero
        valid_counts[valid_counts == 0] = np.nan

        # Calculate the mean of valid values
        return convolved_values / valid_counts

    def mean_filter(self, n):
        """
        Apply mean filter while handling NaN values.

        :param n: Size of the kernel (n x n)
        :return: Filtered grid
        """
        kernel = np.ones((n, n), dtype=float)
        return self.nan_convolution(kernel)

    def median_filter(self, n):
        """
        Apply median filter while handling NaN values.

        :param n: Size of the kernel (n x n)
        :return: Filtered grid
        """

        # Use a sliding window approach with NaN handling
        def nanmedian(values):
            return np.nan if np.isnan(values).all() else np.nanmedian(values)

        output = np.zeros_like(self.grid, dtype=float)
        pad_width = n // 2
        padded_grid = np.pad(
            self.grid, pad_width, mode="constant", constant_values=np.nan
        )

        for i in range(output.shape[0]):
            for j in range(output.shape[1]):
                window = padded_grid[i : i + n, j : j + n]
                output[i, j] = nanmedian(window)

        return output

    def gaussian_filter(self, sigma):
        """
        Apply Gaussian filter while handling NaN values.

        :param sigma: Standard deviation for Gaussian kernel
        :return: Filtered grid
        """
        # Create a Gaussian kernel
        size = int(2 * np.ceil(2 * sigma) + 1)
        x = np.linspace(-size // 2, size // 2, size)
        gaussian_kernel = np.exp(-(x**2 / (2 * sigma**2)))
        gaussian_kernel = np.outer(gaussian_kernel, gaussian_kernel)
        gaussian_kernel /= gaussian_kernel.sum()

        return self.nan_convolution(gaussian_kernel)

    def directional_filter(self, direction, n=3):
        """
        Apply directional filter (NE, N, NW, W, SW, S, SE, E).

        :param direction: Direction of the filter ('NE', 'N', 'NW', 'W', 'SW', 'S', 'SE', 'E')
        :param n: Size of the kernel (n x n, default is 3x3)
        :return: Filtered grid
        """
        direction_kernels = {
            "N": np.array([[-1, -1, -1], [0, 0, 0], [1, 1, 1]]),
            "S": np.array([[1, 1, 1], [0, 0, 0], [-1, -1, -1]]),
            "E": np.array([[-1, 0, 1], [-1, 0, 1], [-1, 0, 1]]),
            "W": np.array([[1, 0, -1], [1, 0, -1], [1, 0, -1]]),
            "NE": np.array([[0, -1, -1], [1, 0, -1], [1, 1, 0]]),
            "NW": np.array([[-1, -1, 0], [-1, 0, 1], [0, 1, 1]]),
            "SE": np.array([[0, 1, 1], [-1, 0, 1], [-1, -1, 0]]),
            "SW": np.array([[1, 1, 0], [1, 0, -1], [0, -1, -1]]),
        }

        if direction not in direction_kernels:
            raise ValueError(
                f"Invalid direction '{direction}'. Must be one of {list(direction_kernels.keys())}."
            )

        kernel = direction_kernels[direction]
        kernel_size = kernel.shape[0]

        if kernel_size != n:
            kernel = np.pad(
                kernel,
                ((n - kernel_size) // 2, (n - kernel_size) // 2),
                mode="constant",
            )

        return convolve(self.grid, kernel, mode="reflect")

    def sun_shading_filter(self, elevation, sun_alt=45.0, sun_az=315.0, resolution=1.0):
        """
        Compute relief shading for a digital elevation model.

        Parameters:
            elevation (numpy.ndarray): 2D array of elevation data (DEM).
            sun_alt (float): Sun altitude in degrees (default is 45.0).
            sun_az (float): Sun azimuth in degrees clockwise from north (default is 315.0).
            resolution (float): Resolution of the grid (default is 1.0).

        Returns:
            numpy.ndarray: 2D array of relief shading values.
        """

        # sometimes doesn't like 90 so offset a bit
        if sun_alt == 90.0:
            sun_alt = 88.0

        # Convert sun altitude and azimuth to radians
        sun_alt_rad = np.radians(sun_alt)
        sun_az_rad = np.radians(sun_az)

        # Compute light source vector
        sun_vec = np.array(
            [
                np.cos(sun_alt_rad) * np.sin(sun_az_rad),  # x component
                np.cos(sun_alt_rad) * np.cos(sun_az_rad),  # y component
                np.sin(sun_alt_rad),  # z component
            ]
        )

        # Calculate gradients using finite differences
        dzdx = (np.roll(elevation, -1, axis=1) - np.roll(elevation, 1, axis=1)) / (
            2 * resolution
        )
        dzdy = (np.roll(elevation, -1, axis=0) - np.roll(elevation, 1, axis=0)) / (
            2 * resolution
        )

        # Compute normal vectors
        norm_x = -dzdx
        norm_y = -dzdy
        norm_z = 1.0

        # Normalize the normal vectors
        norm_length = np.sqrt(norm_x**2 + norm_y**2 + norm_z**2)
        norm_x /= norm_length
        norm_y /= norm_length
        norm_z /= norm_length

        # Dot product with sun vector
        shading = norm_x * sun_vec[0] + norm_y * sun_vec[1] + norm_z * sun_vec[2]

        # Clamp shading values to range [0, 1]
        # shading = np.clip(shading, 0, 1)

        return shading

    def sun_shading_filter_grass(
        self,
        elevation,
        resolution_ns,
        resolution_ew,
        altitude=30.0,
        azimuth=270.0,
        scale=1.0,
        zscale=1.0,
    ):
        """
        Calculate shaded relief from an elevation map following GRASS r.relief implementation.

        Parameters:
        -----------
        elevation : np.ndarray
            Input elevation data (2D array)
        resolution_ns : float
            North-south resolution (vertical) in same units as elevation
        resolution_ew : float
            East-west resolution (horizontal) in same units as elevation
        altitude : float
            Sun altitude in degrees above horizon (default: 30)
        azimuth : float
            Sun azimuth in degrees east of north (default: 270)
        scale : float
            Scale factor for converting meters to elevation units (default: 1.0)
        zscale : float
            Factor for exaggerating relief (default: 1.0)

        Returns:
        --------
        np.ndarray
            Shaded relief array with values 0-255
        """
        # Convert angles to radians
        degrees_to_radians = np.pi / 180.0
        altitude = altitude * degrees_to_radians
        # Correct azimuth to East (GRASS convention)
        azimuth = (azimuth - 90.0) * degrees_to_radians

        # Calculate distances (following GRASS implementation)
        H = resolution_ew * 4 * scale / zscale  # horizontal run for gradient
        V = resolution_ns * 4 * scale / zscale  # vertical run for gradient

        # Pad the elevation array to handle edges
        elev_padded = np.pad(elevation, pad_width=1, mode="edge")

        nrows, ncols = elevation.shape
        output = np.zeros((nrows, ncols), dtype=np.float64)

        # Main computation loop following GRASS implementation
        for row in range(1, nrows + 1):
            for col in range(1, ncols + 1):
                window = elev_padded[row - 1 : row + 2, col - 1 : col + 2]

                # Skip if any null/nan values in 3x3 window
                if np.any(np.isnan(window)):
                    output[row - 1, col - 1] = np.nan
                    continue

                # Calculate gradients using GRASS method
                dx = (
                    window[0, 0]
                    + 2 * window[1, 0]
                    + window[2, 0]
                    - window[0, 2]
                    - 2 * window[1, 2]
                    - window[2, 2]
                ) / H

                dy = (
                    window[0, 0]
                    + 2 * window[0, 1]
                    + window[0, 2]
                    - window[2, 0]
                    - 2 * window[2, 1]
                    - window[2, 2]
                ) / V

                # Calculate slope and aspect
                key = dx * dx + dy * dy
                slope = np.pi / 2.0 - np.arctan(np.sqrt(key))

                # Aspect (GRASS implementation)
                aspect = np.arctan2(dy, dx)
                if np.isnan(aspect):
                    aspect = degrees_to_radians
                if dx != 0 or dy != 0:
                    if aspect == 0:
                        aspect = 2 * np.pi

                # Calculate shaded relief
                cang = np.sin(altitude) * np.sin(slope) + np.cos(altitude) * np.cos(
                    slope
                ) * np.cos(azimuth - aspect)

                # Scale to 0-255 range
                output[row - 1, col - 1] = 255.0 * cang

        # Clip values to 0-255 range
        # output = np.clip(output, 0, 255)
        return output
