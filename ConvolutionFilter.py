import numpy as np
from scipy.ndimage import convolve, median_filter, gaussian_filter
from PyQt5.QtGui import QValidator
from qgis.core import QgsProject


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
    def __init__(self, dx, dy, localGridName):
        """
        Initialize the ConvolutionFilter with a grid.

        :param grid: 2D numpy array representing the input grid
        """
        self.padded_grid = None
        self.dx = dx
        self.dy = dy
        self.localGridName = localGridName

    def apply_padding(self, grid, pad_width):
        """
        Apply reflective padding to the grid.

        :param pad_width: Width of the padding
        :return: Padded grid
        """
        self.padded_grid = np.pad(grid, pad_width, mode="reflect")
        return self.padded_grid

    def mean_filter(self, grid, n):
        """
        Apply mean filter.

        :param n: Size of the kernel (n x n)
        :return: Filtered grid
        """
        kernel = np.ones((n, n)) / (n * n)
        new_grid = convolve(grid, kernel, mode="reflect")

        return new_grid

    def median_filter(self, grid, n):
        """
        Apply median filter.

        :param n: Size of the kernel (n x n)
        :return: Filtered grid
        """
        return median_filter(grid, size=(n, n), mode="reflect")

    def gaussian_filter(self, grid, sigma):
        """
        Apply Gaussian filter.

        :param sigma: Standard deviation for Gaussian kernel
        :return: Filtered grid
        """
        return gaussian_filter(grid, sigma=sigma, mode="reflect")

    def directional_filter(self, grid, direction, n=3):
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

        return convolve(grid, kernel, mode="reflect")

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
