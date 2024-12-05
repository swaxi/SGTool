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
    def __init__(self,dx, dy,localGridName):
        """
        Initialize the ConvolutionFilter with a grid.

        :param grid: 2D numpy array representing the input grid
        """
        self.padded_grid = None
        self.dx = dx
        self.dy = dy
        self.localGridName=localGridName

    def apply_padding(self, grid, pad_width):
        """
        Apply reflective padding to the grid.

        :param pad_width: Width of the padding
        :return: Padded grid
        """
        self.padded_grid = np.pad(grid, pad_width, mode='reflect')
        return self.padded_grid

    def mean_filter(self, grid, n):
        """
        Apply mean filter.

        :param n: Size of the kernel (n x n)
        :return: Filtered grid
        """
        kernel = np.ones((n, n)) / (n * n)
        new_grid=convolve(grid, kernel, mode='reflect')

        return new_grid

    def median_filter(self, grid, n):
        """
        Apply median filter.

        :param n: Size of the kernel (n x n)
        :return: Filtered grid
        """
        return median_filter(grid, size=(n, n), mode='reflect')

    def gaussian_filter(self, grid, sigma):
        """
        Apply Gaussian filter.

        :param sigma: Standard deviation for Gaussian kernel
        :return: Filtered grid
        """
        return gaussian_filter(grid, sigma=sigma, mode='reflect')

    def directional_filter(self,grid, direction, n=3):
        """
        Apply directional filter (NE, N, NW, W, SW, S, SE, E).

        :param direction: Direction of the filter ('NE', 'N', 'NW', 'W', 'SW', 'S', 'SE', 'E')
        :param n: Size of the kernel (n x n, default is 3x3)
        :return: Filtered grid
        """
        direction_kernels = {
            'N': np.array([[-1, -1, -1], [0, 0, 0], [1, 1, 1]]),
            'S': np.array([[1, 1, 1], [0, 0, 0], [-1, -1, -1]]),
            'E': np.array([[-1, 0, 1], [-1, 0, 1], [-1, 0, 1]]),
            'W': np.array([[1, 0, -1], [1, 0, -1], [1, 0, -1]]),
            'NE': np.array([[0, -1, -1], [1, 0, -1], [1, 1, 0]]),
            'NW': np.array([[-1, -1, 0], [-1, 0, 1], [0, 1, 1]]),
            'SE': np.array([[0, 1, 1], [-1, 0, 1], [-1, -1, 0]]),
            'SW': np.array([[1, 1, 0], [1, 0, -1], [0, -1, -1]]),
        }

        if direction not in direction_kernels:
            raise ValueError(f"Invalid direction '{direction}'. Must be one of {list(direction_kernels.keys())}.")

        kernel = direction_kernels[direction]
        kernel_size = kernel.shape[0]

        if kernel_size != n:
            kernel = np.pad(kernel, ((n - kernel_size) // 2, (n - kernel_size) // 2), mode='constant')

        return convolve(grid, kernel, mode='reflect')

    def sun_shading_filter(self, grid, azimuth, zenith):
        """
        Calculate sun-shading based on user-defined azimuth and zenith angles.

        :param azimuth: Sun azimuth angle (degrees, 0 is North and increases clockwise)
        :param zenith: Sun zenith angle (degrees from vertical)
        :return: Shaded grid
        """

        selected_layer=QgsProject.instance().mapLayersByName(self.localGridName)[0]
        if(selected_layer.isValid()):
            crs = selected_layer.crs()
            if (crs.isGeographic()):
                # Convert azimuth and zenith to radians
                avg_spacing = (self.dx + self.dy) / 2
                z_scaling_factor = avg_spacing * 111_000.0  # Approx. conversion factor for degrees to meters

                # Apply scaling to Z values
                grid_scaled = grid * z_scaling_factor
            else:
                grid_scaled = grid
        else:
            print("notValid")

        # Convert azimuth and altitude to radians
        azimuth_rad = np.radians(azimuth)
        altitude_rad = np.radians(zenith)

        # Compute sun shading
        x, y = np.gradient(grid_scaled, self.dx, self.dy)
        slope = np.pi / 2 - np.arctan(np.sqrt(x**2 + y**2))
        aspect = np.arctan2(-x, y)

        shaded = (
            np.sin(altitude_rad) * np.sin(slope)
            + np.cos(altitude_rad) * np.cos(slope) * np.cos(azimuth_rad - aspect)
        )

        return shaded