import numpy as np
from scipy.ndimage import convolve, median_filter, gaussian_filter
from PyQt5.QtGui import QValidator


class SG_Util:
    def __init__(self, grid):
        """
        Initialize the Utils with a grid.

        :param grid: 2D numpy array representing the input grid
        """
        self.grid = np.array(grid, dtype=float)

    def Threshold2Nan(
        self, grid, condition, above_threshold_value, below_threshold_value
    ):
        # Print types for debugging

        # Ensure grid is a floating-point array
        if not np.issubdtype(grid.dtype, np.floating):
            grid = grid.astype(float)

        # Process the raster using numpy
        if condition == "above":
            # Set values above the threshold to NaN
            grid2 = np.where(grid > above_threshold_value, np.nan, grid)
        elif condition == "below":
            # Set values below the threshold to NaN
            grid2 = np.where(grid < below_threshold_value, np.nan, grid)
        elif condition == "between":
            # Set values outside the thresholds to NaN
            grid2 = np.where(
                (grid < above_threshold_value) & (grid > below_threshold_value),
                np.nan,
                grid,
            )
        else:
            raise ValueError("Invalid condition. Use 'above', 'below', or 'both'.")

        return grid2
