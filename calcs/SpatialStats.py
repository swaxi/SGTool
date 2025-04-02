import numpy as np
from scipy import stats
from scipy.ndimage import generic_filter


class SpatialStats:
    def __init__(self, grid):
        """
        Initialize the class SpatialStats filter:

        :param grid: 2D numpy array representing the input grid
        """
        self.grid = np.array(grid, dtype=float)

    def calculate_windowed_stats(self, window_size=3, stat_type="variance"):
        """
        Calculate windowed statistics on a 2D raster grid.

        Parameters:
        -----------
        grid : numpy.ndarray
            2D array representing the raster grid, may contain NaN values
        window_size : int
            Size of the window (square) for calculation, must be odd
        stat_type : str
            Type of statistic to calculate. Options:
            'variance', 'std', 'skewness', 'kurtosis', 'mean', 'median', 'min', 'max'

        Returns:
        --------
        numpy.ndarray
            2D array with calculated statistics, same dimensions as input grid
        """
        if window_size % 2 == 0:
            raise ValueError("Window size must be odd")

        # Create a padded version of the grid to handle edges
        pad_size = window_size // 2
        padded_grid = np.pad(
            self.grid, pad_size, mode="constant", constant_values=np.nan
        )

        # Create a mask for NaN values
        mask = np.isnan(self.grid)

        # Function to calculate statistics on a window, ignoring NaNs
        def window_stat(window):
            # Reshape the window to a 2D array
            window_2d = window.reshape(window_size, window_size)
            # Remove NaN values
            valid_values = window_2d[~np.isnan(window_2d)]

            # Return NaN if not enough valid values in the window
            if (
                len(valid_values) < 3
            ):  # Require at least 3 values for meaningful statistics
                return np.nan

            if stat_type == "variance":
                return np.var(valid_values)
            elif stat_type == "std":
                return np.std(valid_values)
            elif stat_type == "skewness":
                return stats.skew(valid_values)
            elif stat_type == "kurtosis":
                return stats.kurtosis(valid_values)
            elif stat_type == "mean":
                return np.mean(valid_values)
            elif stat_type == "median":
                return np.median(valid_values)
            elif stat_type == "min":
                return np.min(valid_values)
            elif stat_type == "max":
                return np.max(valid_values)
            else:
                return np.nan

        # Apply the function to the padded grid
        result = generic_filter(padded_grid, window_stat, size=window_size)

        # Crop the result to match the original grid dimensions
        result = result[pad_size:-pad_size, pad_size:-pad_size]

        # Make sure NaN areas in the original grid remain NaN in the result
        result[mask] = np.nan

        return result
