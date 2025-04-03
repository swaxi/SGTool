import numpy as np
from scipy import stats
from scipy.ndimage import generic_filter, gaussian_filter


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

    def classify_terrain_with_cell_size(
        self,
        cell_size_x,
        cell_size_y,
        curvature_threshold=0.01,
        slope_threshold=45,
        window_size=3,
        sigma=1.0,
    ):
        """
        Proper wrapper function for terrain classification that incorporates cell size.

        Parameters:
        -----------
        cell_size_x : float
            Cell size in the x direction
        cell_size_y : float
            Cell size in the y direction
        curvature_threshold : float, optional
            Curvature threshold for classification (default: 0.01)
        slope_threshold : float, optional
            Slope threshold in degrees for cliff classification (default: 45)
        window_size : int, optional
            Window size for gradient calculation (default: 3)
        sigma : float, optional
            Sigma for optional smoothing (default: 1.0)

        Returns:
        --------
        numpy.ndarray
            Classified array where:
            -1 = concave up
            0 = flat
            1 = convex up
            2 = cliff
        """
        dem_array = np.float32(self.grid)

        # Apply optional smoothing to reduce noise
        if sigma > 0:
            dem_smooth = gaussian_filter(dem_array, sigma=sigma)
        else:
            dem_smooth = dem_array.copy()

        # Calculate gradients with correct cell sizes
        # Use window_size parameter to control gradient calculation precision
        # Larger window means smoother gradients but less detail
        if window_size > 1:
            # For larger windows, use a Sobel-like approach with window_size
            from scipy.signal import convolve2d

            # Create kernels of the right size
            kernel_y = np.zeros((window_size, window_size))
            kernel_x = np.zeros((window_size, window_size))

            # Fill the kernels for gradient calculation
            half_win = window_size // 2
            for i in range(window_size):
                for j in range(window_size):
                    y_dist = i - half_win
                    x_dist = j - half_win
                    if abs(y_dist) > 0:
                        kernel_y[i, j] = y_dist / (
                            y_dist**2 + x_dist**2 if x_dist != 0 else abs(y_dist)
                        )
                    if abs(x_dist) > 0:
                        kernel_x[i, j] = x_dist / (
                            y_dist**2 + x_dist**2 if y_dist != 0 else abs(x_dist)
                        )

            # Normalize kernels
            kernel_y = (
                kernel_y / np.sum(np.abs(kernel_y))
                if np.sum(np.abs(kernel_y)) > 0
                else kernel_y
            )
            kernel_x = (
                kernel_x / np.sum(np.abs(kernel_x))
                if np.sum(np.abs(kernel_x)) > 0
                else kernel_x
            )

            # Calculate gradients with convolution
            dy = convolve2d(dem_smooth, kernel_y, mode="same") / cell_size_y
            dx = convolve2d(dem_smooth, kernel_x, mode="same") / cell_size_x
        else:
            # Use simple numpy gradient for window_size = 1
            dy, dx = np.gradient(dem_smooth, cell_size_y, cell_size_x)

        # Calculate slope in degrees
        slope = np.degrees(np.arctan(np.sqrt(dx**2 + dy**2)))

        # Calculate second derivatives
        dxx, dxy = np.gradient(dx, cell_size_y, cell_size_x)
        dyx, dyy = np.gradient(dy, cell_size_y, cell_size_x)

        # Calculate profile curvature (in direction of steepest slope)
        p = dx**2 + dy**2
        q = p + 1

        # Avoid division by zero
        p_safe = np.where(p > 0, p, np.finfo(float).eps)

        # Profile curvature formula
        profile_curv = ((dx**2 * dxx) + (2 * dx * dy * dxy) + (dy**2 * dyy)) / (
            p_safe * np.sqrt(q**3)
        )

        # Classify based on thresholds
        classified = np.zeros_like(dem_array, dtype=np.int8)

        # First identify cliffs using direct slope threshold
        # We'll use this to check if linear cliffs are working
        basic_cliff_mask = slope >= slope_threshold

        # Get linear cliffs
        linear_cliffs = self.detect_linear_cliffs(
            slope, min_length=5, slope_threshold=slope_threshold, connectivity=2
        )

        # Use the linear cliffs as our cliff mask
        cliff_mask = linear_cliffs > 0

        # Check if we got any cliffs
        if np.sum(cliff_mask) == 0:
            # Fall back to basic slope threshold if no linear cliffs detected
            cliff_mask = basic_cliff_mask
            print(f"No linear cliffs found. Falling back to basic slope threshold.")
            print(
                f"Basic slope threshold found {np.sum(basic_cliff_mask)} cliff cells."
            )
        else:
            print(f"Found {np.sum(cliff_mask)} linear cliff cells.")

        classified[cliff_mask] = 2  # Cliffs

        # Then classify non-cliff areas based on curvature
        non_cliff_mask = ~cliff_mask

        # Use more appropriate curvature thresholds - these may need adjustment
        # Try a smaller threshold if you're getting no concave/convex features
        concave_mask = non_cliff_mask & (profile_curv > curvature_threshold)
        convex_mask = non_cliff_mask & (profile_curv < -curvature_threshold)
        flat_mask = (
            non_cliff_mask
            & (profile_curv >= -curvature_threshold)
            & (profile_curv <= curvature_threshold)
        )

        classified[concave_mask] = -1  # Concave up
        classified[convex_mask] = 1  # Convex up
        classified[flat_mask] = 0  # Flat areas explicitly set

        # Add debug info
        print(f"Curvature range: {np.min(profile_curv)} to {np.max(profile_curv)}")
        print(f"Number of concave cells: {np.sum(concave_mask)}")
        print(f"Number of convex cells: {np.sum(convex_mask)}")
        print(f"Number of flat cells: {np.sum(flat_mask)}")
        print(f"Number of cliff cells: {np.sum(cliff_mask)}")

        return classified

    def detect_linear_cliffs(
        self, slope_array, min_length=5, slope_threshold=45, connectivity=2
    ):
        """
        Enhanced cliff detection that finds linear cliff features.
        Requires skimage for connected component analysis.

        Parameters:
        -----------
        slope_array : numpy.ndarray
            Pre-calculated slope values in degrees
        min_length : int, optional
            Minimum length of feature to be classified as a cliff (default: 5)
        slope_threshold : float, optional
            Slope threshold in degrees (default: 45)
        connectivity : int, optional
            Connectivity for component analysis (1 or 2, default: 2)
            1 = 4-connectivity, 2 = 8-connectivity in 2D images

        Returns:
        --------
        numpy.ndarray
            Binary array where 1 = linear cliff feature, 0 = not cliff
        """
        try:
            from skimage import measure
            import numpy as np

            # Create binary mask of steep slopes
            steep_mask = slope_array >= slope_threshold

            # Check if we have any steep slopes at all
            if np.sum(steep_mask) == 0:
                print(f"No steep slopes found above threshold {slope_threshold}°")
                return np.zeros_like(slope_array, dtype=np.int8)

            # Find connected components - connectivity must be either 1 or 2 for 2D images
            labeled = measure.label(steep_mask, connectivity=connectivity)
            num_features = np.max(labeled)

            print(f"Found {num_features} connected steep regions")

            # If no regions found, return empty array
            if num_features == 0:
                return np.zeros_like(slope_array, dtype=np.int8)

            # Calculate properties of each feature
            props = measure.regionprops(labeled)

            # Prepare output array
            linear_cliffs = np.zeros_like(slope_array, dtype=np.int8)
            linear_count = 0

            for prop in props:
                # Skip small regions
                if prop.area < min_length:
                    continue

                # Check if feature is linear using eccentricity or axis ratio
                is_linear = False

                # Check for valid eccentricity (avoid potential NaN or division by zero)
                if hasattr(prop, "eccentricity") and prop.eccentricity > 0.8:
                    is_linear = True

                # Check axis ratio if we have both measurements
                if (
                    hasattr(prop, "major_axis_length")
                    and hasattr(prop, "minor_axis_length")
                    and prop.minor_axis_length > 0
                    and prop.major_axis_length / prop.minor_axis_length > 3
                ):
                    is_linear = True

                if is_linear:
                    linear_cliffs[labeled == prop.label] = 1
                    linear_count += 1

            print(
                f"Found {linear_count} linear features out of {num_features} steep regions"
            )
            return linear_cliffs

        except Exception as e:
            # Print the error for debugging
            print(f"Error in detect_linear_cliffs: {str(e)}")
            # Fall back to simple slope-based detection
            return (slope_array >= slope_threshold).astype(np.int8)
