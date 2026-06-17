import numpy as np
from scipy import stats
from scipy.ndimage import generic_filter, gaussian_filter, sobel, uniform_filter

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
        # Create a circular mask
        radius = window_size // 2
        y, x = np.ogrid[-radius : radius + 1, -radius : radius + 1]
        circular_mask = x * x + y * y <= radius * radius

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

            # Remove the circular mask from the window
            # This will create a circular window of valid values
            window_2d[~circular_mask] = np.nan

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
        nan_mask = np.isnan(dem_array)
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

        classified[nan_mask] = -15
        return classified

    def detect_linear_cliffs_skimage(
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
        from skimage import measure
        import numpy as np       
        try:

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

    def detect_linear_cliffs(
        self, slope_array, min_length=5, slope_threshold=45, connectivity=2
    ):
        """
        Enhanced cliff detection that finds linear cliff features.
        Uses scipy.ndimage for connected component analysis.

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
        from scipy import ndimage
        import numpy as np
        from collections import namedtuple
        
        try:
            # Create binary mask of steep slopes
            steep_mask = slope_array >= slope_threshold

            # Check if we have any steep slopes at all
            if np.sum(steep_mask) == 0:
                print(f"No steep slopes found above threshold {slope_threshold}°")
                return np.zeros_like(slope_array, dtype=np.int8)

            # Define structure for connectivity - 4-connectivity or 8-connectivity
            if connectivity == 1:
                # 4-connectivity: only adjacent pixels (not diagonals)
                structure = ndimage.generate_binary_structure(2, 1)
            else:
                # 8-connectivity: adjacent and diagonal pixels
                structure = ndimage.generate_binary_structure(2, 2)

            # Find connected components
            labeled, num_features = ndimage.label(steep_mask, structure=structure)

            print(f"Found {num_features} connected steep regions")

            # If no regions found, return empty array
            if num_features == 0:
                return np.zeros_like(slope_array, dtype=np.int8)

            # Prepare output array
            linear_cliffs = np.zeros_like(slope_array, dtype=np.int8)
            linear_count = 0

            # Process each labeled region
            for label_value in range(1, num_features + 1):
                # Extract region mask
                region_mask = labeled == label_value
                region_area = np.sum(region_mask)
                
                # Skip small regions
                if region_area < min_length:
                    continue
                    
                # Get region coordinates
                coords = np.argwhere(region_mask)
                
                # Check if feature is linear using basic shape analysis
                is_linear = False
                
                # Calculate region properties similar to skimage.measure.regionprops
                # (1) Calculate covariance matrix for the coordinates
                y_coords, x_coords = coords[:, 0], coords[:, 1]
                y_mean, x_mean = np.mean(y_coords), np.mean(x_coords)
                y_centered, x_centered = y_coords - y_mean, x_coords - x_mean
                
                # Covariance matrix
                cov_xx = np.sum(x_centered * x_centered) / region_area
                cov_yy = np.sum(y_centered * y_centered) / region_area
                cov_xy = np.sum(x_centered * y_centered) / region_area
                
                cov_matrix = np.array([[cov_yy, cov_xy], [cov_xy, cov_xx]])
                
                # Get eigenvalues to determine shape
                try:
                    eigenvalues, _ = np.linalg.eigh(cov_matrix)
                    # Ensure eigenvalues are positive and sorted
                    eigenvalues = np.clip(eigenvalues, 1e-10, None)  # Avoid division by zero
                    eigenvalues = np.sort(eigenvalues)[::-1]  # Sort in descending order
                    
                    # Calculate axis lengths (similar to major_axis_length and minor_axis_length)
                    major_axis_length = 4 * np.sqrt(eigenvalues[0])
                    minor_axis_length = 4 * np.sqrt(eigenvalues[1])
                    
                    # Calculate eccentricity
                    if major_axis_length > 0:
                        eccentricity = np.sqrt(1 - (minor_axis_length / major_axis_length) ** 2)
                    else:
                        eccentricity = 0
                        
                    # Check if feature is linear using eccentricity or axis ratio
                    if eccentricity > 0.8:
                        is_linear = True
                        
                    # Check axis ratio
                    if (minor_axis_length > 0 and 
                        major_axis_length / minor_axis_length > 3):
                        is_linear = True
                        
                    if is_linear:
                        linear_cliffs[region_mask] = 1
                        linear_count += 1
                        
                except np.linalg.LinAlgError:
                    # If eigenvalue calculation fails, skip this region
                    continue

            print(
                f"Found {linear_count} linear features out of {num_features} steep regions"
            )
            return linear_cliffs

        except Exception as e:
            # Print the error for debugging
            print(f"Error in detect_linear_cliffs: {str(e)}")
            # Fall back to simple slope-based detection
            return (slope_array >= slope_threshold).astype(np.int8)


    def calculate_local_anisotropy(self, image, window_size=7, window_type='gaussian'):
        """
        Calculates the local anisotropy magnitude and dominant orientation map.
        
        Parameters:
        -----------
        image : ndarray
            The input 2D grayscale image.
        window_size : int
            The diameter/width of the local neighborhood in pixels (e.g., 5, 7, 15).
        window_type : str
            'box' for a uniform square window, or 'gaussian' for a weighted window.
            
        Returns:
        --------
        anisotropy_map : ndarray
            Values from 0.0 (perfectly isotropic/uniform) to 1.0 (highly directional).
        orientation_deg : ndarray
            Angles from -90 to +90 degrees pointing in the direction of maximum intensity change.
        """
        # Preserve NaN locations; replace with 0 so scipy filters don't propagate NaN
        nan_mask = np.isnan(image)
        if nan_mask.any():
            image = np.where(nan_mask, 0.0, image)

        # 1. Compute local image gradients
        Ix = sobel(image, axis=1)
        Iy = sobel(image, axis=0)
        
        # 2. Compute raw tensor components
        Ixx = Ix ** 2
        Ixy = Ix * Iy
        Iyy = Iy ** 2
        
        # 3. Integrate over the specified local window size
        if window_type.lower() == 'box':
            J11 = uniform_filter(Ixx, size=window_size)
            J12 = uniform_filter(Ixy, size=window_size)
            J22 = uniform_filter(Iyy, size=window_size)
            
        elif window_type.lower() == 'gaussian':
            # Standard deviation covering roughly 99% of the requested window size
            sigma_value = window_size / 6.0  
            J11 = gaussian_filter(Ixx, sigma=sigma_value)
            J12 = gaussian_filter(Ixy, sigma=sigma_value)
            J22 = gaussian_filter(Iyy, sigma=sigma_value)
            
        else:
            raise ValueError("window_type must be either 'box' or 'gaussian'")
        
        # 4. Analytic eigenvalues for a 2x2 symmetric matrix
        trace = J11 + J22
        diff = np.sqrt((J11 - J22)**2 + 4 * J12**2)
        
        lambda1 = (trace + diff) / 2.0
        lambda2 = (trace - diff) / 2.0
        
        # 5. Calculate final Anisotropy and Orientation maps
        # Saliency = λ1 - λ2 (raw directional gradient energy).
        # Unlike coherence (λ1-λ2)/(λ1+λ2+ε), saliency is zero for flat areas
        # (both eigenvalues near zero) AND for isotropic texture (λ1≈λ2).
        # Coherence gave ~1 in background because it normalises by total energy,
        # amplifying even tiny directional noise.
        saliency = lambda1 - lambda2  # always >= 0 since λ1 >= λ2
        p99 = float(np.nanpercentile(saliency, 99)) if np.any(~np.isnan(saliency)) else 1.0
        anisotropy_map = np.clip(saliency / (p99 + 1e-10), 0.0, 1.0)
        
        orientation_rad = 0.5 * np.arctan2(2 * J12, J11 - J22)
        orientation_deg = np.degrees(orientation_rad)%180

        if nan_mask.any():
            anisotropy_map[nan_mask] = np.nan
            orientation_deg[nan_mask] = np.nan

        return anisotropy_map, orientation_deg



    def calculate_anisotropy_chain_length(self, anisotropy_map, orientation_deg,
                                          aniso_threshold=0.3, angle_tolerance=22.5,
                                          search_radius=2):
        """
        Scores each active pixel by the size of the connected component it belongs to.

        Two active pixels are linked when they are within search_radius pixels of each
        other (Euclidean disk) AND their line orientations differ by at most
        angle_tolerance degrees.  The direction of the link is unconstrained — pixels
        side-by-side across the width of a multi-pixel-wide feature are linked just as
        readily as pixels ahead/behind along it.  Connected components are found with
        Union-Find; every member receives the same score (component size).

        Parameters
        ----------
        anisotropy_map  : ndarray  Values 0–1
        orientation_deg : ndarray  0–180° line orientation
        aniso_threshold : float    Minimum anisotropy
        angle_tolerance : float    Max orientation difference to accept a link (degrees)
        search_radius   : int      Disk radius in pixels (default 2)

        Returns
        -------
        chain_length_map : ndarray  NaN outside active mask, component size inside.
        """
        nr, nc = anisotropy_map.shape
        orient = orientation_deg % 180

        active = (
            ~np.isnan(anisotropy_map) & ~np.isnan(orient) &
            (anisotropy_map >= aniso_threshold)
        )

        def ang_diff(a, b):
            d = abs(float(a) - float(b)) % 180
            return min(d, 180.0 - d)

        # Pre-compute disk offsets once
        r2max = search_radius * search_radius
        disk = [
            (dr, dc)
            for dr in range(-search_radius, search_radius + 1)
            for dc in range(-search_radius, search_radius + 1)
            if 0 < dr * dr + dc * dc <= r2max
        ]

        active_rows, active_cols = np.where(active)
        active_list = list(zip(active_rows.tolist(), active_cols.tolist()))

        # Build undirected edges: any two active pixels within disk distance with
        # compatible orientation.  Pixels across the width of a wide feature are
        # connected this way, not just pixels along it.
        edges = set()
        for r, c in active_list:
            theta = float(orient[r, c])
            for dr, dc in disk:
                r2, c2 = r + dr, c + dc
                if not (0 <= r2 < nr and 0 <= c2 < nc):
                    continue
                if not active[r2, c2]:
                    continue
                if ang_diff(theta, float(orient[r2, c2])) <= angle_tolerance:
                    edges.add((min((r, c), (r2, c2)), max((r, c), (r2, c2))))

        # Union-Find
        parent = {p: p for p in active_list}
        size   = {p: 1 for p in active_list}

        def find(x):
            root = x
            while parent[root] != root:
                root = parent[root]
            while parent[x] != root:
                parent[x], x = root, parent[x]
            return root

        def union(a, b):
            ra, rb = find(a), find(b)
            if ra == rb:
                return
            if size[ra] < size[rb]:
                ra, rb = rb, ra
            parent[rb] = ra
            size[ra] += size[rb]

        for a, b in edges:
            union(a, b)

        result = np.full((nr, nc), np.nan)
        for r, c in active_list:
            result[r, c] = float(size[find((r, c))])

        return result

    def calculate_anisotropy_streamline_length(self, anisotropy_map, orientation_deg,
                                               aniso_threshold=0.3, angle_tolerance=22.5,
                                               step_size=0.5, max_steps=200):
        """
        For each active pixel traces forward and backward along the orientation field
        using sub-pixel bilinear interpolation.  Accumulates path length (in pixels)
        while anisotropy stays above threshold and the step-to-step orientation change
        stays within angle_tolerance.  Handles curved features naturally.

        Parameters
        ----------
        anisotropy_map : ndarray  Values 0–1 (from calculate_local_anisotropy)
        orientation_deg : ndarray  0–360° (from calculate_local_anisotropy)
        aniso_threshold : float   Active pixel threshold (default 0.3)
        angle_tolerance : float   Max step-to-step orientation change in degrees (default 22.5)
        step_size : float         Sub-pixel step length in pixels (default 0.5)
        max_steps : int           Maximum steps per direction (default 200)

        Returns
        -------
        streamline_length_map : ndarray  Total forward+backward path length in pixels.
        """
        nr, nc = anisotropy_map.shape
        orient = orientation_deg % 180

        aniso_clean = np.where(np.isnan(anisotropy_map), 0.0, anisotropy_map)
        orient_clean = np.where(np.isnan(orient), 0.0, orient)

        active = (
            ~np.isnan(anisotropy_map) & ~np.isnan(orient) &
            (anisotropy_map >= aniso_threshold)
        )

        def _bilinear(arr, r, c):
            """Bilinear interpolation with clamped boundaries."""
            r0 = max(0, min(nr - 2, int(r)))
            c0 = max(0, min(nc - 2, int(c)))
            wr = max(0.0, min(1.0, r - r0))
            wc = max(0.0, min(1.0, c - c0))
            return (arr[r0,   c0]   * (1-wr) * (1-wc) +
                    arr[r0,   c0+1] * (1-wr) *   wc   +
                    arr[r0+1, c0]   *   wr   * (1-wc) +
                    arr[r0+1, c0+1] *   wr   *   wc)

        def _interp_orient(r, c):
            """Orientation interpolation using double-angle trick (handles 0/180 wrap)."""
            r0 = max(0, min(nr - 2, int(r)))
            c0 = max(0, min(nc - 2, int(c)))
            wr = max(0.0, min(1.0, r - r0))
            wc = max(0.0, min(1.0, c - c0))
            # Double angle → interpolate as unit vector → halve back
            t00 = np.radians(orient_clean[r0,   c0]   * 2)
            t01 = np.radians(orient_clean[r0,   c0+1] * 2)
            t10 = np.radians(orient_clean[r0+1, c0]   * 2)
            t11 = np.radians(orient_clean[r0+1, c0+1] * 2)
            s = (np.sin(t00)*(1-wr)*(1-wc) + np.sin(t01)*(1-wr)*wc +
                 np.sin(t10)*wr*(1-wc)      + np.sin(t11)*wr*wc)
            co = (np.cos(t00)*(1-wr)*(1-wc) + np.cos(t01)*(1-wr)*wc +
                  np.cos(t10)*wr*(1-wc)      + np.cos(t11)*wr*wc)
            return (np.degrees(np.arctan2(s, co)) / 2.0) % 180.0

        def _trace(r0, c0, direction):
            r, c = float(r0), float(c0)
            theta_prev = float(orient_clean[r0, c0])
            length = 0.0
            for _ in range(max_steps):
                if not (0.0 <= r < nr and 0.0 <= c < nc):
                    break
                aniso_cur = _bilinear(aniso_clean, r, c)
                if aniso_cur < aniso_threshold:
                    break
                theta_cur = _interp_orient(r, c)
                d = abs(theta_cur - theta_prev) % 180
                if min(d, 180.0 - d) > angle_tolerance:
                    break
                theta_rad = np.radians(theta_cur)
                # Convention: theta=0° → East (dc+), theta=90° → North (dr-)
                r += direction * (-np.sin(theta_rad)) * step_size
                c += direction *   np.cos(theta_rad)  * step_size
                length += step_size
                theta_prev = theta_cur
            return length

        result = np.full((nr, nc), np.nan)
        active_rows, active_cols = np.where(active)
        for r0, c0 in zip(active_rows.tolist(), active_cols.tolist()):
            result[r0, c0] = _trace(r0, c0, 1) + _trace(r0, c0, -1)

        return result

    def calculate_misorientation_cosines(self, map1_deg, map2_deg):
        """
        Calculates misorientation using the arccos of the dot product of direction cosines.
        Accounts for 180-degree line symmetry by using the absolute dot product.
        
        Parameters:
        -----------
        map1_deg : ndarray
            First orientation map in degrees.
        map2_deg : ndarray
            Second independent orientation map in degrees.
            
        Returns:
        --------
        misorientation_deg : ndarray
            The absolute angular difference map strictly bounded between 0° and 90°.
        """
        # 1. Convert input angles from degrees to radians
        rad1 = np.radians(map1_deg)
        rad2 = np.radians(map2_deg)
        
        # 2. Compute 2D unit vector components (Direction Cosines)
        # Vector 1: [u1, v1]
        u1, v1 = np.cos(rad1), np.sin(rad1)
        # Vector 2: [u2, v2]
        u2, v2 = np.cos(rad2), np.sin(rad2)
        
        # 3. Calculate the dot product of the vector fields
        dot_product = u1 * u2 + v1 * v2
        
        # 4. Take absolute value to enforce 180-degree line/structural symmetry
        # Clip to [-1.0, 1.0] to prevent floating-point inaccuracies from breaking arccos
        abs_dot = np.clip(np.abs(dot_product), -1.0, 1.0)
        
        # 5. Compute the final misorientation angle
        misorientation_rad = np.arccos(abs_dot)
        misorientation_deg = np.degrees(misorientation_rad)
        
        return misorientation_deg


