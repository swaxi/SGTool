import numpy as np
import pywt
import matplotlib.pyplot as plt
from scipy import ndimage
from skimage import feature
from skimage.feature import peak_local_max
from skimage.morphology import skeletonize
from scipy import stats


from scipy import signal, ndimage
import networkx as nx
from matplotlib.collections import LineCollection
from matplotlib import cm

class WTMM:
    """Wavelet Transform Modulus Maxima Implementation."""

    def __init__(self):
        pass

    def wtmm_2d(self,image, wavelet='mexh', scales=None, num_scales=5, 
                threshold_rel=0.2, min_distance=5, chaining_threshold=2.0):
        """
        Perform 2D Wavelet Transform Modulus Maxima analysis on an image.
        
        Parameters:
        -----------
        image : numpy.ndarray
            2D input image (grayscale)
        wavelet : str or pywt.Wavelet
            Wavelet to use (default: 'mexh' - Mexican Hat wavelet)
        scales : list or None
            Specific scales to use for the wavelet transform
            If None, generates logarithmically spaced scales
        num_scales : int
            Number of scales to use if scales is None
        threshold_rel : float
            Relative threshold for maxima detection (0.0 to 1.0)
        min_distance : int
            Minimum distance between detected maxima
        chaining_threshold : float
            Maximum distance (in pixels) for chaining maxima across scales
            
        Returns:
        --------
        dict: A dictionary containing:
            'scales': The scales used
            'wavelet_transforms': Wavelet transforms at each scale
            'modulus': Wavelet transform modulus at each scale
            'maxima': Maxima points at each scale
            'maxima_lines': Chained maxima across scales
            'wtmm_image': Visualization of WTMM result
        """
        # Input validation
        if len(image.shape) > 2:
            # Convert RGB to grayscale if needed
            if image.shape[2] == 3 or image.shape[2] == 4:
                image = np.mean(image[:,:,:3], axis=2)
            else:
                raise ValueError("Input must be a 2D image or RGB/RGBA")
        
        # Normalize image to [0, 1]
        image = (image - np.min(image)) / (np.max(image) - np.min(image))
        
        # Generate scales if not provided
        if scales is None:
            # Logarithmically spaced scales
            scales = np.logspace(0, np.log10(min(image.shape) / 4), num_scales)
        
        # Storage for results
        results = {
            'scales': scales,
            'wavelet_transforms': [],
            'modulus': [],
            'phase': [],
            'maxima': [],
            'directions': []
        }
        
        # Apply wavelet transform at each scale
        for scale in scales:
            # Apply stationary wavelet transform (SWT) with proper scaling
            # We use SWT because it's shift-invariant, which is important for edge detection
            sw_coeffs = pywt.swt2(image, wavelet, level=1, start_level=0)[0]
            
            # Extract horizontal and vertical detail coefficients
            # These represent derivatives in x and y directions
            coeffs_h = sw_coeffs[1][0] * scale  # Horizontal details (x derivative)
            coeffs_v = sw_coeffs[1][1] * scale  # Vertical details (y derivative)
            
            # Compute the modulus (gradient magnitude)
            modulus = np.sqrt(coeffs_h**2 + coeffs_v**2)
            
            # Compute the phase (gradient direction)
            phase = np.arctan2(coeffs_v, coeffs_h)
            
            # Find local maxima of the modulus
            maxima_coords = peak_local_max(modulus, 
                                        min_distance=min_distance,
                                        threshold_rel=threshold_rel,
                                        exclude_border=False)
            
            # Store results for this scale
            results['wavelet_transforms'].append((coeffs_h, coeffs_v))
            results['modulus'].append(modulus)
            results['phase'].append(phase)
            results['maxima'].append(maxima_coords)
            
            # Compute gradient direction at each maximum
            directions = []
            for y, x in maxima_coords:
                # Gradient direction at this point (perpendicular to wavelet maxima direction)
                angle = phase[y, x] + np.pi/2
                directions.append(angle)
            
            results['directions'].append(np.array(directions))
        
        # Chain maxima across scales to create maxima lines
        maxima_lines = self.chain_maxima_across_scales(
            results['maxima'], 
            results['directions'], 
            scales,
            chaining_threshold
        )
        
        results['maxima_lines'] = maxima_lines
        
        # Create visualization of WTMM result
        wtmm_image = np.zeros_like(image)
        for line in maxima_lines:
            if len(line) > 1:  # Only consider lines spanning multiple scales
                for point in line:
                    y, x = point[:2]  # Extract coordinates
                    if 0 <= y < image.shape[0] and 0 <= x < image.shape[1]:
                        wtmm_image[int(y), int(x)] = 1
        
        # Skeletonize the result for cleaner lines
        wtmm_image = skeletonize(wtmm_image)
        results['wtmm_image'] = wtmm_image
        
        return results

    def chain_maxima_across_scales(self,maxima_list, directions_list, scales, threshold=2.0):
        """
        Chain maxima points across scales to create maxima lines.
        
        Parameters:
        -----------
        maxima_list : list of arrays
            List of arrays containing maxima coordinates at each scale
        directions_list : list of arrays
            List of arrays containing gradient directions at each maximum
        scales : array-like
            Scales used for wavelet transform
        threshold : float
            Maximum distance for chaining maxima
            
        Returns:
        --------
        list: List of maxima lines, where each line is a list of 
            (y, x, scale_index, direction) tuples
        """
        maxima_lines = []
        
        # Start from the finest scale (smallest structures)
        finest_scale_idx = np.argmin(scales)
        
        # Initialize maxima lines from the finest scale
        for i, (y, x) in enumerate(maxima_list[finest_scale_idx]):
            direction = directions_list[finest_scale_idx][i]
            maxima_lines.append([(y, x, finest_scale_idx, direction)])
        
        # Order of scales from finest to coarsest
        scale_order = np.argsort(scales)
        
        # Chain maxima across consecutive scales
        for i in range(len(scale_order) - 1):
            current_scale_idx = scale_order[i]
            next_scale_idx = scale_order[i + 1]
            
            # Scale factor between current and next scale
            scale_factor = scales[next_scale_idx] / scales[current_scale_idx]
            
            # Adjust search threshold based on scale difference
            adjusted_threshold = threshold * scale_factor
            
            # For each existing line, try to extend to the next scale
            for line in maxima_lines:
                last_point = line[-1]
                last_y, last_x = last_point[0], last_point[1]
                last_direction = last_point[3]
                
                # Find closest maxima in the next scale
                best_distance = float('inf')
                best_match = None
                
                for j, (y, x) in enumerate(maxima_list[next_scale_idx]):
                    distance = np.sqrt((y - last_y)**2 + (x - last_x)**2)
                    
                    # Direction at the candidate point
                    direction = directions_list[next_scale_idx][j]
                    
                    # Check if the direction is consistent (allowing for some variation)
                    direction_diff = np.abs(np.mod(direction - last_direction + np.pi, 2*np.pi) - np.pi)
                    direction_consistent = direction_diff < np.pi/4  # 45 degrees tolerance
                    
                    if distance < adjusted_threshold and distance < best_distance and direction_consistent:
                        best_distance = distance
                        best_match = (y, x, next_scale_idx, direction)
                
                # Add the best match to the line if found
                if best_match is not None:
                    line.append(best_match)
        
        # Filter lines to keep only those spanning multiple scales
        filtered_lines = [line for line in maxima_lines if len(line) > 1]
        
        return filtered_lines


    def wtmm_1d(self,signal_data, wavelet='mexh', scales=None, num_scales=10, 
                threshold_rel=0.1, min_distance=3, chaining_threshold=2.0):
        """
        Perform 1D Wavelet Transform Modulus Maxima analysis.
        
        Parameters:
        -----------
        signal_data : numpy.ndarray
            1D input signal
        wavelet : str or pywt.Wavelet
            Wavelet to use (default: 'mexh' - Mexican Hat wavelet)
        scales : list or None
            Specific scales to use for the wavelet transform
            If None, generates logarithmically spaced scales
        num_scales : int
            Number of scales to use if scales is None
        threshold_rel : float
            Relative threshold for maxima detection (0.0 to 1.0)
        min_distance : int
            Minimum distance between detected maxima
        chaining_threshold : float
            Maximum distance for chaining maxima across scales
            
        Returns:
        --------
        dict: A dictionary containing:
            'scales': The scales used
            'cwt': Continuous wavelet transform coefficients
            'modulus': Absolute values of wavelet coefficients
            'maxima': Maxima points at each scale
            'maxima_lines': Chained maxima across scales
        """
        # Ensure signal is 1D
        signal_data = np.asarray(signal_data).flatten()
        
        # Normalize signal to [0, 1]
        signal_data = (signal_data - np.min(signal_data)) / (np.max(signal_data) - np.min(signal_data))
        
        # Generate scales if not provided
        if scales is None:
            # Logarithmically spaced scales
            scales = np.logspace(0, np.log10(len(signal_data) / 4), num_scales)
        
        # Storage for results
        results = {
            'scales': scales,
            'cwt': [],
            'modulus': [],
            'maxima': [],
            'maxima_values': []
        }
        
        # Apply continuous wavelet transform
        coeffs, frequencies = pywt.cwt(signal_data, scales, wavelet)
        results['cwt'] = coeffs
        
        # Compute the modulus (absolute value of coefficients)
        modulus = np.abs(coeffs)
        results['modulus'] = modulus
        
        # Find local maxima at each scale
        all_maxima = []
        all_maxima_values = []
        
        for i, scale_modulus in enumerate(modulus):
            # Find local maxima
            maxima_indices = peak_local_max(
                scale_modulus, 
                min_distance=min_distance,
                threshold_rel=threshold_rel,
                exclude_border=False,
                #indices=True
            ).flatten()
            
            # Store maxima coordinates and values
            all_maxima.append(maxima_indices)
            all_maxima_values.append(scale_modulus[maxima_indices])
        
        results['maxima'] = all_maxima
        results['maxima_values'] = all_maxima_values
        
        # Chain maxima across scales
        maxima_lines = self.chain_maxima_across_scales_1d(
            all_maxima, 
            all_maxima_values,
            scales,
            chaining_threshold
        )
        
        results['maxima_lines'] = maxima_lines
        
        return results

    def chain_maxima_across_scales_1d(self,maxima_list, values_list, scales, threshold=2.0):
        """
        Chain maxima points across scales to create maxima lines for 1D signal.
        
        Parameters:
        -----------
        maxima_list : list of arrays
            List of arrays containing maxima positions at each scale
        values_list : list of arrays
            List of arrays containing maxima values at each scale
        scales : array-like
            Scales used for wavelet transform
        threshold : float
            Maximum distance for chaining maxima
            
        Returns:
        --------
        list: List of maxima lines, where each line is a list of 
            (position, scale_index, value) tuples
        """
        maxima_lines = []
        
        # Start from the finest scale (smallest structures)
        finest_scale_idx = np.argmin(scales)
        
        # Initialize maxima lines from the finest scale
        for i, pos in enumerate(maxima_list[finest_scale_idx]):
            value = values_list[finest_scale_idx][i]
            maxima_lines.append([(pos, finest_scale_idx, value)])
        
        # Order of scales from finest to coarsest
        scale_order = np.argsort(scales)
        
        # Chain maxima across consecutive scales
        for i in range(len(scale_order) - 1):
            current_scale_idx = scale_order[i]
            next_scale_idx = scale_order[i + 1]
            
            # Scale factor between current and next scale
            scale_factor = scales[next_scale_idx] / scales[current_scale_idx]
            
            # Adjust search threshold based on scale difference
            adjusted_threshold = threshold * scale_factor
            
            # For each existing line, try to extend to the next scale
            for line in maxima_lines:
                last_point = line[-1]
                last_pos = last_point[0]
                
                # Find closest maxima in the next scale
                best_distance = float('inf')
                best_match = None
                
                for j, pos in enumerate(maxima_list[next_scale_idx]):
                    distance = np.abs(pos - last_pos)
                    value = values_list[next_scale_idx][j]
                    
                    if distance < adjusted_threshold and distance < best_distance:
                        best_distance = distance
                        best_match = (pos, next_scale_idx, value)
                
                # Add the best match to the line if found
                if best_match is not None:
                    line.append(best_match)
        
        # Filter lines to keep only those spanning multiple scales
        filtered_lines = [line for line in maxima_lines if len(line) > 1]
        
        return filtered_lines

    def visualize_wtmm_1d(self,signal_data,raster_layer_name, results, save_path=None):
        """
        Visualize the 1D WTMM analysis results.
        
        Parameters:
        -----------
        signal_data : numpy.ndarray
            Original input signal
        results : dict
            Results from wtmm_1d function
        save_path : str or None
            Path to save the visualization, if provided
        """
        # Create figure and axes
        fig = plt.figure(figsize=(12, 10))
        gs = plt.GridSpec(4, 1, height_ratios=[1, 2, 1, 3])
        
        # Plot the original signal
        ax1 = fig.add_subplot(gs[0])
        ax1.plot(signal_data, 'k-')
        ax1.set_title('Original Signal: {}'.format(raster_layer_name))
        ax1.set_ylim(np.min(signal_data) - 0.1, np.max(signal_data) + 0.1)
        ax1.set_xlim(0, len(signal_data))
        ax1.set_ylabel('Amplitude')
        
        # Plot CWT scalogram
        ax2 = fig.add_subplot(gs[1])
        extent = [0, len(signal_data), np.log2(results['scales'][-1]), np.log2(results['scales'][0])]
        ax2.imshow(results['modulus'], aspect='auto', extent=extent, cmap='viridis')
        ax2.set_title('CWT Scalogram')
        ax2.set_ylabel('$\\log_2(a)$')
        
        # Plot WTMM skeleton
        ax3 = fig.add_subplot(gs[2], sharex=ax2)
        ax3.set_ylim(np.log2(results['scales'][-1]), np.log2(results['scales'][0]))
        
        # Plot maxima lines
        for line in results['maxima_lines']:
            line_points = np.array(line)
            positions = line_points[:, 0]
            scale_indices = line_points[:, 1].astype(int)
            log_scales = np.log2(results['scales'][scale_indices])
            
            # Plot the line
            ax3.plot(positions, log_scales, 'r-', linewidth=1, alpha=0.7)
            # Plot points
            ax3.scatter(positions, log_scales, c='r', s=10, alpha=1)
        
        ax3.set_title('WTMM Skeleton')
        ax3.set_ylabel('$\\log_2(a)$')
        
        # Plot the cladogram
        ax4 = fig.add_subplot(gs[3], sharex=ax2)
        self.plot_wtmm_cladogram(signal_data,results, ax4)
        ax4.set_title('WTMM Cladogram')
        ax4.set_xlabel('Position')
        ax4.set_ylabel('$\\log_2(a)$')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300)
        
        return fig

    def plot_wtmm_cladogram(self,signal_data,results, ax=None):
        """
        Plot a cladogram representing the WTMM maxima lines.
        
        Parameters:
        -----------
        results : dict
            Results from wtmm_1d function
        ax : matplotlib.axes or None
            Axes to plot on, or create new if None
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 6))
        
        # Get scales and maxima lines
        scales = results['scales']
        maxima_lines = results['maxima_lines']
        
        # Set y-limits to log scales
        ax.set_ylim(np.log2(scales[-1]), np.log2(scales[0]))
        
        # Plot the cladogram
        max_val = 0
        for line in maxima_lines:
            if len(line) < 2:
                continue
                
            line_points = np.array(line)
            positions = line_points[:, 0]
            scale_indices = line_points[:, 1].astype(int)
            values = line_points[:, 2]
            log_scales = np.log2(scales[scale_indices])
            
            # Get max value for color scaling
            max_val = max(max_val, np.max(values))
            
            # Create line segments
            points = np.array([positions, log_scales]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            
            # Create a line collection with color mapped to the coefficient values
            lc = LineCollection(segments, cmap="viridis", norm=plt.Normalize(0, max_val))
            lc.set_array(values[:-1])
            lc.set_linewidth(3)
            
            # Add the collection to the plot
            line = ax.add_collection(lc)
        
        # Add a colorbar
        cbar = plt.colorbar(line, ax=ax,orientation='horizontal',shrink=0.5)
        cbar.set_label('Wavelet Coefficient Magnitude')
        
        # Create filled triangles at the peak positions
        heights = np.zeros(len(signal_data))
        
        for line in maxima_lines:
            if len(line) < 2:
                continue
                
            line_points = np.array(line)
            positions = line_points[:, 0]
            scale_indices = line_points[:, 1].astype(int)
            values = line_points[:, 2]
            log_scales = np.log2(scales[scale_indices])
            
            # Find the position at the smallest scale
            min_scale_idx = np.argmin(log_scales)
            pos = int(positions[min_scale_idx])
            
            if pos < len(heights):
                # Store the maximum scale reached by this feature
                heights[pos] = max(heights[pos], np.max(log_scales))
        
        # Plot filled triangles at positions where maxima lines start
        for pos, height in enumerate(heights):
            if height > 0:
                # Triangle
                y_top = height
                y_bottom = np.log2(scales[0])
                x_left = pos - 0.5
                x_right = pos + 0.5
                
                triangle = np.array([[pos, y_bottom], [x_left, y_top], [x_right, y_top]])
                ax.fill(triangle[:, 0], triangle[:, 1], 'r', alpha=0.3)
        
        return ax

    def generate_test_signal(self,length=512, num_periods=7):
        """Generate a test signal with multiple sinusoidal components."""
        x = np.linspace(0, num_periods * 2 * np.pi, length)
        
        # Combine signals of different frequencies
        signal = (
            np.sin(x) +                     # Low frequency
            0.5 * np.sin(3 * x) +           # Medium frequency
            0.25 * np.sin(25 * x) +          # Higher frequency
            0.125 * np.sin(10 * x) +        # High frequency detail
            0.0125 * np.sin(3 * x) +        # High frequency detail
            0.0625 * np.sin(66 * x) +        # High frequency detail
            0.00315 * np.sin(215 * x) +        # High frequency detail
            0.00125 * np.sin(44 * x) +        # High frequency detail
            0.006125 * np.sin(140 * x) +        # High frequency detail
            np.random.normal(0, 0.05, size=len(x))  # Noise
        )
        
        return signal

    def financial_cartoon(self,Iterations=10, Multifractal=1, noise_type=False, noise_level=1.0, plot=False):
        if Multifractal:
            turns = ((0.25, 0.5), (0.75, 0.25))
        else:
            turns = ((0.4, 0.6), (0.6, 0.4))
        first_turn, second_turn = turns
        ys = [0, 1]
        ts = [0, 1]

        if not noise_type:
            for i in range(0, Iterations + 1):

                j = 0
                while ts[j] < 1:
                    dt = ts[j + 1] - ts[j]
                    dy = ys[j + 1] - ys[j]

                    ts.insert(j + 1, ts[j] + first_turn[0] * dt)
                    ts.insert(j + 2, ts[j] + second_turn[0] * dt)
                    ys.insert(j + 1, ys[j] + first_turn[1] * dy)
                    ys.insert(j + 2, ys[j] + second_turn[1] * dy)

                    j += 3
        else:
            if noise_type == 'uniform':
                noise = np.random.rand
            elif noise_type == 'normal':
                noise = np.random.randn
            else:
                raise ValueError('Only normal and uniform accepted for noise at this time')

            for i in range(0, Iterations + 1):

                j = 0
                while ts[j] < 1:
                    dt = ts[j + 1] - ts[j]
                    dy = ys[j + 1] - ys[j]

                    # normalize the noise versus the current dt
                    n_a, n_b = (noise(2) * noise_level) * float(dy)

                    ts.insert(j + 1, ts[j] + first_turn[0] * dt)
                    ts.insert(j + 2, ts[j] + second_turn[0] * dt)
                    ys.insert(j + 1, ys[j] + n_a + first_turn[1] * dy)
                    ys.insert(j + 2, ys[j] + n_b + second_turn[1] * dy)

                    j += 3

        if plot:
            fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={'axisbg': '#EEEEEE', 'axisbelow': True})
            ax.grid(color='w', linewidth=2, linestyle='solid')
            ax.plot(ts, ys, color='b', alpha=0.4)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)

        return np.array(ts), np.array(ys)


    def example_data(self):
        data=[
        -2.25, -2.55, -2.44, -2.7, -2.74, -2.71, -2.82, -2.96, -3.12, -3.24, -3.26, -3.42, -3.61, -3.5, -3.85, -3.84, 
        -4.03, -4.14, -4.21, -4.19, -4.31, -4.34, -4.6, -4.75, -4.61, -4.71, -4.62, -4.58, -4.37, -4.37, -4.2, -4.13, 
        -4.15, -4.0, -3.95, -3.85, -3.79, -3.75, -3.53, -3.33, -3.18, -3.25, -3.32, -3.39, -3.39, -3.62, -3.77, -3.86, 
        -3.88, -4.01, -4.09, -4.05, -3.96, -3.74, -3.79, -3.73, -3.6, -3.38, -3.27, -3.06, -3.24, -2.86, -2.86, -2.94, 
        -2.72, -2.49, -2.44, -2.33, -2.42, -2.11, -2.08, -2.03, -1.92, -1.77, -1.63, -1.56, -1.57, -1.34, -1.4, -1.15, 
        -1.07, -1.07, -0.76, -0.65, -0.6, -0.46, -0.46, -0.36, -0.22, -0.2, -0.42, -0.45, -0.6, -0.66, -0.61, -0.83, 
        -0.81, -1.1, -1.14, -1.25, -1.4, -1.33, -1.47, -1.62, -1.84, -1.72, -1.83, -1.87, -2.17, -2.06, -2.16, -2.52, 
        -2.49, -2.7, -2.57, -2.68, -2.86, -3.08, -3.2, -3.16, -3.36, -3.5, -3.63, -3.58, -3.8, -3.73, -3.79, -4.09, 
        -3.95, -4.1, -4.17, -4.38, -4.6, -4.61, -4.63, -4.84, -4.96, -4.98, -5.14, -5.18, -5.0, -4.97, -4.91, -4.86, 
        -4.64, -4.52, -4.54, -4.38, -4.42, -4.32, -4.13, -4.08, -3.88, -3.8, -3.79, -3.45, -3.58, -3.53, -3.33, -3.16, 
        -3.13, -2.98, -3.09, -3.19, -3.19, -3.55, -3.63, -3.71, -3.66, -3.73, -3.97, -3.98, -4.19, -4.34, -4.23, -4.26, 
        -4.48, -4.3, -4.35, -4.31, -3.99, -4.08, -4.0, -3.85, -3.63, -3.53, -3.35, -3.31, -3.41, -3.23, -2.98, -2.96, 
        -2.84, -2.87, -2.65, -2.48, -2.45, -2.35, -2.39, -2.14, -1.96, -1.88, -1.97, -1.66, -1.56, -1.84, -2.0, -1.86, 
        -2.12, -2.33, -2.24, -2.27, -2.44, -2.45, -2.65, -2.86, -2.92, -3.12, -3.04, -3.22, -3.17, -3.31, -3.24, -2.93, 
        -3.0, -2.89, -2.61, -2.56, -2.6, -2.53, -2.36, -2.27, -2.17, -1.92, -1.97, -1.82, -1.56, -1.55, -1.63, -1.6, 
        -1.8, -1.9, -1.8, -1.91, -2.2, -2.08, -2.43, -2.39, -2.41, -2.69, -2.63, -2.85, -2.78, -3.14, -2.96, -3.27, 
        -3.27, -3.43, -3.4, -3.58, -3.78, -3.76, -3.98, -4.11, -4.14, -4.2, -4.36, -4.47, -4.49, -4.37, -4.32, -4.23, 
        -4.06, -4.1, -3.93, -3.94, -3.75, -3.5, -3.49, -3.37, -3.2, -3.07, -3.23, -3.12, -2.91, -2.69, -2.67, -2.6, 
        -2.5, -2.41, -2.24, -2.26, -2.06, -2.14, -1.96, -1.7, -1.77, -1.75, -1.48, -1.34, -1.4, -1.41, -1.61, -1.72, 
        -1.58, -1.68, -1.84, -1.89, -2.18, -2.32, -2.42, -2.26, -2.35
        ]
        data_array=np.array(data)
        return data_array

    # Example usage
    def example_usage_1d(self):
        """
        Demonstrate 1D WTMM on a sample signal.
        """
        # Generate a test signal
        signal_data = self.generate_test_signal(length=1024, num_periods=10)
        
        # Apply WTMM
        results = self.wtmm_1d(
            signal_data, 
            wavelet='mexh',  # Mexican hat wavelet
            num_scales=20,  # More scales for better visualization
            threshold_rel=0.1, 
            min_distance=5,
            chaining_threshold=2.0
        )
        
        # Visualize results
        fig = self.visualize_wtmm_1d(signal_data, results)
        plt.tight_layout()
        plt.show()
        
        return results, fig




    def visualize_wtmm(self,image, results, save_path=None):
        """
        Visualize the WTMM analysis results.
        
        Parameters:
        -----------
        image : numpy.ndarray
            Original input image
        results : dict
            Results from wtmm_2d function
        save_path : str or None
            Path to save the visualization, if provided
        """
        plt.figure(figsize=(15, 10))
        
        # Number of scales
        n_scales = len(results['scales'])
        
        # Plot original image
        plt.subplot(2, 3, 1)
        plt.imshow(image, cmap='gray')
        plt.title('Original Image')
        plt.axis('off')
        
        # Plot wavelet modulus at a few scales
        plot_indices = [0, n_scales//2, n_scales-1]  # First, middle, last
        for i, idx in enumerate(plot_indices):
            plt.subplot(2, 3, 2+i)
            plt.imshow(results['modulus'][idx], cmap='viridis')
            plt.title(f'Wavelet Modulus (Scale: {results["scales"][idx]:.2f})')
            
            # Plot maxima points
            maxima = results['maxima'][idx]
            if len(maxima) > 0:
                plt.scatter(maxima[:, 1], maxima[:, 0], c='r', s=5, alpha=0.7)
            
            plt.axis('off')
        
        # Plot the final WTMM result
        plt.subplot(2, 3, 5)
        plt.imshow(results['wtmm_image'], cmap='gray')
        plt.title('WTMM Skeleton')
        plt.axis('off')
        
        # Overlay WTMM skeleton on original image
        plt.subplot(2, 3, 6)
        plt.imshow(image, cmap='gray')
        masked_wtmm = np.ma.masked_where(results['wtmm_image'] == 0, results['wtmm_image'])
        plt.imshow(masked_wtmm, cmap='autumn', alpha=0.8)
        plt.title('WTMM Overlay')
        plt.axis('off')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path)
        
        plt.show()

    def visualize_wtmm_3d_plotly(self,image, results):
        """
        Create an interactive 3D visualization of the WTMM skeleton using Plotly.
        Offers more interactive features than the matplotlib version.
        
        Requires: pip install plotly
        """
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        
        # Create a figure
        fig = make_subplots(specs=[[{'type': 'surface'}]])
        
        # Add the image as a surface at the lowest scale
        min_log_scale = np.log2(results['scales'][0]) * 0.9  # Slightly below smallest scale
        y, x = np.mgrid[0:image.shape[0], 0:image.shape[1]]
        
        fig.add_trace(go.Surface(
            z=min_log_scale * np.ones_like(image),
            x=x,
            y=y,
            surfacecolor=image,
            colorscale='gray',
            opacity=0.3,
            showscale=False
        ))
        
        # Add all maxima lines as separate traces
        for i, line in enumerate(results['maxima_lines']):
            if len(line) > 1:  # Only consider lines spanning multiple scales
                # Extract coordinates and scales
                line_points = np.array(line)
                x_coords = line_points[:, 1]
                y_coords = line_points[:, 0]
                scale_indices = line_points[:, 2].astype(int)
                
                # Get the actual log2 scale values for the z-coordinate
                log_scales = np.log2(results['scales'][scale_indices])
                
                # Add the line as a 3D scatter trace
                fig.add_trace(go.Scatter3d(
                    x=x_coords,
                    y=y_coords,
                    z=log_scales,
                    mode='lines+markers',
                    line=dict(color='red', width=2),
                    marker=dict(size=3, color='red'),
                    name=f'Line {i}',
                    opacity=0.7
                ))
        
        # Update layout for better visualization
        log_scale_min = np.log2(min(results['scales']))
        log_scale_max = np.log2(max(results['scales']))
        
        fig.update_layout(
            title='WTMM Skeleton (3D)',
            scene=dict(
                xaxis_title='x',
                yaxis_title='y',
                zaxis_title='log₂(a)',
                zaxis=dict(range=[log_scale_min * 0.9, log_scale_max * 1.1]),
                aspectratio=dict(x=1, y=1, z=0.8)
            ),
            width=900,
            height=700,
            showlegend=False
        )
        
        return fig

    # Example usage
    def example_usage(self):
        """
        Demonstrate WTMM on a sample image.
        """
        # Create a sample image or load your own
        try:
            # Try to load a sample image from scikit-image
            from skimage import data
            image = data.camera()
        except:
            # Create a simple test image if loading fails
            image = np.zeros((256, 256))
            # Add some edges and features
            image[64:192, 64:192] = 1
            image = ndimage.gaussian_filter(image, sigma=3)
        
        # Apply WTMM
        results = self.wtmm_2d(
            image, 
            wavelet='haar',  # Using Haar wavelet which works well for edges
            num_scales=6, 
            threshold_rel=0.1, 
            min_distance=3,
            chaining_threshold=3.0
        )
        
        # Visualize results
        #visualize_wtmm(image, results)
        
        # Visualize 3D results
        fig_plotly = self.visualize_wtmm_3d_plotly(image, results)
        print("3D visualization ready. Use fig_plotly.show() to display.")
        fig_plotly.show()    

        return results

    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import stats

    def simplified_multifractal_spectrum(self,results, q_values=None):
        """
        A simplified implementation to compute D(h) vs h spectrum from WTMM results.
        
        Parameters:
        -----------
        results : dict
            Results from wtmm_1d function
        q_values : array-like or None
            Values of q for which to compute the spectrum
            If None, uses a range of values from -5 to 5
            
        Returns:
        --------
        tuple: (h_values, d_values) for plotting D(h) vs h
        """
        # Default q values if none provided
        if q_values is None:
            q_values = np.linspace(-5, 5, 41)
        
        # Extract maxima and scales from results
        scales = results['scales']
        maxima_values_list = results['maxima_values']
        
        # Compute partition function for each q and each scale
        tau_q = np.zeros(len(q_values))
        
        for q_idx, q in enumerate(q_values):
            # Compute log(Z(q,a)) vs log(a)
            log_scales = np.log(scales)
            log_Z_q = np.zeros_like(scales)
            
            for scale_idx, scale_maxima in enumerate(maxima_values_list):
                if len(scale_maxima) > 0:
                    # Sum of |W(a,b)|^q for all maxima at this scale
                    Z_q = np.sum(scale_maxima ** q)
                    log_Z_q[scale_idx] = np.log(Z_q) if Z_q > 0 else -np.inf
                else:
                    log_Z_q[scale_idx] = -np.inf
            
            # Filter out invalid values
            valid_indices = np.isfinite(log_Z_q)
            if np.sum(valid_indices) > 1:
                # Linear regression to find the slope (tau(q))
                slope, _, _, _, _ = stats.linregress(log_scales[valid_indices], log_Z_q[valid_indices])
                tau_q[q_idx] = slope
        
        # Compute h(q) and D(h) using numerical differentiation
        h_values = np.zeros(len(q_values))
        d_values = np.zeros(len(q_values))
        
        # Central differences for interior points
        for i in range(1, len(q_values) - 1):
            h_values[i] = (tau_q[i+1] - tau_q[i-1]) / (q_values[i+1] - q_values[i-1])
            d_values[i] = q_values[i] * h_values[i] - tau_q[i]
        
        # Forward/backward differences for endpoints
        h_values[0] = (tau_q[1] - tau_q[0]) / (q_values[1] - q_values[0])
        d_values[0] = q_values[0] * h_values[0] - tau_q[0]
        
        h_values[-1] = (tau_q[-1] - tau_q[-2]) / (q_values[-1] - q_values[-2])
        d_values[-1] = q_values[-1] * h_values[-1] - tau_q[-1]
        
        # Normalize D(h) to [0,1] range for better plotting
        d_max = np.max(d_values)
        if d_max > 0:
            d_values = d_values / d_max
        
        # Filter out points with invalid or extreme values
        valid_mask = np.isfinite(h_values) & np.isfinite(d_values)
        h_values = h_values[valid_mask]
        d_values = d_values[valid_mask]
        
        return h_values, d_values

    def plot_Dh_vs_h(self,signal_data, raster_layer_name,results, q_values=None, ax=None):
        """
        Generate and plot the D(h) vs h multifractal spectrum.
        
        Parameters:
        -----------
        signal_data : array-like
            Original signal data
        results : dict
            Results from wtmm_1d function
        q_values : array-like or None
            Values of q to use for spectrum
        ax : matplotlib axis or None
            Axis to plot on or None to create new figure
            
        Returns:
        --------
        matplotlib axis
        """
        # Create figure if not provided
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 6))
        
        # Compute spectrum
        h_values, d_values = self.simplified_multifractal_spectrum(results, q_values)
        
        # Plot D(h) vs h
        ax.plot(h_values, d_values, 'bo-', markersize=5)
        ax.set_xlabel('h', fontsize=12)
        ax.set_ylabel('D(h)', fontsize=12)
        ax.set_title('Multifractal Spectrum: {}'.format(raster_layer_name), fontsize=14)
        
        # Add grid
        ax.grid(True, linestyle='--', alpha=0.7)
        
        # Add horizontal line at y=1 for reference
        ax.axhline(y=1, color='r', linestyle='--', alpha=0.5)
        
        # Set axis limits
        h_min, h_max = np.min(h_values), np.max(h_values)
        h_range = h_max - h_min
        ax.set_xlim(h_min - 0.1 * h_range, h_max + 0.1 * h_range)
        ax.set_ylim(-0.1, 1.1)
        
        return ax

    # Example usage
    if __name__ == "__main__":
        # This should be imported from your existing code
        
        # Generate a test signal
        signal_data = generate_test_signal(length=1024, num_periods=20)
        #signal_data=example_data()
        #signal_data = financial_cartoon(Iterations=10, Multifractal=False, plot=False)[1]
        # Apply WTMM
        results = wtmm_1d(
            signal_data, 
            wavelet='mexh',
            num_scales=15,
            threshold_rel=0.05,  # Lower threshold to detect more maxima
            min_distance=3
        )
        
        # Plot the D(h) vs h spectrum
        fig, ax = plt.subplots(figsize=(8, 6))
        plot_Dh_vs_h(signal_data, results, ax=ax)
        plt.tight_layout()
        plt.show()
        visualize_wtmm_1d(signal_data, results, save_path=None)

        plt.show()



