import numpy as np
from scipy.fftpack import fft2, ifft2,fftfreq
from scipy.ndimage import gaussian_filter, sobel
from numpy.polynomial.polynomial import polyvander2d, polyval2d
from scipy.spatial import cKDTree
from scipy.ndimage import uniform_filter



class GeophysicalProcessor:
    def __init__(self, dx, dy,buffer):
        """
        Initialize the processor with grid spacing.
        
        Parameters:
            dx (float): Grid spacing in the x-direction.
            dy (float): Grid spacing in the y-direction.
        """
        self.dx = dx
        self.dy = dy
        self.buffer=buffer

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
        filled_data = np.copy(data)
        filled_data[nan_mask] = np.nanmean(data)
        return filled_data, nan_mask

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
        if method == "mirror":
            return np.pad(data, pad_width=buffer_size, mode="reflect")
        elif method == "zero":
            return np.pad(data, pad_width=buffer_size, mode="constant", constant_values=0)
        else:
            raise ValueError("Invalid buffer method. Choose 'mirror' or 'zero'.")

    def remove_buffer(self, data, buffer_size):
        """
        Remove the buffer from the edges of the grid.
        """
        return data[buffer_size:-buffer_size, buffer_size:-buffer_size]

    # --- Trend Removal ---
    def remove_polynomial_trend(self, x, y, z, grid_x, grid_y, degree=2):
        """
        Remove a regional trend using polynomial fitting.
        """
        X = np.column_stack([x**i * y**j for i in range(degree + 1) for j in range(degree + 1 - i)])
        coeffs, _, _, _ = np.linalg.lstsq(X, z, rcond=None)

        # Evaluate polynomial on grid
        trend = polyval2d(grid_x, grid_y, coeffs.reshape(degree + 1, -1))
        residual = z - polyval2d(x, y, coeffs.reshape(degree + 1, -1))

        return residual, trend

    def remove_regional_trend_fourier(self, data, cutoff_wavelength, buffer_size=10, buffer_method="mirror"):
        """
        Remove the regional trend using a low-pass filter in the Fourier domain.
        """
        return self.low_pass_filter(data, cutoff_wavelength, buffer_size, buffer_method)

    #vertical integration 
    def vertical_integration(self, data, max_wavenumber=None, buffer_size=10, buffer_method="mirror"):
        """
        Perform vertical integration of a field in the frequency domain.

        Parameters:
            data (numpy.ndarray): 2D array of the field data.
            max_wavenumber (float): Optional maximum wavenumber for filtering (to stabilize high-frequency noise).
            buffer_size (int): Buffer size for edge handling.
            buffer_method (str): Buffering method ('mirror' or 'zero').

        Returns:
            numpy.ndarray: Vertically integrated data.
        """
        def filter_function(kx, ky):
            k = np.sqrt(kx**2 + ky**2) + 1e-10  # Avoid division by zero
            vertical_integration_filter = 1 / k

            # Apply max_wavenumber filtering if specified
            if max_wavenumber:
                vertical_integration_filter[k > max_wavenumber] = 0

            return vertical_integration_filter

        # Apply the Fourier filter
        return self._apply_fourier_filter(data, filter_function, buffer_size, buffer_method)

    # --- Continuation ---
    def upward_continuation(self, data, height, buffer_size=10, buffer_method="mirror"):
        """
        Perform upward continuation to attenuate high-frequency signals.
        """
        def filter_function(kx, ky):
            k = np.sqrt(kx**2 + ky**2)
            return np.exp(-k * height)

        return self._apply_fourier_filter(data, filter_function, buffer_size, buffer_method)

    def downward_continuation(self, data, height, buffer_size=10, buffer_method="mirror"):
        """
        Perform downward continuation to enhance high-frequency signals.
        """
        def filter_function(kx, ky):
            k = np.sqrt(kx**2 + ky**2)
            return np.exp(k * height)

        return self._apply_fourier_filter(data, filter_function, buffer_size, buffer_method)

    # --- Reduction Methods ---
    def reduction_to_pole(self, data, inclination, declination, buffer_size=10, buffer_method="mirror"):
        """
        Perform Reduction to Pole (RTP).
        """
        def filter_function(kx, ky):
            k = np.sqrt(kx**2 + ky**2)
            incl = np.radians(inclination)
            decl = np.radians(declination)

            T = np.sin(incl) - 1j * np.cos(incl) * np.cos(decl) * kx / (k + 1e-10)
            return T / np.sqrt(T.real**2 + T.imag**2)

        return self._apply_fourier_filter(data, filter_function, buffer_size, buffer_method)

    def reduction_to_pole(self, data, inclination, declination, buffer_size=10, buffer_method="mirror"):
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

        # Unit vectors for the geomagnetic field
        fe = np.sin(dec) * np.cos(inc)  # Easting component
        fn = np.cos(dec) * np.cos(inc)  # Northing component
        fz = np.sin(inc)                # Downward component

        # Compute Theta_m and Theta_f
        theta_f = fz + 1j * (fe * kx + fn * ky) / k

        # RTP filter
        rtp_filter = 1.0 / theta_f
        rtp_filter[np.abs(rtp_filter) > 1e6] = 0  # Stabilize extreme values

        # Apply RTP filter in the Fourier domain
        data_fft = fft2(buffered_data)
        data_rtp_fft = data_fft * rtp_filter.T

        data_rtp = np.real(ifft2(data_rtp_fft))

        # Remove buffer and restore NaNs
        filtered_data = self.remove_buffer(data_rtp, buffer_size)
        return self.restore_nan(filtered_data, nan_mask)

    def reduction_to_equator(self, data, inclination, declination, buffer_size=10, buffer_method="mirror"):
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

        # Unit vectors for the geomagnetic field
        fe = np.sin(dec) * np.cos(inc)  # Easting component
        fn = np.cos(dec) * np.cos(inc)  # Northing component
        fz = np.sin(inc)                # Downward component

        # Compute Theta_m and Theta_f
        theta_f = fz + 1j * (fe * kx + fn * ky) / k

        # RTP filter
        rtp_filter = theta_f
        rtp_filter[np.abs(rtp_filter) > 1e6] = 0  # Stabilize extreme values

        # Apply RTP filter in the Fourier domain
        data_fft = fft2(buffered_data)
        data_rtp_fft = data_fft * rtp_filter.T

        data_rtp = np.real(ifft2(data_rtp_fft))

        # Remove buffer and restore NaNs
        filtered_data = self.remove_buffer(data_rtp, buffer_size)
        return self.restore_nan(filtered_data, nan_mask)
    

    # --- Derivatives ---
    def compute_derivative(self, data, direction, order=1, buffer_size=10, buffer_method="mirror"):
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

        return self._apply_fourier_filter(data, filter_function, buffer_size, buffer_method)

    def tilt_derivative(self, data, buffer_size=10, buffer_method="mirror"):
        """
        Compute the tilt derivative of the data.
        """
        dfdx = self.compute_derivative(data, "x", buffer_size=buffer_size, buffer_method=buffer_method)
        dfdy = self.compute_derivative(data, "y", buffer_size=buffer_size, buffer_method=buffer_method)
        horizontal_gradient = np.sqrt(dfdx**2 + dfdy**2)
        dfdz = self.compute_derivative(data, "z", buffer_size=buffer_size, buffer_method=buffer_method)
        return np.arctan2(dfdz, horizontal_gradient)

    def analytic_signal(self, data, buffer_size=10, buffer_method="mirror"):
        """
        Compute the analytic signal amplitude.
        """
        dfdx = self.compute_derivative(data, "x", buffer_size=buffer_size, buffer_method=buffer_method)
        dfdy = self.compute_derivative(data, "y", buffer_size=buffer_size, buffer_method=buffer_method)
        dfdz = self.compute_derivative(data, "z", buffer_size=buffer_size, buffer_method=buffer_method)
        return np.sqrt(dfdx**2 + dfdy**2 + dfdz**2)



    def automatic_gain_control(self,grid, window_size):
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
        padded_abs_grid = np.pad(abs_grid, pad_width=int(3 * window_size), mode='reflect')

        # Apply Gaussian smoothing to the padded grid
        smoothed_padded_grid = gaussian_filter(padded_abs_grid, sigma=window_size)

        # Remove padding
        smoothed_grid = smoothed_padded_grid[
            int(3 * window_size):-int(3 * window_size),
            int(3 * window_size):-int(3 * window_size)
        ]

        # Prevent division by zero by adding a small constant
        smoothed_grid[smoothed_grid == 0] = np.finfo(float).eps

        # Apply AGC: normalize the original grid by the smoothed grid
        agc_grid = grid / smoothed_grid

        return  self.restore_nan(agc_grid, nan_mask)

    def minimum_curvature_gridding(self,x, y, z, grid_x, grid_y, max_iterations=1000, tolerance=1e-6):
        """
        Perform minimum curvature gridding to interpolate scattered data onto a regular grid.

        Parameters:
            x, y (1D array): Coordinates of the input data points.
            z (1D array): Values at the input data points.
            grid_x, grid_y (2D arrays): Regular grid coordinates for interpolation.
            max_iterations (int): Maximum number of iterations for convergence.
            tolerance (float): Convergence criterion based on change in grid values.

        Returns:
            2D array: Interpolated grid values.
        """
        # Initialize grid
        grid = np.full(grid_x.shape, np.nan)

        # Assign known values to the grid
        for xi, yi, zi in zip(x, y, z):
            # Find the nearest grid point
            ix = np.argmin(np.abs(grid_x[0, :] - xi))
            iy = np.argmin(np.abs(grid_y[:, 0] - yi))
            grid[iy, ix] = zi

        # Fill initial NaN values with the average of known values
        nan_mask = np.isnan(grid)
        grid[nan_mask] = np.nanmean(grid)

        # Iteratively solve for minimum curvature
        for iteration in range(max_iterations):
            grid_old = grid.copy()

            # Update each grid cell
            for i in range(1, grid.shape[0] - 1):
                for j in range(1, grid.shape[1] - 1):
                    if nan_mask[i, j]:  # Only update NaN (interpolated) cells
                        grid[i, j] = (
                            0.25 * (grid[i + 1, j] + grid[i - 1, j] + grid[i, j + 1] + grid[i, j - 1])
                            - 0.125 * (grid[i + 1, j + 1] + grid[i - 1, j - 1] + grid[i + 1, j - 1] + grid[i - 1, j + 1])
                        )

            # Check for convergence
            max_change = np.nanmax(np.abs(grid - grid_old))
            if max_change < tolerance:
                print(f"Converged after {iteration + 1} iterations.")
                break
        else:
            print("Maximum iterations reached without full convergence.")

        return np.flipud(grid)

    def minimum_curvature_griddingx(self,x, y, z, grid_x, grid_y, max_iterations=1000, tolerance=1e-4):
        """
        Perform minimum curvature gridding on scattered data.

        Parameters:
            points (iterable): Iterable of (x, y, z) tuples.
            grid_x (numpy.ndarray): 2D array of x-coordinates for the output grid.
            grid_y (numpy.ndarray): 2D array of y-coordinates for the output grid.
            max_iterations (int): Maximum number of iterations for the solver.
            tolerance (float): Convergence tolerance.

        Returns:
            numpy.ndarray: Gridded values.
        """
        points=zip(x, y, z)
        # Initialize the grid with a distance-weighted average
        grid = np.zeros_like(grid_x, dtype=float)
        for xi, yi, zi in zip(x, y, z):
            distances = np.sqrt((grid_x - xi)**2 + (grid_y - yi)**2)
            distances[distances == 0] = 1e-10  # Prevent division by zero
            grid += zi / distances
        grid /= len(x)

        # Mask for grid cells that do not have observed values
        mask = np.isnan(grid)

        # Iterative solver
        for iteration in range(max_iterations):
            old_grid = grid.copy()

            # Update the grid using a 5-point Laplacian smoothing operator
            for i in range(1, grid.shape[0] - 1):
                for j in range(1, grid.shape[1] - 1):
                    if mask[i, j]:
                        grid[i, j] = 0.25 * (old_grid[i+1, j] + old_grid[i-1, j] +
                                            old_grid[i, j+1] + old_grid[i, j-1])

            # Check for convergence
            max_diff = np.nanmax(np.abs(grid - old_grid))
            if max_diff < tolerance:
                print(f"Converged after {iteration} iterations.")
                break
        else:
            print("Did not converge within the maximum number of iterations.")

        return np.flipud(grid)


    def low_pass_filter(self, data, cutoff_wavelength, buffer_size=10, buffer_method="mirror"):
        """
        Apply a low-pass filter to suppress high-frequency noise.
        """
        def filter_function(kx, ky):
            k = np.sqrt(kx**2 + ky**2)
            cutoff_k = 2 * np.pi / (cutoff_wavelength+1e-10)
            return k < cutoff_k  # Low-pass filter mask

        return self._apply_fourier_filter(data, filter_function, buffer_size, buffer_method)

    def high_pass_filter(self, data, cutoff_wavelength, buffer_size=10, buffer_method="mirror"):
        """
        Apply a high-pass filter to remove regional trends.
        """
        def filter_function(kx, ky):
            k = np.sqrt(kx**2 + ky**2)
            cutoff_k = 2 * np.pi / (cutoff_wavelength+1e-10)
            return k > cutoff_k  # High-pass filter mask

        return self._apply_fourier_filter(data, filter_function, buffer_size, buffer_method)

    def band_pass_filter(self, data, low_cut, high_cut, buffer_size=10, buffer_method="mirror"):
        """
        Apply a band-pass filter to isolate anomalies within a wavelength range.
        """
        def filter_function(kx, ky):
            k = np.sqrt(kx**2 + ky**2)
            high_cut_k = 2 * np.pi / low_cut
            low_cut_k = 2 * np.pi / high_cut

            # Create a band-pass mask
            filter_mask = (k >= low_cut_k) & (k <= high_cut_k)


            return (filter_mask)  # Band-pass filter mask

        return self._apply_fourier_filter(data, filter_function, buffer_size, buffer_method)
    
    def directional_band_reject_filter(self, data, low_cut, high_cut, direction_angle, buffer_size=10, buffer_method="mirror"):
        """
        Apply a directional band-reject filter to isolate anomalies in a specific frequency band and direction.

        Parameters:
            data (numpy.ndarray): Input 2D grid.
            low_cut (float): Minimum wavelength to reject.
            high_cut (float): Maximum wavelength to reject.
            direction_angle (float): Direction of the filter in degrees (clockwise from North).
            buffer_size (int): Buffer size for edge handling.
            buffer_method (str): Buffering method ('mirror' or 'zero').

        Returns:
            numpy.ndarray: Filtered data with specified directional band removed.
        """
        if low_cut >= high_cut:
            raise ValueError("low_cut must be smaller than high_cut.")

        # Convert direction angle to radians
        theta = np.radians(direction_angle)

        def filter_function(kx, ky):
            k = np.sqrt(kx**2 + ky**2) + 1e-10  # Avoid division by zero
            high_cut_k = 2 * np.pi / low_cut
            low_cut_k = 2 * np.pi / high_cut

            # Directional wavenumber component
            directional_k = kx * np.cos(theta) + ky * np.sin(theta)

            # Create the directional band-reject mask
            filter_mask = ((k >= low_cut_k) & (k <= high_cut_k) & (np.abs(directional_k) > 0))
            # Debugging: Ensure the filter mask is valid


            return filter_mask
        
        return self._apply_fourier_filter(data, filter_function, buffer_size, buffer_method)

    def apply_window(data, window="hann"):
        """
        Apply a tapering window function to reduce edge effects.

        Parameters:
            data (2D array): Input grid.
            window (str): Type of window function ('hann', 'hamming', etc.).

        Returns:
            2D array: Windowed grid.
        """
        ny, nx = data.shape

        if window == "hann":
            wx = np.hanning(nx)
            wy = np.hanning(ny)
        elif window == "hamming":
            wx = np.hamming(nx)
            wy = np.hamming(ny)
        else:
            raise ValueError("Unsupported window type.")

        window_2d = np.outer(wy, wx)
        return data * window_2d
    
    # --- Internal Fourier Filter Application ---
    def _apply_fourier_filter(self, data, filter_function, buffer_size=10, buffer_method="mirror"):
        """
        Apply a Fourier-domain filter with buffering and NaN handling.
        """
        # Handle NaN values
        filled_data, nan_mask = self.fill_nan(data)

        # Add buffer
        buffered_data = self.add_buffer(filled_data, buffer_size, buffer_method)

        #buffered_data=self.apply_window(buffered_data)
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


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    # Create synthetic data
    nx, ny = 100, 100
    dx, dy = 1, 1
    x = np.linspace(0, 99, nx)
    y = np.linspace(0, 99, ny)
    X, Y = np.meshgrid(x, y)
    data = 0.01 * X**2 + 0.02 * Y + np.sin(2 * np.pi * X / 10) + np.cos(2 * np.pi * Y / 15)

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
    residual, trend = processor.remove_polynomial_trend(x, y, z, grid_x, grid_y, degree=2)

    # Plot results
    plt.imshow(trend, cmap="viridis", origin="lower")
    plt.colorbar()
    plt.title("Polynomial Trend")
    plt.show()

    # Remove regional trend using Fourier (low-pass filter)
    regional_removed = processor.remove_regional_trend_fourier(data, cutoff_wavelength=30)

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
    plt.imshow(tilt_result, cmap="seismic", origin="lower", vmin=-np.pi/2, vmax=np.pi/2)
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

