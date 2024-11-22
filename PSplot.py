from PyQt5.QtWidgets import QVBoxLayout, QWidget
from matplotlib.backends.backend_qt5 import FigureCanvas
from scipy.fftpack import fft2, fftshift
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


class PowerSpectrumDock:
    def __init__(self, grid,gridname,dx, dy):
        """
        Initialize the PowerSpectrumDock class to display plots in a docked widget.

        Parameters:
            dock_widget (QDockWidget): The docked widget to display the plot.
            grid (numpy.ndarray): The 2D grid data.
            dx (float): Grid spacing in the x-direction.
            dy (float): Grid spacing in the y-direction.
        """
        self.grid = grid
        self.dx = dx
        self.dy = dy
        self.gridname=gridname


    def plot_grid_and_power_spectrum(self):
        """
        Plot the grid and its radially averaged power spectrum.
        """

        # Compute and plot radially averaged power spectrum
        self.grid = np.nan_to_num(self.grid, nan=0.0)

        #radial_frequencies, power_spectrum = self.radially_averaged_power_spectrum2(self.grid, self.dx, self.dy,apply_window=True, pad_factor=2)
        radial_frequencies, power_spectrum = self.radially_averaged_power_spectrum(self.grid, self.dx, self.dy)
        # Convert radial frequencies to wavelength (wavelength = 1 / frequency)
        wavelengths = 1 / np.array(radial_frequencies)
        valid_indices = ~np.isinf(wavelengths) & ~np.isnan(wavelengths)  # Exclude invalid wavelengths
        wavelengths = wavelengths[valid_indices]
        power_spectrum = power_spectrum[valid_indices]

        """        # Interpolate to finer resolution
                if len(wavelengths) > 1 and len(power_spectrum) > 1:
                    interp_func = interp1d(wavelengths, power_spectrum, kind="cubic", bounds_error=False, fill_value="extrapolate")
                    fine_wavelengths = np.logspace(np.log10(wavelengths.min()), np.log10(wavelengths.max()), num=500)
                    
                    # Ensure the interpolated range is strictly within bounds
                    fine_wavelengths = np.clip(fine_wavelengths, wavelengths.min(), wavelengths.max())
                    fine_power_spectrum = interp_func(fine_wavelengths)
                else:
        """     
        fine_wavelengths = wavelengths
        fine_power_spectrum = power_spectrum

        # Create a standalone figure
        fig, (ax_image, ax_spectrum) = plt.subplots(1, 2, figsize=(12, 6))

        # Plot the grid with the correct aspect ratio
        # Compute 5th and 95th percentiles
        vmin, vmax = np.percentile(self.grid, [5, 95])
        ax_image.imshow(self.grid, cmap="viridis", origin="lower", vmin=vmin,vmax=vmax,aspect=self.dy/self.dx)
        ax_image.set_title(self.gridname)
        ax_image.set_xlabel("X")
        ax_image.set_ylabel("Y")

        # Plot the power spectrum with log-wavelength axis and linear power
        if len(fine_wavelengths) > 0 and len(fine_power_spectrum) > 0:
            ax_spectrum.loglog(fine_wavelengths, fine_power_spectrum,  linestyle="-")
            ax_spectrum.set_xscale("log")  # Logarithmic scale for wavelength
            ax_spectrum.set_xlim(ax_spectrum.get_xlim()[::-1])
            ax_spectrum.set_title("Radially Averaged Power Spectrum")
            ax_spectrum.set_xlabel("Wavelength (log scale)")
            ax_spectrum.set_ylabel("Power Spectrum (log Scale)")
            ax_spectrum.grid(True)
        else:
            ax_spectrum.text(0.5, 0.5, "No data", transform=ax_spectrum.transAxes,
                            ha="center", va="center", fontsize=12)
        
        # Display the plot
        plt.tight_layout()
        plt.show()

    def radially_averaged_power_spectrum(self,grid, dx, dy):
        """
        Compute the radially averaged power spectrum of a 2D grid.

        Parameters:
            grid (numpy.ndarray): Input 2D grid data.
            dx (float): Grid spacing in the x-direction.
            dy (float): Grid spacing in the y-direction.

        Returns:
            tuple: (radial frequencies, power spectrum)
        """
        # Fourier transform and shift the zero frequency to the center
        fft_data = fft2(grid)
        fft_data_shifted = fftshift(fft_data)
        power_spectrum = np.abs(fft_data_shifted) ** 2

        # Grid dimensions and frequency components
        ny, nx = grid.shape
        kx = fftshift(np.fft.fftfreq(nx, dx))
        ky = fftshift(np.fft.fftfreq(ny, dy))
        kx, ky = np.meshgrid(kx, ky)
        radial_distance = np.sqrt(kx**2 + ky**2)

        # Radial bins
        max_radius = np.min([np.max(kx), np.max(ky)])
        radial_bins = np.linspace(0, max_radius, min(nx, ny) // 2)
        radial_bin_centers = 0.5 * (radial_bins[:-1] + radial_bins[1:])

        # Radial averaging
        radial_average = np.zeros_like(radial_bin_centers)
        for i in range(len(radial_bins) - 1):
            mask = (radial_distance >= radial_bins[i]) & (radial_distance < radial_bins[i + 1])
            radial_average[i] = np.mean(power_spectrum[mask]) if np.any(mask) else 0

        return radial_bin_centers, radial_average
    
    def radially_averaged_power_spectrum2(self,grid, dx, dy, apply_window=True, pad_factor=2):
        """
        Compute the radially averaged power spectrum with optional windowing and padding.

        Parameters:
            grid (numpy.ndarray): Input 2D grid data.
            dx (float): Grid spacing in the x-direction.
            dy (float): Grid spacing in the y-direction.
            apply_window (bool): Apply Hanning window to reduce spectral leakage.
            pad_factor (int): Factor by which to pad the grid size for extended frequency range.

        Returns:
            tuple: (radial frequencies, power spectrum)
        """
        # Handle NaNs
        grid = np.nan_to_num(grid, nan=0.0)

        # Apply windowing
        if apply_window:
            from scipy.signal import windows
            window = windows.hann(grid.shape[0])[:, None] * windows.hann(grid.shape[1])[None, :]
            grid = grid * window

        # Pad grid
        if pad_factor > 1:
            pad_width = grid.shape[0] // 2
            grid = np.pad(grid, pad_width=pad_width, mode='constant')

        # Fourier transform and power spectrum
        fft_data = fft2(grid)
        fft_data_shifted = fftshift(fft_data)
        power_spectrum = np.abs(fft_data_shifted) ** 2

        # Frequency components
        ny, nx = grid.shape
        kx = fftshift(np.fft.fftfreq(nx, d=dx))
        ky = fftshift(np.fft.fftfreq(ny, d=dy))
        kx, ky = np.meshgrid(kx, ky)
        radial_distance = np.sqrt(kx**2 + ky**2)

        # Radial binning
        max_radius = np.min([np.max(kx), np.max(ky)])
        radial_bins = np.linspace(0, max_radius, min(nx, ny) // 2)
        radial_bin_centers = 0.5 * (radial_bins[:-1] + radial_bins[1:])
        radial_average = np.zeros_like(radial_bin_centers)

        for i in range(len(radial_bins) - 1):
            mask = (radial_distance >= radial_bins[i]) & (radial_distance < radial_bins[i + 1])
            if np.any(mask):
                radial_average[i] = np.mean(power_spectrum[mask])

        return radial_bin_centers, radial_average
