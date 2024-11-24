from PyQt5.QtWidgets import QVBoxLayout, QWidget
from matplotlib.backends.backend_qt5 import FigureCanvas
from scipy.fftpack import fft2, fftshift
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


class PowerSpectrumDock:
    def __init__(self, grid,gridname,dx, dy,x,y):
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
        self.x = y
        self.y = x


    def plot_grid_and_power_spectrum(self):
        """
        Plot the grid and its radially averaged power spectrum.
        """

        # Compute and plot radially averaged power spectrum
        self.grid = np.nan_to_num(self.grid, nan=0.0)

        kx, ky, pds=self.power_density_spectra(self.x, self.y, self.grid, self.grid.shape)
        wavenumbers, power_spectrum = self.radial_average_spectrum(kx, ky, pds)
        valid_indices = ~np.isinf(wavenumbers) & ~np.isnan(wavenumbers)  # Exclude invalid wavelengths
        wavenumbers = wavenumbers[valid_indices]
        power_spectrum = power_spectrum[valid_indices]



        # Create a standalone figure
        fig, (ax_image, ax_spectrum) = plt.subplots(1, 2, figsize=(12, 6))

        # Plot the grid with the correct aspect ratio
        # Compute 5th and 95th percentiles
        vmin, vmax = np.percentile(self.grid, [5, 95])
        ax_image.imshow(np.flipud(self.grid), cmap="viridis", origin="lower", vmin=vmin,vmax=vmax,aspect=self.dy/self.dx)
        ax_image.set_title(self.gridname)
        ax_image.set_xlabel("X")
        ax_image.set_ylabel("Y")

        # Plot the power spectrum with log-wavelength axis and linear power
        if len(wavenumbers) > 0 and len(power_spectrum) > 0:
            ax_spectrum.plot(wavenumbers*self.dx, np.log(power_spectrum),  linestyle="-")
            ax_spectrum.set_title("Radially Averaged Power Spectrum")
            ax_spectrum.set_xlabel("Wavenumber (linear scale)")
            ax_spectrum.set_ylabel("ln(Power)")
            ax_spectrum.grid(True)
        else:
            ax_spectrum.text(0.5, 0.5, "No data", transform=ax_spectrum.transAxes,
                            ha="center", va="center", fontsize=12)
        
        # Display the plot
        plt.tight_layout()
        plt.show()


    def power_density_spectra(self,x, y, data, shape):
        """
        Calculates the Power Density Spectra of a 2D gridded potential field
        through the FFT:

        .. math::

            \Phi_{\Delta T}(k_x, k_y) = | F\left{\Delta T \right}(k_x, k_y) |^2

        .. note:: Requires gridded data.

        .. note:: x, y, z and height should be in meters.

        Parameters:

        * x, y : 1D-arrays
            The x and y coordinates of the grid points
        * data : 1D-array
            The potential field at the grid points
        * shape : tuple = (nx, ny)
            The shape of the grid

        Returns:

        * kx, ky : 2D-arrays
            The wavenumbers of each Power Density Spectra point
        * pds : 2D-array
            The Power Density Spectra of the data
        """
        kx, ky = self._fftfreqs(x, y, shape, shape)
        pds = abs(np.fft.fft2(np.reshape(data, shape)))**2
        return kx, ky, pds



    def radial_average_spectrum(self,kx, ky, pds, max_radius=None, ring_width=None):
        """
        Calculates the average of the Power Density Spectra points that falls
        inside concentric rings built around the origin of the wavenumber
        coordinate system with constant width.

        The width of the rings and the inner radius of the biggest ring can be
        changed by setting the optional parameters ring_width and max_radius,
        respectively.

        .. note:: To calculate the radially averaged power density spectra
                use the outputs of the function power_density_spectra as
                input of this one.

        Parameters:

        * kx, ky : 2D-arrays
            The x and y coordinates of the grid points
        * data : 1D-array
            The potential field at the grid points
        * shape : tuple = (nx, ny)
            The shape of the grid
        * max_radius : float (optional)
            Inner radius of the biggest ring.
            By default it's set as the minimum of kx.max() and ky.max().
            Making it smaller leaves points outside of the averaging,
            and making it bigger includes points nearer to the boundaries.
        * ring_width : float (optional)
            Width of the rings.
            By default it's set as the largest value of :math:`\Delta k_x` and
            :math:`\Delta k_y`, being them the equidistances of the kx and ky
            arrays.
            Making it bigger gives more populated averages, and
            making it smaller lowers the ammount of points per ring
            (use it carefully).

        Returns:

        * k_radial : 1D-array
            Wavenumbers of each Radially Averaged Power Spectrum point.
            Also, the inner radius of the rings.
        * pds_radial : 1D array
            Radially Averaged Power Spectrum
        """
        nx, ny = pds.shape
        if max_radius is None:
            max_radius = min(kx.max(), ky.max())
        if ring_width is None:
            ring_width = max(kx[1, 0], ky[0, 1])
        k = np.sqrt(kx**2 + ky**2)
        pds_radial = []
        k_radial = []
        radius_i = -1
        while True:
            radius_i += 1
            if radius_i*ring_width > max_radius:
                break
            else:
                if radius_i == 0:
                    inside = k <= 0.5*ring_width
                else:
                    inside = np.logical_and(k > (radius_i - 0.5)*ring_width,
                                            k <= (radius_i + 0.5)*ring_width)
                pds_radial.append(pds[inside].mean())
                k_radial.append(radius_i*ring_width)
        return np.array(k_radial), np.array(pds_radial)



    def _pad_data(self,data, shape):
        n = self._nextpow2(np.max(shape))
        nx, ny = shape
        padx = (n - nx)//2
        pady = (n - ny)//2
        padded = np.pad(data.reshape(shape), ((padx, padx), (pady, pady)),
                        mode='edge')
        return padded, padx, pady


    def _nextpow2(self,i):
        buf = np.ceil(np.log(i)/np.log(2))
        return int(2**buf)


    def _fftfreqs(self,x, y, shape, padshape):
        """
        Get two 2D-arrays with the wave numbers in the x and y directions.
        """
        nx, ny = shape
        dx = (x.max() - x.min())/(nx - 1)
        fx = 2*np.pi*np.fft.fftfreq(padshape[0], dx)
        dy = (y.max() - y.min())/(ny - 1)
        fy = 2*np.pi*np.fft.fftfreq(padshape[1], dy)
        return np.meshgrid(fy, fx)[::-1]