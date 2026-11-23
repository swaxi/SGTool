import numpy as np

# longest side, in cells, of the grid drawn beside the spectrum: a bigger grid is
# block-averaged for the picture (the spectrum always uses the full grid)
DISPLAY_SIDE = 1500


class PowerSpectrumDock:
    def __init__(self, grid, gridname, dx, dy, x, y):
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
        self.gridname = gridname
        self.x = y
        self.y = x

    def plot_grid_and_power_spectrum(self, check=None):
        """
        Plot the grid and its radially averaged power spectrum.
        """
        self.plot(self.compute(check))

    def compute(self, check=None):
        """The radially averaged power spectrum and a picture-sized copy of the
        grid: numpy only, so it can run in a background task. check(percent)
        is called between steps (it raises to cancel)."""
        check = check or (lambda percent=None: None)
        self.grid = np.nan_to_num(self.grid, nan=0.0)
        check(5)
        kx, ky, pds, weights = self.power_density_spectra(
            self.x, self.y, self.grid, self.grid.shape, half=True
        )
        check(60)
        wavenumbers, power_spectrum = self.radial_average_spectrum(
            kx, ky, pds, weights=weights
        )
        del pds
        check(85)
        valid_indices = ~np.isinf(wavenumbers) & ~np.isnan(
            wavenumbers
        )  # Exclude invalid wavelengths
        wavenumbers = wavenumbers[valid_indices]
        power_spectrum = power_spectrum[valid_indices]
        display = self._display_grid(self.grid)
        vmin, vmax = np.percentile(display, [5, 95])
        check(100)
        return display, wavenumbers, power_spectrum, (vmin, vmax)

    def plot(self, computed):
        """Draw what compute() returned (on the main thread)."""
        import matplotlib.pyplot as plt

        display, wavenumbers, power_spectrum, (vmin, vmax) = computed

        # Create a standalone figure
        fig, (ax_image, ax_spectrum) = plt.subplots(1, 2, figsize=(12, 6))

        # Plot the grid with the correct aspect ratio (5th to 95th percentile)
        ax_image.imshow(
            np.flipud(display),
            cmap="viridis",
            origin="lower",
            vmin=vmin,
            vmax=vmax,
            aspect=self.dy / self.dx,
        )
        ax_image.set_title(self.gridname)
        ax_image.set_xlabel("X")
        ax_image.set_ylabel("Y")

        # Plot the power spectrum with log-wavelength axis and linear power
        if len(wavenumbers) > 0 and len(power_spectrum) > 0:
            ax_spectrum.plot(
                wavenumbers * self.dx, np.log(power_spectrum), linestyle="-"
            )
            ax_spectrum.set_title("Radially Averaged Power Spectrum")
            ax_spectrum.set_xlabel("Wavenumber (linear scale)")
            ax_spectrum.set_ylabel("ln(Power)")
            ax_spectrum.grid(True)
        else:
            ax_spectrum.text(
                0.5,
                0.5,
                "No data",
                transform=ax_spectrum.transAxes,
                ha="center",
                va="center",
                fontsize=12,
            )

        # Display the plot
        plt.tight_layout()
        plt.show()

    @staticmethod
    def _display_grid(grid, max_side=DISPLAY_SIDE):
        """The grid block-averaged to at most max_side cells on its longest side
        (drawing millions of cells is slow and shows no more on screen)."""
        f = int(np.ceil(max(grid.shape) / float(max_side)))
        if f <= 1:
            return grid
        rows, cols = (grid.shape[0] // f) * f, (grid.shape[1] // f) * f
        return grid[:rows, :cols].reshape(rows // f, f, cols // f, f).mean(axis=(1, 3))

    def power_density_spectra(self, x, y, data, shape, half=False):
        """
        Calculates the Power Density Spectra of a 2D gridded potential field
        through the FFT:

        .. math::

            Phi_{Delta T}(k_x, k_y) = | F div left{Delta T right}(k_x, k_y) |^2

        .. note:: Requires gridded data.

        .. note:: x, y, z and height should be in meters.

        Parameters:

        * x, y : 1D-arrays
            The x and y coordinates of the grid points
        * data : 1D-array
            The potential field at the grid points
        * shape : tuple = (nx, ny)
            The shape of the grid
        * half : bool
            Return only the half spectrum of the real transform, with weights
            (1 or 2) counting each bin as often as it appears in the full
            spectrum. |F(-k)| = |F(k)| for a real grid and -k is in the same
            ring, so the radial average is the same, at half the time and
            memory. The wavenumbers are then broadcastable (rows x 1, 1 x cols)
            arrays rather than full grids.

        Returns:

        * kx, ky : 2D-arrays
            The wavenumbers of each Power Density Spectra point
        * pds : 2D-array
            The Power Density Spectra of the data
        * weights (only if half) : 1D-array over the columns
        """
        if not half:
            kx, ky = self._fftfreqs(x, y, shape, shape)
            from scipy.fft import fft2

            # the full spectrum; all cores
            pds = abs(fft2(np.reshape(data, shape), workers=-1)) ** 2
            return kx, ky, pds

        from scipy.fft import rfft2

        fx, fy = self._fftfreqs_1d(x, y, shape, shape)
        n1 = shape[1]
        spec = rfft2(np.reshape(data, shape), workers=-1)
        h = spec.shape[1]  # n1 // 2 + 1
        pds = spec.real ** 2 + spec.imag ** 2
        del spec
        weights = np.full(h, 2.0)
        weights[0] = 1.0  # column 0 is its own mirror
        if n1 % 2 == 0:
            weights[-1] = 1.0  # so is the Nyquist column of an even width
        return fx[:, None], fy[None, :h], pds, weights

    def radial_average_spectrum(self, kx, ky, pds, max_radius=None, ring_width=None,
                                weights=None):
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

        * kx, ky : 2D-arrays (or arrays that broadcast to the spectrum's shape)
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
            By default it's set as the largest value of :math:`Delta k_x` and
            :math:`Delta k_y`, being them the equidistances of the kx and ky
            arrays.
            Making it bigger gives more populated averages, and
            making it smaller lowers the ammount of points per ring
            (use it carefully).
        * weights : 1D-array over the columns (optional)
            How often each column's bins count (from power_density_spectra
            with half=True).

        Returns:

        * k_radial : 1D-array
            Wavenumbers of each Radially Averaged Power Spectrum point.
            Also, the inner radius of the rings.
        * pds_radial : 1D array
            Radially Averaged Power Spectrum
        """
        if max_radius is None:
            max_radius = min(kx.max(), ky.max())
        if ring_width is None:
            ring_width = max(kx[1, 0], ky[0, 1])
        w = ring_width
        # rings 0, 1, ... while radius_i * ring_width <= max_radius
        n = int(max_radius // w) + 1
        while n * w <= max_radius:
            n += 1
        while n > 1 and (n - 1) * w > max_radius:
            n -= 1

        # every bin's ring in one pass (not a mask over the whole spectrum per
        # ring): ring 0 is k <= w/2, ring i is (i - 1/2) w < k <= (i + 1/2) w
        k = np.sqrt(kx ** 2 + ky ** 2)
        ring = np.ceil(k / w - 0.5).astype(np.int64)
        np.maximum(ring, 0, out=ring)
        # put bins rounded across a boundary back on the side the comparisons give
        ring[(ring > 0) & (k <= (ring - 0.5) * w)] -= 1
        ring[k > (ring + 0.5) * w] += 1
        del k
        keep = ring < n
        if weights is None:
            wts = np.ones_like(pds)
        else:
            wts = np.broadcast_to(weights, pds.shape)
        sums = np.bincount(ring[keep], weights=(pds * wts)[keep], minlength=n)[:n]
        counts = np.bincount(ring[keep], weights=wts[keep], minlength=n)[:n]
        with np.errstate(invalid="ignore", divide="ignore"):
            pds_radial = sums / counts  # an empty ring gives NaN, as its mean did
        k_radial = np.array([i * w for i in range(n)])
        return k_radial, pds_radial

    def _pad_data(self, data, shape):
        n = self._nextpow2(np.max(shape))
        nx, ny = shape
        padx = (n - nx) // 2
        pady = (n - ny) // 2
        padded = np.pad(data.reshape(shape), ((padx, padx), (pady, pady)), mode="edge")
        return padded, padx, pady

    def _nextpow2(self, i):
        buf = np.ceil(np.log(i) / np.log(2))
        return int(2**buf)

    def _fftfreqs_1d(self, x, y, shape, padshape):
        """The wave numbers along axis 0 and axis 1 (as _fftfreqs, 1-D)."""
        nx, ny = shape
        dx = (x.max() - x.min()) / (nx - 1)
        fx = 2 * np.pi * np.fft.fftfreq(padshape[0], dx)
        dy = (y.max() - y.min()) / (ny - 1)
        fy = 2 * np.pi * np.fft.fftfreq(padshape[1], dy)
        return fx, fy

    def _fftfreqs(self, x, y, shape, padshape):
        """
        Get two 2D-arrays with the wave numbers in the x and y directions.
        """
        fx, fy = self._fftfreqs_1d(x, y, shape, padshape)
        return np.meshgrid(fy, fx)[::-1]
