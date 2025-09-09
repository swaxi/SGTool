"""
Euler deconvolution

A Python program to perform Euler deconvolution on gridded data.

This code is released from the paper:
Reliable Euler deconvolution estimates throughout the
vertical derivatives of the total-field anomaly

The program is under the conditions terms in the file README.txt

authors: Felipe F. Melo and Valeria C.F. Barbosa, 2019
email: felipe146@hotmail.com, valcris@on.br
"""

import numpy as np


"""
Optimized Euler deconvolution with vectorized operations

Key optimizations:
1. Vectorized sliding window operations using stride_tricks
2. Batch matrix operations instead of per-window loops
3. Numba JIT compilation for critical sections
4. Memory-efficient operations
"""

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

"""try:
    from numba import jit, prange
    import numba

    if numba.__version__ < "0.61.0":
        HAS_NUMBA = False
        # print("Numba >= 0.61.0 not available - using pure NumPy (will be slower)")

    else:
        HAS_NUMBA = Trues
except ImportError:
    HAS_NUMBA = False
    print("Numba not available - using pure NumPy (will be slower)")"""
HAS_NUMBA = False

def fft_pad_data(data, mode="edge"):
    """Optimized padding - unchanged from original"""
    n_points = int(2 ** (np.ceil(np.log(np.max(data.shape)) / np.log(2))))
    nx, ny = data.shape
    padx = (n_points - nx) // 2
    pady = (n_points - ny) // 2
    padded_data = np.pad(data, ((padx, padx), (pady, pady)), mode)
    mask = np.zeros_like(padded_data, dtype=bool)
    mask[padx : padx + nx, pady : pady + ny] = True
    fpdat = np.fft.fft2(padded_data)
    return fpdat, mask


def ifft_unpad_data(data_p, mask, shape_dat):
    """Optimized unpadding - unchanged from original"""
    ifft_data = np.real(np.fft.ifft2(data_p))
    data = ifft_data[mask]
    return np.reshape(data, shape_dat)


def fft_wavenumbers(x, y, shape, padshape):
    """Optimized wavenumber computation"""
    nx, ny = shape
    dx = (x.max() - x.min()) / (nx - 1)
    dy = (y.max() - y.min()) / (ny - 1)
    u = 2 * np.pi * np.fft.fftfreq(padshape[0], dx)
    v = 2 * np.pi * np.fft.fftfreq(padshape[1], dy)
    return np.meshgrid(v, u)[::-1]


def deriv(data, shape, area):
    """Optimized derivative computation"""
    anom_FFT, mask = fft_pad_data(data)
    nx, ny = shape
    xa, xb, ya, yb = area
    xs = np.linspace(xa, xb, nx)
    ys = np.linspace(ya, yb, ny)
    # Match original coordinate ordering - Y should be flipped to match data indexing
    Y, X = np.meshgrid(ys, xs, indexing="ij")

    u, v = fft_wavenumbers(X, Y, data.shape, anom_FFT.shape)

    # Vectorized derivative computation
    derivx_ft = anom_FFT * (u * 1j)
    derivy_ft = anom_FFT * (v * 1j)
    derivz_ft = anom_FFT * np.sqrt(u**2 + v**2)

    derivx = ifft_unpad_data(derivx_ft, mask, data.shape)
    derivy = ifft_unpad_data(derivy_ft, mask, data.shape)
    derivz = ifft_unpad_data(derivz_ft, mask, data.shape)

    return derivx, derivy, derivz


def create_sliding_windows(data, window_size):
    """
    Create sliding windows using stride_tricks for vectorized operations
    Much faster than nested loops
    """
    return sliding_window_view(data, (window_size, window_size))


if HAS_NUMBA:

    @jit(nopython=True, parallel=True)
    def solve_euler_systems_numba(
        windows_data,
        windows_dx,
        windows_dy,
        windows_dz,
        windows_x,
        windows_y,
        windows_z,
        SI,
        window_size,
    ):
        """
        Numba-compiled function to solve Euler systems in parallel
        """
        n_windows_y, n_windows_x = windows_data.shape[:2]

        # Pre-allocate output arrays
        estx = np.zeros((n_windows_y, n_windows_x))
        esty = np.zeros((n_windows_y, n_windows_x))
        estz = np.zeros((n_windows_y, n_windows_x))
        estb = np.zeros((n_windows_y, n_windows_x))
        stdzmat = np.zeros((n_windows_y, n_windows_x))

        # Process windows in parallel
        for i in prange(n_windows_y):
            for j in prange(n_windows_x):
                # Extract window data
                win_data = windows_data[i, j].ravel()
                win_dx = windows_dx[i, j].ravel()
                win_dy = windows_dy[i, j].ravel()
                win_dz = windows_dz[i, j].ravel()
                win_x = windows_x[i, j].ravel()
                win_y = windows_y[i, j].ravel()
                win_z = windows_z[i, j].ravel()

                # Build system matrix A
                n_points = window_size * window_size
                A = np.zeros((n_points, 4))
                A[:, 0] = win_dx
                A[:, 1] = win_dy
                A[:, 2] = win_dz
                A[:, 3] = SI

                # Build right-hand side vector
                b = win_dx * win_x + win_dy * win_y + win_dz * win_z + SI * win_data

                # Solve least squares system: A^T A p = A^T b
                ATA = np.dot(A.T, A)
                ATb = np.dot(A.T, b)

                # Simple 4x4 matrix inversion (faster than general case)
                try:
                    p = np.linalg.solve(ATA, ATb)

                    # Calculate standard deviation
                    stdz = np.sqrt(
                        np.sum((A[:, 2] - np.mean(A[:, 2])) ** 2) / (len(A[:, 2]) - 1.0)
                    )

                    estx[i, j] = p[0]
                    esty[i, j] = p[1]
                    estz[i, j] = p[2]
                    estb[i, j] = p[3]
                    stdzmat[i, j] = stdz

                except:
                    # Handle singular matrices
                    estx[i, j] = 0.0
                    esty[i, j] = 0.0
                    estz[i, j] = 0.0
                    estb[i, j] = 0.0
                    stdzmat[i, j] = 1e10

        return estx, esty, estz, estb, stdzmat


def solve_euler_systems_numpy(
    windows_data,
    windows_dx,
    windows_dy,
    windows_dz,
    windows_x,
    windows_y,
    windows_z,
    SI,
    window_size,
):
    """
    Pure NumPy version for when Numba is not available
    Still much faster than the original nested loop approach
    """
    n_windows_y, n_windows_x = windows_data.shape[:2]
    n_points = window_size * window_size

    # Reshape all windows to (n_windows, n_points)
    windows_data_flat = windows_data.reshape(-1, n_points)
    windows_dx_flat = windows_dx.reshape(-1, n_points)
    windows_dy_flat = windows_dy.reshape(-1, n_points)
    windows_dz_flat = windows_dz.reshape(-1, n_points)
    windows_x_flat = windows_x.reshape(-1, n_points)
    windows_y_flat = windows_y.reshape(-1, n_points)
    windows_z_flat = windows_z.reshape(-1, n_points)

    n_windows = windows_data_flat.shape[0]

    # Build system matrices for all windows at once
    A = np.zeros((n_windows, n_points, 4))
    A[:, :, 0] = windows_dx_flat
    A[:, :, 1] = windows_dy_flat
    A[:, :, 2] = windows_dz_flat
    A[:, :, 3] = SI

    # Build right-hand side vectors
    b = (
        windows_dx_flat * windows_x_flat
        + windows_dy_flat * windows_y_flat
        + windows_dz_flat * windows_z_flat
        + SI * windows_data_flat
    )

    # Solve all systems using batch operations
    print("Solving Euler systems using NumPy batch operations...")
    ATA = np.matmul(A.transpose(0, 2, 1), A)  # Shape: (n_windows, 4, 4)
    ATb = np.matmul(
        A.transpose(0, 2, 1), b[..., np.newaxis]
    )  # Shape: (n_windows, 4, 1)

    try:
        print("Computing solutions...")
        # Batch solve
        solutions = np.linalg.solve(ATA, ATb).squeeze()  # Shape: (n_windows, 4)

        # Calculate standard deviations
        dz_means = np.mean(windows_dz_flat, axis=1, keepdims=True)
        stdzmat_flat = np.sqrt(
            np.sum((windows_dz_flat - dz_means) ** 2, axis=1) / (n_points - 1)
        )

    except np.linalg.LinAlgError:
        # Handle potential singular matrices
        solutions = np.zeros((n_windows, 4))
        stdzmat_flat = np.full(n_windows, 1e10)

    # Reshape results back to 2D
    estx = solutions[:, 0].reshape(n_windows_y, n_windows_x)
    esty = solutions[:, 1].reshape(n_windows_y, n_windows_x)
    estz = solutions[:, 2].reshape(n_windows_y, n_windows_x)
    estb = solutions[:, 3].reshape(n_windows_y, n_windows_x)
    stdzmat = stdzmat_flat.reshape(n_windows_y, n_windows_x)

    return estx, esty, estz, estb, stdzmat


def euler_deconv_optimized(data, xi, yi, zi, shape, area, SI, windowSize, filt):
    """
    Optimized Euler deconvolution using vectorized operations

    Expected speedup: 10-100x depending on data size and hardware
    """
    # Compute derivatives (this part is already reasonably optimized)
    print("Computing derivatives...")
    print("Data shape:", data.shape)
    print("Coordinate shapes:", xi.shape, yi.shape, zi.shape)
    print("Area:", area)

    dx, dy, dz = deriv(data, shape, area)
    print("Derivatives computed")
    # Ensure coordinate arrays match data orientation
    # If yi values are inverted, we need to handle this properly
    if xi.shape != data.shape or yi.shape != data.shape:
        raise ValueError("Coordinate arrays must have same shape as data")

    # Create sliding windows for all arrays at once
    print("Creating sliding windows...")
    windows_data = create_sliding_windows(data, windowSize)
    windows_dx = create_sliding_windows(dx, windowSize)
    windows_dy = create_sliding_windows(dy, windowSize)
    windows_dz = create_sliding_windows(dz, windowSize)
    windows_x = create_sliding_windows(xi, windowSize)
    windows_y = create_sliding_windows(yi, windowSize)
    windows_z = create_sliding_windows(zi, windowSize)
    print("Sliding windows created")
    # Solve all Euler systems efficiently
    if HAS_NUMBA:
        print("Using Numba for optimized Euler deconvolution")
        estx, esty, estz, estb, stdzmat = solve_euler_systems_numba(
            windows_data,
            windows_dx,
            windows_dy,
            windows_dz,
            windows_x,
            windows_y,
            windows_z,
            SI,
            windowSize,
        )
    else:
        print("Using NumPy for Euler deconvolution (no Numba available)")
        estx, esty, estz, estb, stdzmat = solve_euler_systems_numpy(
            windows_data,
            windows_dx,
            windows_dy,
            windows_dz,
            windows_x,
            windows_y,
            windows_z,
            SI,
            windowSize,
        )

    # Remove border effects (same as original)
    delta = windowSize // 2
    if delta > 0:
        estx = estx[delta:-delta, delta:-delta]
        esty = esty[delta:-delta, delta:-delta]
        estz = estz[delta:-delta, delta:-delta]
        estb = estb[delta:-delta, delta:-delta]
        stdzmat = stdzmat[delta:-delta, delta:-delta]
        xi_trimmed = xi[delta:-delta, delta:-delta]
        yi_trimmed = yi[delta:-delta, delta:-delta]

    # Group solutions and apply filtering (same as original)
    classic = np.stack(
        (estx.ravel(), esty.ravel(), estz.ravel(), estb.ravel(), stdzmat.ravel()),
        axis=-1,
    )

    # Sort and filter solutions
    sorted_indices = np.argsort(classic[:, -1])[::-1]  # Sort by std (descending)
    n_keep = int(len(classic) * filt)
    classic_est = classic[sorted_indices[:n_keep], :-1]
    print("Euler finished")
    return classic_est


# Wrapper function to maintain compatibility with original API
def euler_deconv_opt(data, xi, yi, zi, shape, area, SI, windowSize, filt):
    """
    Drop-in replacement for the original euler_deconv function
    Uses optimized implementation automatically
    """

    return euler_deconv_optimized(data, xi, yi, zi, shape, area, SI, windowSize, filt)
