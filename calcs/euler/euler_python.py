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
from numpy.lib.stride_tricks import sliding_window_view


def fft_pad_data(data, mode="edge"):
    """
    Pad data and compute the coeficients in Fourier domain
    The data is padded until reach the length of the next higher power
    of two and the values of the pad are the values of the edge

    Parameters:

    * data: 2d-array
        the input data set - gridded

    Returns:

    * fpdat: 2d-array
        the coefficients of the data in Fourier domain
    * mask: 2d-array
        Location of padding points -
             {True: data points.
              False: padded points.}
    """
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
    """
    Computes the inverse Fourier Transform of a padded array and mask
    the data to the original shape.

    Parameters:

    * data_p: 2d-array
        Array with the padded data.
    * mask: 2d-array
        Location of padding points -
             {True: Points to be kept .
              False: Points to be removed.}
    * shape_dat: tube = (ny, nx)
        The number of data points in each direction before padding.

    Returns:

    * data: 2d-array
        The unpadded data in space-domain.
    """
    ifft_data = np.real(np.fft.ifft2(data_p))
    data = ifft_data[mask]
    return np.reshape(data, shape_dat)


def fft_wavenumbers(x, y, shape, padshape):
    """
    Computes the wavenumbers 2d-arrays

    Parameters:

    * x,y: 2d-array
        grid of the coordinates.
    * shape: tuple = (ny, nx)
        the number of data points in each direction before padding.
    * padshape: tuple = (ny, nx)
        the number of data points in each direction after padding.

    Returns:

    * u,v: 2d-array
        wavenumbers in each direction
    """

    nx, ny = shape
    dx = (x.max() - x.min()) / (nx - 1)
    u = 2 * np.pi * np.fft.fftfreq(padshape[0], dx)
    dy = (y.max() - y.min()) / (ny - 1)
    v = 2 * np.pi * np.fft.fftfreq(padshape[1], dy)
    return np.meshgrid(v, u)[::-1]


def deriv(data, shape, area):
    """
    Compute the first derivative of a potential field
    in Fourier domain in the x, y and z directions.

    Parameters:

    * data: 2d-array
        the input data set - gridded
    * shape : tuple = (nx, ny)
        the shape of the grid
    * area : list
        the area of the input data - [south, north, west, east]

    Returns:

    * derivx, derivy, derivz : 2D-array
        derivatives in x-, y- and z-directions
    """

    anom_FFT, mask = fft_pad_data(data)

    nx, ny = shape
    xa, xb, ya, yb = area
    xs = np.linspace(xa, xb, nx)
    ys = np.linspace(ya, yb, ny)
    Y, X = np.meshgrid(ys, xs)

    u, v = fft_wavenumbers(X, Y, data.shape, anom_FFT.shape)

    derivx_ft = anom_FFT * (u * 1j)
    derivy_ft = anom_FFT * (v * 1j)
    derivz_ft = anom_FFT * np.sqrt(u**2 + v**2)
    derivx = ifft_unpad_data(derivx_ft, mask, data.shape)
    derivy = ifft_unpad_data(derivy_ft, mask, data.shape)
    derivz = ifft_unpad_data(derivz_ft, mask, data.shape)

    return derivx, derivy, derivz


def moving_window(data, dx, dy, dz, xi, yi, zi, windowSize):
    """
    Moving data window that selects the data, derivatives and coordinates
    for solve the system of Euler deconvolution.
    For a 2d-array, the window runs from left to right and up to down
    The window moves 1 step for iteration

    Parameters:

    * data : 2d-array
        the input data set - gridded
    * dx, dy, dz : 2d-array
        derivatives in x-, y- and z-directions
    * xi, yi, zi : 2d-array
        grid of coordinates in x-, y- and z-directions
    * windowSize : tuple (x,y)
        size of the window - equal in both directions

    Returns:

    * data : 2d-array
        windowed input data set
    * dx, dy, dz : 2d-array
        windowed derivatives in x-, y- and z-directions
    * xi, yi, zi : 2d-array
        windowed grid of coordinates in x-, y- and z-directions
    """
    for y in range(0, data.shape[0]):
        for x in range(0, data.shape[1]):
            # yield the current window
            yield (
                x,
                y,
                data[y : y + windowSize[1], x : x + windowSize[0]],
                dx[y : y + windowSize[1], x : x + windowSize[0]],
                dy[y : y + windowSize[1], x : x + windowSize[0]],
                dz[y : y + windowSize[1], x : x + windowSize[0]],
                xi[y : y + windowSize[1], x : x + windowSize[0]],
                yi[y : y + windowSize[1], x : x + windowSize[0]],
                zi[y : y + windowSize[1], x : x + windowSize[0]],
            )


def euler_deconv(data, xi, yi, zi, shape, area, SI, windowSize, filt, callback=None):
    """
    Euler deconvolution - solves the system of equations
    for each moving data window

    Parameters:

    * data : 1d-array
        the input data set
    * xi, yi, zi : 1d-array
        grid of coordinates in x-, y- and z-directions
    * shape : tuple = (nx, ny)
        the shape of the grid
    * area : list
        the area of the input data - [south, north, west, east]
    * SI : int
        structural index - 0, 1, 2 or 3
    * windowSize : tuple (dx,dy)
        size of the window - equal in both directions
    * filt : float
        percentage of the solutions that will be keep
    * callback : callable, optional
        called as callback(fraction_done) once per row of windows; it may
        raise calcs.sgt_cancel.OperationCancelled to stop the calculation

    Returns:

    * classic_est : 2d-array
        x, y, z, base-level and std of df/dz of the best estimates kept after
        select a percentage

    * classic : 2d-array
        x, y, z, base-level and standard deviation of all estimates
    """
    data = data.reshape(shape)
    dx, dy, dz = deriv(data, shape, area)

    xi = xi.reshape(shape)
    yi = yi.reshape(shape)
    zi = zi.reshape(shape)

    # the system for each window is solved through its normal equations
    # (A^T A) p = A^T y, whose terms are sums over the window; these sums are
    # computed for all windows at once, in blocks of rows to bound memory
    delta = windowSize // 2
    ny, nx = data.shape
    n_out_y, n_out_x = ny - 2 * delta, nx - 2 * delta
    n = float(windowSize * windowSize)
    vety = dx * xi + dy * yi + dz * zi + SI * data
    # products of the columns of A = [dx, dy, dz, SI] and y, as summed below
    terms = (
        dx * dx, dx * dy, dx * dz, dx,
        dy * dy, dy * dz, dy,
        dz * dz, dz,
        dx * vety, dy * vety, dz * vety, vety,
    )

    def window_sums(a):
        # sums over every windowSize x windowSize window (top-left anchored)
        s = sliding_window_view(a, windowSize, axis=1).sum(axis=-1)
        return sliding_window_view(s, windowSize, axis=0).sum(axis=-1)

    estx = np.empty((n_out_y, n_out_x))
    esty = np.empty_like(estx)
    estz = np.empty_like(estx)
    estb = np.empty_like(estx)
    stdzmat = np.empty_like(estx)

    block = max(1, 2000000 // max(nx, 1))  # windows rows per block
    for r0 in range(0, n_out_y, block):
        if callback is not None:
            callback(r0 / float(n_out_y))
        r1 = min(r0 + block, n_out_y)
        rows = slice(r0, r1 + windowSize - 1)
        (sxx, sxy, sxz, sx, syy, syz, sy, szz, sz,
         sxv, syv, szv, sv) = (window_sums(t[rows])[: r1 - r0, :n_out_x]
                               for t in terms)

        ATA = np.empty(sxx.shape + (4, 4))
        ATA[..., 0, 0] = sxx
        ATA[..., 0, 1] = ATA[..., 1, 0] = sxy
        ATA[..., 0, 2] = ATA[..., 2, 0] = sxz
        ATA[..., 0, 3] = ATA[..., 3, 0] = SI * sx
        ATA[..., 1, 1] = syy
        ATA[..., 1, 2] = ATA[..., 2, 1] = syz
        ATA[..., 1, 3] = ATA[..., 3, 1] = SI * sy
        ATA[..., 2, 2] = szz
        ATA[..., 2, 3] = ATA[..., 3, 2] = SI * sz
        ATA[..., 3, 3] = SI * SI * n
        ATy = np.stack((sxv, syv, szv, SI * sv), axis=-1)[..., None]

        # a singular window (e.g. a flat patch) has no solution: give it zero
        # estimates and a zero std of df/dz, so it is the first to be dropped
        det = np.linalg.det(ATA)
        bad = ~np.isfinite(det) | (det == 0)
        ATA[bad] = np.eye(4)
        p = np.linalg.solve(ATA, ATy)[..., 0]
        p[bad] = 0.0

        # standard deviation of the z derivative (sample, as in the original)
        stdz = np.sqrt(np.maximum(szz - sz * sz / n, 0.0) / (n - 1.0))
        stdz[bad] = 0.0

        estx[r0:r1] = p[..., 0]
        esty[r0:r1] = p[..., 1]
        estz[r0:r1] = p[..., 2]
        estb[r0:r1] = p[..., 3]
        stdzmat[r0:r1] = stdz

    # group the solutions for the classic plot
    classic = np.stack(
        (estx.ravel(), esty.ravel(), estz.ravel(), estb.ravel(), stdzmat.ravel()),
        axis=-1,
    )
    # sort the solutions according to the std of df/dz and filter a percentage;
    # the std of df/dz is kept as the last column for secondary filtering
    order = np.argsort(-classic[:, -1], kind="stable")
    classic_est = classic[order[: int(len(classic) * filt)], :]
    return classic_est
