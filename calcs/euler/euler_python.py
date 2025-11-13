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


def euler_deconv(data, xi, yi, zi, shape, area, SI, windowSize, filt):
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

    Returns:

    * classic_est : 2d-array
        x, y, z and base-level best estimates kept after select a percentage

    * classic : 2d-array
        x, y, z, base-level and standard deviation of all estimates
    """
    data = data.reshape(shape)
    dx, dy, dz = deriv(data, shape, area)

    xi = xi.reshape(shape)
    yi = yi.reshape(shape)
    zi = zi.reshape(shape)

    delta = windowSize // 2
    estx = np.zeros_like(data)
    esty = np.zeros_like(data)
    estz = np.zeros_like(data)
    estb = np.zeros_like(data)
    stdzmat = np.zeros_like(data)

    # run the moving data window and perform the computations
    for east, south, windata, windx, windy, windz, winx, winy, winz in moving_window(
        data, dx, dy, dz, xi, yi, zi, (windowSize, windowSize)
    ):
        # to keep the same size of the window throughout the grid
        if windata.shape[0] != windowSize or windata.shape[1] != windowSize:
            continue
        # system of equations on Euler deconvolution
        A = np.zeros((windowSize * windowSize, 4))
        A[:, 0] = windx.ravel()
        A[:, 1] = windy.ravel()
        A[:, 2] = windz.ravel()
        A[:, 3] = SI * np.ones_like(winx.ravel())

        vety = np.zeros((windowSize * windowSize, 1))
        vety = (
            windx.ravel() * winx.ravel()
            + windy.ravel() * winy.ravel()
            + windz.ravel() * winz.ravel()
            + SI * windata.ravel()
        )
        # compute the estimates
        ATA = np.linalg.inv(np.dot(A.T, A))
        ATy = np.dot(A.T, vety)
        p = np.dot(ATA, ATy)

        # standard deviation of z derivative (for populations population)
        stdz = np.sqrt(
            np.sum(abs(A[:, 2] - A[:, 2].mean()) ** 2) / (len(A[:, 2]) - 1.0)
        )

        estx[south + windowSize // 2][east + windowSize // 2] = p[0]
        esty[south + windowSize // 2][east + windowSize // 2] = p[1]
        estz[south + windowSize // 2][east + windowSize // 2] = p[2]
        estb[south + windowSize // 2][east + windowSize // 2] = p[3]
        stdzmat[south + windowSize // 2][east + windowSize // 2] = stdz

    # get rid of zeros in the border
    estx = estx[delta:-delta, delta:-delta]
    esty = esty[delta:-delta, delta:-delta]
    estz = estz[delta:-delta, delta:-delta]
    estb = estb[delta:-delta, delta:-delta]
    stdzmat = stdzmat[delta:-delta, delta:-delta]
    xi = xi[delta:-delta, delta:-delta]
    yi = yi[delta:-delta, delta:-delta]
    # group the solutions for the classic plot
    classic = np.stack(
        (estx.ravel(), esty.ravel(), estz.ravel(), estb.ravel(), stdzmat.ravel()),
        axis=-1,
    )
    # sort the solutions according to the std of df/dz and filter a percentage
    classic_est = np.array(sorted(classic, key=lambda l: l[-1], reverse=True))[
        : int(len(classic) * filt), :-1
    ]
    return classic_est
