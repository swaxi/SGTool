import numpy as np


def write_rows(f, array, fmt, delimiter=",", chunk=100000):
    """
    Write the rows of a 2D array as delimited text to an open file, formatting
    each chunk of rows with one string operation (much faster than np.savetxt,
    which formats row by row). fmt is one format for every column, e.g. "%s"
    (the same text np.savetxt writes for fmt="%s") or "%.17g" (the shortest
    fixed width that reads back to exactly the same float64 values).
    """
    array = np.asarray(array)
    if array.ndim == 1:
        array = array[:, None]
    row = delimiter.join([fmt] * array.shape[1]) + "\n"
    for i in range(0, len(array), chunk):
        block = array[i : i + chunk]
        f.write((row * len(block)) % tuple(block.ravel().tolist()))


def savetxt(path, array, fmt="%.17g", delimiter=",", header=None):
    """Write a 2D array as delimited text, with an optional header line."""
    with open(path, "w") as f:
        if header:
            f.write(header + "\n")
        write_rows(f, array, fmt, delimiter)
