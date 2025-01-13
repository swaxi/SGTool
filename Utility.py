import numpy as np
from matplotlib import pyplot as plt
import sys
import os

sys.path += [
    "/Users/frank/Documents/Src/Git Stuff/BSDWormer/src/Py3Vtk/pyvtk/build/lib/pyvtk"
]
# import pyvtk as
from osgeo import gdal, ogr, osr

from qgis.core import (
    QgsProject,
    QgsRasterLayer,
)
from qgis.PyQt.QtCore import (
    QFileInfo,
)


def viewRaster(numpy_grid):
    plt.imshow(numpy_grid)
    plt.show()


def writeGDALRasterFromNumpyArray(dst_filename, array, geotransform, proj):
    driver = gdal.GetDriverByName("GTiff")
    y_pixels, x_pixels = array.shape

    dataset = driver.Create(
        dst_filename,
        x_pixels,
        y_pixels,
        1,
        gdal.GDT_Float32,
    )

    dataset.SetGeoTransform(geotransform)
    dataset.SetProjection(proj)
    dataset.GetRasterBand(1).WriteArray(array)
    dataset.FlushCache()  # Write to disk.


"""
def writeVtkImage(filename, image, origin, spacing):
    image_points = [val for val in np.ravel(image)]
    image_coords = [image.shape[1], image.shape[0], 1]
    vtk = PV.VtkData(
        PV.StructuredPoints(image_coords, origin=origin, spacing=spacing),
        "Image",
        PV.PointData(PV.Scalars(image_points, name="field vals")),
        "Image Data",
    )
    vtk.tofile(filename + "_image.vtk")


def writeVtkWorms(filename, points, lines, vals, single_level=True):
    if not single_level:
        # Warning. This code is buggy and produces a messed up vtk file
        # It displays, but there are all kinds of problems in it
        # FIXME: cure this stanza or junk it...
        points_list = [p for v in points.values() for p in v]
        lines_list = [p for v in lines.values() for p in v]
        vals_list = [p for v in vals.values() for p in v]
    else:
        points_list = points
        lines_list = lines
        vals_list = vals

    vtk = PV.VtkData(
        PV.PolyData(points=points_list, lines=lines_list),
        "Worm Segments",
        PV.PointData(PV.Scalars(np.fabs(vals_list), name="mutliscale edge magnitudes")),
        "Edge Magnitudes",
    )
    vtk.tofile(filename + "_worms.vtk")


def writeVtkWormLevels(filename, points, lines, vals):
    assert len(points) == len(lines)
    assert len(lines) == len(vals)
    for level in range(1, len(points) + 1):
        pts = points[level]
        lns = lines[level]
        vs = vals[level]
        name = filename + "_level_{0}".format(level)
        writeVtkWorms(name, pts, lns, vs)
"""


def isclose(a, b, rtol=1.0e-5, atol=1.0e-8, check_invalid=True):
    """Similar to numpy.allclose, but returns a boolean array.
    See numpy.allclose for an explanation of *rtol* and *atol*.

    # A few quick tests...
    >>> assert np.any(isclose(0.300001, np.array([0.1, 0.2, 0.3, 0.4])))

    >>> x = np.array([0.1, np.nan, np.inf, -np.inf])
    >>> y = np.array([0.1000001, np.nan, np.inf, -np.inf])
    >>> assert np.all(isclose(x, y))

    >>> x = np.array([0.1, 0.2, np.inf])
    >>> y = np.array([0.101, np.nan, 0.2])
    >>> assert not np.all(isclose(x, y))
    """

    def within_tol(x, y, atol, rtol):
        return np.less_equal(np.abs(x - y), atol + rtol * np.abs(y))

    x = np.asarray(a)
    y = np.asarray(b)
    if not check_invalid:
        return within_tol(x, y, atol, rtol)
    xfin = np.isfinite(x)
    yfin = np.isfinite(y)
    if np.all(xfin) and np.all(yfin):
        return within_tol(x, y, atol, rtol)
    else:
        # Avoid subtraction with infinite/nan values...
        cond = np.zeros(np.broadcast(x, y).shape, dtype=bool)
        mask = xfin & yfin
        cond[mask] = within_tol(x[mask], y[mask], atol, rtol)
        # Inf and -Inf equality...
        cond[~mask] = x[~mask] == y[~mask]
        # NaN equality...
        cond[np.isnan(x) & np.isnan(y)] = True
        return cond


"""The following code comes from 
<https://gis.stackexchange.com/questions/57834/how-to-get-raster-corner-coordinates-using-python-gdal-bindings>
and supposedly comes from the metageta project <https://code.google.com/p/metageta/>
which is MIT licensed but with the "MetaGETA name substituted for MIT..
I reproduce the text of the MIT license here to (hopefully) remain in 
compliance with the terms of that license.

MetaGETA license

Copyright (c) 2013 Australian Government, Department of the Environment

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""


def GetExtent(gt, cols, rows):
    """Return list of corner coordinates from a geotransform

    @type gt:   C{tuple/list}
    @param gt: geotransform
    @type cols:   C{int}
    @param cols: number of columns in the dataset
    @type rows:   C{int}
    @param rows: number of rows in the dataset
    @rtype:    C{[float,...,float]}
    @return:   coordinates of each corner
    """
    ext = []
    xarr = [0, cols]
    yarr = [0, rows]

    for px in xarr:
        for py in yarr:
            x = gt[0] + (px * gt[1]) + (py * gt[2])
            y = gt[3] + (px * gt[4]) + (py * gt[5])
            ext.append([x, y])
            # print(x, y)
        yarr.reverse()
    return ext


def ReprojectCoords(coords, src_srs, tgt_srs):
    """Reproject a list of x,y coordinates.

    @type geom:     C{tuple/list}
    @param geom:    List of [[x,y],...[x,y]] coordinates
    @type src_srs:  C{osr.SpatialReference}
    @param src_srs: OSR SpatialReference object
    @type tgt_srs:  C{osr.SpatialReference}
    @param tgt_srs: OSR SpatialReference object
    @rtype:         C{tuple/list}
    @return:        List of transformed [[x,y],...[x,y]] coordinates
    """
    trans_coords = []
    transform = osr.CoordinateTransformation(src_srs, tgt_srs)
    for x, y in coords:
        x, y, z = transform.TransformPoint(x, y)
        trans_coords.append([x, y])
    return trans_coords


""" This is converted into a doctest...
    >>> raster=r'somerasterfile.tif'
    >>> ds=gdal.Open(raster)

    >>> gt=ds.GetGeoTransform()
    >>> cols = ds.RasterXSize
    >>> rows = ds.RasterYSize
    >>> ext=GetExtent(gt,cols,rows)

    >>> src_srs=osr.SpatialReference()
    >>> src_srs.ImportFromWkt(ds.GetProjection())
    >>> #tgt_srs=osr.SpatialReference()
    >>> #tgt_srs.ImportFromEPSG(4326)
    >>> tgt_srs = src_srs.CloneGeogCS()

    >>> geo_ext=ReprojectCoords(ext,src_srs,tgt_srs)
"""


def is_layer_loaded(layer_name):
    """
    Check if a layer with the specified name is already loaded in QGIS.

    Parameters:
        layer_name (str): The name of the layer to check.

    Returns:
        bool: True if the layer is loaded, False otherwise.
    """
    for layer in QgsProject.instance().mapLayers().values():
        if layer.name() == layer_name:
            return True
    return False


def loadGrid(diskGridPath):
    fileInfo = QFileInfo(diskGridPath)
    baseName = fileInfo.baseName()

    layer = QgsRasterLayer(diskGridPath, baseName)
    if not is_layer_loaded(baseName):
        QgsProject.instance().addMapLayer(layer)

    dx = layer.rasterUnitsPerPixelX()
    dy = layer.rasterUnitsPerPixelY()
    # Access the raster data provider
    provider = layer.dataProvider()

    # Get raster dimensions
    cols = provider.xSize()  # Number of columns
    rows = provider.ySize()  # Number of rows

    # Read raster data as a block
    band = 1  # Specify the band number (1-based index)
    raster_block = provider.block(band, provider.extent(), cols, rows)

    # Copy the block data into a NumPy array
    extent = layer.extent()
    rows, cols = layer.height(), layer.width()
    raster_block = provider.block(1, extent, cols, rows)  # !!!!!
    raster_array = np.zeros((rows, cols))
    for i in range(rows):
        for j in range(cols):
            raster_array[i, j] = raster_block.value(i, j)

    # Handle NoData values if needed
    no_data_value = provider.sourceNoDataValue(1)  # Band 1

    if no_data_value is not None:
        raster_array[raster_array == no_data_value] = np.nan
    print("mean", np.nanmean(raster_array))
    return raster_array, layer


def numpy_array_to_raster(
    numpy_array,
    raster_path,
    pad_x,
    pad_y,
    dx=None,
    xmin=None,
    ymax=None,
    reference_layer=None,
    no_data_value=np.nan,
):
    """
    Convert a NumPy array to a GeoTIFF raster file.

    Parameters:
        numpy_array (numpy.ndarray): The NumPy array to convert.
        raster_path (str): The path to save the raster file.
        reference_layer (QgsRasterLayer, optional): A reference layer for CRS and geotransform.
        no_data_value: Value to use for no data (default is NaN).
    """

    # Check if the file already exists and remove it
    if os.path.exists(raster_path):
        try:
            os.remove(raster_path)
        except:
            print(
                "Couldn't delete layer, may be open in another program? On windows files on non-C: drive may be hard to delete"
            )
            return -1

    rows, cols = numpy_array.shape
    driver = gdal.GetDriverByName("GTiff")
    output_raster = driver.Create(raster_path, cols, rows, 1, gdal.GDT_Float32)

    # Set geotransform and projection if a reference layer is provided
    if reference_layer:

        provider = reference_layer.dataProvider()
        extent = provider.extent()
        dx = extent.width() / (cols - pad_x * 2)
        dy = -extent.height() / (rows - pad_y * 2)

        geotransform = [
            extent.xMinimum() - (pad_x * dx),
            dx,  # pixel width
            0,
            extent.yMaximum() - (pad_y * dy),
            0,
            dy,  # pixel height (negative)
        ]
        output_raster.SetGeoTransform(geotransform)

        # Set CRS
        srs = osr.SpatialReference()
        srs.ImportFromWkt(reference_layer.crs().toWkt())
        output_raster.SetProjection(srs.ExportToWkt())

    # Write data to raster
    band = output_raster.GetRasterBand(1)
    if no_data_value is not None:
        band.SetNoDataValue(no_data_value)
    numpy_array = np.nan_to_num(
        numpy_array, nan=no_data_value
    )  # Replace NaN with no_data_value
    band.WriteArray(numpy_array)
    band.FlushCache()
    output_raster.FlushCache()
    del output_raster
    return 0


def insert_text_before_extension(file_path, insert_text):
    """
    Insert text at the end of the filename, before the file extension.

    Parameters:
        file_path (str): Full path of the file.
        insert_text (str): Text to insert before the file extension.

    Returns:
        str: The modified file path.
    """
    # Separate the file path into directory, base name, and extension
    dir_name, base_name = os.path.split(file_path)
    file_name, file_ext = os.path.splitext(base_name)

    # Construct the new file name
    new_file_name = f"{file_name}{insert_text}{file_ext}"

    # Combine directory and new file name
    return os.path.join(dir_name, new_file_name)


def fill_nan(data):
    """
    Replace NaN values with the mean of the non-NaN values.
    """
    nan_mask = np.isnan(data)
    filled_data = np.copy(data)
    filled_data[nan_mask] = 0
    return filled_data


if __name__ == "__main__":
    import doctest

    doctest.testmod()
