# -*- coding: utf-8 -*-
# Copyright (c) 2013 Australian Government, Department of the Environment
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

"""Geometry and dataset helper functions

@todo: Does the module name I{really} describe its functions anymore?
       It doesn't just do geometry now...
"""
"""from __future__ import division
from future.builtins import round
from future.builtins import int
from future.builtins import str
from future.builtins import range
from future.builtins import map"""

import os, math, warnings, tempfile, re

# from metageta import utilities
try:
    from osgeo import gdal
    from osgeo import gdalconst
    from osgeo import osr
    from osgeo import ogr
except ImportError:
    import gdal
    import gdalconst
    import osr
    import ogr

# ========================================================================================================
# {String Utilities
# ========================================================================================================
encoding = "utf-8"


def encode(string):
    """Encode a unicode string
    @type     string:  C{unicode}
    @param    string:  Unicode string
    @rtype:   C{str}
    @return:  Encoded string
    """
    # We're in Python3 mode here
    # if type(string) is unicode:return string.encode(encoding)
    # else:return string
    return string.encode(encoding)


def runcmd(cmd, format="s"):
    """Run a command
    @type     cmd:  C{str}
    @param    cmd:  Command (inc arguments) to run
    @rtype:   C{tuple}
    @return:  Returns (exit_code,stdout,stderr)
    """
    import subprocess

    proc = subprocess.Popen(
        cmd,
        shell=True,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if format.lower() == "s":  # string output
        stdout, stderr = proc.communicate()
    # elif format.lower() == 'f': #file object output #doesn't flush IO buffer, causes python to hang
    #    stdout,stderr=proc.stdout,proc.stderr
    elif format.lower() == "l":  # list output
        stdout, stderr = proc.stdout.readlines(), proc.stderr.readlines()
    # else:raise TypeError, "fomat argument must be in ['s','f','l'] (string, file, list)"
    else:
        raise TypeError("fomat argument must be in ['s','l'] (string or list format)")
    exit_code = proc.wait()
    return exit_code, stdout, stderr


debug = False


# ========================================================================================================
# {Custom error class
# ========================================================================================================
class GDALError(Exception):
    """For raising GDAL related errors"""

    def __init__(self, msg=None):
        """For raising GDAL related errors

        @type msg: C{str}
        @param msg: Initial message to include with GDAL the error message
        @return None

        @todo: GDAL now has a UseExceptions method, should this class go away?
        """
        errtype = {
            gdal.CE_None: "None",
            gdal.CE_Debug: "Debug",
            gdal.CE_Warning: "Warning",
            gdal.CE_Failure: "Failure",
            gdal.CE_Fatal: "Fatal",
        }
        self.errmsg = gdal.GetLastErrorMsg().replace("\n", " ")
        self.errnum = gdal.GetLastErrorNo()
        self.errtyp = errtype.get(gdal.GetLastErrorType(), "None")
        gdal.ErrorReset()
        if msg:
            self.errmsg = "\n".join([msg, self.errmsg])
        Exception.__init__(self, self.errmsg, self.errnum, self.errtyp)

    def __str__(self):
        """For printing GDAL related errors

        @rtype:  C{str}
        @return: String representation of the exception
        """

        return "%s\nError Message:%s" % (self.errtyp, self.errmsg)


# ========================================================================================================
# {Dataset Utilities
# ========================================================================================================
def OpenDataset(filepath, mode=gdalconst.GA_ReadOnly):
    """Open & return a gdalDataset object

    @type    filepath: C{str}
    @param   filepath: path to dataset
    @rtype:  C{osgeo.gdal.Dataset}
    @return: GDAL Dataset
    """

    gdal.ErrorReset()
    gdal.PushErrorHandler("CPLQuietErrorHandler")
    gdalDataset = gdal.Open(filepath, mode)
    if not gdalDataset:
        raise GDALError("Unable to open %s" % filepath)
    gdal.PopErrorHandler()
    return gdalDataset


def ParseGDALinfo(filepath):
    """Very basic gdalinfo parser, does not include colo(u)r tables, raster attribute tables.

    @type    filepath: C{str}
    @param   filepath: path to dataset
    @rtype:  C{(dict,list)}
    @return: GDAL info,extent coordinates
    """

    metadata = {}
    extent = []

    cmd = "gdalinfo -noct " + filepath
    exit_code, stdout, stderr = runcmd(cmd)
    if exit_code != 0:
        raise Exception(stderr)

    decimal = r"([-+]?\d+\.?\d+e?[-+]?\d+?)"

    rex = r"Size is (\d+),\s*(\d+)"
    rex = re.compile(rex, re.I)
    rex = rex.findall(stdout)
    if rex:
        metadata["cols"], metadata["rows"] = [int(r) for r in rex[0]]

    rex = r"Coordinate System is:\s*(.*\])"
    rex = re.compile(rex, re.I)
    rex = rex.findall(" ".join([s.strip() for s in stdout.split("\n")]))
    if rex:
        metadata["srs"] = rex[0]

    rex = r"Pixel Size\s*=\s*\(%s,%s\)" % (decimal, decimal)
    rex = re.compile(rex, re.I)
    rex = rex.findall(stdout)
    if rex:
        metadata["cellx"], metadata["celly"] = [float(r) for r in rex[0]]

    rex = r"Upper Left\s*\(\s*%s\s*,\s*%s\)" % (decimal, decimal)
    rex = re.compile(rex, re.I)
    rex = rex.findall(stdout)
    if rex:
        extent.append([float(r) for r in rex[0]])

    rex = r"Upper Right\s*\(\s*%s\s*,\s*%s\)" % (decimal, decimal)
    rex = re.compile(rex, re.I)
    rex = rex.findall(stdout)
    if rex:
        extent.append([float(r) for r in rex[0]])

    rex = r"Lower Left\s*\(\s*%s\s*,\s*%s\)" % (decimal, decimal)
    rex = re.compile(rex, re.I)
    rex = rex.findall(stdout)
    if rex:
        extent.append([float(r) for r in rex[0]])

    rex = r"Lower Right\s*\(\s*%s\s*,\s*%s\)" % (decimal, decimal)
    rex = re.compile(rex, re.I)
    rex = rex.findall(stdout)
    if rex:
        extent.append([float(r) for r in rex[0]])

    rex = r"Band \d"
    rex = re.compile(rex, re.I)
    rex = rex.findall(stdout)
    if rex:
        metadata["nbands"] = len(rex)

    rex = r"Type=(%s)" % "|".join(
        [gdal.GetDataTypeName(dt) for dt in range(0, gdal.GDT_TypeCount)]
    )
    rex = re.compile(rex, re.I)
    rex = rex.findall(stdout)
    if rex:
        metadata["datatype"] = rex[0]

    rex = r"NoData Value=" + decimal
    rex = re.compile(rex, re.I)
    rex = rex.findall(stdout)
    if rex:
        metadata["nodata"] = rex[0]

    return metadata, extent


# ========================================================================================================
# {Coordinate Utilities
# ========================================================================================================
def DMS2DD(dms, format):
    """Convert coordinate in degrees, minutes, seconds to decimal degrees.

    @type  dms: str
    @param dms: degrees, minutes, seconds
    @type  format: str
    @param format: The format parameter is a string of equal length to the DMS string
    and identifies the position of each of the following elements. Any other string in
    the format parameter is ignored.

        - H - hemisphere
        - D - degrees
        - M - minutes
        - S - seconds

    eg:
        - C{'27 45 12 E'   = 'DD MM SS H'}
        - C{'027 45 12 E'  = 'DDD MM SS H'}
        - C{'027-45-12-E'  = 'DDD MM SS H'}
        - C{'27,45,12.3 E' = 'DD MM SSSS H'}
        - C{'-27,45,12.3'  = 'DDD MM SSSS'}

    @rtype: C{float}
    @return: decimal degrees
    """
    if len(dms) != len(format):
        raise ValueError(
            "Format string %s does not match coordinate string %s" % (dms, format)
        )
    H = ""
    D = ""
    M = ""
    S = ""
    for i, f in enumerate(format):
        if f == "H":
            H += dms[i]
        elif f == "D":
            D += dms[i]
        elif f == "M":
            M += dms[i]
        elif f == "S":
            S += dms[i]
    dd = float(D)
    if dd != abs(dd):  # negative value
        dd = abs(dd)
        neg = True
    else:
        neg = False
    if M:
        dd += float(M) / 60.0  # May not have been passed Minutes or Seconds
    if S:
        dd += float(S) / 3600.0
    if neg or (H and H in ["S", "W"]):
        dd = dd * -1.0

    return dd


# ========================================================================================================
# {Geometry Utilities
# ========================================================================================================
def Rotation(gt):
    """Get rotation angle from a geotransform
    @type gt: C{tuple/list}
    @param gt: geotransform
    @rtype: C{float}
    @return: rotation angle
    """
    try:
        return math.degrees(math.tanh(gt[2] / gt[5]))
    except:
        return 0


def CellSize(gt):
    """Get cell size from a geotransform

    @type gt:  C{tuple/list}
    @param gt: geotransform
    @rtype:    C{(float,float)}
    @return:   (x,y) cell size
    """
    cellx = round(math.hypot(gt[1], gt[4]), 7)
    celly = round(math.hypot(gt[2], gt[5]), 7)
    return (cellx, celly)


def SceneCentre(gt, cols, rows):
    """Get scene centre from a geotransform.

    @type gt: C{tuple/list}
    @param gt: geotransform
    @type cols: C{int}
    @param cols: Number of columns in the dataset
    @type rows: C{int}
    @param rows: Number of rows in the dataset
    @rtype:    C{(float,float)}
    @return:   Scene centre coordinates
    """
    px = cols / 2
    py = rows / 2
    x = gt[0] + (px * gt[1]) + (py * gt[2])
    y = gt[3] + (px * gt[4]) + (py * gt[5])
    return x, y


def ExtentToGCPs(ext, cols, rows):
    """Form a gcp list from the 4 corners.

    This function is meant to be used to convert an extent
    to gcp's for use in the gdal.GCPsToGeoTransform function.

    @type ext:   C{tuple/list}
    @param ext:  Extent, must be in order: [[ulx,uly],[urx,ury],[lrx,lry],[llx,lly]]
    @type cols: C{int}
    @param cols: Number of columns in the dataset
    @type rows: C{int}
    @param rows: Number of rows in the dataset
    @rtype:    C{[gcp,...,gcp]}
    @return:   List of GCP objects
    """
    gcp_list = []
    parr = [0, cols]
    larr = [rows, 0]
    id = 0
    if len(ext) == 5:  # Assume ext[0]==ext[4]
        ext = ext[:-1]
    if len(ext) != 4:
        raise ValueError("Extent must be a tuple/list with 4 elements, each an XY pair")

    for px in parr:
        for py in larr:
            cgcp = gdal.GCP()
            cgcp.Id = str(id)
            cgcp.GCPX = ext[id][0]
            cgcp.GCPY = ext[id][1]
            cgcp.GCPZ = 0.0
            cgcp.GCPPixel = px
            cgcp.GCPLine = py
            id += 1
            gcp_list.append(cgcp)
        larr.reverse()

    return gcp_list


def GeoTransformToGCPs(gt, cols, rows):
    """Form a gcp list from a geotransform using the 4 corners.

    This function is meant to be used to convert a geotransform
    to gcp's so that the geocoded information can be reprojected.

    @type gt:   C{tuple/list}
    @param gt: geotransform to convert to gcps
    @type cols:   C{int}
    @param cols: number of columns in the dataset
    @type rows:   C{int}
    @param rows: number of rows in the dataset
    @rtype:    C{[gcp,...,gcp]}
    @return:   List of GCP objects
    """
    ###############################################################################
    # This code is modified from the GeoTransformToGCPs function
    # in the OpenEV module vrtutils.py
    ###############################################################################
    # $Id: vrtutils.py,v 1.17 2005/07/07 21:36:06 gmwalter Exp $
    #
    # Project:  OpenEV
    # Purpose:  Utilities for creating vrt files.
    # Author:   Gillian Walter, gwal...@atlsci.com
    #
    ###############################################################################
    # Copyright (c) 2000, Atlantis Scientific Inc. (www.atlsci.com)
    #
    # This library is free software; you can redistribute it and/or
    # modify it under the terms of the GNU Library General Public
    # License as published by the Free Software Foundation; either
    # version 2 of the License, or (at your option) any later version.
    #
    # This library is distributed in the hope that it will be useful,
    # but WITHOUT ANY WARRANTY; without even the implied warranty of
    # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    # Library General Public License for more details.
    #
    # You should have received a copy of the GNU Library General Public
    # License along with this library; if not, write to the
    # Free Software Foundation, Inc., 59 Temple Place - Suite 330,
    # Boston, MA 02111-1307, USA.
    ###############################################################################

    gcp_list = []
    parr = [0, cols]
    larr = [0, rows]
    id = 0
    for px in parr:
        for py in larr:
            cgcp = gdal.GCP()
            cgcp.Id = str(id)
            cgcp.GCPX = gt[0] + (px * gt[1]) + (py * gt[2])
            cgcp.GCPY = gt[3] + (px * gt[4]) + (py * gt[5])
            cgcp.GCPZ = 0.0
            cgcp.GCPPixel = px
            cgcp.GCPLine = py
            id += 1
            gcp_list.append(cgcp)
        larr.reverse()
    return gcp_list


def GeomFromExtent(ext, srs=None, srs_wkt=None):
    """Get and OGR geometry object from a extent list

    @type ext:  C{tuple/list}
    @param ext: extent coordinates
    @type srs:  C{str}
    @param srs: SRS WKT string
    @rtype:     C{ogr.Geometry}
    @return:    Geometry object
    """
    if type(ext[0]) is list or type(ext[0]) is tuple:  # is it a list of xy pairs
        if ext[0] != ext[-1]:
            ext.append(ext[0])
        wkt = "POLYGON ((%s))" % ",".join(
            map(" ".join, [list(map(str, i)) for i in ext])
        )
    else:  # it's a list of xy values
        xmin, ymin, xmax, ymax = ext
        template = "POLYGON ((%(minx)f %(miny)f, %(minx)f %(maxy)f, %(maxx)f %(maxy)f, %(maxx)f %(miny)f, %(minx)f %(miny)f))"
        r1 = {"minx": xmin, "miny": ymin, "maxx": xmax, "maxy": ymax}
        wkt = template % r1
    if srs_wkt is not None:
        srs = osr.SpatialReference(wkt=srs_wkt)
    geom = ogr.CreateGeometryFromWkt(wkt, srs)
    return geom


def ReprojectGeom(geom, src_srs, tgt_srs):
    """Reproject a geometry object.

    @type geom:     C{ogr.Geometry}
    @param geom:    OGR geometry object
    @type src_srs:  C{osr.SpatialReference}
    @param src_srs: OSR SpatialReference object
    @type tgt_srs:  C{osr.SpatialReference}
    @param tgt_srs: OSR SpatialReference object
    @rtype:         C{ogr.Geometry}
    @return:        OGRGeometry object
    """
    gdal.ErrorReset()
    gdal.PushErrorHandler("CPLQuietErrorHandler")
    geom.AssignSpatialReference(src_srs)
    geom.TransformTo(tgt_srs)
    err = gdal.GetLastErrorMsg()
    if err:
        warnings.warn(err.replace("\n", " "))
    gdal.PopErrorHandler()
    gdal.ErrorReset()
    return geom


def InvGeoTransform(gt_in):
    """
    ************************************************************************
    *                        InvGeoTransform(gt_in)
    ************************************************************************

    **
    * Invert Geotransform.
    *
    * This function will invert a standard 3x2 set of GeoTransform coefficients.
    *
    * @param  gt_in  Input geotransform (six doubles - unaltered).
    * @return gt_out Output geotransform (six doubles - updated) on success,
    *                None if the equation is uninvertable.
    """
    #    ******************************************************************************
    #    * This code ported from GDALInvGeoTransform() in gdaltransformer.cpp
    #    * as it isn't exposed in the python SWIG bindings until GDAL 1.7
    #    * copyright & permission notices included below as per conditions.
    #
    #    ******************************************************************************
    #    * $Id: gdaltransformer.cpp 15024 2008-07-24 19:25:06Z rouault $
    #    *
    #    * Project:  Mapinfo Image Warper
    #    * Purpose:  Implementation of one or more GDALTrasformerFunc types, including
    #    *           the GenImgProj (general image reprojector) transformer.
    #    * Author:   Frank Warmerdam, warmerdam@pobox.com
    #    *
    #    ******************************************************************************
    #    * Copyright (c) 2002, i3 - information integration and imaging
    #    *                          Fort Collin, CO
    #    *
    #    * Permission is hereby granted, free of charge, to any person obtaining a
    #    * copy of this software and associated documentation files (the "Software"),
    #    * to deal in the Software without restriction, including without limitation
    #    * the rights to use, copy, modify, merge, publish, distribute, sublicense,
    #    * and/or sell copies of the Software, and to permit persons to whom the
    #    * Software is furnished to do so, subject to the following conditions:
    #    *
    #    * The above copyright notice and this permission notice shall be included
    #    * in all copies or substantial portions of the Software.
    #    *
    #    * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
    #    * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    #    * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
    #    * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    #    * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    #    * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    #    * DEALINGS IN THE SOFTWARE.
    #    ****************************************************************************

    # we assume a 3rd row that is [1 0 0]

    # Compute determinate
    det = gt_in[1] * gt_in[5] - gt_in[2] * gt_in[4]

    if abs(det) < 0.000000000000001:
        return

    inv_det = 1.0 / det

    # compute adjoint, and divide by determinate
    gt_out = [0, 0, 0, 0, 0, 0]
    gt_out[1] = gt_in[5] * inv_det
    gt_out[4] = -gt_in[4] * inv_det

    gt_out[2] = -gt_in[2] * inv_det
    gt_out[5] = gt_in[1] * inv_det

    gt_out[0] = (gt_in[2] * gt_in[3] - gt_in[0] * gt_in[5]) * inv_det
    gt_out[3] = (-gt_in[1] * gt_in[3] + gt_in[0] * gt_in[4]) * inv_det

    return gt_out


def ApplyGeoTransform(inx, iny, gt):
    """Apply a geotransform
    @param  inx:       Input x coordinate (double)
    @param  iny:       Input y coordinate (double)
    @param  gt:        Input geotransform (six doubles)

    @return: outx,outy Output coordinates (two doubles)
    """
    outx = gt[0] + inx * gt[1] + iny * gt[2]
    outy = gt[3] + inx * gt[4] + iny * gt[5]
    return (outx, outy)


def MapToPixel(mx, my, gt):
    """Convert map to pixel coordinates
    @param  mx:    Input map x coordinate (double)
    @param  my:    Input map y coordinate (double)
    @param  gt:    Input geotransform (six doubles)
    @return: px,py Output coordinates (two ints)

    @change: changed int(p[x,y]+0.5) to int(p[x,y]) as per http://lists.osgeo.org/pipermail/gdal-dev/2010-June/024956.html
    @change: return floats
    @note:   0,0 is UL corner of UL pixel, 0.5,0.5 is centre of UL pixel
    """
    if gt[2] + gt[4] == 0:  # Simple calc, no inversion required
        px = (mx - gt[0]) / gt[1]
        py = (my - gt[3]) / gt[5]
    else:
        px, py = ApplyGeoTransform(mx, my, InvGeoTransform(gt))
    # return int(px),int(py)
    return px, py


def PixelToMap(px, py, gt):
    """Convert pixel to map coordinates
    @param  px:    Input pixel x coordinate (double)
    @param  py:    Input pixel y coordinate (double)
    @param  gt:    Input geotransform (six doubles)
    @return: mx,my Output coordinates (two doubles)

    @note:   0,0 is UL corner of UL pixel, 0.5,0.5 is centre of UL pixel
    """
    mx, my = ApplyGeoTransform(px, py, gt)
    return mx, my


# ========================================================================================================
# {VRT Utilities
# ========================================================================================================
def CreateVRTCopy(ds):
    """Create a VRT copy of a gdal.Dataset object.

    @type ds:  C{gdal.Dataset}
    @param ds: Dataset object
    @rtype:    C{gdal.Dataset}
    @return:   VRT Dataset object
    """
    try:
        vrtdrv = gdal.GetDriverByName("VRT")
        vrtds = vrtdrv.CreateCopy("", ds)
        return vrtds
    except:
        return None


def CreateMosaicedVRT(
    files, bands, srcrects, dstrects, cols, rows, datatype, relativeToVRT=0
):
    """Create a VRT XML string that mosaics datasets.

    For further info on VRT's, see the U{GDAL VRT Tutorial<http://www.gdal.org/gdal_vrttut.html>}

    @type files:          C{[str,...,str]}
    @param files:         List of files to mosaic
    @type bands:          C{[int,...,int]}
    @param bands:         List of band numbers (1 based). Eg. [1,2,3] will mosaic
                          the first band from each file into the 1st band of the output VRT, etc.

    @type srcrects:       C{[SrcRect,...,SrcRect]}
    @param srcrects:      List of SrcRects, one per file, in image not map units. E.g [[0,0,512,512],...]
                          will be output as <SrcRect xOff="0" yOff="0" xSize="512" ySize="512"/>.
                          The SrcRect allows you to subset your input image.
    @type dstrects:       C{[DstRect,...,DstRect]}
    @param dstrects:      List of DstRects, One per file, in image not map units. E.g [[512,512,1024,1024],...]
                          will be output as <DstRect xOff="512" yOff="512" xSize="1024" ySize="1024"/>
                          The DstRect determines the spatial position of the input image in the mosaic.
    @type cols:           C{int}
    @param cols:          The number of columns in the output mosaic
    @type rows:           C{int}
    @param rows:          The number of rows in the output mosaic
    @type datatype:       C{str}
    @param datatype:      GDAL datatype name. Eg. Byte, Int32, UInt16

    @rtype:               C{xml}
    @return:              VRT XML string
    """
    try:
        vrt = []
        for i, band in enumerate(bands):
            vrt.append('  <VRTRasterBand dataType="%s" band="%s">' % (datatype, i + 1))
            for j, f in enumerate(files):
                vrt.append("    <SimpleSource>")
                vrt.append(
                    '      <SourceFilename relativeToVRT="%s">%s</SourceFilename>'
                    % (relativeToVRT, f)
                )
                vrt.append(
                    '      <SourceProperties RasterXSize="%s" RasterYSize="%s" DataType="%s"/>'
                    % (dstrects[j][2], dstrects[j][3], datatype)
                )
                vrt.append("      <SourceBand>%s</SourceBand>" % band)
                vrt.append(
                    '      <SrcRect xOff="%s" yOff="%s" xSize="%s" ySize="%s"/>'
                    % (srcrects[j][0], srcrects[j][1], srcrects[j][2], srcrects[j][3])
                )
                vrt.append(
                    '      <DstRect xOff="%s" yOff="%s" xSize="%s" ySize="%s"/>'
                    % (dstrects[j][0], dstrects[j][1], dstrects[j][2], dstrects[j][3])
                )
                vrt.append("    </SimpleSource>")
            vrt.append("  </VRTRasterBand>")
        return CreateCustomVRT("\n".join(vrt), cols, rows)
    except:
        raise  # return None


def CreateSimpleVRT(bands, cols, rows, datatype, relativeToVRT=0):
    """Create a VRT XML string with a simple source from one or more datasets,
    each dataset will be output as a separate band.

    For further info on VRT's, see the U{GDAL VRT Tutorial<http://www.gdal.org/gdal_vrttut.html>}

    @type bands:          C{[str,...,str]}
    @param bands:         List of files. The first file becomes the first band and so forth.
    @type cols:           C{int}
    @param cols:          The number of columns in the output VRT
    @type rows:           C{int}
    @param rows:          The number of rows in the output VRT
    @type datatype:       C{str}
    @param datatype:      GDAL datatype name. Eg. Byte, Int32, UInt16

    @rtype:               C{xml}
    @return:              VRT XML string
    """
    try:
        vrt = []
        for i, band in enumerate(bands):
            vrt.append('  <VRTRasterBand dataType="%s" band="%s">' % (datatype, i + 1))
            vrt.append("    <SimpleSource>")
            vrt.append(
                '      <SourceFilename relativeToVRT="%s">%s</SourceFilename>'
                % (relativeToVRT, band)
            )
            vrt.append("      <SourceBand>1</SourceBand>")
            vrt.append("    </SimpleSource>")
            vrt.append("  </VRTRasterBand>")
        return CreateCustomVRT("\n".join(vrt), cols, rows)
    except:
        return None


def CreateRawRasterVRT(
    bands,
    cols,
    rows,
    datatype,
    headeroffset=0,
    byteorder=None,
    relativeToVRT=0,
    nodata=None,
):
    """Create RawRaster VRT from one or more _single_ band files

    For further info on VRT's, see the U{GDAL VRT Tutorial<http://www.gdal.org/gdal_vrttut.html>}

    @type bands:          C{[str,...,str]}
    @param bands:         List of files. The first file becomes the first band and so forth.
    @type cols:           C{int}
    @param cols:          The number of columns in the output VRT
    @type rows:           C{int}
    @param rows:          The number of rows in the output VRT
    @type datatype:       C{str}
    @param datatype:      GDAL datatype name. Eg. Byte, Int32, UInt16
    @type headeroffset:   C{int}
    @param headeroffset:  Number of bytes to skip at the start of the file
    @type byteorder:      C{str}
    @param byteorder:     Byte order of the file (MSB or LSB)

    @rtype:               C{xml}
    @return:              VRT XML string
    """

    try:
        vrt = []
        nbits = gdal.GetDataTypeSize(gdal.GetDataTypeByName(datatype))
        for i, band in enumerate(bands):
            vrt.append(
                '  <VRTRasterBand dataType="%s" band="%s" subClass="VRTRawRasterBand">'
                % (datatype, i + 1)
            )
            vrt.append(
                '    <SourceFilename relativeToVRT="%s">%s</SourceFilename>'
                % (relativeToVRT, band)
            )
            vrt.append("    <ImageOffset>%s</ImageOffset>" % (headeroffset))
            vrt.append("    <PixelOffset>%s</PixelOffset>" % (nbits / 8))
            vrt.append("    <LineOffset>%s</LineOffset>" % (nbits / 8 * cols))
            if nodata is not None:
                vrt.append(
                    "    <NoDataValue>%s</NoDataValue>" % (nodata)
                )  # Fix for Issue 17
            if byteorder:
                vrt.append("    <ByteOrder>%s</ByteOrder>" % (byteorder))
            vrt.append("  </VRTRasterBand>")
        return CreateCustomVRT("\n".join(vrt), cols, rows)
    except:
        return None
    return "\n".join(vrt)


def CreateBSQRawRasterVRT(
    filename,
    nbands,
    cols,
    rows,
    datatype,
    nodata=None,
    headeroffset=0,
    byteorder=None,
    relativeToVRT=0,
):
    """Create RawRaster VRT from a BSQ

    BSQ = Band-Sequential or Band-Interleaved

    For further info on VRT's, see the U{GDAL VRT Tutorial<http://www.gdal.org/gdal_vrttut.html>}

    @type filename:       C{str}
    @param filename:      File name
    @type nbands:         C{int}
    @param nbands:        The number of bands in the output VRT (<= nbands in input file)
    @type cols:           C{int}
    @param cols:          The number of columns in the output VRT
    @type rows:           C{int}
    @param rows:          The number of rows in the output VRT
    @type datatype:       C{str}
    @param datatype:      GDAL datatype name. Eg. Byte, Int32, UInt16
    @type nodata:         C{float}
    @param nodata:        No data/Null value
    @type headeroffset:   C{int}
    @param headeroffset:  Number of bytes to skip at the start of the file
    @type byteorder:      C{str}
    @param byteorder:     Byte order of the file (MSB or LSB)

    @rtype:               C{xml}
    @return:              VRT XML string
    """
    try:
        vrt = []
        nbits = gdal.GetDataTypeSize(gdal.GetDataTypeByName(datatype))
        for i in range(nbands):
            vrt.append(
                '  <VRTRasterBand dataType="%s" band="%s" subClass="VRTRawRasterBand">'
                % (datatype, i + 1)
            )
            vrt.append(
                '    <SourceFilename relativeToVRT="%s">%s</SourceFilename>'
                % (relativeToVRT, filename)
            )
            if nodata is not None:
                vrt.append("    <NoDataValue>%s</NoDataValue>" % (nodata))
            vrt.append(
                "    <ImageOffset>%s</ImageOffset>"
                % (headeroffset + nbits / 8 * i * cols * rows)
            )
            vrt.append("    <PixelOffset>%s</PixelOffset>" % (nbits / 8))
            vrt.append("    <LineOffset>%s</LineOffset>" % (nbits / 8 * cols))
            if byteorder:
                vrt.append("    <ByteOrder>%s</ByteOrder>" % (byteorder))
            vrt.append("  </VRTRasterBand>")
        return CreateCustomVRT("\n".join(vrt), cols, rows)
    except:
        return None
    return "\n".join(vrt)


def CreateBILRawRasterVRT(
    filename,
    nbands,
    cols,
    rows,
    datatype,
    nbits,
    nodata=None,
    headeroffset=0,
    byteorder=None,
    relativeToVRT=0,
):
    """Create RawRaster VRT from a BIL

    BIL = Band-Interleaved-by-Line or Row-Interleaved

    For further info on VRT's, see the U{GDAL VRT Tutorial<http://www.gdal.org/gdal_vrttut.html>}

    @type filename:       C{str}
    @param filename:      File name
    @type nbands:         C{int}
    @param nbands:        The number of bands in the output VRT (<= nbands in input file)
    @type cols:           C{int}
    @param cols:          The number of columns in the output VRT
    @type rows:           C{int}
    @param rows:          The number of rows in the output VRT
    @type datatype:       C{str}
    @param datatype:      GDAL datatype name. Eg. Byte, Int32, UInt16
    @type nodata:         C{float}
    @param nodata:        No data/Null value
    @type headeroffset:   C{int}
    @param headeroffset:  Number of bytes to skip at the start of the file
    @type byteorder:      C{str}
    @param byteorder:     Byte order of the file (MSB or LSB)

    @rtype:               C{xml}
    @return:              VRT XML string
    """
    try:
        vrt = []
        nbits = gdal.GetDataTypeSize(gdal.GetDataTypeByName(datatype))
        for i in range(nbands):
            vrt.append(
                '  <VRTRasterBand dataType="%s" band="%s" subClass="VRTRawRasterBand">'
                % (datatype, i + 1)
            )
            vrt.append(
                '    <SourceFilename relativeToVRT="%s">%s</SourceFilename>'
                % (relativeToVRT, filename)
            )
            if nodata is not None:
                vrt.append(
                    "    <NoDataValue>%s</NoDataValue>" % (nodata)
                )  # Fix for Issue 17
            vrt.append(
                "    <ImageOffset>%s</ImageOffset>"
                % (headeroffset + nbits / 8 * i * cols)
            )
            vrt.append("    <PixelOffset>%s</PixelOffset>" % (nbits / 8))
            vrt.append("    <LineOffset>%s</LineOffset>" % (nbits / 8 * cols))
            if byteorder:
                vrt.append("    <ByteOrder>%s</ByteOrder>" % (byteorder))
            vrt.append("  </VRTRasterBand>")
        return CreateCustomVRT("\n".join(vrt), cols, rows)
    except:
        return None
    return "\n".join(vrt)


def CreateBIPRawRasterVRT(
    filename,
    nbands,
    cols,
    rows,
    datatype,
    nbits,
    nodata=None,
    headeroffset=0,
    byteorder=None,
    relativeToVRT=0,
):
    """Create RawRaster VRT from BIP

    BIP = Band-Interleaved-by-Pixel or Pixel-Interleaved

    For further info on VRT's, see the U{GDAL VRT Tutorial<http://www.gdal.org/gdal_vrttut.html>}

    @type filename:       C{str}
    @param filename:      File name
    @type nbands:         C{int}
    @param nbands:        The number of bands in the output VRT (<= nbands in input file)
    @type cols:           C{int}
    @param cols:          The number of columns in the output VRT
    @type rows:           C{int}
    @param rows:          The number of rows in the output VRT
    @type datatype:       C{str}
    @param datatype:      GDAL datatype name. Eg. Byte, Int32, UInt16
    @type nodata:         C{float}
    @param nodata:        No data/Null value
    @type headeroffset:   C{int}
    @param headeroffset:  Number of bytes to skip at the start of the file
    @type byteorder:      C{str}
    @param byteorder:     Byte order of the file (MSB or LSB)

    @rtype:               C{xml}
    @return:              VRT XML string
    """
    try:
        vrt = []
        nbits = gdal.GetDataTypeSize(gdal.GetDataTypeByName(datatype))
        for i in range(nbands):
            vrt.append(
                '  <VRTRasterBand dataType="%s" band="%s" subClass="VRTRawRasterBand">'
                % (datatype, i + 1)
            )
            vrt.append(
                '    <SourceFilename relativeToVRT="%s">%s</SourceFilename>'
                % (relativeToVRT, filename)
            )
            if nodata is not None:
                vrt.append(
                    "    <NoDataValue>%s</NoDataValue>" % (nodata)
                )  # Fix for Issue 17
            vrt.append(
                "    <ImageOffset>%s</ImageOffset>" % (headeroffset + nbits / 8 * i)
            )
            vrt.append("    <PixelOffset>%s</PixelOffset>" % (nbits / 8))
            vrt.append("    <LineOffset>%s</LineOffset>" % (nbits / 8 * cols))
            if byteorder:
                vrt.append("    <ByteOrder>%s</ByteOrder>" % (byteorder))
            vrt.append("  </VRTRasterBand>")
        return CreateCustomVRT("\n".join(vrt), cols, rows)
    except:
        return None
    return "\n".join(vrt)


def CreateCustomVRT(vrtxml, vrtcols, vrtrows):
    try:
        vrt = []
        vrt.append(
            '<VRTDataset rasterXSize="%s" rasterYSize="%s">' % (vrtcols, vrtrows)
        )
        vrt.append("%s" % vrtxml)
        vrt.append("</VRTDataset>")
        return "\n".join(vrt)
    except:
        return None


# ========================================================================================================
# {Shapefile Writer
# ========================================================================================================
class ShapeWriter:
    """A class for writing geometry and fields to ESRI shapefile format"""

    def __init__(self, shapefile, fields={}, srs_wkt=None, update=False):
        """Open the shapefile for writing or appending.

        @type shapefile:  C{gdal.Dataset}
        @param shapefile: Dataset object
        @type fields:     C{dict}
        @param fields:    L{Fields dict<formats.fields>}
        @type srs_wkt:    C{str}
        @param srs_wkt:   Spatial reference system WKT
        @type update:     C{boolean}
        @param update:    Update or overwrite existing shapefile

        @note: Field names can only be <= 10 characters long. Longer names will be silently truncated. This may result in non-unique column names, which will definitely cause problems later.
               Field names can not contain spaces or special characters, except underscores.
               Starting with version 1.7, the OGR Shapefile driver tries to generate unique field names. Successive duplicate field names, including those created by truncation to 10 characters, will be truncated to 8 characters and appended with a serial number from 1 to 99.

        @see: U{http://www.gdal.org/ogr/drv_shapefile.html}
        """
        gdal.ErrorReset()
        ogr.UseExceptions()
        self._driver = ogr.GetDriverByName("ESRI Shapefile")
        self._srs = osr.SpatialReference()
        self._filename = shapefile
        self._fields = fields
        self._srs_wkt = srs_wkt
        self.fields = {}  # Dict to map full names to truncated names if name >10 chars
        try:
            if update and os.path.exists(shapefile):
                self._shape = self.__openshapefile__()
            else:
                self._shape = self.__createshapefile__()
        except Exception as err:
            self.__error__(err)
        ogr.DontUseExceptions()

    def __del__(self):
        """Shutdown and release the lock on the shapefile"""
        try:
            gdal.ErrorReset()
            self._shape.Release()
        except:
            pass

    def __error__(self, err):
        gdalerr = gdal.GetLastErrorMsg()
        gdal.ErrorReset()
        self.__del__()
        errmsg = str(err)
        if gdalerr:
            errmsg += "\n%s" % gdalerr
        raise err.__class__(errmsg)

    def __createshapefile__(self):
        """Open the shapefile for writing"""
        if self._srs_wkt:
            self._srs.ImportFromWkt(self._srs_wkt)
        else:
            self._srs.ImportFromEPSG(4283)  # default=GDA94 Geographic

        if os.path.exists(self._filename):
            self._driver.DeleteDataSource(self._filename)
        shp = self._driver.CreateDataSource(self._filename)
        lyr = os.path.splitext(os.path.split(self._filename)[1])[0]
        lyr = shp.CreateLayer(lyr, geom_type=ogr.wkbPolygon, srs=self._srs)
        i = 0
        fieldnames = sorted(self._fields.keys())
        for f in fieldnames:
            if self._fields[f]:
                # Get field types
                if type(self._fields[f]) in [list, tuple]:
                    ftype = self._fields[f][0]
                    fwidth = self._fields[f][1]
                else:
                    ftype = self._fields[f]
                    fwidth = 0
                if ftype.upper() == "STRING":
                    fld = ogr.FieldDefn(f, ogr.OFTString)
                    fld.SetWidth(fwidth)
                elif ftype.upper() == "INT":
                    fld = ogr.FieldDefn(f, ogr.OFTInteger)
                elif ftype.upper() == "FLOAT":
                    fld = ogr.FieldDefn(f, ogr.OFTReal)
                else:
                    continue
                # raise AttributeError, 'Invalid field definition'
                lyr.CreateField(fld)
                if len(f) > 10:
                    self.fields[f] = lyr.GetLayerDefn().GetFieldDefn(i).GetName()
                else:
                    self.fields[f] = f
                i += 1

        return shp

    def __openshapefile__(
        self,
    ):
        """Open the shapefile for updating/appending"""
        fieldnames = sorted(self._fields.keys())
        shp = self._driver.Open(self._filename, update=1)
        if not shp:
            raise GDALError("Unable to open %s" % self._filename)
        lyr = shp.GetLayer(0)
        lyrdef = lyr.GetLayerDefn()
        if lyrdef.GetFieldCount() == 0:
            del lyrdef
            del lyr
            del shp
            return self.__createshapefile__()
        self._srs = lyr.GetSpatialRef()
        # Loop thru input fields and match to shp fields
        i = 0
        for f in fieldnames:
            if self._fields[f]:  # does it have a field type?
                if len(f) > 10:
                    self.fields[f] = lyrdef.GetFieldDefn(i).GetName()
                else:
                    self.fields[f] = f
                i += 1
        if i == 0:  # No fieldnames...
            for i in range(lyrdef.GetFieldCount()):
                f = lyrdef.GetFieldDefn(i).GetName()
                self.fields[f] = f
        return shp

    def WriteRecord(self, extent, attributes):
        """Write record

        @type extent:      C{list}
        @param extent:     [[ulx,uly],[urx,ury],[lrx,lry],[llx,lly]]
        @type attributes:  C{dict}
        @param attributes: Must match field names passed to __init__()
        """
        try:
            geom = GeomFromExtent(extent, self._srs)
            if (
                self._srs.IsGeographic()
            ):  # basic coordinate bounds test. Can't do for projected though
                srs = osr.SpatialReference()
                srs.ImportFromEPSG(4283)  # 4326)
                valid = GeomFromExtent([-180, -90, 180, 90], srs=srs)
                if not valid.Contains(geom):
                    # raise ValueError, 'Invalid extent coordinates'
                    warnings.warn("Invalid extent coordinates")

            lyr = self._shape.GetLayer(0)

            feat = ogr.Feature(lyr.GetLayerDefn())
            for a in attributes:
                # if a in self.fields:feat.SetField(a, attributes[a])
                if a in self.fields:
                    feat.SetField(self.fields[a], encode(attributes[a]))
            feat.SetGeometryDirectly(geom)
            lyr.CreateFeature(feat)
            lyr.SyncToDisk()

        except Exception as err:
            self.__error__(err)

    def UpdateRecord(self, extent, attributes, where_clause):
        """Update record/s

        @type where_clause:  C{str}
        @param where_clause: Shapefile supported SQL where clause
        @type attributes:    C{dict}
        @param attributes:   Must match field names passed to __init__()
        """
        try:
            geom = GeomFromExtent(extent, self._srs)
            if (
                self._srs.IsGeographic()
            ):  # basic coordinate bounds test. Can't do for projected though
                srs = osr.SpatialReference()
                srs.ImportFromEPSG(4283)  # 4326)
                valid = GeomFromExtent([-180, -90, 180, 90], srs=srs)
                if not valid.Contains(geom):
                    # raise ValueError, 'Invalid extent coordinates'
                    warnings.warn("Invalid extent coordinates")
            lyr = self._shape.GetLayer()
            lyr.SetAttributeFilter(where_clause)
            feat = lyr.GetNextFeature()
            try:
                feat.SetGeometryDirectly(geom)
            except:
                return self.WriteRecord(extent, attributes)
            while feat:
                for a in attributes:
                    if a in self.fields:
                        feat.SetField(self.fields[a], encode(attributes[a]))
                lyr.SetFeature(feat)
                feat = lyr.GetNextFeature()
            lyr.SyncToDisk()
        except Exception as err:
            self.__error__(err)

    def DeleteRecord(self, where_clause):
        """Delete record/s

        @type where_clause:  C{str}
        @param where_clause: Shapefile supported SQL where clause
        """
        try:
            lyr = self._shape.GetLayer()
            lyr.SetAttributeFilter(where_clause)
            feat = lyr.GetNextFeature()
            while feat:
                fid = feat.GetFID()
                lyr.DeleteFeature(fid)
                feat = lyr.GetNextFeature()

            lyr.SyncToDisk()
        except Exception as err:
            self.__error__(err)


# }
