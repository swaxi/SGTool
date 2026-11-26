# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2024-2026 Mark Jessell
"""
The pysgtool package: SGTool without QGIS or ArcGIS.

    pip install .            array calculations (numpy, scipy)
    pip install .[io]        + reading and writing files, and the IGRF for a grid
    pip install .[all]       + PCA / ICA, worm lines, plots and WTMM

pysgtool/ is the API; the calculations are the plugin's own calcs/ folder, installed
inside it as pysgtool.calcs (not copied or moved in the repository). The version is
the plugin's, from metadata.txt.
"""

import atexit
import os
import re
import shutil
import tempfile

from setuptools import setup

HERE = os.path.dirname(os.path.abspath(__file__))


def text(name):
    with open(os.path.join(HERE, name), encoding="utf-8") as f:
        return f.read()


# built in a temporary folder: a file sync (Dropbox, OneDrive) locks the files of a
# build/ folder here, and an old one would leave stale files in the wheel
BUILD = tempfile.mkdtemp(prefix="pysgtool_build_")
atexit.register(shutil.rmtree, BUILD, ignore_errors=True)

SUBPACKAGES = ["euler", "igrf", "worms"]

EXTRAS = {
    # GDAL (osgeo) is used instead of rasterio and pyogrio where it is installed
    "io": ["rasterio", "pyogrio", "pyproj", "shapely>=2"],
    "ml": ["scikit-learn", "shapely>=2"],
    "plot": ["matplotlib"],
    "wtmm": ["PyWavelets", "matplotlib"],
    "test": ["pytest"],
}
EXTRAS["all"] = sorted({r for k, v in EXTRAS.items() if k != "test" for r in v})

setup(
    name="pysgtool",
    version=re.search(r"^version=(.+)$", text("metadata.txt"), re.M).group(1).strip(),
    description="Potential field and radiometric grid processing: the SGTool calculations "
                "without QGIS or ArcGIS",
    long_description=text("docs/pysgtool.md"),
    long_description_content_type="text/markdown",
    author="Mark Jessell",
    url="https://github.com/swaxi/SGTool",
    packages=["pysgtool", "pysgtool.calcs"] + [f"pysgtool.calcs.{p}" for p in SUBPACKAGES],
    package_dir=dict(
        {"pysgtool": "pysgtool", "pysgtool.calcs": "calcs"},
        **{f"pysgtool.calcs.{p}": f"calcs/{p}" for p in SUBPACKAGES},
    ),
    package_data={
        "pysgtool.calcs.igrf": ["SHC_files/*.SHC", "*.txt", "*.md"],
        "pysgtool.calcs.euler": ["*.md"],
        "pysgtool.calcs.worms": ["*.txt"],
    },
    license="GPL-2.0-or-later",
    license_files=["LICENSE"],
    python_requires=">=3.9",
    install_requires=["numpy", "scipy"],
    extras_require=EXTRAS,
    options={"build": {"build_base": BUILD}},
)
