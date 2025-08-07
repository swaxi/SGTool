# This file is your entry point:
# - add you Python files and folder inside this 'flows' folder
# - add your imports
# - just don't change the name of the function 'run()' nor this filename ('sgtool.py')
#   and everything is gonna be ok.
#
# Remember: everything is gonna be ok in the end: if it's not ok, it's not the end.
# Alternatively, ask for help at https://github.com/deeplime-io/onecode/issues

import onecode
from onecode import Logger, text_input, file_input

import rasterio
import numpy as np
from PIL import Image


def geotiff_to_png(input_path, output_path, band=1):
    """
    Convert a single band GeoTIFF to PNG

    Parameters:
    input_path (str): Path to input GeoTIFF file
    output_path (str): Path for output PNG file
    band (int): Band number to read (default: 1)
    """

    # Open the GeoTIFF file
    with rasterio.open(input_path) as src:
        # Read the specified band
        band_data = src.read(band)

        # Normalize data to 0-255 range for PNG
        # Handle different data types and potential nodata values
        if src.nodata is not None:
            # Mask nodata values
            band_data = np.ma.masked_equal(band_data, src.nodata)

        # Normalize to 0-255 range
        data_min = band_data.min()
        data_max = band_data.max()

        if data_max > data_min:
            normalized = ((band_data - data_min) / (data_max - data_min) * 255).astype(
                np.uint8
            )
        else:
            normalized = np.zeros_like(band_data, dtype=np.uint8)

        # Handle masked arrays (convert masked values to 0)
        if hasattr(normalized, "filled"):
            normalized = normalized.filled(0)

    # Convert to PIL Image and save as PNG
    image = Image.fromarray(normalized, mode="L")  # 'L' for grayscale
    image.save(output_path)

    print(f"Successfully converted {input_path} to {output_path}")


def run():
    onecode.Logger.info(
        """
        #####################################################################
        ###> Hello from SGTool!
        ###> Fill in this run() function with something awesome!
        #####################################################################
        """
    )
    input_path = file_input(
        key="my_file",
        value="uploads/my_image.tif",
        label="Select an grid",
        types=["GeoTIFF files (*.tif;*.tiff)"],
    )
    geotiff_to_png(input_path, "./uploads", band=1)
