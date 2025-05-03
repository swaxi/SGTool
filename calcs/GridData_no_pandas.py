import numpy as np
from scipy.interpolate import (
    griddata,
    Akima1DInterpolator,
    CloughTocher2DInterpolator,
    Rbf,
)
from scipy.spatial import ConvexHull, cKDTree

from qgis.core import (
    QgsVectorLayer,
    QgsField,
    QgsFeature,
    QgsPointXY,
    QgsProcessingFeedback,
    QgsRasterLayer,
    QgsGeometry,
    QgsApplication,
)
from qgis.PyQt.QtWidgets import QMessageBox

from qgis.PyQt.QtCore import QVariant
import processing
import numpy as np
import csv


class QGISGridData:
    def __init__(self, iface):
        self.iface = iface

    def launch_r_surf_rst_dialog(self, input, zcolumn, cell_size, mask):
        """
        Launch the v.surf.rst.cvdev dialog from the Processing Toolbox.
        """
        pre_filled_params = {
            "input": input,
            "zcolumn": zcolumn,
            "GRASS_REGION_CELLSIZE_PARAMETER": cell_size,  # cell size from sgtoosl dialog
            # "mask": mask,
        }
        alg_id = "grass7:v.surf.rst.cvdev"
        try:
            # Check if the algorithm exists
            if QgsApplication.processingRegistry().algorithmById(alg_id):
                # Launch the dialog
                processing.execAlgorithmDialog(alg_id, pre_filled_params)
            else:
                self.iface.messageBar().pushMessage(
                    "Error", "GRASS v.surf.rst.cvdev algorithm not found.", level=3
                )
        except Exception as e:
            self.iface.messageBar().pushMessage("Error", str(e), level=3)

    def launch_idw_dialog(self, input, zcolumn, cell_size, mask):
        """
        Launch the v.surf.idw dialog from the Processing Toolbox.
        """
        pre_filled_params = {
            "input": input,
            "column": zcolumn,
            "GRASS_REGION_CELLSIZE_PARAMETER": cell_size,  # cell size from sgtoosl dialog
        }
        alg_id = "grass7:v.surf.idw"
        try:
            # Check if the algorithm exists
            if QgsApplication.processingRegistry().algorithmById(alg_id):
                # Launch the dialog
                processing.execAlgorithmDialog(alg_id, pre_filled_params)
            else:
                self.iface.messageBar().pushMessage(
                    "Error", "GRASS v.surf.rst algorithm not found.", level=3
                )
        except Exception as e:
            self.iface.messageBar().pushMessage("Error", str(e), level=3)

    def list_grass_algorithms(self):
        for alg in QgsApplication.processingRegistry().algorithms():
            if "grass" in alg.id().lower():
                print(alg.id())

    def filter_points(self, layer, decimation_factor):
        filtered_layer = processing.run(
            "native:extractbyexpression",
            {
                "INPUT": layer,
                "EXPRESSION": f"$id % {decimation_factor} = 0",
                "OUTPUT": "memory:",
            },
        )["OUTPUT"]
        return filtered_layer

    def launch_multi_bspline_dialog(self, input, zcolumn, cell_size, mask):

        # Set up the parameters you want pre-filled
        pre_filled_params = {
            "SHAPES": input,  # Reference to your input layer
            "FIELD": zcolumn,  # Z-value field
            "TARGET_USER_SIZE": cell_size,  # cell size
        }

        alg_id = "sagang:multilevelbspline"
        processing.execAlgorithmDialog(alg_id, pre_filled_params)
