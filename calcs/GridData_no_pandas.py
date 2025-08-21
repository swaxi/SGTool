import processing


class QGISGridData:
    def __init__(self, iface):
        self.iface = iface

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

        processing.execAlgorithmDialog(alg_id, pre_filled_params)

    def launch_multi_bspline_dialog(self, input, zcolumn, cell_size, mask):

        # Set up the parameters you want pre-filled
        pre_filled_params = {
            "SHAPES": input,  # Reference to your input layer
            "FIELD": zcolumn,  # Z-value field
            "TARGET_USER_SIZE": cell_size,  # cell size
        }

        alg_id = "sagang:multilevelbspline"
        processing.execAlgorithmDialog(alg_id, pre_filled_params)
