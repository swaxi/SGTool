import numpy as np
from scipy.interpolate import (
    griddata,
    Akima1DInterpolator,
    CloughTocher2DInterpolator,
    Rbf,
)
from scipy.spatial import ConvexHull, cKDTree
from matplotlib.path import Path

from qgis.core import (
    QgsVectorLayer,
    QgsField,
    QgsFeature,
    QgsPointXY,
    QgsProcessingFeedback,
    QgsRasterLayer,
    QgsGeometry,
    QgsApplication,
    QgsApplication,

)

from qgis.PyQt.QtCore import QVariant
import processing
import numpy as np
import csv

class GridData:
    def __init__(self, data, nx, ny, grid_bounds=None, normalize=False):
        self.validate_input(data)
        data = self.clean_data(data)
        self.x, self.y, self.values = data[:, 0], data[:, 1], data[:, 2]
        self.nx = nx
        self.ny = ny

        # Normalize data
        # if normalize:
        self.normalize_data()
        self.interpolator = None

        self.grid_x, self.grid_y = self.create_grid(grid_bounds=grid_bounds)

    @staticmethod
    def validate_input(data):
        if not isinstance(data, np.ndarray) or data.shape[1] != 3:
            raise ValueError(
                "Input data must be a numpy array with three columns: x, y, and value."
            )

    @staticmethod
    def clean_data(data):
        # Remove rows with NaN values
        return data[~np.isnan(data).any(axis=1)]

    def normalize_data(self):
        # Compute normalization factors
        self.x_mean, self.x_std = np.mean(self.x), np.std(self.x)
        self.y_mean, self.y_std = np.mean(self.y), np.std(self.y)
        self.value_mean, self.value_std = np.mean(self.values), np.std(self.values)

        # Normalize x, y, and value
        self.x = (self.x - self.x_mean) / self.x_std
        self.y = (self.y - self.y_mean) / self.y_std
        self.values = (self.values - self.value_mean) / self.value_std

    def denormalize_grid(self, grid_data):
        # Denormalize the interpolated grid values
        return grid_data * self.value_std + self.value_mean

    def create_grid(self, grid_bounds=None):
        if grid_bounds:
            x_min, x_max, y_min, y_max = grid_bounds
        else:
            x_min, x_max = np.min(self.x), np.max(self.x)
            y_min, y_max = np.min(self.y), np.max(self.y)

        grid_x, grid_y = np.meshgrid(
            np.linspace(x_min, x_max, self.nx),
            np.linspace(y_min, y_max, self.ny),
        )

        return grid_x, grid_y

    def interpolate(self, method="akima", **kwargs):
        if method == "akima":
            grid_data = self.akima_interpolation()
        elif method == "clough_tocher":
            grid_data = self.clough_tocher_interpolation()
        elif method == "rbf":
            grid_data = self.rbf_interpolation(kwargs.get("function", "multiquadric"))
        elif method == "minimum_curvature":
            grid_data = self.minimum_curvature_interpolation()
        elif method == "idw":
            grid_data = self.idw_interpolation(kwargs.get("power", 2))
        else:
            raise ValueError(
                "Invalid method. Choose from: 'akima', 'clough_tocher', 'rbf', 'minimum_curvature', 'idw'"
            )

        # Denormalize the interpolated values
        return self.denormalize_grid(grid_data)

    def clough_tocher_interpolation(self):
        interpolator = CloughTocher2DInterpolator(
            list(zip(self.x, self.y)), self.values
        )
        return interpolator(self.grid_x, self.grid_y)

    def rbf_interpolation(self, function="multiquadric"):
        interpolator = Rbf(self.x, self.y, self.values, function=function)
        return interpolator(self.grid_x, self.grid_y)

    def minimum_curvature_interpolation(self):
        interpolated_values = griddata(
            points=(self.x, self.y),  # Scattered points
            values=self.values,  # Values at scattered points
            xi=(self.grid_x, self.grid_y),  # Target grid
            method="cubic",
        )
        return interpolated_values

    def idw_interpolation(self, power=2):
        # Build k-D tree for the input points
        tree = cKDTree(np.c_[self.x, self.y])
        xi = np.c_[self.grid_x.ravel(), self.grid_y.ravel()]

        # Compute distances and weights for IDW
        distances, indices = tree.query(xi, k=10)
        weights = 1 / (distances**power + 1e-12)
        interpolated_values = np.sum(weights * self.values[indices], axis=1) / np.sum(
            weights, axis=1
        )

        # Compute the convex hull of the input data points
        hull = ConvexHull(np.c_[self.x, self.y])
        hull_path = Path(np.c_[self.x[hull.vertices], self.y[hull.vertices]])

        # Check which grid points are inside the convex hull
        inside_hull = hull_path.contains_points(xi).reshape(self.grid_x.shape)

        # Create the interpolated grid and set NaN for points outside the hull
        interpolated_grid = interpolated_values.reshape(self.grid_x.shape)
        interpolated_grid[~inside_hull] = np.nan

        return interpolated_grid

    def save_to_csv(self, filename, grid_data):
        with open(filename, mode="w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(["x", "y", "value"])
            for i in range(self.grid_x.shape[0]):
                for j in range(self.grid_x.shape[1]):
                    writer.writerow(
                        [
                            self.grid_x[i, j] * self.x_std + self.x_mean,
                            self.grid_y[i, j] * self.y_std + self.y_mean,
                            grid_data[i, j],
                        ]
                    )

    def akima_interpolation(self):
        """
        Perform 2D Akima interpolation on the grid.

        Returns:
        - Interpolated values on the grid.
        """
        akima_2d = IterativeAkima2D(
            self.x, self.y, self.values, self.grid_x, self.grid_y, iterations=5
        )
        grid_z = akima_2d.interpolate()

        return grid_z

    def interpolate_with_v_surf_rst(self, output_raster, epsg, cell_size=10):
        """
        Interpolate a 2D grid using GRASS GIS v.surf.rst based on input x, y, and value arrays.

        Parameters:
            x (array-like): Array of x-coordinates.
            y (array-like): Array of y-coordinates.
            values (array-like): Array of values at the (x, y) points.
            cell_size (float): The resolution of the output grid.

        Returns:
            np.ndarray: A 2D numpy array of the interpolated grid.
        """
        print(
            "region",
            f"{min(self.x)},{max(self.x)},{min(self.y)},{max(self.y)},{cell_size},{cell_size}",
            epsg,
            cell_size,
        )
        # Step 1: Create a temporary vector layer from input arrays
        mem_layer = QgsVectorLayer(
            "Point?crs=EPSG:{}".format(epsg), "temp_points", "memory"
        )
        provider = mem_layer.dataProvider()

        # Add fields
        provider.addAttributes([QgsField("value", QVariant.Double)])
        mem_layer.updateFields()

        # Add features
        features = []
        for xi, yi, vi in zip(self.x, self.y, self.values):
            feature = QgsFeature()
            point = QgsPointXY(xi, yi)  # Create a point
            feature.setGeometry(
                QgsGeometry.fromPointXY(point)
            )  # Wrap the point in QgsGeometry
            feature.setAttributes([vi])
            features.append(feature)
        provider.addFeatures(features)
        mem_layer.updateExtents()

        # Step 2: Run GRASS v.surf.rst tool
        # output_raster = filepath
        params = {
            "input": mem_layer,
            "zcolumn": "value",
            "elevation": output_raster,
            "GRASS_REGION_PARAMETER": f"{min(self.x)},{max(self.x)},{min(self.y)},{max(self.y)} [EPSG:{epsg}]",
            "GRASS_REGION_CELLSIZE_PARAMETER": cell_size,
        }

        feedback = QgsProcessingFeedback()
        processing.run("grass7:v.surf.rst", params, feedback=feedback)

        # Step 3: Read the raster as a numpy array
        raster_layer = QgsRasterLayer(output_raster, "interpolated_raster", "gdal")
        if not raster_layer.isValid():
            raise ValueError("Failed to load interpolated raster.")

        # Extract the raster data as a numpy array
        provider = raster_layer.dataProvider()
        block = provider.block(
            1, raster_layer.extent(), raster_layer.width(), raster_layer.height()
        )
        numpy_grid = np.array(block.data(), dtype=float).reshape(
            raster_layer.height(), raster_layer.width()
        )

        return numpy_grid


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

    def launch_bspline_dialog(self, input, zcolumn, cell_size, mask):
        """
        Launch the v.surf.bspline dialog from the Processing Toolbox.
        """
        pre_filled_params = {
            "input": input,
            "column": zcolumn,
            "GRASS_REGION_CELLSIZE_PARAMETER": cell_size,  # cell size from sgtoosl dialog
        }
        alg_id = "grass7:v.surf.bspline"
        try:
            # Check if the algorithm exists
            if QgsApplication.processingRegistry().algorithmById(alg_id):
                # Launch the dialog
                processing.execAlgorithmDialog(alg_id, pre_filled_params)
            else:
                self.iface.messageBar().pushMessage(
                    "Error", "GRASS v.surf.bspline algorithm not found.", level=3
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
            'SHAPES': input,  # Reference to your input layer
            'FIELD': zcolumn,        # Z-value field
            'TARGET_USER_SIZE': cell_size      # cell size
        }

        alg_id = "sagang:multilevelbspline"
        try:
            # Check if the algorithm exists
            if QgsApplication.processingRegistry().algorithmById(alg_id):
                # Launch the dialog
                processing.execAlgorithmDialog(alg_id, pre_filled_params)
            else:
                self.iface.messageBar().pushMessage(
                    "Error", "sagang multilevelbspline algorithm not found.\nTry installing the Plugin: Saga Processing Saga NextGen Provider", level=3
                )
        except Exception as e:
            self.iface.messageBar().pushMessage("Error", str(e), level=3)  
    



class IterativeAkima2D:
    def __init__(self, x, y, values, grid_x, grid_y, iterations=3):
        """
        Iterative 2D Akima Interpolation.

        Parameters:
        - x: 1D array of x-coordinates of data points.
        - y: 1D array of y-coordinates of data points.
        - values: 1D array of values at the data points.
        - grid_x: 2D array of x-coordinates for the output grid.
        - grid_y: 2D array of y-coordinates for the output grid.
        - iterations: Number of iterations to refine the interpolation.
        """
        self.x = x
        self.y = y
        self.values = values
        self.grid_x = grid_x
        self.grid_y = grid_y
        self.iterations = iterations
        self.grid_z = np.full(grid_x.shape, np.nan)

    def interpolate(self):
        """
        Perform iterative 2D Akima interpolation.

        Returns:
        - 2D array of interpolated values.
        """
        # Initialize grid with nearest-neighbor interpolation for better propagation
        self.grid_z = self.nearest_neighbor_init()

        for iteration in range(self.iterations):
            # Copy grid for blending
            row_result = np.copy(self.grid_z)
            col_result = np.copy(self.grid_z)

            # Interpolate row-wise
            for i, yi in enumerate(self.grid_y[:, 0]):
                row_indices = np.where(np.isclose(self.y, yi, atol=1e-6))[0]
                if len(row_indices) < 2:
                    continue

                x_row = self.x[row_indices]
                z_row = self.values[row_indices]

                # Ensure strict ordering and remove duplicates
                sorted_indices = np.argsort(x_row)
                x_row = x_row[sorted_indices]
                z_row = z_row[sorted_indices]
                x_row, unique_indices = np.unique(x_row, return_index=True)
                z_row = z_row[unique_indices]

                if len(x_row) < 2:
                    continue

                akima = Akima1DInterpolator(x_row, z_row)
                row_result[i, :] = akima(self.grid_x[i, :])

            # Interpolate column-wise
            for j, xi in enumerate(self.grid_x[0, :]):
                col_indices = np.where(np.isclose(self.x, xi, atol=1e-6))[0]
                if len(col_indices) < 2:
                    continue

                y_col = self.y[col_indices]
                z_col = self.values[col_indices]

                # Ensure strict ordering and remove duplicates
                sorted_indices = np.argsort(y_col)
                y_col = y_col[sorted_indices]
                z_col = z_col[sorted_indices]
                y_col, unique_indices = np.unique(y_col, return_index=True)
                z_col = z_col[unique_indices]

                if len(y_col) < 2:
                    continue

                akima = Akima1DInterpolator(y_col, z_col)
                col_result[:, j] = akima(self.grid_y[:, j])

            # Blend row-wise and column-wise results
            valid_row = np.isfinite(row_result)
            valid_col = np.isfinite(col_result)
            self.grid_z[valid_row & valid_col] = (
                0.5 * row_result[valid_row & valid_col]
                + 0.5 * col_result[valid_row & valid_col]
            )

            # Fill gaps in the grid from either row-wise or column-wise results
            self.grid_z[valid_row & ~valid_col] = row_result[valid_row & ~valid_col]
            self.grid_z[~valid_row & valid_col] = col_result[~valid_row & valid_col]

        return self.grid_z

    def nearest_neighbor_init(self):
        """
        Initialize the grid using nearest-neighbor interpolation.

        Returns:
        - 2D array of initial grid values.
        """
        from scipy.interpolate import griddata

        return griddata(
            points=np.column_stack((self.x, self.y)),
            values=self.values,
            xi=(self.grid_x, self.grid_y),
            method="nearest",
            fill_value=np.nan,
        )
