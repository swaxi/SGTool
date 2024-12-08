import numpy as np
import pandas as pd
from scipy.interpolate import (
    griddata,
    Akima1DInterpolator,
    CloughTocher2DInterpolator,
    Rbf,
)
from scipy.spatial import cKDTree
from scipy.optimize import minimize


class GridData:
    def __init__(self, df, nx, ny, grid_bounds=None):
        self.validate_input(df)
        df = df.dropna(subset=["value"])
        df = df.apply(pd.to_numeric, errors="coerce")
        self.original_df = df.copy()  # Keep a copy of the original data
        self.x = df["x"].values
        self.y = df["y"].values
        self.values = df["value"].values
        # self.grid_bounds = [-1, 1, -1, 1]
        self.nx = nx
        self.ny = ny

        # Normalize data
        self.normalize_data()
        self.interpolator = None

        self.grid_x, self.grid_y = self.create_grid(grid_bounds=None)

    @staticmethod
    def validate_input(df):
        if not all(col in df.columns for col in ["x", "y", "value"]):
            raise ValueError("DataFrame must contain 'x', 'y', and 'value' columns")
        if df.empty:
            raise ValueError("DataFrame is empty. Provide valid input data.")

    def normalize_data(self):
        # Compute normalization factors
        self.x_mean, self.x_std = self.x.mean(), self.x.std()
        self.y_mean, self.y_std = self.y.mean(), self.y.std()
        self.value_mean, self.value_std = self.values.mean(), self.values.std()

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
            x_min, x_max = self.x.min(), self.x.max()
            y_min, y_max = self.y.min(), self.y.max()

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
        tree = cKDTree(np.c_[self.x, self.y])
        xi = np.c_[self.grid_x.ravel(), self.grid_y.ravel()]
        distances, indices = tree.query(xi, k=10)
        weights = 1 / (distances**power + 1e-12)
        interpolated_values = np.sum(weights * self.values[indices], axis=1) / np.sum(
            weights, axis=1
        )
        return interpolated_values.reshape(self.grid_x.shape)

    def save_to_csv(self, filename, grid_data):
        df = pd.DataFrame(
            {
                "x": self.grid_x.ravel() * self.x_std + self.x_mean,
                "y": self.grid_y.ravel() * self.y_std + self.y_mean,
                "value": grid_data.ravel(),
            }
        )
        df.to_csv(filename, index=False)
        print(f"Grid data saved to {filename}")

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
