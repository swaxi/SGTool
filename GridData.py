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
        self.original_df = df.copy()  # Keep a copy of the original data
        self.x = df["x"].values
        self.y = df["y"].values
        self.values = df["value"].values
        # self.grid_bounds = [-1, 1, -1, 1]
        self.nx = nx
        self.ny = ny

        # Normalize data
        self.normalize_data()

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
            grid_data = self.akima_2d_interpolation()
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

    def akima_2d_interpolation(self):
        """
        Perform Akima interpolation in 2D using 1D interpolation along each axis.

        This function first grids the scattered data (if necessary) and then performs
        1D Akima interpolation along each axis to fill the grid.
        """

        grid_x = self.grid_x
        grid_y = self.grid_y

        # Grid the scattered data using scipy.interpolate.griddata
        grid_values = griddata(
            points=(self.x, self.y),
            values=self.values,  # Scattered points and values
            xi=(grid_x, grid_y),
            method="linear",  # Grid points and interpolation method
        )

        # Handle NaN values in the grid (if any)
        if np.isnan(grid_values).any():
            grid_values = np.nan_to_num(grid_values, nan=np.mean(self.values))

        # Interpolate along the y-axis first (rows)
        intermediate = np.zeros((grid_y.shape[0], grid_x.shape[1]))
        for i in range(grid_x.shape[1]):
            akima = Akima1DInterpolator(
                grid_y[:, 0], grid_values[:, i], method="makima"
            )
            intermediate[:, i] = akima(grid_y[:, 0])

        # Interpolate along the x-axis (columns)
        result = np.zeros((grid_y.shape[0], grid_x.shape[1]))
        for i in range(grid_y.shape[0]):
            akima = Akima1DInterpolator(
                grid_x[0, :], intermediate[i, :], method="makima"
            )
            result[i, :] = akima(grid_x[0, :])

        return result

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
