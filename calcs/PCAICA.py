import numpy as np
from osgeo import gdal
from qgis.core import QgsRasterLayer
import os


class PCAICA:
    def __init__(self, grid):
        """
        Initialize the class SpatialStats filter:

        :param grid: 2D numpy array representing the input grid
        """
        self.grid = np.array(grid, dtype=float)

    def pca_with_nans(self, input_raster_path, output_raster_path, n_components=None):
        """
        Perform PCA on a raster with NaN values and save as a multiband GeoTIFF

        Parameters:
        input_raster_path (str): Path to the input GeoTIFF
        output_raster_path (str): Path to save the output multiband GeoTIFF with PCA components
        n_components (int): Number of principal components to calculate (defaults to all possible)

        Returns:
        tuple: (components, explained_variance_ratio)
        """
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler

        # Load raster using PyQGIS to access metadata
        raster_layer = QgsRasterLayer(input_raster_path, "Input Raster")
        if not raster_layer.isValid():
            raise ValueError(f"Failed to load raster layer: {input_raster_path}")

        # Also open with GDAL for data access
        ds = gdal.Open(input_raster_path)
        if ds is None:
            raise ValueError(f"Failed to open raster with GDAL: {input_raster_path}")

        # Get raster dimensions
        width = ds.RasterXSize
        height = ds.RasterYSize
        bands = ds.RasterCount
        if bands == 1:
            print("Only one band found, PCA not possible")
            return None, None

        # If n_components is not specified, use all possible components
        if n_components is None or n_components == 0:
            n_components = bands
        elif n_components > bands:
            n_components = bands
            print(
                f"Warning: Requested {n_components} components but only {bands} bands available. Using {bands} components."
            )

        # Create a 3D numpy array to hold all band data
        data_array = np.zeros((bands, height, width))

        # Read all bands into the array
        for b in range(bands):
            band = ds.GetRasterBand(b + 1)
            data_array[b] = band.ReadAsArray()

        # Reshape for PCA (bands as features, pixels as samples)
        # Transpose from (bands, height, width) to (height*width, bands)
        X = data_array.reshape(bands, -1).T

        # Find valid pixels (no NaN in any band)
        valid_mask = ~np.isnan(X).any(axis=1)
        X_valid = X[valid_mask]

        # If there are no valid pixels, raise an error
        if X_valid.shape[0] == 0:
            raise ValueError(
                "No valid pixels (all pixels contain NaN values in at least one band)"
            )

        # Standardize the data (important for PCA)
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X_valid)

        # Perform PCA
        pca = PCA(n_components=n_components)
        pca_result = pca.fit_transform(X_scaled)

        # Create output array with same shape as input but with n_components bands
        output_shape = (n_components, height, width)
        pca_full = np.full(output_shape, np.nan)

        # Map the valid pixels back to their original positions
        valid_indices = np.where(valid_mask)[0]
        for i in range(n_components):
            flat_band = np.full(height * width, np.nan)
            flat_band[valid_indices] = pca_result[:, i]
            pca_full[i] = flat_band.reshape(height, width)

        # Create output directory if it doesn't exist
        output_dir = os.path.dirname(output_raster_path)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Create the output multiband GeoTIFF
        driver = gdal.GetDriverByName("GTiff")
        out_ds = driver.Create(
            output_raster_path,
            width,
            height,
            n_components,  # Number of bands = number of components
            gdal.GDT_Float32,
        )

        # Copy projection and geotransform from input
        out_ds.SetProjection(ds.GetProjection())
        out_ds.SetGeoTransform(ds.GetGeoTransform())

        # Write each PCA component as a separate band
        for i in range(n_components):
            out_band = out_ds.GetRasterBand(i + 1)
            out_band.WriteArray(pca_full[i])
            out_band.SetNoDataValue(np.nan)

            # Set band description
            variance_pct = pca.explained_variance_ratio_[i] * 100
            out_band.SetDescription(f"PC{i+1} ({variance_pct:.2f}%)")

        # Write component metadata as dataset description
        metadata = {
            "EXPLAINED_VARIANCE": ",".join(
                [f"{v:.6f}" for v in pca.explained_variance_]
            ),
            "EXPLAINED_VARIANCE_RATIO": ",".join(
                [f"{v:.6f}" for v in pca.explained_variance_ratio_]
            ),
            "LOADINGS": ";".join(
                [
                    ",".join([f"{v:.6f}" for v in component])
                    for component in pca.components_
                ]
            ),
        }

        for key, value in metadata.items():
            out_ds.SetMetadataItem(key, value)

        # Close datasets
        out_ds = None
        ds = None

        # Print information about the PCA
        print("PCA Summary:")
        print(f"Input: {input_raster_path} ({bands} bands)")
        print(f"Output: {output_raster_path} ({n_components} components)")
        print("\nExplained variance ratio by component:")
        for i, var in enumerate(pca.explained_variance_ratio_):
            print(f"PC{i+1}: {var:.4f} ({var*100:.2f}%)")
        print(
            f"\nCumulative explained variance: {np.sum(pca.explained_variance_ratio_)*100:.2f}%"
        )

        return pca.components_, pca.explained_variance_ratio_

    def ica_with_nans(
        self, input_raster_path, output_raster_path, n_components=None, random_state=42
    ):
        """
        Perform ICA on a raster with NaN values and save as a multiband GeoTIFF

        Parameters:
        input_raster_path (str): Path to the input GeoTIFF
        output_raster_path (str): Path to save the output multiband GeoTIFF with ICA components
        n_components (int): Number of independent components to calculate (defaults to all possible)
        random_state (int): Random seed for reproducibility

        Returns:
        tuple: (mixing_matrix, unmixing_matrix)
        """
        from sklearn.decomposition import FastICA
        from sklearn.preprocessing import StandardScaler

        # Load raster using PyQGIS to access metadata
        raster_layer = QgsRasterLayer(input_raster_path, "Input Raster")
        if not raster_layer.isValid():
            raise ValueError(f"Failed to load raster layer: {input_raster_path}")

        # Also open with GDAL for data access
        ds = gdal.Open(input_raster_path)
        if ds is None:
            raise ValueError(f"Failed to open raster with GDAL: {input_raster_path}")

        # Get raster dimensions
        width = ds.RasterXSize
        height = ds.RasterYSize
        bands = ds.RasterCount
        if bands == 1:
            print("Only one band found, PCA not possible")
            return None, None

        # If n_components is not specified, use all possible components
        if n_components is None or n_components == 0:
            n_components = bands
        elif n_components > bands:
            n_components = bands
            print(
                f"Warning: Requested {n_components} components but only {bands} bands available. Using {bands} components."
            )

        # Create a 3D numpy array to hold all band data
        data_array = np.zeros((bands, height, width))

        # Read all bands into the array
        for b in range(bands):
            band = ds.GetRasterBand(b + 1)
            data_array[b] = band.ReadAsArray()

        # Reshape for ICA (bands as features, pixels as samples)
        # Transpose from (bands, height, width) to (height*width, bands)
        X = data_array.reshape(bands, -1).T

        # Find valid pixels (no NaN in any band)
        valid_mask = ~np.isnan(X).any(axis=1)
        X_valid = X[valid_mask]

        # If there are no valid pixels, raise an error
        if X_valid.shape[0] == 0:
            raise ValueError(
                "No valid pixels (all pixels contain NaN values in at least one band)"
            )

        # Standardize the data (important for ICA)
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X_valid)

        # Perform ICA
        ica = FastICA(
            n_components=n_components,
            random_state=random_state,
            max_iter=1000,
            tol=0.0001,
        )
        ica_result = ica.fit_transform(X_scaled)

        # Create output array with same shape as input but with n_components bands
        output_shape = (n_components, height, width)
        ica_full = np.full(output_shape, np.nan)

        # Map the valid pixels back to their original positions
        valid_indices = np.where(valid_mask)[0]
        for i in range(n_components):
            flat_band = np.full(height * width, np.nan)
            flat_band[valid_indices] = ica_result[:, i]
            ica_full[i] = flat_band.reshape(height, width)

        # Create output directory if it doesn't exist
        output_dir = os.path.dirname(output_raster_path)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Create the output multiband GeoTIFF
        driver = gdal.GetDriverByName("GTiff")
        out_ds = driver.Create(
            output_raster_path,
            width,
            height,
            n_components,  # Number of bands = number of components
            gdal.GDT_Float32,
        )

        # Copy projection and geotransform from input
        out_ds.SetProjection(ds.GetProjection())
        out_ds.SetGeoTransform(ds.GetGeoTransform())

        # Write each ICA component as a separate band
        for i in range(n_components):
            out_band = out_ds.GetRasterBand(i + 1)
            out_band.WriteArray(ica_full[i])
            out_band.SetNoDataValue(np.nan)

            # Set band description
            out_band.SetDescription(f"IC{i+1}")

        # Write component metadata as dataset description
        metadata = {
            "MIXING_MATRIX": ";".join(
                [",".join([f"{v:.6f}" for v in row]) for row in ica.mixing_]
            ),
            "UNMIXING_MATRIX": ";".join(
                [",".join([f"{v:.6f}" for v in row]) for row in ica.components_]
            ),
        }

        for key, value in metadata.items():
            out_ds.SetMetadataItem(key, value)

        # Close datasets
        out_ds = None
        ds = None

        # Print information about the ICA
        print("ICA Summary:")
        print(f"Input: {input_raster_path} ({bands} bands)")
        print(f"Output: {output_raster_path} ({n_components} components)")
        print("\nIndependent components have been extracted.")
        print("Note: Unlike PCA, ICA components are not ordered by importance.")

        return ica.mixing_, ica.components_
