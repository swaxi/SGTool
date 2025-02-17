import numpy as np
from scipy.ndimage import convolve, median_filter, gaussian_filter
from PyQt5.QtGui import QValidator
from qgis.core import (
    QgsProject,
    QgsRasterLayer,
    QgsVectorLayer,
    QgsFeature,
    QgsGeometry,
    QgsField,
    QgsFields,
    QgsPointXY,
    QgsVectorFileWriter,
)
from qgis.PyQt.QtCore import QVariant
import numpy as np
from osgeo import gdal
from collections import defaultdict


class SG_Util:
    def __init__(self, grid):
        """
        Initialize the Utils with a grid.

        :param grid: 2D numpy array representing the input grid
        """
        self.grid = np.array(grid, dtype=float)

    def Threshold2Nan(
        self, grid, condition, above_threshold_value, below_threshold_value
    ):
        # Print types for debugging

        # Ensure grid is a floating-point array
        if not np.issubdtype(grid.dtype, np.floating):
            grid = grid.astype(float)

        # Process the raster using numpy
        if condition == "above":
            # Set values above the threshold to NaN
            grid2 = np.where(grid > above_threshold_value, np.nan, grid)
        elif condition == "below":
            # Set values below the threshold to NaN
            grid2 = np.where(grid < below_threshold_value, np.nan, grid)
        elif condition == "between":
            # Set values outside the thresholds to NaN
            grid2 = np.where(
                (grid < above_threshold_value) & (grid > below_threshold_value),
                np.nan,
                grid,
            )
        else:
            raise ValueError("Invalid condition. Use 'above', 'below', or 'both'.")

        return grid2

    def create_data_boundary_lines(self, current_layer, output_path=None):
        """
        Creates a polygon layer that traces the boundaries between data and NaN regions
        in the currently selected raster layer, handling holes as interior rings.
        """

        if not current_layer.isValid():
            print("Current layer not valid")
            return

        # Get the data provider and read the raster
        provider = current_layer.dataProvider()
        extent = current_layer.extent()
        rows = current_layer.height()
        cols = current_layer.width()
        block = provider.block(1, extent, cols, rows)

        # Convert raster data to numpy array
        data = np.zeros((rows, cols))
        for row in range(rows):
            for col in range(cols):
                data[row, col] = block.value(row, col)

        # Get the NoData value if it exists
        nodata_value = provider.sourceNoDataValue(1)

        # Create binary mask (1 for data, 0 for NaN or NoData)
        mask = np.where(np.logical_or(np.isnan(data), data == nodata_value), 0, 1)

        # Get raster georeference information
        transform = [
            extent.xMinimum(),  # x origin
            current_layer.rasterUnitsPerPixelX(),  # pixel width
            0,  # rotation, 0 if image is "north up"
            extent.yMaximum(),  # y origin
            0,  # rotation, 0 if image is "north up"
            -current_layer.rasterUnitsPerPixelY(),  # pixel height
        ]

        def pixel_to_map_coords(row, col):
            """Convert raster coordinates to map coordinates"""
            x = transform[0] + col * transform[1]
            y = transform[3] + row * transform[5]
            return QgsPointXY(x, y)

        def get_boundary_points():
            """Extract boundary points and their connections"""
            points = set()
            edges = []

            # Find boundaries using numpy array operations
            vertical_edges = np.diff(mask, axis=1)
            horizontal_edges = np.diff(mask, axis=0)

            # Process vertical edges
            edge_rows, edge_cols = np.where(vertical_edges != 0)
            for row, col in zip(edge_rows, edge_cols):
                p1 = pixel_to_map_coords(row, col + 1)
                p2 = pixel_to_map_coords(row + 1, col + 1)
                points.add((p1.x(), p1.y()))
                points.add((p2.x(), p2.y()))
                edges.append(((p1.x(), p1.y()), (p2.x(), p2.y())))

            # Process horizontal edges
            edge_rows, edge_cols = np.where(horizontal_edges != 0)
            for row, col in zip(edge_rows, edge_cols):
                p1 = pixel_to_map_coords(row + 1, col)
                p2 = pixel_to_map_coords(row + 1, col + 1)
                points.add((p1.x(), p1.y()))
                points.add((p2.x(), p2.y()))
                edges.append(((p1.x(), p1.y()), (p2.x(), p2.y())))

            # Add raster edge connections
            for row in range(rows):
                if mask[row, 0] == 1:  # Left edge
                    p1 = pixel_to_map_coords(row, 0)
                    p2 = pixel_to_map_coords(row + 1, 0)
                    points.add((p1.x(), p1.y()))
                    points.add((p2.x(), p2.y()))
                    edges.append(((p1.x(), p1.y()), (p2.x(), p2.y())))

                if mask[row, -1] == 1:  # Right edge
                    p1 = pixel_to_map_coords(row, cols)
                    p2 = pixel_to_map_coords(row + 1, cols)
                    points.add((p1.x(), p1.y()))
                    points.add((p2.x(), p2.y()))
                    edges.append(((p1.x(), p1.y()), (p2.x(), p2.y())))

            for col in range(cols):
                if mask[0, col] == 1:  # Top edge
                    p1 = pixel_to_map_coords(0, col)
                    p2 = pixel_to_map_coords(0, col + 1)
                    points.add((p1.x(), p1.y()))
                    points.add((p2.x(), p2.y()))
                    edges.append(((p1.x(), p1.y()), (p2.x(), p2.y())))

                if mask[-1, col] == 1:  # Bottom edge
                    p1 = pixel_to_map_coords(rows, col)
                    p2 = pixel_to_map_coords(rows, col + 1)
                    points.add((p1.x(), p1.y()))
                    points.add((p2.x(), p2.y()))
                    edges.append(((p1.x(), p1.y()), (p2.x(), p2.y())))

            return points, edges

        def points_equal(p1, p2, tolerance=1e-10):
            """Check if two points are equal within tolerance"""
            return abs(p1[0] - p2[0]) < tolerance and abs(p1[1] - p2[1]) < tolerance

        def build_adjacency_graph(edges):
            """Build graph of connected points"""
            graph = defaultdict(list)
            for start, end in edges:
                graph[start].append(end)
                graph[end].append(start)
            return graph

        def trace_boundary(start_point, graph, used_edges):
            """Trace a complete boundary from a starting point"""
            current = start_point
            points = [current]

            while True:
                neighbors = [
                    p
                    for p in graph[current]
                    if (current, p) not in used_edges and (p, current) not in used_edges
                ]

                if not neighbors:
                    break

                next_point = neighbors[0]
                used_edges.add((current, next_point))
                points.append(next_point)
                current = next_point

                if points_equal(points[0], points[-1]):
                    break

            return points

        def is_clockwise(points):
            """Check if a polygon is clockwise"""
            area = 0
            for i in range(len(points) - 1):
                area += (points[i + 1][0] - points[i][0]) * (
                    points[i + 1][1] + points[i][1]
                )
            return area > 0

        def identify_holes(rings):
            """Identify which rings are holes in other rings"""
            from shapely.geometry import Polygon

            # Convert rings to Shapely polygons
            polygons = [Polygon(ring) for ring in rings]

            # Find containment relationships
            holes = []
            exteriors = []

            for i, poly1 in enumerate(polygons):
                is_hole = False
                for j, poly2 in enumerate(polygons):
                    if i != j and poly2.contains(poly1):
                        holes.append((i, j))  # ring i is a hole in ring j
                        is_hole = True
                        break
                if not is_hole:
                    exteriors.append(i)

            return exteriors, holes

        # Create output vector layer
        fields = QgsFields()
        fields.append(QgsField("id", QVariant.Int))
        fields.append(QgsField("type", QVariant.String))

        vector_layer = QgsVectorLayer(
            "MultiPolygon?crs=" + current_layer.crs().authid(),
            "data_boundaries",
            "memory",
        )
        provider = vector_layer.dataProvider()
        provider.addAttributes(fields)
        vector_layer.updateFields()

        # Get boundary points and edges
        points, edges = get_boundary_points()
        graph = build_adjacency_graph(edges)
        used_edges = set()

        # Trace all rings
        rings = []
        for start_point in points:
            if any(
                (start_point, end) in used_edges or (end, start_point) in used_edges
                for end in graph[start_point]
            ):
                continue

            ring_points = trace_boundary(start_point, graph, used_edges)
            if len(ring_points) > 3:  # Need at least 4 points for a valid ring
                rings.append(ring_points)

        # Identify holes and exterior rings
        exteriors, holes = identify_holes(rings)

        # Create features with proper geometry
        features = []
        feature_id = 1

        # Group holes with their exterior rings
        hole_map = defaultdict(list)
        for hole_idx, exterior_idx in holes:
            hole_map[exterior_idx].append(rings[hole_idx])

        # Create multipolygon features
        for exterior_idx in exteriors:
            exterior_ring = rings[exterior_idx]
            if not is_clockwise(exterior_ring):
                exterior_ring = exterior_ring[::-1]

            # Get any holes for this exterior
            interior_rings = hole_map[exterior_idx]
            # Ensure holes are counterclockwise
            interior_rings = [
                ring if is_clockwise(ring) else ring[::-1] for ring in interior_rings
            ]

            # Convert to QgsPointXY
            exterior_points = [QgsPointXY(x, y) for x, y in exterior_ring]
            interior_points = [
                [QgsPointXY(x, y) for x, y in ring] for ring in interior_rings
            ]

            # Create geometry
            feature = QgsFeature()
            geom = QgsGeometry.fromPolygonXY([exterior_points] + interior_points)
            feature.setGeometry(geom)
            feature.setAttributes([feature_id, "polygon_with_holes"])
            features.append(feature)
            feature_id += 1

        # Add features to the layer
        provider.addFeatures(features)

        # Save the layer if output path is provided
        if output_path:
            error = QgsVectorFileWriter.writeAsVectorFormat(
                vector_layer, output_path, "UTF-8", vector_layer.crs(), "ESRI Shapefile"
            )

            if error[0] != QgsVectorFileWriter.NoError:
                print(f"Error saving vector layer: {error}")

        # Add the layer to the project
        QgsProject.instance().addMapLayer(vector_layer)

        return vector_layer
