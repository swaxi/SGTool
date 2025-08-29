import numpy as np
from collections import defaultdict
from osgeo import gdal, ogr, osr
from shapely.geometry import Polygon
import os
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

    def create_data_boundary_lines(
        self, input_raster_path, output_path=None, band_number=1
    ):
        """
        Creates a polygon layer that traces the boundaries between data and NaN regions
        in the specified raster file, handling holes as interior rings.

        Args:
            input_raster_path (str): Path to the input raster file
            output_path (str, optional): Path for output shapefile
            band_number (int): Band number to process (default: 1)

        Returns:
            str: Path to the created shapefile (if output_path provided)
        """

        # Open the raster dataset
        dataset = gdal.Open(input_raster_path, gdal.GA_ReadOnly)
        if not dataset:
            raise ValueError(f"Could not open raster file: {input_raster_path}")

        # Get raster information
        band = dataset.GetRasterBand(band_number)
        cols = dataset.RasterXSize
        rows = dataset.RasterYSize
        geotransform = dataset.GetGeoTransform()
        projection = dataset.GetProjection()

        # Read the raster data
        data = band.ReadAsArray()
        nodata_value = band.GetNoDataValue()

        # Create binary mask (1 for data, 0 for NaN or NoData)
        if nodata_value is not None:
            mask = np.where(np.logical_or(np.isnan(data), data == nodata_value), 0, 1)
        else:
            mask = np.where(np.isnan(data), 0, 1)

        def pixel_to_map_coords(row, col):
            """Convert raster coordinates to map coordinates using geotransform"""
            x = geotransform[0] + col * geotransform[1] + row * geotransform[2]
            y = geotransform[3] + col * geotransform[4] + row * geotransform[5]
            return (x, y)

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
                points.add(p1)
                points.add(p2)
                edges.append((p1, p2))

            # Process horizontal edges
            edge_rows, edge_cols = np.where(horizontal_edges != 0)
            for row, col in zip(edge_rows, edge_cols):
                p1 = pixel_to_map_coords(row + 1, col)
                p2 = pixel_to_map_coords(row + 1, col + 1)
                points.add(p1)
                points.add(p2)
                edges.append((p1, p2))

            # Add raster edge connections
            for row in range(rows):
                if mask[row, 0] == 1:  # Left edge
                    p1 = pixel_to_map_coords(row, 0)
                    p2 = pixel_to_map_coords(row + 1, 0)
                    points.add(p1)
                    points.add(p2)
                    edges.append((p1, p2))

                if mask[row, -1] == 1:  # Right edge
                    p1 = pixel_to_map_coords(row, cols)
                    p2 = pixel_to_map_coords(row + 1, cols)
                    points.add(p1)
                    points.add(p2)
                    edges.append((p1, p2))

            for col in range(cols):
                if mask[0, col] == 1:  # Top edge
                    p1 = pixel_to_map_coords(0, col)
                    p2 = pixel_to_map_coords(0, col + 1)
                    points.add(p1)
                    points.add(p2)
                    edges.append((p1, p2))

                if mask[-1, col] == 1:  # Bottom edge
                    p1 = pixel_to_map_coords(rows, col)
                    p2 = pixel_to_map_coords(rows, col + 1)
                    points.add(p1)
                    points.add(p2)
                    edges.append((p1, p2))

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
            # Convert rings to Shapely polygons
            polygons = [Polygon(ring) for ring in rings if len(ring) > 3]

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

        # Create output shapefile if path provided
        if output_path:
            # Create the output driver
            driver = ogr.GetDriverByName("ESRI Shapefile")

            # Remove existing file if it exists
            if os.path.exists(output_path):
                driver.DeleteDataSource(output_path)

            # Create the data source
            data_source = driver.CreateDataSource(output_path)

            # Create spatial reference from raster projection
            srs = osr.SpatialReference()
            srs.ImportFromWkt(projection)

            # Create the layer
            layer_name = os.path.splitext(os.path.basename(output_path))[0]
            layer = data_source.CreateLayer(layer_name, srs, ogr.wkbMultiPolygon)

            # Add fields
            id_field = ogr.FieldDefn("id", ogr.OFTInteger)
            type_field = ogr.FieldDefn("type", ogr.OFTString)
            type_field.SetWidth(50)

            layer.CreateField(id_field)
            layer.CreateField(type_field)

            # Group holes with their exterior rings
            hole_map = defaultdict(list)
            for hole_idx, exterior_idx in holes:
                hole_map[exterior_idx].append(rings[hole_idx])

            # Create features
            feature_id = 1
            for exterior_idx in exteriors:
                exterior_ring = rings[exterior_idx]
                if not is_clockwise(exterior_ring):
                    exterior_ring = exterior_ring[::-1]

                # Get any holes for this exterior
                interior_rings = hole_map[exterior_idx]
                # Ensure holes are counterclockwise
                interior_rings = [
                    ring if is_clockwise(ring) else ring[::-1]
                    for ring in interior_rings
                ]

                # Create OGR polygon
                polygon = ogr.Geometry(ogr.wkbPolygon)

                # Add exterior ring
                exterior = ogr.Geometry(ogr.wkbLinearRing)
                for x, y in exterior_ring:
                    exterior.AddPoint(x, y)
                polygon.AddGeometry(exterior)

                # Add interior rings (holes)
                for ring in interior_rings:
                    interior = ogr.Geometry(ogr.wkbLinearRing)
                    for x, y in ring:
                        interior.AddPoint(x, y)
                    polygon.AddGeometry(interior)

                # Create feature
                feature = ogr.Feature(layer.GetLayerDefn())
                feature.SetGeometry(polygon)
                feature.SetField("id", feature_id)
                feature.SetField("type", "Multipolygon boundary")

                layer.CreateFeature(feature)
                feature = None  # Clean up
                feature_id += 1

            # Clean up
            data_source = None

            print(f"Boundary shapefile created: {output_path}")
            return output_path

        else:
            # Return the rings data for further processing
            return rings, exteriors, holes
