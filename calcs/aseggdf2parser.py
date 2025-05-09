import os
import re
from qgis.core import (
    QgsVectorLayer,
    QgsProject,
    QgsFeature,
    QgsWkbTypes,
    QgsProject,
    QgsGeometry,
    QgsPointXY,
    QgsVectorFileWriter,
    QgsVectorLayer,
    QgsFields,
    QgsField,
    QgsCoordinateReferenceSystem,
)

from qgis.PyQt.QtCore import (
    QVariant,
)


class AsegGdf2Parser:
    """
    Parser for ASEG-GDF2 format files (DFN, DAT).
    Based on ASEG-GDF2 standard for point located data.
    """

    def __init__(self):
        """Initialize the parser"""
        pass

    def parse_dfn_file(self, dfn_file_path):
        """
        Parse a DFN file according to ASEG-GDF2 standard to extract field definitions

        Args:
            dfn_file_path (str): Path to the .dfn file

        Returns:
            tuple: (header_list, field_defs, epsg)
        """
        print(f"Parsing DFN file: {dfn_file_path}")
        header_list = []
        field_defs = {}  # Will store field name, type, null value, etc.
        epsg = None

        try:
            with open(dfn_file_path, "r") as f:
                lines = f.readlines()
        except Exception as e:
            print(f"Error reading DFN file: {str(e)}")
            return [], {}, None

        # Print the first few lines for debugging
        """print(f"First 5 lines of DFN file:")
        for i, line in enumerate(lines[:5]):
            print(f"Line {i+1}: {line.strip()}")"""

        # Process all DEFN lines
        current_rt = None
        for line in lines:
            line = line.strip()
            if not line or not line.startswith("DEFN"):
                continue

            # Parse the DEFN line according to ASEG-GDF2 standard
            try:
                # Extract RT (Record Type)
                rt_match = re.search(r"RT=([^;]*)", line)
                if rt_match:
                    current_rt = rt_match.group(1).strip()

                # Extract field name and format
                field_match = re.search(r";\s*([^:]*):([^:]*)", line)
                if field_match and "END DEFN" not in line:
                    field_name = field_match.group(1).strip()
                    field_format = field_match.group(2).strip()

                    # Skip record type field itself
                    if field_name.upper() == "RT":
                        continue

                    # Skip DEFN end markers
                    if field_name.upper() == "END DEFN":
                        continue

                    # Add to header list if it's a regular data field
                    if current_rt == "" and field_name not in header_list:
                        header_list.append(field_name)

                    # Extract field type from format
                    field_type = QVariant.String  # Default to string
                    if field_format:
                        if field_format.startswith("I"):  # Integer
                            field_type = QVariant.Int
                        elif field_format.startswith(("F", "E", "D")):  # Float
                            field_type = QVariant.Double

                    # Extract NULL value if provided
                    null_match = re.search(r"NULL=([^,]*)", line)
                    null_value = null_match.group(1) if null_match else None

                    # Extract units if provided
                    unit_match = re.search(r"UNIT=([^,]*)", line)
                    unit = unit_match.group(1) if unit_match else None

                    # Check for EPSG code (for coordinate fields)
                    epsg_match = re.search(r"epsgcode=(\d+)", line)
                    if epsg_match and (
                        field_name.lower()
                        in ["longitude", "latitude", "easting", "northing"]
                    ):
                        epsg = epsg_match.group(1)
                        print(f"Found EPSG code {epsg} for field {field_name}")

                    # Create field definition
                    field_defs[field_name] = {
                        "type": field_type,
                        "format": field_format,
                        "null": null_value,
                        "unit": unit,
                    }
            except Exception as e:
                print(f"Error parsing DEFN line: {line}")
                print(f"Error details: {str(e)}")

        # Check if we found fields
        if header_list:
            print(f"Found {len(header_list)} fields in DFN file: {header_list}")
        else:
            print("No fields found in DFN file")

        return header_list, field_defs, epsg

    def parse_dat_file(self, dat_file_path, header_list, field_defs):
        """
        Parse a DAT file according to ASEG-GDF2 standard

        Args:
            dat_file_path (str): Path to the .dat file
            header_list (list): List of field names
            field_defs (dict): Dictionary of field definitions

        Returns:
            list: List of dictionaries, one per row
        """
        print(f"Parsing DAT file: {dat_file_path}")
        data = []

        try:
            with open(dat_file_path, "r") as f:
                lines = f.readlines()
        except Exception as e:
            print(f"Error reading DAT file: {str(e)}")
            return []

        # Print the first few lines for debugging
        """ print(f"First 5 lines of DAT file:")
        for i, line in enumerate(lines[:5]):
            print(f"Line {i+1}: {line.strip()}")"""

        # Check if records include a record type identifier
        has_rt = False
        if lines and len(lines) > 0:
            first_line = lines[0].strip()
            # ASEG-GDF2 allows one record type to have no identifier
            # Check if the first word is a 4-character RT identifier
            words = first_line.split()
            if words and len(words[0]) <= 4 and not words[0][0].isdigit():
                has_rt = True
                print(
                    f"Data appears to have record type identifiers. First RT: {words[0]}"
                )

        for line_num, line in enumerate(lines):
            line = line.strip()
            if not line:
                continue

            try:
                # Split line into fields
                fields = line.split()

                # Skip if not enough fields
                if len(fields) < len(header_list):
                    print(
                        f"Warning: Line {line_num+1} has {len(fields)} fields, expected {len(header_list)}"
                    )
                    continue

                # Create data row
                row = {}
                for i, field_name in enumerate(header_list):
                    if i < len(fields):
                        value = fields[i]

                        # Convert to appropriate type based on field definition
                        if field_name in field_defs:
                            field_def = field_defs[field_name]

                            # Check for NULL value
                            if field_def["null"] and value == field_def["null"]:
                                row[field_name] = None
                                continue

                            # Convert based on type
                            if field_def["type"] == QVariant.Int:
                                try:
                                    row[field_name] = int(value)
                                except ValueError:
                                    row[field_name] = value
                            elif field_def["type"] == QVariant.Double:
                                try:
                                    row[field_name] = float(value)
                                except ValueError:
                                    row[field_name] = value
                            else:
                                row[field_name] = value
                        else:
                            # Default handling if field not defined
                            try:
                                # Try to convert to numeric if possible
                                if "." in value or "e" in value.lower():
                                    row[field_name] = float(value)
                                else:
                                    row[field_name] = int(value)
                            except ValueError:
                                row[field_name] = value
                    else:
                        row[field_name] = None

                data.append(row)

                # Print first row for debugging
                """if line_num == 0:
                    print(f"First data row: {row}")"""

            except Exception as e:
                print(f"Error parsing line {line_num+1}: {str(e)}")

        print(f"Parsed {len(data)} rows from DAT file")
        return data

    def process_aseg_gdf2_data(
        self, dat_file_path, dfn_file_path=None, x_field_name=None, y_field_name=None
    ):
        """
        Process ASEG-GDF2 format files and return the parsed data

        Args:
            dat_file_path (str): Path to the .dat file
            dfn_file_path (str, optional): Path to the .dfn file (if None, derived from dat_file_path)
            x_field_name (str, optional): Field name for X/longitude coordinate
            y_field_name (str, optional): Field name for Y/latitude coordinate

        Returns:
            tuple: (data, header_list, x_field, y_field, epsg)
        """
        # Input validation
        if not dat_file_path:
            print("Error: No DAT file path provided")
            return None, None, None, None, None

        # Derive paths if not provided
        if not dfn_file_path:
            dfn_file_path = dat_file_path.rsplit(".", 1)[0] + ".dfn"

        print(f"Processing ASEG-GDF2 files:")
        print(f"  DAT file: {dat_file_path}")
        print(f"  DFN file: {dfn_file_path}")

        # Check if files exist
        if not os.path.exists(dat_file_path):
            print(f"Error: DAT file not found: {dat_file_path}")
            return None, None, None, None, None

        if not os.path.exists(dfn_file_path):
            print(f"Error: DFN file not found: {dfn_file_path}")
            return None, None, None, None, None

        # Parse the DFN file to get field definitions
        header_list, field_defs, epsg = self.parse_dfn_file(dfn_file_path)

        # Parse the DAT file to get the data
        data = self.parse_dat_file(dat_file_path, header_list, field_defs)

        # If no data was found, return None
        if not data:
            print("Error: No data found in the file")
            return None, None, None, None, None

        # Determine coordinate field names if not provided
        if x_field_name is None:
            # Try to find X/longitude field in header list
            for field in header_list:
                if (
                    field.lower() == "longitude"
                    or field.lower() == "x"
                    or field.lower() == "easting"
                    or "x_" in field.lower()
                    or field.lower().endswith("_x")
                    or "east" in field.lower()
                ):
                    x_field_name = field
                    print(f"Using {x_field_name} as X coordinate field")
                    break

        if y_field_name is None:
            # Try to find Y/latitude field in header list
            for field in header_list:
                if (
                    field.lower() == "latitude"
                    or field.lower() == "y"
                    or field.lower() == "northing"
                    or "y_" in field.lower()
                    or field.lower().endswith("_y")
                    or "north" in field.lower()
                ):
                    y_field_name = field
                    print(f"Using {y_field_name} as Y coordinate field")
                    break

        # Verify we have coordinate fields
        if x_field_name is None or y_field_name is None:
            print(f"Error: Could not determine coordinate fields in {header_list}")
            return None, None, None, None, None

        # Verify that coordinate fields exist in the data
        if data and (x_field_name not in data[0] or y_field_name not in data[0]):
            print(
                f"Warning: coordinate fields {x_field_name}, {y_field_name} not found in data. Available fields: {list(data[0].keys())}"
            )
            return None, None, None, None, None

        # Print some summary information
        print(f"Data successfully processed:")
        print(f"  Number of fields: {len(header_list)}")
        print(f"  Number of records: {len(data)}")
        print(f"  X coordinate field: {x_field_name}")
        print(f"  Y coordinate field: {y_field_name}")
        if epsg:
            print(f"  EPSG code: {epsg}")

        # Return the processed data
        return data, header_list, x_field_name, y_field_name, epsg

    def export_to_csv(self, data, header_list, output_file):
        """
        Export data to a CSV file

        Args:
            data (list): List of dictionaries containing the data
            header_list (list): List of field names for the header
            output_file (str): Path to the output file

        Returns:
            bool: True if successful, False otherwise
        """
        if not data or not header_list:
            print("Error: No data to export")
            return False

        try:
            import csv

            # Create output directory if it doesn't exist
            output_dir = os.path.dirname(output_file)
            if output_dir and not os.path.exists(output_dir):
                os.makedirs(output_dir)

            with open(output_file, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=header_list)
                writer.writeheader()
                writer.writerows(data)

            print(f"Data exported to CSV: {output_file}")
            print(f"  {len(data)} records written")
            return True
        except Exception as e:
            print(f"Error exporting to CSV: {str(e)}")
            return False

    def export_to_shapefile(
        self, data, header_list, x_field_name, y_field_name, output_file, epsg=None
    ):
        """
        Export data to a shapefile using QGIS API

        Args:
            data (list): List of dictionaries containing the data
            header_list (list): List of field names
            x_field_name (str): Field name for X coordinate
            y_field_name (str): Field name for Y coordinate
            output_file (str): Path to the output shapefile (with or without .shp extension)
            epsg (str, optional): EPSG code for the projection

        Returns:
            QgsVectorLayer: The created layer
        """

        # Check if we have data
        if not data or not header_list:
            print("Error: No data to export")
            return None

        # Get output directory
        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.exists(output_dir):
            try:
                os.makedirs(output_dir)
                print(f"Created output directory: {output_dir}")
            except Exception as e:
                print(f"Error creating output directory: {str(e)}")
                return None

        # Make sure output_file has .shp extension
        if not output_file.lower().endswith(".shp"):
            output_file = output_file + ".shp"

        # Get layer name from the file path
        layer_name = os.path.splitext(os.path.basename(output_file))[0]

        # Create CRS
        crs = QgsCoordinateReferenceSystem("EPSG:4326")  # Default to WGS84
        if epsg:
            crs = QgsCoordinateReferenceSystem(f"EPSG:{epsg}")
            print(f"Using CRS: {crs.description()}")

        # Create fields
        fields = QgsFields()
        for field_name in header_list:
            # Determine field type from data
            field_type = QVariant.String  # Default to string

            for row in data:
                if field_name in row and row[field_name] is not None:
                    value = row[field_name]
                    if isinstance(value, int):
                        field_type = QVariant.Int
                        break
                    elif isinstance(value, float):
                        field_type = QVariant.Double
                        break

            fields.append(QgsField(field_name, field_type))

        print(f"Creating shapefile: {output_file}")

        # Create the shapefile writer
        writer_options = QgsVectorFileWriter.SaveVectorOptions()
        writer_options.driverName = "ESRI Shapefile"
        writer_options.fileEncoding = "UTF-8"

        transform_context = QgsProject.instance().transformContext()

        writer = QgsVectorFileWriter.create(
            output_file,
            fields,
            QgsWkbTypes.Point,
            crs,
            transform_context,
            writer_options,
        )

        # Add features
        feature_count = 0
        error_count = 0

        for i, row in enumerate(data):
            try:
                # Check if coordinate fields have values
                if x_field_name not in row or y_field_name not in row:
                    continue

                if row[x_field_name] is None or row[y_field_name] is None:
                    continue

                # Get coordinates
                x_coord = float(row[x_field_name])
                y_coord = float(row[y_field_name])

                # Create feature
                feat = QgsFeature(fields)

                # Set geometry
                feat.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(x_coord, y_coord)))

                # Set attributes
                for j, field_name in enumerate(header_list):
                    feat.setAttribute(j, row.get(field_name, None))

                # Add to shapefile
                writer.addFeature(feat)
                feature_count += 1

            except Exception as e:
                error_count += 1
                if error_count < 10:  # Limit the number of error messages
                    print(f"Error with row {i}: {e}")

        # Clean up writer
        del writer

        print(f"Features written to shapefile: {feature_count}")
        print(f"Errors encountered: {error_count}")

        # Load the shapefile
        layer = QgsVectorLayer(output_file, layer_name, "ogr")

        if not layer.isValid():
            print(f"Layer is not valid! Path: {output_file}")
            return None

        print(f"Shapefile loaded with {layer.featureCount()} features")

        # Add to project
        QgsProject.instance().addMapLayer(layer)

        return layer
