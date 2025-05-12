import os
import re
from qgis.core import (
    Qgis,
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
                # Critical fix: Changed regex to properly handle fields on the same line as END DEFN
                field_match = re.search(r";\s*([^:]*):([^:;]*)", line)

                if field_match:
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
                        print(f"Added field to header list: {field_name}")

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
                        in [
                            "longitude",
                            "latitude",
                            "easting",
                            "northing",
                            "x_gda94",
                            "y_gda94",
                        ]
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
            print(f"Found {len(header_list)} fields in DFN file:")
            for i, field in enumerate(header_list):
                field_def = field_defs.get(field, {})
                field_type = "String"
                if field_def.get("type") == QVariant.Int:
                    field_type = "Integer"
                elif field_def.get("type") == QVariant.Double:
                    field_type = "Double"
                print(f"  {i+1}. {field} ({field_type})")
        else:
            print("No fields found in DFN file")

        # Debug - check specifically for Y_GDA94
        if "Y_GDA94" in header_list:
            print(f"Y_GDA94 found at position {header_list.index('Y_GDA94') + 1}")
        else:
            print("WARNING: Y_GDA94 field not found in DFN file!")

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

        # Track stats for reporting
        processed_lines = 0
        skipped_lines = 0

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
            if not line:  # Skip empty lines
                skipped_lines += 1
                continue

            try:
                # Split line into fields
                fields = line.split()

                # Skip if not enough fields
                if len(fields) < len(header_list):
                    print(
                        f"Warning: Line {line_num+1} has {len(fields)} fields, expected {len(header_list)}"
                    )
                    skipped_lines += 1
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
                processed_lines += 1

                # Debug first and last row for verification
                if processed_lines == 1 or line_num == len(lines) - 1:
                    field_samples = [
                        f"{field_name}={row.get(field_name, 'N/A')}"
                        for field_name in header_list[-2:]
                    ]
                    print(
                        f"{'First' if processed_lines == 1 else 'Last'} data row processed: FID={row.get('FID', 'N/A')}, {', '.join(field_samples)}"
                    )

            except Exception as e:
                print(f"Error parsing line {line_num+1}: {str(e)}")
                skipped_lines += 1

        print(
            f"Parsed {len(data)} rows from DAT file (Processed: {processed_lines}, Skipped: {skipped_lines})"
        )

        # Additional verification for last record
        if data and len(lines) > 0:
            last_line = lines[-1].strip()
            last_fields = last_line.split()
            last_fid = last_fields[2] if len(last_fields) > 2 else "unknown"

            # Verify the last record was processed
            last_record_processed = any(
                str(row.get("FID", "")).startswith(last_fid) for row in data[-1:]
            )
            if not last_record_processed:
                print(
                    f"WARNING: Last record (FID={last_fid}) may not have been processed!"
                )

        # Verify Y_GDA94 field in all rows
        missing_y_gda94_count = sum(
            1 for row in data if "Y_GDA94" not in row or row["Y_GDA94"] is None
        )
        if missing_y_gda94_count > 0:
            print(
                f"WARNING: {missing_y_gda94_count} rows are missing the Y_GDA94 field value"
            )

        # Report on coordinate field values (for debugging)
        last_row = data[-1] if data else {}
        x_field = "X_GDA94"
        y_field = "Y_GDA94"
        print(
            f"Last row coordinates: {x_field}={last_row.get(x_field, 'N/A')}, {y_field}={last_row.get(y_field, 'N/A')}"
        )

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
                    or field.lower() == "x_gda94"
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
                    or field.lower() == "y_gda94"
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

        # Verify coordinate values in the data
        if data:
            first_row = data[0]
            last_row = data[-1]
            print(
                f"  First row coordinates: {x_field_name}={first_row.get(x_field_name)}, {y_field_name}={first_row.get(y_field_name)}"
            )
            print(
                f"  Last row coordinates: {x_field_name}={last_row.get(x_field_name)}, {y_field_name}={last_row.get(y_field_name)}"
            )

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
            print(f"Adding field to shapefile: {field_name} (Type: {field_type})")

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

        if writer.hasError() != QgsVectorFileWriter.NoError:
            print(f"Error creating shapefile writer: {writer.errorMessage()}")
            return None

        # Add features
        feature_count = 0
        error_count = 0
        missing_coord_count = 0

        for i, row in enumerate(data):
            try:
                # Check if coordinate fields have values
                if x_field_name not in row or y_field_name not in row:
                    print(
                        f"Warning: Row {i+1} missing coordinate field. Available fields: {list(row.keys())}"
                    )
                    missing_coord_count += 1
                    continue

                if row[x_field_name] is None or row[y_field_name] is None:
                    print(
                        f"Warning: Row {i+1} has null coordinate values: {x_field_name}={row[x_field_name]}, {y_field_name}={row[y_field_name]}"
                    )
                    missing_coord_count += 1
                    continue

                # Get coordinates
                try:
                    x_coord = float(row[x_field_name])
                    y_coord = float(row[y_field_name])
                except (ValueError, TypeError) as e:
                    print(
                        f"Warning: Row {i+1} has invalid coordinate values: {x_field_name}={row[x_field_name]}, {y_field_name}={row[y_field_name]}"
                    )
                    print(f"Error: {str(e)}")
                    missing_coord_count += 1
                    continue

                # Create feature
                feat = QgsFeature(fields)

                # Set geometry
                feat.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(x_coord, y_coord)))

                # Set attributes
                for j, field_name in enumerate(header_list):
                    feat.setAttribute(j, row.get(field_name, None))

                # Add to shapefile
                if writer.addFeature(feat):
                    feature_count += 1
                    # Log first and last feature for verification
                    if feature_count == 1 or i == len(data) - 1:
                        coord_msg = f"({x_coord}, {y_coord})"
                        field_samples = [
                            f"{field_name}={row.get(field_name, 'N/A')}"
                            for field_name in header_list[-2:]
                        ]
                        print(
                            f"{'First' if feature_count == 1 else 'Last'} feature written: FID={row.get('FID', 'N/A')}, Coords={coord_msg}, {', '.join(field_samples)}"
                        )
                else:
                    error_count += 1
                    print(f"Error adding feature {i+1}")

            except Exception as e:
                error_count += 1
                if error_count < 10:  # Limit the number of error messages
                    print(f"Error with row {i}: {e}")

        # Clean up writer
        del writer

        print(f"Features written to shapefile: {feature_count}")
        print(f"Errors encountered: {error_count}")
        print(f"Records skipped due to missing coordinates: {missing_coord_count}")

        # Verify the last row was written
        if data and feature_count > 0:
            last_row = data[-1]
            last_fid = last_row.get("FID", "N/A")
            print(f"Last data row FID: {last_fid}")

            if feature_count < len(data):
                print(
                    f"Warning: Not all records were written to the shapefile ({feature_count} of {len(data)})"
                )

        # Load the shapefile
        layer = QgsVectorLayer(output_file, layer_name, "ogr")

        if not layer.isValid():
            print(f"Layer is not valid! Path: {output_file}")
            return None

        print(f"Shapefile loaded with {layer.featureCount()} features")

        # Add to project
        QgsProject.instance().addMapLayer(layer)

        return layer

    def parse_legacy_headers(self, dat_file_path):
        """
        Extract header information, projection data, and actual data from related geophysical data files.
        """
        # Derive paths for the .dfn and .prj files
        base_name = dat_file_path.rsplit(".", 1)[0]
        if dat_file_path.rsplit(".", 1)[1] == ".DAT":
            dfn_file_path = base_name + ".DFN"
            prj_file_path = base_name + ".PRJ"
        else:
            dfn_file_path = base_name + ".dfn"
            prj_file_path = base_name + ".prj"

        if not os.path.exists(dat_file_path):
            self.iface.messageBar().pushMessage(
                "DAT/dat file not found", level=Qgis.Warning, duration=3
            )
            return None, None, None, None

        if not os.path.exists(dfn_file_path):
            self.iface.messageBar().pushMessage(
                "DFN/dfn file not found", level=Qgis.Warning, duration=3
            )
            return None, None, None, None

        # Parse the DFN file to get field definitions
        field_defs = []
        data_format_type = "unknown"

        with open(dfn_file_path, "r") as dfn_file:
            dfn_content = dfn_file.readlines()

            # Determine the format type by checking the first DEFN line
            for line in dfn_content:
                if "ST=RECD,RT=DATA;RT:A4;" in line:
                    data_format_type = "data_as_identifier"
                    break

            # Process all field definitions
            for line in dfn_content:
                line = line.strip()
                if not line:
                    continue

                if line.startswith("DEFN") and (
                    "ST=RECD,RT=;" in line or "ST=RECD,RT=DATA;" in line
                ):
                    parts = line.split(";")
                    field_info = {}

                    # Extract field index if available
                    if " " in line:
                        defn_parts = line.split(" ", 2)
                        if len(defn_parts) > 1 and defn_parts[1].strip().isdigit():
                            field_info["index"] = int(defn_parts[1].strip())

                    # Extract field name and format
                    for part in parts:
                        if ":" in part:
                            name_format = part.split(":", 2)
                            if len(name_format) > 1:
                                field_name = name_format[0].strip()
                                field_format = (
                                    name_format[1].strip()
                                    if len(name_format) > 1
                                    else ""
                                )

                                # Extract field properties - preserve exact field name
                                if "NAME=" in part:
                                    name_part = (
                                        part.split("NAME=")[1].split(",")[0].strip()
                                        if "," in part.split("NAME=")[1]
                                        else part.split("NAME=")[1].strip()
                                    )
                                    field_info["name"] = (
                                        name_part  # Use exact name as defined
                                    )
                                else:
                                    field_info["name"] = field_name

                                # Parse format
                                if field_format and field_format[0] in [
                                    "F",
                                    "I",
                                    "D",
                                    "A",
                                ]:
                                    field_info["type"] = field_format[0]

                                    # Extract width and precision
                                    format_parts = field_format.split(":", 1)[0].strip()
                                    width_parts = "".join(
                                        [
                                            c
                                            for c in format_parts
                                            if c.isdigit() or c == "."
                                        ]
                                    )

                                    if "." in width_parts:
                                        width, precision = width_parts.split(".")
                                        field_info["width"] = (
                                            int(width) if width else 10
                                        )
                                        field_info["precision"] = (
                                            int(precision) if precision else 0
                                        )
                                    else:
                                        field_info["width"] = (
                                            int(width_parts)
                                            if width_parts.isdigit()
                                            else 10
                                        )
                                        field_info["precision"] = 0

                                    # Only add if we have valid field info
                                    if "name" in field_info and "type" in field_info:
                                        field_defs.append(field_info)
                                        break

        # Sort field definitions by index if available
        if field_defs and all("index" in field for field in field_defs):
            field_defs.sort(key=lambda x: x["index"])

        # Extract header list from field definitions (preserving exact names)
        header_list = [
            field.get("name", f"Field_{i}") for i, field in enumerate(field_defs)
        ]

        print(f"Found {len(header_list)} fields: {header_list}")
        print(f"Data format type: {data_format_type}")

        # Determine EPSG code from projection info
        epsg = None

        # First priority: Check PRJ file if it exists
        if os.path.exists(prj_file_path):
            with open(prj_file_path, "r") as prj_file:
                prj_content = prj_file.read()
                if "GDA94" in prj_content and "zone" in prj_content.lower():
                    try:
                        zone_text = (
                            prj_content.lower().split("zone")[1].strip().split()[0]
                        )
                        zone = int("".join(c for c in zone_text if c.isdigit()))
                        if 49 <= zone <= 56:  # Australian zones
                            epsg = 28300 + zone
                            print(f"Found EPSG from PRJ file: {epsg}")
                    except (IndexError, ValueError):
                        pass
                elif "WGS84" in prj_content:
                    epsg = 4326

        # Second priority: Look for projection info in the DFN file
        if epsg is None:
            dfn_content_str = "".join(dfn_content)
            if "GDA94" in dfn_content_str and "zone" in dfn_content_str.lower():
                try:
                    zone_text = (
                        dfn_content_str.lower().split("zone")[1].strip().split()[0]
                    )
                    zone = int("".join(c for c in zone_text if c.isdigit()))
                    if 49 <= zone <= 56:  # Australian zones
                        epsg = 28300 + zone
                        print(f"Found EPSG from DFN content: {epsg}")
                except (IndexError, ValueError):
                    pass
            elif "WGS84" in dfn_content_str:
                epsg = 4326

        print(f"Final EPSG: {epsg}")
        return header_list, field_defs, epsg, data_format_type

    def parse_legacy_data(
        self, dat_file_path, header_list, field_defs, data_format_type
    ):

        # Read the DAT file data
        data = []

        with open(dat_file_path, "r") as dat_file:
            for line_num, line in enumerate(dat_file):
                original_line = line.strip()
                if not original_line:
                    continue

                # Split the line into values
                values = original_line.split()

                # Skip if not enough values
                if len(values) < len(header_list):
                    print(
                        f"Warning: Line {line_num} has fewer values ({len(values)}) than fields ({len(header_list)})"
                    )
                    continue

                # Map values to fields based on the format
                row_data = {}

                # Special handling for files with "DATA" as an identifier
                if data_format_type == "data_as_identifier" and values[0] == "DATA":
                    # Skip the "DATA" identifier and shift all values left by one
                    field_values = values[1:] if len(values) > 1 else []

                    # If we have fewer values than fields, skip this row
                    if len(field_values) < len(header_list):
                        print(
                            f"Warning: Line {line_num} has fewer values ({len(field_values)}) than fields ({len(header_list)}) after removing DATA"
                        )
                        continue

                    # Map each field to its value
                    for i, field_name in enumerate(header_list):
                        if i < len(field_values):
                            value = field_values[i]

                            # Convert to appropriate type
                            if i < len(field_defs):
                                field_type = field_defs[i].get("type", "")

                                if field_type in ["F", "D"]:  # Float
                                    try:
                                        value = float(value)
                                    except (ValueError, TypeError):
                                        pass
                                elif field_type == "I":  # Integer
                                    try:
                                        value = int(value)
                                    except (ValueError, TypeError):
                                        pass

                            row_data[field_name] = value
                else:
                    # For files with DATA prefix to remove
                    if original_line.startswith("DATA"):
                        line_content = original_line[4:].strip()
                        values = line_content.split()
                    else:
                        values = original_line.split()

                    # Map values to fields
                    for i, field_name in enumerate(header_list):
                        if i < len(values):
                            value = values[i]

                            # Convert to appropriate type
                            if i < len(field_defs):
                                field_type = field_defs[i].get("type", "")

                                if field_type in ["F", "D"]:  # Float
                                    try:
                                        value = float(value)
                                    except (ValueError, TypeError):
                                        pass
                                elif field_type == "I":  # Integer
                                    try:
                                        value = int(value)
                                    except (ValueError, TypeError):
                                        pass

                            row_data[field_name] = value

                # Add this row to our data
                data.append(row_data)

                # Print sample row for debugging
                if line_num < 2:
                    print(
                        f"Sample row {line_num}: First 5 items = {list(row_data.items())[:5]}"
                    )
                    print(f"Original line: {original_line}")
                    print(f"Processed values: {values[:10]}")  # Show first 10 values

        print(f"Processed {len(data)} data rows")

        return data
