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
        for line in lines:
            line = line.strip()
            if not line:
                continue

            # Skip lines that don't start with DEFN
            if not line.startswith("DEFN"):
                continue

            # Skip END DEFN lines
            if "END DEFN" in line:
                continue

            # Check different DEFN formats
            defn_pattern = r"DEFN\s*(\d+)?\s*ST=REC[ORD]*"
            defn_match = re.match(defn_pattern, line, re.IGNORECASE)

            if defn_match:
                # Extract field information
                parts = line.split(";")

                for part in parts[1:]:  # Skip the first part (DEFN...)
                    # Look for field:format pattern
                    field_pattern = r"([^:;]+):([^:;,]+)"
                    field_match = re.match(field_pattern, part.strip())

                    if field_match:
                        field_name = field_match.group(1).strip()
                        field_format = field_match.group(2).strip()

                        # Skip RT field
                        if field_name.upper() == "RT":
                            continue

                        # Check for NAME= attribute in the rest of the part
                        name_match = re.search(r"NAME=([^,;]+)", part)
                        display_name = (
                            name_match.group(1).strip() if name_match else field_name
                        )

                        # Check for array format (e.g., 256f5.0)
                        array_match = re.match(r"(\d+)([A-Za-z]+.*)", field_format)
                        if array_match:
                            count = int(array_match.group(1))
                            base_format = array_match.group(2)

                            # Handle large arrays as a single field with a special name
                            if (
                                count > 20
                            ):  # Consider arrays with more than 20 elements as "large"
                                print(
                                    f"Found large array field: {display_name} with {count} elements"
                                )

                                # Use the display name with _array suffix
                                array_field_name = f"{display_name}_array"
                                header_list.append(array_field_name)

                                # Extract field type from format, but use string for the array
                                field_type = (
                                    QVariant.String
                                )  # Store as string for arrays

                                # Extract NULL value if provided
                                null_match = re.search(r"NULL=([^,;]+)", part)
                                null_value = (
                                    null_match.group(1).strip() if null_match else None
                                )

                                # Extract units if provided
                                unit_match = re.search(r"UNIT=([^,;]+)", part)
                                unit = (
                                    unit_match.group(1).strip() if unit_match else None
                                )

                                # Store additional array info
                                field_defs[array_field_name] = {
                                    "type": field_type,
                                    "format": field_format,
                                    "null": null_value,
                                    "unit": unit,
                                    "original_name": field_name,
                                    "is_array": True,
                                    "array_count": count,
                                    "array_base_format": base_format,
                                }
                            else:
                                # Create multiple fields for smaller arrays - use shorter names for shapefile compatibility
                                # Truncate the base name to leave room for the index
                                base_name = (
                                    display_name[:6]
                                    if len(display_name) > 6
                                    else display_name
                                )

                                for i in range(count):
                                    # Create field name that fits in 10 characters
                                    if count > 999:  # Need more compact format
                                        array_field_name = f"{base_name[:4]}_{i:04d}"
                                    elif count > 99:
                                        array_field_name = f"{base_name[:5]}_{i:03d}"
                                    else:
                                        array_field_name = f"{base_name}_{i:02d}"

                                    header_list.append(array_field_name)

                                    # Extract field type from format
                                    field_type = QVariant.String
                                    if base_format.startswith("I"):
                                        field_type = QVariant.Int
                                    elif base_format.startswith(("F", "E", "D")):
                                        field_type = QVariant.Double

                                    field_defs[array_field_name] = {
                                        "type": field_type,
                                        "format": base_format,
                                        "null": None,
                                        "unit": None,
                                        "original_name": field_name,
                                        "is_array": False,  # Not treating as special array
                                    }
                        else:
                            # Regular single field
                            header_list.append(display_name)
                            print(
                                f"Added field to header list: {display_name} (original: {field_name})"
                            )

                            # Extract field type from format
                            field_type = QVariant.String
                            if field_format.startswith("I"):
                                field_type = QVariant.Int
                            elif field_format.startswith(("F", "E", "D")):
                                field_type = QVariant.Double
                            elif field_format.startswith("A"):
                                field_type = QVariant.String

                            # Extract NULL value if provided
                            null_match = re.search(r"NULL=([^,;]+)", part)
                            null_value = (
                                null_match.group(1).strip() if null_match else None
                            )

                            # Extract units if provided
                            unit_match = re.search(r"UNIT=([^,;]+)", part)
                            unit = unit_match.group(1).strip() if unit_match else None

                            # Check for EPSG code
                            epsg_match = re.search(
                                r"epsgcode=(\d+)", part, re.IGNORECASE
                            )
                            if epsg_match and display_name.lower() in [
                                "longitude",
                                "latitude",
                                "easting",
                                "northing",
                                "x_gda94",
                                "y_gda94",
                                "long",
                                "lat",
                                "x",
                                "y",
                            ]:
                                epsg = epsg_match.group(1)
                                print(
                                    f"Found EPSG code {epsg} for field {display_name}"
                                )

                            # Create field definition
                            field_defs[display_name] = {
                                "type": field_type,
                                "format": field_format,
                                "null": null_value,
                                "unit": unit,
                                "original_name": field_name,
                                "is_array": False,
                            }

                        # Only process the first field in this part
                        break

        # Check if we found fields
        if header_list:
            print(f"Found {len(header_list)} fields in DFN file:")
            for i, field in enumerate(header_list[:20]):  # Show first 20 fields
                field_def = field_defs.get(field, {})
                field_type = "String"
                if field_def.get("type") == QVariant.Int:
                    field_type = "Integer"
                elif field_def.get("type") == QVariant.Double:
                    field_type = "Double"

                # Show if it's an array field
                if field_def.get("is_array", False):
                    array_count = field_def.get("array_count", 0)
                    field_type = f"{field_type} Array[{array_count}]"

                print(f"  {i+1}. {field} ({field_type})")
            if len(header_list) > 20:
                print(f"  ... and {len(header_list) - 20} more fields")
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

        # Track stats for reporting
        processed_lines = 0
        skipped_lines = 0

        # Make a working copy of the header list and field definitions
        working_header_list = header_list.copy()
        working_field_defs = field_defs.copy()

        # Find array fields
        array_fields = {}
        for field_name in working_header_list:
            field_def = working_field_defs.get(field_name, {})
            if field_def.get("is_array", False):
                array_fields[field_name] = {
                    "count": field_def.get("array_count", 0),
                    "base_format": field_def.get("array_base_format", ""),
                }
                print(
                    f"Found array field {field_name} with {field_def.get('array_count', 0)} elements"
                )

        # Check if this file has COMMENTS field
        has_comments_field = "COMMENTS" in working_header_list

        # Check the structure of the first few lines to determine format
        if lines:
            # Sample the first few non-empty lines
            sample_lines = []
            for line in lines[:10]:
                line = line.strip()
                if line:
                    sample_lines.append(line)

            if sample_lines:
                # Check if last character of each line is consistently an asterisk (*)
                has_trailing_asterisk = all(line.endswith("*") for line in sample_lines)

                # Check field counts
                field_counts = [len(line.split()) for line in sample_lines]
                most_common_count = max(set(field_counts), key=field_counts.count)

                # Calculate expected field count in data file
                # Start with regular fields (non-array)
                regular_fields = [
                    f
                    for f in working_header_list
                    if f not in array_fields and f != "COMMENTS"
                ]
                expected_count = len(regular_fields)

                # Add count for array fields - each array takes its full number of elements in the data
                for field_name, array_info in array_fields.items():
                    expected_count += array_info["count"]

                # Handle comments field
                if has_comments_field and has_trailing_asterisk:
                    expected_count += 0  # No extra column for COMMENTS with asterisk
                elif has_comments_field and not has_trailing_asterisk:
                    expected_count += 1  # One column for COMMENTS field

                print(
                    f"Header has {len(working_header_list)} fields (expecting {expected_count} columns), data has {most_common_count} fields"
                )

                # Data format determination
                if has_trailing_asterisk:
                    print(f"Detected trailing asterisk (*) at the end of each line")

                    # If there's no COMMENTS field but we need one for the asterisk
                    if not has_comments_field:
                        print(
                            f"Adding COMMENTS field to header list for trailing asterisk"
                        )
                        working_header_list.append("COMMENTS")
                        working_field_defs["COMMENTS"] = {
                            "type": QVariant.String,
                            "format": "A1",
                            "null": None,
                            "unit": None,
                            "is_array": False,
                        }
                        has_comments_field = True

                    # Make sure COMMENTS is the last field
                    if working_header_list[-1] != "COMMENTS":
                        working_header_list.remove("COMMENTS")
                        working_header_list.append("COMMENTS")
                        print(f"Moved COMMENTS field to the end of header list")

                # Missing COMMENTS field in data
                elif has_comments_field and most_common_count < expected_count:
                    # The data doesn't have a COMMENTS field, so remove it from the header
                    if (
                        "COMMENTS" in working_header_list
                        and "COMMENTS" not in array_fields
                    ):
                        working_header_list.remove("COMMENTS")
                        print(
                            f"Removed COMMENTS field from header list as it's not in the data"
                        )
                        has_comments_field = False

        for line_num, line in enumerate(lines):
            line = line.strip()
            if not line:  # Skip empty lines
                skipped_lines += 1
                continue

            try:
                # Split line into fields
                fields = line.split()

                # Handle trailing asterisk if present
                has_trailing_asterisk_in_line = False
                if len(fields) > 0 and fields[-1] == "*":
                    has_trailing_asterisk_in_line = True
                    fields.pop()  # Remove the asterisk from fields

                # Create data row
                row = {}

                # Process normal fields first (non-array fields)
                current_field_index = 0  # Index in the data fields

                for field_name in working_header_list:
                    if field_name == "COMMENTS":
                        # Handle COMMENTS field separately
                        if has_trailing_asterisk_in_line:
                            row["COMMENTS"] = "*"
                        continue

                    field_def = working_field_defs.get(field_name, {})
                    if field_def.get("is_array", False):
                        # This is an array field
                        array_count = field_def.get("array_count", 0)

                        # Make sure we don't go out of bounds
                        if current_field_index + array_count <= len(fields):
                            # Collect array values
                            array_values = fields[
                                current_field_index : current_field_index + array_count
                            ]
                            row[field_name] = ";".join(array_values)
                            current_field_index += array_count
                        else:
                            # Not enough fields for this array
                            print(
                                f"Warning: Line {line_num+1} doesn't have enough fields for array {field_name}"
                            )
                            if current_field_index < len(fields):
                                # Use what we have
                                array_values = fields[current_field_index:]
                                row[field_name] = ";".join(array_values)
                                current_field_index = len(fields)
                            else:
                                row[field_name] = ""  # Empty array
                    else:
                        # Regular field
                        if current_field_index < len(fields):
                            value = fields[current_field_index]
                            current_field_index += 1

                            # Convert to appropriate type based on field definition
                            if field_def.get("null") and value == field_def.get("null"):
                                row[field_name] = None
                            elif field_def.get("type") == QVariant.Int:
                                try:
                                    row[field_name] = int(float(value))
                                except ValueError:
                                    row[field_name] = value
                            elif field_def.get("type") == QVariant.Double:
                                try:
                                    row[field_name] = float(value)
                                except ValueError:
                                    row[field_name] = value
                            else:
                                row[field_name] = value
                        else:
                            # Not enough fields in data
                            row[field_name] = None

                data.append(row)
                processed_lines += 1

                # Debug first few rows
                if processed_lines <= 3:
                    field_samples = []
                    for j, field_name in enumerate(working_header_list[:5]):
                        if field_name in row:
                            value = row[field_name]
                            # For array fields, show just the length or first few elements
                            if field_name in array_fields:
                                if value:
                                    elements = value.split(";")
                                    value = f"[Array with {len(elements)} elements]"
                                else:
                                    value = "[Empty array]"
                            field_samples.append(f"{field_name}={value}")
                    print(f"Row {processed_lines}: {', '.join(field_samples)}...")

            except Exception as e:
                print(f"Error parsing line {line_num+1}: {str(e)}")
                print(f"Line content: {line[:100]}...")
                skipped_lines += 1

        print(
            f"Parsed {len(data)} rows from DAT file (Processed: {processed_lines}, Skipped: {skipped_lines})"
        )

        # Always return just the data (list of dictionaries)
        return data

    def export_to_shapefile(
        self, data, field_defs, output_path, file_name, x_field, y_field, epsg_code=None
    ):
        """
        Export parsed data to shapefile

        Args:
            data (list): List of dictionaries containing the data
            header_list (list): List of field names
            field_defs (dict): Dictionary of field definitions
            output_path (str): Path for output shapefile
            x_field, y_field (str,str):   coordinate fields
            epsg_code (str): EPSG code for the coordinate system (e.g., "4326", "28356")
        """
        print(f"Exporting to shapefile: {output_path}")

        if output_path:
            if os.path.exists(output_path):
                try:
                    os.remove(output_path)
                    # print(f"File '{file_path}' has been deleted")
                    file_deleted = True
                except Exception as e:
                    print(f"Failed to delete '{output_path}': {str(e)}")

        # Ensure data is a list of dictionaries
        if not data or not isinstance(data[0], dict):
            raise ValueError("Data must be a list of dictionaries")

        if not x_field or not y_field:
            raise ValueError("Coordinate fields must be specified")

        # Check if coordinate fields exist in data
        first_row = data[0]
        if x_field not in first_row or y_field not in first_row:
            print(f"Available fields in data: {list(first_row.keys())}")
            raise ValueError(
                f"Coordinate fields '{x_field}' and/or '{y_field}' not found in data"
            )

        # Create CRS from EPSG code
        if epsg_code:
            crs = QgsCoordinateReferenceSystem(f"EPSG:{epsg_code}")
            if not crs.isValid():
                print(f"Warning: Invalid EPSG code {epsg_code}. Using default WGS84.")
                crs = QgsCoordinateReferenceSystem("EPSG:4326")
        else:
            print("No EPSG code provided. Using default WGS84.")
            crs = QgsCoordinateReferenceSystem("EPSG:4326")

        print(f"Using CRS: {crs.authid()} - {crs.description()}")

        # Create fields for the shapefile - use actual fields from data
        qgs_fields = QgsFields()
        field_name_mapping = {}  # Map original names to truncated names

        # Use the actual fields from the data
        for field_name in first_row.keys():
            field_def = field_defs.get(field_name, {})
            field_type = field_def.get("type", QVariant.String)

            # Ensure field name fits shapefile limits (10 characters)
            if len(field_name) > 10:
                # For array fields that were already shortened, keep them as is
                if "_" in field_name and field_name.count("_") == 1:
                    truncated_name = field_name
                else:
                    # Truncate other long field names
                    truncated_name = field_name[:10]

                # Handle duplicates by adding a number
                original_truncated = truncated_name
                counter = 1
                while truncated_name in field_name_mapping.values():
                    truncated_name = f"{original_truncated[:8]}{counter:02d}"
                    counter += 1
            else:
                truncated_name = field_name

            field_name_mapping[field_name] = truncated_name

            qgs_field = QgsField(truncated_name, field_type)
            qgs_fields.append(qgs_field)

        # Create the shapefile writer with the specified CRS
        writer = QgsVectorFileWriter(
            output_path,
            "UTF-8",
            qgs_fields,
            QgsWkbTypes.Point,
            crs,  # Use the CRS from the EPSG code
            "ESRI Shapefile",
        )

        if writer.hasError() != QgsVectorFileWriter.NoError:
            raise Exception(f"Error creating shapefile writer: {writer.errorMessage()}")

        # Write features
        features_written = 0
        errors = 0

        for row_num, row in enumerate(data):
            try:
                # Create feature
                feature = QgsFeature()
                feature.setFields(qgs_fields)

                # Set geometry
                x_val = row.get(x_field)
                y_val = row.get(y_field)

                if x_val is not None and y_val is not None:
                    try:
                        x_coord = float(x_val)
                        y_coord = float(y_val)
                        point = QgsPointXY(x_coord, y_coord)
                        feature.setGeometry(QgsGeometry.fromPointXY(point))
                    except (ValueError, TypeError) as e:
                        print(f"Row {row_num}: Invalid coordinates: {x_val}, {y_val}")
                        errors += 1
                        continue
                else:
                    print(f"Row {row_num}: Missing coordinates")
                    errors += 1
                    continue

                # Set attributes using the mapped field names
                for original_name, truncated_name in field_name_mapping.items():
                    if original_name in row and row[original_name] is not None:
                        value = row[original_name]

                        # Get field index by truncated name
                        field_idx = qgs_fields.indexFromName(truncated_name)
                        if field_idx >= 0:
                            # Convert value to appropriate type
                            field = qgs_fields.at(field_idx)
                            if field.type() == QVariant.Int:
                                try:
                                    value = int(float(value))
                                except (ValueError, TypeError):
                                    value = None
                            elif field.type() == QVariant.Double:
                                try:
                                    value = float(value)
                                except (ValueError, TypeError):
                                    value = None

                            feature.setAttribute(field_idx, value)

                # Add feature to shapefile
                writer.addFeature(feature)
                features_written += 1

            except Exception as e:
                print(f"Error writing row {row_num}: {str(e)}")
                errors += 1

        # Clean up
        del writer

        print(
            f"Shapefile export complete: {features_written} features written, {errors} errors"
        )
        # Load the shapefile
        layer = QgsVectorLayer(output_path, file_name, "ogr")

        if not layer.isValid():
            print(f"Layer is not valid! Path: {output_path}")
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
