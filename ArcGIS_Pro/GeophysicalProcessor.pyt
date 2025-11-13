# GeophysicsTools.pyt - Simplified version

import arcpy
import numpy as np
import os
import sys

# Global variable to store the processor class
GeophysicalProcessor = None


def get_processor_class():
    """Import GeophysicalProcessor with simplified error handling"""
    global GeophysicalProcessor

    if GeophysicalProcessor is not None:
        return GeophysicalProcessor

    try:
        # Get directories
        toolbox_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        calcs_dir = os.path.join(toolbox_dir, "calcs")
        worms_dir = os.path.join(toolbox_dir, "calcs/worms")

        # Verify the GeophysicalProcessor.py file exists
        gp_file = os.path.join(calcs_dir, "GeophysicalProcessor.py")
        if not os.path.exists(gp_file):
            arcpy.AddError(f"GeophysicalProcessor.py not found at: {gp_file}")
            return None

        # Add to Python path - put calcs_dir first to prioritize it
        if calcs_dir not in sys.path:
            sys.path.insert(0, calcs_dir)
        if worms_dir not in sys.path:
            sys.path.insert(0, worms_dir)
        if toolbox_dir not in sys.path:
            sys.path.insert(0, toolbox_dir)

        # Read and modify the file content to handle relative imports
        with open(gp_file, "r", encoding="utf-8") as f:
            content = f.read()

        # Handle multi-line imports from worms.Utility
        lines = content.split("\n")
        modified_lines = []
        skip_next = False

        for line in lines:
            if skip_next and ")" in line:
                skip_next = False
                modified_lines.append(f"# {line}")
                continue

            if skip_next:
                modified_lines.append(f"# {line}")
                continue

            if "from ..worms.Utility import (" in line:
                skip_next = True
                modified_lines.append(f"# {line}")
                continue

            # Comment out any other relative imports
            if line.strip().startswith("from .."):
                modified_lines.append(f"# {line}")
                continue

            modified_lines.append(line)

        modified_content = "\n".join(modified_lines)

        # Create execution namespace with required dependencies
        namespace = {
            "__name__": "GeophysicalProcessor",
            "__file__": gp_file,
            "np": np,
            "numpy": np,
            "os": os,
            "sys": sys,
            "ceil": __import__("math").ceil,
            "fabs": __import__("math").fabs,
        }

        # Add scipy imports
        try:
            from scipy import ndimage, fft

            namespace.update({"ndimage": ndimage, "fft": fft})
        except ImportError:
            arcpy.AddError("Scipy not available")
            return None

        # Add other standard imports
        import csv, glob

        namespace.update({"csv": csv, "glob": glob})

        # Try to add GDAL
        try:
            from osgeo import gdal, ogr, osr

            namespace.update({"gdal": gdal, "ogr": ogr, "osr": osr})
        except ImportError:
            arcpy.AddMessage("GDAL not available - some functions may not work")

        # Add dummy functions for missing worms functions
        def dummy_function(*args, **kwargs):
            pass

        def fill_nan_simple(data):
            """Simple NaN filling for missing fill_nan function"""
            if not np.any(np.isnan(data)):
                return data
            filled = np.copy(data)
            # Simple linear interpolation along rows
            for i in range(data.shape[0]):
                row = filled[i, :]
                valid = ~np.isnan(row)
                if np.any(valid) and not np.all(valid):
                    indices = np.arange(len(row))
                    filled[i, ~valid] = np.interp(
                        indices[~valid], indices[valid], row[valid]
                    )
            # Simple linear interpolation along columns
            for j in range(data.shape[1]):
                col = filled[:, j]
                valid = ~np.isnan(col)
                if np.any(valid) and not np.all(valid):
                    indices = np.arange(len(col))
                    filled[~valid, j] = np.interp(
                        indices[~valid], indices[valid], col[valid]
                    )
            return filled

        # Add dummy implementations for worms utilities
        namespace.update(
            {
                "Wormer": type("DummyWormer", (), {"__init__": lambda self: None}),
                "GetExtent": lambda *args: [[0, 0], [1, 1], [0, 1]],
                "loadGrid": dummy_function,
                "numpy_array_to_raster": dummy_function,
                "insert_text_before_extension": lambda path, text: os.path.splitext(
                    path
                )[0]
                + text
                + os.path.splitext(path)[1],
                "fill_nan": fill_nan_simple,
            }
        )

        # Execute the modified content
        exec(modified_content, namespace)

        # Get the class from the namespace
        GeophysicalProcessor = namespace.get("GeophysicalProcessor")
        if GeophysicalProcessor is None:
            arcpy.AddError("GeophysicalProcessor class not found in namespace")
            return None

        arcpy.AddMessage("✓ GeophysicalProcessor imported successfully")
        return GeophysicalProcessor

    except ImportError as e:
        arcpy.AddError(f"Failed to import GeophysicalProcessor: {e}")
        arcpy.AddError("Make sure the calcs directory contains GeophysicalProcessor.py")
        return None
    except Exception as e:
        arcpy.AddError(f"Unexpected error importing GeophysicalProcessor: {e}")
        import traceback

        arcpy.AddError(traceback.format_exc())
        return None


# Processor wrapper
class ArcGISProcessor:
    def process_raster_tool(
        self, input_raster, output_raster, processing_method, **kwargs
    ):
        ProcessorClass = get_processor_class()

        if ProcessorClass is None:
            arcpy.AddError("GeophysicalProcessor not available")
            return False

        try:
            # Get raster properties
            desc = arcpy.Describe(input_raster)
            extent = desc.extent
            dx = desc.meanCellWidth
            dy = desc.meanCellHeight
            spatial_ref = desc.spatialReference

            # Read raster
            array = arcpy.RasterToNumPyArray(input_raster, nodata_to_value=np.nan)

            # Initialize processor - FIXED: don't remove buffer_size from kwargs
            init_buffer_size = kwargs.get("buffer_size", 100)  # Get without removing
            convert_to_degrees = kwargs.pop(
                "convert_to_degrees", False
            )  # This one can be removed

            processor = ProcessorClass(dx, dy, init_buffer_size)

            # Check if method exists
            if not hasattr(processor, processing_method):
                available = [
                    m
                    for m in dir(processor)
                    if not m.startswith("_") and callable(getattr(processor, m))
                ]
                arcpy.AddError(
                    f"Method '{processing_method}' not found. Available: {available[:5]}"
                )
                return False

            # Call the method - buffer_size is still in kwargs!
            method = getattr(processor, processing_method)
            result = method(array, **kwargs)

            # Convert to degrees if needed
            if convert_to_degrees:
                result = np.degrees(result)

            # Save result with proper spatial reference
            result_clean = np.where(np.isnan(result), -9999, result)
            lower_left = arcpy.Point(extent.XMin, extent.YMin)
            actual_dy = dy if dy > 0 else -dy

            raster = arcpy.NumPyArrayToRaster(
                result_clean, lower_left, dx, actual_dy, -9999
            )
            arcpy.DefineProjection_management(raster, spatial_ref)
            raster.save(output_raster)

            arcpy.AddMessage(f"✓ {processing_method} completed successfully")
            return True

        except Exception as e:
            arcpy.AddError(f"Processing error: {e}")
            import traceback

            arcpy.AddError(traceback.format_exc())
            return False


# Create the processor wrapper instance
processor_wrapper = ArcGISProcessor()


class Toolbox(object):
    def __init__(self):
        self.label = "Geophysical Processing Tools"
        self.alias = "geophysics"
        self.tools = [
            UpwardContinuation,
            DownwardContinuation,
            VerticalIntegration,
            AnalyticSignal,
            TiltAngle,
            ReductionToPole,
            BandPassFilter,
            HighPassFilter,
            LowPassFilter,
            DirectionalButterworthBandPass,
            RemoveRegionalTrend,
            ComputeDerivative,
            TotalHorizontalGradient,
            AutomaticGainControl,
        ]


class UpwardContinuation(object):
    def __init__(self):
        self.label = "Upward Continuation"
        self.description = "Upward continue to attenuate high frequencies"
        self.canRunInBackground = False

    def getParameterInfo(self):
        params = []

        params.append(
            arcpy.Parameter(
                displayName="Input Raster",
                name="input_raster",
                datatype="GPRasterLayer",
                parameterType="Required",
                direction="Input",
            )
        )

        params.append(
            arcpy.Parameter(
                displayName="Continuation Height (meters)",
                name="height",
                datatype="GPDouble",
                parameterType="Required",
                direction="Input",
            )
        )
        params[-1].value = 500.0

        params.append(
            arcpy.Parameter(
                displayName="Buffer Size (pixels)",
                name="buffer_size",
                datatype="GPLong",
                parameterType="Optional",
                direction="Input",
            )
        )
        params[-1].value = 5000

        params.append(
            arcpy.Parameter(
                displayName="Output Continued Raster",
                name="output_raster",
                datatype="DERasterDataset",
                parameterType="Required",
                direction="Output",
            )
        )

        return params

    def execute(self, parameters, messages):
        input_raster = parameters[0].valueAsText
        height = parameters[1].value
        buffer_size = parameters[2].value or 5000
        output_raster = parameters[3].valueAsText

        success = processor_wrapper.process_raster_tool(
            input_raster=input_raster,
            output_raster=output_raster,
            processing_method="upward_continuation",
            height=height,
            buffer_size=buffer_size,
        )

        if not success:
            raise arcpy.ExecuteError("Upward Continuation band pass failed")


class RemoveRegionalTrend(object):
    def __init__(self):
        self.label = "Remove Regional Trend (Polynomial)"
        self.description = "Remove regional trends using polynomial surface fitting"
        self.canRunInBackground = False

    def getParameterInfo(self):
        params = []

        params.append(
            arcpy.Parameter(
                displayName="Input Raster",
                name="input_raster",
                datatype="GPRasterLayer",
                parameterType="Required",
                direction="Input",
            )
        )

        params.append(
            arcpy.Parameter(
                displayName="Polynomial Order",
                name="order",
                datatype="GPString",
                parameterType="Required",
                direction="Input",
            )
        )
        params[-1].filter.type = "ValueList"
        params[-1].filter.list = ["1st Order (Plane)", "2nd Order (Quadratic)"]
        params[-1].value = "2nd Order (Quadratic)"
        params[-1].description = (
            "1st order removes linear gradients, 2nd order removes curved regional trends"
        )

        params.append(
            arcpy.Parameter(
                displayName="Output Residual Raster",
                name="output_raster",
                datatype="DERasterDataset",
                parameterType="Required",
                direction="Output",
            )
        )

        return params

    def execute(self, parameters, messages):
        input_raster = parameters[0].valueAsText
        order = parameters[1].valueAsText
        output_raster = parameters[2].valueAsText

        # Determine which method to use
        if "1st" in order:
            method_name = "remove_gradient"
            order_description = "1st order (plane fitting)"
        else:
            method_name = "remove_2o_gradient"
            order_description = "2nd order (quadratic surface fitting)"

        arcpy.AddMessage("=== Remove Regional Trend (Polynomial) ===")
        arcpy.AddMessage(f"Method: {order_description}")
        arcpy.AddMessage("Purpose: Regional/residual separation")
        arcpy.AddMessage("Output: Residual anomalies (regional trend removed)")

        # Use custom processing since these methods need mask parameter
        success = self._process_gradient_removal(
            input_raster, output_raster, method_name
        )

        if not success:
            raise arcpy.ExecuteError("Regional trend removal failed")

    def _process_gradient_removal(self, input_raster, output_raster, method_name):
        """Custom processing for gradient removal methods that need mask"""
        try:
            # Get processor class
            ProcessorClass = get_processor_class()
            if ProcessorClass is None:
                arcpy.AddError("GeophysicalProcessor not available")
                return False

            # Get raster properties
            desc = arcpy.Describe(input_raster)
            extent = desc.extent
            dx = desc.meanCellWidth
            dy = desc.meanCellHeight
            spatial_ref = desc.spatialReference

            # Read raster
            array = arcpy.RasterToNumPyArray(input_raster, nodata_to_value=np.nan)
            arcpy.AddMessage(f"Processing array of shape: {array.shape}")

            # Create mask for NaN values
            mask = np.isnan(array)

            # Initialize processor
            processor = ProcessorClass(
                dx, dy, 100
            )  # Buffer size not used for these methods

            # Call the appropriate method
            if hasattr(processor, method_name):
                method = getattr(processor, method_name)
                result = method(array, mask)
            else:
                arcpy.AddError(
                    f"Method {method_name} not found in GeophysicalProcessor"
                )
                return False

            # Save result
            result_clean = np.where(np.isnan(result), -9999, result)
            lower_left = arcpy.Point(extent.XMin, extent.YMin)
            actual_dy = dy if dy > 0 else -dy

            raster = arcpy.NumPyArrayToRaster(
                result_clean, lower_left, dx, actual_dy, -9999
            )
            arcpy.DefineProjection_management(raster, spatial_ref)
            raster.save(output_raster)

            arcpy.AddMessage(f"✓ {method_name} completed successfully")
            return True

        except Exception as e:
            arcpy.AddError(f"Processing error: {e}")
            import traceback

            arcpy.AddError(traceback.format_exc())
            return False


class ComputeDerivative(object):
    def __init__(self):
        self.label = "Compute Derivative"
        self.description = "Calculate derivatives in x, y, or z directions"
        self.canRunInBackground = False

    def getParameterInfo(self):
        params = []

        params.append(
            arcpy.Parameter(
                displayName="Input Raster",
                name="input_raster",
                datatype="GPRasterLayer",
                parameterType="Required",
                direction="Input",
            )
        )

        params.append(
            arcpy.Parameter(
                displayName="Direction",
                name="direction",
                datatype="GPString",
                parameterType="Required",
                direction="Input",
            )
        )
        params[-1].filter.type = "ValueList"
        params[-1].filter.list = ["x", "y", "z"]
        params[-1].value = "z"

        params.append(
            arcpy.Parameter(
                displayName="Order",
                name="order",
                datatype="GPLong",
                parameterType="Optional",
                direction="Input",
            )
        )
        params[-1].value = 1

        params.append(
            arcpy.Parameter(
                displayName="Buffer Size (pixels)",
                name="buffer_size",
                datatype="GPLong",
                parameterType="Optional",
                direction="Input",
            )
        )
        params[-1].value = 5000

        params.append(
            arcpy.Parameter(
                displayName="Output Derivative Raster",
                name="output_raster",
                datatype="DERasterDataset",
                parameterType="Required",
                direction="Output",
            )
        )

        return params

    def execute(self, parameters, messages):
        input_raster = parameters[0].valueAsText
        direction = parameters[1].valueAsText
        order = parameters[2].value or 1
        buffer_size = parameters[3].value or 5000
        output_raster = parameters[4].valueAsText

        success = processor_wrapper.process_raster_tool(
            input_raster=input_raster,
            output_raster=output_raster,
            processing_method="compute_derivative",
            direction=direction,
            order=order,
            buffer_size=buffer_size,
        )

        if not success:
            raise arcpy.ExecuteError("Derivative calculation failed")


class TotalHorizontalGradient(object):
    def __init__(self):
        self.label = "Total Horizontal Gradient"
        self.description = "Calculate total horizontal gradient magnitude"
        self.canRunInBackground = False

    def getParameterInfo(self):
        params = []

        params.append(
            arcpy.Parameter(
                displayName="Input Raster",
                name="input_raster",
                datatype="GPRasterLayer",
                parameterType="Required",
                direction="Input",
            )
        )

        params.append(
            arcpy.Parameter(
                displayName="Buffer Size (pixels)",
                name="buffer_size",
                datatype="GPLong",
                parameterType="Optional",
                direction="Input",
            )
        )
        params[-1].value = 5000

        params.append(
            arcpy.Parameter(
                displayName="Output THG Raster",
                name="output_raster",
                datatype="DERasterDataset",
                parameterType="Required",
                direction="Output",
            )
        )

        return params

    def execute(self, parameters, messages):
        input_raster = parameters[0].valueAsText
        buffer_size = parameters[1].value or 5000
        output_raster = parameters[2].valueAsText

        success = processor_wrapper.process_raster_tool(
            input_raster=input_raster,
            output_raster=output_raster,
            processing_method="total_horizontal_gradient",
            buffer_size=buffer_size,
        )

        if not success:
            raise arcpy.ExecuteError("Total horizontal gradient calculation failed")


class AutomaticGainControl(object):
    def __init__(self):
        self.label = "Automatic Gain Control"
        self.description = "Apply AGC normalization to highlight low amplitude features"
        self.canRunInBackground = False

    def getParameterInfo(self):
        params = []

        params.append(
            arcpy.Parameter(
                displayName="Input Raster",
                name="input_raster",
                datatype="GPRasterLayer",
                parameterType="Required",
                direction="Input",
            )
        )

        params.append(
            arcpy.Parameter(
                displayName="Window Size (pixels)",
                name="window_size",
                datatype="GPLong",
                parameterType="Required",
                direction="Input",
            )
        )
        params[-1].value = 10

        params.append(
            arcpy.Parameter(
                displayName="Output AGC Raster",
                name="output_raster",
                datatype="DERasterDataset",
                parameterType="Required",
                direction="Output",
            )
        )

        return params

    def execute(self, parameters, messages):
        input_raster = parameters[0].valueAsText
        window_size = parameters[1].value
        output_raster = parameters[2].valueAsText

        success = processor_wrapper.process_raster_tool(
            input_raster=input_raster,
            output_raster=output_raster,
            processing_method="automatic_gain_control",
            window_size=window_size,
        )

        if not success:
            raise arcpy.ExecuteError("AGC failed")[0].valueAsText


class DownwardContinuation(object):
    def __init__(self):
        self.label = "Downward Continuation"
        self.description = (
            "Enhance high-frequency signals (use with caution - amplifies noise)"
        )
        self.canRunInBackground = False

    def getParameterInfo(self):
        params = []

        params.append(
            arcpy.Parameter(
                displayName="Input Raster",
                name="input_raster",
                datatype="GPRasterLayer",
                parameterType="Required",
                direction="Input",
            )
        )

        params.append(
            arcpy.Parameter(
                displayName="Continuation Height (meters)",
                name="height",
                datatype="GPDouble",
                parameterType="Required",
                direction="Input",
            )
        )
        params[-1].value = 100.0

        params.append(
            arcpy.Parameter(
                displayName="Buffer Size (pixels)",
                name="buffer_size",
                datatype="GPLong",
                parameterType="Optional",
                direction="Input",
            )
        )
        params[-1].value = 5000

        params.append(
            arcpy.Parameter(
                displayName="Output Continued Raster",
                name="output_raster",
                datatype="DERasterDataset",
                parameterType="Required",
                direction="Output",
            )
        )

        return params

    def execute(self, parameters, messages):
        input_raster = parameters[0].valueAsText
        height = parameters[1].value
        buffer_size = parameters[2].value or 5000
        output_raster = parameters[3].valueAsText

        arcpy.AddMessage("⚠️  WARNING: Downward continuation amplifies noise")
        arcpy.AddMessage(f"Downward continuation height: {height} meters")

        success = processor_wrapper.process_raster_tool(
            input_raster=input_raster,
            output_raster=output_raster,
            processing_method="downward_continuation",
            height=height,
            buffer_size=buffer_size,
        )

        if not success:
            raise arcpy.ExecuteError("Downward continuation failed")


class VerticalIntegration(object):
    def __init__(self):
        self.label = "Vertical Integration"
        self.description = (
            "Integrate field vertically to simulate deeper equivalent sources"
        )
        self.canRunInBackground = False

    def getParameterInfo(self):
        params = []

        params.append(
            arcpy.Parameter(
                displayName="Input Raster",
                name="input_raster",
                datatype="GPRasterLayer",
                parameterType="Required",
                direction="Input",
            )
        )

        params.append(
            arcpy.Parameter(
                displayName="Maximum Wavenumber",
                name="max_wavenumber",
                datatype="GPDouble",
                parameterType="Optional",
                direction="Input",
            )
        )
        params[-1].value = 0.1

        params.append(
            arcpy.Parameter(
                displayName="Minimum Wavenumber",
                name="min_wavenumber",
                datatype="GPDouble",
                parameterType="Optional",
                direction="Input",
            )
        )
        params[-1].value = 0.001

        params.append(
            arcpy.Parameter(
                displayName="Buffer Size (pixels)",
                name="buffer_size",
                datatype="GPLong",
                parameterType="Optional",
                direction="Input",
            )
        )
        params[-1].value = 5000

        params.append(
            arcpy.Parameter(
                displayName="Output Integrated Raster",
                name="output_raster",
                datatype="DERasterDataset",
                parameterType="Required",
                direction="Output",
            )
        )

        return params

    def execute(self, parameters, messages):
        input_raster = parameters[0].valueAsText
        max_wavenumber = parameters[1].value
        min_wavenumber = parameters[2].value or 0.001
        buffer_size = parameters[3].value or 5000
        output_raster = parameters[4].valueAsText

        success = processor_wrapper.process_raster_tool(
            input_raster=input_raster,
            output_raster=output_raster,
            processing_method="vertical_integration",
            max_wavenumber=max_wavenumber,
            min_wavenumber=min_wavenumber,
            buffer_size=buffer_size,
        )

        if not success:
            raise arcpy.ExecuteError("Vertical integration failed")


class AnalyticSignal(object):
    def __init__(self):
        self.label = "Analytic Signal"
        self.description = "Calculate analytic signal amplitude"
        self.canRunInBackground = False

    def getParameterInfo(self):
        params = []

        params.append(
            arcpy.Parameter(
                displayName="Input Raster",
                name="input_raster",
                datatype="GPRasterLayer",
                parameterType="Required",
                direction="Input",
            )
        )

        params.append(
            arcpy.Parameter(
                displayName="Buffer Size (pixels)",
                name="buffer_size",
                datatype="GPLong",
                parameterType="Optional",
                direction="Input",
            )
        )
        params[-1].value = 5000

        params.append(
            arcpy.Parameter(
                displayName="Output Analytic Signal Raster",
                name="output_raster",
                datatype="DERasterDataset",
                parameterType="Required",
                direction="Output",
            )
        )

        return params

    def execute(self, parameters, messages):
        input_raster = parameters[0].valueAsText
        buffer_size = parameters[1].value or 5000
        output_raster = parameters[2].valueAsText

        success = processor_wrapper.process_raster_tool(
            input_raster=input_raster,
            output_raster=output_raster,
            processing_method="analytic_signal",
            buffer_size=buffer_size,
        )

        if not success:
            raise arcpy.ExecuteError("Analytic signal calculation failed")


class TiltAngle(object):
    def __init__(self):
        self.label = "Tilt Angle"
        self.description = "Calculate tilt angle derivative (in degrees)"
        self.canRunInBackground = False

    def getParameterInfo(self):
        params = []

        params.append(
            arcpy.Parameter(
                displayName="Input Raster",
                name="input_raster",
                datatype="GPRasterLayer",
                parameterType="Required",
                direction="Input",
            )
        )

        params.append(
            arcpy.Parameter(
                displayName="Buffer Size (pixels)",
                name="buffer_size",
                datatype="GPLong",
                parameterType="Optional",
                direction="Input",
            )
        )
        params[-1].value = 5000

        params.append(
            arcpy.Parameter(
                displayName="Output Tilt Angle Raster (degrees)",
                name="output_raster",
                datatype="DERasterDataset",
                parameterType="Required",
                direction="Output",
            )
        )

        return params

    def execute(self, parameters, messages):
        input_raster = parameters[0].valueAsText
        buffer_size = parameters[1].value or 5000
        output_raster = parameters[2].valueAsText

        success = processor_wrapper.process_raster_tool(
            input_raster=input_raster,
            output_raster=output_raster,
            processing_method="tilt_angle",
            buffer_size=buffer_size,
            convert_to_degrees=True,
        )

        if not success:
            raise arcpy.ExecuteError("Tilt angle calculation failed")


class ReductionToPole(object):
    def __init__(self):
        self.label = "Reduction to Pole"
        self.description = "Transform magnetic data to the pole"
        self.canRunInBackground = False

    def getParameterInfo(self):
        params = []

        params.append(
            arcpy.Parameter(
                displayName="Input Magnetic Raster",
                name="input_raster",
                datatype="GPRasterLayer",
                parameterType="Required",
                direction="Input",
            )
        )

        params.append(
            arcpy.Parameter(
                displayName="Magnetic Inclination (degrees)",
                name="inclination",
                datatype="GPDouble",
                parameterType="Required",
                direction="Input",
            )
        )
        params[-1].value = -60.0

        params.append(
            arcpy.Parameter(
                displayName="Magnetic Declination (degrees)",
                name="declination",
                datatype="GPDouble",
                parameterType="Required",
                direction="Input",
            )
        )
        params[-1].value = 1.0

        params.append(
            arcpy.Parameter(
                displayName="Buffer Size (pixels)",
                name="buffer_size",
                datatype="GPLong",
                parameterType="Optional",
                direction="Input",
            )
        )
        params[-1].value = 5000

        params.append(
            arcpy.Parameter(
                displayName="Output RTP Raster",
                name="output_raster",
                datatype="DERasterDataset",
                parameterType="Required",
                direction="Output",
            )
        )

        return params

    def execute(self, parameters, messages):
        input_raster = parameters[0].valueAsText
        inclination = parameters[1].value
        declination = parameters[2].value
        buffer_size = parameters[3].value or 5000
        output_raster = parameters[4].valueAsText

        success = processor_wrapper.process_raster_tool(
            input_raster=input_raster,
            output_raster=output_raster,
            processing_method="reduction_to_pole",
            inclination=inclination,
            declination=declination,
            buffer_size=buffer_size,
        )

        if not success:
            raise arcpy.ExecuteError("Reduction to pole failed")


class BandPassFilter(object):
    def __init__(self):
        self.label = "Band Pass Filter"
        self.description = "Apply FFT band pass filter with smooth transitions"
        self.canRunInBackground = False

    def getParameterInfo(self):
        params = []

        params.append(
            arcpy.Parameter(
                displayName="Input Raster",
                name="input_raster",
                datatype="GPRasterLayer",
                parameterType="Required",
                direction="Input",
            )
        )

        params.append(
            arcpy.Parameter(
                displayName="Low Cutoff Wavelength (map units)",
                name="low_cut",
                datatype="GPDouble",
                parameterType="Required",
                direction="Input",
            )
        )
        params[-1].value = 50000.0

        params.append(
            arcpy.Parameter(
                displayName="High Cutoff Wavelength (map units)",
                name="high_cut",
                datatype="GPDouble",
                parameterType="Required",
                direction="Input",
            )
        )
        params[-1].value = 5000.0

        params.append(
            arcpy.Parameter(
                displayName="High Transition Width (map units)",
                name="high_transition_width",
                datatype="GPDouble",
                parameterType="Optional",
                direction="Input",
            )
        )

        params.append(
            arcpy.Parameter(
                displayName="Low Transition Width (map units)",
                name="low_transition_width",
                datatype="GPDouble",
                parameterType="Optional",
                direction="Input",
            )
        )

        params.append(
            arcpy.Parameter(
                displayName="Buffer Size (pixels)",
                name="buffer_size",
                datatype="GPLong",
                parameterType="Optional",
                direction="Input",
            )
        )
        params[-1].value = 5000

        params.append(
            arcpy.Parameter(
                displayName="Output Filtered Raster",
                name="output_raster",
                datatype="DERasterDataset",
                parameterType="Required",
                direction="Output",
            )
        )

        return params

    def execute(self, parameters, messages):
        input_raster = parameters[0].valueAsText
        low_cut = parameters[1].value
        high_cut = parameters[2].value
        high_transition_width = parameters[3].value
        low_transition_width = parameters[4].value
        buffer_size = parameters[5].value or 5000
        output_raster = parameters[6].valueAsText

        kwargs = {"low_cut": low_cut, "high_cut": high_cut, "buffer_size": buffer_size}
        if high_transition_width:
            kwargs["high_transition_width"] = high_transition_width
        if low_transition_width:
            kwargs["low_transition_width"] = low_transition_width

        success = processor_wrapper.process_raster_tool(
            input_raster=input_raster,
            output_raster=output_raster,
            processing_method="band_pass_filter",
            **kwargs,
        )

        if not success:
            raise arcpy.ExecuteError("Band pass filter failed")


class HighPassFilter(object):
    def __init__(self):
        self.label = "High Pass Filter"
        self.description = "Remove long wavelength regional trends"
        self.canRunInBackground = False

    def getParameterInfo(self):
        params = []

        params.append(
            arcpy.Parameter(
                displayName="Input Raster",
                name="input_raster",
                datatype="GPRasterLayer",
                parameterType="Required",
                direction="Input",
            )
        )

        params.append(
            arcpy.Parameter(
                displayName="Cutoff Wavelength (map units)",
                name="cutoff_wavelength",
                datatype="GPDouble",
                parameterType="Required",
                direction="Input",
            )
        )
        params[-1].value = 20000.0

        params.append(
            arcpy.Parameter(
                displayName="Transition Width (map units)",
                name="transition_width",
                datatype="GPDouble",
                parameterType="Optional",
                direction="Input",
            )
        )
        params[-1].value = 20000.0

        params.append(
            arcpy.Parameter(
                displayName="Buffer Size (pixels)",
                name="buffer_size",
                datatype="GPLong",
                parameterType="Optional",
                direction="Input",
            )
        )
        params[-1].value = 5000

        params.append(
            arcpy.Parameter(
                displayName="Output High Pass Raster",
                name="output_raster",
                datatype="DERasterDataset",
                parameterType="Required",
                direction="Output",
            )
        )

        return params

    def execute(self, parameters, messages):
        input_raster = parameters[0].valueAsText
        cutoff_wavelength = parameters[1].value
        transition_width = parameters[2].value or 20000.0
        buffer_size = parameters[3].value or 5000
        output_raster = parameters[4].valueAsText

        success = processor_wrapper.process_raster_tool(
            input_raster=input_raster,
            output_raster=output_raster,
            processing_method="high_pass_filter",
            cutoff_wavelength=cutoff_wavelength,
            transition_width=transition_width,
            buffer_size=buffer_size,
        )

        if not success:
            raise arcpy.ExecuteError("High pass filter failed")


class LowPassFilter(object):
    def __init__(self):
        self.label = "Low Pass Filter"
        self.description = "Remove short wavelength noise"
        self.canRunInBackground = False

    def getParameterInfo(self):
        params = []

        params.append(
            arcpy.Parameter(
                displayName="Input Raster",
                name="input_raster",
                datatype="GPRasterLayer",
                parameterType="Required",
                direction="Input",
            )
        )

        params.append(
            arcpy.Parameter(
                displayName="Cutoff Wavelength (map units)",
                name="cutoff_wavelength",
                datatype="GPDouble",
                parameterType="Required",
                direction="Input",
            )
        )
        params[-1].value = 5000.0

        params.append(
            arcpy.Parameter(
                displayName="Transition Width (map units)",
                name="transition_width",
                datatype="GPDouble",
                parameterType="Optional",
                direction="Input",
            )
        )

        params.append(
            arcpy.Parameter(
                displayName="Buffer Size (pixels)",
                name="buffer_size",
                datatype="GPLong",
                parameterType="Optional",
                direction="Input",
            )
        )
        params[-1].value = 5000

        params.append(
            arcpy.Parameter(
                displayName="Output Low Pass Raster",
                name="output_raster",
                datatype="DERasterDataset",
                parameterType="Required",
                direction="Output",
            )
        )

        return params

    def execute(self, parameters, messages):
        input_raster = parameters[0].valueAsText
        cutoff_wavelength = parameters[1].value
        transition_width = parameters[2].value
        buffer_size = parameters[3].value or 5000
        output_raster = parameters[4].valueAsText

        kwargs = {"cutoff_wavelength": cutoff_wavelength, "buffer_size": buffer_size}
        if transition_width:
            kwargs["transition_width"] = transition_width

        success = processor_wrapper.process_raster_tool(
            input_raster=input_raster,
            output_raster=output_raster,
            processing_method="low_pass_filter",
            **kwargs,
        )

        if not success:
            raise arcpy.ExecuteError("Low pass filter failed")


class DirectionalButterworthBandPass(object):
    def __init__(self):
        self.label = "Directional Butterworth Band Pass"
        self.description = "Combined directional and band-pass filtering"
        self.canRunInBackground = False

    def getParameterInfo(self):
        params = []

        params.append(
            arcpy.Parameter(
                displayName="Input Raster",
                name="input_raster",
                datatype="GPRasterLayer",
                parameterType="Required",
                direction="Input",
            )
        )

        params.append(
            arcpy.Parameter(
                displayName="Low Cutoff Wavelength (map units)",
                name="low_cut",
                datatype="GPDouble",
                parameterType="Required",
                direction="Input",
            )
        )
        params[-1].value = 50000.0

        params.append(
            arcpy.Parameter(
                displayName="High Cutoff Wavelength (map units)",
                name="high_cut",
                datatype="GPDouble",
                parameterType="Required",
                direction="Input",
            )
        )
        params[-1].value = 5000.0

        params.append(
            arcpy.Parameter(
                displayName="Direction Angle (degrees from North)",
                name="direction_angle",
                datatype="GPDouble",
                parameterType="Required",
                direction="Input",
            )
        )
        params[-1].value = 0.0

        params.append(
            arcpy.Parameter(
                displayName="Direction Width (degrees)",
                name="direction_width",
                datatype="GPDouble",
                parameterType="Required",
                direction="Input",
            )
        )
        params[-1].value = 45.0

        params.append(
            arcpy.Parameter(
                displayName="Butterworth Order",
                name="order",
                datatype="GPLong",
                parameterType="Optional",
                direction="Input",
            )
        )
        params[-1].value = 4

        params.append(
            arcpy.Parameter(
                displayName="Buffer Size (pixels)",
                name="buffer_size",
                datatype="GPLong",
                parameterType="Optional",
                direction="Input",
            )
        )
        params[-1].value = 5000

        params.append(
            arcpy.Parameter(
                displayName="Output Filtered Raster",
                name="output_raster",
                datatype="DERasterDataset",
                parameterType="Required",
                direction="Output",
            )
        )

        return params

    def execute(self, parameters, messages):
        input_raster = parameters
