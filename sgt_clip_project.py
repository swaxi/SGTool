"""
Clip every layer of the current QGIS project to a polygon, and save a copy of the project that
loads the clipped layers.

Based on the "clip layers to polygon (with styles)" script. For every layer in the project's layer
tree the clipped data is written to the output folder as <name>_clip (GeoPackage for vector layers,
GeoTIFF for rasters) with a <name>_clip.qml holding the layer's style, and a <project>_clip.qgz is
written there with the same groups, order, visibility and styles, pointing at the new files.

Layers that cannot be clipped (web services such as WMS / XYZ tiles, tables without geometry, mesh
and point-cloud layers) are carried over unchanged. A layer that fails to clip is reported and left
out of the new project.
"""

import os

from qgis.core import (
    Qgis,
    QgsLayerTree,
    QgsProcessingContext,
    QgsProcessingFeedback,
    QgsProject,
    QgsRasterLayer,
    QgsReferencedRectangle,
    QgsVectorLayer,
    QgsCoordinateTransform,
)
import processing

SUFFIX = "_clip"


def safe_filename(name):
    """A layer name made safe for use in a file name."""
    keep = (" ", ".", "_", "-")
    text = "".join(c if (c.isalnum() or c in keep) else "_" for c in name)
    return text.strip().replace(" ", "_") or "layer"


def _unique(stem, used):
    candidate, n = stem, 2
    while candidate.lower() in used:
        candidate = f"{stem}_{n}"
        n += 1
    used.add(candidate.lower())
    return candidate


def _geometry_type_is_null(layer):
    null = Qgis.GeometryType.Null if hasattr(Qgis, "GeometryType") else 4
    return layer.geometryType() == null


def _is_float_raster(layer):
    try:
        dtype = layer.dataProvider().dataType(1)
        return "Float" in getattr(dtype, "name", str(dtype))
    except Exception:
        return False


class ClipJob:
    """One run of clipping the project's layers to a polygon layer."""

    def __init__(self, project, clip_layer, out_dir, feedback=None, progress=None,
                 skip_names=("_preview",)):
        """progress(done, total, layer_name) is called before each layer and returns False to
        stop. skip_names: layers whose names end with one of these (the live preview layer) are
        left out."""
        self.project = project
        self.clip_layer = clip_layer
        self.out_dir = out_dir
        self.feedback = feedback or QgsProcessingFeedback()
        self.progress = progress
        self.skip_names = tuple(skip_names)
        self.context = QgsProcessingContext()
        self.context.setProject(project)
        self.log = []                 # (level, text)
        self.clipped = []             # (layer name, path written)
        self.kept = []                # (layer name, why it was carried over unchanged)
        self.skipped = []             # (layer name, why)
        self.errors = []              # (layer name, message)
        self.cancelled = False
        self.project_path = None
        self._masks = {}
        self._base_mask = None

    # ---------------------------------------------------------------- the mask
    def _prepare_mask(self):
        """The clip polygon (the selected polygons, if any are selected), with valid geometry."""
        layer = self.clip_layer
        source = layer
        if layer.selectedFeatureCount() > 0:
            source = processing.run(
                "native:saveselectedfeatures", {"INPUT": layer, "OUTPUT": "memory:"},
                context=self.context, feedback=self.feedback)["OUTPUT"]
            self.log.append((Qgis.Info, f"Clipping to the {layer.selectedFeatureCount()} selected "
                                        f"polygon(s) of '{layer.name()}'"))
        self._base_mask = processing.run(
            "native:fixgeometries", {"INPUT": source, "OUTPUT": "memory:"},
            context=self.context, feedback=self.feedback)["OUTPUT"]

    def _mask_for(self, crs):
        """The mask in the CRS of the layer being clipped."""
        key = crs.toWkt()
        if key not in self._masks:
            if self._base_mask.crs() == crs:
                self._masks[key] = self._base_mask
            else:
                self._masks[key] = processing.run(
                    "native:reprojectlayer",
                    {"INPUT": self._base_mask, "TARGET_CRS": crs, "OUTPUT": "memory:"},
                    context=self.context, feedback=self.feedback)["OUTPUT"]
        return self._masks[key]

    # ---------------------------------------------------------------- one layer
    def _clip_vector(self, layer, path):
        processing.run(
            "native:clip",
            {"INPUT": layer, "OVERLAY": self._mask_for(layer.crs()), "OUTPUT": path},
            context=self.context, feedback=self.feedback)

    def _clip_raster(self, layer, path):
        nodata = None
        if not layer.dataProvider().sourceHasNoDataValue(1) and _is_float_raster(layer):
            nodata = -9999.0  # the cells outside the polygon need a value that means "no data"
        processing.run(
            "gdal:cliprasterbymasklayer",
            {
                "INPUT": layer, "MASK": self._mask_for(layer.crs()),
                "SOURCE_CRS": layer.crs(), "TARGET_CRS": layer.crs(),
                "NODATA": nodata, "ALPHA_BAND": False, "CROP_TO_CUTLINE": True,
                "KEEP_RESOLUTION": True, "SET_RESOLUTION": False,
                "X_RESOLUTION": None, "Y_RESOLUTION": None, "MULTITHREADING": True,
                "OPTIONS": "COMPRESS=LZW", "DATA_TYPE": 0, "EXTRA": "", "OUTPUT": path,
            },
            context=self.context, feedback=self.feedback)

    def _save_copy(self, layer, path):
        processing.run("native:savefeatures", {"INPUT": layer, "OUTPUT": path},
                       context=self.context, feedback=self.feedback)

    def _load_clipped(self, original, path, qml, name):
        """The clipped file as a layer called `name`, in the original's style."""
        if isinstance(original, QgsRasterLayer):
            new = QgsRasterLayer(path, name)
        else:
            new = QgsVectorLayer(path, name, "ogr")
        if not new.isValid():
            raise RuntimeError(f"the clipped file {path} could not be opened")
        new.setName(name)
        message, ok = new.loadNamedStyle(qml)
        if not ok:
            self.log.append((Qgis.Warning, f"Style of '{name}' not applied: {message}"))
        return new

    def _process_layer(self, layer, used):
        """The layer that goes into the new project for `layer` (None if it fails)."""
        name = layer.name()
        stem = _unique(safe_filename(name) + SUFFIX, used)
        qml = os.path.join(self.out_dir, stem + ".qml")

        if layer.id() == self.clip_layer.id():
            path = os.path.join(self.out_dir, stem + ".gpkg")
            self._save_copy(self._base_mask, path)  # the polygon(s) actually used
            layer.saveNamedStyle(qml)
            self.clipped.append((name, path))
            return self._load_clipped(layer, path, qml, name + SUFFIX)

        if isinstance(layer, QgsVectorLayer):
            if _geometry_type_is_null(layer):
                self.kept.append((name, "a table with no geometry cannot be clipped"))
                return layer.clone()
            path = os.path.join(self.out_dir, stem + ".gpkg")
            self._clip_vector(layer, path)
        elif isinstance(layer, QgsRasterLayer):
            source = layer.source().split("|")[0]
            if layer.providerType() != "gdal" or not os.path.exists(source):
                self.kept.append((name, "a web or non-file raster cannot be clipped"))
                return layer.clone()
            path = os.path.join(self.out_dir, stem + ".tif")
            self._clip_raster(layer, path)
        else:
            self.kept.append((name, f"{layer.__class__.__name__} layers are carried over unchanged"))
            return layer.clone()

        message, ok = layer.saveNamedStyle(qml)
        if not ok:
            self.log.append((Qgis.Warning, f"Style of '{name}' not saved: {message}"))
        self.clipped.append((name, path))
        return self._load_clipped(layer, path, qml, name + SUFFIX)

    # ---------------------------------------------------------------- the whole project
    def run(self):
        os.makedirs(self.out_dir, exist_ok=True)
        self._prepare_mask()
        root = self.project.layerTreeRoot()
        nodes = root.findLayers()
        replacement = {}  # original layer id -> layer for the new project
        used = set()
        for i, node in enumerate(nodes):
            layer = node.layer()
            name = node.name()
            if self.progress is not None and self.progress(i, len(nodes), name) is False:
                self.cancelled = True
                return self
            if layer is None or not layer.isValid():
                self.skipped.append((name, "the layer is not valid (broken source)"))
                continue
            if any(layer.name().endswith(s) for s in self.skip_names):
                self.skipped.append((name, "temporary preview layer"))
                continue
            try:
                new_layer = self._process_layer(layer, used)
            except Exception as e:  # one failing layer must not stop the rest
                self.errors.append((name, str(e).strip().splitlines()[0] if str(e).strip() else repr(e)))
                continue
            replacement[layer.id()] = (new_layer, node)
        if self.progress is not None:
            self.progress(len(nodes), len(nodes), "")
        self._write_project(root, replacement)
        return self

    def _write_project(self, root, replacement):
        new = QgsProject()
        src = self.project
        new.setCrs(src.crs())
        new.setTitle(src.title())
        try:
            new.setEllipsoid(src.ellipsoid())
            new.setTransformContext(src.transformContext())
        except Exception:
            pass

        def mirror(src_group, dst_group):
            for child in src_group.children():
                if QgsLayerTree.isGroup(child):
                    group = dst_group.addGroup(child.name())
                    group.setItemVisibilityChecked(child.itemVisibilityChecked())
                    group.setExpanded(child.isExpanded())
                    try:
                        group.setIsMutuallyExclusive(child.isMutuallyExclusive())
                    except Exception:
                        pass
                    mirror(child, group)
                elif QgsLayerTree.isLayer(child):
                    entry = replacement.get(child.layerId())
                    if entry is None:
                        continue
                    new_layer = entry[0]
                    new.addMapLayer(new_layer, False)
                    node = dst_group.addLayer(new_layer)
                    node.setItemVisibilityChecked(child.itemVisibilityChecked())
                    node.setExpanded(child.isExpanded())

        mirror(root, new.layerTreeRoot())

        # open on the clip polygon
        try:
            extent = self._base_mask.extent()
            to_project = QgsCoordinateTransform(self._base_mask.crs(), new.crs(), new.transformContext())
            new.viewSettings().setDefaultViewExtent(
                QgsReferencedRectangle(to_project.transformBoundingBox(extent), new.crs()))
        except Exception:
            pass

        stem = os.path.splitext(os.path.basename(src.fileName()))[0] if src.fileName() else ""
        stem = safe_filename(stem) if stem else "clipped_layers"
        path = os.path.join(self.out_dir, stem + SUFFIX + ".qgz")
        if not new.write(path):
            raise RuntimeError(f"the new project could not be written to {path}")
        self.project_path = path

    # ---------------------------------------------------------------- report
    def summary(self):
        text = (f"{len(self.clipped)} layer(s) clipped, {len(self.kept)} carried over unchanged, "
                f"{len(self.skipped)} skipped, {len(self.errors)} failed")
        if self.project_path:
            text += f"; project: {self.project_path}"
        return text
