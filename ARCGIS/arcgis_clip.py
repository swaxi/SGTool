"""
Clip every layer of an ArcGIS Pro map to a polygon, and save a copy of the project that uses the
clipped layers (the ArcGIS Pro counterpart of the QGIS plugin's Clip All Layers).

For each layer of the map the clipped data is written to the output folder:

* vector layers go into one file geodatabase, <project>_clip.gdb, as <name>_clip feature classes
  (the clip polygon itself too, as <polygon name>_clip);
* grids are written as <name>_clip.tif beside it;
* each original layer's symbology is saved as <name>_clip.lyrx.

<project>_clip.aprx is then saved there: a copy of the project (so groups, order, symbology,
bookmarks and layouts are all kept) in which each clipped layer points at its new data and is
renamed <name>_clip. Layers that cannot be clipped (web layers, annotation, tables) are left as
they are.

Only the chosen map is changed; other maps in the copy still point at the original data.
"""

import os
import copy

import arcpy

SUFFIX = "_clip"


RASTER_EXTENSIONS = (".tif", ".tiff", ".img", ".jp2", ".grd", ".ers", ".asc")


def base_name(name):
    """A layer name without the file extension grids are often named with (Magnetics.tif)."""
    return os.path.splitext(name)[0] if name.lower().endswith(RASTER_EXTENSIONS) else name


def safe_filename(name):
    keep = (".", "_", "-")
    text = "".join(c if (c.isalnum() or c in keep) else "_" for c in name).strip()
    return text or "layer"


class ClipJob:
    """One run of clipping a map's layers to a polygon."""

    def __init__(self, aprx, polygon, out_folder, map_name=None, check_cancel=None, messages=None,
                 skip_suffixes=("_preview",)):
        """aprx: an arcpy.mp.ArcGISProject. polygon: the polygon layer (a layer name in the map, or
        the path of a feature class; if some of its features are selected only those are used).
        check_cancel(): called between steps, may raise to stop (see calcs.sgt_cancel).
        messages: object with addMessage / addWarningMessage, or None."""
        self.aprx = aprx
        self.polygon = polygon
        self.out_folder = out_folder
        self.map_name = map_name
        self.check_cancel = check_cancel or (lambda *a: None)
        self.messages = messages
        self.skip_suffixes = tuple(skip_suffixes)
        self.clipped = []     # (layer name, path of the data written)
        self.kept = []        # (layer name, why it is carried over unchanged)
        self.skipped = []     # (layer name, why)
        self.errors = []      # (layer name, message)
        self.aprx_path = None
        self.gdb = None
        self._used = set()
        self._repoint = {}    # long name of the original layer -> how the new project reaches its data

    def _say(self, text, warning=False):
        if self.messages is not None:
            (self.messages.addWarningMessage if warning else self.messages.addMessage)(text)

    # ---------------------------------------------------------------- names
    def _stem(self, name):
        stem = safe_filename(base_name(name)) + SUFFIX
        candidate, n = stem, 2
        while candidate.lower() in self._used:
            candidate = f"{stem}_{n}"
            n += 1
        self._used.add(candidate.lower())
        return candidate

    def _fc_name(self, stem):
        return arcpy.ValidateTableName(stem, self.gdb)

    # ---------------------------------------------------------------- the map
    def _the_map(self):
        if self.map_name:
            found = self.aprx.listMaps(self.map_name)
            if not found:
                raise ValueError(f"the project has no map called '{self.map_name}'")
            return found[0]
        return self.aprx.activeMap or self.aprx.listMaps()[0]

    # ---------------------------------------------------------------- the mask
    def _prepare_mask(self, stem_polygon):
        """The clip polygon(s) saved in the geodatabase, with valid geometry."""
        name = self._fc_name(stem_polygon)
        self.mask = os.path.join(self.gdb, name)
        arcpy.management.CopyFeatures(self.polygon, self.mask)  # honours a selection
        arcpy.management.RepairGeometry(self.mask, "DELETE_NULL")
        self._used.add(stem_polygon.lower())
        self.mask_name = name
        self._mask_sr = arcpy.Describe(self.mask).spatialReference

    def _mask_in(self, sr):
        """The mask in the given coordinate system (the saved one if it already matches)."""
        if sr is None or sr.name == self._mask_sr.name:
            return self.mask
        key = sr.name
        cache = self.__dict__.setdefault("_mask_cache", {})
        if key not in cache:
            out = os.path.join("memory", "mask_%d" % len(cache))
            arcpy.management.Project(self.mask, out, sr)
            cache[key] = out
        return cache[key]

    # ---------------------------------------------------------------- layers
    def _is_mask_layer(self, layer):
        polygon = str(self.polygon)
        try:
            same_source = os.path.normcase(os.path.normpath(layer.dataSource)) == os.path.normcase(os.path.normpath(polygon))
        except Exception:
            same_source = False
        return layer.isFeatureLayer and (layer.name == polygon or same_source)

    @staticmethod
    def _is_vector(layer):
        if not layer.isFeatureLayer or layer.isWebLayer:
            return False
        try:
            return arcpy.Describe(layer).shapeType in ("Point", "Multipoint", "Polyline", "Polygon")
        except Exception:
            return False

    def _clip_vector(self, layer, stem):
        fc_name = self._fc_name(stem)
        out_fc = os.path.join(self.gdb, fc_name)
        source = layer.dataSource
        query = layer.definitionQuery if hasattr(layer, "definitionQuery") else ""
        src = layer
        if source and arcpy.Exists(source):
            # the data with the layer's filter, but not its selection
            tmp = "sgt_clip_src_%d" % len(self._used)
            src = arcpy.management.MakeFeatureLayer(source, tmp, query or "").getOutput(0)
        arcpy.analysis.Clip(src, self._mask_in(arcpy.Describe(layer).spatialReference), out_fc)
        return {"kind": "fc", "dataset": fc_name, "path": out_fc}

    def _clip_raster(self, layer, stem):
        source = layer.dataSource
        out_tif = os.path.join(self.out_folder, stem + ".tif")
        sr = arcpy.Describe(source).spatialReference
        mask = self._mask_in(sr)
        e = arcpy.Describe(mask).extent
        rectangle = f"{e.XMin} {e.YMin} {e.XMax} {e.YMax}"
        arcpy.management.Clip(source, rectangle, out_tif, mask, "#", "ClippingGeometry",
                              "NO_MAINTAIN_EXTENT")
        return {"kind": "raster", "dataset": stem + ".tif", "path": out_tif}

    def _process(self, layer):
        name = layer.name
        if layer.isGroupLayer:
            return
        if layer.isBroken:
            self.skipped.append((name, "the layer's data source is broken"))
            return
        if any(name.endswith(s) for s in self.skip_suffixes):
            self.skipped.append((name, "temporary preview layer"))
            return
        long_name = layer.longName
        if self._is_mask_layer(layer):
            # the polygon layer: the saved copy of the polygon(s) actually used replaces it
            stem = self.mask_name
            result = {"kind": "fc", "dataset": self.mask_name, "path": self.mask}
        elif self._is_vector(layer):
            stem = self._stem(name)
            result = self._clip_vector(layer, stem)
        elif layer.isRasterLayer and not layer.isWebLayer and layer.dataSource \
                and arcpy.Exists(layer.dataSource) and not str(layer.dataSource).lower().startswith("http"):
            stem = self._stem(name)
            result = self._clip_raster(layer, stem)
        else:
            self.kept.append((name, "a web layer, annotation or table cannot be clipped"))
            return
        try:
            layer.saveACopy(os.path.join(self.out_folder, stem + ".lyrx"))
        except Exception as e:
            self._say(f"The symbology of '{name}' could not be saved as a .lyrx file: {e}", True)
        self._repoint[long_name] = result
        self.clipped.append((name, result["path"]))

    # ---------------------------------------------------------------- the whole project
    def run(self):
        arcpy.env.overwriteOutput = True
        os.makedirs(self.out_folder, exist_ok=True)
        stem_project = os.path.splitext(os.path.basename(self.aprx.filePath))[0] if self.aprx.filePath else ""
        stem_project = safe_filename(stem_project) if stem_project else "clipped_layers"
        gdb_name = stem_project + SUFFIX + ".gdb"
        self.gdb = os.path.join(self.out_folder, gdb_name)
        if not arcpy.Exists(self.gdb):
            arcpy.management.CreateFileGDB(self.out_folder, gdb_name)
        the_map = self._the_map()
        polygon_name = os.path.basename(str(self.polygon)).replace(".shp", "")
        self._prepare_mask(safe_filename(polygon_name) + SUFFIX)

        layers = the_map.listLayers()
        for i, layer in enumerate(layers):
            self.check_cancel(i, len(layers), layer.name)
            self._say(f"Clipping {layer.name} ({i + 1} of {len(layers)})")
            try:
                self._process(layer)
            except Exception as e:  # one failing layer must not stop the rest
                text = str(e).strip().splitlines()
                self.errors.append((layer.name, text[0] if text else repr(e)))
                self._repoint.pop(layer.longName, None)
        self.check_cancel(len(layers), len(layers), "")
        self._write_project(the_map)
        return self

    def _write_project(self, the_map):
        stem = os.path.splitext(os.path.basename(self.aprx.filePath))[0] if self.aprx.filePath else ""
        stem = safe_filename(stem) if stem else "clipped_layers"
        path = os.path.join(self.out_folder, stem + SUFFIX + ".aprx")
        self.aprx.saveACopy(path)
        new = arcpy.mp.ArcGISProject(path)
        new_map = new.listMaps(the_map.name)[0]
        failed_names = {n for n, _ in self.errors}
        failed_names |= {n for n, _ in self.skipped}
        # work out what to do to each layer first: renaming and removing change the listing
        todo = [(layer, self._repoint.get(layer.longName)) for layer in new_map.listLayers()]
        for layer, target in todo:
            if layer.isGroupLayer:
                continue
            if target is None:
                if layer.name in failed_names:  # not clipped: it would still point at the original data
                    new_map.removeLayer(layer)
                continue
            old_name = layer.name
            try:
                if target["kind"] == "fc":
                    props = {"dataset": target["dataset"], "workspace_factory": "File Geodatabase",
                             "connection_info": {"database": self.gdb}}
                else:
                    props = {"dataset": target["dataset"], "workspace_factory": "Raster",
                             "connection_info": {"database": self.out_folder}}
                layer.updateConnectionProperties(copy.deepcopy(layer.connectionProperties), props,
                                                 validate=False)
                layer.name = base_name(old_name) + SUFFIX
            except Exception as e:
                self.errors.append((old_name, "could not point the new project at the clipped data: "
                                              + (str(e).strip().splitlines() or [repr(e)])[0]))
                new_map.removeLayer(layer)
        new.save()
        self.aprx_path = path

    def summary(self):
        text = (f"{len(self.clipped)} layer(s) clipped, {len(self.kept)} carried over unchanged, "
                f"{len(self.skipped)} skipped, {len(self.errors)} failed")
        if self.aprx_path:
            text += f"; project: {self.aprx_path}"
        return text
