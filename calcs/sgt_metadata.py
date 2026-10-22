"""
SGTool provenance metadata

Records, for every file the plugin saves,
  * when it was created and by which SGTool version
  * the source file(s) it was made from, plus any XML metadata those sources
    carry - including their own SGTool provenance, so the whole chain of
    processing steps is preserved
  * the operation and the parameters that were set to create it

Where it is stored
  GeoTIFF outputs   embedded in the file itself, in its own GDAL metadata
                    domain ("SGTOOL", item "provenance_xml"), the same
                    mechanism the Noddy grid import uses. It travels with the
                    file and survives copying/renaming.
  other outputs     shapefiles, csv/txt etc. cannot carry metadata, so they get
                    a small XML sidecar:  <file>.sgt.xml
                    (also used as a fallback if a GeoTIFF cannot be updated,
                    and read for files made by earlier SGTool versions)

Failing to write metadata must never break processing, so write_sgt_metadata()
swallows errors and returns None instead. No QGIS calls here; GDAL is only
imported (lazily) to embed in / read GeoTIFFs.
"""

import os
import html
import datetime
import xml.etree.ElementTree as ET  # nosec B405 - building trees; parsing goes via _safe_fromstring

try:  # hardened parser if available (not bundled with QGIS)
    from defusedxml import ElementTree as _SafeET
except ImportError:
    _SafeET = None



def _safe_fromstring(data):
    """Parse XML text/bytes, refusing DTDs and entity declarations (entity-expansion attacks)."""
    if _SafeET is not None:
        return _SafeET.fromstring(data)
    raw = data.encode("utf-8") if isinstance(data, str) else data
    low = raw.lower()
    if b"<!doctype" in low or b"<!entity" in low:
        raise ET.ParseError("DTD / entity declarations are not allowed")
    return ET.fromstring(raw)  # nosec B314 - DOCTYPE/ENTITY rejected above


def parse_xml_file(path):
    """Root element of an XML file, parsed safely."""
    with open(path, "rb") as f:
        return _safe_fromstring(f.read())


SCHEMA_VERSION = "1"
SIDECAR_SUFFIX = ".sgt.xml"
GDAL_DOMAIN = "SGTOOL"
GDAL_ITEM = "provenance_xml"
MAX_EMBED_BYTES = 2 * 1024 * 1024  # don't embed huge sidecars


# ---------------------------------------------------------------- sidecars
def sidecar_path(output_path):
    return str(output_path) + SIDECAR_SUFFIX


def remove_sgt_metadata(output_path):
    """Delete <output_path>.sgt.xml if present (call when the output itself is
    deleted or replaced, so a stale sidecar is not left describing it).
    Metadata embedded in a GeoTIFF goes with the file, so needs no clean-up.
    Returns True if a sidecar was removed."""
    path = sidecar_path(output_path)
    try:
        if os.path.exists(path):
            os.remove(path)
            return True
    except OSError:
        pass
    return False


def _modified_text(timestamp):
    return datetime.datetime.fromtimestamp(timestamp).astimezone().isoformat(
        timespec="seconds"
    )


def sidecar_status(output_path, tolerance_seconds=2.0):
    """Does a sidecar still describe the file beside it?

    Sidecars record the output's size and modified time when written.
    Returns
      "current"  size and modified time still match
      "stale"    the file has changed since (e.g. overwritten by another tool)
      "missing"  no sidecar, or the output file itself is gone
      "unknown"  sidecar has no size/time record
    (Not needed for GeoTIFFs: their metadata is inside the file.)
    """
    path = sidecar_path(output_path)
    if not os.path.exists(path) or not os.path.exists(output_path):
        return "missing"
    try:
        out = parse_xml_file(path).find("output")
        size = out.findtext("sizeBytes") if out is not None else None
        modified = out.findtext("modified") if out is not None else None
        if size is None or modified is None:
            return "unknown"
        st = os.stat(output_path)
        recorded = datetime.datetime.fromisoformat(modified).timestamp()
        if int(size) == st.st_size and abs(recorded - st.st_mtime) <= tolerance_seconds:
            return "current"
        return "stale"
    except (OSError, ET.ParseError, ValueError):
        return "unknown"


# ---------------------------------------------------------------- helpers
def _clean_path(path):
    """Strip QGIS provider suffixes such as 'file.tif|layername=x'."""
    return os.path.abspath(str(path).split("|")[0])


def _value_text(value):
    if isinstance(value, (list, tuple)):
        return ", ".join(_value_text(v) for v in value)
    return str(value)


def _to_text(root):
    if hasattr(ET, "indent"):  # Python 3.9+
        ET.indent(root, space="  ")
    return ET.tostring(root, encoding="unicode")


def _is_geotiff(path):
    if not str(path).lower().endswith((".tif", ".tiff")):
        return False
    try:
        from osgeo import gdal

        ds = gdal.Open(str(path))
        ok = ds is not None and ds.GetDriver().ShortName == "GTiff"
        ds = None
        return ok
    except Exception:
        return False


def _read_embedded(path):
    """SGTool provenance embedded in a GeoTIFF, or None."""
    try:
        from osgeo import gdal

        ds = gdal.Open(str(path))
        text = ds.GetMetadataItem(GDAL_ITEM, GDAL_DOMAIN) if ds is not None else None
        ds = None
        return _safe_fromstring(text) if text else None
    except Exception:
        return None


def _read_sidecar(path):
    try:
        sc = sidecar_path(path)
        if os.path.isfile(sc) and os.path.getsize(sc) <= MAX_EMBED_BYTES:
            return parse_xml_file(sc)
    except (OSError, ET.ParseError):
        pass
    return None


def _read_with_origin(path):
    """(provenance root, "embedded" | "sidecar") or (None, None)."""
    if not os.path.isfile(path):
        return None, None
    if str(path).lower().endswith((".tif", ".tiff")):
        root = _read_embedded(path)
        if root is not None:
            return root, "embedded"
    root = _read_sidecar(path)
    if root is not None:
        return root, "sidecar"
    return None, None


def read_sgt_metadata(path):
    """The SGTool provenance of a file as an XML Element, or None.

    Looks inside GeoTIFFs first, then for a .sgt.xml sidecar (non-GeoTIFF
    outputs, and files made by earlier SGTool versions)."""
    return _read_with_origin(_clean_path(path))[0]


def _drop_histograms(root):
    """Remove GDAL <Histograms> blocks (bulky, and repeated at every level of
    a processing chain) but keep the band statistics beside them."""
    for parent in list(root.iter()):
        for child in list(parent):
            if child.tag == "Histograms":
                parent.remove(child)


def _embed_xml_file(parent, kind, file_path):
    """Append the contents of an existing XML file under parent, if readable."""
    try:
        if not os.path.isfile(file_path) or os.path.getsize(file_path) > MAX_EMBED_BYTES:
            return
        node = ET.SubElement(parent, "sourceXml", {"type": kind, "file": file_path})
        try:
            embedded = parse_xml_file(file_path)
            if kind == "gdal-pam":
                _drop_histograms(embedded)
            node.append(embedded)
        except ET.ParseError:
            # not well-formed: keep the raw text so nothing is lost
            with open(file_path, "r", encoding="utf-8", errors="replace") as f:
                node.text = f.read()
    except OSError:
        pass


def _add_source(parent, index, source_path):
    src = ET.SubElement(parent, "source", {"index": str(index)})
    path = _clean_path(source_path)
    ET.SubElement(src, "path").text = path
    ET.SubElement(src, "exists").text = str(os.path.exists(path)).lower()

    meta = ET.SubElement(src, "metadata")
    # the source's own SGTool provenance: this is what chains the steps
    provenance, origin = _read_with_origin(path)
    if provenance is not None:
        node = ET.SubElement(
            meta, "sourceXml", {"type": "sgtool", "origin": origin, "file": path}
        )
        node.append(provenance)
    base = os.path.splitext(path)[0]
    seen = set()
    for kind, candidate in (
        ("gdal-pam", path + ".aux.xml"),  # GDAL statistics/CRS/etc.
        ("sidecar", path + ".xml"),  # e.g. ArcGIS / ISO metadata, .grd.xml
        ("sidecar", base + ".xml"),
    ):
        if candidate not in seen:
            seen.add(candidate)
            _embed_xml_file(meta, kind, candidate)
    if len(meta) == 0:
        src.remove(meta)  # nothing to record


def _build(output_path, source_path, operation, parameters, software_version, file_stats):
    root = ET.Element("SGToolMetadata", {"schema": SCHEMA_VERSION})
    ET.SubElement(root, "created").text = (
        datetime.datetime.now().astimezone().isoformat(timespec="seconds")
    )
    software = ET.SubElement(root, "software", {"name": "SGTool"})
    if software_version:
        software.set("version", str(software_version))

    out = ET.SubElement(root, "output")
    ET.SubElement(out, "path").text = os.path.abspath(str(output_path))
    if file_stats:
        try:
            # sidecar only: lets sidecar_status() spot a file changed later
            st = os.stat(output_path)
            ET.SubElement(out, "sizeBytes").text = str(st.st_size)
            ET.SubElement(out, "modified").text = _modified_text(st.st_mtime)
        except OSError:
            pass  # output not on disk (yet): nothing to record
    if operation:
        ET.SubElement(root, "operation").text = str(operation)

    params = ET.SubElement(root, "parameters")
    for name, value in (parameters or {}).items():
        p = ET.SubElement(params, "parameter", {"name": str(name)})
        p.text = _value_text(value)

    sources = ET.SubElement(root, "sources")
    if source_path:
        paths = (
            [s for s in source_path if s]
            if isinstance(source_path, (list, tuple))
            else [source_path]
        )
        for i, s in enumerate(paths):
            _add_source(sources, i, s)
    return root


def _embed_in_geotiff(path, root):
    """Store the provenance in the GeoTIFF's own metadata. True on success."""
    try:
        from osgeo import gdal

        text = _to_text(root)
        ds = gdal.Open(str(path), gdal.GA_Update)
        if ds is None:
            return False
        ds.SetMetadataItem(GDAL_ITEM, text, GDAL_DOMAIN)
        # a few plain items so the essentials show in QGIS layer properties
        ds.SetMetadataItem("SGTOOL_OPERATION", root.findtext("operation") or "")
        ds.SetMetadataItem("SGTOOL_CREATED", root.findtext("created") or "")
        version = root.find("software").get("version", "")
        ds.SetMetadataItem("SGTOOL_VERSION", version)
        ds.FlushCache()
        ds = None
        # make sure it really landed in the file (not an aux.xml side file)
        return _read_embedded(path) is not None
    except Exception:
        return False


# ---------------------------------------------------------------- public
def write_sgt_metadata(
    output_path,
    source_path=None,
    operation=None,
    parameters=None,
    software_version=None,
):
    """Record provenance for output_path; returns where it was stored (the
    output path when embedded in a GeoTIFF, else the sidecar path), or None.

    Parameters
    ----------
    output_path : str
        The file that was saved (must already exist on disk).
    source_path : str or list of str, optional
        File(s) the output was created from.
    operation : str, optional
        Short name of what was done, e.g. "Derivative" or "Euler deconvolution".
    parameters : dict, optional
        Parameters that were set to create the output (flat name -> value).
    software_version : str, optional
        SGTool version string.
    """
    try:
        if _is_geotiff(output_path):
            root = _build(
                output_path, source_path, operation, parameters,
                software_version, file_stats=False,
            )
            if _embed_in_geotiff(output_path, root):
                remove_sgt_metadata(output_path)  # drop any old-style sidecar
                return os.path.abspath(str(output_path))
            # could not update the GeoTIFF (e.g. locked): fall through to sidecar
        root = _build(
            output_path, source_path, operation, parameters,
            software_version, file_stats=True,
        )
        path = sidecar_path(output_path)
        if hasattr(ET, "indent"):
            ET.indent(root, space="  ")
        ET.ElementTree(root).write(path, encoding="utf-8", xml_declaration=True)
        return path
    except Exception as e:  # never let provenance break processing
        print(f"SGTool: could not write metadata for {output_path}: {e}")
        return None


# ---------------------------------------------------------------- display
def _stat_items(source_el):
    """GDAL band statistics recorded for a source (name -> value)."""
    stats = {}
    meta = source_el.find("metadata")
    if meta is None:
        return stats
    for sx in meta.findall("sourceXml"):
        if sx.get("type") == "gdal-pam":
            for mdi in sx.iter("MDI"):
                key = mdi.get("key", "")
                if key.startswith("STATISTICS_"):
                    stats[key[len("STATISTICS_"):].lower()] = (mdi.text or "").strip()
    return stats


def _format_step(root, parts, depth):
    esc = html.escape
    pad = 16 * depth
    parts.append(f'<div style="margin-left:{pad}px">')
    operation = root.findtext("operation") or "Unnamed step"
    parts.append(f"<h3 style='margin-bottom:2px'>{esc(operation)}</h3>")

    created = root.findtext("created") or "unknown time"
    sw = root.find("software")
    version = f" {sw.get('version')}" if sw is not None and sw.get("version") else ""
    parts.append(f"<p style='margin:0'><b>Created:</b> {esc(created)} "
                 f"&nbsp; <b>SGTool</b>{esc(version)}</p>")
    out = root.findtext("output/path")
    if out:
        parts.append(f"<p style='margin:0'><b>File:</b> {esc(out)}</p>")

    params = root.findall("parameters/parameter")
    if params:
        parts.append("<p style='margin:6px 0 0 0'><b>Parameters</b></p><table cellpadding='2'>")
        for p in params:
            parts.append(f"<tr><td>{esc(p.get('name', ''))}</td>"
                         f"<td>{esc((p.text or '').strip())}</td></tr>")
        parts.append("</table>")
    else:
        parts.append("<p style='margin:6px 0 0 0'><i>No parameters</i></p>")

    for src in root.findall("sources/source"):
        parts.append(f"<p style='margin:8px 0 0 0'><b>Made from:</b> "
                     f"{esc(src.findtext('path') or '')}")
        if src.findtext("exists") == "false":
            parts.append(" <i>(file no longer found)</i>")
        parts.append("</p>")
        stats = _stat_items(src)
        if stats:
            parts.append("<p style='margin:0'><i>Source statistics: "
                         + ", ".join(f"{esc(k)} {esc(v)}" for k, v in stats.items())
                         + "</i></p>")
        meta = src.find("metadata")
        nested = []
        others = []
        if meta is not None:
            for sx in meta.findall("sourceXml"):
                if sx.get("type") == "sgtool":
                    nested.extend(list(sx))
                elif sx.get("type") != "gdal-pam":
                    others.append(os.path.basename(sx.get("file", "")))
        if others:
            parts.append("<p style='margin:0'><i>Also carries XML metadata: "
                         + esc(", ".join(others)) + "</i></p>")
        if nested:
            for step in nested:
                _format_step(step, parts, depth + 1)
        else:
            parts.append("<p style='margin:0'><i>No earlier SGTool history "
                         "(original data or made outside SGTool)</i></p>")
    parts.append("</div>")


def format_sgt_metadata_html(root):
    """Readable HTML for a provenance element, newest step first, with each
    earlier step nested beneath the file it produced."""
    parts = []
    _format_step(root, parts, 0)
    return "".join(parts)


def _text_step(root, lines, depth):
    pad = "    " * depth
    lines.append(f"{pad}{root.findtext('operation') or 'Unnamed step'}")
    sw = root.find("software")
    version = f" (SGTool {sw.get('version')})" if sw is not None and sw.get("version") else ""
    lines.append(f"{pad}  Created: {root.findtext('created') or 'unknown time'}{version}")
    out = root.findtext("output/path")
    if out:
        lines.append(f"{pad}  File: {out}")
    params = root.findall("parameters/parameter")
    if params:
        lines.append(f"{pad}  Parameters:")
        for p in params:
            lines.append(f"{pad}    {p.get('name', '')} = {(p.text or '').strip()}")
    else:
        lines.append(f"{pad}  No parameters")
    for src in root.findall("sources/source"):
        note = "  (file no longer found)" if src.findtext("exists") == "false" else ""
        lines.append(f"{pad}  Made from: {src.findtext('path') or ''}{note}")
        stats = _stat_items(src)
        if stats:
            lines.append(
                f"{pad}  Source statistics: "
                + ", ".join(f"{k} {v}" for k, v in stats.items())
            )
        nested, others = [], []
        meta = src.find("metadata")
        if meta is not None:
            for sx in meta.findall("sourceXml"):
                if sx.get("type") == "sgtool":
                    nested.extend(list(sx))
                elif sx.get("type") != "gdal-pam":
                    others.append(os.path.basename(sx.get("file", "")))
        if others:
            lines.append(f"{pad}  Also carries XML metadata: {', '.join(others)}")
        if nested:
            for step in nested:
                lines.append("")
                _text_step(step, lines, depth + 1)
        else:
            lines.append(
                f"{pad}  No earlier SGTool history (original data or made outside SGTool)"
            )


def format_sgt_metadata_text(root):
    """The same history as format_sgt_metadata_html, as plain text (newest step
    first, each earlier step indented beneath the file it produced), for windows
    that cannot show HTML."""
    lines = []
    _text_step(root, lines, 0)
    return "\n".join(lines)


def save_sgt_metadata_xml(root, grid_path):
    """Write a provenance element to <grid_path>.sgt.xml (replacing any existing
    file) and return that path. Raises OSError if it cannot be written."""
    path = sidecar_path(_clean_path(grid_path))
    if hasattr(ET, "indent"):
        ET.indent(root, space="  ")
    ET.ElementTree(root).write(path, encoding="utf-8", xml_declaration=True)
    return path


def sgt_metadata_xml_text(root):
    """Pretty-printed XML text of a provenance element."""
    return _to_text(root)
