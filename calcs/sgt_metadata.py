"""
SGTool provenance metadata

Writes a small XML sidecar next to every file the plugin saves:

    <output file>.sgt.xml        e.g.  grid_DirC.tif  ->  grid_DirC.tif.sgt.xml

The sidecar records
  * when the file was created and by which SGTool version
  * the source file(s) it was made from, plus any XML metadata those sources
    carry (their own .sgt.xml - so the lineage chain is preserved - and
    GDAL .aux.xml / other .xml sidecars)
  * the operation and the parameters that were set to create it

No QGIS or GDAL calls here, so it can be used from any front end.
Failing to write metadata must never break processing, so write_sgt_metadata()
swallows errors and returns None instead.
"""

import os
import datetime
import xml.etree.ElementTree as ET

SCHEMA_VERSION = "1"
SIDECAR_SUFFIX = ".sgt.xml"
MAX_EMBED_BYTES = 2 * 1024 * 1024  # don't embed huge sidecars


def sidecar_path(output_path):
    return str(output_path) + SIDECAR_SUFFIX


def _clean_path(path):
    """Strip QGIS provider suffixes such as 'file.tif|layername=x'."""
    return os.path.abspath(str(path).split("|")[0])


def _value_text(value):
    if isinstance(value, (list, tuple)):
        return ", ".join(_value_text(v) for v in value)
    return str(value)


def _embed_xml_file(parent, kind, file_path):
    """Append the contents of an existing XML file under parent, if readable."""
    try:
        if not os.path.isfile(file_path) or os.path.getsize(file_path) > MAX_EMBED_BYTES:
            return
        node = ET.SubElement(parent, "sourceXml", {"type": kind, "file": file_path})
        try:
            node.append(ET.parse(file_path).getroot())
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
    base = os.path.splitext(path)[0]
    candidates = [
        ("sgtool", path + SIDECAR_SUFFIX),   # lineage of the source
        ("gdal-pam", path + ".aux.xml"),     # GDAL statistics/CRS/etc.
        ("sidecar", path + ".xml"),          # e.g. ArcGIS / ISO metadata
        ("sidecar", base + ".xml"),
    ]
    seen = set()
    for kind, candidate in candidates:
        if candidate not in seen:
            seen.add(candidate)
            _embed_xml_file(meta, kind, candidate)
    if len(meta) == 0:
        src.remove(meta)  # nothing to record


def write_sgt_metadata(
    output_path,
    source_path=None,
    operation=None,
    parameters=None,
    software_version=None,
):
    """Write <output_path>.sgt.xml and return its path (None on failure).

    Parameters
    ----------
    output_path : str
        The file that was saved.
    source_path : str or list of str, optional
        File(s) the output was created from.
    operation : str, optional
        Short name of what was done, e.g. "DirClean" or "Euler deconvolution".
    parameters : dict, optional
        Parameters that were set to create the output (flat name -> value).
    software_version : str, optional
        SGTool version string.
    """
    try:
        root = ET.Element("SGToolMetadata", {"schema": SCHEMA_VERSION})
        ET.SubElement(root, "created").text = (
            datetime.datetime.now().astimezone().isoformat(timespec="seconds")
        )
        software = ET.SubElement(root, "software", {"name": "SGTool"})
        if software_version:
            software.set("version", str(software_version))

        out = ET.SubElement(root, "output")
        ET.SubElement(out, "path").text = os.path.abspath(str(output_path))
        if operation:
            ET.SubElement(root, "operation").text = str(operation)

        params = ET.SubElement(root, "parameters")
        for name, value in (parameters or {}).items():
            p = ET.SubElement(params, "parameter", {"name": str(name)})
            p.text = _value_text(value)

        sources = ET.SubElement(root, "sources")
        if source_path:
            if isinstance(source_path, (list, tuple)):
                paths = [s for s in source_path if s]
            else:
                paths = [source_path]
            for i, s in enumerate(paths):
                _add_source(sources, i, s)

        if hasattr(ET, "indent"):  # Python 3.9+
            ET.indent(root, space="  ")
        path = sidecar_path(output_path)
        ET.ElementTree(root).write(path, encoding="utf-8", xml_declaration=True)
        return path
    except Exception as e:  # never let provenance break processing
        print(f"SGTool: could not write metadata for {output_path}: {e}")
        return None
