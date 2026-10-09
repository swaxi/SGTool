"""
The processing history of a grid, read from its SGTool provenance.

Every grid SGTool saves records the operation and parameters that made it and
embeds the provenance of the grid it was made from, so the whole chain of steps
is there. history_steps() unrolls that chain, oldest step first, so it can be
turned into a recipe and replayed on another grid (see sgt_workflow.py).

Pure Python (no QGIS), so it can be tested on its own.
"""


def _step_from(node):
    params = {}
    for p in node.findall("parameters/parameter"):
        params[p.get("name", "")] = (p.text or "").strip()
    sources = node.findall("sources/source")
    software = node.find("software")
    return {
        "operation": (node.findtext("operation") or "").strip(),
        "parameters": params,
        "created": node.findtext("created") or "",
        "version": software.get("version", "") if software is not None else "",
        "output": node.findtext("output/path") or "",
        "source": (sources[0].findtext("path") or "") if sources else "",
        "n_sources": len(sources),
    }


def _earlier_step(node):
    """The provenance element of the grid this step was made from, or None."""
    sources = node.findall("sources/source")
    if not sources:
        return None
    for sx in sources[0].findall("metadata/sourceXml"):
        if sx.get("type") == "sgtool":
            for child in sx:
                if child.tag == "SGToolMetadata":
                    return child
    return None


def history_steps(root):
    """(original source path, steps) for a provenance element.

    steps are dicts (operation, parameters, created, version, output, source,
    n_sources) in the order they were applied, oldest first. The original
    source path is what the first step was made from (a grid or a points
    layer): the start of the chain.
    """
    steps = []
    node = root
    guard = 0
    while node is not None and guard < 200:
        steps.append(_step_from(node))
        node = _earlier_step(node)
        guard += 1
    steps.reverse()
    origin = steps[0]["source"] if steps else ""
    return origin, steps
