"""
One point per sample from a "long" petrophysics point layer.

The input has one row per measurement (columns SAMPLENO, PROPERTY, Value);
the output is a temporary layer with one point per SAMPLENO and a
<ShortName>_val field for each measurement type. Other columns that have the
same value on every row of a sample (e.g. lithology) are copied across.

Run in the QGIS Python console with the point layer selected in the Layers
panel. If some of its features are selected, only those are used.
"""
import re
from collections import OrderedDict

from qgis.core import QgsFeature, QgsField, QgsGeometry, QgsProject, QgsVectorLayer, QgsWkbTypes

SAMPLE_FIELD = "SAMPLENO"
MEASUREMENT_FIELD = "PROPERTY"
VALUE_FIELD = "Value"

# short names for common measurement types (matched ignoring case and units in brackets)
SHORT_NAMES = {
    "magnetic susceptibility": "MagSus",
    "susceptibility": "MagSus",
    "mag sus": "MagSus",
    "density": "Dens",
    "bulk density": "BulkDens",
    "dry bulk density": "DryDens",
    "wet bulk density": "WetDens",
    "saturated bulk density": "SatDens",
    "grain density": "GrainDens",
    "specific gravity": "SG",
    "porosity": "Poros",
    "apparent porosity": "AppPoros",
    "remanence": "Rem",
    "remanent magnetisation": "Rem",
    "remanent magnetization": "Rem",
    "natural remanent magnetisation": "NRM",
    "natural remanent magnetization": "NRM",
    "koenigsberger ratio": "Q",
    "q ratio": "Q",
    "electrical conductivity": "Cond",
    "conductivity": "Cond",
    "resistivity": "Res",
    "chargeability": "Charg",
    "p-wave velocity": "Vp",
    "s-wave velocity": "Vs",
    "thermal conductivity": "ThCond",
    "potassium": "K",
    "thorium": "Th",
    "uranium": "U",
}

try:  # QGIS 3.38+ and QGIS 4 (Qt6)
    from qgis.PyQt.QtCore import QMetaType
    DOUBLE, STRING = QMetaType.Type.Double, QMetaType.Type.QString
    QgsField("test", DOUBLE)
except (ImportError, AttributeError, TypeError):  # older QGIS
    from qgis.PyQt.QtCore import QVariant
    DOUBLE, STRING = QVariant.Double, QVariant.String


def short_name(measurement):
    key = re.sub(r"\(.*?\)|\[.*?\]", "", measurement).strip().lower()
    key = re.sub(r"\s+", " ", key)
    if key in SHORT_NAMES:
        return SHORT_NAMES[key]
    words = re.findall(r"[A-Za-z0-9]+", key)
    if not words:
        return "Meas"
    return "".join(w[:4].capitalize() for w in words)[:16]


def is_null(value):
    return value is None or str(value) == "NULL"


def as_number(value):
    """A number from a numeric or text value ('2.67', ' 2,67 '), else None."""
    if is_null(value):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        try:
            return float(str(value).strip().replace(",", "."))
        except ValueError:
            return None


def find_field(layer, name):
    for field in layer.fields():
        if field.name().lower() == name.lower():
            return field.name()
    raise ValueError(f"The layer has no '{name}' column")


def pivot_petrophysics(layer):
    if layer is None or not isinstance(layer, QgsVectorLayer):
        raise ValueError("Select the point layer in the Layers panel first")
    sample_f = find_field(layer, SAMPLE_FIELD)
    meas_f = find_field(layer, MEASUREMENT_FIELD)
    value_f = find_field(layer, VALUE_FIELD)
    other_fields = [f for f in layer.fields() if f.name() not in (sample_f, meas_f, value_f)]

    features = (
        layer.getSelectedFeatures() if layer.selectedFeatureCount() else layer.getFeatures()
    )
    samples = OrderedDict()       # SAMPLENO -> its geometry, attributes and values
    measurements = OrderedDict()  # measurement name -> list of all its values
    varying = set()               # other fields that differ within a sample
    for feature in features:
        sample = feature[sample_f]
        if is_null(sample):
            continue
        rec = samples.get(str(sample))
        attrs = {f.name(): feature[f.name()] for f in other_fields}
        if rec is None:
            rec = samples[str(sample)] = {
                "sample": sample,
                "geom": QgsGeometry(feature.geometry()),
                "attrs": attrs,
                "values": {},
            }
        else:
            varying.update(k for k, v in attrs.items() if v != rec["attrs"][k])
        meas = feature[meas_f]
        if is_null(meas) or not str(meas).strip():
            continue
        meas = str(meas).strip()
        value = feature[value_f]
        rec["values"].setdefault(meas, []).append(value)
        measurements.setdefault(meas, []).append(value)
    if not samples:
        raise ValueError("No features with a sample number")

    # output fields: the sample number, the per-sample constant columns, then one
    # <short>_val per measurement type (numeric if all its values are numbers)
    kept = [f for f in other_fields if f.name() not in varying]
    taken = {f.name().lower() for f in [layer.fields().field(sample_f)] + kept}
    columns = OrderedDict()  # measurement -> (field name, numeric?)
    for meas, values in measurements.items():
        base = short_name(meas) + "_val"
        name, n = base, 2
        while name.lower() in taken:
            name, n = f"{short_name(meas)}{n}_val", n + 1
        taken.add(name.lower())
        numeric = all(as_number(v) is not None for v in values if not is_null(v))
        columns[meas] = (name, numeric)

    out = QgsVectorLayer(
        QgsWkbTypes.displayString(layer.wkbType()), f"{layer.name()}_by_sample", "memory"
    )
    out.setCrs(layer.crs())
    provider = out.dataProvider()
    provider.addAttributes(
        [QgsField(layer.fields().field(sample_f))]
        + [QgsField(f) for f in kept]
        + [QgsField(name, DOUBLE if numeric else STRING) for name, numeric in columns.values()]
    )
    out.updateFields()

    combined = 0
    new_features = []
    for rec in samples.values():
        feature = QgsFeature(out.fields())
        feature.setGeometry(rec["geom"])
        feature[sample_f] = rec["sample"]
        for f in kept:
            feature[f.name()] = rec["attrs"][f.name()]
        for meas, values in rec["values"].items():
            name, numeric = columns[meas]
            if len(values) > 1:
                combined += 1
            if numeric:
                numbers = [as_number(v) for v in values if as_number(v) is not None]
                if numbers:
                    feature[name] = sum(numbers) / len(numbers)  # mean of repeats
            else:
                texts = list(OrderedDict.fromkeys(str(v) for v in values if not is_null(v)))
                if texts:
                    feature[name] = "; ".join(texts)
        new_features.append(feature)
    provider.addFeatures(new_features)
    out.updateExtents()
    QgsProject.instance().addMapLayer(out)

    print(f"{out.name()}: {len(new_features)} samples")
    for meas, (name, numeric) in columns.items():
        print(f"  {name:<20} <- {meas}{'' if numeric else '  (text: not all values are numbers)'}")
    if kept:
        print("  copied: " + ", ".join(f.name() for f in kept))
    if varying:
        print("  not copied (differ within a sample): " + ", ".join(sorted(varying)))
    if combined:
        print(f"  {combined} sample/measurement pairs had repeats: numbers averaged, text joined")
    return out


pivot_petrophysics(iface.activeLayer())
