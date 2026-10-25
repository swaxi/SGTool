[← Back to the documentation](README.md)

# Petrophys tab

Apparent susceptibility and density grids, and Henkel plots (density against log
susceptibility) of those grids and of measured samples.

## Apparent Petrophysics

Both calculations model the ground as a layer of vertical prisms whose top is a
depth *z* below the sensor, *t* thick, and divide the layer's response out of the
grid in the Fourier domain (*k* is the radial wavenumber in rad/m). They need a
grid in a metre-based projection, and can be previewed like the FFT filters.

**Apparent susceptibility** (output `<grid>_AppSus`, SI) from a magnetic grid
(nT) already reduced to the pole (or to the equator, when the result is only
approximate), magnetised by induction in an ambient field *F* (nT):

    T(k) = chi(k) * F/2 * exp(-k z) * (1 - exp(-k t))

With an infinite thickness (the default) the last term is dropped. With a finite
thickness the mean level cannot be recovered. As a magnetic anomaly grid (IGRF
removed) has no absolute level either, a background susceptibility (default
0.001 SI) is added, which keeps most values positive for Henkel plots. *F* is filled in
from the IGRF at the centre of the grid on the date set (or from the grid's own
metadata if it records the field intensity), and can be typed in.

**Apparent density** (output `<grid>_AppDens`, g/cc) from a Bouguer gravity grid
(mGal), with *G* = 6.674e-11:

    g(k) = 2 pi G drho(k) * exp(-k z) * (1 - exp(-k t)) / k

At *k* = 0 this is the Bouguer slab 2 pi G drho t, which sets the mean level; the
background density (default 2.67 g/cc) is then added.

Dividing by the response amplifies short wavelengths, so both apply a high-cut:
wavelengths longer than twice the cutoff pass, those at the cutoff or shorter are
removed, with a cosine roll-off between. Blank uses 4 x the cell size; 0 turns it
off.

| Setting | Default |
|---|---|
| Depth to top (m) | 0 (top at the sensor; for airborne data use the terrain clearance) |
| Thickness (m) | infinite (susceptibility), 1000 (density) |
| Background density (g/cc) | 2.67 |
| Background susceptibility (SI) | 0.001 |
| High-cut wavelength (m) | 4 x cell size |

## Henkel plots

Each plot is drawn as soon as its button is pressed, in a matplotlib window (whose
toolbar can save it);
with PNG and/or PDF ticked under Save Plots they are also saved, in the chosen
folder or beside the input layer.

If any susceptibility is zero or negative, all are shifted by -min + 1% of the
range so that the log axis works, and the shift is shown in the axis label.

- **Apparent Henkel Plot**: a heat map (cells per bin, 200 x 200 bins, log colour
  scale) of the chosen apparent density and susceptibility grids. Each cell of the
  coarser grid is paired with the nearest cell of the finer one (reprojected if the
  CRSs differ); cells where either grid has no data are left out. The axes span the
  0.1-99.9 percentiles so FFT edge effects don't squash the plot.
- **Henkel Plot per Category** (Spatial Analysis): the same grids clipped to the
  polygons of each category of a polygon layer (or of its selected polygons), one
  heat map per category with cells where both grids have data, in a grid of
  ceil(sqrt(N)) columns. All panels share their
  axes and colour scale, which shows the fraction of each category's cells per bin.
- **Henkel Plot of Samples** (True Petrophysics): measured density and
  susceptibility fields of a point layer (or of its selected points), each with a
  multiplier to convert it to g/cc or SI (e.g. 0.001 for kg/m3, 0.00001 for
  x10^-5 SI). As points, coloured and shaped by an optional category field (up to
  8 categories; with more, the 7 largest keep their own colour and the rest are
  grouped as Other), or as a heat map.
