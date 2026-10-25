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

With an infinite thickness (the default) the last term is dropped. *F* is filled
in from the IGRF at the centre of the grid on the date set (or from the grid's own
metadata if it records the field intensity), and can be typed in. The input should
be the anomaly (IGRF removed) after RTP, and *F* the full IGRF total intensity.

**Apparent density** (output `<grid>_AppDens`, g/cc) from a Bouguer gravity grid
(in mGal or µm/s², chosen in the dialog), with *G* = 6.674e-11:

    g(k) = 2 pi G drho(k) * exp(-k z) * (1 - exp(-k t)) / k

The absolute level of either grid (its datum, the regional field, the IGRF
removal) carries no density or susceptibility information, so the *k* = 0 term is
not used: each result has a zero-mean variation to which a background is added,
2.67 g/cc for density and 0.001 SI for susceptibility (which keeps most values
positive for Henkel plots).

Dividing by the response amplifies short wavelengths, so both apply a high-cut:
wavelengths longer than twice the cutoff pass, those at the cutoff or shorter are
removed, with a cosine roll-off between. Blank uses 4 x the cell size; 0 turns it
off.

| Setting | Default |
|---|---|
| Depth to top (m) | 0 (top at the sensor; for airborne data use the terrain clearance) |
| Thickness (m) | infinite (susceptibility), 5000 (density) |
| Gravity grid units | mGal |
| Background density (g/cc) | 2.67 |
| Background susceptibility (SI) | 0.001 |
| High-cut wavelength (m) | 4 x cell size |

## Henkel plots

Each plot is drawn as soon as its button is pressed, in a matplotlib window (whose
toolbar can save it);
with PNG and/or PDF ticked under Save Plots they are also saved, in the chosen
folder or beside the input layer.

Susceptibilities of zero or less can't be shown on the log axis, so they are left
out; the plot title says how many.

- **Apparent Henkel Plot**: a heat map (cells per bin, 200 x 200 bins, log colour
  scale) of the chosen apparent density and susceptibility grids. Each cell of the
  coarser grid is paired with the nearest cell of the finer one (reprojected if the
  CRSs differ); cells where either grid has no data are left out. The axes span the
  0.1-99.9 percentiles so FFT edge effects don't squash the plot.
- **Plot per Category** (in the same Apparent Henkel Plot box): the same grids clipped to the
  polygons of each category of a polygon layer (or of its selected polygons), one
  heat map per category with cells where both grids have data, in a grid of
  ceil(sqrt(N)) columns. All panels share their
  axes and colour scale, which shows the fraction of each category's cells per bin.
- **Plot Samples** (True Petrophysics): measured density and
  susceptibility fields of a point layer (or of its selected points), each with a
  multiplier to convert it to g/cc or SI (e.g. 0.001 for kg/m3, 0.00001 for
  x10^-5 SI). As points, coloured and shaped by an optional category field (up to
  8 categories; with more, the 7 largest keep their own colour and the rest are
  grouped as Other), or as a heat map.
