[← Back to the documentation](README.md)

# Petrophysics (Petrophys tab)

The Petrophys tab (QGIS and ArcGIS Pro), between Grid + Wavelets and Utils, turns magnetic and gravity
grids into **apparent susceptibility** and **apparent density** grids, and plots
them against each other as **Henkel plots** (density against log susceptibility),
for the whole area or for each unit of a geology map. Measured samples (any two
properties of a point layer) can be plotted the same way, so apparent and measured
properties can be compared. With only one property chosen, any of the plots
becomes a histogram.

The tab has four boxes:

| Box | What it does |
|---|---|
| **Apparent Petrophysics** | Calculates an apparent susceptibility or apparent density grid from the grid selected at the top of the tab (Apply Processing, with live Preview) |
| **Apparent Henkel Plot / Histogram** | Plots an apparent density grid against an apparent susceptibility grid, for the whole area (Plot Apparent) or for each category of a polygon layer (Plot per Category) |
| **Measured Petrophysics** | Plots any two fields of a point layer of samples against each other (Plot Samples) |
| **Save Plots** | Also saves each plot as PNG and/or PDF |

## Workflow

1. **Magnetics**: reduce the TMI anomaly grid (IGRF removed) to the pole with
   RTE/P on the FFT Filters tab (or to the equator at low latitudes).
2. On the Petrophys tab, select the `_RTP` grid, choose **Apparent
   Susceptibility**, check the field strength (filled in from the IGRF; see
   below), and press **Apply Processing**. The result is saved as
   `<grid>_AppSus.tif`.
3. **Gravity**: select the Bouguer anomaly grid, choose **Apparent Density**,
   set **Grid units** to mGal or µm/s², and press **Apply Processing** to make
   `<grid>_AppDens.tif`.
4. In **Apparent Henkel Plot / Histogram**, pick the two new grids and press
   **Plot Apparent**. Pick a geology polygon layer and its lithology field and
   press **Plot per Category** to compare units.
5. In **Measured Petrophysics**, pick a point layer of samples with density as X
   and susceptibility as Y (Log ticked) and press **Plot Samples**.

Tick **Preview** beside Apply Processing to see either calculation update live as
the settings change.

## Apparent Petrophysics

Apparent Susceptibility and Apparent Density are tick boxes, at most one ticked at
a time (ticking one unticks the other; both start unticked, so arriving on the tab
with Preview on calculates nothing until one is ticked). **Apply Processing** on
this tab runs the ticked one only, and Apply on the other tabs never runs them.
Both need a grid in a metre-based projection.

Both model the ground as a layer of vertical prisms whose top is a depth *z* below
the sensor and which is *t* thick, and divide the layer's response out of the grid
in the Fourier domain (*k* is the radial wavenumber in rad/m).

### Apparent susceptibility

From a magnetic anomaly grid (nT) already reduced to the pole (or to the equator,
when the result is only approximate), magnetised by induction in an ambient field
of strength *F* (nT):

    T(k) = chi(k) * F/2 * exp(-k z) * (1 - exp(-k t))

With an infinite thickness (the default) the last term is dropped. The result is
in SI units.

- **Input**: the anomaly after IGRF removal and RTP. Using the TMI anomaly without
  RTP partly cancels the positive and negative lobes of each anomaly.
- **Field (nT)**: the full IGRF total intensity (about 58,000 nT in Western
  Australia), not anything taken from the anomaly grid. It is filled in from the
  IGRF at the centre of the grid on the date shown whenever a grid is selected, the
  date changes or **IGRF** is pressed, or from the grid's own metadata if it
  records the field intensity. It can also be typed in.

### Apparent density

From a Bouguer gravity grid, with *G* = 6.674e-11:

    g(k) = 2 pi G drho(k) * exp(-k z) * (1 - exp(-k t)) / k

The result is in g/cc.

- **Grid units**: mGal or µm/s² (gravity units; 1 mGal = 10 µm/s²). Geoscience
  Australia's national grids are in µm/s². A Bouguer grid of the Yilgarn spanning
  hundreds of units is almost certainly in µm/s².
- **Thickness**: thin layers exaggerate the density contrasts needed to explain
  long-wavelength anomalies; 5000 m is the default.

### Levels, backgrounds and the high-cut

The absolute level of either grid (its datum, the regional field, the IGRF
removal) carries no density or susceptibility information, so the *k* = 0 term is
not used. Each result is a zero-mean variation to which a background is added:
2.67 g/cc for density and 0.001 SI for susceptibility. The susceptibility
background keeps most values positive for the log axis of Henkel plots; raise it
if a plot reports many values left out.

Dividing by the response amplifies short wavelengths, so both calculations apply
a high-cut: wavelengths longer than twice the cutoff pass, those at the cutoff or
shorter are removed, with a cosine roll-off between. Blank uses 4 x the cell
size; 0 turns it off.

| Setting | Default |
|---|---|
| Depth to top (m) | 0 (top at the sensor; for airborne data use the terrain clearance) |
| Thickness (m) | infinite (susceptibility; tick Thickness for a finite layer), 5000 (density) |
| Gravity grid units | mGal |
| Background susceptibility (SI) | 0.001 |
| Background density (g/cc) | 2.67 |
| High-cut wavelength (m) | 4 x cell size |

Each output records its settings in its provenance metadata (Utils tab, Read
metadata).

## Henkel plots, scatter plots and histograms

Each plot is drawn as soon as its button is pressed, in its own matplotlib window
(whose toolbar can also save it). With PNG and/or PDF ticked under **Save Plots**
it is saved too, in the chosen folder or, if that is blank, beside the input layer.

On a log axis, values of zero or less can't be shown, so they are left out and the
plot title says how many. Heat maps count points per bin on a log colour scale
(viridis), with empty bins left white.

### Apparent Henkel Plot / Histogram

- **Plot Apparent**: a heat map (200 x 200 bins) of the chosen apparent density
  grid (linear x axis) against the apparent susceptibility grid (log y axis). Each
  cell of the coarser grid (usually gravity) is paired with the nearest cell of the
  finer grid, reprojected if the CRSs differ; cells where either grid has no data
  are left out. The axes span the 0.1-99.9 percentiles so FFT edge effects don't
  squash the plot.
- **Plot per Category**: the same grids clipped to the polygons of each category
  of a polygon layer (or of its selected polygons), e.g. the lithologies of a
  geology map. Categories with cells where both grids have data get a heat map
  each, in a grid of ceil(sqrt(N)) columns. All panels share their axes and colour
  scale, which shows the fraction of each category's cells per bin, so small and
  large units compare fairly. With **Colour by position** ticked, each bin is
  coloured by where it is on the plot instead (red to green along the density
  axis, red to blue up the susceptibility axis: low density and susceptibility
  red, dense green, magnetic blue, both cyan), and the fraction per bin is shown
  by opacity, solid where most crowded and faint where sparse.

### Measured Petrophysics

**Plot Samples** plots the **Y field** against the **X field** of a point layer
(or of its selected points). Any measured properties can be used; with density as
X and susceptibility as Y (log) it is a Henkel plot.

- The axes are named after the fields, and the category legend after the category
  field.
- Fields can be numeric, or text holding numbers (`2.67`, ` 2,67 `); values that
  aren't numbers, and NULLs, are skipped.
- Each field has a **multiplier** (e.g. 0.001 for kg/m3 to g/cc, 0.00001 for
  x10^-5 SI to SI) and its own **Log** tickbox (X linear and Y log by default).
- **Points** plots every sample. With a **category field**, each category gets a
  colour and marker: 8 colours x 8 markers give up to 64 categories, assigned by
  sample count (the 8 largest are circles in the 8 colours, the next 8 squares, and
  so on); any further categories are grouped as Other. The legend lists the
  categories by name, with their sample counts.
- **Heat map** shows samples per bin instead, for large datasets (it ignores the
  category).

### One variable only: histograms

Leave either apparent grid, or either measured field, blank to plot a histogram of
the other instead:

- Apparent density has linear bins and apparent susceptibility log bins; a measured
  field follows its Log tickbox.
- With a category (polygons for the apparent grids, a category field for samples)
  there is one histogram per category, in a grid of panels sharing their bins and
  axes, each showing the fraction of its category's values per bin.
- Every histogram marks its **median** with a red line and gives the median, mean
  and standard deviation of its values in the top right corner. The vertical axis
  stops at the tallest bar of all the panels, rounded up.

## Other ways to run it

- **QGIS Processing**: *Apparent susceptibility* and *Apparent density* algorithms in the
  SGTool Petrophysics group, with the same settings (thickness 0 = infinite, high-cut
  0 = 4 x the cell size and -1 = none).
- **ArcGIS Pro**: the same tab in the SGTool GUI, and **Apparent Susceptibility** and
  **Apparent Density** tools in the toolbox's Petrophysics category. In the GUI the
  polygons and sample points are chosen as files, and all their features are used. See
  [ArcGIS Pro differences](arcgis-pro.md).
- **Apply same steps to another grid** repeats both calculations; the field strength of
  apparent susceptibility is used as recorded, not recalculated for the new grid.

### Grouping similar categories

Tick **Group similar** (next to Plot per Category, and next to Plot Samples) to group
categories whose distributions are alike and label them **A, B, C ...**:

- Each category's fraction of cells (or samples) per bin, on the plotted bins and
  lightly smoothed, is compared with every other category's using the
  Jensen-Shannon distance (0 = identical, 1 = no overlap), and the categories are
  grouped by average-linkage hierarchical clustering.
- The box beside it sets the number of groups: **Auto** cuts the clustering at its
  largest jump (one group if all the categories are alike), or choose 2 to 8.
- Groups are lettered in order of their mean position along the x axis (for
  Henkel plots, A is the lightest).
- Per-category panels (heat maps and histograms) are ordered by group, with the
  letter in the title and a frame in the group's colour (histogram bars too). In
  the sample plot the letter is shown in the legend, which is sorted by group.
- Categories with fewer than 30 cells (10 samples) are marked **-** and not grouped.

A similar distribution does not mean the same geology, and apparent values blur
across boundaries, so small polygons tend to group with their neighbours: use the
groups as a prompt for interpretation, not a classification.

## Preparing sample data

Petrophysical databases often hold one row per measurement (sample number,
property, value) rather than one row per sample. The script
[`scripts/pivot_petrophysics.py`](../scripts/pivot_petrophysics.py), run in the
QGIS Python console with the layer selected, makes a temporary layer with one
point per sample and a `<property>_val` field for each property (repeats of a
property on one sample are averaged), copying columns that are the same on every
row of a sample (such as lithology) for use as the category field. The column
names it reads (`SAMPLENO`, `PROPERTY`, `Value`) and the short property names are
set at the top of the script.
