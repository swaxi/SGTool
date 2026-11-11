[← Back to the documentation](README.md)

# SGTool reference

This reference describes every tab, control and filter of the SGTool panel in QGIS, with its settings, defaults and output, and an example of each filter applied to real data. For step-by-step instructions see the [user guide](user-guide.md). For the formulas and the published methods behind each filter, follow the **Method** links to [FFT filters](filters-fft.md), [Convolution filters, spatial statistics, gridding and utilities](filters-convolution-and-statistics.md) and [Petrophysics](petrophysics.md).

**Contents**

- [Example data and figures](#example-data-and-figures)
- [Controls common to all tabs](#controls-common-to-all-tabs)
- [Output names](#output-names)
- [FFT Filters tab](#fft-filters-tab): [Grav/Mag](#gravmag-filters), [Frequency](#frequency-filters), [Gradient](#gradient-filters)
- [Conv + Stats tab](#conv--stats-tab): [Convolution](#convolution-filters), [Spatial statistics](#spatial-statistics), [Multivariate](#multivariate-statistical-analysis), [Euler](#euler-deconvolution)
- [Grid + Wavelets tab](#grid--wavelets-tab)
- [Petrophys tab](#petrophys-tab)
- [Utils tab](#utils-tab)
- [Help tab](#help-tab)
- [Processing algorithm names](#processing-algorithm-names)

---

## Example data and figures

The examples use `ogrm_usgs_mag_tmi.tif`, an aeromagnetic total magnetic intensity grid of part of Mauritania (courtesy of the Mauritanian Government (OMRG) and the USGS), linked from the plugin's Help tab: 924 x 629 cells of 175 m, WGS 84 / UTM zone 29N, -1410 to 5005 nT. See [The example data](user-guide.md#1-the-example-data) in the user guide.

- Each figure shows the input on the left and the result on the right, with the settings in the title. The results were calculated with SGTool's own filter code (`sgt_filters.py`, also used by the dialog and the Processing algorithms).
- Magnetic reductions use the IGRF field for the dialog's default date, 1 January 2000, at the grid centre: inclination 29.2°, declination -5.4°, total field 36,751 nT.
- Unless stated otherwise, filters are applied to the RTP grid, `ogrm_usgs_mag_tmi_RTP`, and window statistics to its first vertical derivative.
- Figures marked *zoomed* show the 38 x 56 km area outlined below, where the effect of filters that act over a few cells can be seen.
- Grids such as TMI and RTP are drawn with a histogram-equalised colour scale. Derivatives and other zero-centred results use a blue-white-red scale centred on zero.

![The example grid](images/examples/example_grid.png)

---

## Controls common to all tabs

| Control | Where | What it does |
|---|---|---|
| **Select Grid** | Top of FFT Filters, Conv + Stats, Petrophys and Utils | The grid the tab's filters are applied to. Each tab has its own. Lists the raster layers in the project. |
| **Load Grid** and **...** | FFT Filters | Loads a grid from a file: GeoTIFF, Geosoft `.grd`, ER Mapper `.ers`, Noddy `.grv` and `.mag`. A `.grd` takes its CRS from its `.xml` file, or EPSG:4326 if there is none, and is saved as a GeoTIFF. Other formats can be loaded with QGIS's own raster loader. |
| **Apply Processing** | Top of each tab | Applies each ticked filter to the selected grid, as a background calculation. Each filter gives its own output; the filters are unticked afterwards. While **Preview** is ticked the button reads **Apply Preview** and saves the previewed filter at full resolution. |
| **Apply to Window** | Beside Apply | While previewing: clips the grid to the map window (as **Clip to Window**), selects the clip and applies the previewed filter to it. |
| **Cancel** | Beside Apply | Stops a running calculation. Nothing is saved. See [Background Calculations](background-calculations.md). |
| **Preview** | Beside Cancel | Shows one filter on the area in the map window in a temporary layer, updated as settings change. Filters that cannot be previewed are disabled while it is ticked. See [Live Preview](live-preview.md). |
| **Max FFT Buffer (pixels)** | Bottom of FFT Filters | Upper limit of the padding added on each side of the grid before an FFT filter. Default 5000. The padding is the grid's smaller dimension, up to this value. |
| **Version / Units / CRS / Pixels / Pixel size** | Bottom of each tab | The plugin version and, for the selected grid, its length units (m or deg), CRS, size in cells and cell size. |
| **●** | Before some filters | Filters marked with a dot should be applied to magnetic data reduced to the pole or equator (or to gravity data). |

---

## Output names

Each result is saved as a GeoTIFF beside the input grid, named after it with a suffix, and added to the project. Apart from the derivative's power and direction, the settings are not in the name; they are recorded in the file's provenance metadata. A name already in use gets `_1`, `_2` and so on.

| Filter | Suffix | Filter | Suffix |
|---|---|---|---|
| Reduction to the pole | `_RTP` | Mean | `_Mn` |
| Low-latitude RTP | `_RTP_LL` | Median | `_Md` |
| Differential RTP | `_DRTP` (+ `_DRTP_inc`, `_DRTP_dec`) | Gaussian | `_Gs` |
| Reduction to the equator | `_RTE` | Directional | `_Dr` |
| Vertical integration | `_VI` | Sun shading | `_Sh` |
| Upward / downward continuation | `_UC` / `_DC` | Min / Max | `_SS_Min` / `_SS_Max` |
| Line-noise removal | `_DirC` (noise estimate: `_DirC_noise`) | Variance / Std deviation | `_SS_Var` / `_SS_StdDev` |
| Remove regional | `_RR` | Skewness / Kurtosis | `_SS_Skew` / `_SS_Kurt` |
| Band pass | `_BP` | Local anisotropy | `_SS_AnisoMag`, `_SS_AnisoOrient` |
| High pass / low pass | `_HP` / `_LP` | MRVBF | `_MRVBF`, `_MRRTF`, `_slope` + composite |
| Automatic gain control | `_AGC` | PCA / ICA | `_PCA` / `_ICA` |
| Derivative | `_d<power><direction>`, e.g. `_d1z`, `_d0.5x` | Euler deconvolution | `_estimates_SI_<n>.txt` |
| Tilt angle | `_TA` | Threshold to NaN | `_Clean` |
| Analytic signal | `_AS` | Clipping polygon | `_boundary` (shapefile) |
| Total horizontal gradient | `_THG` | Clip to window | `_clip` |
| Apparent susceptibility | `_AppSus` | Worms | `_worms.csv` (+ `_worms.shp`) |
| Apparent density | `_AppDens` | LUT to grey scale | `_gray` |

---

## FFT Filters tab

![The FFT Filters tab](images/ui/tab_fft_filters.png)

All the filters on this tab work in the Fourier domain, except Remove Regional and AGC. Before the transform, the grid's median is removed, no-data cells are filled from the nearest valid cell, and the grid is padded and tapered at its edges; afterwards the padding is removed and the no-data cells restored. See [How the FFT filters are applied](filters-fft.md#how-the-fft-filters-are-applied).

### Grav/Mag Filters

#### RTE/P (reduction to the pole or equator)

Tick **RTE/P** and choose the method in the list.

| Setting | Default | Meaning |
|---|---|---|
| Method list | RTP | **RTP**, **Low-lat. RTP**, **Diff. RTP** or **RTE** (below) |
| **IGRF** | | Fills in Inc, Dec and Int from the IGRF model at the centre of the selected grid on the date shown. For a Noddy grid the values are read from the file. |
| **Inc** | 0 | Magnetic inclination, degrees from horizontal (positive down) |
| **Dec** | 0 | Magnetic declination, degrees clockwise from north |
| **Int** | 50000 | Total field strength, nT |
| Date | 1 Jan 2000 | Survey date, 1900 to 2030, used by **IGRF** and by Diff. RTP |
| **Iterations** (Low-lat. RTP only) | blank (auto) | Number of iterations of the equator-to-pole operator; blank chooses the one whose result best matches the analytic signal |
| **Stripe suppression (°)** (Low-lat. RTP only) | 10 | The iterations fade out within this angle of the direction perpendicular to the declination; 0 turns it off |

| Method | When to use | Output |
|---|---|---|
| **RTP** | Inclination steeper than about 20° | `_RTP` |
| **Low-lat. RTP** | Inclination shallower than about 20°: a stable RTP (Hao et al. 2018) | `_RTP_LL` |
| **Diff. RTP** | Large surveys where the field direction varies (Cooper & Cowan 2005). Inc and Dec are calculated from the IGRF at the grid's corners and centre for the date shown, so the Inc and Dec boxes are not used. Cannot be previewed. | `_DRTP`, `_DRTP_inc`, `_DRTP_dec` |
| **RTE** | Low inclinations, as an alternative to low-latitude RTP | `_RTE` |

Input: a TMI anomaly grid (IGRF removed). Method: [Reduction to the Pole, Equator, Low-latitude and Variable RTP](filters-fft.md#gravmag-filters).

![RTP settings](images/ui/step_rtp.png)

**Example: RTP** (Inc 29.2, Dec -5.4). Dipolar anomalies become single anomalies over their sources. At this moderate inclination, RTP also amplifies features that run along the declination, giving the north-south streaks in the south.

![RTP](images/examples/rtp.png)

**Example: low-latitude RTP** (Inc 29.2, Dec -5.4, automatic iterations, stripe suppression 10°). Fewer streaks along the declination than the RTP.

![Low-latitude RTP](images/examples/rtp_ll.png)

**Example: differential RTP.** Compared with the ordinary RTP. The example grid spans only about 1.6° of longitude and 1° of latitude (inclination 28.2° to 30.3° at the corners), so the two are almost the same; differential RTP matters for surveys covering several degrees.

![Differential RTP](images/examples/drtp.png)

**Example: RTE** (Inc 29.2, Dec -5.4). At the equator, an induced anomaly is a low over its source, elongated perpendicular to the declination.

![RTE](images/examples/rte.png)

#### Vertical Integration ●

Integrates the field vertically (divides the spectrum by the wavenumber). Applied to an RTP or RTE grid it gives a result proportional to the **pseudogravity**: the gravity anomaly the sources would produce if their density were proportional to their magnetisation. Also useful for stitched grids with very different line spacings. No settings. Needs a metre-based projection. Output `_VI`. Method: [Vertical Integration](filters-fft.md#gravmag-filters).

**Example** (of the RTP). Anomalies broaden and the long wavelengths dominate.

![Vertical integration](images/examples/vertical_integration.png)

#### Continuation ●

| Setting | Default | Meaning |
|---|---|---|
| **Direction** | up | **up**: the field at a greater height; **down**: at a smaller height |
| **Height (m)** | 500 | Distance to continue, always in metres (converted approximately for a geographic grid) |

Upward continuation smooths the field as if measured higher up, suppressing shallow sources. Downward continuation sharpens it, but amplifies noise very quickly: keep the height small compared with the source depths. Output `_UC` or `_DC`. Method: [Continuation](filters-fft.md#gravmag-filters).

**Example: upward continuation, 2000 m.**

![Upward continuation](images/examples/continuation_up.png)

**Example: downward continuation, 100 m** (zoomed). Even 100 m (less than one cell) sharpens the anomalies noticeably and starts to amplify noise.

![Downward continuation](images/examples/continuation_down.png)

### Frequency Filters

#### Radial Power Spectrum

Opens a window with the selected grid and its radially averaged power spectrum: ln(power) against wavenumber, in radians per cell (a wavenumber *k* is a wavelength of 2π/*k* cells). Straight segments of the curve correspond roughly to groups of sources at similar depths, and the breaks between them suggest cutoffs for the frequency filters. No-data cells are set to 0 first. Needs matplotlib. Method: [Radially averaged power spectrum](filters-fft.md#frequency-filters).

![Radial power spectrum of the example TMI](images/examples/power_spectrum.png)

#### Direction Cos/Butterw. (line-noise removal)

Removes stripes parallel to the flight lines. A noise estimate is made from the wavelengths between 2 times the smallest and 10 times the largest line spacing, within a wedge of directions about the azimuth, scaled and subtracted.

| Setting | Default | Meaning |
|---|---|---|
| **Azimuth** | 0 | Direction of the flight lines (and so of the stripes), degrees clockwise from north |
| **Scale** | 1 | Multiplier applied to the noise estimate before subtracting it. Values of 2 or more over-subtract. |
| **Wedge** | 45 | Half-width of the wedge of directions, degrees. Narrower removes less geology that runs parallel to the lines. |
| **Line Spacing min** | 1000 | Smallest line spacing, in grid units |
| **Line Spacing max** | blank | Largest line spacing; leave blank if the spacing is constant |

Output `_DirC`. The QGIS Processing algorithm *Remove line noise* can output the noise estimate instead (`_DirC_noise`). Method: [Directional line-noise removal](filters-fft.md#frequency-filters).

**Example** (zoomed). The example grid has little line noise, so synthetic north-south line noise (a random level error on lines 700 m apart, plus a slow drift along each line) was added to the TMI. Both grids are shown as their first vertical derivative, where stripes are easiest to see. Settings: Azimuth 0, Line Spacing 700, Scale 1, Wedge 10.

![Line-noise removal](images/examples/line_noise.png)

The noise estimate for the same settings (Output: Noise estimate). It should look like stripes along the lines, not like geology.

![Line-noise estimate](images/examples/line_noise_estimate.png)

#### Remove Regional

| Setting | Default | Meaning |
|---|---|---|
| **1st order fit** / **2nd order fit** | 1st order | Subtract a least-squares plane (1st order) or quadratic surface (2nd order) |

Output `_RR`. Method: [Remove Regional](filters-fft.md#frequency-filters).

**Example: 1st order** (of the TMI).

![Remove regional](images/examples/remove_regional.png)

#### Band Pass ●

| Setting | Default | Meaning |
|---|---|---|
| **Low** | 50000 | The *longest* wavelength kept, in grid units |
| **High (proj units)** | 5000 | The *shortest* wavelength kept, in grid units |
| **Width** | 5000 | Width of the cosine transition at each cutoff; 0 is a sharp cut |

Low is the larger number. Output `_BP`. Method: [Band Pass Filter](filters-fft.md#frequency-filters).

**Example: Low 20000, High 2000, Width 1000.**

![Band pass](images/examples/band_pass.png)

#### High/Low Pass ●

| Setting | Default | Meaning |
|---|---|---|
| **Accept** | Low | **High**: keep wavelengths shorter than the cutoff; **Low**: keep wavelengths longer than the cutoff |
| **Cutoff (proj units)** | 5000 | Cutoff wavelength, in grid units |
| **Width** | 5000 | Width of the cosine transition; 0 is a sharp cut |

Output `_HP` or `_LP`. Method: [High and Low Pass](filters-fft.md#frequency-filters).

**Example: high pass, cutoff 10000, width 2000.** The residual: shallow sources.

![High pass](images/examples/high_pass.png)

**Example: low pass, cutoff 10000, width 2000.** The regional: broad, deeper sources.

![Low pass](images/examples/low_pass.png)

#### Automatic Gain Control ●

| Setting | Default | Meaning |
|---|---|---|
| **Window (pixels)** | 10 | Standard deviation, in cells, of the Gaussian used to measure the local amplitude |

Divides the grid by its local amplitude (the Gaussian-smoothed absolute value), so weak and strong anomalies show equally. The result is a ratio with no physical units. Output `_AGC`. Method: [Automatic Gain Control](filters-fft.md#frequency-filters).

**Example: window 10.**

![AGC](images/examples/agc.png)

### Gradient Filters

#### Derivative ●

| Setting | Default | Meaning |
|---|---|---|
| **Direction** | z | **z** (vertical), **x** (east) or **y** (north) |
| **Power** | 1 | Order of the derivative; fractional values are allowed (e.g. 0.5) |

Output `_d<power><direction>`, for example `_d1z` (first vertical derivative, 1VD). The result is zero-mean. Method: [Derivative](filters-fft.md#gradient-filters).

**Example: first vertical derivative (z, power 1).** Shallow sources and edges stand out; noise is amplified too.

![First vertical derivative](images/examples/derivative_1vd.png)

**Example: x derivative, power 1.** Edges running north-south are enhanced.

![x derivative](images/examples/derivative_x.png)

**Example: vertical derivative, power 0.5.** Between the RTP and the 1VD.

![Half vertical derivative](images/examples/derivative_half_vd.png)

#### Total Horizontal Gradient ●

The amplitude of the horizontal gradient. On RTP or gravity data its maxima lie over the edges of sources. No settings. Output `_THG`. Method: [Total Horizontal Gradient](filters-fft.md#gradient-filters).

![Total horizontal gradient](images/examples/thg.png)

#### Tilt Angle ●

The angle (in degrees, -90 to 90) between the vertical derivative and the total horizontal gradient. It balances strong and weak anomalies; the 0° contour lies roughly over the edges of vertical sources. Tends to over-connect features. No settings. Output `_TA`. Method: [Tilt Angle](filters-fft.md#gradient-filters).

![Tilt angle](images/examples/tilt.png)

#### Analytic Signal

The amplitude of the three gradients. Its maxima lie over the sources, almost independently of the field direction, so it can be used on TMI data at low latitudes. No settings. Output `_AS`. Method: [Analytic Signal](filters-fft.md#gradient-filters).

**Example** (of the TMI).

![Analytic signal](images/examples/analytic_signal.png)

---

## Conv + Stats tab

![The Conv + Stats tab](images/ui/tab_conv_stats.png)

The tab has its own **Select Grid**. Filter sizes and windows are in cells.

### Convolution Filters

No-data cells are left out of each window. Method: [Convolution Filters](filters-convolution-and-statistics.md#convolution-filters).

| Filter | Setting | Default | Output |
|---|---|---|---|
| **Mean** | Filter size (odd number of cells) | 3 | `_Mn` |
| **Median** | Filter size (odd number of cells) | 3 | `_Md` |
| **Gaussian** | Sigma (cells) | 1 | `_Gs` |
| **Directional** | Direction: N, NE, E, SE, S, SW, W, NW | N | `_Dr` |
| **Sun Shading** | Azimuth, Zenith (degrees), relief | -45, 45, ticked | `_Sh` |

**Mean** and **Median** average or take the median of the cells in an n x n window; the median removes spikes and keeps edges. **Gaussian** smooths with a Gaussian kernel. **Directional** applies a 3 x 3 compass-gradient kernel that enhances edges facing the chosen direction. **Sun Shading** lights the grid from the sun's **Azimuth** (degrees clockwise from north) at **Zenith** degrees above the horizon; with **relief** ticked it uses a softer shading translated from GRASS GIS `r.relief`.

**Example: mean, 9 x 9** (zoomed).

![Mean filter](images/examples/mean.png)

**Example: median, 9 x 9** (zoomed). Edges stay sharper than with the mean.

![Median filter](images/examples/median.png)

**Example: Gaussian, sigma 3** (zoomed).

![Gaussian filter](images/examples/gaussian.png)

**Example: directional, NE** (zoomed).

![Directional filter](images/examples/directional.png)

**Example: sun shading, azimuth 45, zenith 45, relief.**

![Sun shading](images/examples/sun_shading.png)

### Spatial Statistics

| Setting | Default | Meaning |
|---|---|---|
| **Window size for all filters** | 3 | Diameter of the circular window, in cells (odd) |
| Window type (beside Local Anisotropy) | Gaussian | Gaussian or box smoothing of the gradient products |

Tick any of the statistics below. Each is saved as its own grid. A window needs at least 3 valid values. Standard deviation and variance are slow on large grids and windows, skewness and kurtosis very slow. Method: [Spatial Statistics](filters-convolution-and-statistics.md#spatial-statistics).

| Statistic | Output | What it shows |
|---|---|---|
| **Min**, **Max** | `_SS_Min`, `_SS_Max` | Smallest and largest value in the window |
| **Variance**, **Standard Deviation** | `_SS_Var`, `_SS_StdDev` | How much the values vary in the window: texture or roughness |
| **Skewness** | `_SS_Skew` | Asymmetry of the values in the window |
| **Kurtosis** | `_SS_Kurt` | How peaked the values are (excess kurtosis, 0 for a normal distribution) |
| **Local Anisotropy** | `_SS_AnisoMag`, `_SS_AnisoOrient` | How linear the grid is (0 to 1) and the strike of the linear features (0 to 180°) |

The examples are of the RTP 1VD, with a window of 9 cells (15 for anisotropy).

**Standard deviation** (zoomed).

![Standard deviation](images/examples/stats_stddev.png)

**Variance** (zoomed).

![Variance](images/examples/stats_var.png)

**Min** (zoomed).

![Min](images/examples/stats_min.png)

**Max** (zoomed).

![Max](images/examples/stats_max.png)

**Skewness** (zoomed).

![Skewness](images/examples/stats_skew.png)

**Kurtosis** (zoomed).

![Kurtosis](images/examples/stats_kurt.png)

**Local anisotropy magnitude**, window 15.

![Anisotropy magnitude](images/examples/anisotropy_magnitude.png)

**Local anisotropy orientation**, window 15. The colour scale is cyclic, as 0° and 180° are the same strike.

![Anisotropy orientation](images/examples/anisotropy_orientation.png)

#### MRVBF / MRRTF / Slope

Multiresolution valley bottom flatness (MRVBF) and ridge top flatness (MRRTF) of a digital elevation model (Gallant & Dowling 2003), translated from SAGA GIS. Writes `_MRVBF`, `_MRRTF` and `_slope` grids and a three-band (MRRTF, MRVBF, slope) composite. Needs a DEM in a metre-based projection. Cannot be previewed. Method: [MRRTF / MRVBF / Slope](filters-convolution-and-statistics.md#spatial-statistics).

| Setting | Default | Meaning |
|---|---|---|
| **Slope Threshold** | 16 | Initial slope threshold (percent) at the finest resolution |
| **Pctl Lowness** | 0.40 | Threshold for the lowness percentile (valley bottoms) |
| **Pctl Upness** | 0.35 | Threshold for the upness percentile (ridge tops) |
| **Slope Shape** | 4 | Shape parameter of the slope transformation |
| **Pctl Shape** | 3 | Shape parameter of the elevation percentile transformation |
| **Max Res %** | 100 | Coarsest resolution, as a percentage of the grid diagonal |

### Multivariate Statistical Analysis

| Filter | Setting | Default | Output |
|---|---|---|---|
| **Principal Component Analysis** | Number of Components (0 = all) | 0 | `_PCA` |
| **Independent Component Analysis** | Number of Components (0 = all) | 0 | `_ICA` |

Both work on the bands of a **multiband** grid (for example K, Th and U radiometrics stacked into one file). The bands are standardised first, and cells that are no-data in any band are left out. Each component is one band of the output. Need scikit-learn. Cannot be previewed. Method: [Multivariate Statistical Analysis](filters-convolution-and-statistics.md#multivariate-statistical-analysis).

### Euler Deconvolution

| Setting | Default | Meaning |
|---|---|---|
| **Euler Deconvolution** | | Tick to run |
| **Stats** | unticked | Also save statistics of the solutions in grid windows (`_window_stats_SI_<n>.txt`, `_window_summary.txt`) |
| **Threshold** | 0.1 | Fraction of the solutions kept (0.1 = the best 10%) |
| **Window (odd)** | 5 | Window size, in cells |

Estimates source positions and depths with the Reliable Euler method of Melo & Barbosa (2020), for structural indices 0 (contact), 1 (dyke or sill), 2 (pipe) and 3 (sphere). The solutions for each are saved beside the grid as `<grid>_estimates_SI_0.txt` to `_SI_3.txt`, with the columns `y_source` (northing), `x_source` (easting), `z_source` (depth), `base_level` and `std_dfdz` (the standard deviation of the vertical derivative in the window). Needs a metre-based projection. Cannot be previewed. Method: [Euler Deconvolution](filters-convolution-and-statistics.md#euler-deconvolution).

![Euler settings](images/ui/step_euler.png)

**Example: SI 1, window 5, threshold 0.1, of the RTP.** The solutions are drawn over the RTP 1VD, coloured by depth (depths between 0 and 5000 m shown).

![Euler solutions](images/examples/euler.png)

---

## Grid + Wavelets tab

![The Grid + Wavelets tab](images/ui/tab_grid_wavelets.png)

**Cancel running calculation** at the top stops gridding or worms.

### Import point or line data

| Setting | Meaning |
|---|---|
| **Point/Line data** and **...** | The file to import: CSV, ASEG-GDF2 `.dat`, or `.xyz` line data |
| CRS | The coordinate system of the data |
| **long_x**, **lat_y** | The columns holding the x (longitude) and y (latitude) coordinates |
| **Load Tie Lines from XYZ files** | Also load the tie lines of XYZ line data |
| **Import Data** | Imports the points as a layer in the project |

### Gridding

| Setting | Default | Meaning |
|---|---|---|
| **Select points layer** | | The point layer to grid |
| **Select Data Field** | | The field to grid |
| **Cell Size (proj Units)** | 100 | Output cell size, in the layer's units; **# Cells** shows the resulting grid size |
| **IDW Gridding** | | Opens QGIS's GRASS `v.surf.idw` (inverse distance weighting) dialog |
| **Threshold Error** | 0.0001 | B-spline: stop adding levels once no point misfits by more than this |
| **Max Levels** | 11 | B-spline: the number of levels of the multilevel B-spline |
| **Ignore values <** | unticked, -99999 | Leave out points with values below this (dummy values) |
| **Auto** | unticked | Choose the cell size and Max Levels from the spacing of the points |
| **BSpline Gridding** | | Grids the field with multilevel B-splines (translated from SAGA GIS) |

The B-spline grid is saved as `<points>_<field>_bspline.tif`. Method: [Gridding](filters-convolution-and-statistics.md#gridding), including how **Auto** chooses the settings.

### BSDWorms

| Setting | Default | Meaning |
|---|---|---|
| **Select Grid** | | The grid to worm |
| **# Levels** | 10 | Number of upward-continuation heights |
| **Base Level** | 1000 | The first height, in metres |
| **Increment** | 1000 | Step between heights, in metres |
| **Also save to shapefile** | unticked | Also save the worms as 3D polylines (PolylineZ, with the height as Z) |
| **Calculate Worms** | | Runs the calculation |

Finds the edges of sources at each height with wavelet transforms of the horizontal gradient. An edge that persists to height *h* comes from a source roughly *h*/2 deep. The points are saved as `<grid>_worms.csv` (and `<grid>_worms.shp`), with the padded grid, beside the input. Needs scikit-learn and networkx. Method: [BSDWorms](filters-convolution-and-statistics.md#wavelets).

---

## Petrophys tab

![The Petrophys tab](images/ui/tab_petrophys.png)

All the settings and the methods are described in [Petrophysics](petrophysics.md). In brief:

### Apparent Petrophysics

Tick one of the two calculations and press **Apply Processing** on this tab. Both need a grid in a metre-based projection and can be previewed.

| Setting | Default | Meaning |
|---|---|---|
| **Apparent Susceptibility** | | From an RTP magnetic anomaly grid (nT); output `_AppSus`, in SI |
| **Field (nT)**, **IGRF**, date | 50000, today | Total field strength; **IGRF** fills it in for the grid centre and date |
| **Depth to top (m)** | 0 | Depth of the top of the magnetic layer below the sensor |
| **Thickness (m)** | unticked (infinite), 1000 | Tick for a layer of finite thickness |
| **Background (SI)** | 0.001 | Added to the result, which is a variation about this value |
| **High-cut (m)** | blank (4 x cell size) | Removes wavelengths shorter than this, which the calculation would amplify; 0 turns it off |
| **Apparent Density** | | From a Bouguer gravity grid; output `_AppDens`, in g/cc |
| **Grid units** | mGal | mGal or µm/s² |
| **Depth to top (m)** | 0 | Depth of the top of the layer below the sensor |
| **Thickness (m)** | 5000 | Thickness of the layer |
| **Background (g/cc)** | 2.67 | Added to the result |
| **High-cut (m)** | blank (4 x cell size) | As for susceptibility |

**Example: apparent susceptibility of the RTP** (Field 36751 nT, depth 0, infinite thickness, background 0.001).

![Apparent susceptibility](images/examples/apparent_susceptibility.png)

### Apparent Henkel Plot / Histogram

| Control | Meaning |
|---|---|
| **App. Susceptibility**, **Log** | The apparent susceptibility grid, on a log axis when ticked |
| **App. Density** | The apparent density grid |
| **Plot Apparent** | Heat map of density against susceptibility for the whole area; with one grid blank, a histogram of the other |
| **Polygons**, **Category field** | A polygon layer (e.g. geology) and its unit field |
| **Plot per Category** | One panel per category |
| **Colour by position** | Colour each bin by its position on the plot instead of by its count |
| **Group similar**, number of groups | Group categories with similar distributions as A, B, C ... (Auto or 2 to 8) |

**Example: histogram of the apparent susceptibility** (App. Density left blank). Because there are negative values, the log axis is shifted (value - minimum + 0.0004), as its label says; the median, mean and SD in the corner are of the unshifted values.

![Histogram](images/examples/appsus_histogram.png)

### Measured Petrophysics

| Control | Meaning |
|---|---|
| **Points** | A point layer of samples |
| **X field**, **Y field**, multiplier, **Log** | The properties to plot, a factor to convert their units, and log axes |
| **Category field** | Colour and marker per category (up to 64) |
| **Plot Samples**, **Points** / **Heat map** | Draw every sample, or samples per bin |
| **Group similar** | As above |

### Save Plots

| Control | Default | Meaning |
|---|---|---|
| **PNG**, **PDF** | unticked | Also save each plot in these formats |
| **Folder** | beside the input layer | Where to save them |

---

## Utils tab

![The Utils tab](images/ui/tab_utils.png)

### SGTool Metadata

| Button | What it does |
|---|---|
| **Read metadata of selected grid** | Shows the history of the selected grid: each processing step, its settings and its source |
| **Save as XML** | Saves that history as `<grid file>.sgt.xml` beside the grid |
| **Apply the same steps to another grid...** | Repeats the steps of the selected grid's history on another grid |

See [Provenance metadata and replaying a processing history](provenance-and-replay.md).

### Threshold to NaN

| Setting | Default | Meaning |
|---|---|---|
| Tick box | unticked | Tick to apply with **Apply Processing** |
| **Above** / **Below** / **Between** | Above | Which values become no-data |
| **Above** | 1 | Upper threshold |
| **Below** | -1 | Lower threshold |

Output `_Clean`. Can be previewed.

**Example: above 1000** (of the TMI). The strongest anomalies become no-data (white).

![Threshold to NaN](images/examples/threshold_to_nan.png)

### Create Clipping Polygon

Tick **Create polygon based on grid** and press **Apply Processing** to outline the parts of the grid that have data, as a shapefile `<grid>_boundary`.

### Clip to Window

Clips the selected grid to the area shown in the map, on cell boundaries, and saves it as `<name>_clip.tif` (`_clip_01`, `_clip_02` ... if that exists).

### Clip All Layers

| Control | Default | Meaning |
|---|---|---|
| **Clip polygon** | Current map window | The map window, or a polygon layer (its selected polygons if any are selected) |
| **Clip all Layers** | | Asks for a folder, clips every layer of the project, and saves them with a copy of the project, `<project>_clip.qgz` |

See [Utilities](filters-convolution-and-statistics.md#utilities).

### Normalise Grids

| Setting | Meaning |
|---|---|
| **Input Directory** | Folder of GeoTIFFs to normalise |
| **Output Directory** | Folder for the results |
| **1st order fit** / **2nd order fit** | Order of the trend removed from each grid |
| **Normalise Grids** | Runs the calculation |

Gives every grid in the folder the same mean and standard deviation, after removing a 1st or 2nd order trend, to reduce steps where they are merged. Overlaps between the grids are not considered.

### Convert LUT to Grayscale

| Setting | Default | Meaning |
|---|---|---|
| **Select RGB Grid** | | A 3-band RGB image of a grid |
| **CSS Colour List or RGB triplets** | | The colour table, from the lowest to the highest value: CSS colour names (`teal, lemonchiffon, red`) or RGB triplets (`45,133,186, 251,248,183, 215,26,29`) |
| **Pick RGB** | | Pick colours from the screen |
| **Full list of CSS Colour Names** | | Opens the list of names |
| **Min**, **Max** | 0, 1000 | The values the first and last colours stand for |
| **Convert Grid** | | Saves the grey-scale grid as `<name>_gray.tif` |

---

## Help tab

Links to the example grid, Noddy models (the Atlas of Structural Geophysics and the Noddy installer), the code repository, the help file, and the people and projects whose code SGTool uses.

![The Help tab](images/ui/tab_help.png)

---

## Processing algorithm names

The filters are also QGIS Processing algorithms (`sgtool:<name>`), with the same settings as the dialog. See [Processing Algorithms](processing-algorithms.md).

| Algorithm | Name | Algorithm | Name |
|---|---|---|---|
| Remove line noise | `directional_line_noise` | Mean filter | `mean_filter` |
| Reduction to the pole | `reduction_to_pole` | Median filter | `median_filter` |
| Reduction to the pole (low latitude) | `reduction_to_pole_low_latitude` | Gaussian filter | `gaussian_filter` |
| Reduction to the equator | `reduction_to_equator` | Directional filter | `directional_filter` |
| Upward / downward continuation | `continuation` | Sun shading | `sun_shading` |
| Vertical integration | `vertical_integration` | Windowed statistic | `windowed_statistic` |
| Apparent susceptibility | `apparent_susceptibility` | Local anisotropy | `local_anisotropy` |
| Apparent density | `apparent_density` | Threshold to NaN | `threshold_to_nan` |
| Remove regional trend | `remove_regional` | Derivative | `derivative` |
| Band pass filter | `band_pass` | Tilt angle | `tilt_angle` |
| High / low pass filter | `high_low_pass` | Analytic signal | `analytic_signal` |
| Automatic gain control | `automatic_gain_control` | Total horizontal gradient | `total_horizontal_gradient` |

Euler deconvolution, PCA, ICA, B-spline gridding, worms and replaying a history are algorithms too. Differential RTP, MRVBF and the live preview are only in the dialog.
