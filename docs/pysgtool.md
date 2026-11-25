# pysgtool: SGTool as a Python package

pysgtool runs the SGTool calculations from Python, without QGIS or ArcGIS. It uses the same code as the plugin, so a filter gives the same result in both. Grids it saves carry the same processing history, and QGIS or ArcGIS Pro can read and replay it.

## Installation

From a copy of the repository:

```
pip install .            # calculations on arrays (numpy and scipy only)
pip install .[io]        # + reading and writing files, and the IGRF for a grid
pip install .[all]       # + PCA / ICA, worm lines, plots and WTMM
```

Use `pip install -e .[all]` to work on the code. GDAL is not needed: files are read and written with rasterio and pyogrio. Where the GDAL Python bindings (`osgeo`) are installed, as in QGIS and many conda environments, they are used instead. `sg.set_io_backend("rasterio")` or the environment variable `PYSGTOOL_IO=rasterio` overrides that.

## Quick start

```python
import pysgtool as sg

mag = sg.read("mag.tif")
rtp = mag.reduction_to_pole(date="2015-06-30")    # inclination and declination from the IGRF
rtp.derivative(direction="z", power=1).save()     # mag_RTP_d1z.tif, beside mag.tif
rtp.save("results/rtp.tif")
```

From an array (row 0 is the north, NaN is no data):

```python
grid = sg.Grid(array, dx=50, dy=50, origin=(xmin, ymax), crs="EPSG:28350", name="mag")
tilt = grid.tilt_angle()
tilt.values        # the result as a numpy array
```

Every method returns a new object and leaves its own unchanged.

A chain in memory keeps full (float64) precision, while the plugin saves each step as float32 and reads it back. One step gives the same result in both, but after several steps they can differ slightly (about 0.2% of the standard deviation in a test of three). To repeat the plugin's numbers exactly, save and re-read each step, or use `replay`.

Fourier filters can also differ slightly between numpy / scipy versions, because the regional plane they remove is fitted iteratively (in a test on a small grid, reduction to the pole by up to 0.06 nT where the standard deviation was 56 nT).

## Grids

| | |
|---|---|
| `sg.read(path, band=1, crs=None)` | A `Grid` from a GeoTIFF, ER Mapper or other GDAL raster, or a Geosoft `.grd` |
| `sg.read_stack(path)` | All the bands of a raster as a `Stack` |
| `grid.values`, `.dx`, `.dy`, `.extent`, `.gt`, `.crs`, `.epsg`, `.geographic`, `.x`, `.y` | The data and its position |
| `grid.save(path=None)` | Save as GeoTIFF with the history; without a path, beside the original with SGTool's usual name |
| `grid.history` | The processing steps so far |
| `grid.resample_to(other)` | The grid on the cells of another, reprojected if need be |

## Filters

The methods below have the names and settings of the QGIS Processing algorithms, which the [reference](reference.md) describes. `sg.filters()` lists them with their defaults, and `help(sg.Grid.band_pass)` shows one. Lengths are in map units.

| Group | Methods |
|---|---|
| Grav/Mag | `reduction_to_pole`, `reduction_to_equator`, `reduction_to_pole_low_latitude`, `continuation`, `vertical_integration` |
| Petrophysics | `apparent_susceptibility`, `apparent_density` |
| Frequency | `band_pass`, `high_low_pass`, `remove_regional`, `automatic_gain_control`, `directional_line_noise` |
| Gradient | `derivative`, `tilt_angle`, `analytic_signal`, `total_horizontal_gradient` |
| Convolution | `mean_filter`, `median_filter`, `gaussian_filter`, `directional_filter`, `sun_shading` |
| Statistics | `windowed_statistic`, `local_anisotropy` |
| Edges and domains | `nstd`, `nav_edge`, `kuwahara`, `snn` |
| Utilities | `threshold_to_nan` |

The reductions to the pole or equator need `inclination` and `declination`, and `apparent_susceptibility` needs `field`. Give them, or give `date=` to take them from the IGRF at the grid's centre. They have no defaults here.

## Other grid methods

| Method | Result |
|---|---|
| `remove_igrf(date)`, `igrf(date)`, `field_direction(date)` | Residual grid; inclination, declination and intensity grids; their values at the centre |
| `differential_rtp(date)` | Reduction to the pole with the field direction of every cell |
| `variable_continuation(observed, output)` | Continuation between uneven surfaces. Each is a height, a `Grid` of heights, or `(DEM Grid, clearance)` |
| `euler_deconvolution(structural_index, window, keep)` | `EulerSolutions`: `x`, `y`, `depth`, `elevation`, `base_level`, `std_dfdz`, `save()` |
| `euler_window_statistics(folder)` | The plugin's window statistics files |
| `worms(levels, base, increment)` | `Worms`: `x`, `y`, `z`, `value`, `save()`, `polylines().save()` |
| `mrvbf()` | MRVBF and MRRTF grids of a DEM |
| `power_spectrum()`, `plot_power_spectrum()` | Radially averaged power spectrum |
| `data_boundary()` | Polygons around the data, `save()` |

## Several bands or grids

```python
rad = sg.read_stack("k_th_u.tif")
pcs = rad.log_ratio("CLR").pca(3)          # also ica()
pcs.info["EXPLAINED_VARIANCE_RATIO"]

k, th, u = rad.band(1), rad.band(2), rad.band(3)
sg.efimov(k, u, th).save("efimov.tif")     # also kd, ud, ternary_gradient
sg.normalise([grid_a, grid_b])             # to the spread of the first
sg.Stack.from_grids([k, th, u])            # on the cells of the first
```

## Points

```python
pts = sg.read_points("gravity.shp")        # or .gpkg, or ASEG-GDF2: read_points("g.dat", x="EAST", y="NORTH")
grav = pts.grid("BOUGUER")                 # cell size and levels from the spacing of the points
sg.grid_points(x, y, z, cell_size=100, crs="EPSG:28350")
```

Gridding is by multilevel B-spline.

## Profiles, plots and history

```python
results = sg.wtmm(profile)                 # 1-D wavelet transform modulus maxima
sg.plot_wtmm(profile, results)

d, s = sg.petrophysics.pair(density, susceptibility)
fig = sg.petrophysics.heatmap(d.values, s.values, "Henkel plot", "Density (g/cc)", "Susceptibility (SI)")

sg.history("mag_RTP_d1z.tif")              # the steps that made it
sg.replay("mag_RTP_d1z.tif", "other.tif", date="2018-03-01")
```

`sg.igrf(lat, lon, date)` gives the IGRF at any position.

## Not in pysgtool yet

IDW and spline gridding, clipping, RGB to greyscale, the Noddy, CSV and XYZ imports, WTMM profiles from layers, and the MRVBF slope and RGB composite. The live preview, layer styling and project clipping belong to the GUIs.
