[← Back to the README](../README.md)   

# Convolution filters, spatial statistics, gridding and utilities   

Each function below gives its primary reference (*Reference*, listed in full at the [end of the page](#references)) and what SGTool's code does with it (*In SGTool*): whether the code is original, translated or adapted from someone else's, and where it differs from the published method.   

## Convolution Filters   
**Mean**
Applies a mean filter using a kernel of size n x n .   
*Reference*: a standard moving-average filter, with no single primary reference.   
*In SGTool*: original code. No-data cells are left out of each average (normalised convolution: the sum of the valid values divided by their number), and the edges are reflected.   
   
**Median**   
Applies a median filter using a kernel of size n x n .   
*Reference*: Tukey (1977).   
*In SGTool*: original code. No-data cells are left out of each median, and cells beyond the edges count as no data.   
   
**Gaussian**   
Applies a Gaussian filter with a specified standard deviation.    
*Reference*: a standard smoothing kernel, with no single primary reference.   
*In SGTool*: original code. The kernel is cut off at about 2 standard deviations (2⌈2σ⌉+1 cells across), normalised to sum to 1, and applied as a normalised convolution that leaves no-data cells out.   
   
**Directional**   
Apply directional filter (NE, N, NW, W, SW, S, SE, E)    
*Reference*: the 3 x 3 compass gradient kernels of Prewitt (1970).   
*In SGTool*: original code. For a kernel size n larger than 3, the 3 x 3 kernel is padded with zeros, so the result is the same as with n = 3.   
   
**Sun Shading**   
Computes relief shading for a digital elevation model (DEM) or other 2D grids.   
*Reference*: Horn (1981).   
*In SGTool*: with Relief ticked, a vectorised Python translation of the GRASS GIS `r.relief` module (GRASS Development Team), with 3 x 3 Sobel-type gradients; otherwise original code that shades the surface normal from central differences with a Lambertian light source. For a geographic grid the cell size is converted roughly to metres.   

## Spatial Statistics   
Calculates 1D statistics in a windowed grid   
*Reference* (Min to Skewness): standard sample statistics, with no single primary reference.   
*In SGTool* (Min to Skewness): original code. The window is a circle whose diameter is the window size, no-data cells are left out, and a window needs at least 3 valid values. Skewness and kurtosis come from `scipy.stats` (biased estimates; kurtosis is the excess kurtosis, 0 for a normal distribution).   
   
**Min**   
Calculate Minimum of values around central pixel for given window size  

**Max**   
Calculate Maximum of values around central pixel for given window size  

**Standard Deviation**   
Calculate Standard Deviation of values around central pixel for given window size  

**Variance**   
Calculate Variance of values around central pixel for given window size  

**Kurtosis**   
Calculate Kurtosis of values around central pixel for given window size  

**Skewness**   
Calculate Skewness of values around central pixel for given window size  

**Local Anisotropy**   
Uses the structure tensor (second-moment matrix) of local Sobel gradients to quantify the strength and orientation of linear structure at each pixel. The gradient products are smoothed over a window (Gaussian or box) of the specified size.   

Structure tensor elements:   
$`J_{11} = \text{smooth}(I_x^2), \quad J_{12} = \text{smooth}(I_x I_y), \quad J_{22} = \text{smooth}(I_y^2)`$   

Eigenvalues:   
$`\lambda_1 = \tfrac{1}{2}\!\left(J_{11}+J_{22}+\sqrt{(J_{11}-J_{22})^2+4J_{12}^2}\right), \quad \lambda_2 = \tfrac{1}{2}\!\left(J_{11}+J_{22}-\sqrt{(J_{11}-J_{22})^2+4J_{12}^2}\right)`$   

Anisotropy magnitude (saliency), normalised to [0, 1] by the 99th-percentile value:   
$`\text{AnisoMag} = \text{clip}\!\left(\tfrac{\lambda_1 - \lambda_2}{p_{99}},\; 0,\; 1\right)`$   

Dominant orientation:   
$`\theta = \tfrac{1}{2}\arctan\!\left(\tfrac{2J_{12}}{J_{11}-J_{22}}\right) \bmod 180°`$   

Returns two layers: `_SS_AnisoMag` (0–1 saliency; 0 = flat or isotropic, 1 = strong linear feature) and `_SS_AnisoOrient` (0–180° dominant strike).   
*Reference*: the structure tensor, Kass & Witkin (1987); Bigün, Granlund & Wiklund (1991).   
*In SGTool*: original code. The gradients are Sobel derivatives, and the normalisation of the saliency by its 99th percentile is SGTool's choice.   

**Chain Length**   
Scores each pixel by the total size of the connected anisotropy component it belongs to. Two active pixels (those exceeding Aniso threshold) are linked if they lie within Search radius pixels of each other and their orientations differ by at most Angle tolerance degrees. The direction of the link is unconstrained so the full width of a multi-pixel-wide lineament is treated as a single entity. Connected components are found with Union-Find; every pixel in a component receives the same score equal to the number of pixels in that component.   
Returns: `_SS_ChainLen`   
*Reference*: the Union-Find (disjoint-set) algorithm, Galler & Fisher (1964).   
*In SGTool*: original method and code.   

**Streamline Length**   
From each active pixel (exceeding Aniso threshold) traces forward and backward along the orientation field using sub-pixel bilinear interpolation and double-angle orientation averaging. The path continues while anisotropy remains above the threshold and the step-to-step orientation change stays within Angle tolerance. The total forward + backward path length in pixels is the score. Unlike Chain Length, Streamline Length follows the curvature of a lineament and produces a continuous (non-integer) distance measure.   
Returns: `_SS_StreamLen`   
*In SGTool*: original method and code.   

**MRRTF / MRVBF / Slope**   
Calculate DTM classification based on hill top & valley bottom curvature and slope  and combines as RGB.Converted from SAGA code  
*Reference*: Gallant & Dowling (2003); slope from Zevenbergen & Thorne (1987).   
*In SGTool*: a statement-by-statement Python translation of SAGA's MRVBF tool (O. Conrad, `ta_morphometry`), which agrees with SAGA to about 32-bit float precision. Combining MRRTF, MRVBF and slope into an RGB image is SGTool's addition.   
   
## Multivariate Statistical Analysis   
**Principal Component Analysis**   
Principal Component Analysis (PCA) transforms correlated variables into orthogonal components that maximize variance, creating a new coordinate system where the first component captures the most variance.   
*Reference*: Pearson (1901); Hotelling (1933).   
*In SGTool*: original code calling scikit-learn's `PCA`. The bands are standardised (zero mean, unit variance) first, and cells that are no-data in any band are left out.   
   
**Independent Component Analysis**   
Independent Component Analysis separates a multivariate signal into additive, statistically independent components by maximizing non-Gaussianity, often used to recover source signals from mixed observations.     
*Reference*: the FastICA algorithm, Hyvärinen (1999).   
*In SGTool*: original code calling scikit-learn's `FastICA`, with the same standardisation and no-data handling as PCA.   
   
## Euler Deconvolution   
**Euler Deconvolution**   
Reliable Euler Deconvolution provides estimates of depth to gravity or magnetic sources based on analysis of gradients. Code derived from Reliable Euler Deconvolution by Felipe F. Melo and Valéria C.F. Barbosa https://github.com/ffigura/Euler-deconvolution-python.   
*Reference*: Euler deconvolution, Thompson (1982) and Reid et al. (1990); the selection of reliable solutions, Melo & Barbosa (2020).   
*In SGTool*: the authors' `euler_python.py`, with the x and y axes swapped to match SGTool's grids (tested against the original code), and the standard deviation of the vertical derivative saved with each solution. The statistics of the solutions in grid windows (`window_stats`) are SGTool's addition to the authors' `estimates_statistics.py`.   
   
   
## Gridding   
**Import points**   
Imports point data in csv, ASEG-GDF2 dat or xyz formats. For xyz line data, tie lines can optionally be loaded as well.   
*In SGTool*: original code.   

**Gridding**   
Grids point data using either BSpline (Converted from SAGA code, in both QGIS and ArcGIS Pro) or IDW built-in gridding algoithms (QGIS dialog; ArcGIS Pro toolbox)   
*Reference*: multilevel B-splines, Lee, Wolberg & Shin (1997); inverse distance weighting, Shepard (1968).   
*In SGTool*: B-spline gridding is a statement-by-statement Python translation of SAGA's Multilevel B-Spline tool (O. Conrad, `grid_spline`), which agrees with SAGA to about 32-bit float precision. IDW is not SGTool's code: the QGIS dialog opens GRASS GIS `v.surf.idw`.   

**Choosing the B-spline settings (Auto)**   
Tick **Auto** next to the B-spline settings (ArcGIS Pro GUI: the **Auto** box beside the cell size; toolbox and Processing algorithm: the *choose automatically* option) to have the cell size and Max Levels worked out from the points. Typing either value yourself switches Auto off.

- **Cell size**: the median distance from a random place inside the survey to the nearest data point is measured. For flight lines it is a quarter of the line spacing, the usual rule of thumb for a gridding cell, so that distance is used; for scattered points it is about half the mean point spacing, so twice it is used. The result is rounded to two significant figures (a suggestion good to 30% does not need more). It works in the units of the layer, including degrees for latitude/longitude data, where east-west distances are shortened by cos(latitude) before measuring and the dialog also shows the cell's size in metres. A cell size you choose yourself is compared with the suggestion, and the dialog says so if it is more than twice it (detail is smoothed away) or less than 0.4 times it (most cells are far from any data).
- **Max Levels**: level L of the multilevel B-spline has a control lattice of range / 2^(L-1) (range being the longer side of the grid). Auto picks the level at which that lattice is half the output cell, `round(2 + log2(range / cell))`. Levels beyond that cannot add detail the grid can show, and on noisy data they start fitting the noise: in a held-out test on 824,000 aeromagnetic points the error was lowest at that level and slightly higher with more, while the run time kept growing with every level.
- **Threshold Error** is not chosen automatically. The calculation only stops early when no point is above it, which noisy data never allows, so Max Levels is what controls the run time.

*In SGTool* (Auto): original method and code.   

   
## Wavelets   
**BSDWorms**   
Worms are calculated by upward continuing the data to a series of heights (levels), finding the edges at each level with wavelet transforms of the horizontal gradient. The source depth of an edge is approximately **half the upward continuation height** of the level it appears at (Jacobsen, 1987), so edges that persist up to a level of 2000 m come from sources at about 1000 m depth. This is a guide, not a precise depth. Levels are set by the base level and the increment (metres) between levels.   

Use wavelet transforms to build multilevel "worms", saves out a single csv file of points (for use in 3D visualisation), and optionally a shapefile (for use in QGIS). Code from Frank Horowitz's bsdwormer  https://bitbucket.org/fghorow/bsdwormer/   
*Reference*: Hornby, Boschetti & Horowitz (1999); depth rule, Jacobsen (1987).   
*In SGTool*: the worm detection is Frank Horowitz's bsdwormer code. SGTool adds the preparation of the grid (mean removed, no-data filled, 20% linear-ramp padding and a Hanning window), the loop over levels, and the CSV and shapefile output.   

    
## Utilities   
*In SGTool*: all the utilities below are original code.   
**Threshold to NaN**   
Define upper or lower bound (or range) for which values will be set to NaN (i.e. excluded from display). Useful when reprojected images produce an unwanted border.      
   
**Create Clipping Polygon**   
Create one or more polygons outlining the available data in the grid. Useful, amongst other things, for clipping worms to grid area..     

**Clip to Window**   
Clips the grid selected in the Utils tab to the area currently shown in the map window. The **Clip to Window** button sits to the left of Clip All Layers. The clip is saved beside the original as `<name>_clip.tif` and loaded into the project with the same styling as the original, and its provenance records the source grid and the cell window that was cut out. All bands of a multiband grid are kept, and the cell size and values are unchanged (the window is cut on cell boundaries, so it can be up to a cell larger than the view). An earlier clip is never overwritten: if `<name>_clip.tif` already exists, or a layer of that name is loaded, the new clip is called `<name>_clip_01`, then `<name>_clip_02` and so on. The window is reprojected to the grid's coordinate system if the map uses a different one, and a window that does not overlap the grid gives a warning.   
**ArcGIS Pro**: the **Clip to Window** button on the GUI's Utils tab does the same for the input grid and the part of the ArcGIS Pro map you are looking at (it uses the map view the toolbox tool writes, so launch the GUI from the toolbox and keep the map open). The result is `<name>_clip.tif` (then `_clip_01`, ...), added to the layer list and loaded into the map with **Load SGTool Outputs**.   

**Clip All Layers**   
Pick a polygon layer in the dropdown (if some of its polygons are selected, only those are used) and press **Clip all Layers**, then choose a folder. Every layer in the project's layer tree is clipped to the polygon and saved there as `<name>_clip` (a GeoPackage for vector layers, a GeoTIFF for grids) with a `<name>_clip.qml` holding its style, and a copy of the project, `<project>_clip.qgz`, is saved beside them. It has the same groups, order, visibility and styles as the current project, loads the clipped files under the names `<name>_clip` with relative paths, and opens on the polygon. The clip polygon itself is saved too. The polygon is reprojected to each layer's coordinate system first. Layers that cannot be clipped (web layers such as WMS or XYZ tiles, tables without geometry) are carried over unchanged, a layer that fails is listed in the SGTool messages log and left out, and the temporary preview layer is skipped. Grids with no no-data value get -9999 outside the polygon. Based on a script by Mark Jessell.   
**ArcGIS Pro**: the toolbox tool **Clip All Layers to Polygon (Current Project)** (Utils) does the same for a map of the open project: choose the polygon layer, an output folder and (optionally) the map. The clipped vector layers go into one file geodatabase, `<project>_clip.gdb`, as `<name>_clip` feature classes (the polygon too), the grids are `<name>_clip.tif`, each layer's symbology is saved as `<name>_clip.lyrx`, and `<project>_clip.aprx`, a copy of the project, has every clipped layer pointing at its new data and renamed `<name>_clip`, so groups, order, symbology, bookmarks and layouts are kept. Only the chosen map is changed. Web layers, annotation and tables are left as they are, a layer's definition query is applied (its selection is not), and the tool stops when you press Cancel. It must be run from inside ArcGIS Pro with the project open.   
   
**Normalise Grids**   
Normalise the means and standard deviations of a series of grids in a directory to minimse mismatches in merged grids. Does not consider overlaps between grids, simply standardises data and removes a first or second order regional.     
   
**Convert LUT to grayscale**   
Takes a 3-band registered RGB image and converts it to a monotonically increasing grayscale image if you provide the correct Look Up Table, either as matplotlib CSS Colour names (https://matplotlib.org/stable/gallery/color/named_colors.html#css-colors) or as a list of RGB triplets.

## References   

- Bigün, J., Granlund, G.H., Wiklund, J., 1991. Multidimensional orientation estimation with applications to texture analysis and optical flow. *IEEE Transactions on Pattern Analysis and Machine Intelligence* 13, 775–790. https://doi.org/10.1109/34.85668   
- Gallant, J.C., Dowling, T.I., 2003. A multiresolution index of valley bottom flatness for mapping depositional areas. *Water Resources Research* 39(12), 1347. https://doi.org/10.1029/2002WR001426   
- Galler, B.A., Fisher, M.J., 1964. An improved equivalence algorithm. *Communications of the ACM* 7, 301–303. https://doi.org/10.1145/364099.364331   
- Hornby, P., Boschetti, F., Horowitz, F.G., 1999. Analysis of potential field data in the wavelet domain. *Geophysical Journal International* 137, 175–196. https://doi.org/10.1046/j.1365-246x.1999.00788.x   
- Horn, B.K.P., 1981. Hill shading and the reflectance map. *Proceedings of the IEEE* 69, 14–47. https://doi.org/10.1109/PROC.1981.11918   
- Hotelling, H., 1933. Analysis of a complex of statistical variables into principal components. *Journal of Educational Psychology* 24, 417–441. https://doi.org/10.1037/h0071325   
- Hyvärinen, A., 1999. Fast and robust fixed-point algorithms for independent component analysis. *IEEE Transactions on Neural Networks* 10, 626–634. https://doi.org/10.1109/72.761722   
- Jacobsen, B.H., 1987. A case for upward continuation as a standard separation filter for potential-field maps. *Geophysics* 52, 1138–1148. https://doi.org/10.1190/1.1442378   
- Kass, M., Witkin, A., 1987. Analyzing oriented patterns. *Computer Vision, Graphics, and Image Processing* 37, 362–385. https://doi.org/10.1016/0734-189X(87)90043-0   
- Lee, S., Wolberg, G., Shin, S.Y., 1997. Scattered data interpolation with multilevel B-splines. *IEEE Transactions on Visualization and Computer Graphics* 3, 228–244. https://doi.org/10.1109/2945.620490   
- Melo, F.F., Barbosa, V.C.F., 2020. Reliable Euler deconvolution estimates throughout the vertical derivatives of the total-field anomaly. *Computers & Geosciences* 138, 104436. https://doi.org/10.1016/j.cageo.2020.104436   
- Pearson, K., 1901. On lines and planes of closest fit to systems of points in space. *Philosophical Magazine* 2, 559–572. https://doi.org/10.1080/14786440109462720   
- Prewitt, J.M.S., 1970. Object enhancement and extraction. In: Lipkin, B.S., Rosenfeld, A. (Eds.), *Picture Processing and Psychopictorics*. Academic Press, New York, 75–149.   
- Reid, A.B., Allsop, J.M., Granser, H., Millett, A.J., Somerton, I.W., 1990. Magnetic interpretation in three dimensions using Euler deconvolution. *Geophysics* 55, 80–91. https://doi.org/10.1190/1.1442774   
- Shepard, D., 1968. A two-dimensional interpolation function for irregularly-spaced data. *Proceedings of the 1968 23rd ACM National Conference*, 517–524. https://doi.org/10.1145/800186.810616   
- Thompson, D.T., 1982. EULDPH: a new technique for making computer-assisted depth estimates from magnetic data. *Geophysics* 47, 31–37. https://doi.org/10.1190/1.1441278   
- Tukey, J.W., 1977. *Exploratory Data Analysis*. Addison-Wesley, Reading, MA.   
- Zevenbergen, L.W., Thorne, C.R., 1987. Quantitative analysis of land surface topography. *Earth Surface Processes and Landforms* 12, 47–56. https://doi.org/10.1002/esp.3290120107   
