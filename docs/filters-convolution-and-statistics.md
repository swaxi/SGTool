[← Back to the README](../README.md)   

# Convolution filters, spatial statistics, gridding and utilities   

## Convolution Filters   
**Mean**
Applies a mean filter using a kernel of size n x n .   
   
**Median**   
Applies a median filter using a kernel of size n x n .   
   
**Gaussian**   
Applies a Gaussian filter with a specified standard deviation.    
   
**Directional**   
Apply directional filter (NE, N, NW, W, SW, S, SE, E)    
   
**Sun Shading**   
Computes relief shading for a digital elevation model (DEM) or other 2D grids.

## Spatial Statistics   
Calculates 1D statistics in a windowed grid
   
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

**Chain Length**   
Scores each pixel by the total size of the connected anisotropy component it belongs to. Two active pixels (those exceeding Aniso threshold) are linked if they lie within Search radius pixels of each other and their orientations differ by at most Angle tolerance degrees. The direction of the link is unconstrained so the full width of a multi-pixel-wide lineament is treated as a single entity. Connected components are found with Union-Find; every pixel in a component receives the same score equal to the number of pixels in that component.   
Returns: `_SS_ChainLen`   

**Streamline Length**   
From each active pixel (exceeding Aniso threshold) traces forward and backward along the orientation field using sub-pixel bilinear interpolation and double-angle orientation averaging. The path continues while anisotropy remains above the threshold and the step-to-step orientation change stays within Angle tolerance. The total forward + backward path length in pixels is the score. Unlike Chain Length, Streamline Length follows the curvature of a lineament and produces a continuous (non-integer) distance measure.   
Returns: `_SS_StreamLen`   

**MRRTF / MRVBF / Slope**   
Calculate DTM classification based on hill top & valley bottom curvature and slope  and combines as RGB.Converted from SAGA code  
   
## Multivariate Statistical Analysis   
**Principal Component Analysis**   
Principal Component Analysis (PCA) transforms correlated variables into orthogonal components that maximize variance, creating a new coordinate system where the first component captures the most variance.   
   
**Independent Component Analysis**   
Independent Component Analysis separates a multivariate signal into additive, statistically independent components by maximizing non-Gaussianity, often used to recover source signals from mixed observations.     
   
## Euler Deconvolution   
**Euler Deconvolution**   
Reliable Euler Deconvolution provides estimates of depth to gravity or magnetic sources based on analysis of gradients. Code derived from Reliable Euler Deconvolution by Felipe F. Melo and Valéria C.F. Barbosa https://github.com/ffigura/Euler-deconvolution-python.   
   
**Independent Component Analysis**   
Independent Component Analysis separates a multivariate signal into additive, statistically independent components by maximizing non-Gaussianity, often used to recover source signals from mixed observations.     
   
## Gridding   
**Import points**   
Imports point data in csv, ASEG-GDF2 dat or xyz formats. For xyz line data, tie lines can optionally be loaded as well.   

**Gridding**   
Grids point data using either BSpline (Converted from SAGA code, in both QGIS and ArcGIS Pro) or IDW built-in gridding algoithms (QGIS dialog; ArcGIS Pro toolbox)   
   
## Wavelets   
**BSDWorms**   
Use wavelet transforms to build multilevel "worms", saves out a single csv file of points (for use in 3D visualisation), and optionally a shapefile (for use in QGIS). Code from Frank Horowitz's bsdwormer  https://bitbucket.org/fghorow/bsdwormer/   
    
## Utilities   
**Threshold to NaN**   
Define upper or lower bound (or range) for which values will be set to NaN (i.e. excluded from display). Useful when reprojected images produce an unwanted border.      
   
**Create Clipping Polygon**   
Create one or more polygons outlining the available data in the grid. Useful, amongst other things, for clipping worms to grid area..     
   
**Normalise Grids**   
Normalise the means and standard deviations of a series of grids in a directory to minimse mismatches in merged grids. Does not consider overlaps between grids, simply standardises data and removes a first or second order regional.     
   
**Convert LUT to grayscale**   
Takes a 3-band registered RGB image and converts it to a monotonically increasing grayscale image if you provide the correct Look Up Table, either as matplotlib CSS Colour names (https://matplotlib.org/stable/gallery/color/named_colors.html#css-colors) or as a list of RGB triplets.   

   
