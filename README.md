# Structural Geophysics Tool v0.3.00
 Simple Potential Field and other Geophysical Grid Calcs to assist WAXI/Agate Structural Geophysics Course    
 https://waxi4.org   and  https://agate-project.org    
    
 This plugin is available directly within QGIS using the Plugin Manager, however the latest ArcGIS Pro version or QGIS version with all new bugs will always be at this site.   

**<a href="https://tectonique.net/sgtools_data/SGTool%20Cheat%20Sheet.pdf">Cheat Sheet (Thanks Chops)</a>**&nbsp;**&nbsp;|&nbsp;&nbsp;&nbsp; <a href="https://tectonique.net/sgtools_data/Structural%20Geophysics%20Tools.pdf">Download Basic Help Document</a>**&nbsp;&nbsp;&nbsp; |&nbsp;&nbsp;&nbsp;  **<a href="https://tectonique.net/sgtools_data/SGTools_large.mp4">Ctrl-click on link to watch demo video</a>**
    
![SGTools image](dialog.png)       

# changelog=0.3.00   
      * Add Spatial stats and convolutions to ArcGIS Pro toolbox
      * Move ArcGIS Pro files to their own directory    
      * Update IGRF to allow up to 2030    
      * Move Euler and IGRF to calcs directory   
      * Compatibility with both QGIS4/QT6 and QGIS3/QT5   
      * Add RGB triplets as LUT definition for RGB to greyscale convertor   
      0.2.18    
      * Fix (again) clipping of csv worms   
      * Fix bug in Filling FFT NaNs in NS direction (Thanks Feargal!)   
      0.2.17    
      * Fix bug saving grid boundary filename   
      * Add ArcPro Toolbox support (very Beta!!)   
      * Remove all QGIS library calls from calc engines   
      * Fix bug in worms spatial filtering   
      * Force processed ers output files to have tif suffix  
      * Check for recognised CRS before processing    
      0.2.16
      * Add Euler Deconvolution from https://github.com/ffigura/Euler-deconvolution-python
      * Add Noddy grid (grv & mag) import
      0.2.15 
      * Use cosine rolloff for high/low and bandpass filters to reduce ringing
      * Change Remove Regional to 1st or 2nd order polynomial
      * New grids use stdev scaling of grid to better highlight features
      0.2.14
      * don't try to load LINE_ID codes if they don't exist in WTMM
      * added PCA & ICA calculations for multiband grids
      * speed up Grass-like relief calc
      * added DAT format points import based on ASEG-GDF2 standard (tested against example datasets only)
      * TA demoted to needing RTP_E first in GUI
      * save all CSV/XYZ/DAT imports as shapefiles
      * fix import matplotlib bug and test for LINE_ID for points in WTMM code
   
Full changelog <a href="https://raw.githubusercontent.com/swaxi/SGTool/refs/heads/main/metadata.txt">Metadata</a>   


# Installation
## QGIS:
1) Either:   
- Download the zip file from the green **<> Code** button and install the zip file in QGIS using the plugin manager for the version on github or   
- Install directly from the QGIS plugin manager from the plugin repository   
2) If you get an warning of the type **The following Python packages are required but not installed: scikit-learn** or any other module name you can install it directly from the QGIS Python Console (Menu Plugins->Python Console) and then type in (for the scikit-learn example):   
   
   **!pip3 install scikit-learn**   
      
   The packages required for specific functions are:   
   **matplotlib** WTMM, Radial Power Spectrum   
   **scikit-learn** BSDWorms, PCA, ICA   
   **PyWavelets** WTMM    

   If you don't use these functions, there is no need to install the extra packages.   
3) For BSpline Gridding you need to install the plugin **Processing Saga NextGen Provider**   
   
## ArcGIS Pro:
1) Download and unzip this respository and store somewhere safe.
2) In the ArcGIS Pro Catalogue area, go to Add toolbox and select the file **GeophysicalProcessor.pyt** in the **ArcGIS_Pro** directory in this repository. Double click on the new Geophysical Processing Toolbox to get the list of functions that can be run (still in Beta so limited to classical geophysical, convolution and spataial stats processing calls for now):   
- # --- Grav/Mag ---
- ReductionToPole,
- UpwardContinuation,
- DownwardContinuation,
- VerticalIntegration,
- # --- Frequency ---
- BandPassFilter,
- HighPassFilter,
- LowPassFilter,
- DirectionalButterworthBandPass,
- RemoveRegionalTrend,            
- AutomaticGainControl,
- # --- Gradient ---
- AnalyticSignal,
- TiltAngle,
- ComputeDerivative,
- TotalHorizontalGradient,
- # --- Spatial Stats ---
- SpatialStatsVariance,
- SpatialStatsStdDev,
- SpatialStatsSkewness,
- SpatialStatsKurtosis,
- SpatialStatsMin,
- SpatialStatsMax,
- # --- Convolution ---
- ConvolutionMeanFilter,
- ConvolutionMedianFilter,
- ConvolutionGaussianFilter,
- ConvolutionDirectionalFilter,
- ConvolutionSunShading,
   
# Inputs   
- QGIS version supports data geotiff, grd, ers and Noddy (grv & mag) grid formats plus any grid format already supported by QGIS. ArcGIS Pro version supports any raster format supported by ArcGIS Pro.
- Supports csv, dat, xyz plus any point format already supported by QGIS
- Existing Noddy mag and grav files can be found at the Atlas of Structural Geophysics: https://tectonique.net/asg/
- New Noddy models can be calculated using the Windows version at https://tectonique.net/noddy/OpenNoddy_installer.exe or a python wrapper at https://github.com/cgre-aachen/pynoddy
   
# Capabilities   

## Grav/Mag Filters   
   
**Reduction to the Pole**    
$`H_{RTP}(k_x, k_y) = \frac{k \cos I \cos D + i k_y \cos I \sin D + k_x \sin I}{k}`$   
Converts magnetic data measured at any inclination and declination to what it would be if measured at the magnetic pole.
Where     
- k<sub>x</sub> and k<sub>y</sub> : The wavenumber components in the x and y directions.
- k = The total wavenumber magnitude = sqrt{k<sub>x</sub><sup>2</sup> + k<sub>y</sub><sup>2</sup>}   
- I : Magnetic inclination (in radians).
- D : Magnetic declination (in radians).
- i : Imaginary unit.


**Reduction to the Equator**    
$`H_{RTE}(k_x, k_y) = \frac{k \cos I \cos D + i k_y \cos I \sin D + k_x \sin I}{k \cos I \cos D - i k_y \cos I \sin D + k_x \sin I}`$     
Converts magnetic data measured at any inclination and declination to what it would be if measured at the magnetic equator.
Where   
- k<sub>x</sub> and k<sub>y</sub> : The wavenumber components in the x and y directions.
- k = The total wavenumber magnitude = sqrt{k<sub>x</sub><sup>2</sup> + k<sub>y</sub><sup>2</sup>}
- I : Magnetic inclination (in radians).
- D : Magnetic declination (in radians).
- i : Imaginary unit. 
   
**Continuation**    
$`H(k) = e^{-k h}`$   
Where   
h > 0 for upward continuation.   
h < 0  for downward continuation.   
   
**Vertical Integration**   
$`H(k_x, k_y) = \frac{1}{k}`$  
When applied to an RTE or RTP image provides the so called Pseudogravity result    
Where    
k = sqrt{k<sub>x</sub><sup>2</sup> + k<sub>y</sub><sup>2</sup>} .   
   
## Frequency Filters   
   
**High Pass Filter**

$$H(k) = \begin{cases} 
0 & \text{if } k \leq k_{low} \\
\frac{1}{2}\left(1 - \cos\left(\pi \frac{k - k_{low}}{k_{high} - k_{low}}\right)\right) & \text{if } k_{low} < k < k_{high} \\
1 & \text{if } k \geq k_{high}
\end{cases}$$

The high-pass filter removes low-frequency components (long wavelengths) while retaining high-frequency components (short wavelengths) with a smooth transition to reduce ringing artifacts.

Where:
- $k$ : Current wavenumber magnitude $\sqrt{k_x^2 + k_y^2}$
- $k_{low} = \frac{2\pi}{\lambda_c + w/2}$ : Lower transition boundary
- $k_{high} = \frac{2\pi}{\lambda_c - w/2}$ : Upper transition boundary  
- $\lambda_c$ : Cutoff wavelength
- $w$ : Transition width (in same units as wavelength)

**Low Pass Filter**

$$H(k) = \begin{cases} 
1 & \text{if } k \leq k_{inner} \\
\frac{1}{2}\left(1 + \cos\left(\pi \frac{k - k_{inner}}{k_{outer} - k_{inner}}\right)\right) & \text{if } k_{inner} < k < k_{outer} \\
0 & \text{if } k \geq k_{outer}
\end{cases}$$

The low-pass filter removes high-frequency components (short wavelengths) while retaining low-frequency components (long wavelengths). Optional smooth transition reduces potential ringing artifacts.

Where:
- $k$ : Current wavenumber magnitude $\sqrt{k_x^2 + k_y^2}$
- $k_{inner} = \frac{2\pi}{\lambda_c + w/2}$ : Inner transition boundary
- $k_{outer} = \frac{2\pi}{\lambda_c - w/2}$ : Outer transition boundary
- $\lambda_c$ : Cutoff wavelength
- $w$ : Transition width (optional, in same units as wavelength)

**Band Pass Filter**

$$H(k) = H_{high}(k) \times H_{low}(k)$$

Where:

***High-pass component:***
$$H_{high}(k) = \begin{cases} 
0 & \text{if } k \leq k_{h,low} \\
\frac{1}{2}\left(1 - \cos\left(\pi \frac{k - k_{h,low}}{k_{h,high} - k_{h,low}}\right)\right) & \text{if } k_{h,low} < k < k_{h,high} \\
1 & \text{if } k \geq k_{h,high}
\end{cases}$$

***Low-pass component:***
$$H_{low}(k) = \begin{cases} 
1 & \text{if } k \leq k_{l,inner} \\
\frac{1}{2}\left(1 + \cos\left(\pi \frac{k - k_{l,inner}}{k_{l,outer} - k_{l,inner}}\right)\right) & \text{if } k_{l,inner} < k < k_{l,outer} \\
0 & \text{if } k \geq k_{l,outer}
\end{cases}$$

The band-pass filter isolates features within a specific wavelength range by combining high-pass and low-pass components with smooth transitions.

Where:
- $k$ : Current wavenumber magnitude $\sqrt{k_x^2 + k_y^2}$
- $k_{h,low} = \frac{2\pi}{\lambda_{low} + w_h/2}$ : High-pass lower transition boundary
- $k_{h,high} = \frac{2\pi}{\lambda_{low} - w_h/2}$ : High-pass upper transition boundary
- $k_{l,inner} = \frac{2\pi}{\lambda_{high} + w_l/2}$ : Low-pass inner transition boundary
- $k_{l,outer} = \frac{2\pi}{\lambda_{high} - w_l/2}$ : Low-pass outer transition boundary
- $\lambda_{low}$ : Low cutoff wavelength (removes longer wavelengths)
- $\lambda_{high}$ : High cutoff wavelength (removes shorter wavelengths)
- $w_h$ : High-pass transition width
- $w_l$ : Low-pass transition width

**Directional Band Pass**   
Removes combined directional and high pass filtered data from original data, with scaling function modify extent of feature suppression.   
   
***Butterworth High-Pass Filter***
$`H(k) = \frac{1}{1 + \left(\frac{k_c}{k}\right)^{2n}}`$    
The Butterworth filter attenuates frequencies below the cutoff k<sub>c</sub> while preserving higher frequencies.    
H(k) : Filter response as a function of wavenumber k.    
k : Wavenumber (spatial frequency).    
k<sub>c</sub> : Cutoff wavenumber, related to the cutoff wavelength by k<sub>c</sub> = \frac{1}{\text{cutoff wavelength}}.    
n : Filter order, determining the sharpness of the transition. Higher \( n \) makes the filter more selective.   
    
***Directional Cosine Filter***    
$`H(k_x, k_y) = \left| \cos(\theta - \theta_c) \right|^p`$   
The Directional Cosine Filter emphasizes or suppresses frequency components along a specific direction.   
H(k<sub>x</sub>, k<sub>y</sub>): Filter response as a function of wavenumber components k<sub>x</sub> and k<sub>y</sub>.   
theta = \arctan\left(\frac{k_y}{k_x}\right) : Angle of the frequency component.   
theta<sub>c</sub> : Center direction (in radians), representing the direction to emphasize.   
p : Degree of the cosine function. Higher \( p \) sharpens the directional emphasis.   

**Remove Regional**   
Remove a 1st order (dipping plane) or 2nd order (parabolic plane) regional from data. 
   
**Automatic Gain Control**    
$`AGC(x, y) = \frac{f(x, y)}{\text{RMS}(f(x, y), w)}`$   
Where    
RMS(f, w)  is the root mean square of the data over a window w.   
   
**Radially averaged power spectrum (but needs testing!)**    
$`P(k) = \frac{1}{N_k} \sum_{(k_x, k_y) \in k} |\text{FFT}(f)|^2`$   
Where    
P(k) is the radially averaged power spectrum, and N<sub>k</sub> is the number of samples in the radial bin.   
   
## Gradient Filters   
   
**Derivative**    
$`\frac{\partial f}{\partial u} = \frac{\partial f}{\partial x} \cos\theta + \frac{\partial f}{\partial y} \sin\theta`$   
Where   
theta is the angle defining the direction of the derivative (x,y or z).   
   
**Total Horizontal Gradient**   
$`THG(x, y) = \sqrt{\left(\frac{\partial f}{\partial x}\right)^2 + \left(\frac{\partial f}{\partial y}\right)^2}`$   
   
**Analytic Signal**    
$`A(x, y) = \sqrt{\left(\frac{\partial f}{\partial x}\right)^2 + \left(\frac{\partial f}{\partial y}\right)^2 + \left(\frac{\partial f}{\partial z}\right)^2}`$   
Computes the total amplitude of the gradients, independent of field inclination or declination.
Useful for locating edges of potential field sources (e.g., faults or contacts).   
       
**Tilt Angle**    
$`T = \tan^{-1}\left(\frac{\frac{\partial f}{\partial z}}{\sqrt{\left(\frac{\partial f}{\partial x}\right)^2 + \left(\frac{\partial f}{\partial y}\right)^2}}\right)`$   
Enhances the contrast of geological features by highlighting gradients relative to the vertical component.
Where   
df/dz : Vertical derivative of the field.
df/dx , df/dy : Horizontal derivatives of the field.   
   
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

**DTM Curvature Classifier**   
Calculate DTM classification based on curvature and slope   
Based on Curvature Threshold, Cliff Threshold, Window Size and Smoothing Parameter   
Classified array where: -1 = concave up; 0 = flat; 1 = convex up and 2 = steep slope  
   
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
Imports point data in csv, ASEG-GDF2 dat or xyz formats   

**Gridding**   
Grids point data using either BSpline or IDW built-in gridding algoithms   
   
## Wavelets   
**BSDWorms**   
Use wavelet transforms to build multilevel "worms", saves out a single csv file of points (for use in 3D visualisation), and optionally a shapefile (for use in QGIS). Code from Frank Horowitz's bsdwormer  https://bitbucket.org/fghorow/bsdwormer/   
    
**WTMM**   
Use wavelet transforms to build multilevel analysis (Wavelet Transform Modulus Maxima) along a selected linestring (polyline) profile extracted from grid or for imported XYZ data. 
    
## Utilities   
**Threshold to NaN**   
Define upper or lower bound (or range) for which values will be set to NaN (i.e. excluded from display). Useful when reprojected images produce an unwanted border.      
   
**Create Clipping Polygon**   
Create one or more polygons outlining the available data in the grid. Useful, amongst other things, for clipping worms to grid area..     
   
**Normalise Grids**   
Normalise the means and standard deviations of a series of grids in a directory to minimse mismatches in merged grids. Does not consider overlaps between grids, simply standardises data and removes a first or second order regional.     
   
**Convert LUT to grayscale**   
Takes a 3-band registered RGB image and converts it to a monotonically increasing grayscale image if you provide the correct Look Up Table, uses matplotlib CSS Colour names: https://matplotlib.org/stable/gallery/color/named_colors.html#css-colors   

   
# How To   
1) Load a raster image from file
- If a GRD grid (Oasis Montaj) is selected, the plugin will attempt to load CRS from the associated xml file, if this is not possible a CRS of EPSG:4326 is assumed. In any case the grid is saved as geotiff.
2) Whatever layer is shown in the layer selector will be the one processed by whatever combination of filters are selected by check boxes, **but must exist as a file, this plugin cannot process grids that are only in memory**. 
- All processed files will be saved as geotiffs, and will be saved in the same directory as the original file, and will have a suffix added describing the processing step.
- If a RTP or RTE calculation is performed, it is possible to define the magnetic field manually or the IGRF mag field parameters can be assigned based on the centroid of grid, plus survey date, or embedded geotiff metadata if the source of the tif was a Noddy grid file..
- If a file exists on disk it will be overwritten, although QGIS plugins don't always like saving to disks other than C: on Windows, and can't overwrite a file if the grid is open in another program.
- Length units are defined by grid properties except for Up/Down Continuation (so Lat/Long wavelengths should be defined in degrees!)
3) If multiple processing steps are required, first apply one process, select the result and then apply subsequent steps.

# Alternatives   
There are several excellent Open Source or at least free alternatives to this plugin if you don't want to use QGIS, or want to do things this plugin can't:   
- Fatiando Harmonica https://www.fatiando.org/harmonica/
- GravMagSuite https://github.com/fcastro25/GravMagSuite
- GSSH https://cires1.colorado.edu/people/jones.craig/GSSH/index.html
- UBC Toolkit https://toolkit.geosci.xyz/content/Demos/SyntheticFilters.html
- gravmag https://github.com/birocoles/gravmag
- Fourpot https://sites.google.com/view/markkussoftware/gravity-and-magnetic-software/fourpot
- GridMerge https://www.gridmerge.com.au/home
- GammaSpec https://www.gammaspec.com.au/ 


# Code development
- You can explore the codebase and functionality at the <a href="https://deepwiki.com/swaxi/SGTool">DeepWiki Code description</a>   
- Calcs ChatGPT, Claude and Mark Jessell
- Plugin construction - Mark Jessell using QGIS Plugin Builder Plugin https://g-sherman.github.io/Qgis-Plugin-Builder/    
- IGRF calculation -  using pyIGRF https://github.com/ciaranbe/pyIGRF
- GRD Loader & Radially averaged power spectrum Fatiando a Terra crew & Mark Jessell https://www.fatiando.org/
- Example geophysics data in image above courtesy of Mauritania Govt. and USGS https://anarpam.mr/en/     
- Worming of grids uses Frank Horowitz's bsdwormer  https://bitbucket.org/fghorow/bsdwormer/
- Wavelet Transform base code - https://github.com/PyWavelets/pywt 
- Multilevel BSpline Gridding piggybacks off SAGA code via the plugin Processing Saga NextGen Provider https://github.com/north-road/qgis-processing-saga-nextgen   
- Euler Deconvolution uses Felipe F. Melo and Valéria C.F. Barbosa's Reliable Euler method https://github.com/ffigura/Euler-deconvolution-python   




