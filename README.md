# SGTool v 0.2
 Simple Potential Field Calcs to assist WAXI/Agate Structural Geophysics Course    
 https://waxi4.org   
 https://agate-project.org    
    
![SGTools image](dialog.png)       
   
**Directional Band Pass**   
Butterworth High-Pass Filter
$`H(k) = \frac{1}{1 + \left(\frac{k_c}{k}\right)^{2n}}`$    
The Butterworth filter attenuates frequencies below the cutoff k<sub>c</sub> while preserving higher frequencies.    
H(k) : Filter response as a function of wavenumber k.    
k : Wavenumber (spatial frequency).    
k<sub>c</sub> : Cutoff wavenumber, related to the cutoff wavelength by k<sub>c</sub> = \frac{1}{\text{cutoff wavelength}}.    
n : Filter order, determining the sharpness of the transition. Higher \( n \) makes the filter more selective.   
    
Directional Cosine Filter   
$`H(k_x, k_y) = \left| \cos(\theta - \theta_c) \right|^p`$   
The Directional Cosine Filter emphasizes or suppresses frequency components along a specific direction.   
H(k<sub>x</sub>, k<sub>y</sub>): Filter response as a function of wavenumber components k<sub>x</sub> and k<sub>y</sub>.   
theta = \arctan\left(\frac{k_y}{k_x}\right) : Angle of the frequency component.   
theta<sub>c</sub> : Center direction (in radians), representing the direction to emphasize.   
p : Degree of the cosine function. Higher \( p \) sharpens the directional emphasis.   

**Band Pass**
$`H(k) = e^{-(k - k_c)^2 / (2 \sigma^2)} - e^{-(k + k_c)^2 / (2 \sigma^2)}`$   
The band-pass filter retains frequencies within a specified range, suppressing both low and high frequencies outside this range.
Where   
k<sub>c</sub> : The central frequency of the band.
sigma : The width of the frequency band.   

**Reduction to the Pole**    
$`H(k_x, k_y) = \frac{1}{\Theta_m \Theta_f}`$   
Converts magnetic data measured at any inclination and declination to what it would be if measured at the magnetic pole.
Where   
Theta<sub>m</sub> : Magnetization vector of the source.
Theta<sub>f/sub> : Geomagnetic field vector.   

**Reduction to the Equator**    
$`H(k_x, k_y) = \frac{1}{f_z + i \frac{f_e k_x + f_n k_y}{|\mathbf{k}|}}`$   
Converts magnetic data measured at any inclination and declination to what it would be if measured at the magnetic equator.
Where   
f<sub>z</sub> : Downward component of the geomagnetic field.
f<sub>e</sub>, f<sub>n</sub> : Easting and northing components of the geomagnetic field.   

**Derivative**    
$`\frac{\partial f}{\partial u} = \frac{\partial f}{\partial x} \cos\theta + \frac{\partial f}{\partial y} \sin\theta`$   
Where   
theta is the angle defining the direction of the derivative (x,y or z).   

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

**High Pass**    
$`H(k) = 1 - e^{-k^2 / (2 k_c^2)}`$   
The high-pass filter removes low-frequency components (long wavelengths) while retaining high-frequency components (short wavelengths).
Where    
k<sub>c</sub> : The cutoff frequency where the filter begins attenuating lower frequencies.   

**Low Pass**    
$`H(k) = e^{-k^2 / (2 k_c^2)}`$   
The low-pass filter suppresses high-frequency components (short wavelengths) while preserving low-frequency components (long wavelengths).
Where: k<sub>c</sub> : The cutoff frequency where the filter begins attenuating higher frequencies.   

**Remove Regional**   
$`H(k) = e^{-k^2 / (2 k_c^2)}`$
The low-pass filter suppresses high-frequency components (short wavelengths) while preserving low-frequency components (long wavelengths).
Where    
k<sub>c</sub> : The cutoff frequency where the filter begins attenuating higher frequencies.   

**Continuation**    
$`H(k) = e^{-k h}`$   
Where   
h > 0 for upward continuation.   
h < 0  for downward continuation.   

**Automatic Gain Control**    
$`AGC(x, y) = \frac{f(x, y)}{\text{RMS}(f(x, y), w)}`$   
Where    
RMS(f, w)  is the root mean square of the data over a window w.   

**Vertical Integration**   
$`H(k_x, k_y) = \frac{1}{k}`$  
When applied to an RTE or RTP image provides the so called Pseudogravity result    
Where    
k = sqrt{k<sub>x</sub><sup>2</sup> + k<sub>y</sub><sup>2</sup>} .   

**Total Horizontal Gradient**   
$`THG(x, y) = \sqrt{\left(\frac{\partial f}{\partial x}\right)^2 + \left(\frac{\partial f}{\partial y}\right)^2}`$   
   
**Radially averaged power spectrum (but needs testing!)**    
$`P(k) = \frac{1}{N_k} \sum_{(k_x, k_y) \in k} |\text{FFT}(f)|^2`$   
Where    
P(k) is the radially averaged power spectrum, and N<sub>k</sub> is the number of samples in the radial bin.   

# Installation
- Dowload the zip file from the green **<> Code** button and install the zip file in QGIS using the plugin manager   
   
# Inputs   
- Supports data geotiff, grd, ers formats

# Tips
- Simple Potential field calculations, mostly FFT-based
- Adds suffix (e.g. _RTP) to input filename and stores results in same directory. Converts grd to tif on loading, but leaves ers as is.
- Calculates IGRF mag field parameters based on centroid of grid, plus date and survey height
- If a layer with a given name already loaded, no calculation is performed
- Processing methods preceded by a dot should be performed on line-direction-filtered and reduced to pole/equator data   
- Length units are defined by grid properties (so Lat/Long wavelengths should be defined in degrees!)

# Code development
- Plugin construction - Mark Jessell using QGIS Plugin Builder Plugin    
- Calcs ChatGPT and Mark Jessell
- IGRF calculation - https://github.com/klaundal/ppigrf  
- GRD Loader Mark Jessell & Fatiando a Terra crew https://www.fatiando.org/
- Example geophysics data in image above courtesy of CPRM http://www.cprm.gov.br/



