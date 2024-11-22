# SGTool v 0.2
 Simple Potential Field Calcs to assist WAXI/Agate Structural Geophysics Course    
 https://waxi4.org   
 https://agate-project.org    
    
- Directional Band Pass 
$`H(k_x, k_y) = \exp\left(-\frac{(k_x \cos\theta + k_y \sin\theta - k_c)^2}{2 \sigma^2}\right)`$
- Band Pass   
$`H(k) = e^{-(k - k_c)^2 / (2 \sigma^2)} - e^{-(k + k_c)^2 / (2 \sigma^2)}`$
- RTP
$`H(k_x, k_y) = \frac{1}{\Theta_m \Theta_f}`$
- RTE
$`H(k_x, k_y) = \frac{1}{f_z + i \frac{f_e k_x + f_n k_y}{|\mathbf{k}|}}`$
- Derivative
- Analytic Signal   
$`A(x, y) = \sqrt{\left(\frac{\partial f}{\partial x}\right)^2 + \left(\frac{\partial f}{\partial y}\right)^2 + \left(\frac{\partial f}{\partial z}\right)^2}
`$
- Tilt Angle
$`T = \tan^{-1}\left(\frac{\frac{\partial f}{\partial z}}{\sqrt{\left(\frac{\partial f}{\partial x}\right)^2 + \left(\frac{\partial f}{\partial y}\right)^2}}\right)`$
- High Pass
$`H(k) = 1 - e^{-k^2 / (2 k_c^2)}`$
- Low Pass
Equation:
$`H(k) = e^{-k^2 / (2 k_c^2)}`$
- Remove Regional
- Continuation
- AGC
- Pseudo Gravity
- Total Horizontal Gradient
- Can now display radially averaged power spectrum (but needs testing!)
   
![SGTools image](dialog.png)    
   
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



