# SGTool v 0.2
 Simple Potential Field Calcs to assist WAXI/Agate Structural Geophysics Course    
 https://waxi4.org   
 https://agate-project.org    
    
- Directional Band Pass
- Band Pass   
- RTP
- RTE
- Derivative
- Analytic Signal   
- Tilt Derivative
- High Pass
- Low Pass
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
- Adds suffix (e.g. _TRP) to input filename and stores results in same directory
- If a layer with a given name already exists, no calc performed
- Processing methods preceded by a dot should be performed on line-direction-filtered and reduced to pole/equator data   
- Length units are defined by grid properties (so Lat/Long wavelengths should be defined in degrees!)

# Code development
- Plugin construction - Mark Jessell using QGIS Plugin Builder Plugin    
- Calcs ChatGPT and Mark Jessell
- IGRF calculation - https://github.com/klaundal/ppigrf  
- GRD Loader Mark Jessell & Fatiando a Terra crew https://www.fatiando.org/



