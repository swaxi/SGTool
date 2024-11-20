# SGTool
 Simple Potential Field Calcs to assist WAXI Structural Geophsyics Course https://waxi4.org 
    
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
   
![SGTools image](dialog.png)    
   
# Inputs   
- Supports data geotiff, grd, ers formats

# Tips
- Simple Potential field calculations, mostly FFT-based
- Adds suffix to input filename and stores results in same directory
- If a layer with a given name already exists, no calc performed
- Processing methods preceded by a dot should be performed on line-direction-filtered and reduced to pole/equator data   

# Code development
- Plugin construction - Mark Jessell using QGIS Plugin Builder Plugin    
- Calcs ChatGPT and Mark Jessell
- IGRF calculation - https://github.com/klaundal/ppigrf  
- GRD Loader Mark Jessell & Fatiando a Terra crew https://www.fatiando.org/



