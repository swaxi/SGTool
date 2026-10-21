[← Back to the README](../README.md)   

# Changelog   

The changelog is also in `metadata.txt`, which QGIS shows in the Plugin Manager.   

## 0.3.8 and earlier
      * BSDWorms is now also a Processing algorithm (Wavelets group): same worms as the dialog, with the points CSV, padded grid and optional polyline shapefile written beside the input and their provenance recorded   
      * Multilevel B-spline gridding: new Auto option (QGIS dialog, Processing algorithm, ArcGIS Pro GUI and toolbox) that chooses the cell size and Max Levels from the spacing of the points; B-spline gridding and MRVBF / MRRTF are several times faster with identical results, and MRVBF works through the grid in bands so a 46-million-cell DEM needs about 3 GB instead of 9.5 GB   
      * Apply the same steps to another grid now repeats local anisotropy, anisotropy chain length and streamline length steps (also new Processing algorithms), and styles the orientation grid; fix the replay dialog failing on histories with steps that cannot be repeated   
      * retain std df/dz for Euler solutions
      * ArcGIS Pro: add a live Preview (FFT Filters, Conv + Stats and Utils/Threshold tabs) in its own window with zoom and pan, opening on the ArcGIS Pro map's extent, following you between tabs, and Apply Preview to save the result at full resolution; fix an error in the power-spectrum error message   
      * ArcGIS Pro: add Multilevel B-Spline gridding (Grid + Wavelets tab and toolbox tool, no Spatial Analyst licence needed); fix the toolbox tools that could not open in current ArcGIS Pro (their parameters were passed to arcpy in the wrong order)   
      * ArcGIS Pro: add MRVBF / MRRTF / slope (Conv + Stats tab and toolbox tool, with the 3-band RGB composite) and remove the DTM curvature classifier (GUI and toolbox) that it supersedes   
      * ArcGIS Pro outputs now use the QGIS file-naming convention (grid_RTP.tif, grid_DirC.tif, grid_d1z.tif ...); the QGIS tilt angle is now saved in degrees (-90 to 90), as in ArcGIS Pro   
      * ArcGIS Pro: grid information (units, CRS, pixels, pixel size) and version in the top bar, Wedge setting and a Remove Line Noise tool, provenance metadata embedded in every output with Read metadata / Save as XML buttons, Apply same steps to another grid (with IGRF recalculated for RTP, RTE and differential RTP), and a Cancel button   
      * Preview: the Keep button is replaced by Apply Preview (the Apply Processing button is renamed Apply Preview while previewing and saves the previewed filter at full resolution); Preview follows you between the FFT Filters, Conv + Stats and Utils tabs, unticking the filters of the tab you leave; the preview layer is always the top layer in the QGIS Layers panel   
      * The layer selected in the QGIS Layers panel and the SGTool Select Grid dropdowns now follow each other   
      * Units, CRS, pixel counts and pixel size of the selected grid, and the version, are shown at the bottom of the first four tabs (cleared when no grids remain)   
      * Cancelling a calculation (task manager X or Cancel) ends the task at once; the dock can be made much narrower and scrolls sideways   
      * Line-noise filter: new Wedge width setting; fix the filter removing too little noise (the opposite spectral lobe was not being filtered); Scale now defaults to 1   
      * Add Apply the same steps to another grid (Utils tab): repeats the processing steps recorded in a grid's history on another grid, recalculating inclination and declination from IGRF for RTP, RTE and differential RTP. Also a Processing algorithm   
      * Improved and debugged Direction Cosine/Butterworth filter: min/max line spacing inputs define the noise band, noise estimate is zero-centred before scaling and subtraction
      * Add live Preview (FFT Filters, Conv + Stats and Utils tabs): one temporary layer that updates as parameters change, for the current map extent or a subsampled grid, with Apply Preview to save the full resolution result
      * Add provenance metadata: how each saved file was made (creation time, source file(s) and their XML, operation and parameters, chained through successive steps) is embedded in saved GeoTIFFs, or in a filename.sgt.xml sidecar for files that cannot carry it
      * Add Read metadata and Save as XML buttons to the Utils tab to show a grid's processing history in a new window or save it as filename.sgt.xml
      * Calculations now run as cancellable QGIS background tasks (Apply Processing, Apply Preview, B-spline gridding, worms, normalise, RGB to grey, grd conversion), with Cancel buttons on each tab
      * Filters exposed as QGIS Processing algorithms (SGTool group) for the modeler, batch processing and scripting
      * Apply Processing moved to a fixed bar at the top of the tabs, next to the Preview controls
      * Fix repeated entries in menus and duplicated signal connections when the plugin is reopened
      * Fix first recalculation failing to overwrite an existing output grid on Windows
      * ArcGIS Pro: Euler CSV includes std df/dz, Euler toolbox tool repaired, directional filter updated to match QGIS
      * Generalise XYZ importer for new variations
      0.3.7
      * Add native python Multilevel BSpline code converted from SAGA
      * Add native python MRVBF/MRRTF/Slope RGB code converted from SAGA
      * Remove DTM Curvature but keep code for now
      * Remove WTMM from GUI but keep code for now
      * Fix clash with modern python for worms calcs
      0.3.6
      * Variable RTP code thanks to Gordon Cooper   
      * Add windowed spatial anisotropy calcs   
      * Add chain length calculations for linear features  
      * Add RGB picker tool to convert LUT picture to grayscale   
      0.3.5   
      * Redo GUI so it removes need for .ui file and is now scrollable
      * redo RTE code following suggestion from google AI   
      * Almost complete ARCGIS Pro version now available with full GUI   
      * Translations of GUI for French, Spanish, Portuguese (Portugal and Brazil) and Mandarin
      0.3.4
      * Fix bug with inclination declination signs for RTE_P now that images always north  
      * Check for commas as decimal separators and replace with decimal points     
      0.3.3
      * change to single decimal versions
      * Don't add mean back for THG calcs
      * Replace raster copy code so it works for negative and positive delta z on linux/macos
      * Add networkx, shapely and geopandas  as requirement so linux is happy
      0.3.01   
      * Fix median removal and resoration for FFT filters   
      0.3.00   
      * Add Spatial stats and convolutions to ArcGIS Pro toolbox
      * Move ArcGIS Pro files to their own directory    
      * Update IGRF to allow up to 2030    
      * Move Euler and IGRF to calcs directory   
      * Compatibility with both QGIS4/QT6 and QGIS3/QT5   
      * Add RGB triplets as LUT definition for RGB to greyscale convertor   

   
Full changelog <a href="https://raw.githubusercontent.com/swaxi/SGTool/refs/heads/main/metadata.txt">Metadata</a>   
