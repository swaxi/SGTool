[← Back to the README](../README.md)   

# Changelog   

The changelog is also in `metadata.txt`, which QGIS shows in the Plugin Manager.   

## 0.4.0
      * New Petrophys tab (between Grid + Wavelets and Utils), described in [docs/petrophysics.md](petrophysics.md)
      * Apparent susceptibility (SI) from an RTP/RTE magnetic anomaly grid, with the field strength from the IGRF at the grid centre, and apparent density (g/cc) from a Bouguer grid in mGal or µm/s², by dividing out the response of a layer of vertical prisms in the Fourier domain; settings for depth to top, thickness, background value and high-cut; the grids' arbitrary absolute level is not used, so each result is a variation about its background (0.001 SI, 2.67 g/cc); at most one calculation ticked at a time, run by that tab's Apply Processing, with live preview and provenance metadata
      * Henkel plots (density against log susceptibility) of two apparent grids as a heat map, pairing each cell of the coarser grid with the nearest cell of the finer one (reprojected if needed), and per category of a polygon layer (e.g. geology units) in a grid of panels with shared axes and colour scale
      * Measured Petrophysics: any two fields of a point layer (numeric, or text holding numbers) plotted against each other, each with a multiplier and its own log option, as points with up to 64 category symbols (8 colours x 8 markers, by sample count) or as a heat map, with axes and legend named after the fields
      * With only one apparent grid or measured field chosen, histograms instead, per category if one is chosen, with the median as a red line, median, mean and standard deviation as text, and the vertical axis clipped to the tallest bar
      * Plots drawn by their own buttons, each in its own window, optionally saved as PNG and/or PDF; values <= 0 left out of log axes (and counted in the title); heat maps in viridis
      * New script scripts/pivot_petrophysics.py: one point per sample from a petrophysics layer with one row per measurement, for Measured Petrophysics
      * ArcGIS Pro: the same Petrophys tab in the GUI (with preview, IGRF field strength and all the plots, polygons and samples chosen as files), and Apparent Susceptibility and Apparent Density toolbox tools
      * Group similar: categories with alike distributions (Jensen-Shannon distance, average-linkage clustering) grouped automatically or into 2 to 8 groups, lettered A, B, C ... on the per-category plots (title, coloured frame) and in the sample plot legend
      * Apparent susceptibility and density are also QGIS Processing algorithms (Petrophysics group), and Apply same steps to another grid repeats them in both versions
      * macOS: the dock uses a smaller font and the compact Fusion style, so its layout is no larger than on Windows (set the QSettings value SGTool/compactMacUI to false to keep the native look)
      * New Low-latitude RTP (RTE/P choice "Low-lat. RTP", output _RTP_LL): stable reduction to the pole at low magnetic latitudes by the Taylor-series iterative method of Hao et al. (2018), reducing to the equator and then to the pole with m iterations; m chosen automatically (0 to 200, best correlation with the analytic signal) or set in Iterations. In the QGIS dialog, Processing, replay, and the ArcGIS Pro GUI and toolbox
      * Low-latitude RTP: Stripe suppression setting (default 10 degrees, 0 = off) fades the iterations out near the direction perpendicular to the declination, where noise, edge effects and levelling errors otherwise show as stripes parallel to the declination; the dialog reports the iterations used
      * Fix reduction to the equator: the operator's denominator was not squared, so at inclination 0 the grid was multiplied by cos^2(D - theta) instead of being left unchanged, and at other inclinations the anomalies were too weak. RTE results will differ from earlier versions
      * RTE and Low-latitude RTP pad the grid by mirroring its edges, which stops stripes from the north and south edges along the declination
      * RTP, RTE and Low-latitude RTP warn and stop again when Inc and Dec are both 0 (the check never fired); in the ArcGIS Pro toolbox the FFT buffer is now passed to RTP and RTE (it was read but not used)
      * Chain Length and Streamline Length are hidden (QGIS dialog, Processing and the ArcGIS Pro GUI) until they work well; the calculation code is kept
      * Dock layout: every tab's groups, grid selectors and file rows share one width (so "..." buttons stay in view), the scroll bars are on the left, and the band pass / high-low pass width and AGC window are on one line with their filters
      * Petrophys: changing the survey date no longer recalculates the field strength; it updates when IGRF is pressed or a grid is selected
      * Fix an error when the plugin is unloaded after QGIS has already removed its Processing provider (seen when uninstalling on macOS)
      * Docs: every documented function gives its primary reference and an "In SGTool" note on how the code relates to the published method, with a reference list per page; equations corrected where they did not match the code (RTP, downward continuation, AGC, derivatives, directional line-noise filter)
      * The grid that was processed stays selected (Layers panel and dropdowns) after Apply Processing, Apply Preview and the other calculations, instead of the new result; the RTE/P dropdown is wide enough to show "Low-lat. RTP"

## 0.3.9
      * Utils tab: new Clip to Window button clips the selected grid to the current map window, saving it beside the original as <name>_clip.tif (then _clip_01, _clip_02 ... rather than overwriting an earlier clip); Clip All Layers is now half width, beside it. Also in the ArcGIS Pro GUI (Utils tab), using the ArcGIS Pro map view
      * Crop to area prior to preview to speed up on large grids
      * IGRF button no longer loads the whole grid just to find its centre
      * Faster NaN filling and multithreaded FFT for Fourier-based filters (RTP etc.) on large grids
      * Apply now reads the grid inside the background task, so QGIS stays responsive and Cancel works while a large grid loads
      * Security and Qt6 compatibility fixes for the plugin repository checks (no silent try/except, safer XML parsing, fully scoped Qt/QGIS enums)

## 0.3.8 and earlier
      * Utils tab: new Clip All Layers tool (QGIS dialog, and the ArcGIS Pro toolbox tool Clip All Layers to Polygon): clips every layer of the project to a chosen polygon layer (or its selected polygons), saving each as <name>_clip with its .qml style in a chosen folder, plus a copy of the project (.qgz) with the same groups, order, visibility and styles that loads the clipped layers   
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
