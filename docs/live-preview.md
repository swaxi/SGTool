[← Back to the README](../README.md)   

# Live Preview (QGIS)   

(ArcGIS Pro has a preview too, in its own window: see [ArcGIS Pro differences](arcgis-pro.md).)   

The **FFT Filters**, **Conv + Stats** and **Utils** (Threshold to NaN) tabs have a fixed bar at the top with **Apply Processing** (which reads **Apply Preview** while previewing) on the left and the preview controls to its right, so they stay visible while you scroll the filter list.   

- **Preview**: tick it to see the selected filter without writing any files. A single temporary layer called `<grid><suffix>_preview` is added to the project and is replaced in place each time a parameter changes, so layers do not accumulate. Updates wait a fraction of a second after the last edit, so typing a value does not trigger a calculation for every digit.   
- **One filter at a time**: while previewing, ticking a second filter unticks the first, and editing a filter's parameters selects that filter. Filters that cannot be previewed (Differential RTP, local anisotropy, MRVBF, PCA, ICA, Euler deconvolution, and the clipping polygon) are disabled until Preview is switched off.   
- **Map extent**: calculates at full resolution for the area currently shown in the map canvas, so panning or zooming updates the preview. Filter edge effects appear along the edges of the view, and very large views are block-averaged to at most 2000 cells on the longest side.   
- **Subsampled grid**: calculates for the whole grid after block-averaging it to at most 600 cells on the longest side. This is quick, but filters set in pixels (convolutions, window statistics, AGC) then act on the coarser cells, and line noise cannot be seen if the preview cell size exceeds the line spacing. Use Map extent for those.   
- **Apply Preview**: while Preview is ticked the **Apply Processing** button is renamed **Apply Preview** and stays available. Pressing it calculates the previewed filter at full resolution, adds it as a normal permanent layer (with the usual name and suffix), switches preview off and unticks the filters.   
- **Unticking Preview** (or closing the plugin) discards the temporary layer. Unticking the last ticked filter clears the image but keeps preview mode on, so ticking a filter brings it straight back. **Changing tab** unticks the filters of the tab you leave and removes the image, and Preview is ticked on the next tab that has one (Grid + Wavelets has none, so Preview waits and returns afterwards); turning Preview off yourself ends it for good. While previewing, the Apply button reads **Apply Preview** and saves the previewed filter at full resolution.   
- **Calculated in the background**: QGIS stays usable while a preview is calculated, and the image appears when it is ready. **Cancel** stops it (useful after zooming out over a big grid), panning, zooming or editing a setting replaces a preview still being calculated, and starting another job replaces it and redraws it afterwards.
- The preview layer is shown with bilinear resampling when zoomed in to reduce the blocky look of a coarse preview.   
