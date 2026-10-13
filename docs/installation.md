[← Back to the README](../README.md)   

# Installation
## QGIS:
1) Either:   
- Download the zip file from the green **<> Code** button and install the zip file in QGIS using the plugin manager for the version on github or   
- Install directly from the QGIS plugin manager from the plugin repository   
   
2) If you get an warning of the type **The following Python packages are required but not installed: scikit-learn** or any other module name the best way to manage this is to install the QGIS Plugin called **qpip** and open it. It will tell you which libraries are missing and allow you to install the correct versions.   
       
   The packages required for specific functions are:   
   **matplotlib** Radial Power Spectrum   
   **scikit-learn** BSDWorms, PCA, ICA   

   If you don't use these functions, there is no need to install the extra packages.   
      

   
## ArcGIS Pro:
1) Download and unzip this respository and store somewhere safe.
2) In order to run python code in ARCGIS Pro you will need to manually add some python libraries so SGTool has everything it needs:
- a) Open Package Manager: Click the Project tab on the ribbon and select Package Manager from the side menu. Access Environment Manager (Gear icon to right of active environemnt): 
  - i) Click the Environment Manager button (top right) to open the management dialog.
  - ii) Clone the Environment: Locate the environment you wish to copy (usually arcgispro-py3) and click the Clone button next to it.
  - iii) Set Destination: In the Clone Environment dialog, provide a name and path for your new environment or leave it as the default.
  - iv) Finalize: Click OK. The cloning process may take several minutes.
- b) Once the cloned environment has been created, search for and install the following packages, **if they are not already installed**:
  - scipy
  - matplotlib
  - scikit-learn
  - pyproj
  - networkx
  - shapely
  - geopandas
- c) In the ArcGIS Pro Catalogue area, go to Add toolbox and select the file **SGTool.pyt** in the **ArcGIS** directory in this repository. 
- d) Double click on the new Geophysical Processing Toolbox to get the list of functions that can be run, and select  **Launch** then **run** to launch the SGTool GUI. Alternatively you can access individual tools from the same Toolbox.    

   
