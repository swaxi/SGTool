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
2) In the ArcGIS Pro Catalog pane, go to **Toolboxes > Add Toolbox** and select the file **SGTool.pyt** in the **ARCGIS** directory in this repository.
3) SGTool needs some Python libraries that ArcGIS Pro does not include. The default environment (arcgispro-py3) cannot be changed, so they go into a clone of it. Either use the setup tool (3A, recommended) or do it by hand (3B).

### 3A) Setup tool
> **Note:** the tool creates a **new** environment, `arcgispro-py3-geotools`, as a copy of `arcgispro-py3`. It does not change `arcgispro-py3` or any other environment you already have. It then makes the new environment the one ArcGIS Pro uses, so **restart ArcGIS Pro** when it has finished. To go back, choose another environment in **Project > Package Manager** and restart ArcGIS Pro. The tool cannot change or delete the environment ArcGIS Pro is using at the time, so to run it again with *Delete the environment first*, switch to another environment and restart first.

- In the SGTool toolbox, open **Setup > Set Up Python Environment** and click **Run**.
- The tool clones arcgispro-py3 into a new environment called **arcgispro-py3-geotools**, installs the missing packages, checks that they all import, and makes the new environment the one ArcGIS Pro uses. It needs an internet connection and takes about 10 minutes. Click **View Details** in the Geoprocessing pane to follow its progress.
- When it has finished, **restart ArcGIS Pro**.
- The same tool is in the **SGTool**, **tomofast_x_q** and **GeoSubSampler** toolboxes and installs what all three need, so you only need to run it once. Running it again only adds what is missing.
- Its options let you choose a different environment name, leave ArcGIS Pro's active environment unchanged (activate it later in Package Manager), or delete the environment and start again from a new clone.
- The tool runs the script `ARCGIS/setup_arcgispro_env.py`, which can also be run from a Windows Command Prompt (change the paths if ArcGIS Pro or this repository are somewhere else):
  ```
  "C:\Program Files\ArcGIS\Pro\bin\Python\envs\arcgispro-py3\python.exe" "C:\path\to\SGTool\ARCGIS\setup_arcgispro_env.py"
  ```

### 3B) By hand
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
  - scikit-image
  - pywavelets
- c) Activate the clone and restart ArcGIS Pro.

4) Double click on the new Geophysical Processing Toolbox to get the list of functions that can be run, and select  **Launch** then **run** to launch the SGTool GUI. Alternatively you can access individual tools from the same Toolbox.    

   
