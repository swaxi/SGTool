# -*- coding: utf-8 -*-
"""
/***************************************************************************
 SGTool
                                 A QGIS plugin
 Simple Potential Field Processing
 Generated by Plugin Builder: http://g-sherman.github.io/Qgis-Plugin-Builder/
                              -------------------
        begin                : 2024-11-17
        git sha              : $Format:%H$
        copyright            : (C) 2024 by Mark Jessell
        email                : mark.jessell@uwa.edu.au
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

from qgis.PyQt.QtCore import QSettings, QTranslator, QCoreApplication, Qt
from qgis.PyQt.QtGui import QIcon
from qgis.PyQt.QtWidgets import QAction
from qgis.core import (
    Qgis,
    QgsCoordinateReferenceSystem,
    QgsCoordinateTransform,
    QgsVectorLayer,
    QgsProject,
    QgsRasterLayer,
    QgsFeature,
    QgsField,
    QgsVectorFileWriter,
    QgsPoint,
)
from qgis.PyQt.QtWidgets import QAction, QFileDialog
from qgis.PyQt.QtCore import (
    QSettings,
    QTranslator,
    QCoreApplication,
    QFileInfo,
    QVariant,
    Qt,
)
from qgis.core import QgsRasterLayer, QgsRectangle

# Initialize Qt resources from file resources.py
from .resources import *
import ntpath

# Import the code for the DockWidget
from .SGTool_dockwidget import SGToolDockWidget
from .GeophysicalProcessor import GeophysicalProcessor
from .geosoft_grid_parser import * 

import os.path
import numpy as np
import pandas as pd
from osgeo import gdal, osr
from .ppigrf import igrf, get_inclination_declination
from .ppigrf import igrf, get_inclination_declination
from datetime import datetime
from pyproj import Transformer

class SGTool:
    """QGIS Plugin Implementation."""

    def __init__(self, iface):
        """Constructor.

        :param iface: An interface instance that will be passed to this class
            which provides the hook by which you can manipulate the QGIS
            application at run time.
        :type iface: QgsInterface
        """
        # Save reference to the QGIS interface
        self.iface = iface

        # initialize plugin directory
        self.plugin_dir = os.path.dirname(__file__)

        # initialize locale
        locale = QSettings().value('locale/userLocale')[0:2]
        locale_path = os.path.join(
            self.plugin_dir,
            'i18n',
            'SGTool_{}.qm'.format(locale))

        if os.path.exists(locale_path):
            self.translator = QTranslator()
            self.translator.load(locale_path)
            QCoreApplication.installTranslator(self.translator)

        # Declare instance attributes
        self.actions = []
        self.menu = self.tr(u'&SGTool')
        # TODO: We are going to let the user set this up in a future iteration
        self.toolbar = self.iface.addToolBar(u'SGTool')
        self.toolbar.setObjectName(u'SGTool')

        #print "** INITIALIZING SGTool"

        self.pluginIsActive = False
        self.dlg = None



    # noinspection PyMethodMayBeStatic
    def tr(self, message):
        """Get the translation for a string using Qt translation API.

        We implement this ourselves since we do not inherit QObject.

        :param message: String for translation.
        :type message: str, QString

        :returns: Translated version of message.
        :rtype: QString
        """
        # noinspection PyTypeChecker,PyArgumentList,PyCallByClass
        return QCoreApplication.translate('SGTool', message)


    def add_action(
        self,
        icon_path,
        text,
        callback,
        enabled_flag=True,
        add_to_menu=True,
        add_to_toolbar=True,
        status_tip=None,
        whats_this=None,
        parent=None):
        """Add a toolbar icon to the toolbar.

        :param icon_path: Path to the icon for this action. Can be a resource
            path (e.g. ':/plugins/foo/bar.png') or a normal file system path.
        :type icon_path: str

        :param text: Text that should be shown in menu items for this action.
        :type text: str

        :param callback: Function to be called when the action is triggered.
        :type callback: function

        :param enabled_flag: A flag indicating if the action should be enabled
            by default. Defaults to True.
        :type enabled_flag: bool

        :param add_to_menu: Flag indicating whether the action should also
            be added to the menu. Defaults to True.
        :type add_to_menu: bool

        :param add_to_toolbar: Flag indicating whether the action should also
            be added to the toolbar. Defaults to True.
        :type add_to_toolbar: bool

        :param status_tip: Optional text to show in a popup when mouse pointer
            hovers over the action.
        :type status_tip: str

        :param parent: Parent widget for the new action. Defaults None.
        :type parent: QWidget

        :param whats_this: Optional text to show in the status bar when the
            mouse pointer hovers over the action.

        :returns: The action that was created. Note that the action is also
            added to self.actions list.
        :rtype: QAction
        """

        icon = QIcon(icon_path)
        action = QAction(icon, text, parent)
        action.triggered.connect(callback)
        action.setEnabled(enabled_flag)

        if status_tip is not None:
            action.setStatusTip(status_tip)

        if whats_this is not None:
            action.setWhatsThis(whats_this)

        if add_to_toolbar:
            self.toolbar.addAction(action)

        if add_to_menu:
            self.iface.addPluginToMenu(
                self.menu,
                action)

        self.actions.append(action)

        return action


    def initGui(self):
        """Create the menu entries and toolbar icons inside the QGIS GUI."""

        icon_path = ':/plugins/SGTool/icon.png'
        self.add_action(
            icon_path,
            text=self.tr(u'SGTool'),
            callback=self.run,
            parent=self.iface.mainWindow())

    #--------------------------------------------------------------------------

    def onClosePlugin(self):
        """Cleanup necessary items here when plugin dockwidget is closed"""

        #print "** CLOSING SGTool"

        # disconnects
        self.dlg.closingPlugin.disconnect(self.onClosePlugin)

        # remove this statement if dockwidget is to remain
        # for reuse if plugin is reopened
        # Commented next statement since it causes QGIS crashe
        # when closing the docked window:
        # self.dlg = None

        self.pluginIsActive = False


    def unload(self):
        """Removes the plugin menu item and icon from QGIS GUI."""

        #print "** UNLOAD SGTool"

        for action in self.actions:
            self.iface.removePluginMenu(
                self.tr(u'&SGTool'),
                action)
            self.iface.removeToolBarIcon(action)
        # remove the toolbar
        del self.toolbar

    def initParams(self):
        self.localGridName=""
        self.diskGridPath=""
        self.diskPointsPath=""
        self.buffer=0
        self.DirClean=False
        self.DC_azimuth=0
        self.DC_lineSpacing=400
        self.RTE_P=False
        self.RTE_P_type="Pole"
        self.RTE_P_inc=0
        self.RTE_P_dec=0
        self.RTE_P_height=0
        self.RTE_P_date=[1,1,2000]
        self.RemRegional=False
        self.remReg_wavelength=5000
        self.Derivative=False
        self.derive_direction="z"
        self.derive_power=1.0
        self.TDR=False
        self.AS=False
        self.Continuation=False
        self.cont_direction="up"
        self.cont_height=500
        self.BandPass=False
        self.band_low=5
        self.band_high=50
        self.AGC=False
        self.agc_window=10
        self.FreqCut=False
        self.FreqCut_type="Low"
        self.FreqCut_cut=1000

    def parseParams(self):
        self.localGridName=self.dlg.mMapLayerComboBox_selectGrid.currentText()

        self.DirClean=self.dlg.checkBox_3_DirClean.isChecked()
        self.DC_azimuth=self.dlg.lineEdit_3_azimuth.text()
        self.DC_lineSpacing=self.dlg.lineEdit_4_lineSpacing.text()

        self.RTE_P=self.dlg.checkBox_4_RTE_P.isChecked()
        self.RTE_P_type=self.dlg.comboBox_3_rte_p_list.currentText()
        self.RTE_P_inc=self.dlg.lineEdit_6_inc.text()
        self.RTE_P_dec=self.dlg.lineEdit_5_dec.text()
        self.RTE_P_height=self.dlg.lineEdit_7_height.text()
        date_text = str(self.dlg.dateEdit.date().toPyDate())
        date_split = date_text.split("-")        
        self.RTE_P_date=[int(date_split[2]),int(date_split[1]),int(date_split[0])]


        self.RemRegional=self.dlg.checkBox_5_regional.isChecked()
        self.remReg_wavelength=self.dlg.lineEdit_9_removeReg_wavelength.text()
        
        self.Derivative=self.dlg.checkBox_6_derivative.isChecked()
        self.derive_direction=self.dlg.comboBox_derivDirection.currentText()
        self.derive_power=self.dlg.lineEdit_9_derivePower.text()

        self.TDR=self.dlg.checkBox_7_tiltDerivative.isChecked()

        self.AS=self.dlg.checkBox_8_analyticSignal.isChecked()
        
        self.Continuation=self.dlg.checkBox_9_continuation.isChecked()
        self.cont_direction=self.dlg.comboBox_2_continuationDirection.currentText()
        self.cont_height=self.dlg.lineEdit_10_continuationHeight.text()

        self.BandPass=self.dlg.checkBox_10_bandPass.isChecked()
        self.band_low=self.dlg.lineEdit_12_bandPassLow.text()
        if(float(self.band_low)<=0.0):
            self.band_low=1e-10
        self.band_high=self.dlg.lineEdit_11_bandPassHigh.text()

        self.AGC=self.dlg.checkBox_11_1vd_agc.isChecked()
        self.agc_window=self.dlg.lineEdit_13_agc_window.text()

        self.FreqCut=self.dlg.checkBox_10_freqCut.isChecked()
        self.FreqCut_type=self.dlg.comboBox_2_FreqCutType.currentText()
        self.FreqCut_cut=self.dlg.lineEdit_12_FreqPass.text()

    def loadGrid(self):
        fileInfo = QFileInfo(self.diskGridPath)
        baseName = fileInfo.baseName()

        self.layer = QgsRasterLayer(self.diskGridPath, baseName)
        if(not self.is_layer_loaded(baseName)):
            QgsProject.instance().addMapLayer(self.layer)

        self.dx = self.layer.rasterUnitsPerPixelX()
        self.dy = self.layer.rasterUnitsPerPixelY()
        # Access the raster data provider
        provider = self.layer.dataProvider()

        # Get raster dimensions
        cols = provider.xSize()  # Number of columns
        rows = provider.ySize()  # Number of rows

        # Read raster data as a block
        band = 1  # Specify the band number (1-based index)
        raster_block = provider.block(band, provider.extent(), cols, rows)

        # Copy the block data into a NumPy array
        extent = self.layer.extent()    
        rows, cols = self.layer.height(), self.layer.width()    
        raster_block = provider.block(1, extent, cols, rows)  # !!!!!  
        self.raster_array = np.zeros((rows, cols))    
        for i in range(rows):    
            for j in range(cols):    
                self.raster_array[i,j] = raster_block.value(i,j)    

        # Handle NoData values if needed
        no_data_value = provider.sourceNoDataValue(1)  # Band 1

        if no_data_value is not None:
            self.raster_array[self.raster_array == no_data_value] = np.nan

    def insert_text_before_extension(self,file_path, insert_text):
        """
        Insert text at the end of the filename, before the file extension.

        Parameters:
            file_path (str): Full path of the file.
            insert_text (str): Text to insert before the file extension.

        Returns:
            str: The modified file path.
        """
        # Separate the file path into directory, base name, and extension
        dir_name, base_name = os.path.split(file_path)
        file_name, file_ext = os.path.splitext(base_name)

        # Construct the new file name
        new_file_name = f"{file_name}{insert_text}{file_ext}"

        # Combine directory and new file name
        return os.path.join(dir_name, new_file_name)

    def procDirClean(self):
        self.new_grid=self.processor.directional_band_reject_filter(self.raster_array, low_cut=0.9*float(self.DC_lineSpacing), high_cut=1.1*float(self.DC_lineSpacing), direction_angle=float(self.DC_azimuth),buffer_size=self.buffer)
        self.new_grid = self.raster_array-self.new_grid
        self.suffix="_DC"

    def procRTP_E(self):
        if(self.RTE_P_type=="Pole"):
            self.new_grid=self.processor.reduction_to_pole(self.raster_array, inclination=float(self.RTE_P_inc), declination=float(self.RTE_P_dec),buffer_size=self.buffer)
            self.suffix="_RTP"
        else:
            self.new_grid=self.processor.reduction_to_equator(self.raster_array, inclination=float(self.RTE_P_inc), declination=float(self.RTE_P_dec),buffer_size=self.buffer)
            self.suffix="_RTE"

    def procRemRegional(self):
        self.new_grid = self.processor.remove_regional_trend_fourier(self.raster_array, cutoff_wavelength=float(self.remReg_wavelength),buffer_size=self.buffer)
        self.new_grid = self.raster_array-self.new_grid
        self.suffix="_RR"+"_"+str(self.remReg_wavelength)

    def procDerivative(self):
        self.new_grid=self.processor.compute_derivative(self.raster_array, direction=self.derive_direction, order=float(self.derive_power),buffer_size=self.buffer)   
        self.suffix="_d"+str(self.derive_power)+self.derive_direction

    def procTiltDerivative(self):
        self.new_grid=self.processor.tilt_derivative(self.raster_array,buffer_size=self.buffer)
        self.suffix="_TDR"

    def procAnalyticSignal(self):
        self.new_grid=self.processor.analytic_signal(self.raster_array,buffer_size=self.buffer)
        self.suffix="_AS"

    def procContinuation(self):
        if(self.cont_direction=="up"):
            self.new_grid=self.processor.upward_continuation(self.raster_array, height=float(self.cont_height),buffer_size=self.buffer)
            self.suffix="_UC"+"_"+str(self.cont_height)
        else:
            self.new_grid=self.processor.downward_continuation(self.raster_array, height=float(self.cont_height),buffer_size=self.buffer)
            self.suffix="_DC"+"_"+str(self.cont_height)

    def procBandPass(self):
        self.new_grid=self.processor.band_pass_filter(self.raster_array, low_cut=float(self.band_low), high_cut=float(self.band_high),buffer_size=self.buffer)
        self.suffix="_BP"+"_"+str(self.band_low)+"_"+str(self.band_high)

    def procAGC(self):
        self.new_grid=self.processor.automatic_gain_control(self.raster_array, window_size=float(self.agc_window))
        self.suffix="_AGC"

    def procFreqCut(self):
        if(self.FreqCut_type=="Low"):
            self.new_grid=self.processor.low_pass_filter(self.raster_array, cutoff_wavelength=float(self.FreqCut_cut),buffer_size=self.buffer)
            self.suffix="_LP"+"_"+str(self.FreqCut_cut)
        else:
            self.new_grid=self.processor.high_pass_filter(self.raster_array, cutoff_wavelength=float(self.FreqCut_cut),buffer_size=self.buffer)
            self.suffix="_HP"+"_"+str(self.FreqCut_cut)

    def gridPoints(self):
        dx=float(self.dlg.lineEdit_13_MC_gridSize.text())
        xmin=self.pointData[self.dlg.comboBox_grid_x.currentText()].min()
        xmax=self.pointData[self.dlg.comboBox_grid_x.currentText()].max()
        ymin=self.pointData[self.dlg.comboBox_grid_y.currentText()].min()
        ymax=self.pointData[self.dlg.comboBox_grid_y.currentText()].max()

        num_points_x=int((xmax-xmin)/dx)
        num_points_y=int((ymax-ymin)/dx)

        grid_x, grid_y = np.meshgrid(np.linspace(xmin, xmax, num_points_x), np.linspace(ymin, ymax, num_points_y))
        self.buffer=min(num_points_y,num_points_x)


        self.processor = GeophysicalProcessor(float(self.dlg.lineEdit_13_MC_gridSize.text()), float(self.dlg.lineEdit_13_MC_gridSize.text()),self.buffer)

        x=self.pointData[self.dlg.comboBox_grid_x.currentText()]
        y=self.pointData[self.dlg.comboBox_grid_y.currentText()]
        z=self.pointData[self.dlg.comboBox_grid_data.currentText()]


        # Ensure inputs are iterable
        self.gridded=self.processor.minimum_curvature_gridding(x,y,z, grid_x,grid_y)

        raster_path=os.path.dirname(self.diskPointsPath)+"/"+os.path.splitext(os.path.basename(self.diskPointsPath))[0]+"_"+self.dlg.comboBox_grid_data.currentText()+".tif"

        self.numpy_array_to_raster(self.gridded, raster_path=raster_path, dx=dx,xmin=xmin,ymax=ymax,reference_layer=None, no_data_value=np.nan)
        gridded_layer = QgsRasterLayer(raster_path, os.path.splitext(os.path.basename(self.diskPointsPath))[0]+"_"+self.dlg.comboBox_grid_data.currentText())
        if gridded_layer.isValid():
            QgsProject.instance().addMapLayer(gridded_layer)
    
    def addNewGrid(self):
        if(not self.is_layer_loaded(self.base_name+self.suffix)):
            self.diskNewGridPath = self.insert_text_before_extension(self.diskGridPath, self.suffix)
            self.numpy_array_to_raster(self.new_grid, self.diskNewGridPath, dx=None,xmin=None,ymax=None,reference_layer=self.layer, no_data_value=np.nan)
            con_raster_layer = QgsRasterLayer(self.diskNewGridPath, self.base_name+self.suffix)
            if con_raster_layer.isValid():
                QgsProject.instance().addMapLayer(con_raster_layer)
    
    def processGeophysics(self):
        self.localGridName=self.dlg.mMapLayerComboBox_selectGrid.currentText()
        process=False
        if(os.path.exists(self.diskGridPath) and self.diskGridPath!=""):
            self.parseParams()
            self.loadGrid()

            paths = os.path.split(self.diskGridPath)
            self.base_name = "".join(paths[1].split(".")[:-1])
            provider = self.layer.dataProvider()

            # Get raster dimensions
            cols = provider.xSize()  # Number of columns
            rows = provider.ySize()  # Number of rows
            process=True

        elif(self.localGridName and self.localGridName!=""):
            self.parseParams()
            self.layer=QgsProject.instance().mapLayersByName(self.localGridName)[0]
            self.base_name = self.localGridName
            self.diskGridPath=self.layer.dataProvider().dataSourceUri()
            self.dx = self.layer.rasterUnitsPerPixelX()
            self.dy = self.layer.rasterUnitsPerPixelY()
            # Access the raster data provider
            provider = self.layer.dataProvider()

            # Get raster dimensions
            cols = provider.xSize()  # Number of columns
            rows = provider.ySize()  # Number of rows

            # Read raster data as a block
            band = 1  # Specify the band number (1-based index)
            raster_block = provider.block(band, provider.extent(), cols, rows)

            # Copy the block data into a NumPy array
            extent = self.layer.extent()    
            rows, cols = self.layer.height(), self.layer.width()    
            raster_block = provider.block(1, extent, cols, rows)  # !!!!!  
            self.raster_array = np.zeros((rows, cols))    
            for i in range(rows):    
                for j in range(cols):    
                    self.raster_array[i,j] = raster_block.value(i,j)    

            # Handle NoData values if needed
            no_data_value = provider.sourceNoDataValue(1)  # Band 1

            if no_data_value is not None:
                self.raster_array[self.raster_array == no_data_value] = np.nan
            process=True

        if(process):
            self.buffer=min(rows,cols)
            self.processor = GeophysicalProcessor(self.dx, self.dy,self.buffer)

            if(self.DirClean):
                self.procDirClean()
                self.addNewGrid()
            if(self.RTE_P):
                self.procRTP_E()
                self.addNewGrid()
            if(self.RemRegional):
                self.procRemRegional()
                self.addNewGrid()
            if(self.Derivative):
                self.procDerivative()
                self.addNewGrid()
            if(self.TDR):
                self.procTiltDerivative()
                self.addNewGrid()
            if(self.AS):
                self.procAnalyticSignal()
                self.addNewGrid()
            if(self.Continuation):
                self.procContinuation()
                self.addNewGrid()
            if(self.BandPass):
                self.procBandPass()
                self.addNewGrid()
            if(self.FreqCut):
                self.procFreqCut()
                self.addNewGrid()
            if(self.AGC):
                self.procAGC()
                self.addNewGrid()

    def is_layer_loaded(self,layer_name):
        """
        Check if a layer with the specified name is already loaded in QGIS.

        Parameters:
            layer_name (str): The name of the layer to check.

        Returns:
            bool: True if the layer is loaded, False otherwise.
        """
        for layer in QgsProject.instance().mapLayers().values():
            if layer.name() == layer_name:
                return True
        return False

    def select_grid_file(self):

        self.diskGridPath, _filter = QFileDialog.getOpenFileName(
            None,
            "Select Data File",
            ".",
            "TIF (*.TIF;*.tif;*.TIFF;*.tiff);;GRD (*.grd;*GRD);;ERS (*.ERS;*.ers)",
        )
        suffix = self.diskGridPath.split(".")[-1].lower()
        epsg="4326"

        if os.path.exists(self.diskGridPath) and self.diskGridPath != "":
            self.dlg.lineEdit_2_loadGridPath.setText(self.diskGridPath)
            self.dlg.pushButton_3_applyProcessing.setEnabled(True)
        if(os.path.exists(self.diskGridPath+'.xml') and suffix=="grd"):
            epsg=extract_proj_str(self.diskGridPath+'.xml')
            if(epsg== None):
                epsg=4326
                self.iface.messageBar().pushMessage("No CRS found in XML, default to 4326", level=Qgis.Warning, duration=15)
            else:
                self.iface.messageBar().pushMessage("CRS Read from XML as "+epsg, level=Qgis.Info, duration=15)
            #self.dlg.mQgsProjectionSelectionWidget.setCrs(QgsCoordinateReferenceSystem('EPSG:'+str(epsg)))
            self.save_a_grid(epsg)
        elif( suffix=="tif"):
            basename =os.path.basename(self.diskGridPath)
            filename_without_extension =os.path.splitext(basename)[0]

            self.layer = QgsRasterLayer(self.diskGridPath, filename_without_extension)
            if(not self.is_layer_loaded(self.diskGridPath)):
                QgsProject.instance().addMapLayer(self.layer)

    #save grd file as geotiff
    def save_a_grid(self,epsg):


        #load grd file and store in memory
        if(self.diskGridPath !=''): 
            if(not os.path.exists(self.diskGridPath)):
                self.iface.messageBar().pushMessage("File: "+self.diskGridPath+" not found", level=Qgis.Warning, duration=3)
            else:    
                grid,header,Gdata_type=load_oasis_montaj_grid(self.diskGridPath)
                directory_path = os.path.dirname(self.diskGridPath)
                basename =os.path.basename(self.diskGridPath)
                filename_without_extension =os.path.splitext(basename)[0]
                self.diskGridPath=directory_path+"/"+filename_without_extension+".tif"

                fn=self.diskGridPath

                driver=gdal.GetDriverByName('GTiff')
                if(header["ordering"]==1):
                    ds = driver.Create(fn,xsize=header["shape_e"],ysize=header["shape_v"],bands=1,eType=Gdata_type)
                else:
                    ds = driver.Create(fn,xsize=header["shape_v"],ysize=header["shape_e"],bands=1,eType=Gdata_type)

                ds.GetRasterBand(1).WriteArray(grid)
                geot=[header["x_origin"]-(header["spacing_e"]/2),
                    header["spacing_e"],
                    0,
                    header["y_origin"]-(header["spacing_v"]/2),
                    0,
                    header["spacing_e"],
                    ]
                ds.SetGeoTransform(geot)
                srs=osr.SpatialReference()
                srs.ImportFromEPSG(int(epsg))
                ds.SetProjection(srs.ExportToWkt())
                ds=None

                self.layer = QgsRasterLayer(self.diskGridPath, filename_without_extension+".tif")
                if(not self.is_layer_loaded(filename_without_extension+".tif")):
                    QgsProject.instance().addMapLayer(self.layer)

                #self.iface.messageBar().pushMessage("GRD saved to file", level=Qgis.Success, duration=5)

        else:
            self.iface.messageBar().pushMessage("You need to select a file first", level=Qgis.Warning, duration=3)

    def select_point_file(self):

        self.diskPointsPath, _filter = QFileDialog.getOpenFileName(
            None,
            "Select Data File",
            ".",
            "CSV (*.csv;*.txt;*.CSV;*.TXT)",
        )
        if os.path.exists(self.diskPointsPath) and self.diskPointsPath != "":

            self.pointData = pd.read_csv(self.diskPointsPath)

            self.dlg.lineEdit_loadPointsPath.setText(self.diskPointsPath)
            self.dlg.comboBox_grid_x.setEnabled(True)
            self.dlg.comboBox_grid_y.setEnabled(True)
            self.dlg.comboBox_grid_data.setEnabled(True)
            self.dlg.comboBox_grid_x.addItems(self.pointData.columns)
            self.dlg.comboBox_grid_x.setCurrentIndex(1)
            self.dlg.comboBox_grid_y.addItems(self.pointData.columns)
            self.dlg.comboBox_grid_y.setCurrentIndex(2)
            self.dlg.comboBox_grid_data.addItems(self.pointData.columns)


    def numpy_array_to_raster(self,numpy_array, raster_path, dx=None,xmin=None,ymax=None,reference_layer=None, no_data_value=np.nan):
        """
        Convert a NumPy array to a GeoTIFF raster file.
        
        Parameters:
            numpy_array (numpy.ndarray): The NumPy array to convert.
            raster_path (str): The path to save the raster file.
            reference_layer (QgsRasterLayer, optional): A reference layer for CRS and geotransform.
            no_data_value: Value to use for no data (default is NaN).
        """

        # Check if the file already exists and remove it
        if os.path.exists(raster_path):
            os.remove(raster_path)

        rows, cols = numpy_array.shape
        driver = gdal.GetDriverByName('GTiff')
        output_raster = driver.Create(raster_path, cols, rows, 1, gdal.GDT_Float32)

        # Set geotransform and projection if a reference layer is provided
        if reference_layer:
            provider = reference_layer.dataProvider()
            extent = provider.extent()
            geotransform = [
                extent.xMinimum(),
                extent.width() / cols,  # pixel width
                0,
                extent.yMaximum(),
                0,
                -extent.height() / rows  # pixel height (negative)
            ]
            output_raster.SetGeoTransform(geotransform)

            # Set CRS
            srs = osr.SpatialReference()
            srs.ImportFromWkt(reference_layer.crs().toWkt())
            output_raster.SetProjection(srs.ExportToWkt())
        else:
            crs=self.dlg.mQgsProjectionSelectionWidget.crs().authid()
            srs = osr.SpatialReference()
            srs.ImportFromEPSG(int(crs.split(":")[1]))
            output_raster.SetProjection(srs.ExportToWkt())
            geotransform = [
                xmin,
                dx,  # pixel width
                0,
                ymax,
                0,
                -dx  # pixel height (negative)
            ]
            output_raster.SetGeoTransform(geotransform)

        # Write data to raster
        band = output_raster.GetRasterBand(1)
        if no_data_value is not None:
            band.SetNoDataValue(no_data_value)
        numpy_array = np.nan_to_num(numpy_array, nan=no_data_value)  # Replace NaN with no_data_value
        band.WriteArray(numpy_array)
        band.FlushCache()
        output_raster = None  # Close the file

    # estimate mag field from centroid of data, date and sensor height
    def update_mag_field(self):
        self.localGridName=self.dlg.mMapLayerComboBox_selectGrid.currentText()
        if(os.path.exists(self.diskGridPath) or  self.localGridName):
            if(os.path.exists(self.diskGridPath)):
                self.loadGrid()

            else:
                self.layer=QgsProject.instance().mapLayersByName(self.localGridName)[0]
                self.base_name = self.localGridName

            # retrieve parameters
            self.magn_SurveyHeight = self.dlg.lineEdit_7_height.text()
            date_text = str(self.dlg.dateEdit.date().toPyDate())

            date_split = date_text.split("-")
            self.magn_SurveyDay = int(date_split[2])
            self.magn_SurveyMonth = int(date_split[1])
            self.magn_SurveyYear = int(date_split[0])
            date = datetime(self.magn_SurveyYear, self.magn_SurveyMonth, self.magn_SurveyDay)

            extent = self.layer.extent()  # Get the extent of the raster layer

            # calculate midpoint of grid
            midx = extent.xMinimum()+(extent.xMaximum()-extent.xMinimum())
            midy = extent.yMinimum()+(extent.yMaximum()-extent.yMinimum())

            # convert midpoint to lat/long
            magn_proj = self.layer.crs().authid().split(":")[1]

            proj = Transformer.from_crs(magn_proj, 4326, always_xy=True)
            x, y = (midx, midy)
            long, lat = proj.transform(x, y)

            # calculate IGRF compnents and  convert to Inc, Dec, Int
            Be, Bn, Bu = igrf(
                long, lat, float(self.magn_SurveyHeight), date
            )  # returns east, north, up
            (
                self.RTE_P_inc,
                self.RTE_P_dec
            ) = get_inclination_declination(Be, Bn, Bu, degrees=True)
            self.forward_magneticField_intensity = np.sqrt(Be**2 + Bn**2 + Bu**2)

            self.RTE_P_inc=self.RTE_P_inc.item()
            self.RTE_P_dec=self.RTE_P_dec.item()

            # update widgets
            self.dlg.lineEdit_5_dec.setText(
                str(round(self.RTE_P_dec,1))
            )
            self.dlg.lineEdit_6_inc.setText(
                str(round(self.RTE_P_inc,1))
            )
            """self.dlg.doubleSpinBox_mag_int.setValue(
                self.forward_magneticField_intensity.item()
            )"""
    
    def update_paths(self):
        self.localGridName=self.dlg.mMapLayerComboBox_selectGrid.currentText()
        self.dlg.lineEdit_2_loadGridPath.setText("")
        self.diskGridPath=""
        self.base_name = self.localGridName

    #--------------------------------------------------------------------------

    def run(self):
        """Run method that loads and starts the plugin"""

        if not self.pluginIsActive:
            self.pluginIsActive = True
            self.initParams()
            #print "** STARTING SGTool"

            # dockwidget may not exist if:
            #    first run of plugin
            #    removed on close (see self.onClosePlugin method)
            if self.dlg == None:
                # Create the dockwidget (after translation) and keep reference
                self.dlg = SGToolDockWidget()

            # connect to provide cleanup on closing of dockwidget
            self.dlg.closingPlugin.connect(self.onClosePlugin)

            # show the dockwidget
            # TODO: fix to allow choice of dock location
            self.iface.addDockWidget(Qt.RightDockWidgetArea, self.dlg)
            self.dlg.show()

            self.deriv_dir_list = []
            self.deriv_dir_list.append("z")
            self.deriv_dir_list.append("x")
            self.deriv_dir_list.append("y")
            self.dlg.comboBox_derivDirection.addItems(self.deriv_dir_list)

            self.contin_dir_list = []
            self.contin_dir_list.append("up")
            self.contin_dir_list.append("down")
            self.dlg.comboBox_2_continuationDirection.addItems(self.contin_dir_list)

            self.ret_p_list = []
            self.ret_p_list.append("Pole")
            self.ret_p_list.append("Eqtr")
            self.dlg.comboBox_3_rte_p_list.addItems(self.ret_p_list)
            
            self.freq_cut_type_list = []
            self.freq_cut_type_list.append("Low")
            self.freq_cut_type_list.append("High")
            self.dlg.comboBox_2_FreqCutType.addItems(self.freq_cut_type_list)
            

            self.dlg.pushButton_4_calcIGRF.clicked.connect(self.update_mag_field)
            self.dlg.pushButton_2_selectGrid.clicked.connect(
                self.select_grid_file
            )
            self.dlg.pushButton_3_applyProcessing.clicked.connect(self.processGeophysics)

            self.dlg.pushButton_selectPoints.clicked.connect(
                self.select_point_file
            )
            self.dlg.pushButton_3_applyGridding.clicked.connect(self.gridPoints)

            self.dlg.mMapLayerComboBox_selectGrid.layerChanged.connect(
                self.update_paths
            )

            self.dlg.mQgsProjectionSelectionWidget.setCrs(
                QgsCoordinateReferenceSystem("EPSG:4326")
            )


