# -*- coding: utf-8 -*-
"""
/***************************************************************************
 raster_tracer
                                 A QGIS plugin
 Raster tracer
                              -------------------
        begin                : 2017-06-08
        git sha              : $Format:%H$
        copyright            : (C) 2017 by Mikhail Kondratyev
        email                : fatune@gmail.com
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
from PyQt4.QtCore import QSettings, QTranslator, qVersion, QCoreApplication, Qt, pyqtSignal
from PyQt4.QtGui import QAction, QIcon, QApplication

from qgis.core import *
from qgis.gui import *
# Initialize Qt resources from file resources.py
import resources

# Import the code for the DockWidget
from raster_tracer_dockwidget import raster_tracerDockWidget
import os.path

from osgeo import gdal 

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage.filters import gaussian_filter
from scipy import ndimage
import implementation
import timeit
from math import atan2, radians, tan
from matplotlib.colors import rgb_to_hsv

#from colorsys import rgb_to_hsv
import _astar
import gc

class raster_tracer:
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
            'raster_tracer_{}.qm'.format(locale))

        if os.path.exists(locale_path):
            self.translator = QTranslator()
            self.translator.load(locale_path)

            if qVersion() > '4.3.3':
                QCoreApplication.installTranslator(self.translator)

        # Declare instance attributes
        self.actions = []
        self.menu = self.tr(u'&raster_tracer')
        # TODO: We are going to let the user set this up in a future iteration
        self.toolbar = self.iface.addToolBar(u'raster_tracer')
        self.toolbar.setObjectName(u'raster_tracer')

        #print "** INITIALIZING raster_tracer"

        self.pluginIsActive = False
        self.dockwidget = None


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
        return QCoreApplication.translate('raster_tracer', message)


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

        icon_path = ':/plugins/raster_tracer/icon.png'
        self.add_action(
            icon_path,
            text=self.tr(u'Raster tracer'),
            callback=self.run,
            parent=self.iface.mainWindow())

    #--------------------------------------------------------------------------

    def onClosePlugin(self):
        """Cleanup necessary items here when plugin dockwidget is closed"""

        #print "** CLOSING raster_tracer"

        # disconnects
        self.dockwidget.closingPlugin.disconnect(self.onClosePlugin)

        # remove this statement if dockwidget is to remain
        # for reuse if plugin is reopened
        # Commented next statement since it causes QGIS crashe
        # when closing the docked window:
        # self.dockwidget = None

        self.pluginIsActive = False


    def unload(self):
        """Removes the plugin menu item and icon from QGIS GUI."""

        #print "** UNLOAD raster_tracer"

        for action in self.actions:
            self.iface.removePluginMenu(
                self.tr(u'&raster_tracer'),
                action)
            self.iface.removeToolBarIcon(action)
        # remove the toolbar
        del self.toolbar

    #--------------------------------------------------------------------------
    def raster_selected(self):
        layer_index = self.dockwidget.comboBox.currentIndex()
        self.map_canvas = self.iface.mapCanvas()
        self.tool_identify2 = PointTool(self.map_canvas, self.layers_list[layer_index],self.iface, self.dockwidget)
        self.map_canvas.setMapTool(self.tool_identify2)
        print "ya!"

    def line_type_selected(self):
        indx = self.dockwidget.comboBoxLinePoly.currentIndex()
        self.linetype = self.line_types[indx]
        self.tool_identify2.linetype = self.linetype
        print "line type selected: %s" % self.linetype

    def run(self):
        """Run method that loads and starts the plugin"""

        if not self.pluginIsActive:
            self.pluginIsActive = True

            #print "** STARTING raster_tracer"

            # dockwidget may not exist if:
            #    first run of plugin
            #    removed on close (see self.onClosePlugin method)
            if self.dockwidget == None:
                # Create the dockwidget (after translation) and keep reference
                self.dockwidget = raster_tracerDockWidget()

            # connect to provide cleanup on closing of dockwidget
            self.dockwidget.closingPlugin.connect(self.onClosePlugin)

            self.dockwidget.comboBox.currentIndexChanged.connect(self.raster_selected)
            self.dockwidget.comboBoxLinePoly.currentIndexChanged.connect(self.line_type_selected)

            # show the dockwidget
            # TODO: fix to allow choice of dock location
            self.iface.addDockWidget(Qt.TopDockWidgetArea, self.dockwidget)
            self.dockwidget.show()

        layers = self.iface.legendInterface().layers()
        self.layer_list = []
        self.layers_list = []
        for layer in layers:
             if not layer.type() == QgsMapLayer.RasterLayer: continue
             self.layer_list.append(layer.name())
             self.layers_list.append(layer)
        self.dockwidget.comboBox.addItems(self.layer_list)
        self.line_types = ["Line","Polygon"]
        self.dockwidget.comboBoxLinePoly.addItems(self.line_types)

        #layer_index = self.dockwidget.comboBox.currentIndex()
        #self.map_canvas = self.iface.mapCanvas()
        #self.tool_identify2 = PointTool(self.map_canvas, self.layers_list[layer_index],self.iface)
        #self.map_canvas.setMapTool(self.tool_identify2)

class PointTool(QgsMapToolEmitPoint):
    canvasDoubleClicked = pyqtSignal(object, object)

    def __init__(self, canvas,layer,iface,dockwidget):
        QApplication.restoreOverrideCursor()
        QApplication.setOverrideCursor(Qt.CrossCursor)
        QgsMapToolEmitPoint.__init__(self, canvas)
        self.layer = layer
        self.start = None
        self.end = None
        self.ends = []
        sample, geo_ref = get_whole_raster(self.layer)
        self.sample = sample
        r = self.sample[0]
        g = self.sample[1]
        b = self.sample[2]
        h,s,v = rgb2hsv(r,g,b)
        self.sample_hsv = h,s,v
        self.geo_ref = geo_ref
        self.iface=iface
        self.layer_v = None
        self.points = []
        self.double = False
        self.dockwidget = dockwidget
        self.linetype = "Line"

    def canvasDoubleClickEvent(self, event):
        self.double = True
    

    def canvasReleaseEvent(self, mouseEvent):
        if self.double:
            self.double = False
            self.layer_v = None
            self.points = []
            self.start = None
            self.end = None
            return

        if mouseEvent.button() == Qt.RightButton: # roll back
            if len(self.points)<2: return
            delete_segment(self.layer_v)
            self.iface.mapCanvas().refresh()
            self.points.pop()
            return


        qgsPoint = self.toMapCoordinates(mouseEvent.pos())
        x,y = qgsPoint.x(), qgsPoint.y()
        i,j = get_indxs_from_raster_coords(self.geo_ref, x,y)
        if self.dockwidget.checkBoxSnap.isChecked():
            i,j = snap_to_neares_pixel(i,j,self.sample_hsv)
            x,y = get_coords_from_raster_indxs(self.geo_ref, i,j)

        if not self.dockwidget.checkBoxTrace.isChecked():
            self.points.append((x,y))
            path = [self.points[-2],self.points[-1]]
            update_segment(path, self.layer_v, self.linetype)
            self.iface.mapCanvas().refresh()
            QApplication.restoreOverrideCursor()
            QApplication.setOverrideCursor(Qt.CrossCursor)
            return

        cvalue = self.layer.dataProvider().identify(QgsPoint(x, y), QgsRaster.IdentifyFormatValue).results()
        self.cval = [cvalue[1], cvalue[2],cvalue[3]]
        if self.cval[0] == None or self.cval[1] == None or self.cval[2] == None: return

        self.points.append((x,y))

        if len(self.points)<2:
            source_crs = self.layer.crs()
            self.layer_v = add_segment(source_crs, linetype = self.linetype)
            QApplication.restoreOverrideCursor()
            QApplication.setOverrideCursor(Qt.CrossCursor)
            return

        QApplication.setOverrideCursor(Qt.WaitCursor)

        self.start = self.points[-2]
        self.end = self.points[-1]

        if not (self.end and self.start): return
        cx = (self.start[0] + x)*0.5
        cy = (self.start[1] + y)*0.5

        time1 = timeit.default_timer()
        print self.cval
        raster = prepare_raster(self.sample, self.cval)


        time3 = timeit.default_timer()
        start_ind = get_indxs_from_raster_coords(self.geo_ref, self.start[0],self.start[1])
        end_ind = get_indxs_from_raster_coords(self.geo_ref, self.end[0],self.end[1])


        ###k = np.median(raster)
        #size=100
        #im = raster[i-size:i+size,j-size:j+size]

        #r = self.sample[0][i-size:i+size,j-size:j+size]
        #g = self.sample[1][i-size:i+size,j-size:j+size]
        #b = self.sample[2][i-size:i+size,j-size:j+size]
        #h,s,v = rgb_to_hsv(r,g,b)
        #h_, s_, v_ = rgb_to_hsv_single(self.cval[0],self.cval[1],self.cval[2])
        #im2 = (h-h_)**2 + (s-s_)**2 + ((v-v_)/255)**2
        #im3 = (v_-v)**2
        #im4 = (h_-h)**2
        #print (h.max(), s.max(), v.max())

        #sx = ndimage.sobel(im, axis=0, mode='constant')
        #sy = ndimage.sobel(im, axis=1, mode='constant')
        #sob = np.hypot(sx, sy)
        ##h,s,v = self.sample_hsv
        ##h = h [i-size:i+size,j-size:j+size]
        ##s = s [i-size:i+size,j-size:j+size]
        ##v = v [i-size:i+size,j-size:j+size]
        #plt.subplot(231)
        #plt.imshow(im,interpolation='None',cmap='gray_r')
        ##plt.scatter(size+1,size+1)
        #plt.subplot(232)
        ##k = np.percentile(raster,25)
        ##im2 = np.zeros_like(im)
        ##im2[im>k] = 1
        ###plt.imshow(im2,interpolation='None',cmap='gray')
        #plt.imshow(im2,interpolation='None',cmap='gray_r')

        #plt.subplot(233)
        #plt.imshow(sob,interpolation='None',cmap='gray_r')

        #plt.subplot(234)
        ##k = np.percentile(raster,50)
        ##im2 = np.zeros_like(im)
        ##im2[im>k] = 1
        ###plt.imshow(im2,interpolation='None',cmap='gray')
        #plt.imshow(h,interpolation='None',cmap='gray_r')
        #plt.subplot(235)
        ##k = np.percentile(raster,75)
        ##im2 = np.zeros_like(im)
        ##im2[im>k] = 1
        ###plt.imshow(im2,interpolation='None',cmap='gray')
        ##plt.imshow(s,interpolation='None',cmap='gray_r')
        ##plt.subplot(236)
        ##plt.imshow(v,interpolation='None',cmap='gray_r')
        ##plt.show()

        ##print raster.min(), np.percentile(raster,25), np.percentile(raster,50), np.percentile(raster,75), raster.max()
        ##k = np.percentile(raster,75)
        ##raster2 = np.zeros_like(raster)
        ##raster2[raster>k] = 1
        ##plt.imshow(raster2, interpolation='None', cmap='gray')
        ##plt.show()
        ##QApplication.restoreOverrideCursor()
        ##return



        print "h"
        cost_so_far = np.zeros_like(raster).astype(np.int)
        j_range, i_range = raster.shape
        try:
            cost_so_far2 = _astar.astar(start_ind[0],start_ind[1],
                                       end_ind[0],end_ind[1],
                                       raster.flatten(),
                                       cost_so_far.flatten(), i_range, j_range)
        except RuntimeError:
            if len(self.points)<2: return
            self.points.pop()
            QApplication.restoreOverrideCursor()
            QApplication.setOverrideCursor(Qt.CrossCursor)
            return


        cost_so_far = cost_so_far2.reshape(raster.shape)
        #plt.imshow(cost_so_far)
        #plt.scatter(start_ind[1], start_ind[0])
        #plt.scatter(end_ind[1], end_ind[0])
        #plt.show()
        print "h2"
        time4 = timeit.default_timer()
        path = restore_path(cost_so_far,end_ind, start_ind, self.geo_ref)
        if not path:
            if len(self.points)<2: return
            self.points.pop()
            QApplication.restoreOverrideCursor()
            QApplication.setOverrideCursor(Qt.CrossCursor)
            return
        path.insert(0,start_ind)
        path = smooth(path,size=5)
        path = simplify(path)

        time5 = timeit.default_timer()
        path2 = [(get_coords_from_raster_indxs(self.geo_ref, i,j)) for i,j in path]

        update_segment(path2, self.layer_v, linetype = self.linetype)
        self.iface.mapCanvas().refresh()
        QApplication.restoreOverrideCursor()
        QApplication.setOverrideCursor(Qt.CrossCursor)
        del raster
        gc.collect()


def restore_path(grid, start, end,geo_ref):
    increments = [(1,0),(1,1),(0,1),(0,-1),(-1,-1),(-1,0),(-1,1),(1,-1)]
    path = []
    current = start
    time_start = timeit.default_timer()
    while True:
        path.append(current)
        neighbours = [(current[0] + dx, current[1]+dy) for dx,dy in increments]
        vals = []
        for x,y in neighbours:
            try:
                val = grid[x,y]
                if val == 0 and not ((x,y) == end): val = 10000000000000
                vals.append(val)
            except IndexError:
                vals.append(10000000000000000000)
        minimum = min(vals)
        current = neighbours[vals.index(minimum)]
        if current == end: break
        #if grid[current] <2: break
        if  timeit.default_timer()-time_start > 0.1: return None
    return path[::-1]



def draw_segment(segment,source_crs):
    s = QSettings()
    # possible values are: prompt, useProject, useGlobal
    default_value = s.value("/Projections/defaultBehaviour")
    s.setValue("/Projections/defaultBehaviour", "useProject")
    
    # here, add your layer to the mapcanvas
    
    layer = QgsVectorLayer('LineString', 'connecting line', "memory") #this will create an in memory layer which is not stored on a HDD
    s.setValue("/Projections/defaultBehaviour", default_value)
    layer.setCrs(source_crs) #we set the CRS to the one provided by the point layer
    pr = layer.dataProvider() #once again we need the provider (should be shapefile/ESRI)
    line = QgsFeature() #let's create an empty Feature which will be the box to put our segment(s) into
    seg = []
    #print segment
    for (x,y) in segment:
        seg.append(QgsPoint(x,y))
    #seg = [QgsPoint(x,y),QgsPoint(x+100,y+100)] #so this is our whole segment
    line.setGeometry(QgsGeometry.fromPolyline(seg)) #we will set the geometry of our so far empty feature from this array of points
    pr.addFeatures([line]) #and add it to the layer
    layer.updateExtents() #update it
    QgsMapLayerRegistry.instance().addMapLayer(layer) #and put it on the map. Store it manually if you need to.

def add_segment(source_crs, linetype = "Line"):
    s = QSettings()
    # possible values are: prompt, useProject, useGlobal
    default_value = s.value("/Projections/defaultBehaviour")
    s.setValue("/Projections/defaultBehaviour", "useProject")
    
    # here, add your layer to the mapcanvas
    
    if linetype == "Line":
        layer = QgsVectorLayer('LineString', 'connecting line', "memory") #this will create an in memory layer which is not stored on a HDD
    elif linetype == "Polygon":
        layer = QgsVectorLayer('Polygon', 'connecting line', "memory") #this will create an in memory layer which is not stored on a HDD
    s.setValue("/Projections/defaultBehaviour", default_value)
    layer.setCrs(source_crs) #we set the CRS to the one provided by the point layer
    return layer

def update_segment (segment, layer, linetype='Line'):
    pr = layer.dataProvider() #once again we need the provider (should be shapefile/ESRI)
    line = QgsFeature() #let's create an empty Feature which will be the box to put our segment(s) into
    seg = []
    #print segment
    for (x,y) in segment:
        seg.append(QgsPoint(x,y))
    #seg = [QgsPoint(x,y),QgsPoint(x+100,y+100)] #so this is our whole segment
    layer.startEditing()
    if linetype == "Line":
        line.setGeometry(QgsGeometry.fromPolyline(seg)) #we will set the geometry of our so far empty feature from this array of points
    elif linetype == "Polygon":
        line.setGeometry(QgsGeometry.fromPolygon([seg])) #we will set the geometry of our so far empty feature from this array of points
    pr.addFeatures([line]) #and add it to the layer
    layer.commitChanges()
    layer.updateExtents() #update it
    
    symbol_layer = QgsSimpleLineSymbolLayerV2.create({'color':'red'})
    layer.rendererV2().symbols()[0].changeSymbolLayer(0, symbol_layer)
    QgsMapLayerRegistry.instance().addMapLayer(layer) #and put it on the map. Store it manually if you need to.

def delete_segment(layer):
    ids = [f.id() for f in layer.getFeatures()]
    layer.startEditing()
    layer.deleteFeature(ids[-1])
    layer.commitChanges()

def get_whole_raster(layer):
    raster_path = layer.source()
    ds = gdal.Open(raster_path) 
    band1 = np.array(ds.GetRasterBand(1).ReadAsArray())
    band2 = np.array(ds.GetRasterBand(2).ReadAsArray())
    band3 = np.array(ds.GetRasterBand(3).ReadAsArray())
    geo_ref  = ds.GetGeoTransform()
    return (band1,band2,band3), geo_ref

def get_indxs_from_raster_coords(geo_ref, x, y):
    top_left_x, we_resolution, _, top_left_y, _, ns_resolution = geo_ref
    
    i = int( (y - top_left_y) /ns_resolution )
    j = int( (x - top_left_x) /we_resolution )
    
    return i,j

def get_coords_from_raster_indxs(geo_ref, i, j):
    top_left_x, we_resolution, _, top_left_y, _, ns_resolution = geo_ref

    y = top_left_y- i*we_resolution
    x = top_left_x - j*ns_resolution

    return x,y

def prepare_raster(sample,value):
    return (sample[0] - value[0])**2 + (sample[1]-value[1])**2 + (sample[2]-value[2])**2
    return np.sqrt((sample[0] - value[0])**2 + (sample[1]-value[1])**2 + (sample[2]-value[2])**2)

def prepare_raster_hsv((h,s,v),value):
    return np.sqrt(((h-value[0])/360)**2 + ((s-value[1]))**2 + ((v-value[2])*10)**2)

def prepare_grey_raster(sample):
    ban0 = ndimage.gaussian_filter(sample[0], sigma=1)
    ban1 = ndimage.gaussian_filter(sample[1], sigma=1)
    ban2 = ndimage.gaussian_filter(sample[2], sigma=1)
    return (ban0+ban1+ban2)/3.0
    return (sample[0] + sample[1] + sample[2]) / 3
    return np.median(sample)

def smooth(path, size=2):
    smoothed = []
    smoothed.append(path[0])
    for i in range(size, len(path)-size):
        xx = [x for (x,y) in path[i-size:i+size]]
        yy = [y for (x,y) in path[i-size:i+size]]
        x = sum(xx)/(len(xx)*1.0)
        y = sum(yy)/(len(yy)*1.0)
        smoothed.append((x,y))
    smoothed.append(path[-1])
    return smoothed

def simplify(path, tolerance = 2):
    previous = None
    previousfactor = None
    nodes_to_delete = []
    tolerance = tan(radians(tolerance))
    for i, node in enumerate(path):
        if not previous: 
            previous = node
            continue
        factor = atan2((node[0]-previous[0]),(node[1]-previous[1]))
        #factor = ((node[0]-previous[0]),(node[1]-previous[1]))
        if not previousfactor:
            previousfactor = factor
            continue
        if abs(factor-previousfactor) < tolerance: nodes_to_delete.append(i-1)
        #print factor, previousfactor, abs(factor-previousfactor), tolerance
        previous = node
        previousfactor = factor

    for i in nodes_to_delete[::-1]:
        path.pop(i)
    return path

def rgb2hsv(r, g, b):
    r, g, b = r/255.0, g/255.0, b/255.0
    mx = np.maximum(np.maximum(r, g), b)
    mn = np.minimum(np.minimum(r, g), b)
    df = mx-mn
    h = np.zeros_like(r)
    s = np.zeros_like(r)
    v = np.zeros_like(r)
    h[mx==r] = (60 * ((g[mx==r]-b[mx==r])/df[mx==r]) + 360) % 360
    h[mx==g] = (60 * ((b[mx==g]-r[mx==g])/df[mx==g]) + 120) % 360
    h[mx==b] = (60 * ((r[mx==b]-g[mx==b])/df[mx==b]) + 240) % 360
    h[mx==mn] = 0
    s = df/mx
    s[mx==0] = 0
    v = mx

    return h, s, v

def rgb_to_hsv(r, g, b):
    maxc = np.maximum(np.maximum(r, g), b)
    minc = np.minimum(np.minimum(r, g), b)
    v = maxc

    h = np.zeros_like(r)
    s = np.zeros_like(r)

    h[minc==maxc] == 0.0
    s[minc==maxc] == 0.0

    s = (maxc-minc) / maxc
    rc = (maxc-r) / (maxc-minc)
    gc = (maxc-g) / (maxc-minc)
    bc = (maxc-b) / (maxc-minc)

    h          = 4.0 + gc          - rc
    h[r==maxc] =       bc[r==maxc] - gc[r==maxc]
    h[g==maxc] = 2.0 + rc[g==maxc] - bc[g==maxc]
    h = (h/6.0) % 1.0

    return h, s, v

def rgb_to_hsv_single(r, g, b):
    maxc = max(r, g, b)
    minc = min(r, g, b)
    v = maxc
    if minc == maxc:
        return 0.0, 0.0, v
    s = (maxc-minc) / maxc
    rc = (maxc-r) / (maxc-minc)
    gc = (maxc-g) / (maxc-minc)
    bc = (maxc-b) / (maxc-minc)
    if r == maxc:
        h = bc-gc
    elif g == maxc:
        h = 2.0+rc-bc
    else:
        h = 4.0+gc-rc
    h = (h/6.0) % 1.0
    return h, s, v

def rgb2hsv_single(r,g,b):
    r, g, b = r/255.0, g/255.0, b/255.0
    mx = max(r, g, b)
    mn = min(r, g, b)
    df = mx-mn
    if mx == mn:
        h = 0
    elif mx == r:
        h = (60 * ((g-b)/df) + 360) % 360
    elif mx == g:
        h = (60 * ((b-r)/df) + 120) % 360
    elif mx == b:
        h = (60 * ((r-g)/df) + 240) % 360
    if mx == 0:
        s = 0
    else:
        s = df/mx
    v = mx
    return h, s, v

def snap_to_neares_pixel(i,j,(h,s,v)):
    size = 100
    im = v[i-size:i+size,j-size:j+size]
    mask = np.zeros_like(im)
    m = np.median(im)
    treshold = m # if background is light
    if m > 0.5: 
        treshold = 1 - (1-m)*5 # if background is light
        mask[im<treshold] = 1
    else: # if backround is dark
        treshold = m*5
        mask[im>treshold] = 1
    mask = ndimage.binary_erosion(mask)

    distance = np.ones(im.shape)*2
    x_ = np.linspace(-1.0,1.0,size*2+1)
    y_ = np.linspace(-1.0,1.0,size*2+1)
    xx,yy = np.meshgrid(x_,y_)
    distance_ = xx**2 + yy**2
    distance[mask==1] = distance_[mask==1]
    indexes = np.where(distance==distance.min())
    i_, j_ = indexes[0][0],indexes[1][0]
    i_, j_ = i+i_-size, j+j_-size
    return i_, j_
