# -*- coding: utf-8 -*-
"""
/***************************************************************************
 GeodesicDensifier
                                 A QGIS plugin
 Adds vertices to geometry along geodesic lines
                              -------------------
        begin                : 2017-10-06
        git sha              : $Format:%H$ test
        copyright            : (C) 2017 by Jonah Sullivan
        email                : jonah.sullivan@ga.gov.au
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
try:
    # use system version of geographiclib
    from geographiclib.geodesic import Geodesic
except ImportError:
    # use version of geographiclib distributed with plugin
    import site
    import os
    # this will get the path for this file and add it to the system PATH
    # so the geographiclib folder can be found
    site.addsitedir(os.path.abspath(os.path.dirname(__file__)))
    from geographiclib.geodesic import Geodesic
import math
from qgis.core import QgsCoordinateReferenceSystem, QgsCoordinateTransform, \
    QGis, QgsFeature, QgsPoint, QgsGeometry, QgsField
from qgis.gui import QgsMessageBar
from PyQt4.QtCore import QSettings, QTranslator, qVersion, QCoreApplication, QVariant
from PyQt4.QtGui import QAction, QIcon
# Initialize Qt resources from file resources.py
import resources
# Import the code for the dialog
from geodesic_densifier_dialog import GeodesicDensifierDialog
import os.path


class GeodesicDensifier:
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

        # Declare instance attributes
        self.actions = []
        self.menu = u'&Geodesic Densifier'
        self.toolbar = self.iface.addToolBar(u'GeodesicDensifier')
        self.toolbar.setObjectName(u'GeodesicDensifier')

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

        # Create the dialog (after translation) and keep reference
        self.dlg = GeodesicDensifierDialog()

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

        icon_path = ':/plugins/GeodesicDensifier/icon.png'
        self.add_action(
            icon_path,
            text=u'Geodesic Densifier',
            callback=self.run,
            parent=self.iface.mainWindow())

    def unload(self):
        """Removes the plugin menu item and icon from QGIS GUI."""
        for action in self.actions:
            self.iface.removePluginMenu(
                u'&Geodesic Densifier', action)
            self.iface.removeToolBarIcon(action)
        # remove the toolbar
        del self.toolbar

    def run(self):
        """Run method that performs all the real work"""

        # show the dialog
        self.dlg.show()

        # set default values
        self.inLayer = self.dlg.mMapLayerComboBox.currentLayer()

        def set_in_layer():
            """ function to set the input layer from the GUI """
            self.inLayer = self.dlg.mMapLayerComboBox.currentLayer()
            if self.inLayer:
                if self.inLayer.crs():
                    self.dlg.messageBox.setText("Input Layer Set: " + str(self.inLayer.name()))
                else:
                    self.dlg.messageBox.setText("Error: Input must have projection defined")

        # listener to set input layer when combo box changes
        self.dlg.mMapLayerComboBox.layerChanged.connect(set_in_layer)

        # clear the ellipsoid combobox
        self.dlg.EllipsoidcomboBox.clear()

        # this is a dictionary of common ellipsoid parameters
        # http://www.ga.gov.au/__data/assets/file/0019/11377/Vincentys-formulae-to-calculate-distance-and-bearing-from-latitude-and-longitude.xls
        ellipsoid_dict = {'165': [6378165.000, 298.3],
                          'ANS': [6378160, 298.25],
                          'CLARKE 1858': [6378293.645, 294.26],
                          'GRS80': [6378137, 298.2572221],
                          'WGS84': [6378137, 298.2572236],
                          'WGS72': [6378135, 298.26],
                          'International 1924': [6378388, 297]}

        # add items to ellipsoid combobox
        for k in ellipsoid_dict.keys():
            self.dlg.EllipsoidcomboBox.addItem(str(k))

        # default ellipsoid is WGS84
        self.ellipsoid_a = 6378137.0
        self.ellipsoid_f = 298.2572236
        self.ellipsoid_name = 'WGS84'

        def set_in_ellipsoid():
            """ This function gets the ellipsoid name from the GUI and sets the parameters """
            in_ellipsoid_name = self.dlg.EllipsoidcomboBox.currentText()
            for k in ellipsoid_dict.keys():
                if k == in_ellipsoid_name:
                    self.ellipsoid_a = ellipsoid_dict[k][0]
                    self.ellipsoid_f = ellipsoid_dict[k][1]
                    self.ellipsoid_name = k
                    self.dlg.messageBox.setText("Ellipsoid set to " + str(k))

        # listener to set input ellipsoid when combo box changes
        self.dlg.EllipsoidcomboBox.currentIndexChanged.connect(set_in_ellipsoid)

        # default point spacing is 900
        self.spacing = 900

        def set_in_spacing():
            """ This function gets the spacing from the GUI and sets the parameter"""
            self.spacing = int(self.dlg.spacingSpinBox.value())
            self.dlg.messageBox.setText("Point spacing set to " + str(self.spacing) + "m")

        # listener to set input point spacing when spin box changes
        self.dlg.spacingSpinBox.valueChanged.connect(set_in_spacing)

        # Run the dialog event loop
        result = self.dlg.exec_()
        # See if OK was pressed
        if result:

            # set the input layer
            self.inLayer = self.dlg.mMapLayerComboBox.currentLayer()

            # get the field list
            fields = self.inLayer.pendingFields()

            # handle layers that aren't WGS84 (EPSG:4326)
            wgs84crs = QgsCoordinateReferenceSystem("EPSG:4326")
            if self.inLayer.crs() != wgs84crs:
                transtowgs84 = QgsCoordinateTransform(self.inLayer.crs(), wgs84crs)
                transfromwgs84 = QgsCoordinateTransform(wgs84crs, self.inLayer.crs())

            # get input geometry type
            self.inType = 'Unknown'
            if self.inLayer.wkbType() == QGis.WKBPoint:
                self.inType = 'Point'           # works

            elif self.inLayer.wkbType() == QGis.WKBMultiPoint:
                self.inType = 'MultiPoint'          # multipoint won't be implemented

            elif self.inLayer.wkbType() == QGis.WKBLineString:
                self.inType = 'LineString'          # works

            elif self.inLayer.wkbType() == QGis.WKBMultiLineString:
                self.inType = 'MultiLineString'     # not working

            elif self.inLayer.wkbType() == QGis.WKBPolygon:
                self.inType = 'Polygon'             # works

            elif self.inLayer.wkbType() == QGis.WKBMultiPolygon:
                self.inType = 'MultiPolygon'        # multipolygon detection isn't working

            # setup output layers
            if self.inType == 'Point':
                self.create_point = True
                # create and add to map canvas a point memory layer
                layer_name = "Densified Point " + str(self.ellipsoid_name) + " " + str(self.spacing) + "m"
                out_point_layer = self.iface.addVectorLayer("Point?crs={}".format(self.inLayer.crs().authid()),
                                                            layer_name,
                                                            "memory")
                # set data provider
                pointPr = out_point_layer.dataProvider()
                # add attribute fields
                pointPr.addAttributes(fields)
                # add field to store "original" or "densified" point attribute
                pointTypeField = ''
                for fieldName in ["pointType", "pntType", "pntTyp"]:
                if fieldName not in [field.name() for field in fields]:
                    pointTypeField = fieldName
                pointPr.addAttributes([QgsField(pointTypeField, QVariant.String)])
                out_point_layer.updateFields()
            else:
                self.create_point = False

            if self.inType == 'LineString':
                self.create_polyline = True
                # create and add to map canvas a polyline memory layer
                layer_name = "Densified Line " + str(self.ellipsoid_name) + " " + str(self.spacing) + "m"
                out_line_layer = self.iface.addVectorLayer("LineString?crs={}".format(self.inLayer.crs().authid()),
                                                           layer_name,
                                                           "memory")
                # set data provider
                linePr = out_line_layer.dataProvider()
                # add attribute fields
                linePr.addAttributes(fields)
                out_line_layer.updateFields()
            else:
                self.create_polyline = False

            if self.inType == 'Polygon':
                self.create_polygon = True
                # create and add to map canvas a polyline memory layer
                layer_name = "Densified Polygon " + str(self.ellipsoid_name) + " " + str(self.spacing) + "m"
                out_poly_layer = self.iface.addVectorLayer("Polygon?crs={}".format(self.inLayer.crs().authid()),
                                                           layer_name,
                                                           "memory")
                # set data provider
                polyPr = out_poly_layer.dataProvider()
                # add attribute fields
                polyPr.addAttributes(fields)
                out_poly_layer.updateFields()
            else:
                self.create_polygon = False

            # Create a geographiclib Geodesic object
            self.geod = Geodesic(self.ellipsoid_a, 1 / self.ellipsoid_f)

            def densifyPoint(inLayer, pr):
                """ This function densifies the input point layer and writes it to the output provider"""
                # iterator to read input layer
                iterator = inLayer.getFeatures()
                # counter to mark first point as "original"
                counter = 0
                # empty feature used to store temporary data
                currentFeature = QgsFeature()
                # counter to report features that don't work
                badGeom = 0
                for feature in iterator:
                    try:
                        if counter == 0:
                            # this is only for the first point
                            geom = feature.geometry().asPoint()
                            currentFeature.setGeometry(QgsGeometry.fromPoint(geom))
                            attr = feature.attributes()
                            attr.append("Original")
                            currentFeature.setAttributes(attr)
                            pr.addFeatures([currentFeature])
                        elif counter > 0:
                            startPt = currentFeature.geometry().asPoint()
                            endPt = feature.geometry().asPoint()
                            if self.inLayer.crs() != wgs84crs:
                                startPt = transtowgs84.transform(startPt)
                                endPt = transtowgs84.transform(endPt)
                            # create a geographiclib line object
                            lineObject = self.geod.InverseLine(startPt.y(), startPt.x(), endPt.y(), endPt.x())
                            # determine how many densified segments there will be
                            n = int(math.ceil(lineObject.s13 / self.spacing))
                            # adjust the spacing distance
                            seglen = lineObject.s13 / n
                            # create densified points along the line object
                            for i in range(1, n):
                                if i > 0:
                                    s = seglen * i
                                    g = lineObject.Position(s,
                                                            Geodesic.LATITUDE |
                                                            Geodesic.LONGITUDE |
                                                            Geodesic.LONG_UNROLL)
                                    geom = QgsPoint(g['lon2'], g['lat2'])
                                    attr = feature.attributes()
                                    attr.append("Densified")
                                    currentFeature.setAttributes(attr)
                                    if self.inLayer.crs() != wgs84crs:  # Convert each point back to the output CRS
                                        geom = transfromwgs84.transform(geom)
                                    currentFeature.setGeometry(QgsGeometry.fromPoint(geom))
                                    # write the point
                                    pr.addFeatures([currentFeature])
                            # write the last point
                            geom = feature.geometry().asPoint()
                            currentFeature.setGeometry(QgsGeometry.fromPoint(geom))
                            attr = feature.attributes()
                            attr.append("Original")
                            currentFeature.setAttributes(attr)
                            pr.addFeatures([currentFeature])
                        counter += 1
                    except:
                        badGeom += 1
                if badGeom > 0:
                    # report number of features that didn't work
                    self.iface.messageBar().pushMessage("", "{} features failed".format(badGeom),
                                                        level=QgsMessageBar.WARNING, duration=5)

            def densifyLine(inLayer, pr):
                # counter to report features that don't work
                badGeom = 0
                # iterator to read input layer
                iterator = inLayer.getFeatures()
                # create empty feature to write to
                newLine = QgsFeature()
                for feature in iterator:
                    try:
                        # convert lineString into list of lists like a multiLineString so they are both processed the same way
                        if self.inType == 'LineString':
                            segments = [feature.geometry().asPolyline()]
                        elif self.inType == 'MultiLineString':
                            segments = feature.geometry().asMultiPolyline()
                        else:
                            badGeom += 1
                    except:
                        badGeom += 1
                    if self.inType == 'LineString':
                        line = segments[0]
                        pointCount = len(line)
                        startPt = QgsPoint(line[0][0], line[0][1])
                        if self.inLayer.crs() != wgs84crs:
                            startPt = transtowgs84.transform(startPt)
                        pointList = [startPt]
                        for i in range(1, pointCount):
                            endPt = QgsPoint(line[i][0], line[i][1])
                            if self.inLayer.crs() != wgs84crs:
                                endPt = transtowgs84.transform(endPt)
                            # create a geographiclib line object
                            lineObject = self.geod.InverseLine(startPt.y(), startPt.x(), endPt.y(), endPt.x())
                            # determine how many densified segments there will be
                            n = int(math.ceil(lineObject.s13 / self.spacing))
                            if lineObject.s13 > self.spacing:
                                seglen = lineObject.s13 / n
                                for j in range(1, n):
                                    s = seglen * j
                                    g = lineObject.Position(s,
                                                            Geodesic.LATITUDE |
                                                            Geodesic.LONGITUDE |
                                                            Geodesic.LONG_UNROLL)
                                    pointList.append(QgsPoint(g['lon2'], g['lat2']))
                            pointList.append(endPt)
                            startPt = endPt
                        if self.inLayer.crs() != wgs84crs:  # Convert each point back to the output CRS
                            for x, pt in enumerate(pointList):
                                pointList[x] = transfromwgs84.transform(pt)
                        newLine.setGeometry(QgsGeometry.fromPolyline(pointList))
                    elif self.inType == 'MultiLineString':
                        outsegment = []
                        for line in segments:
                            pointCount = len(line)
                            startPt = QgsPoint(line[0][0], line[0][1])
                            if self.inLayer.crs() != wgs84crs:  # Convert to 4326
                                startPt = transtowgs84.transform(startPt)
                            pts = [startPt]
                            for x in range(1, pointCount):
                                endPt = QgsPoint(line[x][0], line[x][1])
                                if self.inLayer.crs() != wgs84crs:  # Convert to 4326
                                    endPt = transtowgs84.transform(endPt)
                                lineObject = self.geod.InverseLine(startPt.y(), startPt.x(), endPt.y(), endPt.x())
                                n = int(math.ceil(lineObject.s13 / self.spacing))
                                if lineObject.s13 > self.spacing:
                                    seglen = lineObject.s13 / n
                                    for i in range(1, n):
                                        s = seglen * i
                                        g = lineObject.Position(s,
                                                                Geodesic.LATITUDE |
                                                                Geodesic.LONGITUDE |
                                                                Geodesic.LONG_UNROLL)
                                        pts.append(QgsPoint(g['lon2'], g['lat2']))
                                pts.append(endPt)
                                startPt = endPt

                            if self.inLayer.crs() != wgs84crs:  # Convert each point back to the output CRS
                                for x, pt in enumerate(pts):
                                    pts[x] = transfromwgs84.transform(pt)
                            outsegment.append(pts)

                        newLine.setGeometry(QgsGeometry.fromMultiPolyline(outsegment))

                    newLine.setAttributes(feature.attributes())
                    pr.addFeatures([newLine])
                if badGeom > 0:
                    self.iface.messageBar().pushMessage("", "{} features failed".format(badGeom),
                                                        level=QgsMessageBar.WARNING, duration=5)

            def densifyPolygon(inLayer, pr):
                badGeom = 0
                iterator = inLayer.getFeatures()
                # create empty feature to write to
                for feature in iterator:
                    try:
                        if self.inType == 'Polygon':
                            polygon = feature.geometry().asPolygon()
                            for points in polygon:
                                pointCount = len(points)
                                startPt = QgsPoint(points[0][0], points[0][1])
                                if self.inLayer.crs() != wgs84crs:
                                    startPt = transtowgs84.transform(startPt)
                                polyPointList = [startPt]
                                for i in range(1, pointCount):
                                    endPt = QgsPoint(points[i][0], points[i][1])
                                    if self.inLayer.crs() != wgs84crs:  # Convert to 4326
                                        endPt = transtowgs84.transform(endPt)
                                    lineObject = self.geod.InverseLine(startPt.y(), startPt.x(), endPt.y(), endPt.x())
                                    n = int(math.ceil(lineObject.s13 / self.spacing))
                                    seglen = lineObject.s13 / n
                                    for j in range(1, n):
                                        s = seglen * j
                                        g = lineObject.Position(s,
                                                                Geodesic.LATITUDE |
                                                                Geodesic.LONGITUDE |
                                                                Geodesic.LONG_UNROLL)
                                        polyPointList.append(QgsPoint(g['lon2'], g['lat2']))
                                    polyPointList.append(endPt)
                                    startPt = endPt

                                if self.inLayer.crs() != wgs84crs:
                                    for x, pt in enumerate(polyPointList):
                                        polyPointList[x] = transfromwgs84.transform(pt)

                            if len(polyPointList) > 0:
                                outPolygon = QgsFeature()
                                outPolygon.setGeometry(QgsGeometry.fromPolygon([polyPointList]))
                                outPolygon.setAttributes(feature.attributes())
                                pr.addFeatures([outPolygon])

                        else:
                            print "multipoly"
                            multipolygon = feature.geometry().asMultiPolygon()
                            multiPointList = []
                            for polygon in multipolygon:
                                polyPointList = []
                                for points in polygon:
                                    print 'points:', points
                                    pointCount = len(points)
                                    startPt = QgsPoint(points[0][0], points[0][1])
                                    print 'startPt:', startPt
                                    if self.inLayer.crs() != wgs84crs:
                                        startPt = transtowgs84.transform(startPt)
                                    polyPointList = [startPt]
                                    for i in range(1, pointCount):
                                        endPt = QgsPoint(points[i][0], points[i][1])
                                        if self.inLayer.crs() != wgs84crs:  # Convert to 4326
                                            endPt = transtowgs84.transform(endPt)
                                        lineObject = self.geod.InverseLine(startPt.y(), startPt.x(), endPt.y(),
                                                                           endPt.x())
                                        n = int(math.ceil(lineObject.s13 / self.spacing))
                                        seglen = lineObject.s13 / n
                                        for j in range(1, n):
                                            s = seglen * j
                                            g = lineObject.Position(s,
                                                                    Geodesic.LATITUDE |
                                                                    Geodesic.LONGITUDE |
                                                                    Geodesic.LONG_UNROLL)
                                            polyPointList.append(QgsPoint(g['lon2'], g['lat2']))
                                        polyPointList.append(endPt)
                                        startPt = endPt

                                    if self.inLayer.crs() != wgs84crs:
                                        for x, pt in enumerate(polyPointList):
                                            polyPointList[x] = transfromwgs84.transform(pt)
                                    print 'polypointlist', polyPointList
                                multiPointList.append(polyPointList)

                            outMultiPolygon = QgsFeature()
                            outMultiPolygon.setGeometry(QgsGeometry.fromMultiPolygon(multiPointList))
                            outMultiPolygon.setAttributes(feature.attributes())
                            pr.addFeatures([outMultiPolygon])
                    except:
                        badGeom += 1
                if badGeom > 0:
                    self.iface.messageBar().pushMessage("", "{} features failed".format(badGeom),
                                                        level=QgsMessageBar.WARNING, duration=5)

            if self.create_point:
                densifyPoint(self.inLayer, pointPr)
                out_point_layer.updateExtents()

            elif self.create_polyline:
                densifyLine(self.inLayer, linePr)
                out_line_layer.updateExtents()

            elif self.create_polygon:
                densifyPolygon(self.inLayer, polyPr)
                out_poly_layer.updateExtents()

            else:
                self.iface.messageBar().pushMessage("", "geometry type not recognized",
                                                    level=QgsMessageBar.WARNING, duration=5)
