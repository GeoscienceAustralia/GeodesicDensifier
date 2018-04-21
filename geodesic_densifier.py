# -*- coding: utf-8 -*-
"""
/***************************************************************************
 GeodesicDensifier
                                 A QGIS plugin
 Adds vertices to geometry along geodesic lines
                              -------------------
        begin                : 2018-02-21
        git sha              : $Format:%H$
        copyright            : (C) 2018 by Jonah Sullivan
        email                : jonah.sullivan@ga.gov.au
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the Creative Commons 4.0 International License. *
 *                                                                         *
 ***************************************************************************/
"""
try:
    # use system version of geographiclib
    from geographiclib.geodesic import Geodesic
except ImportError:
    # use version of geographiclib distributed with plugin
    import site
    from os.path import abspath, dirname
    from inspect import getsourcefile
    # this will get the path for this file and add it to the system PATH
    # so the geographiclib folder can be found
    site.addsitedir(dirname(abspath(getsourcefile(lambda: 0))))
    from geographiclib.geodesic import Geodesic
import math
from qgis.core import (QgsCoordinateReferenceSystem,
                       QgsCoordinateTransform,
                       QgsWkbTypes,
                       QgsFeature,
                       QgsPointXY,
                       QgsGeometry,
                       QgsField,
                       QgsProject,
                       QgsMapLayerProxyModel)
from PyQt5.QtCore import (QSettings,
                          QTranslator,
                          qVersion,
                          QCoreApplication,
                          QVariant)
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import QAction

# Initialize Qt resources from file resources.py
from .resources import *
# Import the code for the dialog
from .geodesic_densifier_dialog import GeodesicDensifierDialog
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

        # Create the dialog (after translation) and keep reference
        self.dlg = GeodesicDensifierDialog()
        self.dlg.mMapLayerComboBox.setFilters(QgsMapLayerProxyModel.LineLayer |
                                              QgsMapLayerProxyModel.PolygonLayer |
                                              QgsMapLayerProxyModel.PointLayer)

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
            self.iface.removePluginMenu(u'&Geodesic Densifier', action)
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

        ellipsoid_dict = {'165': [6378165.000, 298.3],
                          'ANS': [6378160, 298.25],
                          'CLARKE 1858': [6378293.645, 294.26],
                          'GRS80': [6378137, 298.2572221],
                          'WGS72': [6378135, 298.26],
                          'International 1924': [6378388, 297],
                          'WGS84': [6378137, 298.2572236]}

        # add items to ellipsoid combobox
        for k in list(ellipsoid_dict.keys()):
            self.dlg.EllipsoidcomboBox.addItem(str(k))

        # default ellipsoid is WGS84
        self.ellipsoid_a = 6378137.0
        self.ellipsoid_f = 298.2572236
        self.ellipsoid_name = 'WGS84'
        self.dlg.EllipsoidcomboBox.setCurrentText(self.ellipsoid_name)

        def set_in_ellipsoid():
            in_ellipsoid_name = self.dlg.EllipsoidcomboBox.currentText()
            for k in list(ellipsoid_dict.keys()):
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
            self.spacing = int(self.dlg.spacingSpinBox.value())
            self.dlg.messageBox.setText("Point spacing set to " + str(self.spacing) + "m")

        # listener to set input point spacing when spin box changes
        self.dlg.spacingSpinBox.valueChanged.connect(set_in_spacing)

        # Run the dialog event loop
        result = self.dlg.exec_()
        # See if OK was pressed
        if result:

            # check input parameters are valid
            self.dlg.mMapLayerComboBox.currentLayer() != None

            # set the input layer
            self.inLayer = self.dlg.mMapLayerComboBox.currentLayer()

            # get the field list
            fields = self.inLayer.fields()

            # handle layers that aren't WGS84 (EPSG:4326)
            wgs84crs = QgsCoordinateReferenceSystem("EPSG:4326")
            if self.inLayer.crs() != wgs84crs:
                transtowgs84 = QgsCoordinateTransform(self.inLayer.crs(), wgs84crs, QgsProject.instance())
                transfromwgs84 = QgsCoordinateTransform(wgs84crs, self.inLayer.crs(), QgsProject.instance())

            # get input geometry type
            if self.inLayer.geometryType() == QgsWkbTypes.PointGeometry:
                self.inType = 'Point'           # works

            if self.inLayer.geometryType() == QgsWkbTypes.LineGeometry:
                self.inType = 'LineString'      # works

            if self.inLayer.geometryType() == QgsWkbTypes.PolygonGeometry:
                self.inType = 'Polygon'         # works

            # setup output layers
            if self.inType == 'Point':
                self.create_point = True
                # create and add to map canvas a point memory layer
                layer_name = "Densified Point " + str(self.ellipsoid_name) + " " + str(self.spacing) + "m"
                out_point_layer = self.iface.addVectorLayer("Point?crs={}".format(self.inLayer.crs().authid()),
                                                            layer_name,
                                                            "memory")
                # set data provider
                provider = out_point_layer.dataProvider()
                # add attribute fields
                provider.addAttributes(fields)
                pointTypeField = ''
                if "pointType" not in [field.name() for field in fields]:
                    pointTypeField = "pointType"
                elif "pntType" not in [field.name() for field in fields]:
                    pointTypeField = "pntType"
                elif "pntTyp" not in [field.name() for field in fields]:
                    pointTypeField = "pntTyp"
                provider.addAttributes([QgsField(pointTypeField, QVariant.String)])
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
                provider = out_line_layer.dataProvider()
                # add attribute fields
                provider.addAttributes(fields)
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
                provider = out_poly_layer.dataProvider()
                # add attribute fields
                provider.addAttributes(fields)
                out_poly_layer.updateFields()
            else:
                self.create_polygon = False

            # Create a geographiclib Geodesic object
            self.geod = Geodesic(self.ellipsoid_a, 1 / self.ellipsoid_f)

            def densifyPoint(inLayer, pr):
                pointTypeFieldIdx = pr.fieldNameIndex(pointTypeField)
                iterator = inLayer.getFeatures()
                featureCount = pr.featureCount()
                counter = 0
                currentFeature = QgsFeature()
                badGeom = 0
                for feature in iterator:
                    if feature.geometry().wkbType() == QgsWkbTypes.Point:
                        try:
                            if counter == 0:
                                pointxy = feature.geometry().asPoint()
                                currentFeature.setGeometry(QgsGeometry.fromPointXY(pointxy))
                                attr = feature.attributes()
                                attr.append("Original")
                                currentFeature.setAttributes(attr)
                                pr.addFeatures([currentFeature])
                            else:
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
                                for i in range(1, n):
                                    if i > 0:
                                        s = seglen * i
                                        g = lineObject.Position(s, Geodesic.LATITUDE | Geodesic.LONGITUDE | Geodesic.LONG_UNROLL)
                                        geom = QgsPointXY(g['lon2'], g['lat2'])
                                        attr = feature.attributes()
                                        attr.append("Densified")
                                        currentFeature.setAttributes(attr)
                                        if self.inLayer.crs() != wgs84crs:  # Convert each point back to the output CRS
                                            geom = transfromwgs84.transform(geom)
                                        currentFeature.setGeometry(QgsGeometry.fromPointXY(geom))
                                        pr.addFeatures([currentFeature])
                                geom = feature.geometry().asPoint()
                                currentFeature.setGeometry(QgsGeometry.fromPointXY(geom))
                                attr = feature.attributes()
                                attr.append("Original")
                                currentFeature.setAttributes(attr)
                                pr.addFeatures([currentFeature])
                            counter += 1
                        except:
                            badGeom += 1
                            counter += 1
                    else:
                        badGeom += 1
                        self.iface.messageBar().pushWarning("multipoint geometries will not be densified")
                if badGeom > 0:
                    self.iface.messageBar().pushWarning("", "{} features failed".format(badGeom))

            def densifyPoly(inLayer, pr):
                badGeom = 0
                iterator = inLayer.getFeatures()
                for feature in iterator:
                    try:
                        if feature.geometry().wkbType() ==  QgsWkbTypes.LineString:
                            line_geom = feature.geometry().asPolyline()
                            geomType = "LineString"
                        elif feature.geometry().wkbType() ==  QgsWkbTypes.MultiLineString:
                            multiline_geom = feature.geometry().asMultiPolyline()
                            geomType = "MultiLineString"
                        elif feature.geometry().wkbType() ==  QgsWkbTypes.Polygon:
                            poly_geom = feature.geometry().asPolygon()
                            geomType = "Polygon"
                        elif feature.geometry().wkbType() == QgsWkbTypes.MultiPolygon:
                            multipoly_geom = feature.geometry().asMultiPolygon()
                            geomType = "MultiPolygon"
                        else:
                            badGeom += 1
                    except:
                        badGeom += 1

                    if geomType == "LineString":
                        dense_points = []
                        pointCount = len(line_geom)
                        startPt = QgsPointXY(line_geom[0][0], line_geom[0][1])
                        dense_points.append(startPt)
                        if self.inLayer.crs() != wgs84crs:
                            startPt = transtowgs84.transform(startPt)
                        for j in range(1, pointCount):
                            endPt = QgsPointXY(line_geom[j][0], line_geom[j][1])
                            if self.inLayer.crs() != wgs84crs:
                                endPt = transtowgs84.transform(endPt)
                            # create a geographiclib line object
                            lineObject = self.geod.InverseLine(startPt.y(), startPt.x(), endPt.y(), endPt.x())
                            # determine how many densified segments there will be
                            n = int(math.ceil(lineObject.s13 / self.spacing))
                            if lineObject.s13 > self.spacing:
                                seglen = lineObject.s13 / n
                                for k in range(1, n):
                                    s = seglen * k
                                    g = lineObject.Position(s, Geodesic.LATITUDE | Geodesic.LONGITUDE | Geodesic.LONG_UNROLL)
                                    waypoint = QgsPointXY(g['lon2'], g['lat2'])
                                    if self.inLayer.crs() != wgs84crs:
                                        waypoint = transfromwgs84.transform(waypoint)
                                    dense_points.append(waypoint)
                                if self.inLayer.crs() != wgs84crs:
                                    endPt = transfromwgs84.transform(endPt)
                                dense_points.append(endPt)
                            startPt = endPt

                    elif geomType == "MultiLineString":
                        dense_features = []
                        for i in range(len(multiline_geom)):
                            dense_points = []
                            line = multiline_geom[i]
                            pointCount = len(line)
                            startPt = QgsPointXY(line[0][0], line[0][1])
                            dense_points.append(startPt)
                            for j in range(1, pointCount):
                                endPt = QgsPointXY(line[j][0], line[j][1])
                                if self.inLayer.crs() != wgs84crs:
                                    startPt = transtowgs84.transform(startPt)
                                    endPt = transtowgs84.transform(endPt)
                                # create a geographiclib line object
                                lineObject = self.geod.InverseLine(startPt.y(), startPt.x(), endPt.y(), endPt.x())
                                # determine how many densified segments there will be
                                n = int(math.ceil(lineObject.s13 / self.spacing))
                                if lineObject.s13 > self.spacing:
                                    seglen = lineObject.s13 / n
                                    for k in range(1, n):
                                        s = seglen * k
                                        g = lineObject.Position(s, Geodesic.LATITUDE | Geodesic.LONGITUDE | Geodesic.LONG_UNROLL)
                                        waypoint = QgsPointXY(g['lon2'], g['lat2'])
                                        if self.inLayer.crs() != wgs84crs:
                                            waypoint = transfromwgs84.transform(waypoint)
                                        dense_points.append(waypoint)
                                    if self.inLayer.crs() != wgs84crs:
                                        endPt = transfromwgs84.transform(endPt)
                                    dense_points.append(endPt)
                                startPt = endPt
                            dense_features.append(dense_points)

                    elif geomType == "Polygon":
                        for poly in poly_geom:
                            dense_points = []
                            pointCount = len(poly)
                            startPt = QgsPointXY(poly[0][0], poly[0][1])
                            dense_points.append(startPt)
                            for j in range(1,pointCount):
                                endPt = QgsPointXY(poly[j][0], poly[j][1])
                                if self.inLayer.crs() != wgs84crs:
                                    endPt = transtowgs84.transform(endPt)
                                    startPt= transtowgs84.transform(startPt)
                                # create a geographiclib line object
                                lineObject = self.geod.InverseLine(startPt.y(), startPt.x(), endPt.y(), endPt.x())
                                # determine how many densified segments there will be
                                n = int(math.ceil(lineObject.s13 / self.spacing))
                                if lineObject.s13 > self.spacing:
                                    seglen = lineObject.s13 / n
                                    for k in range(1, n):
                                        s = seglen * k
                                        g = lineObject.Position(s, Geodesic.LATITUDE | Geodesic.LONGITUDE | Geodesic.LONG_UNROLL)
                                        waypoint = QgsPointXY(g['lon2'], g['lat2'])
                                        if self.inLayer.crs() != wgs84crs:
                                            waypoint = transfromwgs84.transform(waypoint)
                                        dense_points.append(waypoint)
                                    if self.inLayer.crs() != wgs84crs:
                                        endPt = transfromwgs84.transform(endPt)
                                    dense_points.append(endPt)
                                startPt = endPt

                    if geomType == "MultiPolygon":
                        dense_features = []
                        for i in range(len(multipoly_geom)):
                            dense_points = []
                            poly = multipoly_geom[i][0]
                            pointCount = len(poly)
                            startPt = QgsPointXY(poly[0][0], poly[0][1])
                            dense_points.append(startPt)
                            for j in range(1, pointCount):
                                endPt = QgsPointXY(poly[j][0], poly[j][1])
                                if self.inLayer.crs() != wgs84crs:
                                    startPt = transtowgs84.transform(startPt)
                                    endPt = transtowgs84.transform(endPt)
                                # create a geographiclib line object
                                lineObject = self.geod.InverseLine(startPt.y(), startPt.x(), endPt.y(), endPt.x())
                                # determine how many densified segments there will be
                                n = int(math.ceil(lineObject.s13 / self.spacing))
                                if lineObject.s13 > self.spacing:
                                    seglen = lineObject.s13 / n
                                    for k in range(1, n):
                                        s = seglen * k
                                        g = lineObject.Position(s, Geodesic.LATITUDE | Geodesic.LONGITUDE | Geodesic.LONG_UNROLL)
                                        waypoint = QgsPointXY(g['lon2'], g['lat2'])
                                        if self.inLayer.crs() != wgs84crs:
                                            waypoint = transfromwgs84.transform(waypoint)
                                        dense_points.append(waypoint)
                                    if self.inLayer.crs() != wgs84crs:
                                        endPt = transfromwgs84.transform(endPt)
                                    dense_points.append(endPt)
                                startPt = endPt
                            dense_features.append(dense_points)

                    new_poly = QgsFeature()
                    if geomType == "LineString":
                        new_poly.setGeometry(QgsGeometry.fromPolylineXY(dense_points))
                    elif geomType == "MultiLineString":
                        new_poly.setGeometry(QgsGeometry.fromMultiPolylineXY(dense_features))
                    elif geomType == "Polygon":
                        new_poly.setGeometry(QgsGeometry.fromPolygonXY([dense_points]))
                    elif geomType == "MultiPolygon":
                        new_poly.setGeometry(QgsGeometry.fromMultiPolygonXY([dense_features]))
                    new_poly.setAttributes(feature.attributes())
                    provider.addFeatures([new_poly])
                if badGeom > 0:
                    self.iface.messageBar().pushWarning("", "{} features failed".format(badGeom))

            if self.create_point:
                densifyPoint(self.inLayer, provider)
                out_point_layer.reload()

            if self.create_polyline:
                densifyPoly(self.inLayer, provider)
                out_line_layer.reload()

            if self.create_polygon:
                densifyPoly(self.inLayer, provider)
                out_poly_layer.reload()
