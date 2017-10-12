# -*- coding: utf-8 -*-
"""
/***************************************************************************
 GeodesicDensifier
                                 A QGIS plugin
 Adds vertices to geometry along geodesic lines
                              -------------------
        begin                : 2017-10-06
        git sha              : $Format:%H$
        copyright            : (C) 2017 by Jonah Sullivan
        email                : jonahsullivan79@gmail.com
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
    from geographiclib.geodesic import Geodesic
except ImportError:
    import sys
    import inspect
    import os

    sys.path.append(os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda: 0))))
    from geographiclib.geodesic import Geodesic
import math
from qgis.core import *
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
        # initialize locale
        locale = QSettings().value('locale/userLocale')[0:2]
        locale_path = os.path.join(
            self.plugin_dir,
            'i18n',
            'GeodesicDensifier_{}.qm'.format(locale))

        if os.path.exists(locale_path):
            self.translator = QTranslator()
            self.translator.load(locale_path)

            if qVersion() > '4.3.3':
                QCoreApplication.installTranslator(self.translator)

        # Declare instance attributes
        self.actions = []
        self.menu = self.tr(u'&Geodesic Densifier')
        self.toolbar = self.iface.addToolBar(u'GeodesicDensifier')
        self.toolbar.setObjectName(u'GeodesicDensifier')

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
        return QCoreApplication.translate('GeodesicDensifier', message)

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
            text=self.tr(u'Geodesic Densifier'),
            callback=self.run,
            parent=self.iface.mainWindow())

    def unload(self):
        """Removes the plugin menu item and icon from QGIS GUI."""
        for action in self.actions:
            self.iface.removePluginMenu(
                self.tr(u'&Geodesic Densifier'),
                action)
            self.iface.removeToolBarIcon(action)
        # remove the toolbar
        del self.toolbar

    def run(self):
        """Run method that performs all the real work"""

        # find all of the layers in the map
        layers = []
        for i in range(self.iface.mapCanvas().layerCount()):
            layer = self.iface.mapCanvas().layer(i)
            if layer.type() == layer.VectorLayer:
                layers.append(layer)

        # list of layer names
        lyr_name_list = [layer.name() for layer in layers]

        # clear the layer combo box
        self.dlg.inLayerComboBox.clear()

        # add layer names to layer combo box
        self.dlg.inLayerComboBox.addItem('')
        self.dlg.inLayerComboBox.addItems(lyr_name_list)

        # show the dialog
        self.dlg.show()

        # set the layer to process
        # create an empty layer object
        self.inLayer = QgsVectorLayer()

        def set_in_layer():
            in_layer_name = self.dlg.inLayerComboBox.currentText()
            for i in range(self.iface.mapCanvas().layerCount()):
                layer = self.iface.mapCanvas().layer(i)
                if layer.name() == in_layer_name:
                    if layer.crs().geographicFlag():
                        self.inLayer = layer
                        self.dlg.messageBox.setText("Input Layer Set: " + str(in_layer_name))
                        for field in layer.pendingFields():
                            self.dlg.uidFieldComboBox.addItem(field.name())
                    else:
                        self.dlg.messageBox.setText("Error: Input must be in Geographic coordinates")

        # listener to set input layer when combo box changes
        self.dlg.inLayerComboBox.currentIndexChanged.connect(set_in_layer)

        # set the uid field
        self.in_uid_field = QgsField()

        # add field
        def set_in_field():
            for field in self.inLayer.pendingFields():
                if field.name == self.dlg.uidFieldComboBox.currentText():
                    self.in_uid_field = field

        # listener to set input uid field when combo box changes
        self.dlg.uidFieldComboBox.currentIndexChanged.connect(set_in_field)

        # clear the ellipsoid combobox
        self.dlg.EllipsoidcomboBox.clear()

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
            self.spacing = int(self.dlg.spacingSpinBox.value())
            self.dlg.messageBox.setText("Point spacing set to " + str(self.spacing) + "m")

        # listener to set input point spacing when spin box changes
        self.dlg.spacingSpinBox.valueChanged.connect(set_in_spacing)

        # Run the dialog event loop
        result = self.dlg.exec_()
        # See if OK was pressed
        if result:

            # set output types
            if self.dlg.pointCheckBox.isChecked():
                self.create_point = True
            else:
                self.create_point = False
            if self.dlg.polygonCheckBox.isChecked():
                self.create_polygon = True
            else:
                self.create_polygon = False
            if self.dlg.polylineCheckBox.isChecked():
                self.create_polyline = True
            else:
                self.create_polyline = False

            # Create a geographiclib Geodesic object
            self.geod = Geodesic(self.ellipsoid_a, 1 / self.ellipsoid_f)

            # create empty lists to hold output points
            self.dens_point_list = []
            self.dens_line_list = []

            def densify_points(lat1, lon1, lat2, lon2, ds, segment):
                # create a geographiclib line object
                line_object = self.geod.InverseLine(lat1, lon1, lat2, lon2)
                # determine how many densified segments there will be
                n = int(math.ceil(line_object.s13 / ds)) + 1
                # adjust the spacing distance
                ds = line_object.s13 / n
                # this variable gives each point an ID
                point_id = 0
                # this variable tracks how far along the line we are
                dist = 0.0
                # this variable tracks whether it is an original or densified point
                # loop through all of the densified segments
                for i in range(n + 1):
                    if i == 0 or i == n:
                        point_type = "Original"
                    else:
                        point_type = "Densified"
                    g = line_object.Position(dist, Geodesic.STANDARD)
                    if i == n:
                        dist += ds
                        point_id += 1
                    else:
                        # add points to the line
                        self.dens_point_list.append([segment, point_id, g['lon2'], g['lat2'], point_type])
                        dist += ds
                        point_id += 1

            def densify_line(lat1, lon1, lat2, lon2, ds, segment):
                # create a geographiclib line object
                line_object = self.geod.InverseLine(lat1, lon1, lat2, lon2)
                # determine how many densified segments there will be
                n = int(math.ceil(line_object.s13 / ds)) + 1
                # adjust the spacing distance
                ds = line_object.s13 / n
                # this variable tracks how far along the line we are
                dist = 0.0
                # loop through all of the densified segments
                seg = []
                for i in range(n + 1):
                    g = line_object.Position(dist, Geodesic.STANDARD)
                    seg.append([g['lon2'], g['lat2']])
                    dist += ds
                self.dens_line_list.append([segment, seg])

            # get geometry of input point layer
            in_point_list = []
            layer_iter = self.inLayer.getFeatures()
            for feature in layer_iter:
                geom = feature.geometry()
                uid_idx = self.inLayer.fieldNameIndex(self.in_uid_field.name())
                uid = feature.attributes()[uid_idx]
                if geom.type() == QGis.Point:
                    x = geom.asPoint()
                    in_point_list.append([uid, x])

            # get input projection
            in_crs = self.inLayer.crs().authid()

            if self.create_point:

                if len(in_point_list) <= 2:  # densify a pair of points
                    from_point = in_point_list[0]
                    to_point = in_point_list[1]
                    segment = str(in_point_list[0][0])
                    densify_points(from_point[1][1],
                                   from_point[1][0],
                                   to_point[1][1],
                                   to_point[1][0],
                                   self.spacing,
                                   segment)
                    # add last point
                    self.dens_point_list.append([segment,
                                                 len(self.dens_point_list) + 1,
                                                 to_point[1][0],
                                                 to_point[1][1],
                                                 'Original'])

                else:  # densify more than two points and go back to the start
                    # create list of pairs as tuples
                    pair_list = zip(*[in_point_list[i:] + in_point_list[:i] for i in range(2)])
                    for pair in pair_list:
                        from_point = pair[0]
                        to_point = pair[1]
                        segment = str(from_point[0])
                        densify_points(from_point[1][1],
                                       from_point[1][0],
                                       to_point[1][1],
                                       to_point[1][0],
                                       self.spacing,
                                       segment)

                # create and add to map canvas a point memory layer
                layer_name = "Densified Point " + str(self.ellipsoid_name) + " " + str(self.spacing) + "m"
                out_point_layer = self.iface.addVectorLayer("Point?crs={0}".format(in_crs),
                                                            layer_name,
                                                            "memory")

                # set data provider
                pr = out_point_layer.dataProvider()
                # add attribute fields
                pr.addAttributes([QgsField("Segment", QVariant.String),
                                  QgsField("ID", QVariant.String),
                                  QgsField("LAT", QVariant.Double),
                                  QgsField("LON", QVariant.Double),
                                  QgsField("PntType", QVariant.String),
                                  QgsField("DensT", QVariant.String)])
                out_point_layer.updateFields()

                # loop through points adding geometry and attributes
                for point_object in self.dens_point_list:
                    # create a feature
                    feat = QgsFeature(out_point_layer.pendingFields())
                    # set geometry to the feature
                    x_value = point_object[2]
                    y_value = point_object[3]
                    feat.setGeometry(QgsGeometry.fromPoint(QgsPoint(x_value, y_value)))
                    # set attribute fields
                    feat.setAttribute("Segment", str(point_object[0]))
                    feat.setAttribute("ID", point_object[1])
                    feat.setAttribute("LAT", float(y_value))
                    feat.setAttribute("LON", float(x_value))
                    feat.setAttribute("PntType", point_object[4])
                    feat.setAttribute("DensT", "G")
                    out_point_layer.dataProvider().addFeatures([feat])

            if self.create_polyline:

                if len(in_point_list) <= 2:  # densify one line segment
                    from_point = in_point_list[0]
                    to_point = in_point_list[1]
                    segment = str(in_point_list[0][0])
                    densify_line(from_point[1][1],
                                 from_point[1][0],
                                 to_point[1][1],
                                 to_point[1][0],
                                 self.spacing,
                                 segment)

                else:  # densify more than one line segment and back to the start
                    # create list of pairs as tuples
                    pair_list = zip(*[in_point_list[i:] + in_point_list[:i] for i in range(2)])
                    for pair in pair_list:
                        from_point = pair[0]
                        to_point = pair[1]
                        segment = str(from_point[0])
                        densify_line(from_point[1][1],
                                     from_point[1][0],
                                     to_point[1][1],
                                     to_point[1][0],
                                     self.spacing,
                                     segment)

                # create and add to map canvas a polyline memory layer
                layer_name = "Densified Line " + str(self.ellipsoid_name) + " " + str(self.spacing) + "m"
                out_line_layer = self.iface.addVectorLayer("LineString?crs={0}".format(in_crs),
                                                           layer_name,
                                                           "memory")
                # set data provider
                pr = out_line_layer.dataProvider()
                # add attribute fields
                pr.addAttributes([QgsField("Segment", QVariant.String),
                                  QgsField("DensT", QVariant.String)])
                out_line_layer.updateFields()

                for segmentObject in self.dens_line_list:
                    feat = QgsFeature(out_line_layer.pendingFields())
                    qgs_point_list = []
                    for i in range(len(segmentObject[1])):
                        point_object = segmentObject[1][i]
                        qgs_point_list.append(QgsPoint(point_object[0], point_object[1]))
                    feat.setGeometry(QgsGeometry.fromPolyline(qgs_point_list))
                    # set attribute fields
                    feat.setAttribute("Segment", str(segmentObject[0]))
                    feat.setAttribute("DensT", "G")
                    out_line_layer.dataProvider().addFeatures([feat])

            if self.create_polygon:

                if len(in_point_list) > 2: # need three points to make a polygon

                    # create list of pairs as tuples
                    pair_list = zip(*[in_point_list[i:] + in_point_list[:i] for i in range(2)])
                    for pair in pair_list:
                        from_point = pair[0]
                        to_point = pair[1]
                        segment = str(from_point[0])
                        densify_line(from_point[1][1],
                                     from_point[1][0],
                                     to_point[1][1],
                                     to_point[1][0],
                                     self.spacing,
                                     segment)

                    # create and add to map canvas a polyline memory layer
                    layer_name = "Densified Polygon " + str(self.ellipsoid_name) + " " + str(self.spacing) + "m"
                    out_poly_layer = self.iface.addVectorLayer("Polygon?crs={0}".format(in_crs),
                                                               layer_name,
                                                               "memory")
                    # set data provider
                    pr = out_poly_layer.dataProvider()
                    # add attribute fields
                    pr.addAttributes([QgsField("DensT", QVariant.String)])
                    out_poly_layer.updateFields()

                    qgs_point_list = []
                    feat = QgsFeature(out_poly_layer.pendingFields())
                    for segmentObject in self.dens_line_list:
                        for i in range(len(segmentObject[1])):
                            point_object = segmentObject[1][i]
                            qgs_point_list.append(QgsPoint(point_object[0], point_object[1]))
                    feat.setGeometry(QgsGeometry.fromPolygon([qgs_point_list]))
                    # set attribute fields
                    feat.setAttribute("DensT", "G")
                    out_poly_layer.dataProvider().addFeatures([feat])
