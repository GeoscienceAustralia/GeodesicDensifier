# -*- coding: utf-8 -*-
"""
/***************************************************************************
 GeodesicDensifierDialog
                                 A QGIS plugin
 Adds vertices to geometry along geodesic lines
                             -------------------
        begin                : 2017-10-06
        git sha              : $Format:%H$
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

import os

from PyQt4 import QtGui, uic

FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'geodesic_densifier_dialog_base.ui'))


class GeodesicDensifierDialog(QtGui.QDialog, FORM_CLASS):
    def __init__(self, parent=None):
        """Constructor."""
        super(GeodesicDensifierDialog, self).__init__(parent)
        self.setupUi(self)
