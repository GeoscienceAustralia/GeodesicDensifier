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
 *   it under the terms of the Apache 2.0 License.                         *
 *                                                                         *
 ***************************************************************************/
"""

import os

from PyQt5 import uic
from PyQt5 import QtWidgets

FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'geodesic_densifier_dialog_base.ui'))


class GeodesicDensifierDialog(QtWidgets.QDialog, FORM_CLASS):
    def __init__(self, parent=None):
        """Constructor."""
        super(GeodesicDensifierDialog, self).__init__(parent)
        self.setupUi(self)
