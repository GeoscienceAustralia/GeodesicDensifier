# -*- coding: utf-8 -*-
"""
/***************************************************************************
 GeodesicDensifier
                                 A QGIS plugin
 Adds vertices to geometry along geodesic lines
                             -------------------
        begin                : 2017-10-06
        copyright            : (C) 2017 by Jonah Sullivan
        email                : jonah.sullivan@ga.gov.au
        git sha              : $Format:%H$
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the Apache 2.0 license                          *
 *                                                                         *
 ***************************************************************************/
 This script initializes the plugin, making it known to QGIS.
"""


# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load GeodesicDensifier class from file GeodesicDensifier.

    :param iface: A QGIS interface instance.
    :type iface: QgsInterface
    """
    #
    from .geodesic_densifier import GeodesicDensifier
    return GeodesicDensifier(iface)
