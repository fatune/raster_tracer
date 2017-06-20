# -*- coding: utf-8 -*-
"""
/***************************************************************************
 raster_tracer
                                 A QGIS plugin
 Raster tracer
                             -------------------
        begin                : 2017-06-08
        copyright            : (C) 2017 by Mikhail Kondratyev
        email                : fatune@gmail.com
        git sha              : $Format:%H$
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 This script initializes the plugin, making it known to QGIS.
"""


# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load raster_tracer class from file raster_tracer.

    :param iface: A QGIS interface instance.
    :type iface: QgsInterface
    """
    #
    from .raster_tracer import raster_tracer
    return raster_tracer(iface)
