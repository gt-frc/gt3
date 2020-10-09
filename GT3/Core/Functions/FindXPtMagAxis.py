#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point, MultiPoint, LineString

def find_xpt_mag_axis(core, R, Z, psi):
    # type: (GT3.core, np.array, np.array, np.array) -> (np.array, np.array, np.array)
    """

    :type core: .core
    """
    xpt_u = None
    xpt_l = None
    m_axis = None
    # find x-point location
    dpsidR = np.gradient(psi, R[0, :], axis=1)
    dpsidZ = np.gradient(psi, Z[:, 0], axis=0)
    d2psidR2 = np.gradient(dpsidR, R[0, :], axis=1)
    d2psidZ2 = np.gradient(dpsidZ, Z[:, 0], axis=0)

    # find line(s) where dpsidR=0
    csR = plt.contour(R, Z, dpsidR, [0])
    csZ = plt.contour(R, Z, dpsidZ, [0])

    dpsidR_0 = csR.collections[0].get_paths()
    # dpsidR_0 = cntr.contour(R, Z, dpsidR).trace(0.0)

    # find line(s) where dpsidZ=0
    # dpsidZ_0 = cntr.contour(R, Z, dpsidZ).trace(0.0)
    dpsidZ_0 = csZ.collections[0].get_paths()
    for i, path1 in enumerate(dpsidR_0):
        for j, path2 in enumerate(dpsidZ_0):
            try:
                # find intersection points between curves for dpsidR=0 and dpsidZ=0
                ints = LineString(path1.vertices).intersection(LineString(path2.vertices))
                # if there is only one intersection ('Point'), then we're probably not
                # dealing with irrelevant noise in psi
                if isinstance(ints, Point):
                    # check if local maximum or minimum
                    # If the point is within the walls and less than 0.5m from the vessel centroid, this is
                    # very likely the magnetic axis
                    if core.wall_line.convex_hull.contains(ints) and core.wall_line.centroid.distance(ints) < 0.5:
                        # we've found the magnetic axis
                        m_axis = np.array([ints.x, ints.y])
                    # If the point is not inside the walls, it's a magnet
                    elif not core.wall_line.convex_hull.contains(ints):
                        # we've found a magnet. Do nothing.
                        pass
                    elif core.wall_line.convex_hull.contains(ints) and ints.y < 0:
                        # The point is within the walls and not near the magnetic axis. This is likely the
                        # lower x-point.
                        # TODO: Provides lower x-point here, but upper x-point can be built out
                        xpt = np.array([ints.x, ints.y])
                        xpt_l = xpt
                    elif core.wall_line.convex_hull.contains(ints) and ints.y > 0:
                        # This is an upper x-point. This functionality needs to be built in but is here for
                        # later
                        xpt_u = np.array([ints.x, ints.y])

                    # uncomment this line when debugging
                    # print list(ints.coords), d2psidR2(ints.x, ints.y), d2psidZ2(ints.x, ints.y)

                # If multiple points are found, the flux surfaces may have come back around onto each other.
                if isinstance(ints, MultiPoint):
                    for point in ints:
                        # check if local maximum or minimum
                        # If the point is within the walls and less than 0.5m from the vessel centroid, this is
                        # very likely the magnetic axis
                        if core.wall_line.convex_hull.contains(point) and core.wall_line.centroid.distance(
                                point) < 0.5:
                            # we've found the magnetic axis
                            m_axis = np.array([point.x, point.y])
                        # If the point is not inside the walls, it's a magnet
                        elif not core.wall_line.convex_hull.contains(point):
                            # we've found a magnet. Do nothing.
                            pass
                        elif core.wall_line.convex_hull.contains(point) and point.y < 0:
                            # The point is within the walls and not near the magnetic axis. This is likely the
                            # lower x-point.
                            # TODO: Provides lower x-point here, but upper x-point can be built out
                            xpt = np.array([point.x, point.y])
                            xpt_l = xpt
                        elif core.wall_line.convex_hull.contains(point) and point.y > 0:
                            # This is an upper x-point. This functionality needs to be built in but is here for
                            # later
                            xpt_u = np.array([point.x, point.y])

            except:
                pass
    return xpt_l, xpt_u, m_axis
