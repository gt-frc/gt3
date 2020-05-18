#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from contours.quad import QuadContourGenerator
from scipy.interpolate import griddata
from shapely.geometry import Point, MultiPoint, LineString


def find_xpt_mag_axis(R, Z, psi):

    # calculate gradients of psi in R and Z directions
    dpsidR = np.gradient(psi, R[0, :], axis=1)
    dpsidZ = np.gradient(psi, Z[:, 0], axis=0)

    # find line(s) where dpsidR=0
    c_dpsidR = QuadContourGenerator.from_rectilinear(R[0], Z[:, 0], dpsidR)
    c_dpsidZ = QuadContourGenerator.from_rectilinear(R[0], Z[:, 0], dpsidZ)
    dpsidR_0 = c_dpsidR.contour(0.0)
    dpsidZ_0 = c_dpsidZ.contour(0.0)

    # populate list all intersection points (as Shapely points)
    int_pts = []
    for i, path1 in enumerate(dpsidR_0):
        for j, path2 in enumerate(dpsidZ_0):
            try:
                # find intersection points between curves for dpsidR=0 and dpsidZ=0
                # these correspond to local minima, local maxima, or saddle-points in psi
                ints = LineString(path1).intersection(LineString(path2))
                if ints.type == 'Point':
                    int_pts.append(ints)
                elif ints.type == 'MultiPoint':
                    for pt in ints:
                        int_pts.append(pt)
            except:
                pass

    # magnetic axis is the intersection point with the minimum value of psi
    psi_int_pts = griddata(np.column_stack((R.flatten(), Z.flatten())),
                           psi.flatten(),
                           np.asarray(MultiPoint(int_pts)),
                           method='cubic')
    mag_axis = np.asarray(MultiPoint(int_pts))[np.argmin(psi_int_pts)]

    # find the dpsidR_0 = 0 contour that contains the magnetic axis.
    # Of the other intersection points along the same dpsidR_0 = 0 contour, the x-point
    # will be the highest (in Z) of those points that are below the magnetic axis.

    # create list of potential x-points (i.e. points that also fall on the same dpsidR_0 = 0 contour
    # as the magnetic axis and whose y-values are lower than that of the magnetic axis.
    pot_xpts = []
    for i, path in enumerate(dpsidR_0):
        line = LineString(path)
        if line.distance(Point(mag_axis)) < 1E-8:
            # then we've found the contour that includes the magnetic axis
            for pt in int_pts:
                if line.distance(pt) < 1E-8 and pt.y < Point(mag_axis).y:
                    pot_xpts.append(pt)

    # of the potential x-points, take the one with the largest y value (there might be only one potential x-point)
    pts_xpts_array = np.asarray(MultiPoint(pot_xpts))
    xpt = pts_xpts_array[pts_xpts_array[:, 1].argsort()][-1]

    # plt.contourf(R, Z, psi)
    # for pt in pts_xpts_array:
    #     plt.scatter(pt[0], pt[1], color='yellow')
    # plt.scatter(xpt[0], xpt[1],marker='+',color='black')
    # plt.show()
    # sys.exit()
    return xpt, mag_axis
