#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
import sys
from math import cos, sin
from shapely.geometry import Point, LineString
from GT3.Core.Functions.DrawCoreLine import draw_core_line


def calc_RZ(rho, theta, theta_xpt, pts, psi_data, psi_norm, lines):
    # get parameters that depend on both rho and theta

    sep_pts = np.asarray(lines.sep_closed.coords)

    R = np.zeros(rho.shape)
    Z = np.zeros(rho.shape)

    line_length = 5.0

    if pts.xpt[1] is None or pts.xpt[0] is None:
        # Grab the actual x-point
        if pts.xpt[0] is not None:
            xpt_temp = pts.xpt[0]
        else:
            xpt_temp = pts.xpt[1]
    else:
        # THere are 2 xpoints. Use the bottom one.
        xpt_temp = pts.xpt[0]

    for i, psi_norm_val in enumerate(psi_norm[:, 0]):
        if i == 0:
            R[i] = pts.axis.mag[0]
            Z[i] = pts.axis.mag[1]
        else:
            # attempt to draw flux surface line through that point
            # (may not work for flux surfaces very close to the magnetic axis)
            fs_line, fs_pts = draw_core_line(psi_data.R, psi_data.Z, psi_data.psi_norm, psi_norm_val, sep_pts)
            if fs_line is not None and fs_pts.axis is None:
                # then draw_core_line didn't find any contours, which probably means it was trying
                # to draw a surface closer to psi_norm=0 than the contours package would cooperate with.
                # When this happens, the best thing to do for now is decrease the radial resolution in
                # the vicinity of the magnetic axis. Ideas for the fixing this in the future:
                #   1) resample the raw psi data (maybe an Rbf interpolation) onto a finer mesh. May or may not work.
                #   2) some kind of interpolation based on the flux surfaces it can draw.  # TODO
                raise RuntimeError("""GT3 had trouble drawing a contour line when getting the R and Z points. This ' \
                      'is most likely due to an overly fine radial mesh in the vicnity of the magnetic ' \
                      'axis. Try reducing your number of radial meshes in the core and try again. This ' \
                      'will hopefully be fixed in a future update.""")

            for j, thetaval in enumerate(theta[0]):
                if psi_norm_val < 1.0:
                    thetaline = LineString([Point(fs_pts.axis),
                                            Point([line_length * cos(thetaval) + fs_pts.axis[0],
                                                   line_length * sin(thetaval) + fs_pts.axis[1]])])
                    int_pt = fs_line.intersection(thetaline)
                else:
                    if thetaval == theta_xpt:
                        int_pt = Point(xpt_temp)

                    else:
                        thetaline = LineString([Point(pts.axis.geo),
                                                Point([line_length * cos(thetaval) + pts.axis.geo[0],
                                                       line_length * sin(thetaval) + pts.axis.geo[1]])])
                        int_pt = LineString(sep_pts).intersection(thetaline)
                if isinstance(int_pt, Point):
                    R[i, j] = int_pt.x
                    Z[i, j] = int_pt.y
                else:
                    raise
    return R, Z