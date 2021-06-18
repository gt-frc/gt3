#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from Core.Functions.DrawCounterLine import draw_contour_line
from contours.quad import QuadContourGenerator
from shapely.geometry import LineString
import numpy as np


def draw_core_line(R, Z, psi, psi_val, sep_pts):
    # create contour generator
    c = QuadContourGenerator.from_rectilinear(R[0], Z[:, 0], psi)

    # draw contours with psi_val
    contours = c.contour(psi_val)

    if len(contours) == 1:
        # then we're definitely dealing with a surface inside the seperatrix
        x, y = draw_contour_line(R, Z, psi, psi_val, 0)
    else:
        # we need to find which of the surfaces is inside the seperatrix
        for j, line in enumerate(contours):
            x, y = draw_contour_line(R, Z, psi, psi_val, j)

            if (np.amax(x) < np.amax(sep_pts[:, 0]) and
                np.amin(x) > np.amin(sep_pts[:, 0]) and
                np.amax(y) < np.amax(sep_pts[:, 1]) and
                np.amin(y) > np.amin(sep_pts[:, 1])):
                # then it's an internal flux surface
                break
    pts = np.column_stack((x, y))
    line = LineString(pts)
    out_pt = pts[np.argmax(pts, axis=0)[0]]
    in_pt = pts[np.argmin(pts, axis=0)[0]]
    top_pt = pts[np.argmax(pts, axis=0)[1]]
    bot_pt = pts[np.argmin(pts, axis=0)[1]]
    fs_axis = [(out_pt[0]+in_pt[0])/2, (out_pt[1]+in_pt[1])/2]
    return line, fs_axis