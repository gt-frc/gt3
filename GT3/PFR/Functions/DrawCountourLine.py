#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from contours.quad import QuadContourGenerator

def draw_contour_line(R, Z, array, val, pathnum):
    c = QuadContourGenerator.from_rectilinear(R[0], Z[:, 0], array)

    res = c.contour(val)[pathnum]
    x = res[:, 0]
    y = res[:, 1]
    return x, y