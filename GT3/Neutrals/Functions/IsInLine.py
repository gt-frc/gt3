#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from shapely.geometry import Point

def isinline(pt, line):
    pt_s = Point(pt)
    dist = line.distance(pt_s)
    if dist < 1E-6:
        return True
    else:
        return False