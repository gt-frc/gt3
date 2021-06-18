#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from shapely.geometry import  Point
from math import pi

def getangle(p1, p2):
    if isinstance(p1, Point) and isinstance(p2, Point):
        p1 = [p1.coords.xy[0][0], p1.coords.xy[1][0]]
        p2 = [p2.coords.xy[0][0], p2.coords.xy[1][0]]
    p1 = np.asarray(p1)
    p1 = np.reshape(p1, (-1, 2))
    p2 = np.asarray(p2)
    p2 = np.reshape(p2, (-1, 2))
    theta = np.arctan2(p1[:, 1]-p2[:, 1], p1[:, 0]-p2[:, 0])
    theta_mod = np.where(theta<0, theta+pi, theta)  # makes it so the angle is always measured counterclockwise from the horizontal
    return theta
