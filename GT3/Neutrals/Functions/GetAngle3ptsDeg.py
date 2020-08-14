#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from math import degrees, sqrt, acos

def getangle3ptsdeg(p1, p2, p3):
    a = sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2)
    b = sqrt((p2[0]-p3[0])**2+(p2[1]-p3[1])**2)
    c = sqrt((p1[0]-p3[0])**2+(p1[1]-p3[1])**2)
    theta = degrees(acos((c**2 - a**2 - b**2)/(-2*a*b)))  # returns degree in radians
    return theta
