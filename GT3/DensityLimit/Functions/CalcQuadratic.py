#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from math import sqrt

def calc_quadratic(a, b, c):
    y1 = -b * (1 + sqrt(1 - 4 * (a * c / b ** 2))) / (2 * a)
    y2 = -b * (1 - sqrt(1 - 4 * (a * c / b ** 2))) / (2 * a)
    return y1, y2