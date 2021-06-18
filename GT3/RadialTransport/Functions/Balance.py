#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from scipy.interpolate import UnivariateSpline


def balance(gamma, interm, r, sa, x):
    return r[x], UnivariateSpline(r, interm, k=3, s=0).integral(0., r[x]), gamma[x]*sa(r[x])
