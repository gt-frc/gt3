#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np

def calc_cxcool(core, n, T):
    """Calculates charge exchange cooling with the slow neutrals"""
    slowNeutDens = np.array([n.n.s[-1] * x ** 15 for x in core.r[:, 0] / core.a])
    totalNeutDens = np.array([n.n.tot[-1] * x ** 15 for x in core.r[:, 0] / core.a])
    #result = 1.5 * n.i * T.i.J * slowNeutDens* ((core.sv.el.st[:,0] + core.sv.el.st[:,0]* (totalNeutDens/n.i)) + core.sv.cx.st[:,0])
    #elastW = 10**(-14 + np.log10(T.i.ev.val))
    elastW = 0.
    result = 1.5 * n.i.val * T.i.J.val * slowNeutDens * ((core.sv.el.st[:, 0] +
                                                          core.sv.cx.st[:, 0] +
                                                           + elastW * slowNeutDens / n.i.val)) / 1E4
            # TODO: Work this out for n-n elastic scattering
    return result