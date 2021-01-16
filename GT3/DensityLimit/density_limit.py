#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 20:24:53 2018

@author: max
"""

import numpy as np
from scipy.special import jv
import sys
from math import sqrt
from .Functions.CalcQuadratic import calc_quadratic


class DensityLimit:
    def __init__(self, core, nbi):
        sys.dont_write_bytecode = True
        self.r = core.r[:, 0]
        self.a = core.a
        self.eq44(core, nbi)
        pass

    def eq44(self, core, nbi):
        def av(X, p=self):
            numerator = np.sum(p.r * X * jv(0, 5.5 * p.r / p.a) * p.a / (len(p.r) - 1))
            denominator = np.sum(p.r * jv(0, 5.5 * p.r / p.a) * p.a / (len(p.r) - 1))
            return numerator / denominator

        ni = core.n.i[:, 0]
        Ti = core.T.i.J.to1D()
        n0 = ni[0]
        n_av = av(ni)
        f = n0 / n_av

        chi = core.chi.bohm[:, 0]
        chi_hat = av(ni * chi) / n_av
        D = 0.5
        D_hat = av(ni * D) / n_av

        a = core.a
        g = ni / n0

        fz = 0.05
        Lz = core.Lz.t[:, 0]
        dLzdT = core.Lz.ddT.t[:, 0]

        sv_fus = core.sv.fus.dd[:, 0]
        dsv_fus_dT = core.sv.fus.d_dT.dd[:, 0]
        Ua = 0

        H_aux = nbi.beams.D1.dPdV.v1D.W

        H_ohm = 0
        dHdT = 3 / (2 * Ti) * (H_ohm - H_aux)
        dHdn = 0

        nn = core.n.n.tot[:, 0]
        sv_ion = core.sv.ion.st[:, 0]
        dsv_ion_dT = core.sv.ion.d_dT.st[:, 0]

        dSdn = nn * sv_ion
        dSdT = ni * nn * dsv_ion_dT

        if False:
            print('sv_ion = ', sv_ion)
            print('nn = ', nn)
            print('ni = ', ni)
            print('Ti = ', Ti)
            print('a = ', a)
            print('dSdn = ', dSdn)
            print('D_hat = ', D_hat)
            print('dHdn = ', dHdn)
            print('Ua = ', Ua)
            print('sv_fus = ', sv_fus)
            print('fz = ', fz)
            print('Lz = ', Lz)

        ya = 3 * av(Ti) * (av(dSdn) - (5.5 / a) ** 2 * D_hat) - av(dHdn + 2 * ni * (1 / 4 * Ua * sv_fus - fz * Lz))
        yb = 3 * av(ni) * (av(dSdn) - (5.5 / a) ** 2 * D_hat) + 3 * av(Ti) * av(dSdT) - av(
            dHdT + ni ** 2 * (1 / 4 * Ua * dsv_fus_dT + fz * (-dLzdT)))
        yc = 3 * av(ni) * av(dSdT)

        if False:
            print('ya = ', ya)
            print('yb = ', yb)
            print('yc = ', yc)
            print()

        y1, y2 = calc_quadratic(ya, yb, yc)

        if False:
            print()
            print('y1 = ', y1)
            print('y2 = ', y2)

        t1y1 = chi_hat * (5.5 / a) ** 2 * av(g) + 2 * y1 * (fz * av(g * Lz) - av(1 / 4 * Ua * g * sv_fus))
        t2y1 = 2 * av(g ** 2 * (1 / 4 * Ua * dsv_fus_dT + fz * (-dLzdT)))
        t3y1 = 4 * (av(-dHdT) - y1 * av(dHdn)) * av(g ** 2 * (1 / 4 * Ua * dsv_fus_dT + fz * (-dLzdT)))
        t4y1 = chi_hat * (5.5 / a) ** 2 * av(g) + 2 * y1 * (fz * av(g * Lz) - av(1 / 4 * Ua * g * sv_fus)) ** 2

        t1y2 = chi_hat * (5.5 / a) ** 2 * av(g) + 2 * y2 * (fz * av(g * Lz) - av(1 / 4 * Ua * g * sv_fus))
        t2y2 = 2 * av(g ** 2 * (1 / 4 * Ua * dsv_fus_dT + fz * (-dLzdT)))
        t3y2 = 4 * (av(-dHdT) - y2 * av(dHdn)) * av(g ** 2 * (1 / 4 * Ua * dsv_fus_dT + fz * (-dLzdT)))
        t4y2 = chi_hat * (5.5 / a) ** 2 * av(g) + 2 * y2 * (fz * av(g * Lz) - av(1 / 4 * Ua * g * sv_fus)) ** 2

        if False:
            print()
            print('t1y1 = ', t1y1)
            print('t2y1 = ', t2y1)
            print('t3y1 = ', t3y1)
            print('t4y1 = ', t4y1)
            print('t1y2 = ', t1y2)
            print('t2y2 = ', t2y2)
            print('t3y2 = ', t3y2)
            print('t4y2 = ', t4y2)

        nlim1 = (1 / f) * (t1y1 / t2y1) * (1 + sqrt(1 + t3y1 / t4y1))
        nlim2 = (1 / f) * (t1y2 / t2y2) * (1 + sqrt(1 + t3y2 / t4y2))
        nlim3 = (1 / f) * (t1y1 / t2y1) * (1 - sqrt(1 + t3y1 / t4y1))
        nlim4 = (1 / f) * (t1y2 / t2y2) * (1 - sqrt(1 + t3y2 / t4y2))

        if False:
            print()
            print('nlim1 = ', nlim1)
            print('nlim2 = ', nlim2)
            print('nlim3 = ', nlim3)
            print('nlim4 = ', nlim4)
            print()

        nlim = min(i for i in [nlim1, nlim2, nlim3, nlim4] if i > 0)

        if n_av > nlim:
            print('n_av = ', n_av)
            print('nlim = ', nlim)
            print('Disruption predicted: YES')
        else:
            print('n_av = ', n_av)
            print('nlim = ', nlim)
            print('Disruption predicted: NO')
        print()
