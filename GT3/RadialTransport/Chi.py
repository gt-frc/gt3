#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from scipy.interpolate import UnivariateSpline
import numpy as np
from collections import namedtuple
from scipy import constants
from math import sqrt
from deprecation import deprecated



m_d = constants.physical_constants['deuteron mass'][0]
e = constants.elementary_charge
z_d = 1  # atomic number of deuterium
z_c = 6  # atomic number of carbon
ch_d = e * z_d  # charge of deuterium
ch_c = e * z_c  # charge of carbon
from GT3 import Core


@deprecated(deprecated_in="0.0.3", removed_in="0.0.4", details="Chi has been folded back into RadialTransport")
class Chi:
    """
    Chi class claculates various chis and provides source terms as necessary
    """
    def __init__(self, data, core, n, T, L, reInterp=False):

        self.Qi = data.Qi_diff
        self.Qe = data.Qe_diff
        self.conv15 = UnivariateSpline(core.r[:,0], .5 * ch_d * data.gamma_diff_D * T.i.ev, k=3, s=0)(core.r[:,0])
        self.conv25 = UnivariateSpline(core.r[:,0], .5 * ch_d * data.gamma_diff_D * T.i.ev, k=3, s=0)(core.r[:,0])
        #self.heatvisc = self.viscCalc(data, core, n, T)
        self.heatvisc = np.zeros(self.conv25.shape)
        self.heatin = UnivariateSpline(core.r[:,0], .5 * data.gamma_diff_D * m_d * (data.vtor_D_total**2 + data.vpol_D**2), k=3, s=0)(core.r[:,0]) # TODO: Provide logic that uses vtor_D_intrin/fluid depending on IOL Switch, currently too small to matter
        if reInterp:
            L_T_i = UnivariateSpline(core.r[:,0], L.T.i, k=3, s=0)(core.r[:,0])
            T_i_ev = UnivariateSpline(core.r[:,0], T.i.ev, k=3, s=0)(core.r[:,0])
            n_i = UnivariateSpline(core.r[:,0], n.i, k=3, s=0)(core.r[:,0])
            self.chi = namedtuple('chi','i e')(
                namedtuple('i','chi1 chi2 chi3 chi4')(
                    UnivariateSpline(core.r[:,0], (self.Qi) * L_T_i / (n_i * T_i_ev * ch_d), k=1, s=5)(core.r[:,0]),
                    UnivariateSpline(core.r[:,0], (self.Qi - self.conv25) * L_T_i / (n_i * T_i_ev * ch_d), k=1, s=5)(core.r[:,0]),
                    UnivariateSpline(core.r[:,0], ((self.Qi - self.conv25 - self.heatin) * L_T_i / (n_i * T_i_ev * ch_d)), k=1, s=5)(core.r[:,0]),
                    UnivariateSpline(core.r[:,0], ((self.Qi - self.conv25 - self.heatin * self.heatvisc) * L_T_i / (n_i * T_i_ev * ch_d)), k=1, s=5)(core.r[:,0])
                ),self.calc_chi_e(data, n, L, T)
            )
        else:
            self.chi = namedtuple('chi','i e')(
                namedtuple('i','chi1 chi2 chi3 chi4')(
                    (self.Qi) * L.T.i / (n.i * T.i.ev * ch_d),
                    (self.Qi - self.conv25) * L.T.i / (n.i * T.i.ev * ch_d),
                    (self.Qi - self.conv25 - self.heatin) * L.T.i / (n.i * T.i.ev * ch_d),
                    (self.Qi - self.conv25 - self.heatin * self.heatvisc) * L.T.i / (n.i * T.i.ev * ch_d)
                ),self.calc_chi_e(data, n, L, T)
            )


    def calc_chi_e(self, data, n, L, T):
        gameltemp = 1.0* data.gamma_diff_D + 6.0 * data.gamma_C   #OHHHH gamma electron!!! HA HA HA HA WE DON'T KNOW WHAT THE FUCK GAMMA_C IS

        return L.T.e * ((self.Qe / (ch_d * n.e * T.e.ev)) - 2.5 * gameltemp / n.e)

    def viscCalc(self, data, core: Core, n, T):
        f1 = 1 # Must determine later, appears to be geometeric factor in Eq(4) in POP052504(2010)
        fp = core.B.pol.fsa/core.B.tor.fsa
        vth = [sqrt(2. * a  * ch_d / m_d) for a in T.i.ev]
        eta0 = [a * m_d * b * c * core.R0_a * f1 for a,b,c in zip(n.i, vth, core.q[:,0])]
        eta4 = [a * m_d * c * ch_d / (ch_d * abs(b)) for a, b, c in zip(n.i, core.B.tor.fsa, T.i.ev)]
        vrad = data.gamma_diff_D / n.i

        # a = vtor    b = fp   c = eta0
        # d = vrad    f = vthet g = eta 4

        return [a * (b * c * d - .5 * g * (4.0 * a + f)) - .5 * f * (c * d + g * (a + .5 * f)) for a, b,c,d,f,g in zip(data.vtor_D_total, fp, eta0, vrad, data.vpol_D, eta4)]

