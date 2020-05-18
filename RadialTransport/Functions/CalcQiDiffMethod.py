#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from scipy.interpolate import UnivariateSpline
from scipy.integrate import ode

def calc_Qi_diff_method(r, a, cxcool, Qie, en_src_nbi_kept, en_src_nbi_lost, dVdrho, iol_adjusted=False, E_orb=None):
    dVdr = UnivariateSpline(r, dVdrho(r / a) / a, k=3, s=0)
    dE_orb = UnivariateSpline(r, E_orb, k=3, s=0).derivative()
    en_src_nbi_keptint = UnivariateSpline(r, en_src_nbi_kept, k=3, s=0)
    en_src_nbi_lostint = UnivariateSpline(r, en_src_nbi_lost, k=3, s=0)
    Qie_int = UnivariateSpline(r, Qie, k=3, s=0)
    cxcoolint = UnivariateSpline(r, cxcool, k=3, s=0)

    def f(t, flux, cxcool, Qie, Q_i_nbi, Qnbi_loss, dEdr, iolFlag):
        S = Q_i_nbi(t) + cxcool(t) + Qie(t)
        if iolFlag:
            return S - Qnbi_loss(t) - (flux * (dEdr(t) + 1 ) / (t + 0.003))
        else:
            return S - (flux * (dEdr(t) + 1 ) / (t + 0.003))


    flux = ode(f).set_integrator('vode', with_jacobian=False)
    flux.set_initial_value(0., 0.).set_f_params(cxcoolint, Qie_int, en_src_nbi_keptint, en_src_nbi_lostint, dE_orb,
                                                iol_adjusted)
    dt = a / len(r)
    x, y = [], []
    while flux.successful() and flux.t < a:
        x.append(flux.t + dt)
        y.append(flux.integrate(flux.t + dt))
    flux = UnivariateSpline(x, y, k=3, s=0)

    print "Total volume in Qi_diff calc: " + str(UnivariateSpline(r, dVdr(r), k=3, s=0).integral(0., a))
    print "Total nbi ion energy: " + str(UnivariateSpline(r, (en_src_nbi_keptint(r) + en_src_nbi_lostint(r)) * dVdr(r), k=3, s=0).integral(0., 1.)/(1E6))+" MW"
    return flux(r)