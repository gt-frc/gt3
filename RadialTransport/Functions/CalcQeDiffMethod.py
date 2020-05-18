#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from scipy.interpolate import UnivariateSpline


def calc_Qe_diff_method(r, a, en_src_nbi_e, cool_rate, Qie):
    en_src_nbi_eint = UnivariateSpline(r, en_src_nbi_e, k=3, s=0)
    cool_rateint = UnivariateSpline(r, cool_rate, k=3, s=0)
    Qie_int = UnivariateSpline(r, Qie, k=3, s=0)

    def f(t, flux, Qie, Q_e_nbi, cool_rate):
        S = Q_e_nbi(t) - Qie(t) - cool_rate(t)
        return S - flux * (1 / (t + 0.003))

    from scipy.integrate import ode

    flux = ode(f).set_integrator('vode', with_jacobian=False)
    flux.set_initial_value(0., 0.).set_f_params(Qie_int, en_src_nbi_eint, cool_rateint)

    dt = a / len(r)
    x, y = [], []
    while flux.successful() and flux.t < a:
        x.append(flux.t + dt)
        y.append(flux.integrate(flux.t + dt))
    flux = UnivariateSpline(x, y, k=3, s=0)

    # print "Total volume in Qi_diff calc: " + str(UnivariateSpline(r, dVdr(r), k=3, s=0).integral(0., a))
    # print "Total nbi ion energy: " + str(UnivariateSpline(r, (en_src_nbi_keptint(r) + en_src_nbi_lostint(r)) * dVdr(r), k=3, s=0).integral(0., 1.)/(1E6))+" MW"
    return flux(r)
