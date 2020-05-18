#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from scipy.interpolate import UnivariateSpline

def calc_gamma_diff_method(r, a, part_src_nbi_tot, part_src_nbi_lost, izn_rate, dVdrho, iol_adjusted=False, F_orb=None):
    dVdr = UnivariateSpline(r, dVdrho(r / a) / a, k=3, s=0)
    dF_orb = UnivariateSpline(r, F_orb, k=3, s=0).derivative()
    izn_rateint = UnivariateSpline(r, izn_rate, k=3, s=0)
    part_src_nbi_totint = UnivariateSpline(r, part_src_nbi_tot, k=3, s=0)
    part_src_nbi_lostint = UnivariateSpline(r, part_src_nbi_lost, k=3, s=0)

    def f(t, gamma, sion, snbi, snbi_loss, dFdr, iolFlag):
        S = snbi(t) + sion(t)
        if iolFlag:
            return S - snbi_loss(t) - gamma * (dFdr(t) + 1) / (t + 0.003)
        else:
            return S - gamma * (1 / (t + 0.003))

    from scipy.integrate import ode

    gamma = ode(f).set_integrator('vode', with_jacobian=False)
    gamma.set_initial_value(0., 0.).set_f_params(izn_rateint, part_src_nbi_totint, part_src_nbi_lostint, dF_orb, iol_adjusted)
    dt = a / len(r)
    x, y = [], []
    while gamma.successful() and gamma.t < a:
        x.append(gamma.t+dt)
        y.append(gamma.integrate(gamma.t+dt))
    gamma = UnivariateSpline(x, y, k=3, s=0)

    print "Total volume in gamma_diff calc: " + str(UnivariateSpline(r, dVdr(r), k=3, s=0).integral(0., a))
    print "Total nbi source: " + str(UnivariateSpline(r, part_src_nbi_totint(r) * dVdr(r), k=3, s=0).integral(0., 1.))
    return gamma(r)