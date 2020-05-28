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

def calc_gamma_int_method(r, part_src_nbi_tot, part_src_nbi_lost, izn_rate, iol_adjusted=False, F_orb=None):
    # Piper Changes: Added cylindrical integral method as a separate function. This will be set to a separate variable in the main code.
    gamma = np.zeros(r.shape)

    # Boundary condition at magnetic axis. Needs to be in units of ions/m^3.
    # Only has the second term, since it's the center value. Also uses delta_r of the next point to avoid indexing into a non-existant location.
    # If not adjusted for IOL, part_src_nbi_lost = 0 anyway, so no need for an IOL check.
    gamma[0] = (part_src_nbi_tot[0] - 2*part_src_nbi_lost[0] + izn_rate[0])*(r[1] - r[0])
    gamma[1] = gamma[0] + (part_src_nbi_tot[1] - 2*part_src_nbi_lost[1] + izn_rate[1])*(r[1] - r[0])

    # You'll  prolly want to change this since it uses a dreaded for loop.
    for n in range(2,len(r)): # Pretty sure len() is still valid for multidimensional arrays.
        # The 2*part_src_nbi_lost is the factor of 2 in the fast IOL.
        if iol_adjusted:
            # Imported the exp() function for the thermal IOL attenuation.
            gamma[n] = (r[n-1]/r[n]) * gamma[n-1] * exp(-2*(F_orb[n]-F_orb[n-1])) + (part_src_nbi_tot[n] - 2*part_src_nbi_lost[n] + izn_rate[n])*(r[n] - r[n-1])
        else:
            gamma[n] = (r[n-1]/r[n]) * gamma[n-1] + (part_src_nbi_tot[n] + izn_rate[n])*(r[n] - r[n-1])

    return gamma