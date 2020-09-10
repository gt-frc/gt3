#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from scipy.interpolate import UnivariateSpline, interp1d
import numpy as np
from math import exp
import matplotlib.pyplot as plt

def calc_gamma_diff_method(r, a, part_src_nbi_tot, part_src_nbi_lost, izn_rate, dVdrho, iol_adjusted=False, F_orb=None,
                           verbose=False):
    dVdr = UnivariateSpline(r, dVdrho(r / a) / a, k=2, s=0)
    dF_orb = UnivariateSpline(r, F_orb, k=3, s=0).derivative()
    izn_rateint = UnivariateSpline(r, izn_rate, k=2, s=0)
    part_src_nbi_totint = UnivariateSpline(r, part_src_nbi_tot, k=2, s=0)
    part_src_nbi_lostint = UnivariateSpline(r, part_src_nbi_lost, k=2, s=0)
    iolPeak = np.where(dF_orb(r) == dF_orb(r).max())

    def f(t, gamma, sion, snbi, snbi_loss, dFdr, iolFlag, peak):
        S = snbi(t) + sion(t)
        # Physically, if the IOL peak has occured, everything radially outward should be dFdr = 0.0 since F(r)
        # should equal 0.5 until r=1.0.66
        dFdrval = dFdr(t)
        if t >= peak:
            dFdrval = 0.0
        if iolFlag:
            return S - snbi_loss(t) - gamma * (dFdrval + 1) / (t + 0.003)
        else:
            return S - gamma * (1 / (t + 0.003))

    from scipy.integrate import ode

    gamma = ode(f).set_integrator('vode', with_jacobian=False)
    gamma.set_initial_value(0., 0.).set_f_params(izn_rateint, part_src_nbi_totint, part_src_nbi_lostint, dF_orb, iol_adjusted, r[iolPeak])
    dt = a / len(r)
    x, y = [], []
    while gamma.successful() and gamma.t < a:
        x.append(gamma.t+dt)
        y.append(gamma.integrate(gamma.t+dt))
    #gamma = UnivariateSpline(x, y, k=3, s=0)
    gamma = interp1d(x, np.array([float(b) for b in y]), kind="linear", fill_value="extrapolate")

    if verbose:
        plot = plt.figure()
        fig1 = plot.add_subplot(311)
        fig2 = plot.add_subplot(312)
        fig3 = plot.add_subplot(313)

        fig1.scatter(r, izn_rateint(r), color="green")
        fig1.scatter(r, part_src_nbi_totint(r), color="red")
        fig1.scatter(r, part_src_nbi_lostint(r), color="black")
        fig1.legend([r"$S_{ion}$", r"$S_{nbi,tot}$", r"$S_{nbi,lost}$"])
        fig1.set_xlim(0.85 * a, a)

        fig2.scatter(r, dF_orb(r), color="red")
        fig2.set_xlim(0.85 * a, a)

        fig3.scatter(r, gamma(r))
        fig3.set_xlim(0.85 * a, a)
        plt.show()

        return fig1, fig2, fig3

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
            gamma[n] = (r[n-1]/r[n]) * gamma[n-1] * exp(-2*(F_orb[n]-F_orb[n-1])) + (part_src_nbi_tot[n] - 2.*part_src_nbi_lost[n] + izn_rate[n])*(r[n] - r[n-1])
        else:
            gamma[n] = (r[n-1]/r[n]) * gamma[n-1] + (part_src_nbi_tot[n] + izn_rate[n])*(r[n] - r[n-1])

    return gamma