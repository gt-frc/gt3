#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from scipy.interpolate import UnivariateSpline
from math import sqrt, pi, exp
from GT3.RadialTransport.Functions.CalcCoulLog import calc_coul_log
from scipy import constants
import numpy as np
import matplotlib.pyplot as plt

eps_0 = constants.epsilon_0
e = constants.elementary_charge
m_e = constants.electron_mass
m_d = constants.physical_constants['deuteron mass'][0]
m_c = 12 / constants.N_A / 1E3  # in kg


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

def calc_Qe_int_method(r, n, T, en_src_nbi_e_tot, cool_rate):  # Piper Changes: Same as Qi changes.

    Qe = np.zeros(r.shape)
    qie = calc_qie(n, T, ion_species='D')

    Qe[0] = (en_src_nbi_e_tot[0] - qie[0] - cool_rate[0]) * (r[1] - r[0])
    Qe[1] = Qe[0] + (en_src_nbi_e_tot[1] - qie[1] - cool_rate[1]) * (r[1] - r[0])

    # Integral cylindrical form of the energy balance equation.
    # Identical in form to the continuity equation, but different source term.
    for n in range(2, len(r)):
        Qe[n] = (r[n - 1] / r[n]) * Qe[n - 1] + (en_src_nbi_e_tot[n] - qie[n] - cool_rate[n]) * (r[n] - r[n - 1])

    return Qe


def calc_Qi_diff_method(r, a, cxcool, Qie, en_src_nbi_tot, en_src_nbi_lost, dVdrho, iol_adjusted=False, E_orb=None, verbose=False):
    dE_orb = UnivariateSpline(r, E_orb, k=3, s=0).derivative()
    en_src_nbi_totint = UnivariateSpline(r, en_src_nbi_tot, k=3, s=0)
    en_src_nbi_lostint = UnivariateSpline(r, en_src_nbi_lost, k=3, s=0)
    Qie_int = UnivariateSpline(r, Qie, k=3, s=0)
    cxcoolint = UnivariateSpline(r, cxcool, k=3, s=0)
    iolPeak = np.where(dE_orb(r) == dE_orb(r).max())

    def f(t, flux, cxcool, Qie, Q_i_nbi, Qnbi_loss, dEdr, iolFlag, peak):
        S = Q_i_nbi(t) - cxcool(t) + Qie(t)
        # Physically, if the IOL peak has occured, everything radially outward should be dFdr = 0.0 since F(r)
        # should equal 0.5 until r=1.0.66
        dEdrval = dEdr(t)
        if t >= peak:
            dEdrval = 0.0
        if iolFlag:
            return S - Qnbi_loss(t) - (flux * (dEdrval + 1 ) / (t + 0.003))
        else:
            return S - (flux * (dEdr(t) + 1 ) / (t + 0.003))

    from scipy.integrate import ode

    flux = ode(f).set_integrator('vode', with_jacobian=False)
    flux.set_initial_value(0., 0.).set_f_params(cxcoolint, Qie_int, en_src_nbi_totint, en_src_nbi_lostint, dE_orb,
                                                iol_adjusted, r[iolPeak])
    dt = a / len(r)
    x, y = [], []
    while flux.successful() and flux.t < a:
        x.append(flux.t + dt)
        y.append(flux.integrate(flux.t + dt))
    flux = UnivariateSpline(x, y, k=3, s=0)

    if verbose:
        plot = plt.figure()
        fig1 = plot.add_subplot(311)
        fig2 = plot.add_subplot(312)
        fig3 = plot.add_subplot(313)

        fig1.scatter(r, Qie_int(r), color="green")
        fig1.scatter(r, cxcoolint(r), color="yellow")
        fig1.scatter(r, en_src_nbi_totint(r), color="red")
        fig1.scatter(r, en_src_nbi_lostint(r), color="black")
        fig1.legend([r"$Q_{ie}$", r"$Q_{cx}$", r"$Q_{nbi,kept}$", r"$Q_{nbi,lost}$"])
        fig1.set_xlim(0.85 * a, a)

        fig2.scatter(r, dE_orb(r), color="red")
        fig2.set_xlim(0.85 * a, a)

        fig3.scatter(r, flux(r))
        fig3.set_xlim(0.85 * a, a)
        plt.show()

    return flux(r)

def calc_Qi_int_method(r, n, T, qie, en_src_nbi_i_kept, cool_rate, iol_adjusted=False, E_orb=None):  # formerly qheat

    Qi = np.zeros(r.shape)

    # Boundary condition at the magnetic axis.
    # Only has the second term, since it's the center value. Also uses delta_r of the next point.
    # If not adjusted for IOL, en_src_nbi_kept = en_src_nbi_tot, so no need for an IOL check.
    Qi[0] = (en_src_nbi_i_kept[0] - cool_rate[0] + qie[0]) * (r[1] - r[0])
    Qi[1] = Qi[0] + (en_src_nbi_i_kept[1] - cool_rate[1] + qie[1]) * (r[1] - r[0])

    # Integral cylindrical form of the energy balance equation.
    # Identical in form to the particle continuity equation, but different source term.
    for n in range(2, len(r)):
        if iol_adjusted:
            # Imported the exp() function for the thermal IOL attenuation.
            Qi[n] = (r[n - 1] / r[n]) * Qi[n - 1] * exp(-(E_orb[n] - E_orb[n - 1])) + (
                        en_src_nbi_i_kept[n] - cool_rate[n] + qie[n]) * (r[n] - r[n - 1])
        else:
            Qi[n] = (r[n - 1] / r[n]) * Qi[n - 1] + (en_src_nbi_i_kept[n] - cool_rate[n] - qie[n]) * (r[n] - r[n - 1])

    return Qi


def calc_qie(n, T, ion_species='D'):
    """Calculates collisional energy transfer between an ion species and electrons
    Reference Equation 4.90 in Stacey's Fusion Plasma Physics Book

    :param n:
    :param T:
    :param ion_species:
    :return:
    """

    if ion_species == 'C':
        zi = 1
        Ti = T.i.J  # assumed carbon is at the same temperature as main ion species
        mi = m_c
    else:
        zi = 1
        Ti = T.i.J
        mi = m_d

    coul_log = calc_coul_log(zi, 1, T.i.J, n.i)  # TODO: check these inputs

    qie = n.e * (zi * e**2)**2 * m_e * coul_log * (1 - Ti / T.e.J) / \
          (2*pi * eps_0**2 * np.sqrt(2*pi*m_e*T.e.J) * mi * (1 + 4*sqrt(pi)/3 * (3*m_e*Ti/(2*mi*T.e.J))**1.5))

    return qie * n.i #multiply this shit by density
