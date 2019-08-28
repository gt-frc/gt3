#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 08:58:46 2018

@author: max
"""
from __future__ import division
import numpy as np
from math import pi, cos, sqrt, exp
import sys
import matplotlib.pyplot as plt
from collections import namedtuple
from scipy import constants
from scipy.integrate import odeint
from scipy.interpolate import interp1d, UnivariateSpline
from core import calc_fsa

eps_0 = constants.epsilon_0
e = constants.elementary_charge
m_e = constants.electron_mass
m_d = constants.physical_constants['deuteron mass'][0]
m_c = 12 / constants.N_A / 1E3  # in kg

z_d = 1  # atomic number of deuterium
z_c = 6  # atomic number of carbon
ch_d = e * z_d  # charge of deuterium
ch_c = e * z_c  # charge of carbon

E_phi = 0.04  # toroidal electrostatic potential

def neutPatch(core):
    """

    Re-runs FSA of quantities needed in rad transport that are not available because neutrals update is run after
    these values, which depend on neutral data, are calculated in core
    :param core:
    :return:
    """
    core.izn_rate_fsa_s = calc_fsa(core.izn_rate.s, core.R, core.Z)
    core.izn_rate_fsa_t = calc_fsa(core.izn_rate.t, core.R, core.Z)
    core.izn_rate_fsa = calc_fsa(core.izn_rate.s + core.izn_rate.t, core.R, core.Z)

    #Piper chages: Use carbon Lz for the cool_rate calculation
    core.cool_rate_fsa = calc_fsa(core.n.e * core.n.C * np.nan_to_num(core.Lz_C.s) + \
    +                            core.n.e * core.n.C * np.nan_to_num(core.Lz_C.t), \
                                  core.R, core.Z)

    """
    Calculate the neutrals now
    """

    core.dn_dr_fsa = calc_fsa(core.dn_dr, core.R, core.Z)

    """
    Cuts off neutrals calculations at 0.8 to give realistic neutral shit
    """

    core.izn_rate_fsa = np.array(map(lambda x: core.izn_rate_fsa[-1] * x**10, core.r[:,0]/core.a))
    core.cool_rate_fsa = np.array(map(lambda x: core.cool_rate_fsa[-1] * x**10, core.r[:,0]/core.a))
    core.dn_dr_fsa = np.array(map(lambda x: core.dn_dr_fsa[-1] * x**10, core.r[:,0]/core.a))


    core.n_fsa = namedtuple('n', 'i e n C')(
        calc_fsa(core.n.i, core.R, core.Z),
        calc_fsa(core.n.e, core.R, core.Z),
        namedtuple('n', 's t tot')(
            np.array(map(lambda x: core.n_fsa.n.s[-1] * x ** 10, core.r[:, 0] / core.a))/1E3,  # slow
            np.array(map(lambda x: core.n_fsa.n.t[-1] * x ** 10, core.r[:, 0] / core.a))/1E3,  # thermal
            np.array(map(lambda x: core.n_fsa.n.tot[-1] * x ** 10, core.r[:, 0] / core.a))/1E3),  # total
        calc_fsa(core.n.C, core.R, core.Z))


def corePatch(core, neutFlag=True):
    """
    Updates deuterium ion density and zeff. D density is invalid if no D density file is given because it will be
    set to 0. This screws up subsequent calculations in non-obvious ways. As an example, z_eff calculation will
    include D density as 0s but still have a carbon density from the fracz input file, giving an invalid result

    This is addressed by recreating the n_fsa namedtuple and z_eff_fsa in full, as namedtuples cannot be changed piecemeal

    :param core:
    :return:
    """


    if neutFlag:
        core.n = namedtuple('n', 'i e n C')(core.n.e/(1.+.025*6.0), core.n.e, core.n.n, 0.025 * core.n.e/(1.+.025*6.0))   # TODO: Update 0.025 and 6.0 to actual fracz and zbar2
    else:
        core.n = namedtuple('n', 'i e n C')(core.n.e / (1. + .025 * 6.0), core.n.e,
                                   namedtuple('n', 's t tot')(
                                       np.zeros(core.n.i.shape),  # slow
                                       np.zeros(core.n.i.shape),  # thermal
                                       np.zeros(core.n.i.shape)),  # total
                                   0.025 * core.n.e / (1. + .025 * 6.0))
    core.n_fsa = namedtuple('n', 'i e n C')(
        calc_fsa(core.n.i, core.R, core.Z),
        calc_fsa(core.n.e, core.R, core.Z),
        namedtuple('n', 's t tot')(
            calc_fsa(core.n.n.s, core.R, core.Z),  # slow
            calc_fsa(core.n.n.t, core.R, core.Z),  # thermal
            calc_fsa(core.n.n.tot, core.R, core.Z)  # total
        ),
        calc_fsa(core.n.C, core.R, core.Z))

    core.z_eff_fsa = calc_fsa((core.n.i * (1.**2) + core.n.C * (6.0**2))/(core.n.i * 1.0 + core.n.C * 6.0), core.R, core.Z) #TODO Similar updates (1.0 = atnum, 6.0 = zbar2)


def calc_part_src_nbi(beam, iol_adjusted=False, F_orb_nbi=None):
    """
    """
    # # old suspicious calculation
    # snbi = (0.624E25 * beam.dPdV.v1D.W / (vol * beam.E.J)) * beam.num * beam.dp

    # Piper Changes: Changed function to calculate particle source and source density.
    # tot, lost, and kept are all densities. The first return variable is a total source.
    # particle source needs to be split into a total, amount kept, and amount lost to calculate the return current.

    # These equations assume:
    # beam.dPdV.v1D.W = H(rho) * Pbeam * pwrfrac / Volp_nbi
    # beam.dPdr.v1D.W = beam.dPdV.v1D.W * dV

    # The factor of 2 can be applied to F_orb in the continuity equation later.

    # nbi source in particles/sec. Needed to calculate gamma with sources.
    part_src_nbi = beam.dPdr.v1D.W / beam.E.J

    # nbi source in # of particles/(m^3 * sec). Needed to calculate gamma with source densities.
    part_src_nbi_tot = beam.dPdV.v1D.W / beam.E.J

    part_src_nbi = beam.dPdV.v1D.W / beam.E.J #* beam.dp
    if iol_adjusted:
        part_src_nbi = part_src_nbi * (1 - F_orb_nbi)
        part_src_nbi_lost = part_src_nbi_tot * (F_orb_nbi)
        part_src_nbi_kept = part_src_nbi_tot * (1 - F_orb_nbi)
    else:
        part_src_nbi_lost = np.zeros(part_src_nbi_tot.shape)
        part_src_nbi_kept = part_src_nbi_tot

    return part_src_nbi, part_src_nbi_tot, part_src_nbi_lost, part_src_nbi_kept


def calc_en_src_nbi(beam, iol_adjusted=False, E_orb_nbi=None):
    """
    """
    
    # Piper Changes: verify this calculation. Assumes the same things as the particle source eq.
    
    # nbi energy source in Joules/sec. Surprisingly, no need to actually calculate anything here.
    en_src_nbi = beam.dPdr.v1D.W
    
    # nbi energy source in Joules/(m^3 * sec)
    en_src_nbi_tot = beam.dPdV.v1D.W
    
    if iol_adjusted:
        en_src_nbi = en_src_nbi * (1 - E_orb_nbi)
        en_src_nbi_lost = en_src_nbi_tot * (E_orb_nbi)
        en_src_nbi_kept = en_src_nbi_tot * (1 - E_orb_nbi)
    else:
        en_src_nbi_lost = np.zeros(en_src_nbi_tot.shape)
        en_src_nbi_kept = en_src_nbi_tot

    return en_src_nbi, en_src_nbi_tot, en_src_nbi_lost, en_src_nbi_kept


def calc_mom_src_nbi(beam,  n, z_eff, R0_a, fforb):
    """Calculates toroidal momentum input from a neutral beam (or beam component)

    :param beam:
    :param n:
    :param z_eff:
    :param rpts:
    :param R0_a:
    :return:
    """

    def calc_atten(beam, n, z_eff, rpts=200):  # TODO: Figure out what this is. has to do with beam attenuation
        delma = 1 / (rpts - 1)
        alphain = 0.6475  # TODO: no idea what this is
        xlam = 5.5E17 * (beam.E.J / 1) / (2.0 * 0.5 * n.i * z_eff ** 0.5)
        atten = 1 - np.exp(-delma / cos(alphain) / xlam)
        return atten

    atten = calc_atten(beam, n, z_eff)
    unatten = 1-atten  # might need to flip atten...
    torque = calc_torque(beam,fforb)

    tor_mom = unatten*atten*torque/R0_a
    return tor_mom


def calc_torque(beam,fforb):
    """ Calculates torque from a neutral beam (or beam component)

    torque = F * r_tan = (P/v) * r_tan = (P/sqrt(2E/m)) * r_tan = P * sqrt(m/(2E)) * r_tan

    :param beam: beam object with attributes z, m, a, en, pwr, rtan
    :return: torque
    """

    power = beam.P.W
    energy = beam.E.J
    mass = beam.m
    rtan = beam.rtan

    torque = power * np.sqrt(0.5 * mass / energy) * rtan * (1.0-fforb) # Piper Changes: Included fast ion losses.
    return torque


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

    from scipy.integrate import ode

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


def calc_Qi_int_method(r, n, T, en_src_nbi_i_kept, cool_rate, iol_adjusted=False, E_orb=None):  # formerly qheat

    Qi = np.zeros(r.shape)

    qie = calc_qie(n, T, ion_species='D')

    # Boundary condition at the magnetic axis.
    # Only has the second term, since it's the center value. Also uses delta_r of the next point.
    # If not adjusted for IOL, en_src_nbi_kept = en_src_nbi_tot, so no need for an IOL check.
    Qi[0] = (en_src_nbi_i_kept[0] - cool_rate[0] - qie[0]) * (r[1] - r[0])
    Qi[1] = Qi[0] + (en_src_nbi_i_kept[1] - cool_rate[1] - qie[1]) * (r[1] - r[0])

    # Integral cylindrical form of the energy balance equation.
    # Identical in form to the particle continuity equation, but different source term.
    for n in range(2, len(r)):
        if iol_adjusted:
            # Imported the exp() function for the thermal IOL attenuation.
            Qi[n] = (r[n - 1] / r[n]) * Qi[n - 1] * exp(-(E_orb[n] - E_orb[n - 1])) + (
                        en_src_nbi_i_kept[n] - cool_rate[n] - qie[n]) * (r[n] - r[n - 1])
        else:
            Qi[n] = (r[n - 1] / r[n]) * Qi[n - 1] + (en_src_nbi_i_kept[n] - cool_rate[n] - qie[n]) * (r[n] - r[n - 1])

    return Qi, qie



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

    #print "Total volume in Qi_diff calc: " + str(UnivariateSpline(r, dVdr(r), k=3, s=0).integral(0., a))
    #print "Total nbi ion energy: " + str(UnivariateSpline(r, (en_src_nbi_keptint(r) + en_src_nbi_lostint(r)) * dVdr(r), k=3, s=0).integral(0., 1.)/(1E6))+" MW"
    return flux(r)
    
def calc_Qe_int_method(r, n, T, en_src_nbi_e_tot, cool_rate): # Piper Changes: Same as Qi changes.

    Qe = np.zeros(r.shape)
    qie = calc_qie(n,T,ion_species='D')
    
    Qe[0] = (en_src_nbi_e_tot[0] - qie[0] - cool_rate[0])*(r[1] - r[0])
    Qe[1] = Qe[0] + (en_src_nbi_e_tot[1] - qie[1] - cool_rate[1])*(r[1] - r[0])
    
    # Integral cylindrical form of the energy balance equation.
    # Identical in form to the continuity equation, but different source term.
    for n in range(2,len(r)):
        Qe[n] = (r[n-1]/r[n]) * Qe[n-1] + (en_src_nbi_e_tot[n] - qie[n] - cool_rate[n])*(r[n] - r[n-1])

    return Qe


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

def calc_cxcool(core, n, T):
    """Calculates charge exchange cooling with the slow neutrals"""
    slowNeutDens = np.array(map(lambda x: n.n.s[-1] * x ** 15, core.r[:, 0] / core.a))
    totalNeutDens = np.array(map(lambda x: n.n.tot[-1] * x ** 15, core.r[:, 0] / core.a))
    #result = 1.5 * n.i * T.i.J * slowNeutDens* ((core.sv.el.st[:,0] + core.sv.el.st[:,0]* (totalNeutDens/n.i)) + core.sv.cx.st[:,0])
    result = 1.5 * n.i * T.i.J * slowNeutDens * ((core.sv.el.st[:, 0] + core.sv.cx.st[:, 0])) / 1E4 # TODO: Work this out for n-n elastic scattering
    return result

def calc_coul_log(z1, z2, T_J, n2):
    """Calculate the coulomb logarithm.

    Reference Equation 1.36 in Dr. Stacey's Fusion Plasma Physics book

    :param z1:
    :param z2:
    :param T_J:
    :param n2:
    :return:
    """

    coul_log = np.log(12 * pi * np.sqrt((eps_0 * T_J)**3 / (n2 * (z2*e)**4 * (z1*e)**2)))
    return coul_log

def calc_coul_log_j_k(z_j, z_k, T_j, n_k):
    # Coulomb logarithm calculation in GTEDGE PARAM subroutine.
    coul_log = np.log(12 * pi * (T_j**1.5) * ((eps_0/e)**1.5) / (np.sqrt(n_k) * (z_k**2.0) * z_j))
    return coul_log

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



def calc_return_cur(r, part_src_nbi_lost, gamma, izn_rate, ch_d, iol_adjusted=False, F_orb=None):

    # Piper Changes: This calculates the return current from fast and thermal IOL.
    # This is a new quantity that i want to calculate. Big part of my thesis.
    
    if iol_adjusted:
        diff_F_orb = UnivariateSpline(r, F_orb, k=1, s=0).derivative()(r) # Need the derivative of F_orb with respect to r.
    else:
        diff_F_orb = np.zeros(r.shape)
    
    Jr_iol = np.zeros(r.shape)
    Jr_iol[0] = (part_src_nbi_lost[0] + (gamma[0] * diff_F_orb[0])) * (r[1]-r[0]) * ch_d
    Jr_iol[1] = Jr_iol[0] + (part_src_nbi_lost[1] + (gamma[1] * diff_F_orb[1]))*(r[1] - r[0]) * ch_d
    
    # Integral cylindrical form of the continuity of (current?) equation.
    for n in range(2,len(r)):
        Jr_iol[n] = (r[n-1]/r[n]) * Jr_iol[n-1] + (part_src_nbi_lost[n] + (gamma[n] * diff_F_orb[n]))*(r[n] - r[n-1]) * ch_d # No factor of 2, because only 1 direction.

    return Jr_iol

# Calculates the radial electric field including the J x B force from IOL
def calc_Er_iol(n_i, n_e, m_i, n_n, Bphi, Lp_i, dp_dr, e_i, T_i, dn_dr, izn_rate, Jr_iol):
    Jr_visc = 0.0
    Jr_neut = -Jr_iol - Jr_visc
    iol_term = (Jr_neut * Bphi**2.0) / (m_i * izn_rate) # Both the ion and neutral densities are built into izn_rate
    diamag_term = -1.0 * (Lp_i * T_i.J) / e_i # Pressure gradient is negative, so multiply by -1
    diamag_term_orig = -1.0 * (dp_dr / (e_i * n_i)) # Original form of diamagnetic term
    neut_dens_term = -1.0 * (T_i.J * dn_dr) / (e_i * n_n.t)
    E_r_iol = iol_term + diamag_term + neut_dens_term
    return E_r_iol, iol_term, diamag_term, diamag_term_orig, neut_dens_term


def calc_mbal_rhs(mom_src_ext, z, n, B_p, gamma):
    """ The function formerly known as y11 and y22
    """

    mbal_rhs = mom_src_ext + (z * e) * (n * E_phi + B_p * gamma)

    return mbal_rhs


def calc_intrin_rot(M_orb, T_J, m):

    intrin_rot = 2 / sqrt(pi) * M_orb * np.sqrt(2 * T_J / m)
    return intrin_rot


def calc_vtor_d_pert(vtor_C_total, vtor_C_intrin, vtor_D_intrin, mom_src_ext, z, n, T, B_p, gamma_D, gamma_C):
    """
    """
    # Piper Changes: Changed n to n.i and n.C.
    # Changed z for mbal_rhs_C to 6.
    mbal_rhs_D = calc_mbal_rhs(mom_src_ext, z, n.i, B_p, gamma_D)
    mbal_rhs_C = calc_mbal_rhs(mom_src_ext, z, n.C, B_p, gamma_C) # TODO: This z should be zbar.

    nu_c_DC = 1 / calc_t90(m_d, m_c, 1, 6, n.C, T.i.J)

    #vtor_C_total = vtor_fluid_C + vtor_intrin_C
    del_v0 = vtor_D_intrin - vtor_C_intrin

    nu_drag_approx = (mbal_rhs_D + mbal_rhs_C) / ((n.i * m_d + n.C * m_c) * vtor_C_total + n.i * m_d * del_v0)

    del_v1 = (mbal_rhs_D - n.i * m_d * nu_drag_approx * vtor_C_total) / (n.i * m_d * (nu_c_DC + nu_drag_approx))

    vtorDPert = vtor_C_total + del_v1
    return vtorDPert


def calc_nu_drag(n_j, m_j, v_tor_j, v_tor_k, mbal_rhs, nu_c):
    nu_drag = (mbal_rhs + (n_j * m_j * nu_c * v_tor_k)) / (v_tor_j * n_j * m_j) - nu_c # Piper Changes: Multiplied V_tor_k in the numerator by n_j*m_j. The algebra was wrong.
    return nu_drag

def calc_vpol(Er, vphi_j, Lp, T, n, z_d, B_t, B_p, vphi_k, vpol_k, z_k):
    vpol = (1.0/B_t) * (1.0/(e*z_d) * -Lp.i * T.i.J + vphi_j * B_p - Er)
    vpol_assum = vpol_k - (T.C.J / (e * B_t)) * (Lp.i - (Lp.C / z_c)) + (B_p / B_t) * (vphi_j - vphi_k)
    vpol_alt = vpol_k - 1.0 / (e * B_t) * (T.i.J * Lp.i - (T.C.J * Lp.C / z_c)) + (B_p / B_t) * (vphi_j - vphi_k)
    return vpol, vpol_assum, vpol_alt


# # TODO: these equations don't match Equations 15.8, 15.9, etc.
# # def calc_xnudrag_old(n):
# #     """
# #     """
# #     # PROBLEM: When nudrag is negative, it can come close to xnuc12 in magnitude
# #     #          and blow up the pertrubation theory.
# #     torv = vtorChat + intrin_C
# #     delv0 = intrin_D - intrin_C
# #
# #     nu_drag_eff1 = (y11 + y22) / \
# #                   ((n.i * m_D + n.C * m_C) * torv +  n.i * m_D * delv0)
# #
# #     delv1 = (y11 - n.i * m_d * xnudrageff1 * torv) / \
# #             (n.i * m_D * (nu_c_DC + nu_drag_eff1))
# #
# #     xnudrageff2 = (y22 + n.C * m_C * nu_c_CD * delv1) / \
# #                   (n.C * m_C * torv)
# #
# #     nu_drag_tot1 = (y11 + y22 - n.i * m_D * xnudrageff1 * delv1) / \
# #                ((n.i * m_D + n.C * m_C) * torv)
# #
# #     delv2 = (y11 - n.i * m_D * xnudtot1 * torv) / \
# #             (n.i * m_D * (nu_c_DC * xnudtot1))
# #
# #     xnudrag_D = (y11 - n.i * m_D * nu_c_DC * delv2) / \
# #                      (n.i * m_D * (torv + delv2))
# #
# #     xnudrag_C = (y22 + n.C * m_C * nu_c_CD * delv2) / \
# #                      (n.C * m_C * torv)
# #
# #     return xnudrag_D, xnudrag_C


def calc_t90(m1, m2, z1, z2, n2, T_J):
    """

    :param m1:
    :param m2:
    :param z1:
    :param z2:
    :param n2:
    :param T_J:
    :return:
    """
    m_r = calc_reduced_mass(m1, m2)
    coul_log = calc_coul_log(z1, z2, T_J, n2)
    t_90 = 2*pi*sqrt(m_r) * eps_0**2 * (3*T_J)**1.5 / \
           (n2 * ((z1*e) * (z2 * e))**2 * coul_log) # Piper Changes: Put the denominator in parenthises.
    return t_90

def calc_nu_j_k(m_j, m_k, z_j, z_k, T_j, n_k):
    m_r = calc_reduced_mass(m_j, m_k)
    coul_log_j_k = calc_coul_log_j_k(z_j, z_k, T_j, n_k)
    C1 = 1/((((4.8E-10)/(1.6E-12))**1.5)*((4.8E-10)**2.5))
    nu_j_k = 3.34*(coul_log_j_k * (z_j**2) * (z_k**2) * 1E-6 * n_k) / \
            (C1*np.sqrt(m_r*1E3)*(T_j**1.5))
    #nu_j_k = (1.2734E14) * (z_j**2) * (z_k**2) * n_k * coul_log_j_k / \
    #        (np.sqrt(m_r * 1E3) * (T_j**1.5))
    return nu_j_k


def calc_reduced_mass(m1, m2):
    return m1 * (1 + (m1 / m2))  # * 1E3  # TODO: not sure why this 1E3 is here




def calc_vr_pinch():
    vr_pinch = 0  # TODO: Add the equation for vr_pinch
    return vr_pinch


def calc_nustar(nu_c, q95, R0_a, vpol):
    nustar = nu_c * abs(q95) * R0_a / vpol
    return nustar


# def calc_gamma_diff_method(r, a, part_src_nbi_tot, part_src_nbi_lost, izn_rate, dVdrho, iol_adjusted=False, F_orb=None):
#     dVdr = UnivariateSpline(r, dVdrho(r / a) / a, k=3, s=0)
#     dF_orb = UnivariateSpline(r, F_orb, k=3, s=0).derivative()
#     izn_rateint = UnivariateSpline(r, izn_rate, k=3, s=0)
#     part_src_nbi_totint = UnivariateSpline(r, part_src_nbi_tot, k=3, s=0)
#     part_src_nbi_lostint = UnivariateSpline(r, part_src_nbi_lost, k=3, s=0)
#
#     def f(t, gamma, sion, snbi, snbi_loss, dFdr, iolFlag):
#         S = snbi(t) + sion(t)
#         if iolFlag:
#             return S - snbi_loss(t) - gamma * (dFdr(t) + 1) / (t + 0.003)
#         else:
#             return S - gamma * (1 / (t + 0.003))
#
#     from scipy.integrate import ode
#
#     gamma = ode(f).set_integrator('vode', with_jacobian=False)
#     gamma.set_initial_value(0., 0.).set_f_params(izn_rateint, part_src_nbi_totint, part_src_nbi_lostint, dF_orb, iol_adjusted)
#     dt = a / len(r)
#     x, y = [], []
#     while gamma.successful() and gamma.t < a:
#         x.append(gamma.t+dt)
#         y.append(gamma.integrate(gamma.t+dt))
#     gamma = UnivariateSpline(x, y, k=3, s=0)
#
#     print "Total volume in gamma_diff calc: " + str(UnivariateSpline(r, dVdr(r), k=3, s=0).integral(0., a))
#     print "Total nbi source: " + str(UnivariateSpline(r, part_src_nbi_totint(r) * dVdr(r), k=3, s=0).integral(0., 1.))
#     return gamma(r)

def calc_external_term(M_phi, n_j, ch_j, B_p):
    ext_term = (-M_phi - (n_j * ch_j * E_phi)) / (n_j * ch_j * B_p)
    return ext_term

def calc_poloidal_term(n_j, m_j, ch_j, nu_jk, nu_dj, B_t, B_p, v_pol_j):
    pol_term = (n_j * m_j * (nu_jk + nu_dj) * B_t * v_pol_j) / (n_j * ch_j * (B_p**2.0))
    return pol_term
    
def calc_radial_E_field_term(n_j, m_j, ch_j, nu_jk, nu_dj, Er, B_p):
    Er_term = (n_j * m_j * (nu_jk + nu_dj) * Er) / (n_j * ch_j * (B_p**2.0))
    return Er_term
    
def calc_toroidal_term(n_j, m_j, ch_j, nu_jk, B_p, v_tor_k):
    tor_term = (-n_j * m_j * nu_jk * v_tor_k) / (n_j * ch_j * B_p)
    return tor_term
    
def calc_pinch_velocity(ext_term, pol_term, Er_term, tor_term):
    vr_pinch = ext_term + pol_term + Er_term + tor_term
    return vr_pinch


def calc_nustar(nu_c, q95, R0_a, vpol):
    nustar = nu_c * abs(q95) * R0_a / vpol
    return nustar
def calc_Er_mom_bal(n, e_z, dp_dr, T, Lp, v_tor, v_pol, B_tor, B_pol):
    pres_term = -1.0 / (n * e_z) * dp_dr
    #pres_term_simp = -1.0 * (Lp * T.J) / e_z # alternative method using grad scale length.
    vxb_term = -1.0 * (v_pol * B_tor - v_tor * B_pol)
    E_r = pres_term + vxb_term
    #E_r_simp = pres_term_simp + vxb_term
    return E_r, pres_term, vxb_term


# def calc_xpartdot_old(beam, n, z_eff, rpts):  # TODO: wtf is xpartdot?
#     atten = calc_atten(beam, n, z_eff, rpts)
#     unatten = 1-atten  # might need to flip atten...
#     volm =
#     xpartdot = unatten * atten * beam.P * 1E6 / \
#     (beam.E.J / comp * 1E3 * e) / volm
#     return xpartdot


# def calc_qnbi_old(beam, n, z_eff, rpts):
#     xpartdot1 = calc_xpartdot(beam.e, beam.Pwr, pwr_frac, 1, beam.a, n, z_eff, rpts)
#     xpartdot2 = calc_xpartdot(beam.e, beam.Pwr, pwr_frac, 2, beam.a, n, z_eff, rpts)
#     xpartdot3 = calc_xpartdot(beam.e, beam.Pwr, pwr_frac, 3, beam.a, n, z_eff, rpts)
#     qnbi = (xpartdot1 + xpartdot2/2 + xpartdot3/3) * 1E3 * e * e_beam
#     return qnbi


class Chi:
    """
    Chi class claculates various chis and provides source terms as necessary
    """
    def __init__(self, data, core, n, T, L, nustar, reInterp=False):

        self.Qi = data.Qi_diff
        self.Qe = data.Qe_diff
        self.conv15 = UnivariateSpline(core.r[:,0], .5 * ch_d * data.gamma_diff_D * T.i.ev, k=3, s=0)(core.r[:,0])
        self.conv25 = UnivariateSpline(core.r[:,0], .5 * ch_d * data.gamma_diff_D * T.i.ev, k=3, s=0)(core.r[:,0])
        self.heatvisc = self.viscCalc(data, core, n, T)
        #self.heatvisc = np.zeros(self.conv25.shape)
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

    def viscCalc(self, data, core, n, T):
        f1 = 1 # Must determine later, appears to be geometeric factor in Eq(4) in POP052504(2010)
        fp = core.B_p_fsa/core.B_t_fsa
        vth = [sqrt(2. * a  * ch_d / m_d) for a in T.i.ev]
        eta0 = [a * m_d * b * c * core.R0_a * f1 for a,b,c in zip(n.i, vth, core.q[:,0])]
        eta4 = [a * m_d * c * ch_d / (ch_d * abs(b)) for a, b, c in zip(n.i, core.B_t_fsa, T.i.ev)]
        vrad = data.gamma_diff_D / n.i

        # a = vtor    b = fp   c = eta0
        # d = vrad    f = vthet g = eta 4

        return [a * (b * c * d - .5 * g * (4.0 * a + f)) - .5 * f * (c * d + g * (a + .5 * f)) for a,b,c,d,f,g in zip(data.vtor_D_total, fp, eta0, vrad, data.vpol_D, eta4)]


class RadialTransport(Chi):
    
    def __init__(self, inp, core, iol, nbi, iolFlag=True, neutFlag=True, debugFlag=False ):
        sys.dont_write_bytecode = True





        ##############################################################
        # prep quantities for 1D transport analysis
        ##############################################################

        #prepare beams object
        beam_D = nbi.beams.D1 # Piper changes: Changed beam object references to match beamdep.
        beam_D2 = nbi.beams.D2
        beam_D3 = nbi.beams.D3
        corePatch(core, neutFlag) # Patch to update values not brought in via ffiles (ni, zeff)
        neutPatch(core)
        dn_dr = core.dn_dr_fsa


        # prepare core and iol quantities
        r = core.r.T[0]    # TODO: Should this be a flux surface average?
        self.rhor = core.r[:,0]/core.a

        if neutFlag:
            izn_rate = core.izn_rate_fsa  # TODO: Should this be a flux surface average or a flux surface total?
        else:
            izn_rate = np.zeros(core.izn_rate_fsa.shape)
        cool_rate = core.cool_rate_fsa  # TODO: Should this be a flux surface average or a flux surface total?
        n = core.n_fsa

        T = core.T_fsa
        L = core.L_fsa
        z_eff = core.z_eff_fsa
        R0_a = core.R0_a
        vol_P = core.vol
        B_p = core.B_p_fsa
        B_t = core.B_t_fsa
        Er = core.E_r_fsa #* 1000.0 # Piper Changes: Convert input Er from kV/m to V/m
        Lp = core.Lp_fsa
        dp_dr = core.dp_dr_fsa

        #prepare iol quantities
        F_orb_d = iol.forb_d_therm_1D
        M_orb_d = iol.morb_d_therm_1D
        E_orb_d = iol.eorb_d_therm_1D

        F_orb_c = iol.forb_c_therm_1D
        M_orb_c = iol.morb_c_therm_1D
        E_orb_c = iol.eorb_c_therm_1D

        F_orb_t = iol.forb_t_therm_1D
        M_orb_t = iol.morb_t_therm_1D
        E_orb_t = iol.eorb_t_therm_1D
        
        # prepare fast iol quantities
        F_orb_d_nbi = iol.forb_d_nbi_1D
        M_orb_d_nbi = iol.morb_d_nbi_1D
        E_orb_d_nbi = iol.eorb_d_nbi_1D

        ##############################################################
        # particle balance
        ##############################################################

        # Piper Changes: Function now has 4 outputs, so we set 4 variables. The first variable is for the differential balance methods, and the other 3 are for the integral methods.
        part_src_nbi_D, part_src_nbi_D_tot, part_src_nbi_D_lost, part_src_nbi_D_kept = calc_part_src_nbi(beam_D, iol_adjusted=iolFlag, F_orb_nbi=F_orb_d_nbi)
        part_src_nbi_D2, part_src_nbi_D2_tot, part_src_nbi_D2_lost, part_src_nbi_D2_kept = calc_part_src_nbi(beam_D2, iol_adjusted=iolFlag, F_orb_nbi=F_orb_d_nbi)
        part_src_nbi_D3, part_src_nbi_D3_tot, part_src_nbi_D3_lost, part_src_nbi_D3_kept = calc_part_src_nbi(beam_D3, iol_adjusted=iolFlag, F_orb_nbi=F_orb_d_nbi)

        # Piper changes: Calculate beam totals for sources, sinks, and totals.
        part_src_nbi = part_src_nbi_D + part_src_nbi_D2 + part_src_nbi_D3
        part_src_nbi_tot = part_src_nbi_D_tot + part_src_nbi_D2_tot + part_src_nbi_D3_tot
        part_src_nbi_lost = part_src_nbi_D_lost + part_src_nbi_D2_lost + part_src_nbi_D3_lost
        part_src_nbi_kept = part_src_nbi_D_kept + part_src_nbi_D2_kept + part_src_nbi_D3_kept

        # Piper changes: Changed names of particle and heat flux so it's easier to tell what method is used.
        self.gamma_diff_D = calc_gamma_diff_method(r, core.a, part_src_nbi, part_src_nbi_lost, izn_rate, core.dVdrho, iol_adjusted=iolFlag, F_orb=F_orb_d) # Differential Cylindrical Method
        self.gamma_int_D = calc_gamma_int_method(r, part_src_nbi_tot, part_src_nbi_lost, izn_rate, iol_adjusted=iolFlag, F_orb=F_orb_d) # Integral Cylindrical Method
        self.gamma_C = np.zeros(self.gamma_int_D.shape)
        
        # Piper changes: Calculate radial return current (Uses integral cylindrical gamma)
        self.jr_iol = calc_return_cur(r, part_src_nbi_lost, self.gamma_int_D, izn_rate, ch_d, iol_adjusted=iolFlag, F_orb=F_orb_d)
        self.Er_iol, self.iol_term, self.diamag_term, self.diamag_term_orig, self.neut_dens_term = calc_Er_iol(n.i, n.e, m_d, n.n, B_t, Lp.i, dp_dr.i, e * z_d, T.i, dn_dr, izn_rate, self.jr_iol)

        ##############################################################
        # momentum balance
        ##############################################################

        # calculate toroidal momentum source rates from beams
        mom_src_nbi_D = calc_mom_src_nbi(beam_D, n, z_eff, R0_a, F_orb_d_nbi)
        mom_src_nbi_D2 = calc_mom_src_nbi(beam_D2, n, z_eff, R0_a, F_orb_d_nbi)
        mom_src_nbi_D3 = calc_mom_src_nbi(beam_D3, n, z_eff, R0_a, F_orb_d_nbi)
        self.mom_src_nbi = (mom_src_nbi_D + mom_src_nbi_D2 + mom_src_nbi_D3) / vol_P

        # calculate momentum source from anomalous torque
        self.mom_src_anom = np.zeros(r.shape)  # TODO: Anomolous torque

        frac = n.i / (n.i + n.C)
        self.mom_src_tor_D_tot = (1 - frac) * (self.mom_src_nbi + self.mom_src_anom)
        self.mom_src_tor_C_tot = frac * (self.mom_src_nbi + self.mom_src_anom)

        ##############################################################
        # rotation
        ##############################################################

        # calculate carbon toroidal rotation
        self.vtor_C_intrin = calc_intrin_rot(M_orb_c, T.i.J, m_c)
        self.vtor_C_total = core.v_1D.tor.C
        self.vtor_C_fluid = self.vtor_C_total - self.vtor_C_intrin

        # calculate deuterium toroidal rotation
        self.vtor_D_intrin = calc_intrin_rot(M_orb_d, T.i.J, m_d)

        # Piper Changes: Changed core.v_1D.tor.C.any() to core.v_1D.tor.D.any(). Carbon velocity should be a given.
        if not core.v_1D.tor.D.any():  # if array is all zeros, then no input. Use perturbation theory.
            self.vtor_D_total = calc_vtor_d_pert(self.vtor_C_total,
                                                 self.vtor_C_intrin,
                                                 self.vtor_D_intrin,
                                                 self.mom_src_tor_D_tot,
                                                 1,
                                                 n,
                                                 T,
                                                 B_p,
                                                 self.gamma_int_D,
                                                 self.gamma_C) # Piper Changes: Uses integral cylindrical gamma
            # Piper changes: added a message to let the user know the D velocity was calculated.
            print 'Deuterium toroidal velocity calculated from perturbation theory.'
        else:
            # Piper Changes: For some reason this used to set D velocity to C velocity,
            # which overwrote the input D velocity.
            self.vtor_D_total = core.v_1D.tor.D

        self.vtor_fluid_D = self.vtor_D_total - self.vtor_D_intrin

        # calculate carbon and deuterium poloidal rotation
        try:
            self.vpol_C = core.v_1D.pol.C
            self.vpol_D, self.vpol_D_assum, self.vpol_D_alt = calc_vpol(Er, self.vtor_D_total, Lp, T, n, z_d, B_t, B_p, self.vtor_C_total, self.vpol_C, z_c)
        except:
            self.vpol_D = self.vpol_C / 0.4
            print 'could not calculate deuterium poloidal rotation'
            pass
        
        # Nick Changes: TEMPORARY - Calculate Er using pressure gradient vs. scale length.
        self.Er_calc_D, self.Er_pres_term_D, self.Er_vxb_term_D = calc_Er_mom_bal(n.i, e * z_d, dp_dr.i, T.i, Lp.i, self.vtor_D_total, self.vpol_D, B_t, B_p)
        self.Er_calc_C, self.Er_pres_term_C, self.Er_vxb_term_C = calc_Er_mom_bal(n.C, e * z_c, dp_dr.C, T.C, Lp.C, self.vtor_C_total, self.vpol_C, B_t, B_p)

        # calculate nu_drags
        mbal_rhs_D = calc_mbal_rhs(self.mom_src_tor_D_tot, z_d, n.i, B_p, self.gamma_int_D) # Piper Changes: Uses integral cylindrical gamma
        mbal_rhs_C = calc_mbal_rhs(self.mom_src_tor_C_tot, z_c, n.C, B_p, self.gamma_C)

        nu_c_DC = 1 / calc_t90(m_d, m_c, z_d, z_c, n.C, T.i.J)
        nu_c_CD = 1 / calc_t90(m_c, m_d, z_c, z_d, n.i, T.i.J)
        
        #Piper changes: added alternate collision frequency calculation for comparison.
        nu_c_j_k = calc_nu_j_k(m_d, m_c, z_d, z_c, T.i.ev, n.C)
        nu_c_k_j = calc_nu_j_k(m_c, m_d, z_c, z_d, T.C.ev, n.i)

        self.nu_drag_D = calc_nu_drag(n.i, m_d, self.vtor_D_total, self.vtor_C_total, mbal_rhs_D, nu_c_DC)
        self.nu_drag_C = calc_nu_drag(n.i, m_d, self.vtor_D_total, self.vtor_C_total, mbal_rhs_C, nu_c_CD)
        
        ##############################################################
        # Pinch Velocity
        ##############################################################
        
        # Piper Changes: Added pinch velocity section and calculations.
        
        self.vrpinch_ext_term = calc_external_term(self.mom_src_tor_D_tot, n.i, ch_d, B_p)
        self.vrpinch_poloidal_term = calc_poloidal_term(n.i, m_d, ch_d, nu_c_DC, self.nu_drag_D, B_t, B_p, self.vpol_D)
        self.vrpinch_Er_term = calc_radial_E_field_term(n.i, m_d, ch_d, nu_c_DC, self.nu_drag_D, Er, B_p)
        self.vrpinch_toroidal_term = calc_toroidal_term(n.i, m_d, ch_d, nu_c_DC, B_p, self.vtor_C_total)
        self.vrpinch = calc_pinch_velocity(self.vrpinch_ext_term, self.vrpinch_poloidal_term, self.vrpinch_Er_term, self.vrpinch_toroidal_term)

        ##############################################################
        # energy balance
        ##############################################################

        # Piper Changes: Same as changes made to particle balance.
        en_src_nbi_i_D, en_src_nbi_i_D_tot, en_src_nbi_i_D_lost, en_src_nbi_i_D_kept = calc_en_src_nbi(beam_D, iol_adjusted=iolFlag, E_orb_nbi=E_orb_d_nbi)
        en_src_nbi_i_D2, en_src_nbi_i_D2_tot, en_src_nbi_i_D2_lost, en_src_nbi_i_D2_kept = calc_en_src_nbi(beam_D2, iol_adjusted=iolFlag, E_orb_nbi=E_orb_d_nbi)
        en_src_nbi_i_D3, en_src_nbi_i_D3_tot, en_src_nbi_i_D3_lost, en_src_nbi_i_D3_kept = calc_en_src_nbi(beam_D3, iol_adjusted=iolFlag, E_orb_nbi=E_orb_d_nbi)


        ################################################################################################################
        #
        #   NBI energy split - Currently 50:50 split ions and electrons
        #
        #   TODO: Implement accurate split
        #
        ################################################################################################################

        # calculate beam totals
        en_src_nbi_i = 0.5 * (en_src_nbi_i_D + en_src_nbi_i_D2 + en_src_nbi_i_D3)
        en_src_nbi_i_tot = .5 * (en_src_nbi_i_D_tot + en_src_nbi_i_D2_tot + en_src_nbi_i_D3_tot)
        en_src_nbi_i_lost = .5 * (en_src_nbi_i_D_lost + en_src_nbi_i_D2_lost + en_src_nbi_i_D3_lost)
        en_src_nbi_i_kept = .5 * (en_src_nbi_i_D_kept + en_src_nbi_i_D2_kept + en_src_nbi_i_D3_kept)

        en_src_nbi_e_D = np.zeros(en_src_nbi_i_D.shape)  # TODO: This isn't correct.
        en_src_nbi_e_D2 = np.zeros(en_src_nbi_i_D2.shape)  # TODO: This isn't correct.
        en_src_nbi_e_D3 = np.zeros(en_src_nbi_i_D3.shape)  # TODO: This isn't correct.

        #en_src_nbi_e = en_src_nbi_e_D + en_src_nbi_e_D2 + en_src_nbi_e_D3
        en_src_nbi_e = en_src_nbi_i_kept

        # calculate radial heat flux. Piper changes: Separated heat flux equations into differential and integral cylindrical methods.
        self.Qi_diff = calc_Qi_diff_method(r, core.a, calc_cxcool(core, n, T), calc_qie(n, T), en_src_nbi_i_kept, en_src_nbi_i_lost, core.dVdrho , iol_adjusted=iolFlag, E_orb=E_orb_d)  # previously called qheat. Differential Method.
        self.Qi_int, self.qie = calc_Qi_int_method(r, n, T, en_src_nbi_i_kept, cool_rate, iol_adjusted=iolFlag, E_orb=E_orb_d) # Integral method.
        self.Qe_diff = calc_Qe_diff_method(r, core.a, en_src_nbi_e, cool_rate, calc_qie(n, T)) # Differential Method.
        self.Qe_int = calc_Qe_int_method(r, n, T, en_src_nbi_e, cool_rate) # Integral method.

        self.chi = Chi(self, core, n, T, L, calc_nustar(nu_c_DC, core.q95, core.R0_a, self.vpol_C), reInterp=True)

        nn = namedtuple('nn', 's t tot')(
            core.n.n.s,  # slow
            core.n.n.t,  # thermal
            core.n.n.tot  # total
        )
        if debugFlag:
            """Beam source debugging"""
            nbiDebug(self.rhor, core.a, en_src_nbi_i_kept, en_src_nbi_e, part_src_nbi_kept, self.mom_src_nbi, core.dVdrho)
            """Particle flux debuggin"""
            gammaDebug(self.rhor, core.a, core.r2sa, core.dVdrho, part_src_nbi, part_src_nbi_kept, izn_rate, self.gamma_int_D, self.gamma_diff_D)
            """Energy flux debugging"""
            QDebug(self.rhor, core.a, core.r2sa, core.dVdrho, calc_qie(n,T,ion_species='D') , en_src_nbi_i, en_src_nbi_i_kept, cool_rate, E_orb_d, self.chi.Qi, T)

        ##############################################################
        # Release profiles used for rtransport calculations
        ##############################################################

        self.profiles = namedtuple('profiles', 'n T L nn')(
            namedtuple('n', 'i e')(
                core.n_fsa.i,
                core.n_fsa.e
            ),
            namedtuple('T', 'i e')(
                core.T_fsa.i,
                core.T_fsa.e
            ),
            namedtuple('L', 'n T')(
                namedtuple('n', 'i e')(
                    core.L_fsa.n.i,
                    core.L_fsa.n.e
                ),
                namedtuple('T', 'i e')(
                    core.L_fsa.T.i,
                    core.L_fsa.T.e)
            ),
            namedtuple('nn',' s t tot')(
                core.n_fsa.n.s,
                core.n_fsa.n.t,
                core.n_fsa.n.tot
            )
        )


def nbiDebug(rho, a, Qi, Qe, Si, Mi, dVdrho):
    """ """
    dVdr = dVdrho(rho)/a

    fig = plt.figure(figsize=(12, 8))
    fig.tight_layout()
    fig.suptitle(r'Radial Transport Debug Info: NBI')
    ax1 = fig.add_subplot(221)
    ax1.set_title(r'$Q^{nbi}_i$', fontsize=16)
    ax1.set_ylabel(r"""$Q^{nbi}_i$
                    $\left[\frac{W}{m^3}\right]$""", fontsize=16, rotation=0, ha='right')
    ax1.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax1.plot(rho, Qi)

    ax2 = fig.add_subplot(222)
    ax2.set_title(r'$Q^{nbi}_e$', fontsize=16)
    ax2.set_ylabel(r"""$Q^{nbi}_e$
                    $\left[J/{m^3}\right]$""", fontsize=16, rotation=0, ha='right')
    ax2.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax2.plot(rho, Qe)

    ax3 = fig.add_subplot(223)
    ax3.set_title(r'$S^{nbi}_i$', fontsize=16)
    ax3.set_ylabel(r"""$S^{nbi}_i$
                    $\left[\frac{\#}{m^3}\right]$""", fontsize=16, rotation=0, ha='right')
    ax3.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax3.plot(rho, Si)

    ax4 = fig.add_subplot(224)
    ax4.set_title(r'$M^{nbi}_i$', fontsize=16)
    ax4.set_ylabel(r"""$M^{nbi}_i$
                   $\left[\frac{N s}{m^3}\right]$""", fontsize=16, rotation=0, ha='right')
    ax4.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax4.plot(rho, Mi)

    enerTot = UnivariateSpline(rho*a, Qi * dVdr, k=3, s=0).integral(0., a)
    partTot = UnivariateSpline(rho*a, Si * dVdr, k=3, s=0).integral(0., a)
    volTot = UnivariateSpline(rho*a, dVdr, k=3, s=0).integral(0., a)

    print r"""NBI Debug info:
    
            {:<16}          {}
            {:<16}          {} MW
            {:<16}          {} $m^3$
            """.format("Total # particles", str(partTot),
                       "Total energy", str(enerTot/(1E6)),
                       "Total volume", str(volTot)
                       )

    plt.show(block=False)

def balance(gamma, interm, r, sa, x):
    return r[x], UnivariateSpline(r, interm, k=3, s=0).integral(0., r[x]), gamma[x]*sa(r[x])

def gammaDebug(rho, a, r2sa, dVdrho, Snbi_d_i, Snbi_kept_i, Sizn, gamma_i_D, gamma_D):

    dVdr = dVdrho(rho)/a

    fig = plt.figure(figsize=(12, 8))
    fig.tight_layout()
    fig.suptitle(r'Radial Transport Debug Info: Deuterium Radial Particle Flux')

    ax1 = fig.add_subplot(231)
    ax1.set_title(r'$S^{nbi}_{diff,i}$', fontsize=16)
    ax1.set_ylabel(r"""$S_i$
                       $\left[\frac{\#}{m^3 s}\right]$""", fontsize=16, rotation=0, ha='right')
    ax1.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax1.plot(rho, Snbi_d_i)

    ax2 = fig.add_subplot(232)
    ax2.set_title(r'$S^{nbi}_{int,i}$', fontsize=16)
    ax2.set_ylabel(r"""$S_i$
                       $\left[\frac{\#}{m^3 s}\right]$""", fontsize=16, rotation=0, ha='right')
    ax2.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax2.plot(rho, Snbi_kept_i)

    ax3 = fig.add_subplot(233)
    ax3.set_title(r'$S_{ion}$', fontsize=16)
    ax3.set_ylabel(r"""$S_i$
                       $\left[\frac{\#}{m^3 s}\right]$""", fontsize=16, rotation=0, ha='right')
    ax3.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax3.plot(rho, Sizn)

    ax4 = fig.add_subplot(234)
    ax4.set_title(r'$\Gamma_{r,int}$', fontsize=16)
    ax4.set_ylabel(r"""$\Gamma_{r}$
                        $\left[\frac{\#}{m^2 s}\right]$""", fontsize=16, rotation=0, ha='right')
    ax4.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax4.plot(rho, gamma_i_D)

    ax5 = fig.add_subplot(235)
    ax5.set_title(r'$\Gamma_{r,diff}$', fontsize=16)
    ax5.set_ylabel(r"""$\Gamma_{r}$
                        $\left[\frac{\#}{m^2 s}\right]$""", fontsize=16, rotation=0, ha='right')
    ax5.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax5.plot(rho, gamma_D)

    plt.show(block=False)

    partIn = UnivariateSpline(rho*a, (Snbi_kept_i + Sizn) * dVdr, k=3, s=0).integral(0., a)
    partOut = gamma_i_D[-1]*r2sa(a)
    totVol = UnivariateSpline(rho*a, dVdr, k=3, s=0).integral(0., a)

    print r"""Radial Particle Flux Debug info:

            {:<16}          {}
            {:<16}          {}
            {:<16}          {}
            {:<16}          {}
            """.format("Total # particles in", str(partIn),
                       "Total # particles out", str(partOut),
                       "Total volume", str(totVol),
                       "Surface area at LCFS", str(r2sa(a))
                       )
    print r"""r          Particles in      Particles out"""
    # for x in range(len(rho)):
    #     strList=balance(gamma_i_D, (Snbi_kept_i + Sizn) * dVdr, rho*a, r2sa, x)
    #     strDiff = 100.*np.abs(strList[1]-strList[2])/np.average(np.array(strList))
    #     print """{:<15}    {:<15}     {:<15}      {:<15}%""".format(strList[0], strList[1], strList[2], str(strDiff))

def QDebug(rho, a, r2sa, dVdrho, Qie, Qnbi_d_i, Qnbi_kept_i, coolrate, Eorb, Qi, T):
    # TODO: FINISH THIS
    dVdr = dVdrho(rho)/a

    fig = plt.figure(figsize=(12, 8))
    fig.tight_layout()
    fig.suptitle(r'Radial Transport Debug Info: Deuterium Radial Heat Flux')

    ax1 = fig.add_subplot(231)
    ax1.set_title(r'$Q^{nbi}_{diff,i}$', fontsize=16)
    ax1.set_ylabel(r"""$Q_i$
                       $\left[\frac{J}{m^3 s}\right]$""", fontsize=16, rotation=0, ha='right')
    ax1.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax1.plot(rho, Qnbi_d_i)

    ax2 = fig.add_subplot(232)
    ax2.set_title(r'$Q^{nbi}_{int,i}$', fontsize=16)
    ax2.set_ylabel(r"""$Q_i$
                       $\left[\frac{J}{m^3 s}\right]$""", fontsize=16, rotation=0, ha='right')
    ax2.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax2.plot(rho, Qnbi_kept_i)

    ax3 = fig.add_subplot(233)
    ax3.set_title(r'$Q_{ioniz}$', fontsize=16)
    ax3.set_ylabel(r"""$Q_i$
                       $\left[\frac{J}{m^3 s}\right]$""", fontsize=16, rotation=0, ha='right')
    ax3.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax3.plot(rho, coolrate)

    ax4 = fig.add_subplot(234)
    ax4.set_title(r'$Q_{ie}$', fontsize=16)
    ax4.set_ylabel(r"""$Q_{r,ie}$
                        $\left[\frac{J}{m^2 s}\right]$""", fontsize=16, rotation=0, ha='right')
    ax4.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax4.plot(rho, Qie)

    ax5 = fig.add_subplot(235)
    ax5.set_title(r'$dEdr$', fontsize=16)
    ax5.set_ylabel(r"""$dEdr{r}$
                        $\left[\frac{1}{m}\right]$""", fontsize=16, rotation=0, ha='right')
    ax5.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax5.plot(rho, Eorb)

    ax6 = fig.add_subplot(236)
    ax6.set_title(r'Temperatures (red = ion)', fontsize=16)
    ax6.set_ylabel(r"""T{r}$
                        $\left[eV\right]$""", fontsize=16, rotation=0, ha='right')
    ax6.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax6.plot(rho, T.i.ev, color='red')
    ax6.plot(rho, T.e.ev, color='blue')

    plt.show(block=False)

    energyIn = UnivariateSpline(rho*a, (Qnbi_d_i - coolrate) * dVdr, k=3, s=0).integral(0., a)
    energyOut = Qi[-1]*r2sa(a)
    totVol = UnivariateSpline(rho*a, dVdr, k=3, s=0).integral(0., a)

    print r"""Radial Energy Flux Debug info:

            {:<16}          {}  MW
            {:<16}          {}  MW
            {:<16}          {}
            {:<16}          {}
            """.format("Total energy in", str(energyIn / (1E6)),
                       "Total energy out", str(energyOut / (1E6)),
                       "Total volume", str(totVol),
                       "Surface area at LCFS", str(r2sa(a))
                       )
    # print r"""r          Energy in      Energy out"""
    # for x in range(len(rho)):
    #     strList=balance(Qi, (Qnbi_kept_i - coolrate) * dVdr, rho*a, r2sa, x)
    #     strDiff = 100.*np.abs(strList[1]-strList[2])/np.average(np.array(strList))
    #     print """{:<15}    {:<15}     {:<15}      {:<15}%""".format(strList[0], strList[1], strList[2], str(strDiff))

