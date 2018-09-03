#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 08:58:46 2018

@author: max
"""
from __future__ import division
import numpy as np
from math import pi, cos, sqrt, exp # Piper Changes: Imported exp for solving balance equations.
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

    core.cool_rate_fsa = calc_fsa(core.n.e * core.n.C * np.nan_to_num(core.Lz.s) +
                                  core.n.e * core.n.C * np.nan_to_num(core.Lz.t),
                                  core.R, core.Z)

    """
    Cuts off neutrals calculations at 0.8 to give realistic neutral shit
    """

    core.izn_rate_fsa = np.array(map(lambda x: core.izn_rate_fsa[-1] * x**10, core.r[:,0]/core.a))
    core.cool_rate_fsa = np.array(map(lambda x: core.cool_rate_fsa[-1] * x**10, core.r[:,0]/core.a))

def corePatch(core):
    """
    Updates deuterium ion density and zeff. D density is invalid if no D density file is given because it will be
    set to 0. This screws up subsequent calculations in non-obvious ways. As an example, z_eff calculation will
    include D density as 0s but still have a carbon density from the fracz input file, giving an invalid result

    This is addressed by recreating the n_fsa namedtuple and z_eff_fsa in full, as namedtuples cannot be changed piecemeal

    :param core:
    :return:
    """

    core.n = namedtuple('n', 'i e n C')(core.n.e/(1.+.025*6.0), core.n.e, core.n.n, 0.025 * core.n.e/(1.+.025*6.0))   # TODO: Update 0.025 and 6.0 to actual fracz and zbar2

    core.n_fsa = namedtuple('n', 'i e n C')(
        calc_fsa(core.n.i, core.R, core.Z),
        calc_fsa(core.n.e, core.R, core.Z),
        namedtuple('nn', 's t tot')(
            calc_fsa(core.n.n.s, core.R, core.Z),  # slow
            calc_fsa(core.n.n.t, core.R, core.Z),  # thermal
            calc_fsa(core.n.n.tot, core.R, core.Z)  # total
        ),
        calc_fsa(core.n.C, core.R, core.Z))

    core.z_eff_fsa = calc_fsa((core.n.i * (1.**2) + core.n.C * (6.0**2))/(core.n.i * 1.0 + core.n.C * 6.0), core.R, core.Z) #TODO Similar updates (1.0 = atnum, 6.0 = zbar2)


def calc_part_src_nbi(beam, iol_adjusted=False, F_orb_nbi=None):
    """
    """

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


def calc_mom_src_nbi(beam, n, z_eff, R0_a,fforb):
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
        # Piper Changes: changed beam.E.J to beam.E.kev
        xlam = 5.5E17 * beam.E.kev / (2.0 * 0.5 * n.i * z_eff ** 0.5) # Piper Changes: beam.a doesn't exist anymore. changed to a hard coded 2.
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

    # Piper Changes: Changed function name to denote differential method. If the name sux, feel free to change it.
def calc_Qi_diff_method(r, n, T, en_src_nbi_i, cool_rate, r2sa, iol_adjusted=False, E_orb=None):  # formerly qheat

    if iol_adjusted:
        diff_E_orb = UnivariateSpline(r, E_orb, k=1, s=0).derivative()(r)
    else:
        diff_E_orb = np.zeros(r.shape)

    qie = calc_qie(n,T,ion_species='D')

    dQi_dr_interp = interp1d(r, (en_src_nbi_i - qie - cool_rate)*(1-diff_E_orb), fill_value='extrapolate')

    # ion energy balance equation
    def dQi_dr(Qi_0, r):
        return dQi_dr_interp(r)

    # boundary condition at magnetic axis
    Qi_0 = en_src_nbi_i[0]  # or something like this

    # solve ODE and divide by surface area
    # note that r2sa is an interpolation function from core that calculates the
    # surface area for a flux surface based on its corresponding value of r
    surf_area = r2sa(r)
    surf_area[0] = surf_area[1]  # this is just to prevent divide by zero errors at rho=0. Will be replaced later anyway

    Qi = odeint(dQi_dr, Qi_0, r, hmax=0.1).T[0] / surf_area

    # set the center value of Qi equal to the next point radially outward.
    # Qi(0) is mathematically infinite, so this just gives us a reasonable representation.
    Qi[0] = Qi[1]

    return Qi
    
    # Piper Changes: Added cylindrical integral method as a separate function. This will be called separately in the main code.
def calc_Qi_int_method(r, n, T, en_src_nbi_i_kept, cool_rate, iol_adjusted=False, E_orb=None):  # formerly qheat

    Qi = np.zeros(r.shape)
    
    qie = calc_qie(n,T,ion_species='D')
    
    # Boundary condition at the magnetic axis.
    # Only has the second term, since it's the center value. Also uses delta_r of the next point.
    # If not adjusted for IOL, en_src_nbi_kept = en_src_nbi_tot, so no need for an IOL check.
    Qi[0] = (en_src_nbi_i_kept[0] - qie[0])*(r[1] - r[0])
    Qi[1] = Qi[0] + (en_src_nbi_i_kept[1] - qie[1])*(r[1] - r[0])
    
    # Integral cylindrical form of the energy balance equation.
    # Identical in form to the particle continuity equation, but different source term.
    for n in range(2,len(r)):
        if iol_adjusted:
            # Imported the exp() function for the thermal IOL attenuation. 
            Qi[n] = (r[n]/r[n-1]) * Qi[n-1] * exp(-(E_orb[n]-E_orb[n-1])) + (en_src_nbi_i_kept[n] - qie[n])*(r[n] - r[n-1])
        else:
            Qi[n] = (r[n]/r[n-1]) * Qi[n-1] + (en_src_nbi_i_kept[n] - qie[n])*(r[n] - r[n-1])
            
    return Qi

# OLD VERSION
#
# rhovals = rho[:, 0]
#
# qie = calc_qie(n, T, ion_species='D')
#
# qheat = np.zeros(rho.shape)
# for i, rho in enumerate(rhovals):
#     if rho > 0:
#
#         if iol_adjusted:
#             xponq = np.exp(-1 * (E_orb[i] - E_orb[i-1]))
#         else:
#             xponq = 1
#         cxcool = 1.5 * (n.i[i] + n.i[i-1]) / 2 * (T.i.kev[i, :] + T.i.kev[i-1]) / 2 * e * xnuati
#         srprimq = qnbi[i] - cxcool - qie[i]
#         qheat[i] = (rho[i-1] / rho[i]) * qheat[i-1] * xponq + srprimq * delma
# return qheat


def calc_Qe_diff_method(r, n, T, en_src_nbi_e, cool_rate, r2sa): # Piper Changes: Same as Qi changes.

    qie = calc_qie(n, T, ion_species='D')
    dQe_dr_interp = interp1d(r, en_src_nbi_e - qie - cool_rate, fill_value='extrapolate')

    # electron energy balance equation
    def dQe_dr(Qe_0, r):
        return dQe_dr_interp(r)

    # boundary condition at magnetic axis
    Qe_0 = en_src_nbi_e[0]  # or something like this

    # solve ODE and divide by surface area
    # note that r2sa is an interpolation function from core that calculates the
    # surface area for a flux surface based on its corresponding value of r

    surf_area = r2sa(r)
    surf_area[0] = surf_area[1]  # this is just to prevent divide by zero errors at rho=0. Will be replaced later anyway

    Qe = odeint(dQe_dr, Qe_0, r, hmax=0.1).T[0] / surf_area

    # set the center value of Qi equal to the next point radially outward.
    # Qe(0) is mathematically infinite, so this just gives us a reasonable representation.
    Qe[0] = Qe[1]

    return Qe
    
def calc_Qe_int_method(r, n, T, en_src_nbi_e_tot, cool_rate): # Piper Changes: Same as Qi changes.

    Qe = np.zeros(r.shape)
    qie = calc_qie(n,T,ion_species='D')
    
    Qe[0] = (en_src_nbi_e_tot[0] + qie[0] - cool_rate[0])*(r[1] - r[0])
    Qe[1] = Qe[0] + (en_src_nbi_e_tot[1] + qie[1] - cool_rate[1])*(r[1] - r[0])
    
    # Integral cylindrical form of the energy balance equation.
    # Identical in form to the continuity equation, but different source term.
    for n in range(2,len(r)):
        Qe[n] = (r[n]/r[n-1]) * Qe[n-1] + (en_src_nbi_e_tot[n] + qie[n] - cool_rate[n])*(r[n] - r[n-1])

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

    
def calc_gamma_diff_method(r, part_src_nbi, izn_rate, r2sa, iol_adjusted=False, F_orb=None):
    # Piper Changes: Changed function to denote differential method.
    if iol_adjusted:
        diff_F_orb = UnivariateSpline(r, F_orb, k=1, s=0).derivative()(r)
    else:
        diff_F_orb = np.zeros(r.shape)

    # the two is to account for the inward radial current, assumed to be deuterium as well
    dgamma_dr_interp = interp1d(r, (izn_rate + part_src_nbi)*(1-2*diff_F_orb), fill_value='extrapolate')

    # continuity equation
    def dgamma_dr(gamma_0, r):
        return dgamma_dr_interp(r)

    # boundary condition at magnetic axis

    gamma_0 = part_src_nbi[0]  # or something like this

    # solve ODE and divide by surface area
    # note that r2sa is an interpolation function from core that calculates the
    # surface area for a flux surface based on its corresponding value of r

    surf_area = r2sa(r)
    surf_area[0] = surf_area[1]  # this is just to prevent divide by zero errors at rho=0. Will be replaced later anyway

    gamma = odeint(dgamma_dr, gamma_0, r).T[0] / surf_area

    # set the center value of gamma equal to the next point radially outward.
    # Gamma(0) is mathematically infinite, so this just gives us a reasonable representation.
    gamma[0] = gamma[1]

    return gamma
    
    
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
            gamma[n] = (r[n]/r[n-1]) * gamma[n-1] * exp(-2*(F_orb[n]-F_orb[n-1])) + (part_src_nbi_tot[n] - 2*part_src_nbi_lost[n] + izn_rate[n])*(r[n] - r[n-1])
        else:
            gamma[n] = (r[n]/r[n-1]) * gamma[n-1] + (part_src_nbi_tot[n] + izn_rate[n])*(r[n] - r[n-1])

    return gamma

# OLD VERSION
#
# rhovals = rho[:, 0]
#
# xpon = np.zeros(rho.shape)
# gamma = np.zeros(rho.shape)
# for i, rho in enumerate(rhovals):
#     if rho > 0:
#
#         if iol_adjusted:
#             xpon[i] = np.exp(-2.0 * (F_orb[i] - F_orb[i-1]))
#         else:
#             xpon[i] = 1
#
#         srprim = src_nbi[i] + 0.5 * (n.i[i] + n.i[i-1]) * xnuioni * (1 + fracz[i] * zbar2[i])
#         gamma[i] = rho[i-1] / rho[i] * gamma[i-1] * xpon[i] + srprim * delma
# return gamma

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
        Jr_iol[n] = (r[n]/r[n-1]) * Jr_iol[n-1] + (part_src_nbi_lost[n] + (gamma[n] * diff_F_orb[n]))*(r[n] - r[n-1]) * ch_d # No factor of 2, because only 1 direction.

    return Jr_iol


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

def calc_vpol(Er, vphi, Lp, T, n, z, B_t, B_p,):
    vpol = (1.0/B_t) * (1.0/(n.i*e*z) * Lp.i * T.i.kev + vphi * B_p - Er)
    return vpol


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


def calc_reduced_mass(m1, m2):
    return m1 * (1 + (m1 / m2))  # * 1E3  # TODO: not sure why this 1E3 is here


def calc_chi(r, n, T, L, cool_rate, gamma, conv_mult=3/2, visc=False, pressure=False, iol_adjusted=False, E_orb=None):

    q_cond = calc_Qi(r, n, T, en_src_nbi_i, cool_rate, iol_adjusted=True, E_orb=E_orb) - conv_mult * gamma * T.i.J

    if pressure == True:
        q_cond = q_cond - calc_heatin()

    if visc == True:
        q_cond = q_cond - calc_visc_heat()

    chi = L.T.i / (n.i * T.i.J) * q_cond
    return chi


def calc_heatin(gamma, vtor_D, vpol_D):  # TODO: What quantity are we calculating here?
    # Looks like 1/2 mV^2 * Gamma. Is that the energy flux? Why is it considered heat in and not heat out?
    heatin = gamma * 0.5 * m_d * (vtor_D ** 2 + vpol_D ** 2)
    return heatin


def calc_visc_heat(a, R0_a, kappa, n, T, q95, vpol_D, gammahat, vtor_Dhat, B_p, B_t, q, nu_c_DD):
    # TODO: This function is a train wreck.
    se = np.sqrt(0.5 * (1 + kappa ** 2))  # TODO: WTF are we calcuating here?
    ep = a * se / R0_a  # TODO: WTF are we calcuating here?

    #  f varies with theta, not sure about the best way to make 1D - MH
    fp = B_p / B_t  # TODO: This is already calculated in core.

    xnustar11 = 0.

    # QUESTIONABLE CALCULATION COMMENTED OUT - JRo
    #        for a, b in zip(data.xnuc[0, 1], data.vpolD):
    #            xnustar11=xnustar11+a*abs(data.q95)*data.rmajor/b

    xnustar11 = calc_nustar(nu_c_DD, q95, R0_a, vpol_D)

    eff = xnustar11 / ((1 + xnustar11) * (ep ** 1.5 + xnustar11))
    vrad1 = gammahat / n.i

    eta0 = n.i * m_d * vpol_D * q * R0_a * eff
    eta4 = n.i * m_d * T.i.kev * e / (ch_D * abs(B_t))
    #       Calculate viscous heating:  a=vtord, b=fp, c = eta0, d=vrad1, f = eta4, g= vpold
    #   TODO: THIs does not match POP17-052504
    asymR = 0.1 / R0_a

    visc_heat = asymR * vtor_Dhat * (fp * eta0 * vrad1 - 0.5 * eta4 * (4. * vtor_Dhat + vpol_D)) - \
           0.5 * vpol_D * (eta0 * vrad1 + eta4 * (vtor_Dhat + 0.5 * vpol_D))
    return visc_heat


#    Piper Changes: Added calculation of pinch velocity components and the pinch velocity itself.
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

def calc_chi_i(self, args):

    return

class Chi:
    """
    Chi class claculates various chis and provides source terms as necessary
    """
    def __init__(self, data, n, T, L, nustar):

        self.Qi_i = data.Qi_i
        self.Qe_i = data.Qe_i
        self.conv15 = 1.5 * ch_d * data.gamma_i_D * T.i.ev
        self.conv25 = 2.5 * ch_d * data.gamma_i_D * T.i.ev
        self.heatvisc = np.zeros(data.gamma_i_D.shape)
        self.heatin = 0.5 * data.gamma_i_D * m_d * (data.vtor_D_total**2 + data.vpol_D**2) # TODO: Provide logic that uses vtor_D_intrin/fluid depending on IOL Switch, currently too small to matter

        self.chi = namedtuple('chi','i e')(
            namedtuple('i','chi1 chi2 chi3 chi4')(
                (self.Qi_i) * L.T.i / (n.i * T.i.ev * ch_d),
                (self.Qi_i - self.conv25) * L.T.i / (n.i * T.i.ev * ch_d),
                (self.Qi_i - self.conv25 - self.heatin) * L.T.i / (n.i * T.i.ev * ch_d),
                (self.Qi_i - self.conv25 - self.heatin * self.heatvisc) * L.T.i / (n.i * T.i.ev * ch_d)
            ),self.calc_chi_e(data, n, L, T)
        )

    def calc_chi_e(self, data, n, L, T):
        gameltemp = 1.0* data.gamma_i_D + 6.0 * data.gamma_C   #OHHHH gamma electron!!! HA HA HA HA WE DON'T KNOW WHAT THE FUCK GAMMA_C IS

        return L.T.e * ((self.Qe_i / (ch_d * n.e * T.e.ev)) - 1.5 * gameltemp / n.e)


    def calc_heatvisc(self, data):
        se = sqrt(0.5 * (1 + ((data.kappa_vals[0] + data.kappa_vals[1])/2)**2))
        ep = data.a * se / data.R0_a
        fp = data.B_p_fsa / data.B_t_fsa

        # TODO : finish this clusterfuck

        return None


class RadialTransport(Chi):
    
    def __init__(self, inp, core, iol, nbi, ntrl=False, iolFlag=True):
        sys.dont_write_bytecode = True

        ##############################################################
        # prep quantities for 1D transport analysis
        ##############################################################

        #prepare beams object
        beam_D = nbi.beams.D1 # Piper changes: Changed beam object references to match beamdep.
        beam_D2 = nbi.beams.D2
        beam_D3 = nbi.beams.D3

        if ntrl: neutPatch(core) # Patch to update neutrals-related quantities
        corePatch(core) # Patch to update values not brought in via files (ni, zeff)

        # prepare core and iol quantities
        r = core.r.T[0]    # TODO: Should this be a flux surface average?
        izn_rate = core.izn_rate_fsa  # TODO: Should this be a flux surface average or a flux surface total?
        cool_rate = core.cool_rate_fsa  # TODO: Should this be a flux surface average or a flux surface total?
        n = core.n_fsa
        T = core.T_fsa
        L = core.L_fsa
        z_eff = core.z_eff_fsa
        R0_a = core.R0_a
        vol_P = core.vol
        B_p = core.B_p_fsa
        B_t = core.B_t_fsa
        Er = core.E_r_fsa
        #Lp = core.Lp_fsa

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
        
        # Piper changes: Calculate gamma using the differential and integral methods, respectively.
        self.gamma_d_D = calc_gamma_diff_method(r, part_src_nbi, izn_rate, core.r2sa, iol_adjusted=iolFlag, F_orb=F_orb_d) # Differential Cylindrical Method
        self.gamma_i_D = calc_gamma_int_method(r, part_src_nbi_tot, part_src_nbi_lost, izn_rate, iol_adjusted=iolFlag, F_orb=F_orb_d) # Integral Cylindrical Method
        self.gamma_C = np.zeros(self.gamma_i_D.shape)
        
        # Piper changes: Calculate radial return current (Uses integral cylindrical gamma)
        self.jr_iol = calc_return_cur(r, part_src_nbi_lost, self.gamma_i_D, izn_rate, ch_d, iol_adjusted=iolFlag, F_orb=F_orb_d)

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
                                                 self.gamma_i_D,
                                                 self.gamma_C) # Piper Changes: Uses integral cylindrical gamma
        else:
            self.vtor_D_total = core.v_1D.tor.C

        self.vtor_fluid_D = self.vtor_D_total - self.vtor_D_intrin

        # calculate carbon and deuterium poloidal rotation
        try:
            self.vpol_C = core.v_1D.pol.C
            #self.vpol_D = calc_vpol(Er, self.vtor_D_total, Lp, T, n, z_d, B_t, B_p,)
            self.vpol_D = self.vpol_C / 0.4
        except:
            self.vpol_D = self.vpol_C / 0.4
            print 'could not calculate deuterium poloidal rotation'
            pass

        # calculate nu_drags
        mbal_rhs_D = calc_mbal_rhs(self.mom_src_tor_D_tot, z_d, n.i, B_p, self.gamma_i_D) # Piper Changes: Uses integral cylindrical gamma
        mbal_rhs_C = calc_mbal_rhs(self.mom_src_tor_C_tot, z_c, n.C, B_p, self.gamma_C)

        nu_c_DC = 1 / calc_t90(m_d, m_c, z_d, z_c, n.C, T.i.J)
        nu_c_CD = 1 / calc_t90(m_c, m_d, z_c, z_d, n.i, T.i.J)

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
        self.Qi_d = calc_Qi_diff_method(r, n, T, en_src_nbi_i, cool_rate, core.r2sa, iol_adjusted=iolFlag, E_orb=E_orb_d)  # previously called qheat. Differential Method.
        self.Qi_i = calc_Qi_int_method(r, n, T, en_src_nbi_i_kept, cool_rate, iol_adjusted=iolFlag, E_orb=E_orb_d) # Integral method.
        self.Qe_d = calc_Qe_diff_method(r, n, T, en_src_nbi_e, cool_rate, core.r2sa) # Differential Method.
        self.Qe_i = calc_Qe_int_method(r, n, T, en_src_nbi_e, cool_rate) # Integral method.

        self.chi = Chi(self, n, T, L, nustar=None)

        nn = namedtuple('nn', 's t tot')(
            core.n.n.s,  # slow
            core.n.n.t,  # thermal
            core.n.n.tot  # total
        )
        plt.plot(r, self.gamma_i_D)
        plt.show()
