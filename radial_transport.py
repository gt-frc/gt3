#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 08:58:46 2018

@author: max
"""
from __future__ import division
import numpy as np
from math import pi, cos, sqrt
import sys
import matplotlib.pyplot as plt
from collections import namedtuple
from scipy import constants
from scipy.integrate import odeint
from scipy.interpolate import interp1d

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


def calc_part_src_nbi(beam, iol_adjusted=False, F_orb_nbi=None):
    """
    """
    # # old suspicious calculation
    # snbi = (0.624E25 * beam.p.W / (vol * beam.e.J)) * beam.num * beam.dp

    # new shiny calculation
    # TODO: verify this calculation

    part_src_nbi = beam.p.W / beam.e.J * beam.dp
    if iol_adjusted:
        part_src_nbi = part_src_nbi * (1 - F_orb_nbi)

    return part_src_nbi


def calc_en_src_nbi(beam, iol_adjusted=False, F_orb_nbi=None):
    """
    """
    # # old suspicious calculation
    # snbi = (0.624E25 * beam.p.W / (vol * beam.e.J)) * beam.num * beam.dp

    # new shiny calculation
    # TODO: verify this calculation
    part_src_nbi = beam.p.W / beam.e.J * beam.dp
    if iol_adjusted:
        part_src_nbi = part_src_nbi * (1 - F_orb_nbi)

    return part_src_nbi


def calc_mom_src_nbi(beam, n, z_eff, R0_a):
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
        xlam = 5.5E17 * (beam.e.J / 1) / (beam.a * 0.5 * n.i * z_eff ** 0.5)
        atten = 1 - np.exp(-delma / cos(alphain) / xlam)
        return atten

    atten = calc_atten(beam, n, z_eff)
    unatten = 1-atten  # might need to flip atten...
    torque = calc_torque(beam)

    tor_mom = unatten*atten*torque/R0_a
    return tor_mom


def calc_torque(beam):
    """ Calculates torque from a neutral beam (or beam component)

    torque = F * r_tan = (P/v) * r_tan = (P/sqrt(2E/m)) * r_tan = P * sqrt(m/(2E)) * r_tan

    :param beam: beam object with attributes z, m, a, en, pwr, rtan
    :return: torque
    """

    power = beam.p.W
    energy = beam.e.J
    mass = beam.m
    rtan = beam.rtan

    torque = power * np.sqrt(0.5 * mass / energy) * rtan
    return torque


def calc_Qi(r, n, T, en_src_nbi_i, cool_rate, iol_adjusted=False, E_orb=None):  # formerly qheat

    qie = calc_qie(n,T,ion_species='D')

    dQi_dr_interp = interp1d(r, en_src_nbi_i - qie - cool_rate, fill_value='extrapolate')

    # continuity equation
    def dQi_dr(Qi_0, r):
        return dQi_dr_interp(r)

    # boundary condition at magnetic axis
    Qi_0 = en_src_nbi_i[0]  # or something like this

    # solve ODE
    Qi = odeint(dQi_dr, Qi_0, r, hmax=0.1)

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


def calc_Qe(r, n, T, en_src_nbi_e, cool_rate):

    qie = calc_qie(n, T, ion_species='D')
    dQe_dr_interp = interp1d(r, en_src_nbi_e - qie - cool_rate, fill_value='extrapolate')

    # continuity equation
    def dQe_dr(Qe_0, r):
        return dQe_dr_interp(r)

    # boundary condition at magnetic axis
    Qe_0 = en_src_nbi_e[0]  # or something like this

    # solve ODE
    Qe = odeint(dQe_dr, Qe_0, r, hmax=0.1)

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

    return qie


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


def calc_gamma(r, part_src_nbi, izn_rate, iol_adjusted=False, F_orb=None):

    dgamma_dr_interp = interp1d(r, (izn_rate + part_src_nbi), fill_value='extrapolate')

    # continuity equation
    def dgamma_dr(gamma_0, r):
        return dgamma_dr_interp(r)

    # boundary condition at magnetic axis
    gamma_0 = part_src_nbi[0]  # or something like this

    # solve ODE
    gamma = odeint(dgamma_dr, gamma_0, r)

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

def calc_mbal_rhs(mom_src_ext, z, n, B_p, gamma):
    """ The function formerly known as y11 and y22
    """

    mbal_rhs = mom_src_ext + (z * e) * (n * E_phi + B_p * gamma)

    return mbal_rhs


def calc_intrin_rot(M_orb, T_J, m):

    intrin_rot = 2 / sqrt(pi) * M_orb * np.sqrt(2 * T_J / m)
    return intrin_rot


def calc_vtor_d_pert(vtor_C_total, vtor_C_intrin, vtor_D_intrin, mom_src_ext, z, n, T, B_p, gamma):
    """
    """
    mbal_rhs_D = calc_mbal_rhs(mom_src_ext, z, n, B_p, gamma)
    mbal_rhs_C = calc_mbal_rhs(mom_src_ext, z, n, B_p, gamma)

    nu_c_DC = 1 / calc_t90(m_d, m_c, 1, 6, n.C, T.i.J)

    #vtor_C_total = vtor_fluid_C + vtor_intrin_C
    del_v0 = vtor_D_intrin - vtor_C_intrin

    nu_drag_approx = (mbal_rhs_D + mbal_rhs_C) / ((n.i * m_d + n.C * m_c) * vtor_C_total + n.i * m_d * del_v0)

    del_v1 = (mbal_rhs_D - n.i * m_d * nu_drag_approx * vtor_C_total) / (n.i * m_d * (nu_c_DC + nu_drag_approx))

    vtorDPert = vtor_C_total + del_v1
    return vtorDPert


def calc_nu_drag(n_j, m_j, v_tor_j, v_tor_k, mbal_rhs, nu_c):
    nu_drag = (mbal_rhs + nu_c * v_tor_k) / (v_tor_j * n_j * m_j) - nu_c
    return nu_drag


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
           n2 * ((z1*e) * (z2 * e))**2 * coul_log
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


def calc_nustar(nu_c, q95, R0_a, vpol):
    nustar = nu_c * abs(q95) * R0_a / vpol
    return nustar

# def calc_xpartdot_old(beam, n, z_eff, rpts):  # TODO: wtf is xpartdot?
#     atten = calc_atten(beam, n, z_eff, rpts)
#     unatten = 1-atten  # might need to flip atten...
#     volm =
#     xpartdot = unatten * atten * beam.p * 1E6 / \
#     (beam.e.J / comp * 1E3 * e) / volm
#     return xpartdot


# def calc_qnbi_old(beam, n, z_eff, rpts):
#     xpartdot1 = calc_xpartdot(beam.e, beam.pwr, pwr_frac, 1, beam.a, n, z_eff, rpts)
#     xpartdot2 = calc_xpartdot(beam.e, beam.pwr, pwr_frac, 2, beam.a, n, z_eff, rpts)
#     xpartdot3 = calc_xpartdot(beam.e, beam.pwr, pwr_frac, 3, beam.a, n, z_eff, rpts)
#     qnbi = (xpartdot1 + xpartdot2/2 + xpartdot3/3) * 1E3 * e * e_beam
#     return qnbi


class RadialTransport:
    
    def __init__(self, inp, core, iol, nbi):
        sys.dont_write_bytecode = True

        ##############################################################
        # prep quantities for 1D transport analysis
        ##############################################################

        #prepare beams object
        beam_D = nbi.beams_1D.D
        beam_D2 = nbi.beams_1D.D2
        beam_D3 = nbi.beams_1D.D3

        # prepare core and iol quantities
        r = core.r[:, 0]    # TODO: Should this be a flux surface average?
        izn_rate = core.izn_rate_fsa  # TODO: Should this be a flux surface average or a flux surface total?
        cool_rate = core.cool_rate_fsa  # TODO: Should this be a flux surface average or a flux surface total?
        n = core.n_fsa
        T = core.T_fsa
        z_eff = core.z_eff_fsa
        R0_a = core.R0_a
        vol_P = core.vol
        B_p = core.B_p_fsa

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

        ##############################################################
        # particle balance
        ##############################################################

        part_src_nbi_D = calc_part_src_nbi(beam_D, iol_adjusted=False, F_orb_nbi=None)
        part_src_nbi_D2 = calc_part_src_nbi(beam_D2, iol_adjusted=False, F_orb_nbi=None)
        part_src_nbi_D3 = calc_part_src_nbi(beam_D3, iol_adjusted=False, F_orb_nbi=None)
        part_src_nbi = part_src_nbi_D + part_src_nbi_D2 + part_src_nbi_D3

        # calculate radial heat flux
        self.gamma_D = calc_gamma(r, part_src_nbi, izn_rate, iol_adjusted=False, F_orb=None).T[0]
        self.gamma_C = np.zeros(self.gamma_D.shape)

        ##############################################################
        # momentum balance
        ##############################################################

        # calculate toroidal momentum source rates from beams
        mom_src_nbi_D = calc_mom_src_nbi(beam_D, n, z_eff, R0_a)
        mom_src_nbi_D2 = calc_mom_src_nbi(beam_D2, n, z_eff, R0_a)
        mom_src_nbi_D3 = calc_mom_src_nbi(beam_D3, n, z_eff, R0_a)
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
        try:
            self.vtor_D_total = core.v_1D.tor.C
        except:
            self.vtor_D_total = calc_vtor_d_pert(self.vtor_C_fluid,
                                                 self.vtor_C_intrin,
                                                 self.vtor_D_intrin,
                                                 self.mom_src_tor_D_tot,
                                                 1,
                                                 n,
                                                 B_p,
                                                 self.gamma_D)

        self.vtor_fluid_D = self.vtor_D_total - self.vtor_D_intrin

        # calculate carbon and deuterium poloidal rotation
        try:
            self.vpol_C = core.v_1D.pol.C
            self.vpol_D = self.vpol_C / 0.4
        except:
            print 'no poloidal rotation data available'
            pass

        # calculate nu_drags
        mbal_rhs_D = calc_mbal_rhs(self.mom_src_tor_D_tot, z_d, n.i, B_p, self.gamma_D)
        mbal_rhs_C = calc_mbal_rhs(self.mom_src_tor_C_tot, z_c, n.C, B_p, self.gamma_C)

        nu_c_DC = 1 / calc_t90(m_d, m_c, z_d, z_c, n.C, T.i.J)
        nu_c_CD = 1 / calc_t90(m_c, m_d, z_c, z_d, n.i, T.i.J)

        self.nu_drag_D = calc_nu_drag(n.i, m_d, self.vtor_D_total, self.vtor_C_total, mbal_rhs_D, nu_c_DC)
        self.nu_drag_C = calc_nu_drag(n.i, m_d, self.vtor_D_total, self.vtor_C_total, mbal_rhs_C, nu_c_CD)

        ##############################################################
        # energy balance
        ##############################################################

        en_src_nbi_i_D = part_src_nbi_D * beam_D.e.J
        en_src_nbi_i_D2 = part_src_nbi_D * beam_D.e.J
        en_src_nbi_i_D3 = part_src_nbi_D * beam_D.e.J
        en_src_nbi_i = en_src_nbi_i_D + en_src_nbi_i_D2 + en_src_nbi_i_D3

        en_src_nbi_e_D = np.zeros(en_src_nbi_i_D.shape)  # TODO: This isn't correct.
        en_src_nbi_e_D2 = np.zeros(en_src_nbi_i_D2.shape)  # TODO: This isn't correct.
        en_src_nbi_e_D3 = np.zeros(en_src_nbi_i_D3.shape)  # TODO: This isn't correct.
        en_src_nbi_e = en_src_nbi_e_D + en_src_nbi_e_D2 + en_src_nbi_e_D3

        # calculate radial heat flux
        self.Qi = calc_Qi(r, n, T, en_src_nbi_i, cool_rate, iol_adjusted=False, E_orb=None)  # previously called qheat
        self.Qe = calc_Qe(r, n, T, en_src_nbi_e, cool_rate)

        # # calculate chi
        # chi = calc_chi(r,
        #                n,
        #                T,
        #                L,
        #                cool_rate,
        #                self.gamma,
        #                conv_mult=3/2,
        #                visc=False,
        #                pressure=False,
        #                iol_adjusted=False,
        #                E_orb=None)
        plt.plot(r, self.gamma_D)
        plt.show()
