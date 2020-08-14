#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 08:58:46 2018

@author: max
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
from collections import namedtuple
from scipy import constants
from Chi import Chi
from GT3.RadialTransport.Functions.CorePatch import corePatch
from GT3.RadialTransport.Functions.CalcPartSrcNBI import calc_part_src_nbi
from GT3.RadialTransport.Functions.CalcMomSrcNBI import calc_mom_src_nbi
from GT3.RadialTransport.Functions.CalcReturnCur import calc_return_cur
from GT3.RadialTransport.Functions.CalcNu import calc_nu_j_k, calc_nustar, calc_nu_drag
from GT3.RadialTransport.Functions.CalcMbalRHS import calc_mbal_rhs
from GT3.RadialTransport.Functions.CalcT90 import calc_t90
from GT3.RadialTransport.Functions.CalcQ import calc_Qe_diff_method, calc_qie, calc_Qe_int_method, calc_Qi_int_method, \
    calc_Qi_diff_method
from GT3.RadialTransport.Test.QDebug import QDebug
from GT3.RadialTransport.Test.GammaDebug import gammaDebug
from GT3.RadialTransport.Test.NBIDebug import nbiDebug
from GT3.RadialTransport.Test.IOLDebug import IOLDebug
from GT3.RadialTransport.Functions.NeutPatch import neutPatch
from GT3.RadialTransport.Functions.CalcGamma import calc_gamma_diff_method, calc_gamma_int_method
from GT3.RadialTransport.Functions.CalcIntrinRot import calc_intrin_rot
from GT3.RadialTransport.Functions.CalcVTorDPert import calc_vtor_d_pert
from GT3.RadialTransport.Functions.CalcErMomBal import calc_Er_mom_bal
from GT3.RadialTransport.Functions.CalcEnSrcNBI import calc_en_src_nbi
from GT3.RadialTransport.Functions.CalcErIOL import calc_Er_iol
from GT3.RadialTransport.Functions.CalcVpol import calc_vpol
from GT3.RadialTransport.Functions.CalcCXCool import calc_cxcool

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
    pol_term = (n_j * m_j * (nu_jk + nu_dj) * B_t * v_pol_j) / (n_j * ch_j * (B_p ** 2.0))
    return pol_term


def calc_radial_E_field_term(n_j, m_j, ch_j, nu_jk, nu_dj, Er, B_p):
    Er_term = (n_j * m_j * (nu_jk + nu_dj) * Er) / (n_j * ch_j * (B_p ** 2.0))
    return Er_term


def calc_toroidal_term(n_j, m_j, ch_j, nu_jk, B_p, v_tor_k):
    tor_term = (-n_j * m_j * nu_jk * v_tor_k) / (n_j * ch_j * B_p)
    return tor_term


def calc_pinch_velocity(ext_term, pol_term, Er_term, tor_term):
    vr_pinch = ext_term + pol_term + Er_term + tor_term
    return vr_pinch


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


class RadialTransport(Chi):

    def __init__(self, core, iol, nbi, iolFlag=True, neutFlag=True, debugFlag=False):
        sys.dont_write_bytecode = True

        ##############################################################
        # prep quantities for 1D transport analysis
        ##############################################################

        # prepare beams object
        beam_D = nbi.beams.D1  # Piper changes: Changed beam object references to match beamdep.
        beam_D2 = nbi.beams.D2
        beam_D3 = nbi.beams.D3
        corePatch(core, neutFlag)  # Patch to update values not brought in via ffiles (ni, zeff)
        neutPatch(core)
        dn_dr = core.dn_dr_fsa

        # prepare core and iol quantities
        r = core.r.T[0]  # TODO: Should this be a flux surface average?
        self.rhor = core.r[:, 0] / core.a

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
        Er = core.E_r_fsa  # * 1000.0 # Piper Changes: Convert input Er from kV/m to V/m
        Lp = core.Lp_fsa
        dp_dr = core.dp_dr_fsa

        # prepare iol quantities
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
        part_src_nbi_D, part_src_nbi_D_tot, part_src_nbi_D_lost, part_src_nbi_D_kept = calc_part_src_nbi(beam_D,
                                                                                                         iol_adjusted=iolFlag,
                                                                                                         F_orb_nbi=F_orb_d_nbi)
        part_src_nbi_D2, part_src_nbi_D2_tot, part_src_nbi_D2_lost, part_src_nbi_D2_kept = calc_part_src_nbi(beam_D2,
                                                                                                             iol_adjusted=iolFlag,
                                                                                                             F_orb_nbi=F_orb_d_nbi)
        part_src_nbi_D3, part_src_nbi_D3_tot, part_src_nbi_D3_lost, part_src_nbi_D3_kept = calc_part_src_nbi(beam_D3,
                                                                                                             iol_adjusted=iolFlag,
                                                                                                             F_orb_nbi=F_orb_d_nbi)

        # Piper changes: Calculate beam totals for sources, sinks, and totals.
        part_src_nbi = part_src_nbi_D + part_src_nbi_D2 + part_src_nbi_D3
        part_src_nbi_tot = part_src_nbi_D_tot + part_src_nbi_D2_tot + part_src_nbi_D3_tot
        part_src_nbi_lost = part_src_nbi_D_lost + part_src_nbi_D2_lost + part_src_nbi_D3_lost
        part_src_nbi_kept = part_src_nbi_D_kept + part_src_nbi_D2_kept + part_src_nbi_D3_kept

        # Piper changes: Changed names of particle and heat flux so it's easier to tell what method is used.
        self.gamma_diff_D = calc_gamma_diff_method(r, core.a, part_src_nbi, part_src_nbi_lost, izn_rate, core.dVdrho,
                                                   iol_adjusted=iolFlag,
                                                   F_orb=F_orb_d)  # Differential Cylindrical Method
        self.gamma_int_D = calc_gamma_int_method(r, part_src_nbi_tot, part_src_nbi_lost, izn_rate, iol_adjusted=iolFlag,
                                                 F_orb=F_orb_d)  # Integral Cylindrical Method
        self.gamma_C = np.zeros(self.gamma_int_D.shape)

        # Piper changes: Calculate radial return current (Uses integral cylindrical gamma)
        self.jr_iol = calc_return_cur(r, part_src_nbi_lost, self.gamma_int_D, izn_rate, ch_d, iol_adjusted=iolFlag,
                                      F_orb=F_orb_d)
        self.Er_iol, self.iol_term, self.diamag_term, self.diamag_term_orig, self.neut_dens_term = calc_Er_iol(n.i, n.e,
                                                                                                               m_d, n.n,
                                                                                                               B_t,
                                                                                                               Lp.i,
                                                                                                               dp_dr.i,
                                                                                                               e * z_d,
                                                                                                               T.i,
                                                                                                               dn_dr,
                                                                                                               izn_rate,
                                                                                                               self.jr_iol)

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
                                                 self.gamma_C)  # Piper Changes: Uses integral cylindrical gamma
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
            self.vpol_D, self.vpol_D_assum, self.vpol_D_alt = calc_vpol(Er, self.vtor_D_total, Lp, T, n, z_d, B_t, B_p,
                                                                        self.vtor_C_total, self.vpol_C, z_c)
        except:
            self.vpol_D = self.vpol_C / 0.4
            print 'could not calculate deuterium poloidal rotation'
            pass

        # Nick Changes: TEMPORARY - Calculate Er using pressure gradient vs. scale length.
        self.Er_calc_D, self.Er_pres_term_D, self.Er_vxb_term_D = calc_Er_mom_bal(n.i, e * z_d, dp_dr.i, T.i, Lp.i,
                                                                                  self.vtor_D_total, self.vpol_D, B_t,
                                                                                  B_p)
        self.Er_calc_C, self.Er_pres_term_C, self.Er_vxb_term_C = calc_Er_mom_bal(n.C, e * z_c, dp_dr.C, T.C, Lp.C,
                                                                                  self.vtor_C_total, self.vpol_C, B_t,
                                                                                  B_p)

        # calculate nu_drags
        mbal_rhs_D = calc_mbal_rhs(self.mom_src_tor_D_tot, z_d, n.i, B_p,
                                   self.gamma_int_D)  # Piper Changes: Uses integral cylindrical gamma
        mbal_rhs_C = calc_mbal_rhs(self.mom_src_tor_C_tot, z_c, n.C, B_p, self.gamma_C)

        nu_c_DC = 1 / calc_t90(m_d, m_c, z_d, z_c, n.C, T.i.J)
        nu_c_CD = 1 / calc_t90(m_c, m_d, z_c, z_d, n.i, T.i.J)

        # Piper changes: added alternate collision frequency calculation for comparison.
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
        self.vrpinch = calc_pinch_velocity(self.vrpinch_ext_term, self.vrpinch_poloidal_term, self.vrpinch_Er_term,
                                           self.vrpinch_toroidal_term)

        ##############################################################
        # energy balance
        ##############################################################

        # Piper Changes: Same as changes made to particle balance.
        en_src_nbi_i_D, en_src_nbi_i_D_tot, en_src_nbi_i_D_lost, en_src_nbi_i_D_kept = calc_en_src_nbi(beam_D,
                                                                                                       iol_adjusted=iolFlag,
                                                                                                       E_orb_nbi=E_orb_d_nbi)
        en_src_nbi_i_D2, en_src_nbi_i_D2_tot, en_src_nbi_i_D2_lost, en_src_nbi_i_D2_kept = calc_en_src_nbi(beam_D2,
                                                                                                           iol_adjusted=iolFlag,
                                                                                                           E_orb_nbi=E_orb_d_nbi)
        en_src_nbi_i_D3, en_src_nbi_i_D3_tot, en_src_nbi_i_D3_lost, en_src_nbi_i_D3_kept = calc_en_src_nbi(beam_D3,
                                                                                                           iol_adjusted=iolFlag,
                                                                                                           E_orb_nbi=E_orb_d_nbi)

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

        # en_src_nbi_e = en_src_nbi_e_D + en_src_nbi_e_D2 + en_src_nbi_e_D3
        en_src_nbi_e = en_src_nbi_i_kept

        # calculate radial heat flux. Piper changes: Separated heat flux equations into differential and integral cylindrical methods.
        self.Qi_diff = calc_Qi_diff_method(r, core.a, calc_cxcool(core, n, T), calc_qie(n, T), en_src_nbi_i_kept,
                                           en_src_nbi_i_lost, core.dVdrho, iol_adjusted=iolFlag,
                                           E_orb=E_orb_d)  # previously called qheat. Differential Method.
        self.Qi_int, self.qie = calc_Qi_int_method(r, n, T, en_src_nbi_i_kept, cool_rate, iol_adjusted=iolFlag,
                                                   E_orb=E_orb_d)  # Integral method.
        self.Qe_diff = calc_Qe_diff_method(r, core.a, en_src_nbi_e, cool_rate, calc_qie(n, T))  # Differential Method.
        self.Qe_int = calc_Qe_int_method(r, n, T, en_src_nbi_e, cool_rate)  # Integral method.

        self.chi = Chi(self, core, n, T, L, calc_nustar(nu_c_DC, core.q95, core.R0_a, self.vpol_C), reInterp=True)

        nn = namedtuple('nn', 's t tot')(
            core.n.n.s,  # slow
            core.n.n.t,  # thermal
            core.n.n.tot  # total
        )

    def plot_Q_e_int(self):

        plot = plt.figure()
        fig = plot.add_subplot(111)
        fig.set_xlabel(r'$\rho$', fontsize=20)
        fig.set_ylabel(r'$Q_e^{int} [W/m^2]$', fontsize=20)
        fig.set_title('GT3.RT electron heat flux')
        fig.scatter(self.rhor, self.Qe_int, marker='o', color='red')
        plt.figtext(.5, .01, "Electron heat flux using the integral method", wrap=True)

        #
        #
        # if debugFlag:
        #     """Beam source debugging"""
        #     nbiDebug(self.rhor, core.a, en_src_nbi_i_kept, en_src_nbi_e, part_src_nbi_kept, self.mom_src_nbi, core.dVdrho)
        #     """Particle flux debuggin"""
        #     gammaDebug(self.rhor, core.a, core.r2sa, core.dVdrho, part_src_nbi, part_src_nbi_kept, izn_rate, self.gamma_int_D, self.gamma_diff_D)
        #     """Energy flux debugging"""
        #     QDebug(self.rhor, core.a, core.r2sa, core.dVdrho, calc_qie(n, T, ion_species='D'), en_src_nbi_i, en_src_nbi_i_kept, cool_rate, E_orb_d, self.chi.Qi, T)
        #     """IOL Debugging"""
        #     IOLDebug(self.rhor, F_orb_d, E_orb_d, M_orb_d)
