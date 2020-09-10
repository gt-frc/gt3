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
from GT3 import Core, BeamDeposition
from collections import namedtuple
from scipy import constants
from math import sqrt
from scipy.interpolate import UnivariateSpline, interp1d
from GT3.RadialTransport.Functions.CorePatch import corePatch
from GT3.RadialTransport.Functions.CalcReturnCur import calc_return_cur
from GT3.RadialTransport.Functions.CalcNu import calc_nu_j_k, calc_nu_drag, calc_nustar
from GT3.RadialTransport.Functions.CalcMbalRHS import calc_mbal_rhs
from GT3.RadialTransport.Functions.CalcT90 import calc_t90
from GT3.RadialTransport.Functions.CalcQ import calc_Qe_diff_method, calc_qie, calc_Qe_int_method, calc_Qi_int_method, \
    calc_Qi_diff_method
from GT3.RadialTransport.Functions.CalcGamma import calc_gamma_diff_method, calc_gamma_int_method
from GT3.RadialTransport.Functions.CalcIntrinRot import calc_intrin_rot
from GT3.RadialTransport.Functions.CalcVTorDPert import calc_vtor_d_pert
from GT3.RadialTransport.Functions.CalcErMomBal import calc_Er_mom_bal
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


def viscCalc(data, core):
    """

    :type core: Core
    :type data: RadialTransport
    """
    fp = core.B_p_fsa / core.B_t_fsa
    ni = core.n_fsa.i
    Ti = core.T_fsa.i.J
    q = core.q[:, 0]
    R0 = core.R0_a
    vth = data.vpol_D[0]
    vtor = data.vtor_D_total
    vpol = data.vpol_D
    eps = core.a / core.R0_a
    nustar = data.nustar[0]
    geom = (eps ** (-3./2.) * nustar) / ((1 + eps ** (-3. / 2.) * nustar) * (1 + nustar))
    vtorS = 0.1
    vpolS = 0.1
    #eta0 = [a * m_d * b * c * core.R0_a * f1 for a, b, c in zip(n.i, vth, core.q[:, 0])]
    eta0 = ni * m_d * vth * q * R0 * geom
    #eta4 = [a * m_d * c * ch_d / (ch_d * abs(b)) for a, b, c in zip(n.i, core.B_t_fsa, T.i.ev)]
    eta4 = ni * m_d * Ti / (ch_d * abs(core.B_t_fsa))
    vrad = data.gamma_diff_D / ni

    # a = vtor    b = fp   c = eta0
    # d = vrad    f = vthet g = eta 4

    #return [a * (b * c * d - .5 * g * (4.0 * a + f)) - .5 * f * (c * d + g * (a + .5 * f)) for a, b, c, d, f, g in
    #        zip(data.vtor_D_total, fp, eta0, vrad, data.vpol_D, eta4)]

    res = vtor * vtorS * (eta0 * fp * vrad - eta4 * (2. * vtor + .5 * vpol))
    res = res - 0.5 * vpol * vpolS * (eta0 * vrad + eta4 * (vtor + .5 * vpol))
    return res / R0

def calc_chi_e(Qe, gamma_diff_D, gamma_C, n, L, T):
    gameltemp = 1.0 * gamma_diff_D + 6.0 * gamma_C   #OHHHH gamma electron!!! HA HA HA HA WE DON'T KNOW WHAT THE FUCK GAMMA_C IS

    return L.T.e * ((Qe / (ch_d * n.e * T.e.ev)) - 2.5 * gameltemp / n.e)

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


class RadialTransport:

    def __init__(self, core, iol, nbi, iolFlag=True, neutFlag=True):
        """

        :type nbi: BeamDeposition
        :type core: Core
        """
        sys.dont_write_bytecode = True

        ##############################################################
        # prep quantities for 1D transport analysis
        ##############################################################

        # prepare beams object
        corePatch(core, neutFlag)  # Patch to update values not brought in via ffiles (ni, zeff)
        #neutPatch(core)
        try:
            dn_dr = core.dn_dr_fsa
        except:
            dn_dr = None

        # prepare core and iol quantities
        r = core.r.T[0]  # TODO: Should this be a flux surface average?
        self.rhor = core.r[:, 0] / core.a
        """The rho vector"""
        self.core = core
        """The utilized GT3 core background plasma"""
        self.nbi = nbi
        """The utilized GT3 NBI module data"""
        self.iol = iol
        """The utilized GT3 IOL module data"""
        self.iolFlag = iolFlag


        self.izn_rate = core.izn_rate_fsa  # TODO: Should this be a flux surface average or a flux surface total?
        self.cool_rate = core.cool_rate_fsa  # TODO: Should this be a flux surface average or a flux surface total?
        n = core.n_fsa

        T = core.T_fsa
        L = core.L_fsa
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

        self.part_src_nbi = np.array(UnivariateSpline(nbi.beams_space, nbi.combined_beam_src_dens_total.Snbi)(self.rhor))
        self.part_src_nbi_tot = np.array(UnivariateSpline(nbi.beams_space, nbi.combined_beam_src_dens_total.Snbi)(self.rhor))
        self.part_src_nbi_lost = np.array(interp1d(nbi.beams_space, nbi.combined_beam_src_dens_lost.Snbi, kind="linear")(self.rhor))
        self.part_src_nbi_kept = np.array(UnivariateSpline(nbi.beams_space, nbi.combined_beam_src_dens_kept.Snbi)(self.rhor))

        # Piper changes: Changed names of particle and heat flux so it's easier to tell what method is used.
        self.gamma_diff_D = calc_gamma_diff_method(r, core.a, self.part_src_nbi, self.part_src_nbi_lost, self.izn_rate, core.dVdrho,
                                                   iol_adjusted=iolFlag,
                                                   F_orb=F_orb_d)  # Differential Cylindrical Method
        self.gamma_int_D = calc_gamma_int_method(r, self.part_src_nbi_tot, self.part_src_nbi_lost, self.izn_rate, iol_adjusted=iolFlag,
                                                 F_orb=F_orb_d)  # Integral Cylindrical Method
        self.gamma_C = np.zeros(self.gamma_int_D.shape)

        # Piper changes: Calculate radial return current (Uses integral cylindrical gamma)
        self.jr_iol = calc_return_cur(r, self.part_src_nbi_lost, self.gamma_int_D, self.izn_rate, ch_d, iol_adjusted=iolFlag,
                                      F_orb=F_orb_d)
        self.Er_iol, self.iol_term, self.diamag_term, self.diamag_term_orig, self.neut_dens_term = calc_Er_iol(n.i, n.e,
                                                                                                               m_d, n.n,
                                                                                                               B_t,
                                                                                                               Lp.i,
                                                                                                               dp_dr.i,
                                                                                                               e * z_d,
                                                                                                               T.i,
                                                                                                               dn_dr,
                                                                                                               self.izn_rate,
                                                                                                               self.jr_iol)

        ##############################################################
        # momentum balance
        ##############################################################
        self.mom_src_nbi = np.array(UnivariateSpline(nbi.beams_space, nbi.combined_beam_src_dens_total.Mnbi)(self.rhor))
        self.mom_src_nbi_tot = np.array(UnivariateSpline(nbi.beams_space, nbi.combined_beam_src_dens_total.Mnbi)(self.rhor))
        self.mom_src_nbi_lost = np.array(interp1d(nbi.beams_space, nbi.combined_beam_src_dens_lost.Mnbi, kind="linear")(self.rhor))
        self.mom_src_nbi_kept = np.array(UnivariateSpline(nbi.beams_space, nbi.combined_beam_src_dens_kept.Mnbi)(self.rhor))

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
        self.nu_c_j_k = calc_nu_j_k(m_d, m_c, z_d, z_c, T.i.ev, n.C)
        self.nu_c_k_j = calc_nu_j_k(m_c, m_d, z_c, z_d, T.C.ev, n.i)
        self.nu_c_j_j = calc_nu_j_k(m_d, m_d, z_d, z_d, T.i.ev, n.i)

        self.nu_drag_D = calc_nu_drag(n.i, m_d, self.vtor_D_total, self.vtor_C_total, mbal_rhs_D, nu_c_DC)
        self.nu_drag_C = calc_nu_drag(n.i, m_d, self.vtor_D_total, self.vtor_C_total, mbal_rhs_C, nu_c_CD)
        self.nustar = calc_nustar(self.nu_c_j_j, core.q[:, 0], core.R0_a, self.vpol_D)

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

        # I don't honestly understand why anything but densities are used anywhere. Changed to just use densities.

        part_src_nbi = np.array(UnivariateSpline(nbi.beams_space, nbi.combined_beam_src_total.Snbi)(self.rhor))
        part_src_nbi_tot = np.array(UnivariateSpline(nbi.beams_space, nbi.combined_beam_src_dens_total.Snbi)(self.rhor))
        part_src_nbi_lost = np.array(UnivariateSpline(nbi.beams_space, nbi.combined_beam_src_dens_lost.Snbi)(self.rhor))
        part_src_nbi_kept = np.array(UnivariateSpline(nbi.beams_space, nbi.combined_beam_src_dens_kept.Snbi)(self.rhor))

        self.en_src_nbi_i = 0.5 * np.array(UnivariateSpline(nbi.beams_space, nbi.combined_beam_src_dens_total.Qnbi)(self.rhor))
        self.en_src_nbi_i_tot = 0.5 * np.array(UnivariateSpline(nbi.beams_space, nbi.combined_beam_src_dens_total.Qnbi)(self.rhor))
        self.en_src_nbi_i_lost = 0.5 * np.array(interp1d(nbi.beams_space, nbi.combined_beam_src_dens_lost.Qnbi, kind="linear")(self.rhor))
        self.en_src_nbi_i_kept = 0.5 * np.array(UnivariateSpline(nbi.beams_space, nbi.combined_beam_src_dens_kept.Qnbi)(self.rhor))

        ################################################################################################################
        #
        #   NBI energy split - Currently 50:50 split ions and electrons
        #
        #   TODO: Implement accurate split
        #
        ################################################################################################################

        self.en_src_nbi_e = self.en_src_nbi_i_kept
        self.cxcool = calc_cxcool(core, n, T)
        self.qie = calc_qie(n, T, ion_species='D')

        # calculate radial heat flux. Piper changes: Separated heat flux equations into differential and integral cylindrical methods.
        self.Qi_diff = calc_Qi_diff_method(r, core.a, self.cxcool, self.qie, self.en_src_nbi_i_tot,
                                           self.en_src_nbi_i_lost, core.dVdrho, iol_adjusted=iolFlag,
                                           E_orb=E_orb_d)  # previously called qheat. Differential Method.

        self.Qi_int = calc_Qi_int_method(r, n, T, self.qie, self.en_src_nbi_i_kept, self.cool_rate, iol_adjusted=iolFlag,
                                                   E_orb=E_orb_d)  # Integral method.
        self.Qe_diff = calc_Qe_diff_method(r, core.a, self.en_src_nbi_e, self.cool_rate, calc_qie(n, T))  # Differential Method.

        self.Qe_int = calc_Qe_int_method(r, n, T, self.en_src_nbi_e, self.cool_rate)  # Integral method.

        self.conv15 = UnivariateSpline(core.r[:, 0], 3. * .5 * ch_d * self.gamma_diff_D * T.i.ev, k=3, s=0)(core.r[:, 0])
        self.conv25 = UnivariateSpline(core.r[:, 0], 5. * .5 * ch_d * self.gamma_diff_D * T.i.ev, k=3, s=0)(core.r[:, 0])

        self.heatvisc = viscCalc(self, core)

        self.heatin = UnivariateSpline(core.r[:, 0],
                                       .5 * self.gamma_diff_D * m_d * (self.vtor_D_total ** 2 + self.vpol_D ** 2), k=3,
                                       s=0)(core.r[:, 0])  # TODO: Provide logic that uses vtor_D_intrin/fluid depending on IOL Switch, currently too small to matter

        self.chi = namedtuple('chi', 'i e')(
            namedtuple('i', 'chi1 chi2 chi3 chi4')(
                (self.Qi_diff) * L.T.i / (n.i * T.i.ev * ch_d),
                (self.Qi_diff - self.conv25) * L.T.i / (n.i * T.i.ev * ch_d),
                (self.Qi_diff - self.conv25 - self.heatin) * L.T.i / (n.i * T.i.ev * ch_d),
                (self.Qi_diff - self.conv25 - self.heatin - self.heatvisc) * L.T.i / (n.i * T.i.ev * ch_d)
            ), calc_chi_e(self.Qe_diff, self.gamma_diff_D, self.gamma_C, n, L, T)
        )

        self.D_i = m_d * core.T_fsa.i.J * (self.nu_c_j_k * (1. - ch_d / ch_c) + self.nu_drag_D) / ((ch_d * core.B_p_fsa)**2)

    def plot_chi_terms(self, edge=True):
        fig = self._plot_base(self.conv25, title="", yLabel="q[W/m^2]", edge=edge)
        fig.scatter(self.rhor, self.heatin, color="blue", s=6)
        fig.scatter(self.rhor, self.heatvisc, color="purple", s=6)
        fig.scatter(self.rhor, self.Qi_diff, color="black", s=6)
        #fig.legend([r"$q^{conv}$", r"$q^{heatin}$", r"$q^{tot}$"])
        fig.legend([r"$q^{conv}$", r"$q^{heatin}$", r"$q^{visc}$", r"$q^{tot}$"])
        return fig

    def plot_gamma(self, edge=True):
        fig = self._plot_base(self.gamma_int_D, title="Radial particle flux", yLabel=r"$\Gamma[W/m^3]$", edge=edge)
        fig.scatter(self.rhor, self.gamma_diff_D, color="blue", s=6)
        fig.legend([r"$\Gamma_{r, int}$", r"$\Gamma_{r, diff}$"])
        return fig

    def plot_gamma_diff(self, edge=True):
        fig = self._plot_base(self.gamma_diff_D, title="Radial particle flux", yLabel=r"$\Gamma[W/m^3]$", edge=edge)
        return fig

    def plot_gamma_int(self, edge=True):
        fig = self._plot_base(self.gamma_int_D, title="Radial particle flux", yLabel=r"$\Gamma[W/m^3]$", edge=edge)
        return fig

    def plot_gamma_diff_calc(self):
        calc_gamma_diff_method(self.rhor * self.core.a,
                               self.core.a,
                               self.part_src_nbi,
                               self.part_src_nbi_lost,
                               self.izn_rate,
                               self.core.dVdrho,
                               iol_adjusted=self.iolFlag,
                               F_orb=self.iol.forb_d_therm_1D,
                               verbose=True)

    def plot_Q_i(self, edge=True):
        fig = self._plot_base(self.Qi_int, yLabel=r'$Q_i [W/m^2]$', title="Ion heat flux", edge=edge)
        fig.scatter(self.rhor, self.Qi_diff, color="black")
        fig.legend([r"$Q^{int}_i$", r"$Q^{diff}_i$"])
        return fig

    def plot_Q_e_int(self, edge=True):
        return self._plot_base(self.Qe_int, yLabel=r'$Q_e^{int} [W/m^2]$', title="Electron heat flux", edge=edge)

    def plot_q_diff_calc(self):
        calc_Qi_diff_method(self.rhor * self.core.a,
                            self.core.a,
                            self.cxcool,
                            self.qie,
                            self.en_src_nbi_i_tot,
                            self.en_src_nbi_i_lost,
                            self.core.dVdrho,
                            iol_adjusted=self.iolFlag,
                            E_orb=self.iol.eorb_d_therm_1D,
                            verbose=True)

    def _plot_base(self, val, xLabel=r'$\rho$', yLabel="Value", title="Title", color='red', edge=False):

        plot = plt.figure()
        fig = plot.add_subplot(111)
        fig.set_xlabel(xLabel, fontsize=30)
        fig.set_ylabel(yLabel, fontsize=30)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        fig.set_title(title)
        if edge:
            fig.set_xlim(0.85, 1.0)
        fig.scatter(self.rhor, val, color=color, s=8)
        plt.show()
        return fig

    def plot_ni(self, edge=True):
        return self._plot_base(self.core.n_fsa.i, yLabel=r'$n_i[\#/m^3]$', title="Ion Density", edge=edge)

    def plot_ne(self, edge=True):
        return self._plot_base(self.core.n_fsa.e, yLabel=r'$n_e[\#/m^3]$', title="Electron Density", edge=edge)

    def plot_n(self, edge=True):
        fig = self._plot_base(self.core.n_fsa.e, yLabel=r'$n_e[\#/m^3]$', title="Plasma Densities", edge=edge)
        fig.scatter(self.rhor, self.core.n_fsa.i, color="blue", s=4)
        return fig

    def plot_nn_therm(self, edge=True):
        return self._plot_base(self.core.n_fsa.n.t, yLabel=r'$n_{n,therm}[\#/m^3]$', title="Thermal Neutral Density", edge=edge)

    def plot_nn_cold(self, edge=True):
        return self._plot_base(self.core.n_fsa.n.s, yLabel=r'$n_{n,slow}[\#/m^3]$', title="Cold Neutral Density", edge=edge)

    def plot_nn_total(self, edge=True):
        return self._plot_base(self.core.n_fsa.n.tot, yLabel=r'$n_{n,total}[\#/m^3]$', title="Total Neutral Density", edge=edge)

    def plot_Ti(self, edge=True):
        return self._plot_base(self.core.T_fsa.i.kev, yLabel=r'$T_i[keV]$', title="Ion Temperature", edge=edge)

    def plot_Te(self, edge=True):
        return self._plot_base(self.core.T_fsa.e.kev, yLabel=r'$T_e[keV]$', title="Electron Temperature", edge=edge)

    def plot_T(self, edge=True):
        fig = self._plot_base(self.core.T_fsa.e.kev, yLabel=r'$T[keV]$', title="Plasma Temperature", edge=edge)
        fig.scatter(self.rhor, self.core.T_fsa.i.kev,  color="blue", s=4)
        plt.show()
        return fig

    def plot_Er(self, edge=True):
        return self._plot_base(self.Er_calc_C, yLabel=r'$E_r[V/m]$', title="Radial Electric Field", edge=edge)

    def plot_S_sources(self, edge=True, logPlot=True):
        fig = self._plot_base(self.part_src_nbi, yLabel=r'$S_r[#/m^3s]$', title="Radial sources", edge=edge)
        fig.scatter(self.rhor, self.izn_rate, color="green", s=4)
        if logPlot:
            fig.set_yscale("log")

        fig.legend([r"$S_{nbi}$", r"$S_{izn}$"])
        plt.show()
        return fig

    def plot_Q_sources(self, edge=True, logPlot=False):
        fig = self._plot_base(self.en_src_nbi_i, yLabel=r'$Q_r[W/m^3]$', title="Radial sources", edge=edge)
        fig.scatter(self.rhor, self.cool_rate, color="green", s=4)
        fig.scatter(self.rhor, self.qie, color="black", s=4)
        if logPlot:
            fig.set_yscale("log")
        plt.show()
        fig.legend(["Q_{nbi}", "Q_{cxcool}", "Q_{ie}"])
        return fig

    def plot_Chi_i_comp(self, edge=True):
        fig = self._plot_base(self.chi.i.chi1, yLabel=r'$\chi_{r,i}$', title="", edge=edge)
        fig.scatter(self.rhor, self.chi.i.chi2, color="blue", s=8)
        fig.scatter(self.rhor, self.chi.i.chi3, color="green", s=8)
        fig.scatter(self.rhor, self.chi.i.chi4, color="purple", s=4)
        fig.legend([r"$q^{cond} = q^{tot}$",
                    r"$q^{cond} = q^{tot}-q^{conv}$",
                    r"$q^{cond} = q^{tot}-q^{conv}-q^{heatin}$",
                    r"$q^{cond} = q^{tot}-q^{conv}-q^{heatin}-q^{visc}$"],  prop={'size': 20})
        return fig

    def plot_D(self, edge=True):
        fig = self._plot_base(self.D_i, yLabel=r'$D_{r, i} [m^2/s]$', title='', edge=edge)
        return fig