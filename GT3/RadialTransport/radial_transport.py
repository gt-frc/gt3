#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 08:58:46 2018

@author: max
"""


import numpy as np
import matplotlib.pyplot as plt
import sys
from math import exp
from collections import namedtuple
from scipy.interpolate import UnivariateSpline, interp1d

from GT3.RadialTransport.Functions.CalcReturnCur import calc_return_cur
from GT3.RadialTransport.Functions.CalcNu import calc_nu_j_k, calc_nu_drag, calc_nustar
from GT3.RadialTransport.Functions.CalcMbalRHS import calc_mbal_rhs
from GT3.RadialTransport.Functions.CalcT90 import calc_t90
from GT3.RadialTransport.Functions.CalcQ import calc_qie
from GT3.RadialTransport.Functions.CalcIntrinRot import calc_intrin_rot
from GT3.RadialTransport.Functions.CalcVTorDPert import calc_vtor_d_pert
from GT3.RadialTransport.Functions.CalcErMomBal import calc_Er_mom_bal
from GT3.RadialTransport.Functions.CalcErIOL import calc_Er_iol
from GT3.RadialTransport.Functions.CalcVpol import calc_vpol
from GT3.RadialTransport.Functions.CalcCXCool import calc_cxcool
import GT3.constants as constants
from GT3.Core.Functions.ProfileClasses import OneDProfile, TemperatureProfiles, DensityProfiles, PressureProfiles, Flux
from GT3.utilities.PlotBase import PlotBase
from GT3 import Core, BeamDeposition

eps_0 = constants.epsilon_0
e = constants.elementary_charge
m_e = constants.electron_mass
m_d = constants.deuteron_mass
m_c = constants.carbon_mass

z_d = 1  # atomic number of deuterium
z_c = 6  # atomic number of carbon
ch_d = e * z_d  # charge of deuterium
ch_c = e * z_c  # charge of carbon

E_phi = 0.04  # toroidal electrostatic potential

class RadTransTemperatureProfiles(TemperatureProfiles):

    def __init__(self, core: Core, *args, **kwargs):
        super().__init__(core.psi, core.R, core.Z, ProfileType=OneDProfile, *args, **kwargs)

class RadTransDensityProfiles(DensityProfiles):
    def __init__(self, core: Core, *args, **kwargs):
        super().__init__(core.psi, core.R, core.Z, ProfileType=OneDProfile, *args, **kwargs)

class RadTransPressureProfiles(PressureProfiles):
    def __init__(self, core: Core, *args, **kwargs):
        super().__init__(core.psi, core.R, core.Z, ProfileType=OneDProfile, *args, **kwargs)



def calc_chi_e(Qe, gamma_diff_D, gamma_C, n, T):
    gameltemp = 1.0 * gamma_diff_D + 6.0 * gamma_C   #OHHHH gamma electron!!! HA HA HA HA WE DON'T KNOW WHAT THE FUCK GAMMA_C IS
    return T.e.J.L * ((Qe / (ch_d * n.e * T.e.ev)) - 2.5 * gameltemp / n.e.val)

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

class RadialTransport(PlotBase):

    def __init__(self, core, iol, nbi: BeamDeposition, iolFlag=True, neutFlag=True):
        """

        :type nbi: BeamDeposition
        :type core: Core
        """
        super().__init__()
        sys.dont_write_bytecode = True

        ##############################################################
        # prep quantities for 1D transport analysis
        ##############################################################

        # prepare beams object
        #corePatch(core, neutFlag)  # Patch to update values not brought in via ffiles (ni, zeff)

        # prepare core and iol quantities
        r = core.r.T[0]  # TODO: Should this be a flux surface average?
        self.rhor = core.r[:, 0] / core.a
        self.set_plot_rho1d(self.rhor)
        """The rho vector"""
        self.core = core
        """The utilized GT3 core background plasma"""
        self.nbi = nbi
        """The utilized GT3 NBI module data"""
        self.iol = iol
        """The utilized GT3 IOL module data"""
        self.iolFlag = iolFlag


        self.izn_rate = core.izn_rate.tot.fsa  # TODO: Should this be a flux surface average or a flux surface total?
        self.cool_rate = core.cool_rate.fsa  # TODO: Should this be a flux surface average or a flux surface total?

        n = RadTransDensityProfiles(self.core,
                                    i=self.core.n.i.fsa,
                                    e=self.core.n.e.fsa,
                                    C=self.core.n.C.fsa,
                                    ns=self.core.n.n.s.fsa,
                                    nt=self.core.n.n.t.fsa)


        T = RadTransTemperatureProfiles(self.core,
                                        i=self.core.T.i.kev.fsa,
                                        e=self.core.T.e.kev.fsa,
                                        C=self.core.T.C.kev.fsa)

        p = RadTransPressureProfiles(self.core,
                                     i=self.core.p.i.fsa,
                                     e=self.core.p.e.fsa,
                                     C=self.core.p.C.fsa)


        B_p = core.B.pol.fsa
        B_t = core.B.tor.fsa
        Er = core.E_r.fsa  # * 1000.0 # Piper Changes: Convert input Er from kV/m to V/m

        # Put some information in instances for debugging purposes. The use of n, T, p, etc. makes writing equations
        # easier.
        self._n = n
        self._T = T
        self._p = p
        self._Bp = B_p
        self._Bt = B_t
        self._Er = Er


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

        self.part_src_nbi = nbi.combined_beam_src_dens_total.Snbi
        self.part_src_nbi_tot = nbi.combined_beam_src_dens_total.Snbi
        self.part_src_nbi_lost = nbi.combined_beam_src_dens_lost.Snbi
        self.part_src_nbi_kept = nbi.combined_beam_src_dens_kept.Snbi

        # Piper changes: Changed names of particle and heat flux so it's easier to tell what method is used.
        gamma_diff_D = self._calc_gamma_diff_method(iol_adjusted=iolFlag, F_orb=F_orb_d,
                                                    neutFlag=neutFlag)  # Differential Cylindrical Method
        gamma_int_D = self._calc_gamma_int_method(r, iol_adjusted=iolFlag, F_orb=F_orb_d,
                                                  neutFlag=neutFlag)  # Integral Cylindrical Method

        self.gamma = Flux(core,
                          D_diff=gamma_diff_D,
                          D_int=gamma_int_D,
                          C_diff=np.zeros(gamma_int_D.shape),
                          C_int=np.zeros(gamma_int_D.shape))

        # Piper changes: Calculate radial return current (Uses integral cylindrical gamma)
        self.jr_iol = calc_return_cur(r, self.part_src_nbi_lost, self.gamma.D.int, ch_d, iol_adjusted=iolFlag,
                                      F_orb=F_orb_d)
        self.Er_iol, self.iol_term, self.diamag_term, self.diamag_term_orig, self.neut_dens_term = calc_Er_iol(n.i, n.e,
                                                                                                               m_d, n.n,
                                                                                                               B_t,
                                                                                                               p,
                                                                                                               e * z_d,
                                                                                                               T.i,
                                                                                                               n.n.tot.derivative(),
                                                                                                               self.izn_rate,
                                                                                                               self.jr_iol)

        ##############################################################
        # momentum balance
        ##############################################################
        self.mom_src_nbi = nbi.combined_beam_src_dens_total.Mnbi
        self.mom_src_nbi_tot = nbi.combined_beam_src_dens_total.Mnbi
        self.mom_src_nbi_lost = nbi.combined_beam_src_dens_lost.Mnbi
        self.mom_src_nbi_kept = nbi.combined_beam_src_dens_kept.Mnbi

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
        self.vtor_C_total = core.v.C.tor.fsa
        self.vtor_C_fluid = self.vtor_C_total - self.vtor_C_intrin

        # calculate deuterium toroidal rotation
        self.vtor_D_intrin = calc_intrin_rot(M_orb_d, T.i.J, m_d)

        # Piper Changes: Changed core.v_1D.tor.C.any() to core.v_1D.tor.D.any(). Carbon velocity should be a given.
        if not core.v.i.tor.isNonZero():  # if array is all zeros, then no input. Use perturbation theory.
            vtor_D_total = calc_vtor_d_pert(self.vtor_C_total,
                                                 self.vtor_C_intrin,
                                                 self.vtor_D_intrin,
                                                 self.mom_src_tor_D_tot,
                                                 1,
                                                 n,
                                                 T,
                                                 B_p,
                                                 self.gamma.D.int,
                                                 self.gamma.C.int)  # Piper Changes: Uses integral cylindrical gamma
            self.core.v.D.tor.update_from_1D(vtor_D_total)

            # Broadcast to 2D before replacing
            vtor_D_prep = np.broadcast_to(vtor_D_total, (self.core.rho.shape[1], len(vtor_D_total))).T

            # Update TwoDProfile
            self.core.v.update_D(tor = vtor_D_prep)
            # Piper changes: added a message to let the user know the D velocity was calculated.
            print('Deuterium toroidal velocity calculated from perturbation theory.')
            self.vtor_D_total = OneDProfile(core.psi, vtor_D_total, core.R, core.Z)
        else:
            # Piper Changes: For some reason this used to set D velocity to C velocity,
            # which overwrote the input D velocity.
            self.vtor_D_total = core.v.i.tor.fsa

        self.vtor_fluid_D = self.vtor_D_total - self.vtor_D_intrin

        # calculate carbon and deuterium poloidal rotation
        try:
            self.vpol_C = core.v.C.pol.fsa
            vpol_D, vpol_D_assum, vpol_D_alt = calc_vpol(Er, self.vtor_D_total, p, T, n, z_d, B_t, B_p,
                                                                        self.vtor_C_total, self.vpol_C, z_c)
            self.vpol_D = OneDProfile(self.core.psi, vpol_D, self.core.R, self.core.Z)
            self.vpol_D_assum = vpol_D_assum
            self.vpol_D_alt = vpol_D_alt
        except:
            self.vpol_D = OneDProfile(self.core.psi, self.vpol_C.val / 0.4, self.core.R, self.core.Z)
            print('could not calculate deuterium poloidal rotation')
            pass

        # Update core velocities
        # Broadcast to 2D before replacing
        vpol_D_prep = np.broadcast_to(self.vpol_D, (self.core.rho.shape[1], len(self.vpol_D))).T

        # Update TwoDProfile
        self.core.v.update_D(pol=vpol_D_prep)

        # Nick Changes: TEMPORARY - Calculate Er using pressure gradient vs. scale length.
        self.Er_calc_D, self.Er_pres_term_D, self.Er_vxb_term_D = calc_Er_mom_bal(n.i, e * z_d, p.i.derivative(),
                                                                                  self.vtor_D_total, self.vpol_D, B_t,
                                                                                  B_p)
        self.Er_calc_C, self.Er_pres_term_C, self.Er_vxb_term_C = calc_Er_mom_bal(n.C, e * z_c, p.C.derivative(),
                                                                                  self.vtor_C_total, self.vpol_C, B_t,
                                                                                  B_p)

        # calculate nu_drags
        #mbal_rhs_D = calc_mbal_rhs(self.mom_src_tor_D_tot, z_d, n.i, B_p,
        #                           self.gamma_int_D)  # Piper Changes: Uses integral cylindrical gamma
        mbal_rhs_D = calc_mbal_rhs(self.mom_src_tor_D_tot, z_d, n.i, B_p,
                                   self.gamma.D.diff)  # Piper Changes: Uses integral cylindrical gamma
        mbal_rhs_C = calc_mbal_rhs(self.mom_src_tor_C_tot, z_c, n.C, B_p, self.gamma.C.int)

        nu_c_DC = 1 / calc_t90(m_d, m_c, z_d, z_c, n.C, T.i.J)
        nu_c_CD = 1 / calc_t90(m_c, m_d, z_c, z_d, n.i, T.i.J)

        # Piper changes: added alternate collision frequency calculation for comparison.
        self.nu_c_j_k = calc_nu_j_k(m_d, m_c, z_d, z_c, T.i.ev, n.C)
        self.nu_c_k_j = calc_nu_j_k(m_c, m_d, z_c, z_d, T.C.ev, n.i)
        self.nu_c_j_j = calc_nu_j_k(m_d, m_d, z_d, z_d, T.i.ev, n.i)
        self.nu_c_j_e = calc_nu_j_k(m_d, m_e, z_d, z_d, T.i.ev, n.e)
        self.nu_c_e_j = calc_nu_j_k(m_e, m_d, z_d, z_d, T.e.ev, n.i)
        self.nu_c_e_e = calc_nu_j_k(m_e, m_e, z_d, z_d, T.e.ev, n.e)

        self.nu_drag_D = calc_nu_drag(n.i, m_d, self.vtor_D_total, self.vtor_C_total, mbal_rhs_D, nu_c_DC)
        self.nu_drag_C = calc_nu_drag(n.i, m_d, self.vtor_D_total, self.vtor_C_total, mbal_rhs_C, nu_c_CD)
        self.nustar = calc_nustar(self.nu_c_j_j, core.q.fsa, core.R0_a, self.vpol_D)

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

        self.en_src_nbi_i = OneDProfile(core.psi, 0.5 * nbi.combined_beam_src_dens_total.Qnbi, core.R, core.Z)
        self.en_src_nbi_i_tot = OneDProfile(core.psi, 0.5 * nbi.combined_beam_src_dens_total.Qnbi, core.R, core.Z)
        self.en_src_nbi_i_lost = OneDProfile(core.psi, 0.5 * nbi.combined_beam_src_dens_lost.Qnbi, core.R, core.Z)
        self.en_src_nbi_i_kept = OneDProfile(core.psi, 0.5 * nbi.combined_beam_src_dens_kept.Qnbi, core.R, core.Z)

        ################################################################################################################
        #
        #   NBI energy split - Currently 50:50 split ions and electrons
        #
        #   TODO: Implement accurate split
        #
        ################################################################################################################

        self.en_src_nbi_e = OneDProfile(core.psi, self.en_src_nbi_i_kept, core.R, core.Z)
        self.cxcool = OneDProfile(core.psi, calc_cxcool(core, n, T), core.R, core.Z)
        self.qie = OneDProfile(core.psi, calc_qie(n, T, ion_species='D'), core.R, core.Z)

        # calculate radial heat flux. Piper changes: Separated heat flux equations into differential and integral cylindrical methods.
        Qi_diff = self._calc_Qi_diff_method(iol_adjusted=iolFlag, E_orb=E_orb_d)  # previously called qheat. Differential Method.

        Qi_int = self._calc_Qi_int_method(iol_adjusted=iolFlag, E_orb=E_orb_d)  # Integral method.
        Qe_diff = self._calc_Qe_diff_method(self.cool_rate, calc_qie(n, T))  # Differential Method.

        Qe_int = self._calc_Qe_int_method()  # Integral method.

        self.Q = Flux(core, label=r"$Q_r",
                      D_int=Qi_int,
                      D_diff=Qi_diff,
                      e_int=Qe_int,
                      e_diff=Qe_diff)

        conv15 = 3. * .5 * ch_d * self.gamma.D.diff * T.i.ev
        conv25 = 5. * .5 * ch_d * self.gamma.D.diff * T.i.ev
        hvisc = self._calc_visc_heat()
        heatin = .5 * self.gamma.D.diff * m_d * (self.vtor_D_total ** 2 + self.vpol_D ** 2) # TODO: Provide logic that uses vtor_D_intrin/fluid depending on IOL Switch, currently too small to matter

        self.conv15 = OneDProfile(self.core.psi, conv15, self.core.R, self.core.Z)
        self.conv25 = OneDProfile(self.core.psi, conv25, self.core.R, self.core.Z)
        self.heatvisc = OneDProfile(self.core.psi, hvisc, self.core.R, self.core.Z)
        self.heatin = OneDProfile(self.core.psi, heatin, self.core.R, self.core.Z)

        self.chi = namedtuple('chi', 'i e')(
            namedtuple('i', 'chi1 chi2 chi3 chi4')(
                (self.Q.D.diff) * T.i.J.L / (n.i * T.i.ev * ch_d),
                (self.Q.D.diff - self.conv25) * T.i.kev.L / (n.i * T.i.ev * ch_d),
                (self.Q.D.diff - self.conv25 - self.heatin) * T.i.J.L / (n.i * T.i.ev * ch_d),
                self._calc_chi_i_visc()
            ), calc_chi_e(self.Q.e.diff, self.gamma.D.diff, self.gamma.C.diff, n, T)
        )

        D_i = m_d * T.i.J * (self.nu_c_j_k * (1. - ch_d / ch_c) + self.nu_drag_D) / ((ch_d * core.B.pol.fsa)**2)
        self.D_i = OneDProfile(self.core.psi, D_i, self.core.R, self.core.Z)

    def _calc_chi_i_visc(self, vtorS=0.1, vpolS=0.1):
        heatvis = OneDProfile(self.core.psi, self._calc_visc_heat(vtorS, vpolS), self.core.R, self.core.Z)
        return (self.Q.D.diff - self.conv25 - self.heatin - heatvis) * self._T.i.J.L / (self._n.i * self._T.i.ev * ch_d)

    def _calc_gamma_diff_method(self, iol_adjusted=False, F_orb=None, neutFlag=True, verbose=False, *args, **kwargs):
        a = self.core.a
        r = self.rhor * self.core.a
        dF_orb = UnivariateSpline(r, F_orb, k=3, s=0).derivative()
        izn_rateint = UnivariateSpline(r, self.izn_rate.val, k=2, s=0)
        part_src_nbi_totint = UnivariateSpline(r, self.part_src_nbi_tot.val, k=2, s=0)
        part_src_nbi_lostint = UnivariateSpline(r, self.part_src_nbi_lost.val, k=2, s=0)
        iolPeak = np.where(dF_orb(r) == dF_orb(r).max())

        def f(t, gamma, sion, snbi, snbi_loss, dFdr, iolFlag, peak):
            if neutFlag:
                # if t/a >= 0.95:
                # S = snbi(t) + sion(t)
                # else:
                # S = snbi(t)
                S = snbi(t) + sion(t)
            else:
                S = snbi(t)
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
        gamma.set_initial_value(0., 0.).set_f_params(izn_rateint, part_src_nbi_totint, part_src_nbi_lostint, dF_orb,
                                                     iol_adjusted, r[iolPeak])
        dt = a / len(r)
        x, y = [], []
        while gamma.successful() and gamma.t < a:
            x.append(gamma.t + dt)
            y.append(gamma.integrate(gamma.t + dt))
        # gamma = UnivariateSpline(x, y, k=3, s=0)
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

        if kwargs.get("splineVerify"):
            plot1 = plt.figure()
            #plot1.set_title("Spline fit verification")
            fig1 = plot1.add_subplot(411)
            fig2 = plot1.add_subplot(412)
            fig3 = plot1.add_subplot(413)
            fig4 = plot1.add_subplot(414)

            fig1.scatter(r, self.izn_rate.val, color="red", marker="x")
            fig1.plot(r, izn_rateint(r))
            fig1.set_title("izn_rate")
            fig2.scatter(r, self.part_src_nbi_tot.val, color="red", marker="x")
            fig2.plot(r, part_src_nbi_totint(r))
            fig2.set_title("nbi_tot")
            fig3.scatter(r, self.part_src_nbi_lost.val, color="red", marker="x")
            fig3.plot(r, part_src_nbi_lostint(r))
            fig3.set_title("nbi_list")
            fig4.scatter(x, y, color="red", marker="x")
            fig4.plot(r, gamma(r))
            fig4.set_title("gamma")

        return gamma(r)

    def _calc_gamma_int_method(self, r, iol_adjusted=False, F_orb=None, neutFlag=True):
        # Piper Changes: Added cylindrical integral method as a separate function. This will be set to a separate variable in the main code.
        gamma = np.zeros(r.shape)

        if not neutFlag:
            izn_rate = [0.] * len(self.izn_rate.val)
        else:
            izn_rate = self.izn_rate

        # Boundary condition at magnetic axis. Needs to be in units of ions/m^3.
        # Only has the second term, since it's the center value. Also uses delta_r of the next point to avoid indexing into a non-existant location.
        # If not adjusted for IOL, part_src_nbi_lost = 0 anyway, so no need for an IOL check.
        gamma[0] = (self.part_src_nbi_tot[0] - 2 * self.part_src_nbi_lost[0] + izn_rate[0]) * (r[1] - r[0])
        gamma[1] = gamma[0] + (self.part_src_nbi_tot[1] - 2 * self.part_src_nbi_lost[1] + izn_rate[1]) * (r[1] - r[0])

        # You'll  prolly want to change this since it uses a dreaded for loop.
        for n in range(2, len(r)):  # Pretty sure len() is still valid for multidimensional arrays.
            # The 2*part_src_nbi_lost is the factor of 2 in the fast IOL.
            if iol_adjusted:
                # Imported the exp() function for the thermal IOL attenuation.
                gamma[n] = (r[n - 1] / r[n]) * gamma[n - 1] * exp(-2 * (F_orb[n] - F_orb[n - 1])) + (
                            self.part_src_nbi_tot[n] - 2. * self.part_src_nbi_lost[n] + izn_rate[n]) * (r[n] - r[n - 1])
            else:
                gamma[n] = (r[n - 1] / r[n]) * gamma[n - 1] + (self.part_src_nbi_tot[n] + izn_rate[n]) * (r[n] - r[n - 1])

        return gamma

    def _calc_Qe_diff_method(self, cool_rate, Qie):
        a = self.core.a
        r = self.rhor * a
        en_src_nbi_e = self.en_src_nbi_e
        en_src_nbi_eint = UnivariateSpline(r, en_src_nbi_e, k=3, s=0)
        cool_rateint = UnivariateSpline(r, cool_rate.val, k=3, s=0)
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

    def _calc_Qe_int_method(self):  # Piper Changes: Same as Qi changes.
        r = self.rhor * self.core.a
        n = self._n
        T = self._T
        cool_rate = self.cool_rate
        en_src_nbi_e_tot = self.en_src_nbi_e
        Qe = np.zeros(r.shape)
        qie = calc_qie(n, T, ion_species='D')

        Qe[0] = (en_src_nbi_e_tot[0] - qie[0] - cool_rate[0]) * (r[1] - r[0])
        Qe[1] = Qe[0] + (en_src_nbi_e_tot[1] - qie[1] - cool_rate[1]) * (r[1] - r[0])

        # Integral cylindrical form of the energy balance equation.
        # Identical in form to the continuity equation, but different source term.
        for n in range(2, len(r)):
            Qe[n] = (r[n - 1] / r[n]) * Qe[n - 1] + (en_src_nbi_e_tot[n] - qie[n] - cool_rate[n]) * (r[n] - r[n - 1])

        return Qe

    def _calc_Qi_diff_method(self, iol_adjusted=False, E_orb=None, verbose=False, *args, **kwargs):
        a = self.core.a
        r = self.rhor * a
        en_src_nbi_tot = self.en_src_nbi_i_tot
        en_src_nbi_lost = self.en_src_nbi_i_lost
        Qie = self.qie
        cxcool = self.cxcool
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
                return S - Qnbi_loss(t) - (flux * (dEdrval + 1) / (t + 0.003))
            else:
                return S - (flux * (1) / (t + 0.003))

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

        if kwargs.get("splineVerify"):
            plot1 = plt.figure()
            #plot1.set_title("Spline fit verification")
            fig1 = plot1.add_subplot(321)
            fig2 = plot1.add_subplot(322)
            fig3 = plot1.add_subplot(323)
            fig4 = plot1.add_subplot(324)
            fig5 = plot1.add_subplot(325)

            fig1.scatter(r, en_src_nbi_tot, color="red", marker="x")
            fig1.plot(r, en_src_nbi_totint(r))
            fig1.set_title("nbi_src_tot")
            fig2.scatter(r, en_src_nbi_lost, color="red", marker="x")
            fig2.plot(r, en_src_nbi_lostint(r))
            fig2.set_title("nbi_src_lost")
            fig3.scatter(r, Qie, color="red", marker="x")
            fig3.plot(r, Qie_int(r))
            fig3.set_title("Qie")
            fig4.scatter(r, cxcool, color="red", marker="x")
            fig4.plot(r, cxcoolint(r))
            fig4.set_title("CX cooling")
            fig5.scatter(x, y, color="red", marker="x")
            fig5.plot(r, flux(r))
            fig5.set_title("flux")

        return flux(r)

    def _calc_Qi_int_method(self, iol_adjusted=False, E_orb=None):  # formerly qheat
        r = self.rhor * self.core.a
        en_src_nbi_i_kept = self.en_src_nbi_i_kept
        cool_rate = self.cool_rate
        qie = self.qie
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
                Qi[n] = (r[n - 1] / r[n]) * Qi[n - 1] + (en_src_nbi_i_kept[n] - cool_rate[n] - qie[n]) * (
                            r[n] - r[n - 1])

        return Qi

    def _calc_visc_heat(self, vtorS=0.1, vpolS=0.1):
        """

        :type core: Core
        :type self: RadialTransport
        """
        fp = self.core.B.pol.fsa / self.core.B.tor.fsa
        ni = self._n.i
        Ti = self._T.i.J
        q = self.core.q.fsa
        R0 = self.core.R0_a
        vtor = self.core.v.D.tor.fsa
        vpol = self.core.v.D.pol.fsa
        vth = self.core.v.D.tot.fsa
        eps = self.core.a / self.core.R0_a
        nustar = self.nustar
        geom = (eps ** (-3. / 2.) * nustar) / ((1 + eps ** (-3. / 2.) * nustar) * (1 + nustar))
        # eta0 = [a * m_d * b * c * core.R0_a * f1 for a, b, c in zip(n.i, vth, core.q[:, 0])]
        eta0 = ni * m_d * vth * q * R0 * geom
        # eta4 = [a * m_d * c * ch_d / (ch_d * abs(b)) for a, b, c in zip(n.i, core.B_t_fsa, T.i.ev)]
        eta4 = ni * m_d * Ti / (ch_d * abs(self.core.B.tor.fsa.val))
        vrad = OneDProfile(self.core.psi, self.gamma.D.diff / ni, self.core.r, self.core.Z)

        # a = vtor    b = fp   c = eta0
        # d = vrad    f = vthet g = eta 4

        # return [a * (b * c * d - .5 * g * (4.0 * a + f)) - .5 * f * (c * d + g * (a + .5 * f)) for a, b, c, d, f, g in
        #        zip(data.vtor_D_total, fp, eta0, vrad, data.vpol_D, eta4)]

        res = vtor * vtorS * (eta0 * fp * vrad - eta4 * (2. * vtor + .5 * vpol))
        res = res - 0.5 * vpol * vpolS * (eta0 * vrad + eta4 * (vtor + .5 * vpol))
        return res / R0

    def plot_chi_terms(self, edge=True):

        fig = self._plot_base(self.conv25, title="", yLabel="q[W/m^2]", edge=edge)
        fig.scatter(self.rhor, self.heatin, color="blue", s=self._markerSize)
        fig.scatter(self.rhor, self.heatvisc, color="purple", s=self._markerSize)
        fig.scatter(self.rhor, self.Q.D.diff, color="black", s=self._markerSize)
        #fig.legend([r"$q^{conv}$", r"$q^{heatin}$", r"$q^{tot}$"])
        fig.legend([r"$q^{conv}$", r"$q^{heatin}$", r"$q^{visc}$", r"$q^{tot}$"])
        return fig


    def plot_gamma_diff_calc(self):
        self._calc_gamma_diff_method(iol_adjusted=self.iolFlag, F_orb=self.iol.forb_d_therm_1D, verbose=True)


    def plot_nu_jk(self, edge=True):
        return self._plot_base(self.nu_c_j_k, yLabel=r'$\nu_{j,k}$', title="Ion-Impurity Collision frequency", edge=edge)

    def plot_nu_kj(self, edge=True):
        return self._plot_base(self.nu_c_k_j, yLabel=r'$\nu_{k,j}$', title="Impurity-Ion Collision frequency", edge=edge)

    def plot_nu_jj(self, edge=True):
        return self._plot_base(self.nu_c_j_j, yLabel=r'$\nu_{j,j}$', title="Ion-Ion Collision frequency", edge=edge)

    def plot_nu_ee(self, edge=True):
        return self._plot_base(self.nu_c_e_e, yLabel=r'$\nu_{e,e}$', title="Electron-Electron Collision frequency", edge=edge)

    def plot_nu_je(self, edge=True):
        return self._plot_base(self.nu_c_j_e, yLabel=r'$\nu_{j,e}$', title="Ion-Electron Collision frequency", edge=edge)

    def plot_nu_ej(self, edge=True):
        return self._plot_base(self.nu_c_e_j, yLabel=r'$\nu_{e,j}$', title="Electron-Ion Collision frequency", edge=edge)

    def plot_Er(self, edge=True):
        return self._plot_base(self.Er_calc_C, yLabel=r'$E_r[V/m]$', title="Radial Electric Field", edge=edge)

    def plot_S_sources(self, edge=True, logPlot=True):
        fig = self._plot_base(self.part_src_nbi, yLabel=r'$S_r[#/m^3s]$', title="Radial sources", edge=edge)
        fig.scatter(self.rhor, self.izn_rate, color="green", s=self._markerSize)
        if logPlot:
            fig.set_yscale("log")

        fig.legend([r"$S_{nbi}$", r"$S_{izn}$"])
        plt.show()
        return fig

    def plot_Q_sources(self, edge=True, logPlot=False):
        fig = self._plot_base(self.en_src_nbi_i_kept, yLabel=r'$Q_r[W/m^3]$', title="", edge=edge)
        fig.scatter(self.rhor, self.cxcool, color="green", s=self._markerSize)
        fig.scatter(self.rhor, abs(self.qie.val), color="black", s=self._markerSize)
        if logPlot:
            fig.set_yscale("log")
        plt.show()
        fig.legend([r"$Q_{nbi}$", r"$Q_{cxcool}$", r"$|Q_{ie}|$"], prop={'size': 30}, markerscale=1.5)
        return fig

    def plot_Chi_i_comp(self, edge=True, marker=False):
        fig = self._plot_base(self.chi.i.chi1, yLabel=r'$\chi_{r,i}$', title="", edge=edge)
        if marker:
            fig.scatter(self.rhor, self.chi.i.chi2, color="blue", s=self._markerSize, marker="x")
            fig.scatter(self.rhor, self.chi.i.chi3, color="green", s=self._markerSize, marker="o", facecolors="None")
            fig.scatter(self.rhor, self.chi.i.chi4, color="purple", s=self._markerSize, marker="^")
        else:
            fig.scatter(self.rhor, self.chi.i.chi2, color="blue", s=self._markerSize)
            fig.scatter(self.rhor, self.chi.i.chi3, color="green", s=self._markerSize)
            fig.scatter(self.rhor, self.chi.i.chi4, color="purple", s=self._markerSize)
        fig.legend([r"$q^{cond} = q^{tot}$",
                    r"$q^{cond} = q^{tot}-q^{conv}$",
                    r"$q^{cond} = q^{tot}-q^{conv}-q^{heatin}$",
                    r"$q^{cond} = q^{tot}-q^{conv}-q^{heatin}-q^{visc}$"],  prop={'size': 20}, markerscale=1.5)
        return fig

    def plot_D(self, edge=True):
        fig = self._plot_base(self.D_i, yLabel=r'$D_{r, i} [m^2/s]$', title='', edge=edge)
        return fig