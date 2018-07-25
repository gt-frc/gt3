#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  5 16:05:08 2017

@author: max
"""
from __future__ import division
from math import pi, sin, acos, ceil
import numpy as np
from scipy.constants import mu_0, elementary_charge, k
from scipy.interpolate import UnivariateSpline, interp1d, interp2d, griddata
from scipy.integrate import quad
import matplotlib.pyplot as plt
import sys
from collections import namedtuple


def PolyArea(x, y):
    """finds the area of a polygon"""
    return 0.5*np.abs(np.dot(x, np.roll(y, 1))-np.dot(y, np.roll(x, 1)))


def calc_radial_prof(rho, v0, vp, vs, nu, ped_loc=0.9):
    if ped_loc == 1:  # then there is no pedestal
        prof = v0 * (1 - rho ** 2) ** nu
    else:
        prof = np.where(rho < ped_loc,
                        (v0 - vp) * (1 - rho ** 2) ** nu + vp,
                        (vs - vp) / (1 - ped_loc) * (rho - 1) + vs)
    return prof


def xmiller(kappa, rho, theta, inp):
    """
    """

    # PART 1: SELECT CONVOLUTION KERNEL. WE'LL USE A STANDARD BUMP FUNCTION.
    def bump(x, epsilon):
        """
        """

        # define eta0
        def eta0(x2):
            """
            """
            eta0 = np.where(np.logical_and((x2 > -1), (x2 < 1)), np.exp(-1.0 / (1.0 - np.power(np.abs(x2), 2.0))),
                            0.)
            return eta0

        # calculate normalization coefficient
        C = quad(eta0, -1, 1)[0]
        # calculate eta_eps
        eta_eps = 1 / epsilon / C * eta0(x / epsilon)
        return eta_eps

        # PART 2: DEFINE SEPERATRIX FUNCTION

    def kappa_sep(x, gamma1, gamma2):
        """
        """
        # Necessary kappa to get the z-value of the x-point
        kappa_x = (inp.xpt[1] - inp.Z0) / (inp.a * sin(3. * pi / 2.))
        # Amount of "extra" kappa needed at the x-point, i.e. kappa_tilda
        delta_kappa = kappa_x - inp.kappa_lo

        kappa_sep = np.piecewise(
            x,
            [x <= pi,  # condition 1
             np.logical_and((pi < x), (x <= 3. * pi / 2.)),  # condition 2
             np.logical_and((3. * pi / 2. < x), (x <= 2. * pi)),  # condition 3
             x > 2 * pi],  # condition 4
            [lambda x: 0,  # function for condition 1
             lambda x: delta_kappa * (1. - abs(2. * x / pi - 3.) ** (1.0)) ** gamma1,  # function for condition 2
             lambda x: delta_kappa * (1. - abs(2. * x / pi - 3.) ** (1.0)) ** gamma2,  # function for condition 3
             lambda x: 0]  # function for condition 4
        )
        return kappa_sep

        # PART 2: RESCALE THETA (OPTIONAL) (HOW QUICKLY TO DO YOU LIMIT THETA RANGE OF CONVLUTION, IF AT ALL)

    def rescale_theta(theta, epsilon):
        """
        """
        return (theta - 3 * pi / 2) / (1 - epsilon) + 3 * pi / 2

    # PART 3: DEFINE Y-SCALE (HOW QUICKLY DO FLUX SURFACESAPPROACH THE SEPERATRIX AS A FUNCTION OF r)
    def yscale(r, a):
        """
        """
        nu = 5.0  # nu=10ish works well
        return np.power(r / a, nu)  # * (xpt[1] / (a*sin(3.*pi/2.)) - kappa_lo)

    # PART 4: DEFINE EPSILON AS A RUNCTION OF r (HOW QUICKLY DO FLUX SURFACES GET 'POINTY' AS THEY APPROACH THE SEPERATRIX)
    def calc_eps(r, a):
        """
        """
        D = 2.0  # What is this?
        nu = 3.0
        return D * (1 - np.power(r / a, nu))

    # PART 5: POST TRANSFORM (SUBTRACT ENDPOINTS, ETC.)
    def posttransform(x, f):  # here x is from pi to 2pi
        """
        """
        f_pt = f - ((f[-1] - f[0]) / (2 * pi - pi) * (x - pi) + f[0])
        return f_pt

    # INITIALIZE KAPPA_TILDA ARRAY
    kappa_tilda = np.zeros(kappa.shape)

    # DEFINE POINTS FOR THE CONVOLUTION. THIS IS SEPARATE FROM THETA POINTS.
    xnum = 1001

    # DEFINE SEPERATRIX FUNCTION
    gamma1 = 5.0  # 3.8
    gamma2 = 2.0  # 1.0
    k_sep = kappa_sep(np.linspace(pi, 2 * pi, xnum), gamma1, gamma2)

    # For each flux surface (i.e. r value)
    for i, rhoval in enumerate(rho.T[0]):
        if rhoval < 1:
            # Calculate epsilon
            epsilon = calc_eps(rhoval, inp.a)
            # OPTIONALLY - modify domain to ensure a constant domain of compact support based on current epsilon
            scale_theta = 0
            if scale_theta == 1:
                thetamod = rescale_theta(theta, epsilon)

            # define convolution kernel for the flux surface, eta_epsilon (eta for short)
            eta = bump(np.linspace(-1, 1, xnum), epsilon)

            # scale eta. The convolution operation doesn't
            scaled_eta = eta * yscale(rhoval, inp.a)

            # convolve seperatrix function and bump function
            kappa_tilda_pre = np.convolve(k_sep, scaled_eta, 'same') / (( xnum - 1) / 2)  # Still don't understand why we need to divide by this, but we definitely need to.

            # post processing
            kappa_tilda_post = posttransform(np.linspace(pi, 2 * pi, xnum), kappa_tilda_pre)

            # create a 1D interpolation function for everywhere except for the seperatrix
            kappa_tilda_f = interp1d(np.linspace(pi, 2 * pi, xnum), kappa_tilda_post, kind="linear")
            for j in range(0, kappa_tilda.shape[1]):
                if theta[i, j] > pi and theta[i, j] <= 2 * pi:
                    kappa_tilda[i, j] = kappa_tilda_f(theta[i, j])
        else:
            kappa_tilda_f = interp1d(np.linspace(pi, 2 * pi, xnum), k_sep, kind="linear")
            for j in range(0, kappa_tilda.shape[1]):
                if theta[i, j] > pi and theta[i, j] <= 2 * pi:
                    kappa_tilda[i, j] = kappa_tilda_f(theta[i, j])

    kappa = kappa + kappa_tilda
    return kappa


def calc_sigv_fus(T, mode='dd'):
    def sigv(Ti, mode):  # function takes T in kev
        if mode == 'dt':
            B_G = 34.3827
            m_rc2 = 1124656

            C1 = 1.17302E-9
            C2 = 1.51361E-2
            C3 = 7.51886E-2
            C4 = 4.60643E-3
            C5 = 1.35000E-2
            C6 = -1.06750E-4
            C7 = 1.36600E-5

            theta = Ti / (1.0 - (Ti * (C2 + Ti * (C4 + Ti * C6))) / (1.0 + Ti * (C3 + Ti * (C5 + Ti * C7))))
            xi = (B_G ** 2.0 / (4.0 * theta)) ** (1.0 / 3.0)
            sigv = C1 * theta * np.sqrt(xi / (m_rc2 * Ti ** 3.0)) * np.exp(-3.0 * xi)
            sigv = sigv / 1.0E6  # convert from cm^3/s to m^3/s

        elif mode == 'dd':

            B_G = 31.3970
            m_rc2 = 937814

            # first for the D(d, p)T reaction
            C1_1 = 5.65718E-12
            C2_1 = 3.41267E-3
            C3_1 = 1.99167E-3
            C4_1 = 0.0
            C5_1 = 1.05060E-5
            C6_1 = 0.0
            C7_1 = 0.0

            theta_1 = Ti / (
                        1.0 - (Ti * (C2_1 + Ti * (C4_1 + Ti * C6_1))) / (1.0 + Ti * (C3_1 + Ti * (C5_1 + Ti * C7_1))))
            xi_1 = (B_G ** 2.0 / (4.0 * theta_1)) ** (1.0 / 3.0)
            sigv_1 = C1_1 * theta_1 * np.sqrt(xi_1 / (m_rc2 * Ti ** 3.0)) * np.exp(-3.0 * xi_1)

            # then for the D(d, n)He3 reaction

            C1_2 = 5.43360E-12
            C2_2 = 5.85778E-3
            C3_2 = 7.68222E-3
            C4_2 = 0.0
            C5_2 = -2.96400E-6
            C6_2 = 0.0
            C7_2 = 0.0

            theta_2 = Ti / (
                        1.0 - (Ti * (C2_2 + Ti * (C4_2 + Ti * C6_2))) / (1.0 + Ti * (C3_2 + Ti * (C5_2 + Ti * C7_2))))
            xi_2 = (B_G ** 2.0 / (4.0 * theta_2)) ** (1.0 / 3.0)
            sigv_2 = C1_2 * theta_2 * np.sqrt(xi_2 / (m_rc2 * Ti ** 3.0)) * np.exp(-3.0 * xi_2)

            sigv = (0.5 * sigv_1 + 0.5 * sigv_2) / 1.0E6  # convert from cm^3/s to m^3/s
        return sigv

    # create logspace over the relevant temperature range
    # (bosch hale technically only valid over 0.2 - 100 kev)
    Ti_range = np.logspace(-1, 2, 1000)  # values in kev
    sigv_fus_range = sigv(Ti_range, mode=mode)  # in m^3/s
    sigv_fus_interp = UnivariateSpline(Ti_range * 1.0E3 * 1.6021E-19, sigv_fus_range, s=0)  # converted to Joules

    sv_fus = sigv_fus_interp(T.i.J)
    dsv_fus_dT = sigv_fus_interp.derivative()(T.i.J)

    return sv_fus, dsv_fus_dT


def calc_sigv_ion(T):
    # TODO: configure so it can use any of the cross section libraries
    # currently using the Stacey-Thomas cross sections
    T_exps_fit = np.array([-1, 0, 1, 2, 3, 4, 5])
    sigv_exps_fit = np.array([-2.8523E+01, -1.7745E+01, -1.3620E+01,
                              -1.3097E+01, -1.3301E+01, -1.3301E+01, -1.3301E+01])
    interp1 = UnivariateSpline(T_exps_fit, sigv_exps_fit, s=0)

    T_exps_range = np.linspace(-1, 5, 1000)
    sigv_vals_range = 10.0 ** interp1(T_exps_range)  # in m^3/s

    T_vals_range = np.logspace(-1, 5, 1000) * 1.6021E-19  # in joules
    interp2 = UnivariateSpline(T_vals_range, sigv_vals_range, s=0)

    sv_ion = interp2(T.i.J)
    dsv_ion_dT = interp2.derivative()(T.i.J)

    return sv_ion, dsv_ion_dT


def calc_svel(T):
    tint = np.array([-1, 0, 1, 2, 3])
    tnnt = np.array([0, 1, 2])

    elast = np.array([[-1.3569E+01, -1.3337E+01, -1.3036E+01, -1.3569E+01, -1.3337E+01],
                      [-1.3036E+01, -1.3337E+01, -1.3167E+01, -1.3046E+01, -1.3036E+01],
                      [-1.3046E+01, -1.2796E+01, -1.3036E+01, -1.3046E+01, -1.2796E+01]])

    interp1 = interp2d(tint, tnnt, elast)

    Ti_exps = np.linspace(-1, 3, 100)
    Tn_exps = np.linspace(0, 2, 100)
    svel_vals = 10.0 ** (interp1(Ti_exps, Tn_exps))  # in m^3/s

    Ti_vals = np.logspace(-1, 3, 100) * 1.6021E-19  # in joules
    Tn_vals = np.logspace(0, 2, 100) * 1.6021E-19  # in joules

    dsvel_dTi_vals = np.gradient(svel_vals, Ti_vals, axis=0)

    Ti_vals2d, Tn_vals2d = np.meshgrid(Ti_vals, Tn_vals)

    Ti_mod = np.where(T.i.ev > 1E3, 1E3 * 1.6021E-19, T.i.ev * 1.6021E-19)
    Tn_mod = np.zeros(Ti_mod.shape) + 2.0 * 1.6021E-19

    sv_el = griddata(np.column_stack((Ti_vals2d.flatten(), Tn_vals2d.flatten())),
                     svel_vals.flatten(),
                     (Ti_mod, Tn_mod),
                     method='linear', rescale=False)

    dsv_el_dT = griddata(np.column_stack((Ti_vals2d.flatten(), Tn_vals2d.flatten())),
                         dsvel_dTi_vals.flatten(),
                         (Ti_mod, Tn_mod),
                         method='linear', rescale=False)

    return sv_el, dsv_el_dT


def calc_svcx_st(T):
    tint = np.array([-1, 0, 1, 2, 3])
    tnnt = np.array([0, 1, 2])

    cx = np.array([[-1.4097E+01, -1.3921E+01, -1.3553E+01, -1.4097E+01, -1.3921E+01],
                   [-1.3553E+01, -1.3921E+01, -1.3824E+01, -1.3538E+01, -1.3553E+01],
                   [-1.3538E+01, -1.3432E+01, -1.3553E+01, -1.3538E+01, -1.3432E+01]])

    interp1 = interp2d(tint, tnnt, cx)

    Ti_exps = np.linspace(-1, 3, 100)
    Tn_exps = np.linspace(0, 2, 100)
    svcx_vals = 10.0 ** (interp1(Ti_exps, Tn_exps))  # in m^3/s

    Ti_vals = np.logspace(-1, 3, 100) * 1.6021E-19  # in joules
    Tn_vals = np.logspace(0, 2, 100) * 1.6021E-19  # in joules

    dsvcx_dTi_vals = np.gradient(svcx_vals, Ti_vals, axis=0)

    Ti_vals2d, Tn_vals2d = np.meshgrid(Ti_vals, Tn_vals)

    Ti_mod = np.where(T.i.ev > 1E3, 1E3 * 1.6021E-19, T.i.ev * 1.6021E-19)
    Tn_mod = np.zeros(Ti_mod.shape) + 2.0 * 1.6021E-19

    sv_cx = griddata(np.column_stack((Ti_vals2d.flatten(), Tn_vals2d.flatten())),
                     svcx_vals.flatten(),
                     (Ti_mod, Tn_mod),
                     method='linear', rescale=False)

    dsv_cx_dT = griddata(np.column_stack((Ti_vals2d.flatten(), Tn_vals2d.flatten())),
                         dsvcx_dTi_vals.flatten(),
                         (Ti_mod, Tn_mod),
                         method='linear', rescale=False)

    return sv_cx, dsv_cx_dT


def calc_svrec_st(n, T):
    # TODO: check this calculation. -MH
    znint = np.array([16, 18, 20, 21, 22])
    Tint = np.array([-1, 0, 1, 2, 3])

    rec = np.array([[-1.7523E+01, -1.6745E+01, -1.5155E+01, -1.4222E+01, -1.3301E+01],
                    [-1.8409E+01, -1.8398E+01, -1.8398E+01, -1.7886E+01, -1.7000E+01],
                    [-1.9398E+01, -1.9398E+01, -1.9398E+01, -1.9398E+01, -1.9398E+01],
                    [-2.0155E+01, -2.0155E+01, -2.0155E+01, -2.0155E+01, -2.0155E+01],
                    [-2.1000E+01, -2.1000E+01, -2.1000E+01, -2.1000E+01, -2.1000E+01]])

    interp1 = interp2d(znint, Tint, rec)

    zni_exps = np.linspace(16, 22, 100)
    Ti_exps = np.linspace(-1, 3, 100)
    svrec_vals = 10.0 ** (interp1(zni_exps, Ti_exps))  # in m^3/s

    zni_vals = np.logspace(16, 22, 100)
    Ti_vals = np.logspace(-1, 3, 100) * 1.6021E-19  # in joules

    dsvrec_dTi_vals = np.gradient(svrec_vals, Ti_vals, axis=0)

    zni_vals2d, Ti_vals2d = np.meshgrid(zni_vals, Ti_vals)

    zni_mod = np.where(n.i > 1E22, 1E22, n.i)
    zni_mod = np.where(n.i < 1E16, 1E16, zni_mod)
    Ti_mod = np.where(T.i.ev > 1E3, 1E3 * 1.6021E-19, T.i.ev * 1.6021E-19)
    Ti_mod = np.where(T.i.ev < 1E-1, 1E-1 * 1.6021E-19, Ti_mod)

    sv_rec = griddata(np.column_stack((zni_vals2d.flatten(), Ti_vals2d.flatten())),
                      svrec_vals.flatten(),
                      (zni_mod, Ti_mod),
                      method='linear',
                      rescale=False)

    dsv_rec_dT = griddata(np.column_stack((zni_vals2d.flatten(), Ti_vals2d.flatten())),
                          dsvrec_dTi_vals.flatten(),
                          (zni_mod, Ti_mod),
                          method='linear',
                          rescale=False)

    return sv_rec, dsv_rec_dT


class MilCoreBrnd():
    """Calculates various plasma properties using a modified Miller geometry
    
    Methods:
        createbackround
        xmiller
    
    Attributes:
        r
        theta
        rho
        ni
        ne
        Ti_kev
        Ti_K
        Ti_J
        Te_kev
        Te_K
        Te_J
        nC
        E_pot
        pressure
        j_r
        kappa
        tri
        R
        Z
        diff_vol
        IP
        B_phi
        psi
        psi_norm
        B_p
        B_tot
        f_phi
    """
    def __init__(self, inp):
        sys.dont_write_bytecode = True

        # specify rho values
        try:
            rho1d = np.concatenate((np.linspace(0, inp.edge_rho, inp.rhopts_core, endpoint=False),
                                    np.linspace(inp.edge_rho, 1, inp.rhopts_edge)), axis=0)
        except AttributeError:
            try:
                rho1d = np.linspace(0, 1, inp.rhopts)
            except AttributeError:
                raise AttributeError("You haven't specified the number of radial points.")
        self.rhopts = len(rho1d)

        # calculate thetapts and specify theta values
        self.thetapts = int(4 * ceil(float(inp.thetapts_approx)/4))+1
        theta1d = np.linspace(0, 2*pi, self.thetapts)

        #create rho, theta, and r matrices
        # TODO: This isn't technically correct because rho is not a constant function of normalized psi
        self.theta, self.rho = np.meshgrid(theta1d, rho1d)
        self.r = self.rho * inp.a
        self.a = inp.a
        self.R0_a = inp.R0_a
        # populate radial density and temperature profiles
        try:
            ni = UnivariateSpline(inp.ni_data[:, 0], inp.ni_data[:, 1], k=3, s=2.0)(self.rho)
        except AttributeError:
            ni = calc_radial_prof(self.rho, inp.ni0, inp.ni9, inp.ni_sep, inp.nu_ni, ped_loc=0.9)

        try:
            ne = UnivariateSpline(inp.ne_data[:, 0], inp.ne_data[:, 1], k=3, s=2.0)(self.rho)
        except AttributeError:
            ne = calc_radial_prof(self.rho, inp.ne0, inp.ne9, inp.ne_sep, inp.nu_ne, ped_loc=0.9)

        try:
            fracz = UnivariateSpline(inp.fracz_data[:, 0], inp.fracz_data[:, 1], k=5, s=2.0)(self.rho)
        except AttributeError:
            fracz = np.full(self.rho.shape, 0.025)
        nC = ne * fracz

        try:
            Ti_kev = UnivariateSpline(inp.Ti_data[:, 0], inp.Ti_data[:, 1], k=5, s=2.0)(self.rho)
        except AttributeError:
            Ti_kev = calc_radial_prof(self.rho, inp.Ti0, inp.Ti9, inp.Ti_sep, inp.nu_Ti, ped_loc=0.9)

        try:
            Te_kev = UnivariateSpline(inp.Te_data[:, 0], inp.Te_data[:, 1], k=5, s=2.0)(self.rho)
        except AttributeError:
            Te_kev = calc_radial_prof(self.rho, inp.Te0, inp.Te9, inp.Te_sep, inp.nu_Te, ped_loc=0.9)

        Tns_kev = np.full(self.rho.shape, 0.002)
        Tnt_kev = Ti_kev
        TC_kev = Ti_kev
        Ti_K = Ti_kev * 1.159E7

        # populate main density namedtuple
        self.n = namedtuple('n', 'i e n C')(
            ni,
            ne,
            namedtuple('nn', 's t tot')(
                np.zeros(self.rho.shape),  # slow
                np.zeros(self.rho.shape),  # thermal
                np.zeros(self.rho.shape)  # total
            ),
            ne * fracz)

        # populate main temperature namedtuple
        self.T = namedtuple('T', 'i e n C')(
            namedtuple('T', 'kev ev J')(Ti_kev, Ti_kev * 1E3, Ti_kev * 1E3 * 1.6021E-19),
            namedtuple('T', 'kev ev J')(Te_kev, Te_kev * 1E3, Te_kev * 1E3 * 1.6021E-19),
            namedtuple('Tn', 's t')(
                namedtuple('T', 'kev ev J')(Tns_kev, Tns_kev * 1E3, Tns_kev * 1E3 * 1.6021E-19),
                namedtuple('T', 'kev ev J')(Tnt_kev, Tnt_kev * 1E3, Tnt_kev * 1E3 * 1.6021E-19)
            ),
            namedtuple('T', 'kev ev J')(TC_kev, TC_kev * 1E3, TC_kev * 1E3 * 1.6021E-19))

        # calculate pressure
        self.pressure = ni * k * Ti_K

        # calculate some impurity related quantities
        self.z_0 = nC*6.0**2 / ni
        self.z_eff = (ni*1.0**2 + nC*6.0**2) / ne

        # populate E_r and calculate the radial electrostatic potential
        try:
            E_r_fit = UnivariateSpline(inp.er_data[:, 0], inp.er_data[:, 1], k=5, s=2.0)
            self.E_r = E_r_fit(self.rho)
            self.E_pot = np.zeros(self.rho.shape)
            for i, rhoval in enumerate(rho1d):
                self.E_pot[i] = E_r_fit.integral(rhoval, 1.0)
        except AttributeError:
            raise AttributeError("You need E_r data")

        # populate the radial current density profile
        try:
            self.j_r = calc_radial_prof(self.rho, inp.j0, 0, 0, inp.nu_j, ped_loc=1.0)
        except AttributeError:
            raise AttributeError("You haven't specified a current distribution.")

        # populate rotation velocities
        try:
            self.vpolC = UnivariateSpline(inp.vpolC_data[:, 0], inp.vpolC_data[:, 1], k=5, s=2.0)(self.rho)
        except AttributeError:
            self.vpolC = np.zeros(self.rho.shape)

        try:
            self.vtorC = UnivariateSpline(inp.vtorC_data[:, 0], inp.vtorC_data[:, 1], k=5, s=2.0)(self.rho)
        except AttributeError:
            self.vtorC = np.zeros(self.rho.shape)

        try:
            self.vpolD = UnivariateSpline(inp.vpolD_data[:, 0], inp.vpolD_data[:, 1], k=5, s=2.0)(self.rho)
        except AttributeError:
            self.vpolD = np.zeros(self.rho.shape)

        try:
            self.vtorD = UnivariateSpline(inp.vtorD_data[:, 0], inp.vtorD_data[:, 1], k=5, s=2.0)(self.rho)
        except AttributeError:
            self.vtorD = np.zeros(self.rho.shape)

        # populate experimental q-profile, if present
        # TODO: decide whether to even keep this
        try:
            self.q = UnivariateSpline(inp.q_data[:, 0], inp.q_data[:, 1], k=5, s=2.0)(self.rho)
        except AttributeError:
            self.q = np.zeros(self.rho.shape)  # will calculated later with the other miller stuff

        # create kappa and triangularity matrices
        upperhalf = (self.theta >= 0) & (self.theta < pi)
        kappa = np.where(upperhalf,
                         inp.kappa_up / (inp.a**inp.s_k_up) * self.r**inp.s_k_up,
                         inp.kappa_lo / (inp.a**inp.s_k_lo) * self.r**inp.s_k_lo)

        if inp.xmil==1:
            xpt = np.array([inp.xpt_R, inp.xpt_Z])
            kappa = xmiller(kappa, inp)
            tri_lo = sin(3*pi/2 - acos((xpt[0]-inp.R0_a)/inp.a))
            tri_up = inp.tri_up
        else:
            tri_lo = inp.tri_lo
            tri_up = inp.tri_up

        self.kappa_vals = namedtuple('kappa', 'axis sep')(0, inp.kappa_up)
        self.tri_vals = namedtuple('tri', 'axis sep')(0, inp.tri_up)

        tri = np.where(upperhalf,
                         tri_up * self.rho,
                         tri_lo * self.rho)

        # CALCULATE INITIAL R, Z WITH NO SHAFRANOV SHIFT
        # (NECESSARY TO GET ESTIMATES OF L_r WHEN CALCULATING SHAFRANOV SHIFT)
        R0 = np.ones(self.rho.shape) * inp.R0_a
        R = R0 + self.r * np.cos(self.theta + np.arcsin(tri * np.sin(self.theta)))
        Z = kappa * self.r * np.sin(self.theta)

        # THIS CALCULATES A MATRIX OF THE LENGTHS OF EACH SECTION OF EACH FLUX
        # SURFACE AND THEN SUMS THEM TO GET THE PERIMETER IN 2D OF EACH FLUX
        # SURFACE (VALUE OF r).
        L_seg = np.sqrt((Z - np.roll(Z, -1, axis=1)) ** 2 + (R - np.roll(R, -1, axis=1)) ** 2)
        L_seg[:, -1] = 0
        L_r = np.tile(np.sum(L_seg, axis=1), (self.thetapts, 1)).T

        # CALCULATE CROSS-SECTIONAL AREA CORRESPONDING TO EACH r AND ASSOCIATED
        # DIFFERENTIAL AREAS
        area = np.zeros(self.rho.shape)
        for i in range(0, len(self.rho)):
            area[i, :] = PolyArea(R[i, :], Z[i, :])

        diff_area = area - np.roll(area, 1, axis=0)
        diff_area[0, :] = 0

        diff_vol = diff_area * 2 * pi * inp.R0_a  # approx because it uses R0_a instead of shifted R0
        vol = np.cumsum(diff_vol, axis=0)

        # Calculate each differential I and sum to get cumulative I
        j_r_ave = np.roll((self.j_r + np.roll(self.j_r, -1, axis=0)) / 2.0, 1, axis=0)
        j_r_ave[0, :] = 0
        diff_I = diff_area * j_r_ave
        self.I = np.cumsum(diff_I, axis=0)
        self.IP = self.I[-1, 0]

        # Calculate B_p_bar
        B_p_bar = mu_0 * self.I / L_r
        B_p_bar[0, :] = 0

        # Calculate li
        li = (np.cumsum(B_p_bar ** 2 * diff_vol, axis=0) / vol) / (2 * B_p_bar ** 2)
        li[0, :] = 0

        # Calculate beta_p
        beta_p = 2 * mu_0 * (np.cumsum(self.pressure * diff_vol, axis=0) / vol - self.pressure) / B_p_bar ** 2

        # Calculate dR0dr
        dR0dr = np.zeros(self.rho.shape)
        R0 = np.zeros(self.rho.shape)

        f = 2 * (kappa ** 2 + 1) / (3 * kappa ** 2 + 1) * (beta_p + li / 2) + 1 / 2 * (
                    kappa ** 2 - 1) / (3 * kappa ** 2 + 1)
        f[0, :] = f[1, :]  # TODO: NEED TO REVISIT, SHOULD EXTRAPOLATE SOMEHOW

        dR0dr[-1, :] = -2.0 * inp.a * f[-1, :] / inp.R0_a
        R0[-1, :] = inp.R0_a

        for i in range(len(self.rho) - 2, -1, -1):
            R0[i, :] = dR0dr[i + 1, :] * (self.r[i, :] - self.r[i + 1, :]) + R0[i + 1, :]
            dR0dr[i, :] = -2.0 * self.r[i, :] * f[i, :] / R0[i, :]

        # NOW USE UPDATED R0 AND dR0dr to get new R, Z.
        R = R0 + self.r * np.cos(self.theta + np.arcsin(tri * np.sin(self.theta)))
        Z = kappa * self.r * np.sin(self.theta) + inp.Z0

        # RECALCULATE L_seg and L_r
        L_seg = np.sqrt((Z - np.roll(Z, -1, axis=1)) ** 2 + (R - np.roll(R, -1, axis=1)) ** 2)
        L_seg[:, -1] = 0
        L_r = np.tile(np.sum(L_seg, axis=1), (self.thetapts, 1)).T

        # Update areas and volumes
        area = np.zeros(self.rho.shape)
        for i in range(0, len(self.rho)):
            area[i, :] = PolyArea(R[i, :], Z[i, :])

        self.diff_area = area - np.roll(area, 1, axis=0)
        self.diff_area[0, :] = 0

        # TODO: fix this to use the correct R0 for each value of rho
        self.diff_vol = diff_area * 2 * pi * inp.R0_a  # approx because it uses R0_a instead of shifted R0
        self.vol = np.cumsum(diff_vol, axis=0)

        # RECALCULATE GRAD-r
        dkappa_dtheta = np.gradient(kappa, edge_order=1)[1] * self.thetapts / (2 * pi)
        dkappa_dr = np.gradient(kappa, edge_order=1)[0] * self.rhopts / inp.a

        dkappa_dtheta[-1] = dkappa_dtheta[-2]
        dkappa_dr[-1] = dkappa_dr[-2]

        dZ_dtheta = np.gradient(Z, edge_order=2)[1] * self.thetapts / (2 * pi)  # self.r*(self.kappa*np.cos(self.theta)+dkappa_dtheta*np.sin(self.theta))
        dZ_dr = np.gradient(Z, edge_order=2)[0] * self.rhopts / inp.a  # np.sin(self.theta)*(self.r*dkappa_dr + self.kappa)
        dR_dr = np.gradient(R, edge_order=2)[0] * self.rhopts / inp.a  # dR0dr - np.sin(self.theta + np.sin(self.theta)*np.arcsin(tri))*(np.sin(self.theta)*s_tri) + np.cos(self.theta+np.sin(self.theta)*np.arcsin(tri))
        dR_dtheta = np.gradient(R, edge_order=2)[1] * self.thetapts / ( 2 * pi)  # -self.r*np.sin(self.theta+np.sin(self.theta)*np.arcsin(tri))*(1+np.cos(self.theta)*np.arcsin(tri))

        abs_grad_r = np.sqrt(dZ_dtheta ** 2 + dR_dtheta ** 2) / np.abs(dR_dr * dZ_dtheta - dR_dtheta * dZ_dr)

        # we now have the final R and Z matrices, and can calculate quantities that depend on R
        self.R = R
        self.Z = Z
        self.BT = inp.BT0 * self.R[0, 0] / self.R
        self.xpt_loc = np.where(self.Z == np.amin(self.Z))

        self.shaf_shift = (self.R[0, 0] - inp.R0_a)/self.a
        # We want to calculate the poloidal field strength everywhere. The problem is that we've got two equations
        # in three unknowns. However, if we assume that the poloidal integral of the flux surface average of the
        # poloidal magnetic field is approximately the same as the poloidal integral of the actual poloidal magnetic
        # field, then we can calculate the q profile.

        # Calculate initial crappy guess on q
        # Equation 16 in the miller paper. The last term is how I'm doing a flux surface average
        q_mil = inp.BT0 * self.R[0, 0] / (2 * pi * B_p_bar) * np.tile(np.sum(L_seg / self.R ** 2, axis=1),
                                                                      (self.thetapts, 1)).T
        q_mil[0, :] = q_mil[1, :]

        dpsidr = (inp.BT0 * self.R[0, 0]) / \
                 (2 * pi * q_mil) * np.tile(np.sum(L_seg / (self.R * abs_grad_r), axis=1),
                                            (self.thetapts, 1)).T

        self.psi = np.zeros(self.r.shape)
        for index, row in enumerate(self.r):
            if index >= 1:
                self.psi[index] = dpsidr[index] * (self.r[index, 0] - self.r[index - 1, 0]) + self.psi[index - 1]
        self.psi_norm = self.psi / self.psi[-1, 0]

        self.B_p = dpsidr * 1 / self.R * abs_grad_r
        self.B_p[0, :] = 0

        self.B_t = inp.BT0 * self.R[0, 0] / self.R
        self.B_tot = np.sqrt(self.B_p ** 2 + self.B_t ** 2)
        self.f_phi = self.B_t / self.B_tot


        # calculate cross sections
        self.sv_fus, self.dsv_fus_dT = calc_sigv_fus(self.T, mode='dt')
        self.sv_ion, self.dsv_ion_dT = calc_sigv_ion(self.T)
        self.sv_el, self.dsv_el_dT = calc_svel(self.T)
        self.sv_cx, self.dsv_cx_dT = calc_svcx_st(self.T)


    def miller(self, p):

        # All we're doing with kappa in this next part is making the derivative between upper and lower
        # elongation continuous by "smoothing out" the "step function"
        # using f(x) = tanh(B*sin(x)), where be controlls how smooth or squre the function is.
        # Plot that function and you'll see what we're doing. This is necessary
        # to prevent shafranov shift from producing ugly pictures with high poloidal
        # resolution. It also makes Richard's stuff easier. Just deal with it
        # and don't put this in any papers. It's just a bandaid. We do the same
        # thing with triangularity. - MH
        
        # B_kappa = 0.0
        # self.kappa  = (((p.kappa_up / (p.a**p.s_k_up) * self.r**p.s_k_up) - (p.kappa_lo / (p.a**p.s_k_lo) * self.r**p.s_k_lo))/2.0
        #        * np.tanh(B_kappa*np.sin(self.theta))
        #        + ((p.kappa_up / (p.a**p.s_k_up) * self.r**p.s_k_up) + (p.kappa_lo / (p.a**p.s_k_lo) * self.r**p.s_k_lo))/2.0)

        # s_tri = np.where(upperhalf,
        #                 self.r*p.tri_up/(p.a*np.sqrt(1-tri)),
        #                 self.r*tri_lo/(p.a*np.sqrt(1-tri)))
        pass