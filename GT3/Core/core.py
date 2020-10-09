#!/usr/bin/env python2
# !/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 18 13:22:31 2018

@author: max
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata, UnivariateSpline
from math import pi
from collections import namedtuple
from scipy import constants
from shapely.geometry import LineString, Point, MultiPoint
from matplotlib.path import Path
from Functions.FindXPtMagAxis import find_xpt_mag_axis
from Functions.CalcPsiNorm import calc_psi_norm
from Functions.CalcRho2PsiInterp import calc_rho2psi_interp
from Functions.CalcGrad import calc_grad
from Functions.CalcFSA import calc_fsa
from Functions.CalcSV import calc_svion_st, calc_svel_st, calc_svrec_st, calc_svcx_st, calc_svfus
from Functions.CalcFsPerimInt import calc_fs_perim_int
from Functions.CreateSurfAreaInterp import create_surf_area_interp
from Functions.CalcTheta1D import calc_theta1d
from Functions.CalcPtsLines import calc_pts_lines
from Functions.CalcRZ import calc_RZ
from Functions.CalcKappaElong import calc_kappa_elong
from Functions.CreateVolInterp import create_vol_interp
from Functions.CalcChiJet import calc_chi_jet
from scipy.interpolate import interp1d
import GT3

e = constants.elementary_charge
u_0 = constants.mu_0
m_d = 3.343583719e-27
m_t = 5.006e-27
m_c = 1.9926467e-26
m_a = 6.643e-27


# def calc_imp_rad(nz, ne, rho, Lz, calc_dVdrho, rho_start=0, rho_end=1):
#     nz_fit = UnivariateSpline(rho[:,0], nz[:,0], k=3, s=0)
#     ne_fit = UnivariateSpline(rho[:,0], ne[:,0], k=3, s=0)
#
#     return


# noinspection SpellCheckingInspection
class Core:

    sep_val_overridden = False  # type: bool
    """Has the separatrix value been overridden?"""
    sep_val = 1.0  # type: float
    """The value of the separatrix. If this is set in the input file as <1.0, the X-point calculations are omitted"""
    r2vol = None  # type: UnivariateSpline
    """The V(r) interpolator"""
    rho2vol = None  # type: UnivariateSpline
    """The V(rho) interpolator"""
    psinorm2vol = None  # type: UnivariateSpline
    """The V(psi) interpolator on [0., 1.]"""

    def __init__(self, inp):
        # type: (GT3.ReadInfile) -> GT3.Core
        # Hold the input as a class attribute for use in functions

        self.inp = inp
        self.wall_line = inp.wall_line
        self.sep_val = inp.sep_val
        if abs(1.0 - float(self.sep_val)) > 0.0001:
            self.sep_val_overridden = True
            print "The separatrix value is being overwritten. GT3.SOL will be omitted."
        # this assumes that psi is given as a square array
        psi_shape = int(np.sqrt(inp.psirz_exp[:, 0].shape[0]))  # type: int

        raw_psi_R = inp.psirz_exp[:, 0].reshape(-1, psi_shape)  # type: np.ndarray
        raw_psi_Z = inp.psirz_exp[:, 1].reshape(-1, psi_shape)  # type: np.ndarray
        try:
            raw_psi = inp.psirz_exp[:, 2].reshape(-1, psi_shape) * inp.psi_scale  # type: np.ndarray
        except AttributeError:
            raw_psi = inp.psirz_exp[:, 2].reshape(-1, psi_shape)  # type: np.ndarray

        xpt_l, xpt_u, mag_axis = find_xpt_mag_axis(self, raw_psi_R, raw_psi_Z, raw_psi)
        xpt = [xpt_l, xpt_u]
        raw_psi_norm = calc_psi_norm(raw_psi_R, raw_psi_Z, raw_psi, xpt, mag_axis)

        raw_dpsidR = np.abs(np.gradient(raw_psi_norm, raw_psi_R[0, :], axis=1))
        raw_1_R_dpsidR = raw_dpsidR / raw_psi_R
        raw_d_dR_1_R_dpsidR = np.abs(np.gradient(raw_1_R_dpsidR, raw_psi_R[0, :], axis=1))
        raw_dpsidZ = np.abs(np.gradient(raw_psi_norm, raw_psi_Z[:, 0], axis=0))
        raw_d2psidZ2 = np.abs(np.gradient(raw_dpsidZ, raw_psi_Z[:, 0], axis=0))
        raw_dpsidr = raw_dpsidR + raw_dpsidZ
        raw_j = -(raw_d_dR_1_R_dpsidR * raw_psi_R + raw_d2psidZ2) / (raw_psi_R * u_0)

        PsiData = namedtuple('PsiData', 'R Z psi psi_norm dpsidR dpsidZ dpsidr j')
        """NamedTuple containing psi data, including R, Z, psi, psi_norm, dpsidR, dpsidZ, dpsidr, j"""

        self.psi_data = PsiData(
            raw_psi_R,
            raw_psi_Z,
            raw_psi,
            raw_psi_norm,
            raw_dpsidR,
            raw_dpsidZ,
            raw_dpsidr,
            raw_j
        )

        # calculate some important points and lines
        self.pts, self.lines = calc_pts_lines(self.psi_data, xpt, self.wall_line, mag_axis, self.sep_val, core=self)

        # create interpolation functions to convert rho to psi and vice versa
        interp_fns = calc_rho2psi_interp(self.pts, self.psi_data, self.sep_val)

        self.rho2psi = interp_fns[0]  # type: interp1d
        self.rho2psinorm = interp_fns[1]  # type: interp1d
        self.psi2rho = interp_fns[2]  # type: interp1d
        self.psi2psinorm = interp_fns[3]  # type: interp1d
        self.psinorm2rho = interp_fns[4]  # type: interp1d
        self.psinorm2psi = interp_fns[5]  # type: interp1d

        self.R0_a = self.pts.axis.mag[0]
        self.R0_g = self.pts.axis.geo[0]


        # specify rho values
        try:
            rho1d = np.concatenate((np.linspace(0, inp.edge_rho, inp.rhopts_core, endpoint=False),
                                    np.linspace(inp.edge_rho, self.sep_val, inp.rhopts_edge)), axis=0)
        except AttributeError:
            try:
                rho1d = np.linspace(0, self.sep_val, inp.rhopts)
            except AttributeError:
                raise AttributeError("You haven't specified the number of radial points.")
        self.rhopts = len(rho1d)

        theta1d, theta_markers = calc_theta1d(self.pts, inp.thetapts_approx)
        theta_xpt = theta_markers[3]
        self.thetapts = len(theta1d)

        # create rho and theta arrays (initializing the main computational grid)
        self.theta, self.rho = np.meshgrid(theta1d, rho1d)


        self.psi = self.rho2psi(self.rho)
        self.psi_norm = self.rho2psinorm(self.rho)

        # estimate elongation and triangularity
        self.kappa_vals, self.tri_vals = calc_kappa_elong(self.psi_data, np.asarray(self.lines.sep_closed.coords))

        # Calculate plasma radius considering potential for wall limiter
        # Inboard limiter
        if self.pts.ibmp[0] < self.inp.wall_line.bounds[0]:
            self.a = np.average((self.pts.obmp[0] - self.pts.axis.mag[0], self.pts.axis.mag[0] - self.inp.wall_line.bounds[0]))
        # Outboard limiter
        elif self.pts.obmp[0] >= self.inp.wall_line.bounds[2]:
            self.a = np.average((self.inp.wall_line.bounds[2] - self.pts.axis.mag[0], self.pts.axis.mag[0] - self.pts.ibmp[0]))
        else:
            self.a = np.average((self.pts.obmp[0] - self.pts.axis.mag[0], self.pts.axis.mag[0] - self.pts.ibmp[0]))

        self.shaf_shift = (self.pts.axis.mag[0] - self.pts.axis.geo[0]) / self.a
        self.r = self.rho * self.a

        self.R, self.Z = calc_RZ(self, self.rho, self.theta, theta_xpt, self.pts, self.psi_data, self.psi_norm, self.lines)

        self.xpt_loc = np.where(self.Z == np.amin(self.Z))
        self.obmp_loc = np.where(self.R == np.amax(self.R))
        self.ibmp_loc = np.where(self.R == np.amin(self.R))
        self.top_loc = np.where(self.Z == np.amax(self.Z))



        # create interpolation functions to obtain the flux surface surface area for any value of r, rho, or psi_norm
        self.r2sa, self.rho2sa, self.psinorm2sa = create_surf_area_interp(self.rho2psinorm,
                                                                          self.psinorm2rho,
                                                                          self.psi_data,
                                                                          np.asarray(self.lines.sep_closed.coords),
                                                                          self.R0_a,
                                                                          self.a,
                                                                          self.sep_val)

        # create interpolation functions to obtain the plasma volume for any value of r, rho, or psi_norm
        self.r2vol, self.rho2vol, self.psinorm2vol = create_vol_interp(self.rho2psinorm,
                                                                       self.psinorm2rho,
                                                                       self.psi_data,
                                                                       np.asarray(self.lines.sep_closed.coords),
                                                                       self.R0_a,
                                                                       self.a,
                                                                       self.sep_val)

        np.savetxt('psi2rho.txt', np.column_stack((np.linspace(0, self.sep_val, 100),
                                                   self.psi2rho(np.linspace(0, self.sep_val, 100)))))

        # initialize volume as a function of rho
        self.vol_rho = self.rho2vol(self.rho)
        """The volume as a function of rho V(rho)"""

        self.dVdrho = UnivariateSpline(np.linspace(0, self.sep_val, 100), self.rho2vol(np.linspace(0, self.sep_val, 100)), k=3,
                                       s=0).derivative()
        """The dV/drho(rho) interpolator interpolated on [0.,1.]"""

        self.dVdr = UnivariateSpline(np.linspace(0, self.sep_val, 100), self.dVdrho(np.linspace(0, self.sep_val, 100)) / self.a, k=3, s=0)
        """The dV/r(rho) interpolator interpolated on [0.,1.]"""

        self.vol = self.rho2vol(self.sep_val)
        """THe total volume of the plasma"""

        # initialize ionization rate arrays with zero
        self.izn_rate = namedtuple('izn_rate', 's t tot')(
            np.zeros(self.rho.shape),  # slow
            np.zeros(self.rho.shape),  # thermal
            np.zeros(self.rho.shape)  # total
        )

        # initialize cooling rate array with zero
        self.cool_rate = np.zeros(self.rho.shape)

        # initialize main density object

        # deuterium density
        try:
            nD = UnivariateSpline(inp.nD_data[:, 0], inp.nD_data[:, 1], k=3, s=0)(self.rho)
        except (AttributeError, TypeError):
            nD = np.zeros(self.rho.shape)

        # tritium density
        try:
            nT = UnivariateSpline(inp.nT_data[:, 0], inp.nT_data[:, 1], k=3, s=0)(self.rho)
        except (AttributeError, TypeError):
            nT = np.zeros(self.rho.shape)

        ni = nD + nT

        # electron density
        try:
            ne = UnivariateSpline(inp.ne_data[:, 0], inp.ne_data[:, 1], k=3, s=0)(self.rho)
        except (AttributeError, TypeError):
            ne = np.zeros(self.rho.shape)

        # carbon density
        try:
            nC = UnivariateSpline(inp.nC_data[:, 0], inp.nC_data[:, 1], k=3, s=0)(self.rho)
        except (AttributeError, TypeError):
            try:
                frac_C = UnivariateSpline(inp.frac_C_data[:, 0], inp.frac_C_data[:, 1], k=3, s=0)(self.rho)
                nC = ne * frac_C
            except (AttributeError, TypeError):
                nC = np.zeros(self.rho.shape)

        # tungsten density
        try:
            nW = UnivariateSpline(inp.nW_data[:, 0], inp.nW_data[:, 1], k=3, s=0)(self.rho)
        except (AttributeError, TypeError):
            try:
                frac_W = UnivariateSpline(inp.frac_W_data[:, 0], inp.frac_W_data[:, 1], k=3, s=0)(self.rho)
                nW = ne * frac_W
            except (AttributeError, TypeError):
                nW = np.zeros(self.rho.shape)

        # beryllium density
        try:
            nBe = UnivariateSpline(inp.nBe_data[:, 0], inp.nBe_data[:, 1], k=3, s=0)(self.rho)
        except (AttributeError, TypeError):
            try:
                frac_Be = UnivariateSpline(inp.frac_Be_data[:, 0], inp.frac_Be_data[:, 1], k=3, s=0)(self.rho)
                nBe = ne * frac_Be
            except (AttributeError, TypeError):
                nBe = np.zeros(self.rho.shape)

        # neon density
        try:
            nNe = UnivariateSpline(inp.nNe_data[:, 0], inp.nNe_data[:, 1], k=3, s=0)(self.rho)
        except (AttributeError, TypeError):
            try:
                frac_Ne = UnivariateSpline(inp.frac_Ne_data[:, 0], inp.frac_Ne_data[:, 1], k=3, s=0)(self.rho)
                nNe = ne * frac_Ne
            except (AttributeError, TypeError):
                nNe = np.zeros(self.rho.shape)

        # krypton density
        try:
            nKr = UnivariateSpline(inp.nKr_data[:, 0], inp.nKr_data[:, 1], k=3, s=0)(self.rho)
        except (AttributeError, TypeError):
            try:
                frac_Kr = UnivariateSpline(inp.frac_Kr_data[:, 0], inp.frac_Kr_data[:, 1], k=3, s=0)(self.rho)
                nKr = ne * frac_Kr
            except (AttributeError, TypeError):
                nKr = np.zeros(self.rho.shape)

        # argon density
        try:
            nAr = UnivariateSpline(inp.nAr_data[:, 0], inp.nAr_data[:, 1], k=3, s=0)(self.rho)
        except (AttributeError, TypeError):
            try:
                frac_Ar = UnivariateSpline(inp.frac_Ar_data[:, 0], inp.frac_Ar_data[:, 1], k=3, s=0)(self.rho)
                nAr = ne * frac_Ar
            except (AttributeError, TypeError):
                nAr = np.zeros(self.rho.shape)

        # alpha density
        try:
            na = UnivariateSpline(inp.na_data[:, 0], inp.na_data[:, 1], k=3, s=0)(self.rho)
        except (AttributeError, TypeError):
            try:
                frac_a = UnivariateSpline(inp.frac_a_data[:, 0], inp.frac_a_data[:, 1], k=3, s=0)(self.rho)
                na = ne * frac_a
            except (AttributeError, TypeError):
                na = np.zeros(self.rho.shape)

        nn = namedtuple('nn', 's t tot')(
            ne * 1E-7,  # slow
            ne * 1E-7,  # thermal
            2 * ne * 1E-7  # total
        )
        self.n = namedtuple('n', 'D T i e C W Be Ne Ar Kr a n')(nD, nT, ni, ne, nC, nW, nBe, nNe, nAr, nKr, na, nn)

        self.z_0 = nC * 6.0 ** 2 / (nD + nT)  # TODO: Update this calculation
        self.z_eff = ((nD + nT) * 1.0 ** 2 + nC * 6.0 ** 2) / ne  # TODO: Update this calculation

        # populate temperature namedtuples, including the main T object
        Ti_kev_fit = UnivariateSpline(
            inp.Ti_data[:, 0],
            inp.Ti_data[:, 1],
            k=3,
            s=0.0
        )

        Ti_kev = Ti_kev_fit(self.rho)
        Te_kev = UnivariateSpline(inp.Te_data[:, 0], inp.Te_data[:, 1], k=3, s=0.0)(self.rho)
        Tns_kev = np.full(self.rho.shape, 0.002)
        Tnt_kev = Ti_kev
        try:
            TC_kev = UnivariateSpline(inp.TC_data[:, 0], inp.TC_data[:, 1], k=3, s=0.0)(self.rho)
        except:
            TC_kev = Ti_kev

        self.T = namedtuple('T', 'i e n C')(
            namedtuple('Ti', 'kev ev J')(
                Ti_kev,
                Ti_kev * 1E3,
                Ti_kev * 1E3 * e),
            namedtuple('eT', 'kev ev J')(
                Te_kev,
                Te_kev * 1E3,
                Te_kev * 1E3 * e),
            namedtuple('Tn', 's t')(
                namedtuple('Tn_s', 'kev ev J')(
                    Tns_kev,
                    Tns_kev * 1E3,
                    Tns_kev * 1E3 * e),
                namedtuple('Tn_t', 'kev ev J')(
                    Tnt_kev,
                    Tnt_kev * 1E3,
                    Tnt_kev * 1E3 * e)
            ),
            namedtuple('TC', 'kev ev J')(
                TC_kev,
                TC_kev * 1E3,
                TC_kev * 1E3 * e))

        # initialize pressures
        self.p = namedtuple('p', 'i e C')(
            self.n.i * self.T.i.J,
            self.n.e * self.T.e.J,
            self.n.C * self.T.C.J)

        # initialize spatial gradients and gradient scale lengths
        self.dni_dr = calc_grad(self.rho, self.n.i, self.psi_norm, self.R, self.Z, self.psi_data)
        # plt.contourf(self.R, self.Z, dni_dr, 500)
        # plt.axis('equal')
        # plt.colorbar()
        # plt.show()
        # sys.exit()
        self.dne_dr = calc_grad(self.rho, self.n.e, self.psi_norm, self.R, self.Z, self.psi_data)
        self.dnC_dr = calc_grad(self.rho, self.n.C, self.psi_norm, self.R, self.Z, self.psi_data)
        self.dTi_kev_dr = calc_grad(self.rho, self.T.i.kev, self.psi_norm, self.R, self.Z, self.psi_data)
        self.dTi_ev_dr = self.dTi_kev_dr * 1E3
        self.dTi_J_dr = self.dTi_kev_dr * 1E3 * e
        self.dTe_kev_dr = calc_grad(self.rho, self.T.i.kev, self.psi_norm, self.R, self.Z, self.psi_data)
        self.dTe_ev_dr = self.dTe_kev_dr * 1E3
        self.dTe_J_dr = self.dTe_kev_dr * 1E3 * e

        self.dpi_dr = calc_grad(self.rho, self.n.i * self.T.i.J, self.psi_norm, self.R, self.Z, self.psi_data)
        # Nick: Added Carbon Pressure gradient
        self.dpC_dr = calc_grad(self.rho, self.n.C * self.T.C.J, self.psi_norm, self.R, self.Z, self.psi_data)
        self.dpe_dr = calc_grad(self.rho, self.n.e * self.T.e.J, self.psi_norm, self.R, self.Z, self.psi_data)

        self.L = namedtuple('L', 'n T p')(
            namedtuple('Ln', 'i e')(
                -self.n.i / self.dni_dr,
                -self.n.e / self.dne_dr
            ),
            namedtuple('LT', 'i e')(
                -self.T.i.kev / self.dTi_kev_dr,  # note: independent of units,
                -self.T.e.kev / self.dTe_kev_dr  # note: independent of units
            ),
            namedtuple('Lp', 'i e')(
                -self.p.i / self.dpi_dr,  # note: independent of units,
                -self.p.e / self.dpe_dr  # note: independent of units
            )
        )

        try:
            try:
                E_r_fit = UnivariateSpline(inp.er_data[:, 0], inp.er_data[:, 1] * inp.Er_scale, k=3, s=0)
            except (AttributeError, TypeError):
                E_r_fit = UnivariateSpline(inp.er_data[:, 0], inp.er_data[:, 1], k=3, s=0)
            self.E_r = E_r_fit(self.rho)
            self.E_pot = np.zeros(self.rho.shape)
            try:
                for i, rhoval in enumerate(rho1d):
                    self.E_pot[i] = E_r_fit.integral(rhoval, self.sep_val)
            except:
                print "Error in E_pot integration"
        except:
            print 'Er data not supplied. Setting E_r and E_pot to zero.'
            self.E_r = np.zeros(self.rho.shape)
            self.E_pot = np.zeros(self.rho.shape)
        # initialize E_r and the corresponding electric potential

        # initialize rotation velocities from data
        try:
            vpolD = UnivariateSpline(inp.vpolD_data[:, 0], inp.vpolD_data[:, 1], k=3, s=0)(self.rho)
        except (AttributeError, TypeError):
            vpolD = np.zeros(self.rho.shape)
        try:
            vpolC = UnivariateSpline(inp.vpolC_data[:, 0], inp.vpolC_data[:, 1], k=3, s=0)(self.rho)
        except (AttributeError, TypeError):
            vpolC = np.zeros(self.rho.shape)
        try:
            vtorD = UnivariateSpline(inp.vtorD_data[:, 0], inp.vtorD_data[:, 1], k=3, s=0)(self.rho)
        except (AttributeError, TypeError):
            vtorD = np.zeros(self.rho.shape)
        try:
            vtorC = UnivariateSpline(inp.vtorC_data[:, 0], inp.vtorC_data[:, 1], k=3, s=0)(self.rho)
        except (AttributeError, TypeError):
            vtorC = np.zeros(self.rho.shape)

        self.v = namedtuple('v_comp', 'pol tor')(
            namedtuple('v_spec', 'D C')(vpolD, vpolC),
            namedtuple('v_spec', 'D C')(vtorD, vtorC),
        )

        self.v_1D = namedtuple('v_comp', 'pol tor')(
            namedtuple('v_spec', 'D C')(vpolD[:, 0], vpolC[:, 0]),
            namedtuple('v_spec', 'D C')(vtorD[:, 0], vtorC[:, 0]),
        )

        # initialize magnetic field-related quantities
        B_pol_raw = np.sqrt((np.gradient(self.psi_data.psi, self.psi_data.R[0], axis=1) / self.psi_data.R) ** 2 +
                            (np.gradient(self.psi_data.psi, self.psi_data.Z.T[0], axis=0) / self.psi_data.R) ** 2)

        self.B_p = griddata(np.column_stack((raw_psi_R.flatten(), raw_psi_Z.flatten())),
                            B_pol_raw.flatten(),
                            (self.R, self.Z),
                            method='linear')

        self.B_t = inp.BT0 * self.pts.axis.mag[0] / self.R
        self.B_tot = np.sqrt(self.B_p ** 2 + self.B_t ** 2)
        self.f_phi = self.B_t / self.B_tot

        try:
            # use input q profile if given
            q_interp = UnivariateSpline(inp.q_data[:, 0], inp.q_data[:, 1], k=3, s=1.0)
            self.q = q_interp(self.rho)
            self.q_1D = q_interp(self.rho.T[0])
        except:
            # otherwise calculate q-profile from psi data
            self.q_1D = inp.BT0 * self.pts.axis.mag[0] / (2 * pi) * calc_fs_perim_int(1.0 / (self.R ** 2 * self.B_p),
                                                                                      self.R, self.Z)
            self.q_1D[0] = self.q_1D[1]
            self.q = np.repeat(self.q_1D[np.newaxis, :], self.rho.shape[1], axis=0).T
        self.q0 = self.q_1D[0]
        self.q95 = UnivariateSpline(self.psi_norm[:, 0], self.q_1D, k=1, s=0)(0.95)

        # create Lz-related variables. These will remain zero unless set by the ImpRad module
        self.Lz_C = namedtuple('Lz_C', 's t ddT')(
            np.zeros(self.rho.shape),
            np.zeros(self.rho.shape),
            namedtuple('Lz_ddT', 's t')(
                np.zeros(self.rho.shape),
                np.zeros(self.rho.shape)
            )
        )

        self.Lz_Be = namedtuple('Lz_Be', 's t ddT')(
            np.zeros(self.rho.shape),
            np.zeros(self.rho.shape),
            namedtuple('Lz_ddT', 's t')(
                np.zeros(self.rho.shape),
                np.zeros(self.rho.shape)
            )
        )

        self.Lz_W = namedtuple('Lz_W', 's t ddT')(
            np.zeros(self.rho.shape),
            np.zeros(self.rho.shape),
            namedtuple('Lz_ddT', 's t')(
                np.zeros(self.rho.shape),
                np.zeros(self.rho.shape)
            )
        )

        self.Lz_Ne = namedtuple('Lz_Ne', 's t ddT')(
            np.zeros(self.rho.shape),
            np.zeros(self.rho.shape),
            namedtuple('Lz_ddT', 's t')(
                np.zeros(self.rho.shape),
                np.zeros(self.rho.shape)
            )
        )

        self.Lz_Ar = namedtuple('Lz_Ar', 's t ddT')(
            np.zeros(self.rho.shape),
            np.zeros(self.rho.shape),
            namedtuple('Lz_ddT', 's t')(
                np.zeros(self.rho.shape),
                np.zeros(self.rho.shape)
            )
        )

        self.Lz_Kr = namedtuple('Lz_Kr', 's t ddT')(
            np.zeros(self.rho.shape),
            np.zeros(self.rho.shape),
            namedtuple('Lz_ddT', 's t')(
                np.zeros(self.rho.shape),
                np.zeros(self.rho.shape)
            )
        )

        # calculate cross sections on the main computational grid
        svfus_dd, svfus_dd_ddT = calc_svfus(self.T, mode='dd')
        svfus_dt, svfus_dt_ddT = calc_svfus(self.T, mode='dt')
        svrec_st, svrec_st_ddT = calc_svrec_st(self.n, self.T)
        svcx_st, svcx_st_ddT = calc_svcx_st(self.T)
        svion_st, svion_st_ddT = calc_svion_st(self.T)
        svel_st, svel_st_ddT = calc_svel_st(self.T)

        self.sv = namedtuple('sv', 'fus rec cx ion el')(
            namedtuple('sv_fus', 'dd dt d_dT')(
                svfus_dd,
                svfus_dt,
                namedtuple('sv_fus_d_dT', 'dd dt')(
                    svfus_dd_ddT,
                    svfus_dt_ddT)),
            namedtuple('sv_rec', 'st d_dT')(
                svrec_st,
                namedtuple('sv_rec_d_dT', 'st')(
                    svrec_st_ddT)),
            namedtuple('sv_cx', 'st d_dT')(
                svcx_st,
                namedtuple('sv_cx_d_dT', 'st')(
                    svcx_st_ddT)),
            namedtuple('sv_ion', 'st d_dT')(
                svion_st,
                namedtuple('sv_ion_d_dT', 'st')(
                    svion_st_ddT)),
            namedtuple('sv_el', 'st d_dT')(
                svel_st,
                namedtuple('sv_el_d_dT', 'st')(
                    svel_st_ddT)))

        # Calculate some 1D Flux-surface averaged quantities
        self.izn_rate_fsa_s = calc_fsa(self.izn_rate.s, self.R, self.Z)
        self.izn_rate_fsa_t = calc_fsa(self.izn_rate.t, self.R, self.Z)
        self.izn_rate_fsa = calc_fsa(self.izn_rate.s + self.izn_rate.t, self.R, self.Z)
        self.cool_rate_fsa = calc_fsa(self.cool_rate, self.R, self.Z)
        self.dn_dr = calc_grad(self.rho, self.n.n.s + self.n.n.t, self.psi_norm, self.R, self.Z, self.psi_data)
        self.z_eff_fsa = calc_fsa(self.z_eff, self.R, self.Z)
        self.B_p_fsa = calc_fsa(self.B_p, self.R, self.Z)
        self.B_t_fsa = calc_fsa(self.B_t, self.R, self.Z)
        self.E_r_fsa = calc_fsa(self.E_r, self.R, self.Z)  # Piper changes: Also get FSA'd Er.

        self.n_fsa = namedtuple('n', 'i e n C')(
            calc_fsa(self.n.i, self.R, self.Z),
            calc_fsa(self.n.e, self.R, self.Z),
            namedtuple('nn', 's t tot')(
                calc_fsa(self.n.n.s, self.R, self.Z),  # slow
                calc_fsa(self.n.n.t, self.R, self.Z),  # thermal
                calc_fsa(self.n.n.tot, self.R, self.Z)  # total
            ),
            calc_fsa(self.n.C, self.R, self.Z))

        Ti_kev_fsa = calc_fsa(self.T.i.kev, self.R, self.Z)
        Te_kev_fsa = calc_fsa(self.T.e.kev, self.R, self.Z)
        Tns_kev_fsa = calc_fsa(self.T.n.s.kev, self.R, self.Z)
        Tnt_kev_fsa = calc_fsa(self.T.n.t.kev, self.R, self.Z)
        TC_kev_fsa = calc_fsa(self.T.C.kev, self.R, self.Z)

        self.T_fsa = namedtuple('T', 'i e n C')(
            namedtuple('Ti', 'kev ev J')(
                Ti_kev_fsa,
                Ti_kev_fsa * 1E3,
                Ti_kev_fsa * 1E3 * e),
            namedtuple('Te', 'kev ev J')(
                Te_kev_fsa,
                Te_kev_fsa * 1E3,
                Te_kev_fsa * 1E3 * e),
            namedtuple('Tn', 's t')(
                namedtuple('Tn_s', 'kev ev J')(
                    Tns_kev_fsa,
                    Tns_kev_fsa * 1E3,
                    Tns_kev_fsa * 1E3 * e),
                namedtuple('Tn_t', 'kev ev J')(
                    Tnt_kev_fsa,
                    Tnt_kev_fsa * 1E3,
                    Tnt_kev_fsa * 1E3 * e)
            ),
            namedtuple('TC', 'kev ev J')(
                TC_kev_fsa,
                TC_kev_fsa * 1E3,
                TC_kev_fsa * 1E3 * e))
        # Piper Changes: Calculate FSA of pressures and make namedtuple.
        pi_fsa = calc_fsa(self.n.i * self.T.i.J, self.R, self.Z)
        pC_fsa = calc_fsa(self.n.C * self.T.C.J, self.R, self.Z)
        pe_fsa = calc_fsa(self.n.e * self.T.e.J, self.R, self.Z)

        self.p_fsa = namedtuple('p_fsa', 'i C e')(
            -pi_fsa,
            -pC_fsa,
            -pe_fsa)

        # Piper Changes: Calculate FSA of pressure gradient.
        dpi_dr_fsa = calc_fsa(self.dpi_dr, self.R, self.Z)
        dpC_dr_fsa = calc_fsa(self.dpC_dr, self.R, self.Z)
        dpe_dr_fsa = calc_fsa(self.dpe_dr, self.R, self.Z)

        self.dp_dr_fsa = namedtuple('dp_dr_fsa', 'i C e')(
            -dpi_dr_fsa,
            -dpC_dr_fsa,
            -dpe_dr_fsa)

        # FSA of gradient scale length.
        self.Lp_fsa = namedtuple('Lp_fsa', 'i C e')(
            -dpi_dr_fsa / pi_fsa,
            -dpC_dr_fsa / pC_fsa,
            -dpe_dr_fsa / pe_fsa)

        self.L_fsa = namedtuple('L', 'T n p')(
            namedtuple('T', 'i e')(
                calc_fsa(self.L.T.i, self.R, self.Z),
                calc_fsa(self.L.T.e, self.R, self.Z)),
            namedtuple('n', 'i e')(
                calc_fsa(self.L.n.i, self.R, self.Z),
                calc_fsa(self.L.n.e, self.R, self.Z)),
            namedtuple('p', 'i e')(
                calc_fsa(self.L.p.i, self.R, self.Z),
                calc_fsa(self.L.p.e, self.R, self.Z))
        )
        # initialize chi_r using Bohm diffusion. This might get updated later
        chi_bohm = 5 / 32 * self.T.i.J / (e * self.B_tot)

        # calculate chi using the JET Bohm-GyroBohm model
        chi_i_jet, chi_e_jet = calc_chi_jet(self.T, self.L, self.a, self.q, self.B_t, m_d, self.rho)

        # chi hybrid
        # this chi uses the JET chi up to rho = 0.9 and then linearly interpolates to the value of chi_bohm
        # at the seperatrix

        self.chi = namedtuple('chi', 'bohm jet')(
            chi_bohm,
            namedtuple('jet', 'i e')(
                chi_i_jet,
                chi_e_jet
            )
        )
        # TODO: Verify this fus_rate calculation. I think it's wrong.
        self.fus_rate = 1 / 4 * self.n.D * self.n.D * self.sv.fus.dd + 1 / 4 * self.n.D * self.n.T * self.sv.fus.dt

    def plot_te(self, edge=False):
        """
        Plots the GT3.Core electron temperature
        """
        fig = self._plot_base(self.T_fsa.e.kev,yLabel=r'$T_e [keV]$', title="GT3.Core Electron Temperature", edge=edge)
        return fig

    def _plot_base(self, val, xLabel=r'$\rho$', yLabel="Value", title="Title", color='red', edge=False):

        plot = plt.figure()
        fig = plot.add_subplot(111)
        fig.set_xlabel(xLabel, fontsize=20)
        fig.set_ylabel(yLabel, fontsize=20)
        fig.set_title(title)
        if edge:
            fig.set_xlim(0.85, self.sep_val)
        fig.scatter(self.rho[:, 0], val, color=color, s=4)
        plt.show()
        return fig

    def _plot_with_wall(self):
        """
        Generates a Matplotlib plot with the wall pre-plotted.

        :return: An Axis object with the wall line plotted
        """
        fig_width = 6.0
        # Check to see if self.Z/R have been defind yet to generate a figure height, as they aren't immediately
        # calculated
        try:
            fig_height = (np.amax(self.Z) - np.amin(self.Z)) / (np.amax(self.R) - np.amin(self.R)) * fig_width
        except:
            fig_height = 9.0

        fig = plt.figure(figsize=(0.975*fig_width, fig_height))
        ax1 = fig.add_subplot(1, 1, 1)
        ax1.axis('equal')
        ax1.plot(np.asarray(self.wall_line.xy)[0], np.asarray(self.wall_line.xy)[1], color='black', lw=1.5)
        return ax1

    def _shapely_obj_plot_hanlder(self, obj, ax):
        if isinstance(obj, Point):
            ax.scatter(obj.x, obj.y, color='red', marker='o')
            return ax
        if isinstance(obj, MultiPoint):
            for p in obj:
                ax.scatter(p.x, p.y, color='red', marker='o')
            return ax
        if isinstance(obj, Path):
            ls = LineString(obj.vertices)
            ax.plot(*ls.xy, color='red', marker='o')
            return ax
        if isinstance(obj, LineString):
            ax.plot(obj.xy, color='red', marker='o')
            return ax

    def _unknown_data_plot_helper(self, obj, ax):
        # Is it a shapely or similar object?
        if isinstance(obj, (Point, MultiPoint, Path, LineString)):
            return self._shapely_obj_plot_hanlder(obj, ax)

        # Does it have the vertices property to convert to a LineString?
        try:
            ls = LineString(obj.vertices)
            return self._shapely_obj_plot_hanlder(ls, ax)
        except:
            pass

        # Does it have an xy property?

        try:
            ax.plot(*obj.xy)
            return ax
        except:
            pass

        # Can I turn it into a point?
        try:
            p = Point(obj)
            return self._shapely_obj_plot_hanlder(p, ax)
        except:
            pass

        # Can I turn it into a LineString directly?
        try:
            ls = LineString(obj)
            return self._shapely_obj_plot_hanlder(ls, ax)
        except:
            pass

        # No dice. Raising error.
        raise

    def plot_with_wall(self, obj=None):
        ax = self._plot_with_wall()

        try:
            return self._unknown_data_plot_helper(obj, ax)
        except:
            pass

        # Is it an iterable?
        try:
            obj.__iter__
            for p in obj:
                try:
                    ax = self._unknown_data_plot_helper(p, ax)
                except:
                    raise
            return ax
        except:
            pass

        print "Could not plot given data"

    def plot_exp_psi(self, res=50):
        ax = self._plot_with_wall()
        try:
            ax.contour(self.R, self.Z, self.psi, res)
            return ax
        except NameError:
            print "Psi not defined"
            pass

    def plot_L_t_i(self, edge=False):
        """
        Plots the GT3.Core ion temperature scale length
        """
        fig = self._plot_base(self.L.T.i[:,0], yLabel=r'$L_{T,i} [m]$', title="GT3.Core Ion Temperature Scale Length", edge=edge)
        return fig

    def plot_L_t_e(self, edge=False):
        """
        Plots the GT3.Core electron temperature scale length
        """
        fig = self._plot_base(self.L.T.e[:,0], yLabel=r'$L_{T,e} [m]$', title="GT3.Core Electron Temperature Scale Length", edge=edge)
        return fig

    def plot_L_n_i(self, edge=False):
        """
        Plots the GT3.Core ion density scale length
        """
        fig = self._plot_base(self.L.n.i[:,0], yLabel=r'$L_{n,i} [m]$', title="GT3.Core Ion Density Scale Length", edge=edge)
        return fig

    def plot_L_n_e(self, edge=False):
        """
        Plots the GT3.Core electron density scale length
        """
        fig = self._plot_base(self.L.n.e[:,0], yLabel=r'$L_{n,e} [m]$', title="GT3.Core Electron Density Scale Length", edge=edge)
        return fig

    def plot_ti(self, edge=False):
        """
        Plots the GT3.Core ion temperature
        """
        fig = self._plot_base(self.T_fsa.e.kev, yLabel=r'$T_i [keV]$', title="GT3.Core Ion Temperature", edge=edge)
        return fig

    def plot_t(self, edge=False):
        """
        Plots the 1D GT3.Core ion and electron temperatures
        """
        fig = self._plot_base(self.T_fsa.e.kev, yLabel=r'$T [keV]$', title="GT3.Core ion/electron temperatures", edge=edge)
        fig.scatter(self.rho[:, 0], self.T_fsa.i.kev, color="blue")
        fig.legend([r'$T_e$', r'$T_i$'])
        return fig

    def plot_ni(self, edge=False):
        """
        Plots the 1D GT3.Core ion density
        """
        fig = self._plot_base(self.n_fsa.i, yLabel=r'$n_i [#/m^3]$', title="GT3.Core Ion density", edge=edge)
        return fig

    def plot_ne(self, edge=False):
        """
        Plots the 1D GT3.Core electron density
        """
        fig = self._plot_base(self.n_fsa.e, yLabel=r'$n_e [#/m^3]$', title="GT3.Core electron density", edge=edge)
        return fig

    def plot_n(self, edge=False):
        """
        Plots the 1D GT3.Core ion and electron densities
        """
        fig = self._plot_base(self.n_fsa.e, yLabel=r'$n [#/m^3]$', title="GT3.Core ion/electron densities", edge=edge)
        fig.scatter(self.rho[:, 0], self.n_fsa.i, color="blue")
        fig.legend([r'$n_e$', r'$n_i$'])
        return fig

    def plot_ntrl_s(self, edge=False):
        """
        Plots the 1D GT3.Core slow neutral density
        """
        fig = self._plot_base(self.n_fsa.n.s, yLabel=r'$n_{o,s} [#/m^3]$', title="GT3.Core Slow Neutral Density", edge=edge)
        return fig

    def plot_ntrl_s(self, edge=False):
        """
        Plots the 1D GT3.Core thermal neutral density
        """
        fig = self._plot_base(self.n_fsa.n.t, yLabel=r'$n_{o,t} [#/m^3]$', title="GT3.Core Thermal Neutral Density", edge=edge)
        return fig

    def plot_ntrl_total(self, edge=False):
        """
        Plots the 1D GT3.Core total neutral density
        """
        fig = self._plot_base(self.n_fsa.n.tot, yLabel=r'$n_{o,tot} [#/m^3]$', title="GT3.Core Total Neutral Density", edge=edge)
        return fig

    def plot_er(self, edge=False):
        """
        Plots the 1D GT3.Core electric field
        """
        fig = self._plot_base(self.E_r_fsa, yLabel=r'$E_r [V/m]$', title="GT3.Core Radial Electric Field",
                              edge=edge)
        return fig

    def plot_vpol_C(self, edge=False):
        """
        Plots the 1D GT3.Core carbon poloidal velocity
        """

        fig = self._plot_base(self.v_1D.pol.C, yLabel=r'$V_{C,\theta} [km/s]$',
                              title="GT3.Core Carbon Poloidal Rotation Velocity",
                              edge=edge)
        return fig

    def plot_vpol_D(self, edge=False):
        """
        Plots the 1D GT3.Core deuterium poloidal velocity
        """

        fig = self._plot_base(self.v_1D.pol.D, yLabel=r'$V_{C,\theta} [km/s]$',
                              title="GT3.Core Deuterium Poloidal Rotation Velocity",
                              edge=edge)
        return fig

    def plot_vtor_C(self, edge=False):
        """
        Plots the 1D GT3.Core carbon toroidal velocity
        """
        fig = self._plot_base(self.v_1D.tor.C, yLabel=r'$V_{c,\phi} [km/s]$',
                              title="GT3.Core Carbon toroidal velocity",
                              edge=edge)
        return fig

    def plot_vtor_D(self, edge=False):
        """
        Plots the 1D GT3.Core deuterium toroidal velocity
        """

        fig = self._plot_base(self.v_1D.tor.D, yLabel=r'$V_{c,\phi} [km/s]$',
                              title="GT3.Core deuterium toroidal velocity",
                              edge=edge)
        return fig

    def plot_pressure_C(self, edge=False):
        """
        Plots the 1D GT3.Core carbon pressure
        """

        fig = self._plot_base(self.p_fsa.C, yLabel=r'$P_C [??]$',
                              title="GT3.Core carbon pressure",
                              edge=edge)
        return fig

    def plot_pressure_D(self, edge=False):
        """
        Plots the 1D GT3.Core deuterium pressure
        """

        fig = self._plot_base(self.p_fsa.i, yLabel=r'$P_i [??]$',
                              title="GT3.Core deuterium pressure",
                              edge=edge)
        return fig

    def plot_pressure_e(self, edge=False):
        """
        Plots the 1D GT3.Core electron pressure
        """

        fig = self._plot_base(self.p_fsa.e, yLabel=r'$P_e [??]$',
                              title="GT3.Core electron pressure",
                              edge=edge)
        return fig

    def plot_izn_rate_s(self, edge=False):
        """
        Plots the 1D GT3.Core slow neutral ionization rate
        """

        fig = self._plot_base(self.izn_rate_fsa_s, yLabel=r'$<\sigma v>_{slow} [??]$',
                              title="GT3.Core slow neutral ionization rate",
                              edge=edge)
        return fig

    def plot_izn_rate_t(self, edge=False):
        """
        Plots the 1D GT3.Core thermal neutral ionization rate
        """

        fig = self._plot_base(self.izn_rate_fsa_t, yLabel=r'$<\sigma v>_{thermal} [??]$',
                              title="GT3.Core thermal neutral ionization rate",
                              edge=edge)
        return fig

    def plot_izn_rate_total(self, edge=False):
        """
        Plots the 1D GT3.Core total neutral ionization rate
        """

        fig = self._plot_base(self.izn_rate_fsa, yLabel=r'$<\sigma v>_{tot} [??]$',
                              title="GT3.Core total neutral ionization rate",
                              edge=edge)
        return fig

    def plot_q(self, edge=False):
        """
        Plots the 1D GT3.Core safety factor
        """

        fig = self._plot_base(self.q_1D, yLabel=r'$q$',
                              title="GT3.Core safety factor",
                              edge=edge)
        return fig

    def update_ntrl_data(self, data):
        try:
            n_n_s = griddata(np.column_stack((data.R, data.Z)),
                             data.n_n_slow,
                             (self.R, self.Z),
                             method='linear')
        except:
            n_n_s = self.n.n.s

        try:
            n_n_t = griddata(np.column_stack((data.R, data.Z)),
                             data.n_n_thermal,
                             (self.R, self.Z),
                             method='linear')
        except:
            n_n_t = self.n.n.t

        try:
            izn_rate_s = griddata(np.column_stack((data.R, data.Z)),
                                  data.izn_rate_slow,
                                  (self.R, self.Z),
                                  method='linear')
        except:
            izn_rate_s = self.izn_rate.s

        try:
            izn_rate_t = griddata(np.column_stack((data.R, data.Z)),
                                  data.izn_rate_thermal,
                                  (self.R, self.Z),
                                  method='linear')
        except:
            izn_rate_t = self.izn_rate.t

        nn = namedtuple('nn', 's t tot')(
            n_n_s,  # slow
            n_n_t,  # thermal
            n_n_s + n_n_t  # total
        )

        self.n = namedtuple('n', 'i e n C')(self.n.i, self.n.e, nn, self.n.C)

        self.n_fsa = namedtuple('n', 'i e n C')(
            calc_fsa(self.n.i, self.R, self.Z),
            calc_fsa(self.n.e, self.R, self.Z),
            namedtuple('nn', 's t tot')(
                calc_fsa(self.n.n.s, self.R, self.Z),  # slow
                calc_fsa(self.n.n.t, self.R, self.Z),  # thermal
                calc_fsa(self.n.n.tot, self.R, self.Z)  # total
            ),
            calc_fsa(self.n.C, self.R, self.Z))

        self.izn_rate = namedtuple('izn_rate', 's t tot')(
            izn_rate_s,  # slow
            izn_rate_t,  # thermal
            izn_rate_s + izn_rate_t  # total
        )

        self.izn_rate_fsa_s = calc_fsa(self.izn_rate.s, self.R, self.Z)
        self.izn_rate_fsa_t = calc_fsa(self.izn_rate.t, self.R, self.Z)
        self.izn_rate_fsa = calc_fsa(self.izn_rate.s + self.izn_rate.t, self.R, self.Z)

        self.cool_rate_fsa = calc_fsa(self.n.e * self.n.C * np.nan_to_num(self.Lz_C.s) + \
                                      + self.n.e * self.n.C * np.nan_to_num(self.Lz_C.t),
                                      self.R, self.Z)
        self.dn_dr_fsa = calc_fsa(self.dn_dr, self.R, self.Z)

        # self.izn_rate_fsa = np.array(map(lambda x: self.izn_rate_fsa[-1] * x ** 10, self.r[:, 0] / self.a))
        # self.cool_rate_fsa = np.array(map(lambda x: self.cool_rate_fsa[-1] * x ** 10, self.r[:, 0] / self.a))
        # self.dn_dr_fsa = np.array(map(lambda x: self.dn_dr_fsa[-1] * x ** 10, self.r[:, 0] / self.a))

    def update_Lz_data(self, z, Lz, dLzdT):

        Lz_slow = Lz(np.log10(self.T.n.s.kev),
                     np.log10(self.n.n.s / self.n.e),
                     np.log10(self.T.e.kev))

        dLzdT_slow = dLzdT(np.log10(self.T.n.s.kev),
                           np.log10(self.n.n.s / self.n.e),
                           np.log10(self.T.e.kev))

        Lz_thermal = Lz(np.log10(self.T.n.t.kev),
                        np.log10(self.n.n.t / self.n.e),
                        np.log10(self.T.e.kev))

        dLzdT_thermal = dLzdT(np.log10(self.T.n.t.kev),
                              np.log10(self.n.n.t / self.n.e),
                              np.log10(self.T.e.kev))

        if z == 6:
            self.Lz_C = namedtuple('Lz_C', 's t ddT')(
                Lz_slow,
                Lz_thermal,
                namedtuple('Lz_ddT', 's t')(
                    dLzdT_slow,
                    dLzdT_thermal
                )
            )

        if z == 4:
            self.Lz_Be = namedtuple('Lz_Be', 's t ddT')(
                Lz_slow,
                Lz_thermal,
                namedtuple('Lz_ddT', 's t')(
                    dLzdT_slow,
                    dLzdT_thermal
                )
            )

        if z == 74:
            self.Lz_W = namedtuple('Lz_W', 's t ddT')(
                Lz_slow,
                Lz_thermal,
                namedtuple('Lz_ddT', 's t')(
                    dLzdT_slow,
                    dLzdT_thermal
                )
            )

        if z == 10:
            self.Lz_Ne = namedtuple('Lz_Ne', 's t ddT')(
                Lz_slow,
                Lz_thermal,
                namedtuple('Lz_ddT', 's t')(
                    dLzdT_slow,
                    dLzdT_thermal
                )
            )

        if z == 18:
            self.Lz_Ar = namedtuple('Lz_Ar', 's t ddT')(
                Lz_slow,
                Lz_thermal,
                namedtuple('Lz_ddT', 's t')(
                    dLzdT_slow,
                    dLzdT_thermal
                )
            )

        if z == 36:
            self.Lz_Kr = namedtuple('Lz_Kr', 's t ddT')(
                Lz_slow,
                Lz_thermal,
                namedtuple('Lz_ddT', 's t')(
                    dLzdT_slow,
                    dLzdT_thermal
                )
            )

        self.cool_rate_fsa = calc_fsa(self.n.e * self.n.C * np.nan_to_num(self.Lz_C.s) + \
                                      + self.n.e * self.n.C * np.nan_to_num(self.Lz_C.t),
                                      self.R, self.Z)  # TODO: Figure out how to not run this for all Z in a normal shot
