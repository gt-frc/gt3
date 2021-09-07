#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 18 13:22:31 2018

@author: max
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata, UnivariateSpline
from math import pi
from collections import namedtuple
from .Functions.CalcSV import calc_svion_st, calc_svel_st, calc_svrec_st, calc_svcx_st, calc_svfus
from .Functions.CalcFsPerimInt import calc_fs_perim_int
from .Functions.CreateSurfAreaInterp import create_surf_area_interp
from .Functions.CalcTheta1D import calc_theta1d
from .Functions.CalcPtsLines import calc_pts_lines
from .Functions.CalcRZ import calc_RZ
from .Functions.CalcKappaElong import calc_kappa_elong
from .Functions.CreateVolInterp import create_vol_interp
from .Functions.CalcChiJet import calc_chi_jet
from scipy.interpolate import interp1d
from GT3.Psi import Psi
import GT3.constants as constants
from GT3.utilities import PlotBase
from warnings import warn
from .Functions.ProfileClasses import SlowFastSplit, TwoDProfile, ImpurityProfiles, TemperatureProfiles, \
    DensityProfiles, PressureProfiles, VectorialProfiles, VectorialBase, TwoDProfileWithHM

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
class Core(PlotBase.PlotBase):

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
    num_xpt = None

    def __init__(self, inp, *args, **kwargs):
        # Hold the input as a class attribute for use in functions

        super().__init__()
        if kwargs.get("debug"):
            self.debug = True
        else:
            self.debug = False

        self.inp = inp
        self.wall_line = inp.wall_line
        self.sep_val = inp.sep_val
        if abs(1.0 - float(self.sep_val)) > 0.0001:
            self.sep_val_overridden = True
            print("The separatrix value is being overwritten. GT3.SOL will be omitted.")

        # Calculate PsiData information
        self._set_psiData(inp, **kwargs)
        self.a = self.psi.a


        self.R0_a = self.psi.get_mag_axis()[0]
        self.R0_g = self.psi.get_geom_axis()[0]


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

        theta1d, theta_markers = calc_theta1d(self.psi, inp.thetapts_approx)
        theta_xpt = theta_markers[3]
        self.thetapts = len(theta1d)

        # create rho and theta arrays (initializing the main computational grid)
        self.theta, self.rho = np.meshgrid(theta1d, rho1d)
        self.psi.set_rho_mesh(self.rho)

        self.set_plot_rho1d(self.rho[:, 0])

        # estimate elongation and triangularity
        self.kappa_vals, self.tri_vals = calc_kappa_elong(self.psi, np.asarray(self.psi.lcfs_line_closed.coords))



        # Create Psi object containing lots of important information

        self.shaf_shift = (self.psi.mag_axis[0] - self.psi.geo_axis[0]) / self.a
        self.r = self.rho * self.a

        self.R, self.Z = calc_RZ(self.rho, self.theta, theta_xpt, self.psi, self.psi.psi_norm)
        self.set_plot_RZ(self.R, self.Z)


        self.xpt_loc = np.where(self.Z == np.amin(self.Z))
        self.obmp_loc = np.where(self.R == np.amax(self.R))
        self.ibmp_loc = np.where(self.R == np.amin(self.R))
        self.top_loc = np.where(self.Z == np.amax(self.Z))



        # create interpolation functions to obtain the flux surface surface area for any value of r, rho, or psi_norm
        self.r2sa, self.rho2sa, self.psinorm2sa = create_surf_area_interp(self.psi.rho2psinorm,
                                                                          self.psi.psinorm2rho,
                                                                          self.psi,
                                                                          np.asarray(self.psi.get_lcfs_closed().coords),
                                                                          self.R0_a,
                                                                          self.a,
                                                                          self.sep_val)

        # create interpolation functions to obtain the plasma volume for any value of r, rho, or psi_norm
        self.r2vol, self.rho2vol, self.psinorm2vol = create_vol_interp(self.psi.rho2psinorm,
                                                                       self.psi.psinorm2rho,
                                                                       self.psi,
                                                                       np.asarray(self.psi.get_lcfs_closed().coords),
                                                                       self.R0_a,
                                                                       self.a,
                                                                       self.sep_val)

        #np.savetxt('psi2rho.txt', np.column_stack((np.linspace(0, self.sep_val, 100),
        #                                           self.psi.psi2rho(np.linspace(0, self.sep_val, 100)))))

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
        self.izn_rate = SlowFastSplit(TwoDProfileWithHM(self.psi, np.zeros(self.rho.shape), self.R, self.Z),
                                      TwoDProfileWithHM(self.psi, np.zeros(self.rho.shape), self.R, self.Z),
                                      TwoDProfileWithHM(self.psi, np.zeros(self.rho.shape), self.R, self.Z))

        self.izn_rate.tot.set_wall(self.wall_line)
        self.izn_rate.s.set_wall(self.wall_line)
        self.izn_rate.t.set_wall(self.wall_line)

        # initialize cooling rate array with zero
        self.cool_rate = TwoDProfile(self.psi, np.zeros(self.rho.shape), self.R, self.Z, wall=self.wall_line)

        # Initialize densities
        self._set_densities(inp)

        # Initialize temperatures
        self._set_temperatures(inp)

        # Initialize pressures
        self.p = PressureProfiles(self.psi, self.R, self.Z,
                                  i=self.n.i * self.T.i.J,
                                  e=self.n.e * self.T.e.J,
                                  C=self.n.C * self.T.C.J,
                                  wall=self.wall_line)

        self._set_efield(inp)

        self._set_velocities(inp)

        self._set_bfields(inp)


        try:
            # use input q profile if given
            q_interp = UnivariateSpline(inp.q_data[:, 0], inp.q_data[:, 1], k=3, s=1.0)
            q = q_interp(self.rho)
            self.q = TwoDProfile(self.psi, q, self.R, self.Z, wall=self.wall_line)
        except:
            # otherwise calculate q-profile from psi data
            q_1D = inp.BT0 * self.psi.mag_axis[0] / (2 * pi) * calc_fs_perim_int(1.0 / (self.R ** 2 * self.B.pol),
                                                                                 self.R, self.Z)
            q_1D[0] = q_1D[1]
            q = np.repeat(q_1D[np.newaxis, :], self.rho.shape[1], axis=0).T
            self.q = TwoDProfile(self.psi, q, self.R, self.Z, wall=self.wall_line)
        self.q0 = self.q.fsa.val[0]
        self.q95 = self.q.fsa.Spline(0.95*self.a)

        # create Lz-related variables. These are initiated as np.zeros until updated by ImpRad module
        self.Lz = ImpurityProfiles(core=self)

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


        # initialize chi_r using Bohm diffusion. This might get updated later
        chi_bohm = 5 / 32 * self.T.i.J / (e * self.B.tot)

        # calculate chi using the JET Bohm-GyroBohm model
        chi_i_jet, chi_e_jet = calc_chi_jet(self.T, self.p, self.a, self.q, self.B.tor, m_d, self.rho)

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
        self.fus_rate = 1 / 4 * self.n.i * self.n.i * self.sv.fus.dd + 1 / 4 * self.n.i * self.n.T * self.sv.fus.dt
        plt.clf()
        plt.close()

    def _set_densities(self, inp):
        # initialize main density object
        raw_data = {}
        # electron density
        try:
            ne = UnivariateSpline(inp.ne_data[:, 0], inp.ne_data[:, 1], k=3, s=0)(self.rho)
            raw_data['e'] = [inp.ne_data[:, 0], inp.ne_data[:, 1]]
        except (AttributeError, TypeError):
            ne = np.zeros(self.rho.shape)

        # carbon density
        try:
            nC = UnivariateSpline(inp.nC_data[:, 0], inp.nC_data[:, 1], k=3, s=0)(self.rho)
            raw_data['C'] = [inp.nC_data[:, 0], inp.nC_data[:, 1]]
        except (AttributeError, TypeError):
            try:
                frac_C = UnivariateSpline(inp.frac_C_data[:, 0], inp.frac_C_data[:, 1], k=3, s=0)(self.rho)
                nC = ne * frac_C
            except (AttributeError, TypeError):
                # Set the carbon density to 2.5% across the board
                nC = 0.025 * ne / (1. + .025 * 6.0)
                # nC = np.zeros(self.rho.shape)

            # deuterium density
        try:
            nD = UnivariateSpline(inp.nD_data[:, 0], inp.nD_data[:, 1], k=3, s=0)(self.rho)
            raw_data['D'] = [inp.nD_data[:, 0], inp.nD_data[:, 1]]
        except (AttributeError, TypeError):
            nD = np.zeros(self.rho.shape)

            # tritium density
        try:
            nT = UnivariateSpline(inp.nT_data[:, 0], inp.nT_data[:, 1], k=3, s=0)(self.rho)
            raw_data['T'] = [inp.nT_data[:, 0], inp.nT_data[:, 1]]
        except (AttributeError, TypeError):
            nT = np.zeros(self.rho.shape)

        # If no ion density data are given
        if not (nD + nT).any():
            nD = ne / (1. + .025 * 6.0)

        # tungsten density
        try:
            nW = UnivariateSpline(inp.nW_data[:, 0], inp.nW_data[:, 1], k=3, s=0)(self.rho)
            raw_data['W'] = [inp.nW_data[:, 0], inp.nW_data[:, 1]]
        except (AttributeError, TypeError):
            try:
                frac_W = UnivariateSpline(inp.frac_W_data[:, 0], inp.frac_W_data[:, 1], k=3, s=0)(self.rho)
                nW = ne * frac_W
            except (AttributeError, TypeError):
                nW = np.zeros(self.rho.shape)

        # beryllium density
        try:
            nBe = UnivariateSpline(inp.nBe_data[:, 0], inp.nBe_data[:, 1], k=3, s=0)(self.rho)
            raw_data['Be'] = [inp.nBe_data[:, 0], inp.nBe_data[:, 1]]
        except (AttributeError, TypeError):
            try:
                frac_Be = UnivariateSpline(inp.frac_Be_data[:, 0], inp.frac_Be_data[:, 1], k=3, s=0)(self.rho)
                nBe = ne * frac_Be
            except (AttributeError, TypeError):
                nBe = np.zeros(self.rho.shape)

        # neon density
        try:
            nNe = UnivariateSpline(inp.nNe_data[:, 0], inp.nNe_data[:, 1], k=3, s=0)(self.rho)
            raw_data['Ne'] = [inp.nNe_data[:, 0], inp.nNe_data[:, 1]]
        except (AttributeError, TypeError):
            try:
                frac_Ne = UnivariateSpline(inp.frac_Ne_data[:, 0], inp.frac_Ne_data[:, 1], k=3, s=0)(self.rho)
                nNe = ne * frac_Ne
            except (AttributeError, TypeError):
                nNe = np.zeros(self.rho.shape)

        # krypton density
        try:
            nKr = UnivariateSpline(inp.nKr_data[:, 0], inp.nKr_data[:, 1], k=3, s=0)(self.rho)
            raw_data['Kr'] = [inp.nKr_data[:, 0], inp.nKr_data[:, 1]]
        except (AttributeError, TypeError):
            try:
                frac_Kr = UnivariateSpline(inp.frac_Kr_data[:, 0], inp.frac_Kr_data[:, 1], k=3, s=0)(self.rho)
                nKr = ne * frac_Kr
            except (AttributeError, TypeError):
                nKr = np.zeros(self.rho.shape)

        # argon density
        try:
            nAr = UnivariateSpline(inp.nAr_data[:, 0], inp.nAr_data[:, 1], k=3, s=0)(self.rho)
            raw_data['Ar'] = [inp.nAr_data[:, 0], inp.nAr_data[:, 1]]
        except (AttributeError, TypeError):
            try:
                frac_Ar = UnivariateSpline(inp.frac_Ar_data[:, 0], inp.frac_Ar_data[:, 1], k=3, s=0)(self.rho)
                nAr = ne * frac_Ar
            except (AttributeError, TypeError):
                nAr = np.zeros(self.rho.shape)

        # alpha density
        try:
            na = UnivariateSpline(inp.na_data[:, 0], inp.na_data[:, 1], k=3, s=0)(self.rho)
            raw_data['alpha'] = [inp.na_data[:, 0], inp.na_data[:, 1]]
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
        self.z_0 = nC * 6.0 ** 2 / (nD + nT)  # TODO: Update this calculation
        zeff = ((nD + nT) * 1.0 ** 2 + nC * 6.0 ** 2) / ne
        self.z_eff = TwoDProfile(self.psi, zeff, self.R, self.Z)  # TODO: Update this calculation

        self.n = DensityProfiles(self.psi, self.R, self.Z,
                                 i=nD,
                                 e=ne,
                                 C=nC,
                                 T=nT,
                                 W=nW,
                                 Be=nBe,
                                 Ne=nNe,
                                 Ar=nAr,
                                 Kr=nKr,
                                 alpha=na,
                                 ns=nn.s,
                                 nt=nn.t,
                                 wall=self.wall_line, raw=raw_data)

    def _set_temperatures(self, inp):

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

        raw_data = {
            'i': [inp.Ti_data[:, 0], inp.Ti_data[:, 1]],
            'e': [inp.Te_data[:, 0], inp.Te_data[:, 1]],
        }

        try:
            TC_kev = UnivariateSpline(inp.TC_data[:, 0], inp.TC_data[:, 1], k=3, s=0.0)(self.rho)
            raw_data['C'] = [inp.TC_data[:, 0], inp.TC_data[:, 1]]
        except:
            TC_kev = Ti_kev


        self.T = TemperatureProfiles(self.psi, self.R, self.Z,
                                     i=Ti_kev,
                                     e=Te_kev,
                                     C=TC_kev,
                                     ns=Tns_kev,
                                     nt=Tnt_kev,
                                     wall=self.wall_line, raw=raw_data)

    def _set_psiData(self, inp, **kwargs):
        # this assumes that psi is given as a square array
        psi_shape = int(np.sqrt(inp.psirz_exp[:, 0].shape[0]))  # type: int
        try:
            psi_scale = inp.psi_scale
        except AttributeError:
            psi_scale = 1.

        self.psi = Psi(inp.psirz_exp[:, 0].reshape(-1, psi_shape),
                        inp.psirz_exp[:, 1].reshape(-1, psi_shape),
                        inp.psirz_exp[:, 2].reshape(-1, psi_shape) * psi_scale,
                        self.wall_line, **kwargs)


        xpt_l = self.psi.get_xpt_lower()
        xpt_u = self.psi.get_xpt_upper()
        self.mag_axis = self.psi.get_mag_axis()
        if self.debug:
            print("Lower X-Point:" + str(xpt_l))
            print("Upper X-Point:" + str(xpt_u))
        self.xpt = [xpt_l, xpt_u]


    def _set_efield(self, inp):
        try:
            try:
                E_r_fit = UnivariateSpline(inp.er_data[:, 0], inp.er_data[:, 1] * inp.Er_scale, k=3, s=0)
            except (AttributeError, TypeError):
                E_r_fit = UnivariateSpline(inp.er_data[:, 0], inp.er_data[:, 1], k=3, s=0)
            self.E_r = TwoDProfile(self.psi, E_r_fit(self.rho), self.R, self.Z, wall=self.wall_line, units=r"$V/m")
            if np.average(np.abs(self.E_r.fsa.val)) <= 100.:
                warn("Electric field average is <100V. Verify E_r data and any multipliers are correct.", UserWarning)
                print("I AM WRNIANG")
            E_pot = np.zeros(self.rho.shape)
            try:
                for i, rhoval in enumerate(self.rho[:, 0]):
                    E_pot[i] = E_r_fit.integral(rhoval, self.sep_val)
                self.E_pot = TwoDProfile(self.psi, E_pot, self.R, self.Z, wall=self.wall_line, units=r"$V/m")
            except:
                print("Error in E_pot integration. Setting to zeros")
                self.E_pot = TwoDProfile(self.psi, np.zeros(self.rho.shape), self.R, self.Z, wall=self.wall_line, units=r"$V/m")
        except:
            print('Er data not supplied. Setting E_r and E_pot to zero.')
            self.E_r = TwoDProfile(self.psi, np.zeros(self.rho.shape), self.R, self.Z, wall=self.wall_line, units=r"$V/m")
            self.E_pot = TwoDProfile(self.psi, np.zeros(self.rho.shape), self.R, self.Z, wall=self.wall_line, units=r"$V")
        # initialize E_r and the corresponding electric potential

    def _set_velocities(self, inp):
        # initialize rotation velocities from data
        try:
            vpolD = UnivariateSpline(inp.vpolD_data[:, 0], inp.vpolD_data[:, 1], k=3, s=0)(self.rho)
        except (AttributeError, TypeError):
            # Check if frequencies are given instead
            try:
                vpolD = UnivariateSpline(inp.omegDpol_data[:, 0], inp.omegDpol_data[:, 1] * self.psi.mag_axis[0], k=3, s=0)(self.rho)
            except (AttributeError, TypeError):
                vpolD = np.zeros(self.rho.shape)

        try:
            vpolC = UnivariateSpline(inp.vpolC_data[:, 0], inp.vpolC_data[:, 1], k=3, s=0)(self.rho)
        except (AttributeError, TypeError):
            # Check if frequencies are given instead
            try:
                vpolC = UnivariateSpline(inp.omegCpol_data[:, 0], inp.omegCpol_data[:, 1] * self.psi.mag_axis[0], k=3, s=0)(self.rho)
            except (AttributeError, TypeError):
                vpolC = np.zeros(self.rho.shape)

        try:
            vtorD = UnivariateSpline(inp.vtorD_data[:, 0], inp.vtorD_data[:, 1], k=3, s=0)(self.rho)
        except (AttributeError, TypeError):
            # Check if frequencies are given instead
            try:
                vtorD = UnivariateSpline(inp.omegDtor_data[:, 0], inp.omegDtor_data[:, 1] * self.psi.mag_axis[0], k=3, s=0)(self.rho)
            except (AttributeError, TypeError):
                vtorD = np.zeros(self.rho.shape)

        try:
            vtorC = UnivariateSpline(inp.vtorC_data[:, 0], inp.vtorC_data[:, 1], k=3, s=0)(self.rho)
        except (AttributeError, TypeError):
            # Check if frequencies are given instead
            try:
                vtorC = UnivariateSpline(inp.omegCtor_data[:, 0], inp.omegCtor_data[:, 1] * self.psi.mag_axis[0], k=3, s=0)(self.rho)
            except (AttributeError, TypeError):
                vtorC = np.zeros(self.rho.shape)

        self.v = VectorialProfiles(self.psi, self.R, self.Z, wall=self.wall_line, ProfileType=TwoDProfile,
                                   pol_D=vpolD,
                                   tor_D=vtorD,
                                   pol_C=vpolC,
                                   tor_C=vtorC)

    def _set_bfields(self, inp):
        # initialize magnetic field-related quantities
        B_pol_raw = np.sqrt((np.gradient(self.psi.psi_exp, self.psi.R[0], axis=1) / self.psi.R) ** 2 +
                            (np.gradient(self.psi.psi_exp, self.psi.Z.T[0], axis=0) / self.psi.R) ** 2)

        B_p = griddata(np.column_stack((self.psi.R.flatten(), self.psi.Z.flatten())),
                       B_pol_raw.flatten(),
                       (self.R, self.Z),
                       method='linear')

        B_t = inp.BT0 * self.psi.mag_axis[0] / self.R


        self.B = VectorialBase(B_t, B_p, self.psi, self.R, self.Z, self.wall_line, profileType=TwoDProfile)
        self.f_phi = self.B.tor / self.B.tot

    def update_ntrl_data(self, data):
        try:
            n_n_s = griddata(np.column_stack((data.R, data.Z)),
                             data.n_n_slow,
                             (self.R, self.Z),
                             method='linear')
            if np.isnan(n_n_s).any():
                n_n_s = griddata(np.column_stack((data.R, data.Z)), np.array(data.n_n_slow),
                                (self.R, self.Z), method='nearest')
                if np.isnan(n_n_s).any():
                    raise
        except:
            n_n_s = self.n.n.s

        try:
            n_n_t = griddata(np.column_stack((data.R, data.Z)),
                             data.n_n_thermal,
                             (self.R, self.Z),
                             method='linear')
            if np.isnan(n_n_t).any():
                n_n_t = griddata(np.column_stack((data.R, data.Z)), np.array(data.n_n_thermal),
                                  (self.R, self.Z), method='nearest')
                if np.isnan(n_n_t).any():
                    raise
        except:
            n_n_t = self.n.n.t

        try:
            izn_rate_s = griddata(np.column_stack((data.R, data.Z)),
                                  data.izn_rate_slow,
                                  (self.R, self.Z),
                                  method='linear')
            if np.isnan(izn_rate_s).any():
                izn_rate_s = griddata(np.column_stack((data.R, data.Z)), np.array(data.izn_rate_slow),
                                     (self.R, self.Z), method='nearest')
                if np.isnan(izn_rate_s).any():
                    raise
        except:
            izn_rate_s = self.izn_rate.s

        try:
            izn_rate_t = griddata(np.column_stack((data.R, data.Z)),
                                  data.izn_rate_thermal,
                                  (self.R, self.Z),
                                  method='linear')
            if np.isnan(izn_rate_t).any():
                izn_rate_t = griddata(np.column_stack((data.R, data.Z)), np.array(data.izn_rate_thermal),
                                  (self.R, self.Z), method='nearest')
                if np.isnan(izn_rate_t).any():
                    raise
        except:
            izn_rate_t = self.izn_rate.t

        self.n.update_neutrals(n_n_s, n_n_t)

        self.izn_rate.s.update(izn_rate_s)
        self.izn_rate.t.update(izn_rate_t)
        self.izn_rate.tot.update(izn_rate_s + izn_rate_t)

        # self.izn_rate_fsa = np.array(map(lambda x: self.izn_rate_fsa[-1] * x ** 10, self.r[:, 0] / self.a))
        # self.cool_rate_fsa = np.array(map(lambda x: self.cool_rate_fsa[-1] * x ** 10, self.r[:, 0] / self.a))
        # self.dn_dr_fsa = np.array(map(lambda x: self.dn_dr_fsa[-1] * x ** 10, self.r[:, 0] / self.a))

    def num_xpts(self):
        if self.num_xpt == 1:
            return 1
        elif self.num_xpt == 2:
            return 2
        else:
            raise Exception("0 or >2 Xpoints found")
