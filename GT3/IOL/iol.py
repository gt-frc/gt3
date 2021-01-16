# -*- coding: utf-8 -*-
"""
Created on Sat Aug  5 19:48:34 2017

@author: max
"""

import numpy as np
import matplotlib.pyplot as plt
from math import pi, sqrt
from collections import namedtuple
from .Functions.CalcIOLMaxwellian import calc_iol_maxwellian
from .Functions.CalcIOLMonoEn import calc_iol_mono_en
from .Functions.CalcIOLBeams import calc_iol_beams
from GT3.utilities.PlotBase import PlotBase
import sys
import GT3.constants as constants
from GT3 import Core

m_d = constants.deuteron_mass
m_t = constants.triton_mass
m_c = constants.carbon_mass
m_a = constants.alpha_mass


class IOL(PlotBase):
    """
    """

    def __init__(self, inp, core: Core):
        super().__init__()
        sys.dont_write_bytecode = True
        np.warnings.filterwarnings('ignore')

        numcos = inp.numcos  # TODO: Clean this up
        angles = np.linspace(-1, 1, inp.numcos + 1)
        self.coslist = ((angles + np.roll(angles, -1)) / 2)[:-1]
        self.rho = core.rho
        self.set_plot_rho1d(self.rho[:, 0])
        self.calc_iol_beams = calc_iol_beams
        polpts = len(core.rho[-1])
        radpts = len(core.rho.T[-1])

        # THE FOLLOWING ARRAYS ARE 4-DIMENSIONAL ARRAYS
        # [ LAUNCH THETA POSITION , LAUNCH ANGLE COSINE, LAUNCH r  , EXIT THETA POSITION  ]

        # NOTE TO FUTURE DEVELOPERS: IF YOU TRY TO LOOP OVER THE PLASMA POINTS, LAUNCH ANGLES, AND
        # EXIT LOCATIONS IN PYTHON THE WAY YOU WOULD IN C OR FORTRAN, IT'S GOING TO TAKE FOREVER.
        # ALTHOUGH THESE ARRAYS TAKE MORE MEMORY THAT I'D LIKE, IT'S CURRENTLY NECESSARY TO DO IT THIS WAY.
        # MAYBE SOMETHING TO IMPROVE ON IN THE FUTURE.

        # CREATE ARRAYS FOR LAUNCH POINTS IN THE PLASMA
        r0 = np.broadcast_arrays(np.ones(polpts)[:, None, None, None],
                                 np.ones(inp.numcos)[:, None, None],
                                 core.rho)[-1]

        B0 = np.broadcast_arrays(np.ones(polpts)[:, None, None, None],
                                 np.ones(inp.numcos)[:, None, None],
                                 core.B.tot.val)[-1]

        f0 = np.broadcast_arrays(np.ones(polpts)[:, None, None, None],
                                 np.ones(inp.numcos)[:, None, None],
                                 core.f_phi)[-1]

        psi0 = np.broadcast_arrays(np.ones(polpts)[:, None, None, None],
                                   np.ones(inp.numcos)[:, None, None],
                                   core.psi.psi)[-1]

        phi0 = np.broadcast_arrays(np.ones(polpts)[:, None, None, None],
                                   np.ones(inp.numcos)[:, None, None],
                                   core.E_pot.val)[-1]  # * 1E3  # now in volts

        zeta0 = np.broadcast_arrays(np.ones(polpts)[:, None, None, None],
                                    self.coslist[:, None, None],
                                    np.ones(core.R.shape))[1]

        # CREATE ARRAYS FOR DESTINATION POINTS ALONG THE SEPERATRIX
        R1 = np.broadcast_arrays(np.ones(polpts)[:, None, None, None],
                                 np.ones(inp.numcos)[:, None, None],
                                 np.ones(radpts)[:, None],
                                 np.ones(polpts)[:],
                                 core.R[-1][:, None, None, None])[-1]

        f1 = np.broadcast_arrays(np.ones(polpts)[:, None, None, None],
                                 np.ones(inp.numcos)[:, None, None],
                                 np.ones(radpts)[:, None], np.ones(polpts)[:],
                                 core.f_phi[-1][:, None, None, None])[-1]

        B1 = np.broadcast_arrays(np.ones(polpts)[:, None, None, None],
                                 np.ones(inp.numcos)[:, None, None],
                                 np.ones(radpts)[:, None],
                                 np.ones(polpts)[:],
                                 core.B.tot.val[-1][:, None, None, None])[-1]

        psi1 = np.broadcast_arrays(np.ones(polpts)[:, None, None, None],
                                   np.ones(inp.numcos)[:, None, None],
                                   np.ones(radpts)[:, None],
                                   np.ones(polpts)[:],
                                   core.psi.psi[-1][:, None, None, None])[-1]

        phi1 = np.broadcast_arrays(np.ones(polpts)[:, None, None, None],
                                   np.ones(inp.numcos)[:, None, None],
                                   np.ones(radpts)[:, None],
                                   np.ones(polpts)[:],
                                   core.E_pot.val[-1][:, None, None, None])[-1]  # * 1E3  # now in volts

        Tprofile = namedtuple('Tprofile', 'i C')(
            core.T.i.kev.val.T[0],
            core.T.C.kev.val.T[0]
        )

        # iol_params = {}
        # iol_params['r0'] = r0
        # iol_params['B0'] = B0
        # iol_params['f0'] = f0
        # iol_params['psi0'] = psi0
        # iol_params['phi0'] = phi0
        # iol_params['zeta0'] = zeta0
        # iol_params['R1'] = R1
        # iol_params['f1'] = f1
        # iol_params['B1'] = B1
        # iol_params['psi1'] = psi1
        # iol_params['phi1'] = phi1
        # iol_params['Tprofile'] = Tprofile

        # convert iol_params to namedtuple so individual parameters can be accessed as normal attributes
        self.iol_p = namedtuple('iol_p', 'r0 B0 f0 psi0 phi0 zeta0 R1 f1 B1 psi1 phi1')(
            r0,
            B0,
            f0,
            psi0,
            phi0,
            zeta0,
            R1,
            f1,
            B1,
            psi1,
            phi1
        )

        # Calculate IOL for thermal deuterium
        forb_d_therm, morb_d_therm, eorb_d_therm = calc_iol_maxwellian(1,
                                                                       m_d,
                                                                       self.iol_p,
                                                                       core.thetapts,
                                                                       Tprofile.i,
                                                                       self.coslist,
                                                                       numcos)
        self.forb_d_therm = inp.R_loss * forb_d_therm
        self.morb_d_therm = inp.R_loss * morb_d_therm
        self.eorb_d_therm = inp.R_loss * eorb_d_therm
        self.forb_d_therm_1D = self.forb_d_therm[:, 0]
        self.morb_d_therm_1D = self.morb_d_therm[:, 0]
        self.eorb_d_therm_1D = self.eorb_d_therm[:, 0]

        # Calculate IOL for thermal tritium
        forb_t_therm, morb_t_therm, eorb_t_therm = calc_iol_maxwellian(1,
                                                                       m_t,
                                                                       self.iol_p,
                                                                       core.thetapts,
                                                                       Tprofile.i,
                                                                       self.coslist,
                                                                       numcos)
        self.forb_t_therm = inp.R_loss * forb_t_therm
        self.morb_t_therm = inp.R_loss * morb_t_therm
        self.eorb_t_therm = inp.R_loss * eorb_t_therm
        self.forb_t_therm_1D = self.forb_t_therm[:, 0]
        self.morb_t_therm_1D = self.morb_t_therm[:, 0]
        self.eorb_t_therm_1D = self.eorb_t_therm[:, 0]

        # Calculate IOL for thermal carbon
        forb_c_therm, morb_c_therm, eorb_c_therm = calc_iol_maxwellian(6,
                                                                       m_c,
                                                                       self.iol_p,
                                                                       core.thetapts,
                                                                       Tprofile.C,
                                                                       self.coslist,
                                                                       numcos)
        self.forb_c_therm = inp.R_loss * forb_c_therm
        self.morb_c_therm = inp.R_loss * morb_c_therm
        self.eorb_c_therm = inp.R_loss * eorb_c_therm
        self.forb_c_therm_1D = self.forb_c_therm[:, 0]
        self.morb_c_therm_1D = self.morb_c_therm[:, 0]
        self.eorb_c_therm_1D = self.eorb_c_therm[:, 0]

        # Calculate IOL for thermal alphas
        forb_a_therm, morb_a_therm, eorb_a_therm = calc_iol_maxwellian(2,
                                                                       m_a,
                                                                       self.iol_p,
                                                                       core.thetapts,
                                                                       Tprofile.i,
                                                                       self.coslist,
                                                                       numcos)
        self.forb_a_therm = inp.R_loss * forb_a_therm
        self.morb_a_therm = inp.R_loss * morb_a_therm
        self.eorb_a_therm = inp.R_loss * eorb_a_therm
        self.forb_a_therm_1D = self.forb_a_therm[:, 0]
        self.morb_a_therm_1D = self.morb_a_therm[:, 0]
        self.eorb_a_therm_1D = self.eorb_a_therm[:, 0]

        # Calculate IOL for fast, monoenergetic alphas
        v_alpha = sqrt(2 * 3.5E6 * 1.6021E-19 / m_a)
        forb_a_fast, morb_a_fast, eorb_a_fast = calc_iol_mono_en(2,
                                                                 m_a,
                                                                 self.iol_p,
                                                                 core.thetapts,
                                                                 v_alpha,
                                                                 self.coslist,
                                                                 numcos)

        # currently applying R_loss to fast alphas as well as thermal, although I'm skeptical of this. -MH
        self.forb_a_fast = inp.R_loss * forb_a_fast
        self.morb_a_fast = inp.R_loss * morb_a_fast
        self.eorb_a_fast = inp.R_loss * eorb_a_fast
        self.forb_a_fast_1D = self.forb_a_fast[:, 0]
        self.morb_a_fast_1D = self.morb_a_fast[:, 0]
        self.eorb_a_fast_1D = self.eorb_a_fast[:, 0]

        # Calculate IOL for neutral deuterium beams
        v_beam = sqrt(2 * 80.0E3 * 1.6021E-19 / m_d)
        zeta_beam = -0.96
        forb_d_nbi, morb_d_nbi, eorb_d_nbi = self.calc_iol_beams(1,
                                                                 m_d,
                                                                 self.iol_p,
                                                                 core.thetapts,
                                                                 v_beam,
                                                                 zeta_beam,
                                                                 self.coslist)
        # currently applying R_loss to fast alphas as well as thermal, although I'm skeptical of this. -MH
        self.forb_d_nbi = inp.R_loss * forb_d_nbi
        self.morb_d_nbi = inp.R_loss * morb_d_nbi
        self.eorb_d_nbi = inp.R_loss * eorb_d_nbi
        self.forb_d_nbi_1D = self.forb_d_nbi[:, 0]
        self.morb_d_nbi_1D = self.morb_d_nbi[:, 0]
        self.eorb_d_nbi_1D = self.eorb_d_nbi[:, 0]

    def plot_F_i(self, edge=True):
        """
        Plots the 1D GT3.IOL differential ion number loss fraction
        """
        fig = self._plot_base(self.forb_d_therm_1D, yLabel=r'$F(\rho)$',
                              title="GT3.IOL differential ion number loss fraction", edge=edge, color='blue')
        return fig

    def plot_M_i(self, edge=True):
        """
        Plots the 1D GT3.IOL differential ion momentum loss fraction
        """
        fig = self._plot_base(self.morb_d_therm_1D, yLabel=r'$M(\rho)$',
                              title="GT3.IOL differential ion momentum loss fraction", edge=edge, color='green')
        return fig

    def plot_E_i(self, edge=True):
        """
        Plots the 1D GT3.IOL differential ion energy loss fraction
        """
        fig = self._plot_base(self.eorb_d_therm_1D, yLabel=r'$E(\rho)$',
                              title="GT3.IOL differential ion energy loss fraction", edge=edge, color='red')
        return fig

    def plot_all_i(self):
        plot = plt.figure()
        fig1 = plot.add_subplot(111)
        fig1.set_xlabel(r'$\rho$', fontsize=20)
        fig1.set_ylabel(r'$\frac{\partial F}{\partial r}$', fontsize=25)
        fig1.set_title('Ion number loss fraction')
        fig1.scatter(self.rho, self.forb_d_therm_1D, marker='o', color='blue')

        fig2 = plot.add_subplot(121)
        fig2.set_xlabel(r'$\rho$', fontsize=20)
        fig2.set_ylabel(r'$\frac{\partial E}{\partial r}$', fontsize=25)
        fig2.set_title('Ion energy loss fraction')
        fig2.scatter(self.rho, self.eorb_d_therm_1D, marker='o', color='red')

        fig3 = plot.add_subplot(131)
        fig3.set_xlabel(r'$\rho$', fontsize=20)
        fig3.set_ylabel(r'$\frac{\partial M}{\partial r}$', fontsize=25)
        fig3.set_title('Ion momentum loss fraction')
        fig3.scatter(self.rho, self.morb_d_therm_1D, marker='o', color='green')

        return plot

    def plot_F_C(self, edge=True):
        """
        Plots the 1D GT3.IOL differential carbon number loss fraction
        """
        fig = self._plot_base(self.forb_c_therm_1D, yLabel=r'$\frac{\partial F}{\partial r}$',
                              title="GT3.IOL differential carbon number loss fraction", edge=edge, color='blue')
        return fig

    def plot_M_C(self, edge=True):
        """
        Plots the 1D GT3.IOL differential carbon momentum loss fraction
        """
        fig = self._plot_base(self.morb_c_therm_1D, yLabel=r'$\frac{\partial M}{\partial r}$',
                              title="GT3.IOL differential carbon momentum loss fraction", edge=edge, color='green')
        return fig

    def plot_E_C(self, edge=True):
        """
        Plots the 1D GT3.IOL differential carbon energy loss fraction
        """
        fig = self._plot_base(self.eorb_c_therm_1D, yLabel=r'$\frac{\partial E}{\partial r}$',
                              title="GT3.IOL differential carbon energy loss fraction", edge=edge, color='green')
        return fig


    def plot_all_C(self):
        plot = plt.figure()
        fig1 = plot.add_subplot(111)
        fig1.set_xlabel(r'$\rho$', fontsize=20)
        fig1.set_ylabel(r'$\frac{\partial F}{\partial r}$', fontsize=25)
        fig1.set_title('Carbon number loss fraction')
        fig1.scatter(self.rho, self.forb_c_therm_1D, marker='o', color='blue')

        fig2 = plot.add_subplot(121)
        fig2.set_xlabel(r'$\rho$', fontsize=20)
        fig2.set_ylabel(r'$\frac{\partial E}{\partial r}$', fontsize=25)
        fig2.set_title('Carbon energy loss fraction')
        fig2.scatter(self.rho, self.eorb_c_therm_1D, marker='o', color='red')

        fig3 = plot.add_subplot(131)
        fig3.set_xlabel(r'$\rho$', fontsize=20)
        fig3.set_ylabel(r'$\frac{\partial M}{\partial r}$', fontsize=25)
        fig3.set_title('Carbon momentum loss fraction')
        fig3.scatter(self.rho, self.morb_c_therm_1D, marker='o', color='green')

        return plot

    def plot_F_i_fast(self, edge=True):
        """
        Plots the 1D GT3.IOL differential fast ion number loss fraction
        """
        fig = self._plot_base(self.forb_d_nbi_1D, yLabel=r'$\frac{\partial F}{\partial r}$',
                              title="GT3.IOL differential fast ion number loss fraction", edge=edge, color='blue')
        return fig

    def plot_M_i_fast(self,edge=True):
        """
        Plots the 1D GT3.IOL differential fast ion momentum loss fraction
        """
        fig = self._plot_base(self.morb_d_nbi_1D, yLabel=r'$\frac{\partial M}{\partial r}$',
                              title="GT3.IOL differential fast ion momentum loss fraction", edge=edge, color='green')
        return fig

    def plot_E_i_fast(self, edge=True):
        """
        Plots the 1D GT3.IOL differential fast ion energy loss fraction
        """
        fig = self._plot_base(self.eorb_d_nbi_1D, yLabel=r'$\frac{\partial E}{\partial r}$',
                              title="GT3.IOL differential fast ion energy loss fraction", edge=edge, color='red')
        return fig


    def plot_all_i_fast(self):
        plot = plt.figure()
        fig1 = plot.add_subplot(111)
        fig1.set_xlabel(r'$\rho$', fontsize=20)
        fig1.set_ylabel(r'$\frac{\partial F}{\partial r}$', fontsize=25)
        fig1.set_title('Fast ion number loss fraction')
        fig1.scatter(self.rho, self.forb_d_nbi_1D, marker='o', color='blue')

        fig2 = plot.add_subplot(121)
        fig2.set_xlabel(r'$\rho$', fontsize=20)
        fig2.set_ylabel(r'$\frac{\partial E}{\partial r}$', fontsize=25)
        fig2.set_title('Fast ion energy loss fraction')
        fig2.scatter(self.rho, self.eorb_d_nbi_1D, marker='o', color='red')

        fig3 = plot.add_subplot(131)
        fig3.set_xlabel(r'$\rho$', fontsize=20)
        fig3.set_ylabel(r'$\frac{\partial M}{\partial r}$', fontsize=25)
        fig3.set_title('Fast ion momentum loss fraction')
        fig3.scatter(self.rho, self.morb_d_nbi_1D, marker='o', color='green')

        return plot