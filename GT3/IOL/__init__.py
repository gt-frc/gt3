# -*- coding: utf-8 -*-
"""
Created on Sat Aug  5 19:48:34 2017

@author: max
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from math import pi, sqrt
from collections import namedtuple
from Functions.CalcIOLMaxwellian import calc_iol_maxwellian
from Functions.CalcIOLMonoEn import calc_iol_mono_en
from Functions.CalcIOLBeams import calc_iol_beams
import sys

m_d = 3.343583719e-27
m_t = 5.006e-27
m_c = 1.9926467e-26
m_a = 6.643e-27


class IOL:
    """
    """
    def __init__(self, inp, core):
        sys.dont_write_bytecode = True
        np.warnings.filterwarnings('ignore')

        numcos = inp.numcos  # TODO: Clean this up
        angles = np.linspace(-1, 1, inp.numcos + 1)
        self.coslist = ((angles + np.roll(angles, -1)) / 2)[:-1]
        self.rho = core.rho
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
                                 core.B_tot)[-1]

        f0 = np.broadcast_arrays(np.ones(polpts)[:, None, None, None],
                                 np.ones(inp.numcos)[:, None, None],
                                 core.f_phi)[-1]

        psi0 = np.broadcast_arrays(np.ones(polpts)[:, None, None, None],
                                   np.ones(inp.numcos)[:, None, None],
                                   core.psi)[-1]

        phi0 = np.broadcast_arrays(np.ones(polpts)[:, None, None, None],
                                   np.ones(inp.numcos)[:, None, None],
                                   core.E_pot)[-1]# * 1E3  # now in volts

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
                                 core.B_tot[-1][:, None, None, None])[-1]

        psi1 = np.broadcast_arrays(np.ones(polpts)[:, None, None, None],
                                   np.ones(inp.numcos)[:, None, None],
                                   np.ones(radpts)[:, None],
                                   np.ones(polpts)[:],
                                   core.psi[-1][:, None, None, None])[-1]

        phi1 = np.broadcast_arrays(np.ones(polpts)[:, None, None, None],
                                   np.ones(inp.numcos)[:, None, None],
                                   np.ones(radpts)[:, None],
                                   np.ones(polpts)[:],
                                   core.E_pot[-1][:, None, None, None])[-1] #* 1E3  # now in volts

        Tprofile = namedtuple('Tprofile', 'i C')(
            core.T.i.kev.T[0],
            core.T.C.kev.T[0]
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
        v_alpha = sqrt(2*3.5E6*1.6021E-19/m_a)
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
        v_beam = sqrt(2*80.0E3*1.6021E-19/m_d)
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

    def plot_F_i(self):
        """
        Plots the 1D GT3.IOL differential ion number loss fraction
        """

        plot = plt.figure()
        fig1 = plot.add_subplot(111)
        fig1.set_xlabel(r'$\rho$', fontsize=20)
        fig1.set_ylabel(r'$\frac{\partial F}{\partial r}$', fontsize=25)
        fig1.set_title('GT3.IOL differential ion number loss fraction')
        fig1.scatter(self.rho[:, 0], self.forb_d_therm_1D, marker='o', color='blue')

    def plot_M_i(self):
        """
        Plots the 1D GT3.IOL differential ion momentum loss fraction
        """

        plot = plt.figure()
        fig1 = plot.add_subplot(111)
        fig1.set_xlabel(r'$\rho$', fontsize=20)
        fig1.set_ylabel(r'$\frac{\partial M}{\partial r}$', fontsize=25)
        fig1.set_title('GT3.IOL differential ion momentum loss fraction')
        fig1.scatter(self.rho[:, 0], self.morb_d_therm_1D, marker='o', color='green')

    def plot_E_i(self):
        """
        Plots the 1D GT3.IOL differential ion energy loss fraction
        """

        plot = plt.figure()
        fig1 = plot.add_subplot(111)
        fig1.set_xlabel(r'$\rho$', fontsize=20)
        fig1.set_ylabel(r'$\frac{\partial E}{\partial r}$', fontsize=25)
        fig1.set_title('GT3.IOL differential ion energy loss fraction')
        fig1.scatter(self.rho[:, 0], self.eorb_d_therm_1D, marker='o', color='red')

    def plot_all_i(self):

        plot = plt.figure()
        fig1 = plot.add_subplot(131)
        fig1.set_xlabel(r'$\rho$', fontsize=20)
        fig1.set_ylabel(r'$\frac{\partial F}{\partial r}$', fontsize=25)
        fig1.set_title('Ion number loss fraction')
        fig1.scatter(self.rho[:, 0], self.forb_d_therm_1D, marker='o', color='blue')

        fig2 = plot.add_subplot(132)
        fig2.set_xlabel(r'$\rho$', fontsize=20)
        fig2.set_ylabel(r'$\frac{\partial E}{\partial r}$', fontsize=25)
        fig2.set_title('Ion energy loss fraction')
        fig2.scatter(self.rho[:, 0], self.eorb_d_therm_1D, marker='o', color='red')

        fig3 = plot.add_subplot(133)
        fig3.set_xlabel(r'$\rho$', fontsize=20)
        fig3.set_ylabel(r'$\frac{\partial M}{\partial r}$', fontsize=25)
        fig3.set_title('Ion momentum loss fraction')
        fig3.scatter(self.rho[:, 0], self.morb_d_therm_1D, marker='o', color='green')

    def plot_F_C(self):
        """
        Plots the 1D GT3.IOL differential carbon number loss fraction
        """

        plot = plt.figure()
        fig1 = plot.add_subplot(111)
        fig1.set_xlabel(r'$\rho$', fontsize=20)
        fig1.set_ylabel(r'$\frac{\partial F}{\partial r}$', fontsize=25)
        fig1.set_title('GT3.IOL differential carbon number loss fraction')
        fig1.scatter(self.rho[:, 0], self.forb_c_therm_1D, marker='o', color='blue')

    def plot_M_C(self):
        """
        Plots the 1D GT3.IOL differential carbon momentum loss fraction
        """

        plot = plt.figure()
        fig1 = plot.add_subplot(111)
        fig1.set_xlabel(r'$\rho$', fontsize=20)
        fig1.set_ylabel(r'$\frac{\partial M}{\partial r}$', fontsize=25)
        fig1.set_title('GT3.IOL differential carbon momentum loss fraction')
        fig1.scatter(self.rho[:, 0], self.morb_c_therm_1D, marker='o', color='green')

    def plot_E_C(self):
        """
        Plots the 1D GT3.IOL differential carbon energy loss fraction
        """

        plot = plt.figure()
        fig1 = plot.add_subplot(111)
        fig1.set_xlabel(r'$\rho$', fontsize=20)
        fig1.set_ylabel(r'$\frac{\partial E}{\partial r}$', fontsize=25)
        fig1.set_title('GT3.IOL differential carbon energy loss fraction')
        fig1.scatter(self.rho[:, 0], self.eorb_c_therm_1D, marker='o', color='red')

    def plot_all_C(self):

        plot = plt.figure()
        fig1 = plot.add_subplot(111)
        fig1.set_xlabel(r'$\rho$', fontsize=20)
        fig1.set_ylabel(r'$\frac{\partial F}{\partial r}$', fontsize=25)
        fig1.set_title('Carbon number loss fraction')
        fig1.scatter(self.rho[:, 0], self.forb_c_therm_1D, marker='o', color='blue')

        fig2 = plot.add_subplot(121)
        fig2.set_xlabel(r'$\rho$', fontsize=20)
        fig2.set_ylabel(r'$\frac{\partial E}{\partial r}$', fontsize=25)
        fig2.set_title('Carbon energy loss fraction')
        fig2.scatter(self.rho[:, 0], self.eorb_c_therm_1D, marker='o', color='red')

        fig3 = plot.add_subplot(131)
        fig3.set_xlabel(r'$\rho$', fontsize=20)
        fig3.set_ylabel(r'$\frac{\partial M}{\partial r}$', fontsize=25)
        fig3.set_title('Carbon momentum loss fraction')
        fig3.scatter(self.rho[:, 0], self.morb_c_therm_1D, marker='o', color='green')

    def plot_F_i_fast(self):
        """
        Plots the 1D GT3.IOL differential fast ion number loss fraction
        """

        plot = plt.figure()
        fig1 = plot.add_subplot(111)
        fig1.set_xlabel(r'$\rho$', fontsize=20)
        fig1.set_ylabel(r'$\frac{\partial F}{\partial r}$', fontsize=25)
        fig1.set_title('GT3.IOL differential fast ion number loss fraction')
        fig1.scatter(self.rho[:, 0], self.forb_d_nbi_1D, marker='o', color='blue')

    def plot_M_i_fast(self):
        """
        Plots the 1D GT3.IOL differential fast ion momentum loss fraction
        """

        plot = plt.figure()
        fig1 = plot.add_subplot(111)
        fig1.set_xlabel(r'$\rho$', fontsize=20)
        fig1.set_ylabel(r'$\frac{\partial M}{\partial r}$', fontsize=25)
        fig1.set_title('GT3.IOL differential fast ion momentum loss fraction')
        fig1.scatter(self.rho[:, 0], self.morb_d_nbi_1D, marker='o', color='green')

    def plot_E_i_fast(self):
        """
        Plots the 1D GT3.IOL differential fast ion energy loss fraction
        """

        plot = plt.figure()
        fig1 = plot.add_subplot(111)
        fig1.set_xlabel(r'$\rho$', fontsize=20)
        fig1.set_ylabel(r'$\frac{\partial E}{\partial r}$', fontsize=25)
        fig1.set_title('GT3.IOL differential fast ion energy loss fraction')
        fig1.scatter(self.rho[:, 0], self.eorb_d_nbi_1D, marker='o', color='red')

    def plot_all_i_fast(self):

        plot = plt.figure()
        fig1 = plot.add_subplot(111)
        fig1.set_xlabel(r'$\rho$', fontsize=20)
        fig1.set_ylabel(r'$\frac{\partial F}{\partial r}$', fontsize=25)
        fig1.set_title('Fast ion number loss fraction')
        fig1.scatter(self.rho[:, 0], self.forb_d_nbi_1D, marker='o', color='blue')

        fig2 = plot.add_subplot(121)
        fig2.set_xlabel(r'$\rho$', fontsize=20)
        fig2.set_ylabel(r'$\frac{\partial E}{\partial r}$', fontsize=25)
        fig2.set_title('Fast ion energy loss fraction')
        fig2.scatter(self.rho[:, 0], self.eorb_d_nbi_1D, marker='o', color='red')

        fig3 = plot.add_subplot(131)
        fig3.set_xlabel(r'$\rho$', fontsize=20)
        fig3.set_ylabel(r'$\frac{\partial M}{\partial r}$', fontsize=25)
        fig3.set_title('Fast ion momentum loss fraction')
        fig3.scatter(self.rho[:, 0], self.morb_d_nbi_1D, marker='o', color='green')


# iolplot=1
# if iolplot==1:
#    fig = plt.figure(figsize=(5, 5))
#    fig.suptitle('IOL a, b, c in DIII-D with cos:{}'.format(coslist[-2]), fontsize=15)
#    ax1 = fig.add_subplot(221)
#    ax1.set_title(r'$v_{sep-1}$', fontsize=12)
#    ax1.set_ylabel(r'R', fontsize=12)
#    ax1.set_xlabel(r'Z', fontsize=12)
#    ax1.axis('equal')
#    ax1.grid(b=True, which='both', axis='both')
#    CS = ax1.contourf(brnd.R, brnd.Z, v_sep_1[:, -2, :, 0].T, 500)
#    plt.colorbar(CS)
#
#    ax2 = fig.add_subplot(222)
#    ax2.set_title(r'$v_{sep-2}$', fontsize=12)
#    ax2.set_ylabel(r'R', fontsize=12)
#    ax2.set_xlabel(r'Z', fontsize=12)
#    ax2.axis('equal')
#    ax2.grid(b=True, which='both', axis='both')
#    CS = ax2.contourf(brnd.R, brnd.Z, v_sep_2[:, -2, :, 0].T, 500)
#    plt.colorbar(CS)
#
#    ax3 = fig.add_subplot(223)
#    ax3.set_title(r'$v_{sep}$', fontsize=12)
#    ax3.set_ylabel(r'R', fontsize=12)
#    ax3.set_xlabel(r'Z', fontsize=12)
#    ax3.axis('equal')
#    ax3.grid(b=True, which='both', axis='both')
#    #CS = ax3.pcolor(brnd.R, brnd.Z, v_sep[:, 10, :, 0].T, vmin=0, vmax=1E8)
#    CS = ax3.contourf(brnd.R, brnd.Z, v_sep[:, -2, :, 0].T, 500)
#    plt.colorbar(CS)

#    plt.tight_layout()
#    fig.subplots_adjust(top=0.84)
# if iolplot==1:
#     fig = plt.figure(figsize=(6, 8))
#     fig.suptitle(r'IOL v_{} with cos:{}'.format('esc', coslist[0]), fontsize=15)
#     ax1 = fig.add_subplot(111)
#     #ax1.set_title(r'$v_{sep}$', fontsize=12)
#     ax1.set_ylabel(r'Z', fontsize=12)
#     ax1.set_xlabel(r'R', fontsize=12)
#     ax1.axis('equal')
#     ax1.grid(b=True, which='both', axis='both')
#
#     test = np.nanmin(v_sep, axis=3).T
#
#     CS = ax1.contourf(brnd.R, brnd.Z, np.log10((0.5*m*test[:, 0, :]**2)/(1.6021E-19 * 1.0E3)), 500)
#     plt.colorbar(CS)
#     #plt.tight_layout()
#     fig.subplots_adjust(top=0.95)
#     #fig.subplots_adjust(top=0.95)
# return
