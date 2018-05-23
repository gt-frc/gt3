#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 10:21:23 2018

@author: max
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys


class marfe:
    """
    """
    def __init__(self, inp, core, nbi, imp, ntrl):
        self.calc_marfe_denlim(inp, core, nbi, imp)
    
    def calc_marfe_denlim(self, inp, core, nbi, imp):
        """
        """
        # specified quantities
        chi_r = 2.0
        nu = 5.0/2.0
        sv_ion = core.sv_ion
        sv_cx = core.sv_cx
        sv_el = core.sv_el
        L_T = core.L_Ti_J
        L_n = core.L_ni
        E_ion = 15.466 * 1.6021E-19  # ionization energy of deuterium in Joules
        T = core.Ti_J
        fz = core.nC / core.ni  # ne?
        f0 = core.n_n_total / core.ni
        f0c = core.n_n_slow / core.ni
        Lz = imp.core_emissivity
        dLzdT = imp.core_dEmiss_dT
        dsv_ion_dT = core.dsv_ion_dT
        dsv_cxel_dT = core.dsv_cx_dT + core.dsv_el_dT
        
        # TODO: Generalize this to non carbon impurities. Refer to equations 14.39 in Stacey's book.
        Ci2 = 1.56*(1.0+np.sqrt(2)*core.z_0)*(1.0+0.52*core.z_0) / \
            ((1.0+2.65*core.z_0)*(1.0+0.285*core.z_0)*(core.z_0 + np.sqrt(0.5*(1.0+(1.0/6.0)))))
        Ce2 = 1.5*(1.0 - 0.6934/1.3167**core.z_eff)
        C2 = Ce2 - core.z_0 * Ci2
        
        t1 = chi_r * (nu * L_T**-2 - (1.0 - C2) * L_T**-1 * L_n**-1)
        t2 = fz*((nu + 1 - C2)*Lz/T - dLzdT)
        t3 = f0 * (E_ion * sv_ion / T * (nu - T / sv_ion * dsv_ion_dT))
        t4 = f0c * (3.0/2.0*(sv_cx + sv_el) * (nu-1.0-T*dsv_cxel_dT/(sv_cx + sv_el)))

        n_marfe = t1 / (t2 + t3 + t4)
        n_marfe_average = np.nanmean(n_marfe)
        print 'n_marfe_average = ', n_marfe_average
        n_marfe_edge = np.where(core.psi_norm > 0.5, n_marfe, np.nan)
        n_marfe_met = np.where((core.ni > n_marfe) & (n_marfe > 0) & (core.psi_norm > 0.5), core.ni-n_marfe, np.nan)
        
        marfe_fig1 = plt.figure(figsize=(13, 6))
        ax1 = marfe_fig1.add_subplot(1, 3, 1)
        ax1.axis('equal')
        ax1.set_title(r'$n_{MARFE}$')
        cs1 = ax1.contourf(core.R, core.Z, np.log10(n_marfe_edge), 500)
        ax1.plot(core.R[-1, :], core.Z[-1, :], lw=1, color='red')
        ax1.plot(inp.wall_exp[:, 0], inp.wall_exp[:, 1], lw=1, color='black')
        marfe_fig1.colorbar(cs1, ax=ax1)

        ax2 = marfe_fig1.add_subplot(1, 3, 2)
        ax2.axis('equal')
        ax2.set_title(r'$n_i$')
        cs2 = ax2.contourf(core.R, core.Z, np.log10(core.ni), 500)
        # ax1.plot(core.R[-1, :], core.Z[-1, :], lw=1, color='red')
        ax2.plot(inp.wall_exp[:, 0], inp.wall_exp[:, 1], lw=1, color='black')
        marfe_fig1.colorbar(cs2, ax=ax2)

        ax3 = marfe_fig1.add_subplot(1, 3, 3)
        ax3.axis('equal')
        ax3.set_title('$n_i-n_{MARFE}$ where greater')
        cs3 = ax3.contourf(core.R, core.Z, n_marfe_met, 500)
        ax3.plot(core.R[-1, :], core.Z[-1, :], lw=1, color='red')
        ax3.plot(inp.wall_exp[:, 0], inp.wall_exp[:, 1], lw=1, color='black')
        marfe_fig1.colorbar(cs3, ax=ax3)
        plt.tight_layout()
