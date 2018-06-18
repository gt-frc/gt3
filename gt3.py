#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
"""
import sys
import numpy as np
from collections import namedtuple
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from exp_neutpy_prep import Neutrals
from thermaliol import iol_calc
from read_infile import read_infile
from imp_rad import ImpRad
from exp_core_brnd import exp_core_brnd
from exp_sol_brnd import exp_sol_brnd
from exp_pfr_brnd import exp_pfr_brnd
from mil_core_brnd import mil_core_brnd
from mil_sol_brnd import mil_sol_brnd
from mil_pfr_brnd import mil_pfr_brnd
from beamdep import beamdep
from dens_lim import dens_lim
from marfe import marfe
from rad_trans import rad_trans


class gt3:
    """GT3 calculates various tokamak-related quantities

    Methods:
        allthethings
        justneutrals
        plotstuff

    Attributes:

    External Dependencies:
        Triangle:       Used to create the triangular mesh for neutpy. Source
                        code and documentation can be found at
                        https://www.cs.cmu.edu/~quake/triangle.html. Can be
                        installed from Ubuntu repositories as well.
        Neutpy:
        nbeams:
        adpack:
    """
    def __init__(self, shotlabel=None, mode='coreonly'):
        sys.dont_write_bytecode = True 
        # Create shotlabel as an attribute of plasma class
        self.shotlabel = shotlabel
        self.inp = read_infile(self.shotlabel)
        self.core = exp_core_brnd(self.inp) if self.inp.exp_inp else mil_core_brnd(self.inp)

        if mode == 'coreonly':
            pass
        elif mode == 'thermaliol':
            self.iol = iol_calc(self.inp, self.core)
        elif mode == 'fulliol':
            self.nbi = beamdep(self.inp, self.core)
            self.iol = iol_calc(self.inp, self.core)
        elif mode == 'imp':
            self.imp = ImpRad(self.inp, self.core)
        elif mode == 'ntrls':
            self.ntrl = Neutrals(self.inp, self.core)
        elif mode == 'ntrlsandiol':
            self.ntrl = Neutrals(self.inp, self.core)
            self.iol = iol_calc(self.inp, self.core)
        elif mode == 'nbi':
            self.nbi = beamdep(self.inp, self.core)
        elif mode == 'therm_instab':
            self.nbi = beamdep(self.inp, self.core)
            self.ntrl = Neutrals(self.inp, self.core)
            self.imp = ImpRad(self.inp, self.core)
            self.dl = dens_lim(self.inp, self.core, self.nbi, self.imp, self.ntrl)
            self.mar = marfe(self.inp, self.core, self.imp)
        elif mode == 'allthethings':
            self.nbi = beamdep(self.inp, self.core)
            self.iol = iol_calc(self.inp, self.core)
            self.ntrl = Neutrals(self.inp, self.core)
            self.imp = ImpRad(self.inp, self.core)
            self.dl = dens_lim(self.inp, self.core, self.nbi, self.imp, self.ntrl)
            self.mar = marfe(self.inp, self.core, self.imp)
        elif mode == 'radialtrans':
            self.rtrans = rad_trans(self.inp,self.core,self.iol,self.ntrl,self.nbi)

if __name__ == "__main__":
    #d3d_144977_3000 = gt3('144977_3000/togt3_d3d_144977_3000_mil')
    d3d_144977_3000 = gt3('144977_3000/togt3_d3d_144977_3000', mode='ntrlsandiol')
    #iter_shot = gt3('iter/togt3_iter_mil')
    # myshot.coreonly()
    #iter_shot.coreandiol()
    # myshot.therm_instab()
    # myshot.ntrlsonly()
    # myshot.coreandimp()
    # plt.axis('equal')
    # plt.contourf(myshot.core.R, myshot.core.Z, np.log10(myshot.core.n_n_total), 500)
    # plt.colorbar()
    # sys.exit()
    
    # fig1 = plt.figure(figsize=(6, 8))
    # ax1 = fig1.add_subplot(1, 1, 1)
    # ax1.axis('equal')
    # ax1.contour(myshot.core.R, myshot.core.Z, myshot.core.rho, 10)
    # ax1.plot(myshot.core.R[-1, :], myshot.core.Z[-1, :])
    # ax1.plot(myshot.inp.lim_vertex_closed[:, 0], myshot.inp.lim_vertex_closed[:, 1])
    # ax1.plot(myshot.inp.sep_exp_closed[:, 0], myshot.inp.sep_exp_closed[:, 1])

    myshot = d3d_144977_3000

    fontsize = 12
    fig2 = plt.figure(figsize=(7, 10))
    rows = 4
    cols = 3
    num = 1
    try:
        ax1 = fig2.add_subplot(rows, cols, num)
        ax1.set_title(r'$n_i$', fontsize=fontsize)
        ax1.plot(myshot.core.rho[:, 0], myshot.core.ni[:, 0], lw=2, color='black')
        num += 1
    except:
        pass
    
    try:
        ax2 = fig2.add_subplot(rows, cols, num)
        ax2.set_title(r'$n_e$', fontsize=fontsize)
        ax2.plot(myshot.core.rho[:, 0], myshot.core.ne[:, 0], lw=2, color='black')
        num += 1
    except:
        pass
        
    try:
        ax3 = fig2.add_subplot(rows, cols, num)
        ax3.set_title(r'$T_i$', fontsize=fontsize)
        ax3.plot(myshot.core.rho[:, 0], myshot.core.Ti_kev[:, 0], lw=2, color='black')
        num += 1
    except:
        pass
        
    try:
        ax4 = fig2.add_subplot(rows, cols, num)
        ax4.set_title(r'$T_e$', fontsize=fontsize)
        ax4.plot(myshot.core.rho[:, 0], myshot.core.Te_kev[:, 0], lw=2, color='black')
        num += 1
    except:
        pass
        
    try:
        ax5 = fig2.add_subplot(rows, cols, num)
        ax5.set_title(r'$j_r$', fontsize=fontsize)
        ax5.plot(myshot.core.rho[:, 0], myshot.core.j_r[:, 0], lw=2, color='black')
        num += 1
    except:
        pass
        
    try:
        ax6 = fig2.add_subplot(rows, cols, num)
        ax6.set_title(r'$E_r$', fontsize=fontsize)
        ax6.plot(myshot.core.rho[:, 0], myshot.core.E_r[:, 0], lw=2, color='black')
        num += 1
    except:
        pass
        
    try:
        ax7 = fig2.add_subplot(rows, cols, num)
        ax7.set_title(r'$fracz$', fontsize=fontsize)
        ax7.plot(myshot.core.rho[:, 0], myshot.core.fracz[:, 0], lw=2, color='black')
        num += 1
    except:
        pass
        
    try:
        ax8 = fig2.add_subplot(rows, cols, num)
        ax8.set_title(r'$v_{\theta, C}$', fontsize=fontsize)
        ax8.plot(myshot.core.rho[:, 0], myshot.core.vpolC[:, 0], lw=2, color='black')
        num += 1
    except:
        pass
        
    try:
        ax9 = fig2.add_subplot(rows, cols, num)
        ax9.set_title(r'$v_{\phi, C}$', fontsize=fontsize)
        ax9.plot(myshot.core.rho[:, 0], myshot.core.vtorC[:, 0], lw=2, color='black')
        num += 1
    except:
        pass
        
    try:
        ax10 = fig2.add_subplot(rows, cols, num)
        ax10.set_title(r'$F_{orb}$', fontsize=fontsize)
        ax10.set_xlim(0.9, 1.0)
        ax10.plot(myshot.core.rho[:, 0], myshot.iol.F_orb_1D, lw=2, color='black')
        num += 1
    except:
        pass
            
    try:
        ax11 = fig2.add_subplot(rows, cols, num)
        ax11.set_title(r'$M_{orb}$', fontsize=fontsize)
        ax11.set_xlim(0.9, 1.0)
        ax11.plot(myshot.core.rho[:, 0], myshot.iol.M_orb_1D, lw=2, color='black')
        num += 1
    except:
        pass
            
    try:
        ax12 = fig2.add_subplot(rows, cols, num)
        ax12.set_title(r'$E_{orb}$', fontsize=fontsize)
        ax12.set_xlim(0.9, 1.0)
        ax12.plot(myshot.core.rho[:, 0], myshot.iol.E_orb_1D, lw=2, color='black')
        num += 1
    except:
        pass
            
    try:
        ax13 = fig2.add_subplot(rows, cols, num)
        ax13.set_title(r'$I_{cum}$', fontsize=fontsize)
        ax13.plot(myshot.core.rho[:, 0], myshot.core.I[:, 0], lw=2, color='black')
        num += 1
    except:
        pass
    
    try:
        ax14 = fig2.add_subplot(rows, cols, num)
        ax14.set_title(r'$Z_{eff}$', fontsize=fontsize)
        ax14.plot(myshot.core.rho[:, 0], myshot.core.z_eff[:, 0], lw=2, color='black')
        num += 1
    except:
        pass
    
    try:
        ax15 = fig2.add_subplot(rows, cols, num)
        ax15.set_title(r'$NBI Dep. Prof.$', fontsize=fontsize)
        ax15.plot(myshot.core.rho[:, 0], myshot.nbi.pNB_tot[:, 0], lw=2, color='black')
        num += 1
    except:
        pass
    plt.tight_layout()


    fontsize = 12
    fig3 = plt.figure(figsize=(6, 4))
    num = 1
    try:
        ax1 = fig3.add_subplot(1, 1, 1)
        ax1.set_title(r'$\psi(\rho)$', fontsize=fontsize)
        ax1.plot(myshot.core.rho[:, 0], myshot.core.psi[:, 0], lw=2, color='black')
        num += 1
    except:
        pass


    fontsize = 12
    fiol_d3d_v_iter = plt.figure(figsize=(6, 6))
    try:
        ax1 = fiol_d3d_v_iter.add_subplot(1, 1, 1)
        ax1.set_title(r'$F_{orb}$ DIII-D vs ITER', fontsize=fontsize)
        ax1.set_xlim(0.8,1.0)
        ax1.set_xlabel(r'normalized radius ($\rho$)')
        ax1.plot(d3d_144977_3000.core.rho[:, 0], d3d_144977_3000.iol.forb_d_therm[:, 0], lw=2, color='black',label='DIII-D')
        ax1.plot(iter_shot.core.rho[:, 0], iter_shot.iol.forb_d_therm[:, 0], lw=2, color='red',label='ITER')
        ax1.legend()
    except:
       pass

    miol_d3d_v_iter = plt.figure(figsize=(6, 6))
    try:
        ax1 = miol_d3d_v_iter.add_subplot(1, 1, 1)
        ax1.set_title(r'$M_{orb}$ DIII-D vs ITER', fontsize=fontsize)
        ax1.set_xlim(0.8,1.0)
        ax1.set_xlabel(r'normalized radius ($\rho$)')
        ax1.plot(d3d_144977_3000.core.rho[:, 0], d3d_144977_3000.iol.morb_d_therm[:, 0], lw=2, color='black',label='DIII-D')
        ax1.plot(iter_shot.core.rho[:, 0], iter_shot.iol.morb_d_therm[:, 0], lw=2, color='red',label='ITER')
        ax1.legend()
    except:
        pass

    eiol_d3d_v_iter = plt.figure(figsize=(6, 6))
    try:
        ax1 = eiol_d3d_v_iter.add_subplot(1, 1, 1)
        ax1.set_title(r'$E_{orb}$ DIII-D vs ITER', fontsize=fontsize)
        ax1.set_xlim(0.8,1.0)
        ax1.set_xlabel(r'normalized radius ($\rho$)')
        ax1.plot(d3d_144977_3000.core.rho[:, 0], d3d_144977_3000.iol.eorb_d_therm[:, 0], lw=2, color='black',label='DIII-D')
        ax1.plot(iter_shot.core.rho[:, 0], iter_shot.iol.eorb_d_therm[:, 0], lw=2, color='red',label='ITER')
        ax1.legend()
    except:
        pass


    forb_species_d3d = plt.figure(figsize=(6, 6))
    try:
        ax1 = forb_species_d3d.add_subplot(1, 1, 1)
        ax1.set_title(r'$F_{orb}$ on DIII-D for various species', fontsize=fontsize)
        ax1.set_xlim(0, 1.0)
        ax1.set_xlabel(r'normalized radius ($\rho$)')
        ax1.plot(d3d_144977_3000.core.rho[:, 0], d3d_144977_3000.iol.forb_d_therm[:, 0], lw=2, label='Thermal Deuterium')
        ax1.plot(d3d_144977_3000.core.rho[:, 0], d3d_144977_3000.iol.forb_t_therm[:, 0], lw=2, label='Thermal Tritium')
        ax1.plot(d3d_144977_3000.core.rho[:, 0], d3d_144977_3000.iol.forb_c_therm[:, 0], lw=2, label='Thermal Carbon')
        ax1.plot(d3d_144977_3000.core.rho[:, 0], d3d_144977_3000.iol.forb_a_therm[:, 0], lw=2, label='Thermal Alphas')
        ax1.plot(d3d_144977_3000.core.rho[:, 0], d3d_144977_3000.iol.forb_a_fast[:, 0], lw=2, label='Fast Alphas')
        ax1.legend()
    except:
        pass

    morb_species_d3d = plt.figure(figsize=(6, 6))
    try:
        ax1 = morb_species_d3d.add_subplot(1, 1, 1)
        ax1.set_title(r'$M_{orb}$ on DIII-D for various species', fontsize=fontsize)
        ax1.set_xlim(0, 1.0)
        ax1.set_xlabel(r'normalized radius ($\rho$)')
        ax1.plot(d3d_144977_3000.core.rho[:, 0], d3d_144977_3000.iol.morb_d_therm[:, 0], lw=2, label='Thermal Deuterium')
        ax1.plot(d3d_144977_3000.core.rho[:, 0], d3d_144977_3000.iol.morb_t_therm[:, 0], lw=2, label='Thermal Tritium')
        ax1.plot(d3d_144977_3000.core.rho[:, 0], d3d_144977_3000.iol.morb_c_therm[:, 0], lw=2, label='Thermal Carbon')
        ax1.plot(d3d_144977_3000.core.rho[:, 0], d3d_144977_3000.iol.morb_a_therm[:, 0], lw=2, label='Thermal Alphas')
        ax1.plot(d3d_144977_3000.core.rho[:, 0], d3d_144977_3000.iol.morb_a_fast[:, 0], lw=2, label='Fast Alphas')
        ax1.legend()
    except:
        pass

    eorb_species_d3d = plt.figure(figsize=(6, 6))
    try:
        ax1 = eorb_species_d3d.add_subplot(1, 1, 1)
        ax1.set_title(r'$E_{orb}$ on DIII-D for various species', fontsize=fontsize)
        ax1.set_xlim(0, 1.0)
        ax1.set_xlabel(r'normalized radius ($\rho$)')
        ax1.plot(d3d_144977_3000.core.rho[:, 0], d3d_144977_3000.iol.eorb_d_therm[:, 0], lw=2, label='Thermal Deuterium')
        ax1.plot(d3d_144977_3000.core.rho[:, 0], d3d_144977_3000.iol.eorb_t_therm[:, 0], lw=2, label='Thermal Tritium')
        ax1.plot(d3d_144977_3000.core.rho[:, 0], d3d_144977_3000.iol.eorb_c_therm[:, 0], lw=2, label='Thermal Carbon')
        ax1.plot(d3d_144977_3000.core.rho[:, 0], d3d_144977_3000.iol.eorb_a_therm[:, 0], lw=2, label='Thermal Alphas')
        ax1.plot(d3d_144977_3000.core.rho[:, 0], d3d_144977_3000.iol.eorb_a_fast[:, 0], lw=2, label='Fast Alphas')
        ax1.legend()
    except:
        pass

    forb_species_iter = plt.figure(figsize=(6, 6))
    try:
        ax1 = forb_species_iter.add_subplot(1, 1, 1)
        ax1.set_title(r'$F_{orb}$ on ITER for various species', fontsize=fontsize)
        ax1.set_xlim(0, 1.0)
        ax1.set_xlabel(r'normalized radius ($\rho$)')
        ax1.plot(iter_shot.core.rho[:, 0], iter_shot.iol.forb_d_therm[:, 0], lw=2, label='Thermal Deuterium')
        ax1.plot(iter_shot.core.rho[:, 0], iter_shot.iol.forb_t_therm[:, 0], lw=2, label='Thermal Tritium')
        ax1.plot(iter_shot.core.rho[:, 0], iter_shot.iol.forb_c_therm[:, 0], lw=2, label='Thermal Carbon')
        ax1.plot(iter_shot.core.rho[:, 0], iter_shot.iol.forb_a_therm[:, 0], lw=2, label='Thermal Alphas')
        ax1.plot(iter_shot.core.rho[:, 0], iter_shot.iol.forb_a_fast[:, 0], lw=2, label='Fast Alphas')
        ax1.legend()
    except:
        pass

    morb_species_iter = plt.figure(figsize=(6, 6))
    try:
        ax1 = morb_species_iter.add_subplot(1, 1, 1)
        ax1.set_title(r'$M_{orb}$ on ITER for various species', fontsize=fontsize)
        ax1.set_xlim(0, 1.0)
        ax1.set_xlabel(r'normalized radius ($\rho$)')
        ax1.plot(iter_shot.core.rho[:, 0], iter_shot.iol.morb_d_therm[:, 0], lw=2, label='Thermal Deuterium')
        ax1.plot(iter_shot.core.rho[:, 0], iter_shot.iol.morb_t_therm[:, 0], lw=2, label='Thermal Tritium')
        ax1.plot(iter_shot.core.rho[:, 0], iter_shot.iol.morb_c_therm[:, 0], lw=2, label='Thermal Carbon')
        ax1.plot(iter_shot.core.rho[:, 0], iter_shot.iol.morb_a_therm[:, 0], lw=2, label='Thermal Alphas')
        ax1.plot(iter_shot.core.rho[:, 0], iter_shot.iol.morb_a_fast[:, 0], lw=2, label='Fast Alphas')
        ax1.legend()
    except:
        pass

    eorb_species_iter = plt.figure(figsize=(6, 6))
    try:
        ax1 = eorb_species_iter.add_subplot(1, 1, 1)
        ax1.set_title(r'$E_{orb}$ on ITER for various species', fontsize=fontsize)
        ax1.set_xlim(0, 1.0)
        ax1.set_xlabel(r'normalized radius ($\rho$)')
        ax1.plot(iter_shot.core.rho[:, 0], iter_shot.iol.eorb_d_therm[:, 0], lw=2, label='Thermal Deuterium')
        ax1.plot(iter_shot.core.rho[:, 0], iter_shot.iol.eorb_t_therm[:, 0], lw=2, label='Thermal Tritium')
        ax1.plot(iter_shot.core.rho[:, 0], iter_shot.iol.eorb_c_therm[:, 0], lw=2, label='Thermal Carbon')
        ax1.plot(iter_shot.core.rho[:, 0], iter_shot.iol.eorb_a_therm[:, 0], lw=2, label='Thermal Alphas')
        ax1.plot(iter_shot.core.rho[:, 0], iter_shot.iol.eorb_a_fast[:, 0], lw=2, label='Fast Alphas')
        ax1.legend()
    except:
        pass

    #plt.xlim(0.8,1.0)
    # plt.plot(brnd.rho[:,0],self.fiol_d_therm[:,0],label='Deuterium')
    # plt.plot(brnd.rho[:,0],self.fiol_t_therm[:,0],label='Tritium')
    # plt.plot(brnd.rho[:,0],self.fiol_c_therm[:,0],label='Carbon')
    # plt.plot(brnd.rho[:,0],self.fiol_a_therm[:,0],label='Thermal Alphas')
    # plt.plot(brnd.rho[:,0],self.fiol_a_fast[:,0],label='Fast Alphas')
    # plt.plot(brnd.rho[:,0],self.miol_d_therm[:,0],label='Deuterium')
    # plt.plot(brnd.rho[:,0],self.miol_t_therm[:,0],label='Tritium')
    # plt.plot(brnd.rho[:,0],self.miol_c_therm[:,0],label='Carbon')
    # plt.plot(brnd.rho[:,0],self.miol_a_therm[:,0],label='Thermal Alphas')
    # plt.plot(brnd.rho[:,0],self.miol_a_fast[:,0],label='Fast Alphas')
    # plt.legend()
    # plt.show()
    # sys.exit()
    plt.show()
