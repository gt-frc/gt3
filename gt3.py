#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
"""


class gt3():
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
    def __init__(self, shotlabel=None):
        sys.dont_write_bytecode = True 
        # Create shotlabel as an attribute of plasma class
        self.shotlabel = shotlabel
        self.inp = read_infile(self.shotlabel)
        self.ntrl_switch = self.inp.ntrl_switch

    def coreonly(self):
        ntrl_switch = 0
        self.core = exp_core_brnd(self.inp, ntrl_switch) if self.inp.exp_inp else mil_core_brnd(self.inp, ntrl_switch)

    def coreandiol(self):
        ntrl_switch = 0
        self.core = exp_core_brnd(self.inp, ntrl_switch) if self.inp.exp_inp else mil_core_brnd(self.inp, ntrl_switch)
        self.tiol = thermaliol(self.inp, self.core)
        self.fiol = fastiol(self.inp, self.core)
        
    def coreandimp(self):
        ntrl_switch = 0
        self.core = exp_core_brnd(self.inp, ntrl_switch) if self.inp.exp_inp else mil_core_brnd(self.inp, ntrl_switch)
        self.imp = imp_rad(self.inp, self.core)
        
    def ntrlsonly(self):
        ntrl_switch = self.ntrl_switch
        if ntrl_switch == 1:  # neutrals output file already exists. Read it in. No need to run neutpy.
            self.core = exp_core_brnd(self.inp, ntrl_switch) if self.inp.exp_inp else mil_core_brnd(self.inp, ntrl_switch)
            self.ntrl = read_ntrl_data(self.inp, self.core)
        elif ntrl_switch == 2:  # need to run neutpy
            self.core = exp_core_brnd(self.inp, ntrl_switch) if self.inp.exp_inp else mil_core_brnd(self.inp,ntrl_switch)
            self.sol = exp_sol_brnd(self.inp, self.core) if self.inp.exp_inp else mil_sol_brnd(self.inp)
            self.pfr = exp_pfr_brnd(self.inp, self.core) if self.inp.exp_inp else mil_pfr_brnd(self.inp)
            self.ntrl = exp_neutpy_prep(self.inp, self.core, self.sol, self.pfr)
        
    def coreandnbi(self):
        ntrl_switch = 0
        self.core = exp_core_brnd(self.inp, ntrl_switch) if self.inp.exp_inp else mil_core_brnd(self.inp, ntrl_switch)
        self.nbi = beamdep(self.inp, self.core)
        
    def therm_instab(self):
        ntrl_switch = 1
        if ntrl_switch == 1:  # neutrals output file already exists. Read it in. No need to run neutpy.
            self.core = exp_core_brnd(self.inp, ntrl_switch) if self.inp.exp_inp else mil_core_brnd(self.inp, ntrl_switch)
            self.ntrl = read_ntrl_data(self.inp, self.core)
        elif ntrl_switch == 2:  # need to run neutpy
            self.core = exp_core_brnd(self.inp, ntrl_switch) if self.inp.exp_inp else mil_core_brnd(self.inp, ntrl_switch)
            self.sol = exp_sol_brnd(self.inp, self.core) if self.inp.exp_inp else mil_sol_brnd(self.inp)
            self.pfr = exp_pfr_brnd(self.inp, self.core) if self.inp.exp_inp else mil_pfr_brnd(self.inp)
            self.ntrl = exp_neutpy_prep(self.inp, self.core, self.sol, self.pfr)
        self.nbi = beamdep(self.inp, self.core)
        self.imp = imp_rad(self.inp, self.core)
        # self.rtrn = rad_trans(self.inp,self.core,self.tiol,self.fiol,self.ntrl,self.nbi)
        # self.ti = thermal_inst(self.inp,self.core,self.nbi,self.imp,self.ntrl)
        self.dl = dens_lim(self.inp, self.core, self.nbi, self.imp, self.ntrl)
        self.mar = marfe(self.inp, self.core, self.nbi, self.imp, self.ntrl)

    def allthethings(self):
        ntrl_switch = 1
        self.core = exp_core_brnd(self.inp, ntrl_switch) if self.inp.exp_inp else mil_core_brnd(self.inp, ntrl_switch)
        self.sol = exp_sol_brnd(self.inp, self.core) if self.inp.exp_inp else mil_sol_brnd(self.inp)
        self.pfr = exp_pfr_brnd(self.inp, self.core) if self.inp.exp_inp else mil_pfr_brnd(self.inp)
        self.ntrl = exp_neutpy_prep(self.inp, self.core, self.sol, self.pfr)
        self.imp = imp_rad(self.inp, self.core)
        self.tiol = thermaliol(self.inp, self.core)
        self.fiol = fastiol(self.inp, self.core)
        self.nbi = beamdep(self.inp, self.core)
        self.rtrn = rad_trans(self.inp, self.core, self.tiol, self.fiol, self.ntrl, self.nbi)
        # self.ti = thermal_inst(self.inp,self.core,self.nbi,self.imp,self.rtrn)

    def plotstuff(self):
        # self.plots = gt3plots(self)
        pass


class read_ntrl_data():
    def __init__(self,inp,core):
        print 'reading ntrl data'
        ntrl_data = np.loadtxt(inp.neutfile_loc,skiprows=1)
        self.ntrl_R = ntrl_data[:, 0]
        self.ntrl_Z = ntrl_data[:, 1]
        self.n_n_slow = ntrl_data[:, 2]
        self.n_n_thermal = ntrl_data[:, 3]
        self.n_n_total = ntrl_data[:, 4]
        self.izn_rate_slow = ntrl_data[:, 5]
        self.izn_rate_thermal = ntrl_data[:, 6]
        self.izn_rate_total = ntrl_data[:, 7]
        
        n_n_slow = griddata(np.column_stack((self.ntrl_R,self.ntrl_Z)),
                            self.n_n_slow,
                            (core.R, core.Z),
                            method='linear')
        n_n_thermal = griddata(np.column_stack((self.ntrl_R,self.ntrl_Z)),
                               self.n_n_thermal,
                               (core.R, core.Z),
                               method='linear')
        izn_rate_slow = griddata(np.column_stack((self.ntrl_R,self.ntrl_Z)),
                                 self.izn_rate_slow,
                                 (core.R, core.Z),
                                 method='linear')
        izn_rate_thermal = griddata(np.column_stack((self.ntrl_R,self.ntrl_Z)),
                                    self.izn_rate_thermal,
                                    (core.R, core.Z),
                                    method='linear')
        
        core.update_ntrl_data(n_n_slow, n_n_thermal, izn_rate_slow, izn_rate_thermal)


if __name__ == "__main__":
    myshot = gt3('144977_3000/togt3_d3d_144977_3000')
    # myshot.coreonly()
    # myshot.coreandiol()
    myshot.therm_instab()
    # myshot.ntrlsonly()
    # plt.axis('equal')
    # plt.contourf(myshot.core.R,myshot.core.Z,np.log10(myshot.core.n_n_total),500)
    # plt.colorbar()
    # sys.exit()
    
    # fig1 = plt.figure(figsize=(6,8))
    # ax1 = fig1.add_subplot(1,1,1)
    # ax1.axis('equal')
    # ax1.contour(myshot.core.R,myshot.core.Z,myshot.core.rho,10)
    # ax1.plot(myshot.core.R[-1,:],myshot.core.Z[-1,:])
    # ax1.plot(myshot.inp.lim_vertex_closed[:,0],myshot.inp.lim_vertex_closed[:,1])
    # ax1.plot(myshot.inp.sep_exp_closed[:,0],myshot.inp.sep_exp_closed[:,1])
    
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
        ax8.set_title(r'$v_{\theta,C}$', fontsize=fontsize)
        ax8.plot(myshot.core.rho[:, 0], myshot.core.vpolC[:, 0], lw=2, color='black')
        num += 1
    except:
        pass
        
    try:
        ax9 = fig2.add_subplot(rows, cols, num)
        ax9.set_title(r'$v_{\phi,C}$', fontsize=fontsize)
        ax9.plot(myshot.core.rho[:, 0], myshot.core.vtorC[:, 0], lw=2, color='black')
        num += 1
    except:
        pass
        
    try:
        ax10 = fig2.add_subplot(rows, cols, num)
        ax10.set_title(r'$F_{orb}$', fontsize=fontsize)
        ax10.set_xlim(0.9, 1.0)
        ax10.plot(myshot.core.rho[:, 0], myshot.tiol.F_orb_1D, lw=2, color='black')
        num += 1
    except:
        pass
            
    try:
        ax11 = fig2.add_subplot(rows, cols, num)
        ax11.set_title(r'$M_{orb}$', fontsize=fontsize)
        ax11.set_xlim(0.9, 1.0)
        ax11.plot(myshot.core.rho[:, 0], myshot.tiol.M_orb_1D, lw=2, color='black')
        num += 1
    except:
        pass
            
    try:
        ax12 = fig2.add_subplot(rows, cols, num)
        ax12.set_title(r'$E_{orb}$', fontsize=fontsize)
        ax12.set_xlim(0.9, 1.0)
        ax12.plot(myshot.core.rho[:, 0], myshot.tiol.E_orb_1D, lw=2, color='black')
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