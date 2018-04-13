# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 11:56:27 2017

@author: Max Hill
"""

from read_infile import read_infile
from background import background
from imp_rad import imp_rad
from thermaliol import thermaliol
from fastiol import fastiol
from neutprep import neutprep
from beamdep import beamdep
from jrostuff import jrostuff
#from thermal_inst import thermal_inst
from dens_lim import dens_lim
from marfe import marfe
import matplotlib.pyplot as plt
#from gt3plots import gt3plots
import sys

class gt3():
    """GT3 calculates various tokamak-related quantities

    Methods:
        allthethings
        justneutrals
        plotstuff

    Attributes:
        inp (obj):    Contains as attributes those quantities that are 
                        specified in input files
        brnd (obj):     Contains as attributes those quantities related
                        to the background plasma (i.e. ni, Te, j_r, etc.)
        tiol (obj):     Contains as attributes those quantities related
                        to the thermal ion-orbit loss calculation
        ntrl (obj):     Contains as attributes those quantities related
                        to the neutrals calculation performed by neutpy
        plots (obj):    Contains as attributes things related to plotting
        
    External Dependencies:
        Triangle:       Used to create the triangular mesh for neutpy. Source
                        code and documentation can be found at 
                        https://www.cs.cmu.edu/~quake/triangle.html. Can be
                        installed from Ubuntu repositories as well.
        Neutpy:         
        nbeams:
        adpack:
                            6
    """
    
    def __init__(self, shotlabel=None):
        sys.dont_write_bytecode = True 
        #Create shotlabel as an attribute of plasma class
        self.shotlabel = shotlabel

    def brndonly(self):
        self.inp = read_infile(self.shotlabel)
        self.brnd  = background(self.inp)
        
    def brndandiol(self):
        self.inp = read_infile(self.shotlabel)
        self.brnd  = background(self.inp)
        self.tiol  = thermaliol(self.inp,self.brnd)
        self.fiol  = fastiol(self.inp,self.brnd)
        
    def brndandimp(self):
        self.inp = read_infile(self.shotlabel)
        self.brnd  = background(self.inp)
        self.imp   = imp_rad(self.inp,self.brnd)
        
    def brndandntrls(self):
        self.inp = read_infile(self.shotlabel)
        self.brnd  = background(self.inp)
        self.ntrl  = neutprep(self.inp,self.brnd)
        
    def brndandnbi(self):
        self.inp = read_infile(self.shotlabel)
        self.brnd  = background(self.inp)
        self.nbi   = beamdep(self.inp,self.brnd)

    def noneutrals(self):
        self.inp = read_infile(self.shotlabel)
        self.brnd  = background(self.inp)
        self.tiol  = thermaliol(self.inp,self.brnd)
        self.fiol  = fastiol(self.inp,self.brnd)
        #self.ntrl  = neutprep(self.inp,self.brnd)
        self.ntrl = 0
        self.nbi   = beamdep(self.inp,self.brnd)
        self.jro   = jrostuff(self.inp,self.brnd,self.tiol,self.fiol,self.ntrl,self.nbi)
        
    def therm_instab(self):
        self.inp = read_infile(self.shotlabel)
        self.brnd  = background(self.inp)
        self.nbi   = beamdep(self.inp,self.brnd)        
        self.imp   = imp_rad(self.inp,self.brnd)
        self.ntrl  = neutprep(self.inp,self.brnd)
        #self.jro   = jrostuff(self.inp,self.brnd,self.tiol,self.fiol,self.ntrl,self.nbi)
        #self.ti    = thermal_inst(self.inp,self.brnd,self.nbi,self.imp,self.ntrl)
        self.dl    = dens_lim(self.inp,self.brnd,self.nbi,self.imp,self.ntrl)
        self.mar   = marfe(self.inp,self.brnd,self.nbi,self.imp,self.ntrl)

    def allthethings(self):
        self.inp = read_infile(self.shotlabel)
        self.brnd  = background(self.inp)
        self.imp   = imp_rad(self.inp,self.brnd)
        self.tiol  = thermaliol(self.inp,self.brnd)
        self.fiol  = fastiol(self.inp,self.brnd)
        self.ntrl  = neutprep(self.inp,self.brnd)
        self.nbi   = beamdep(self.inp,self.brnd)
        self.jro   = jrostuff(self.inp,self.brnd,self.tiol,self.fiol,self.ntrl,self.nbi)
        #self.ti    = thermal_inst(self.inp,self.brnd,self.nbi,self.imp,self.jro)

    def plotstuff(self):
        #self.plots = gt3plots(self)
        pass

if __name__ == "__main__":
    myshot = gt3('togt3_d3d_118888_1525')
    #myshot = gt3('togt3_dens_lim_test')
    #myshot.therm_instab()
    myshot.brndandiol()
    
    fig1 = plt.figure(figsize=(6,8))
    ax1 = fig1.add_subplot(1,1,1)
    ax1.axis('equal')
    ax1.contour(myshot.brnd.R,myshot.brnd.Z,myshot.brnd.r,25)
    ax1.plot(myshot.brnd.R[-1,:],myshot.brnd.Z[-1,:])
    ax1.plot(myshot.inp.lim_vertex[:,0],myshot.inp.lim_vertex[:,1])
    
    fontsize=12
    fig2 = plt.figure(figsize=(7,10))
    rows = 4
    cols = 3
    num = 1
    try:
        ax1 = fig2.add_subplot(rows,cols,num)
        ax1.set_title(r'$n_i$', fontsize=fontsize)
        ax1.plot(myshot.brnd.rho[:,0],myshot.brnd.ni[:,0],lw=2,color='black')
        num +=1
    except:
        pass
    
    try:
        ax2 = fig2.add_subplot(rows,cols,num)
        ax2.set_title(r'$n_e$', fontsize=fontsize)
        ax2.plot(myshot.brnd.rho[:,0],myshot.brnd.ne[:,0],lw=2,color='black')
        num +=1
    except:
        pass
        
    try:
        ax3 = fig2.add_subplot(rows,cols,num)
        ax3.set_title(r'$T_i$', fontsize=fontsize)
        ax3.plot(myshot.brnd.rho[:,0],myshot.brnd.Ti_kev[:,0],lw=2,color='black')
        num +=1
    except:
        pass
        
    try:
        ax4 = fig2.add_subplot(rows,cols,num)
        ax4.set_title(r'$T_e$', fontsize=fontsize)
        ax4.plot(myshot.brnd.rho[:,0],myshot.brnd.Te_kev[:,0],lw=2,color='black')
        num +=1
    except:
        pass
        
    try:
        ax5 = fig2.add_subplot(rows,cols,num)
        ax5.set_title(r'$j_r$', fontsize=fontsize)
        ax5.plot(myshot.brnd.rho[:,0],myshot.brnd.j_r[:,0],lw=2,color='black')
        num +=1
    except:
        pass
        
    try:
        ax6 = fig2.add_subplot(rows,cols,num)
        ax6.set_title(r'$E_r$', fontsize=fontsize)
        ax6.plot(myshot.brnd.rho[:,0],myshot.brnd.E_r[:,0],lw=2,color='black')
        num +=1
    except:
        pass
        
    try:
        ax7 = fig2.add_subplot(rows,cols,num)
        ax7.set_title(r'$fracz$', fontsize=fontsize)
        ax7.plot(myshot.brnd.rho[:,0],myshot.brnd.fracz[:,0],lw=2,color='black')
        num +=1
    except:
        pass
        
    try:
        ax8 = fig2.add_subplot(rows,cols,num)
        ax8.set_title(r'$v_{\theta,C}$', fontsize=fontsize)
        ax8.plot(myshot.brnd.rho[:,0],myshot.brnd.vpolC[:,0],lw=2,color='black')
        num +=1
    except:
        pass
        
    try:
        ax9 = fig2.add_subplot(rows,cols,num)
        ax9.set_title(r'$v_{\phi,C}$', fontsize=fontsize)
        ax9.plot(myshot.brnd.rho[:,0],myshot.brnd.vtorC[:,0],lw=2,color='black')
        num +=1
    except:
        pass
        
    try:
        ax10 = fig2.add_subplot(rows,cols,num)
        ax10.set_title(r'$F_{orb}$', fontsize=fontsize)
        ax10.set_xlim(0.9,1.0)
        ax10.plot(myshot.brnd.rho[:,0],myshot.tiol.F_orb_1D,lw=2,color='black')    
        num +=1
    except:
        pass
            
    try:
        ax11 = fig2.add_subplot(rows,cols,num)
        ax11.set_title(r'$M_{orb}$', fontsize=fontsize)
        ax11.set_xlim(0.9,1.0)
        ax11.plot(myshot.brnd.rho[:,0],myshot.tiol.M_orb_1D,lw=2,color='black')    
        num +=1
    except:
        pass
            
    try:
        ax12 = fig2.add_subplot(rows,cols,num)
        ax12.set_title(r'$E_{orb}$', fontsize=fontsize)
        ax12.set_xlim(0.9,1.0)
        ax12.plot(myshot.brnd.rho[:,0],myshot.tiol.E_orb_1D,lw=2,color='black')    
        num +=1
    except:
        pass
            
    try:
        ax13 = fig2.add_subplot(rows,cols,num)
        ax13.set_title(r'$I_{cum}$', fontsize=fontsize)
        ax13.plot(myshot.brnd.rho[:,0],myshot.brnd.I[:,0],lw=2,color='black')    
        num +=1
    except:
        pass
    
    try:
        ax14 = fig2.add_subplot(rows,cols,num)
        ax14.set_title(r'$Z_{eff}$', fontsize=fontsize)
        ax14.plot(myshot.brnd.rho[:,0],myshot.brnd.z_eff[:,0],lw=2,color='black')    
        num +=1
    except:
        pass
    
    try:
        ax15 = fig2.add_subplot(rows,cols,num)
        ax15.set_title(r'$NBI Dep. Prof.$', fontsize=fontsize)
        ax15.plot(myshot.brnd.rho[:,0],myshot.nbi.pNB_tot[:,0],lw=2,color='black')    
        num +=1
    except:
        pass
    plt.tight_layout()
    #plt.plot(myshot.brnd.rho[-1,:],myshot.brnd.rho[-1,:],lw=2,color='black')
    
    #ax1.plot(myshot.ntrl.sol_lim_pts[:,0],myshot.ntrl.sol_lim_pts[:,1],'o',color='red')
    #ax1.plot(sep[:,0],sep[:,1], color='green',lw=3)
    #CS = ax1.contourf(myshot.brnd.R,myshot.brnd.Z,myshot.brnd.B_p,500) #plot something calculated by miller
    #plt.colorbar(CS)