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
from thermal_inst import thermal_inst
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
        input (obj):    Contains as attributes those quantities that are 
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
                            
    """
    def __init__(self, shotlabel=None):
        sys.dont_write_bytecode = True 
        #Create shotlabel as an attribute of plasma class
        self.shotlabel = shotlabel

    def brndonly(self):
        self.input = read_infile(self.shotlabel)
        self.brnd  = background(self.input)
        
    def brndandiol(self):
        self.input = read_infile(self.shotlabel)
        self.brnd  = background(self.input)
        self.tiol  = thermaliol(self.input,self.brnd)
        self.fiol  = fastiol(self.input,self.brnd)
        
    def brndandimp(self):
        self.input = read_infile(self.shotlabel)
        self.brnd  = background(self.input)
        self.imp   = imp_rad(self.input,self.brnd)
        
    def brndandntrls(self):
        self.input = read_infile(self.shotlabel)
        self.brnd  = background(self.input)
        self.ntrl  = neutprep(self.brnd,self.input)
        
    def brndandnbi(self):
        self.input = read_infile(self.shotlabel)
        self.brnd  = background(self.input)
        self.nbi   = beamdep(self.input,self.brnd)

    def noneutrals(self):
        self.input = read_infile(self.shotlabel)
        self.brnd  = background(self.input)
        self.tiol  = thermaliol(self.input,self.brnd)
        self.fiol  = fastiol(self.input,self.brnd)
        #self.ntrl  = neutprep(self.brnd,self.input)
        self.ntrl = 0
        self.nbi   = beamdep(self.input,self.brnd)
        self.jro   = jrostuff(self.input,self.brnd,self.tiol,self.fiol,self.ntrl,self.nbi)

    def allthethings(self):
        self.input = read_infile(self.shotlabel)
        self.brnd  = background(self.input)
        self.imp   = imp_rad(self.input,self.brnd)
        self.tiol  = thermaliol(self.input,self.brnd)
        self.fiol  = fastiol(self.input,self.brnd)
        self.ntrl  = neutprep(self.brnd,self.input)
        self.nbi   = beamdep(self.input,self.brnd)
        self.jro   = jrostuff(self.input,self.brnd,self.tiol,self.fiol,self.ntrl,self.nbi)
        self.ti    = thermal_inst(self.inp,self.brnd,self.nbi,self.imp,self.jro)

    def plotstuff(self):
        #self.plots = gt3plots(self)
        pass

if __name__ == "__main__":
    myshot = gt3('togt3_d3d_118888_1525')
    myshot.brndandimp()
    
    fig1 = plt.figure(figsize=(8,12))
    ax1 = fig1.add_subplot(1,1,1)
    ax1.axis('equal')
    ax1.contour(myshot.brnd.R,myshot.brnd.Z,myshot.brnd.r,25)
    ax1.plot(myshot.brnd.R[-1,:],myshot.brnd.Z[-1,:])
    ax1.plot(myshot.input.lim_vertex[:,0],myshot.input.lim_vertex[:,1])
    
    fig2 = plt.figure(figsize=(12,16))
    ax1 = fig2.add_subplot(4,3,1)
    ax1.set_title(r'$n_i$', fontsize=18)
    ax1.plot(myshot.brnd.rho[:,0],myshot.brnd.ni[:,0],lw=2,color='black')
    
    ax2 = fig2.add_subplot(4,3,2)
    ax2.set_title(r'$n_e$', fontsize=18)
    ax2.plot(myshot.brnd.rho[:,0],myshot.brnd.ne[:,0],lw=2,color='black')
    
    ax3 = fig2.add_subplot(4,3,3)
    ax3.set_title(r'$T_i$', fontsize=18)
    ax3.plot(myshot.brnd.rho[:,0],myshot.brnd.Ti_kev[:,0],lw=2,color='black')
    
    ax4 = fig2.add_subplot(4,3,4)
    ax4.set_title(r'$T_e$', fontsize=18)
    ax4.plot(myshot.brnd.rho[:,0],myshot.brnd.Te_kev[:,0],lw=2,color='black')
    
    ax5 = fig2.add_subplot(4,3,5)
    ax5.set_title(r'$j_r$', fontsize=18)
    ax5.plot(myshot.brnd.rho[:,0],myshot.brnd.j_r[:,0],lw=2,color='black')
    
    ax6 = fig2.add_subplot(4,3,6)
    ax6.set_title(r'$E_r$', fontsize=18)
    ax6.plot(myshot.brnd.rho[:,0],myshot.brnd.E_r[:,0],lw=2,color='black')
    
    ax7 = fig2.add_subplot(4,3,7)
    ax7.set_title(r'$fracz$', fontsize=18)
    ax7.plot(myshot.brnd.rho[:,0],myshot.brnd.fracz[:,0],lw=2,color='black')
    
    ax8 = fig2.add_subplot(4,3,8)
    ax8.set_title(r'$v_{\theta,C}$', fontsize=18)
    ax8.plot(myshot.brnd.rho[:,0],myshot.brnd.vpolC[:,0],lw=2,color='black')
    
    ax9 = fig2.add_subplot(4,3,9)
    ax9.set_title(r'$v_{\phi,C}$', fontsize=18)
    ax9.plot(myshot.brnd.rho[:,0],myshot.brnd.vtorC[:,0],lw=2,color='black')
    
    ax10 = fig2.add_subplot(4,3,10)
    ax10.set_title(r'$F_{orb}$', fontsize=18)
    ax10.plot(myshot.brnd.rho[:,0],myshot.tiol.F_orb,lw=2,color='black')    
    plt.tight_layout()
        
    ax11 = fig2.add_subplot(4,3,11)
    ax11.set_title(r'$M_{orb}$', fontsize=18)
    ax11.plot(myshot.brnd.rho[:,0],myshot.tiol.M_orb,lw=2,color='black')    
    plt.tight_layout()
        
    ax12 = fig2.add_subplot(4,3,12)
    ax12.set_title(r'$E_{orb}$', fontsize=18)
    ax12.plot(myshot.brnd.rho[:,0],myshot.tiol.E_orb,lw=2,color='black')    
    plt.tight_layout()
    #for i,v in enumerate(zip(myshot.brnd.ni[:,0],myshot.brnd.rho[-1])):
    #    print i,v
    #plt.plot(myshot.brnd.rho[-1,:],myshot.brnd.rho[-1,:],lw=2,color='black')
    
    #ax1.plot(myshot.ntrl.sol_lim_pts[:,0],myshot.ntrl.sol_lim_pts[:,1],'o',color='red')
    #ax1.plot(sep[:,0],sep[:,1], color='green',lw=3)
    #CS = ax1.contourf(myshot.brnd.R,myshot.brnd.Z,myshot.brnd.B_p,500) #plot something calculated by miller
    #plt.colorbar(CS)