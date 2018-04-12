# -*- coding: utf-8 -*-
"""
Created on Sat Aug  5 19:48:34 2017

@author: max
"""
from __future__ import division
import numpy as np
from math import pi, sqrt
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.constants import e
import matplotlib as mpl
import matplotlib.ticker as mtick
from scipy.stats import maxwell
from scipy.special import gamma, gammaincc
from scipy.interpolate import InterpolatedUnivariateSpline, griddata
from scipy.integrate import quad, quadrature
import sys
import math

m_d = 3.343583719e-27
m_t = 5.006e-27

m = m_t
class thermaliol():
    '''background plasma stuff.'''
    def __init__(self,inp,brnd):
        sys.dont_write_bytecode = True 
        self.calctiol(inp,brnd)

    def calctiol(self,inp,brnd):
        '''Calculate thermal Ion Orbit Loss.'''
        numcos = 40
        coslist = np.linspace(-1,1,num=numcos)
        Tprofile = brnd.Ti_kev.T[0]

        polpts = len(brnd.r[-1])
        radpts = len(brnd.r.T[-1])
        
        #THE FOLLOWING ARRAYS ARE 4-DIMENSIONAL ARRAYS
        #[ LAUNCH THETA POSITION , LAUNCH ANGLE COSINE,  LAUNCH r  , EXIT THETA POSITION  ]
        #NOTE TO FUTURE DEVELOPERS: IF YOU TRY TO LOOP OVER THE PLASMA POINTS, LAUNCH ANGLES, AND
        #EXIT LOCATIONS IN PYTHON THE WAY YOU WOULD IN C OR FORTRAN, IT'S GOING TO TAKE FOREVER. 
        #ALTHOUGH THESE ARRAYS TAKE MORE MEMORY THAT I'D LIKE, IT'S CURRENTLY NECESSARY TO DO IT THIS WAY.
        #MAYBE SOMETHING TO IMPROVE ON IN THE FUTURE.
        
        #Launch point values
        r0          = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.ones(numcos)[:,None,None],brnd.r)[-1]
        B0          = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.ones(numcos)[:,None,None],brnd.B_tot)[-1]
        f0          = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.ones(numcos)[:,None,None],brnd.f_phi)[-1]
        Psi0        = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.ones(numcos)[:,None,None],brnd.Psi)[-1]
        phi0        = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.ones(numcos)[:,None,None],brnd.E_pot)[-1]*1.0E3 #now in volts
        xi0         = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.linspace(-1,1,num=numcos)[:,None,None],np.ones(brnd.R.shape))[1]

        #Destination Point Values
        R1          = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.ones(numcos)[:,None,None],np.ones(radpts)[:,None],np.ones(polpts)[:],brnd.R[-1][:,None,None,None])[-1]
        f1          = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.ones(numcos)[:,None,None],np.ones(radpts)[:,None],np.ones(polpts)[:],brnd.f_phi[-1][:,None,None,None])[-1]
        B1          = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.ones(numcos)[:,None,None],np.ones(radpts)[:,None],np.ones(polpts)[:],brnd.B_tot[-1][:,None,None,None])[-1]
        Psi1        = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.ones(numcos)[:,None,None],np.ones(radpts)[:,None],np.ones(polpts)[:],brnd.Psi[-1][:,None,None,None])[-1]
        phi1        = np.broadcast_arrays(np.ones(polpts)[:,None,None,None],np.ones(numcos)[:,None,None],np.ones(radpts)[:,None],np.ones(polpts)[:],brnd.E_pot[-1][:,None,None,None])[-1]*1.0E3 #now in volts
        
        a = (np.abs(B1/B0)*f0/f1)**2 - 1 + (1 - xi0**2)*np.abs(B1/B0)
        b = 2*e*(Psi0-Psi1)/(R1*m*f1) * np.abs(B1/B0)*f0/f1*xi0
        c = (e*(Psi0-Psi1)/(R1*m*f1))**2 - 2*e*(phi0-phi1)/m

        v_sep_1 = (-b + np.sqrt(b**2 - 4*a*c))/(2*a)
        v_sep_2 = (-b - np.sqrt(b**2 - 4*a*c))/(2*a)
        
        v_sep = np.zeros(r0.shape)
        v_sep = np.where(
                np.logical_and(
                        np.logical_or(v_sep_1<=0,np.isnan(v_sep_1)),
                        np.logical_or(v_sep_2<=0,np.isnan(v_sep_2))
                        ),
                np.nan,v_sep)
        v_sep = np.where(
                np.logical_or(
                        np.logical_and(v_sep_1>0,np.logical_or(v_sep_2<=0,np.isnan(v_sep_2))),
                        np.logical_and(np.logical_and(v_sep_1>0,v_sep_2>0),v_sep_1<v_sep_2)
                        )
                ,v_sep_1,v_sep)
        v_sep = np.where(
                np.logical_or(
                        np.logical_and(v_sep_2>0,np.logical_or(v_sep_1<=0,np.isnan(v_sep_1))),
                        np.logical_and(np.logical_and(v_sep_1>0,v_sep_2>0),v_sep_2<=v_sep_1)
                        )
                ,v_sep_2,v_sep)
        
        ## PREP FOR FLOSS, MLOSS, AND ELOSS CALCUALTIONS
        v_sep_min = np.nanmin(np.nanmin(v_sep,axis=0),axis=2).T
        v_sep_min[-1] = 0
        T_matrix = np.zeros(v_sep_min.shape)
        for indx,row in enumerate(T_matrix):
            T_matrix[indx,:] = Tprofile[indx]
        zeta_matrix = np.zeros(v_sep_min.shape)
        for indx,column in enumerate(zeta_matrix.T):
            zeta_matrix[:,indx] = coslist[indx]
        
        eps_min = m * v_sep_min**2 / (2*T_matrix*1E3*1.6021E-19)
        
        ## F_orb calculation
        integrand = gammaincc(3/2,eps_min)
        self.F_orb_1D = np.sum(integrand,axis=1)*(2./(numcos-1.))/(2.*gamma(3./2.))
        self.F_orb = np.repeat(self.F_orb_1D.reshape(-1,1),inp.thetapts,axis=1)
        self.F_orb_C = np.zeros(self.F_orb.shape) #TODO:

        ## M_orb calculation
        integrand = zeta_matrix*gammaincc(2.,eps_min)
        self.M_orb_1D = np.sum(integrand,axis=1)*(2./(numcos-1.))/(2.*gamma(2.))
        self.M_orb = np.repeat(self.M_orb_1D.reshape(-1,1),inp.thetapts,axis=1)
        self.M_orb_C = np.zeros(self.F_orb.shape) #TODO:

        ## E_orb calculation
        integrand = gammaincc(5/2,eps_min)
        self.E_orb_1D = np.sum(integrand,axis=1)*(2./(numcos-1.))/(5./2.*gamma(2.))
        self.E_orb = np.repeat(self.E_orb_1D.reshape(-1,1),inp.thetapts,axis=1)
        self.E_orb_C = np.zeros(self.F_orb.shape) #TODO:
        
        #Temporary fix for nan values that mess things up downstream
        #TODO: FIX THIS
        self.F_orb = np.nan_to_num(self.F_orb)
        self.F_orb_C = np.nan_to_num(self.F_orb_C)
        self.M_orb = np.nan_to_num(self.M_orb)
        self.M_orb_C = np.nan_to_num(self.M_orb_C)
        self.E_orb = np.nan_to_num(self.E_orb)
        self.E_orb_C = np.nan_to_num(self.E_orb_C)
        
        iolplot=0
        if iolplot==1:
            fig = plt.figure(figsize=(12,8))
            fig.suptitle('IOL a,b,c in DIII-D with cos:-10', fontsize=24)
            ax1 = fig.add_subplot(131)
            ax1.set_title(r'$v_{sep-1}$',fontsize=20)
            ax1.set_ylabel(r'R',fontsize=20)
            ax1.set_xlabel(r'Z',fontsize=20)
            ax1.grid(b=True,which='both',axis='both')
            CS = ax1.contourf(brnd.R,brnd.Z,v_sep_1[:,-10,:,0].T, 500)
            plt.colorbar(CS)
            
            ax2 = fig.add_subplot(132)
            ax2.set_title(r'$v_{sep-2}$',fontsize=20)
            ax2.set_ylabel(r'R',fontsize=20)
            ax2.set_xlabel(r'Z',fontsize=20)
            ax2.grid(b=True,which='both',axis='both')
            CS = ax2.contourf(brnd.R,brnd.Z,v_sep_2[:,-10,:,0].T, 500)
            plt.colorbar(CS)
            
            ax3 = fig.add_subplot(133)
            ax3.set_title(r'$v_{sep}$',fontsize=20)
            ax3.set_ylabel(r'R',fontsize=20)
            ax3.set_xlabel(r'Z',fontsize=20)
            ax3.grid(b=True,which='both',axis='both')
            #CS = ax3.pcolor(brnd.R,brnd.Z,v_sep[:,10,:,0].T,vmin=0, vmax=1E8)
            CS = ax3.contourf(brnd.R,brnd.Z,v_sep[:,-10,:,0].T, 500)
            plt.colorbar(CS)
    
            plt.tight_layout()
            fig.subplots_adjust(top=0.84)
        return