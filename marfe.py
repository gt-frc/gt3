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

class marfe():
    """
    """
    def __init__(self,inp,brnd,nbi,imp,ntrl):
        self.calc_marfe_denlim(inp,brnd,nbi,imp,ntrl)
    
    def calc_marfe_denlim(self,inp,brnd,nbi,imp,ntrl):
        """
        """
        chi_r           = 1.0E3 / (100.0*100.0) #in m^3
        nu              = 5/2 #??
        L_T             = 0.1*100.0 #brnd.L_Ti_J
        C2              = 0.0 #??
        L_n             = 0.01 * 100.0 #brnd.L_ni
        fz              = brnd.fracz
        Lz              = imp.brnd_emissivity
        dLzdT           = 100.0*imp.brnd_dEmiss_dT
        f0              = ntrl.nn / brnd.ni
        f0c             = ntrl.nn_s / brnd.ni
        E_ion           = 15.0 * 1.6021E-19 #13.59844 * 1.6021E-19
        sigv_ion        = 5.0E-14 #brnd.sigv_ion
        T               = brnd.Ti_J
        dsigv_ion_dT    = brnd.dsigv_ion_dT
        sigv_cx         = 5.0E-14 #brnd.svcx
        sigv_el         = 5.0E-14 #brnd.svel
        dsigv_cxel_dT   = brnd.dsvel_dT + brnd.dsvcx_dT
        
        #t1 = f_cond * Q_perp * (nu * L_T**-1 + (C2-1)*L_n**-1)
        #t2 = fz*((nu + 1 - C2)*Lz/T - dLzdT)
        #t3 = f0  * (E_ion * sigv_ion / T * (nu - T / sigv_ion * dsigv_ion_dT))
        #t4 = f0c * (3/2*(sigv_cx + sigv_el) * (nu-1-T*dsigv_cxel_dT/(sigv_cx + sigv_el)))
        #n_marfe = t1 / (T*(t2 + t3 + t4))

        #print 'f0c'
        #print f0c
        #print 
        #print 'sigv_cx'
        #print sigv_cx
        #print
        #print 'sigv_el'
        #print sigv_el

        t1 = chi_r * (nu * L_T**-2 - (1 - C2) * L_T**-1 * L_n**-1)
        t2 = fz*((nu + 1 - C2)*Lz/T - dLzdT)
        t3 = f0  * (E_ion * sigv_ion / T * (nu - T / sigv_ion * dsigv_ion_dT))
        t4 = f0c * (3/2*(sigv_cx + sigv_el) * (nu-1-T*dsigv_cxel_dT/(sigv_cx + sigv_el)))
        n_marfe = t1 / (T*(t2 + t3 + t4))

        n_marfe = np.where(brnd.rho<0.8,np.nan,n_marfe)
        #print 't1 = ',t1
        #print 't2 = ',t2
        #print 't3 = ',t3
        #print 't4 = ',t4

        print 'n_marfe = ',n_marfe
        
        #marfe_fig = plt.figure(figsize=(6,4))
        #ax1 = marfe_fig.add_subplot(1,1,1)
        #ax1.set_xlim(4.5,5.5)
        #ax1.set_yscale('symlog')
        #ax1.set_xscale('log')
        #ax1.set_ylim(-1E34,1E34)
        #ax1.set_xlim(0.8,1)
        #ax1.loglog(new_T_kev,self.emiss_tot_interp2.derivative()(new_T_J))
        #ax1.set_xlabel(r'Electron Temperature ($keV$)')
        #ax1.set_ylabel(r'Carbon Radiation Emissivity ($W*m^3$)')
        #ax1.set_yscale('symlog')
        #ax1.plot(brnd.rho[:,0],t2[:,0],label='t2')
        #ax1.plot(brnd.rho[:,0],t3[:,0],label='t3')
        #ax1.plot(brnd.rho[:,0],t1[:,0],label='t4')
        #for i in range(32):
        #    ax1.plot(brnd.rho[:,0],n_marfe[:,i],label=i)
        #ax1.legend()
        
        #marfe_fig2 = plt.figure(figsize=(6,4))
        #ax1 = marfe_fig2.add_subplot(1,1,1)
        #ax1.axis('equal')
        #ax1.contourf(brnd.R,brnd.Z,n_marfe,500)
        
        #for i,(v1,v2,v3) in enumerate(zip(t2[:,0],t3[:,0],t4[:,0])):
        #    print v1,v2,v3
        #sys.exit()