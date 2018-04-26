#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 20:24:53 2018

@author: max
"""
from __future__ import division
import numpy as np
from scipy.special import jv
import sys
from math import sqrt
class dens_lim():
    #def __init__(self,inp,brnd,nbi,imp,jro):
    def __init__(self,inp,brnd,nbi,imp,ntrl):
        sys.dont_write_bytecode = True
        self.r = brnd.r[:,0]
        self.a = inp.a
        self.eq44(inp,brnd,nbi,imp,ntrl)
        pass

    def eq44(self,inp,brnd,nbi,imp,ntrl):
        def av(X,p=self):
            numerator   = np.sum(p.r*X*jv(0,5.5*p.r/p.a)*p.a/(len(p.r)-1))
            denominator = np.sum(p.r*  jv(0,5.5*p.r/p.a)*p.a/(len(p.r)-1))
            return numerator / denominator

        ni = brnd.ni[:,0]
        Ti = brnd.Ti_J[:,0]
        n0 = ni[0]
        n_av = av(ni)
        f = n0/n_av
        
        chi = 1.0
        chi_hat = av(ni*chi)/n_av
        D = 0.5
        D_hat = av(ni*D)/n_av
        
        a = self.a
        g = ni/n0
        
        fz = 0.05
        Lz = imp.brnd_emissivity[:,0] * 1E1
        dLzdT = imp.brnd_dEmiss_dT[:,0] * 1E1
        
        sigv_fus = brnd.sigv_fus[:,0]
        dsigv_fus_dT = brnd.dsigv_fus_dT[:,0]
        Ua = 0
        
        H_aux = nbi.pNB_tot[:,0] * 1E6
        H_ohm = 0
        dHdT = 3/(2*Ti)*(H_ohm - H_aux)
        dHdn = 0
        
        nn = ntrl.nn[:,0]
        sigv_ion = brnd.sigv_ion[:,0]
        dsigv_ion_dT = brnd.dsigv_ion_dT[:,0]
        dSdn = nn*sigv_ion
        dSdT = ni*nn*dsigv_ion_dT
        
        ya = 3*av(Ti)*(av(dSdn) - (5.5/a)**2*D_hat) - av(dHdn + 2*ni*(1/4*Ua*sigv_fus - fz*Lz))
        yb = 3*av(ni)*(av(dSdn) - (5.5/a)**2*D_hat) + 3*av(Ti)*av(dSdT) - av(dHdT + ni**2 * (1/4*Ua*dsigv_fus_dT + fz*(-dLzdT)))
        yc = 3*av(ni)*av(dSdT)
        
        y1 = -yb*(1+sqrt(1-4*(ya*yc/yb**2)))/(2*ya)
        y2 = -yb*(1-sqrt(1-4*(ya*yc/yb**2)))/(2*ya)
        
        t1y1 = chi_hat * (5.5/a)**2 * av(g) + 2 * y1 * (fz*av(g*Lz) - av(1/4*Ua*g*sigv_fus))
        t2y1 = 2 * av(g**2 * (1/4*Ua*dsigv_fus_dT + fz*(-dLzdT)))
        t3y1 = 4*(av(-dHdT)-y1*av(dHdn)) * av(g**2*(1/4*Ua*dsigv_fus_dT + fz*(-dLzdT)))
        t4y1 = chi_hat * (5.5/a)**2*av(g) + 2*y1*(fz*av(g*Lz) - av(1/4*Ua*g*sigv_fus))**2
        
        t1y2 = chi_hat * (5.5/a)**2 * av(g) + 2 * y2 * (fz*av(g*Lz) - av(1/4*Ua*g*sigv_fus))
        t2y2 = 2 * av(g**2 * (1/4*Ua*dsigv_fus_dT + fz*(-dLzdT)))
        t3y2 = 4*(av(-dHdT)-y2*av(dHdn)) * av(g**2*(1/4*Ua*dsigv_fus_dT + fz*(-dLzdT)))
        t4y2 = chi_hat * (5.5/a)**2*av(g) + 2*y2*(fz*av(g*Lz) - av(1/4*Ua*g*sigv_fus))**2
        
        nlim1 = (1/f) * (t1y1/t2y1) * (1 + sqrt(1 + t3y1 / t4y1))        
        nlim2 = (1/f) * (t1y2/t2y2) * (1 + sqrt(1 + t3y2 / t4y2))        
        nlim3 = (1/f) * (t1y1/t2y1) * (1 - sqrt(1 + t3y1 / t4y1))        
        nlim4 = (1/f) * (t1y2/t2y2) * (1 - sqrt(1 + t3y2 / t4y2))        
        #print 'nlim1 = ',nlim1
        #print 'nlim2 = ',nlim2
        #print 'nlim3 = ',nlim3
        #print 'nlim4 = ',nlim4
        #print 'n_av = ',n_av
        #sys.exit()