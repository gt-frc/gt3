#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 10:21:23 2018

@author: max
"""

class marfe():
    """
    """
    def __init__(self,inp,brnd,nbi,imp,ntrls):
        pass
    
    def calc_marfe_denlim(self,inp,brnd,nbi,imp,ntrls):
        """
        """
        f_cond          = 0
        Q_perp          = 0
        nu              = 0
        L_T             = 0
        C2              = 0
        L_n             = 0
        fz              = 0
        Lz              = 0
        dLzdT           = 0
        f0              = 0
        f0c             = 0
        E_ion           = 0
        sigv_ion        = 0
        T               = 0
        dsigv_ion_dT    = 0
        sigv_cx         = 0
        sigv_el         = 0
        dsigv_cxel_dT   = 0
        
        
        t1 = f_cond * Q_perp * (nu * L_T**-2 + (C2-1)*L_n**-1)
        t2 = fz*((nu + 1 - C2)*Lz/T - dLzdT)
        t3 = f0  * (E_ion * sigv_ion / T * (nu - T / sigv_ion * dsigv_ion_dT))
        t4 = f0c * (3/2*(sigv_cx + sigv_el) * (nu-1-T*dsigv_cxel_dT/(sigv_cx + sigv_el)))
        
        n_marfe = t1 / (T*(t2 + t3 + t4))
        
        print n_marfe