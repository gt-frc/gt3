#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 14:15:35 2018

@author: max
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
from math import sqrt,exp
from scipy.special import jv


class thermal_inst():
    #def __init__(self,inp,brnd,nbi,imp,jro):
    def __init__(self,inp,brnd,nbi,imp):
        sys.dont_write_bytecode = True
        self.therm_collapse(inp,brnd,nbi,imp)
        pass
    
    def J0_weight(X, brnd=brnd, inp=inp):
        r = brnd.r[:,0]
        a = inp.a
        numerator = np.sum(r*X[:,0]*jv(0,5.5*r/a)*a/(len(r)-1.0))
        denominator = np.sum(r*jv(0,5.5*r/a)*a/(len(r)-1.0))
        return numerator / denominator
    
    def therm_collapse_eq14(self,inp,brnd,nbi,imp):

        def therm_cond(inp,brnd):
            taue = 1.0
            chi = inp.a**2.0 / taue
            return np.average(brnd.ni)*chi

        def dHdT(inp,brnd,nbi,imp):
            H_ohm = 0.0
            H_aux = (nbi.pNBi_dvol + nbi.pNBe_dvol) * 1.0E6 # in J * m^-3 * s^-1
            return 3.0/2.0 * (1.0 / brnd.Ti_J) * ( H_ohm - H_aux )

        def dhdT(inp,brnd,nbi,imp):
            return dHdT(inp,brnd,nbi,imp) / np.average(brnd.ni)

        fusion_on = 0
        if fusion_on:
            U_alpha = 1.0E3
        else:
            U_alpha = 0.0

        #Equation 14
        numerator1 = therm_cond(inp,brnd) * (5.5/inp.a)**2 #+ J0_weight(-dHdT(inp,brnd,nbi,imp))
        denominator1 = J0_weight(0.25 * U_alpha * brnd.dsigv_fus_dT + brnd.fracz*(-imp.brnd_dEmiss_dT))
        self.n_therm_collapse1 = np.sqrt(numerator2 / denominator2)

        print 'n_collapse = ',self.n_therm_collapse1

    def therm_collapse_eq44(self,inp,brnd,nbi,imp):
        def chi():
        

        