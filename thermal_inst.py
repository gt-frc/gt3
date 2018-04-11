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
from scipy.interpolate import interp1d, UnivariateSpline


class thermal_inst():
    #def __init__(self,inp,brnd,nbi,imp,jro):
    def __init__(self,inp,brnd,nbi,imp,ntrl):
        sys.dont_write_bytecode = True
        self.therm_collapse_eq44(inp,brnd,nbi,imp,ntrl)
        pass

    def therm_collapse_eq44(self,inp,brnd,nbi,imp,ntrl):
        def J0_weight(X, r, a):
            numerator = np.sum(r*X*jv(0,5.5*r/a)*a/(len(r)-1.0))
            denominator = np.sum(r*jv(0,5.5*r/a)*a/(len(r)-1.0))
            return numerator / denominator
        
        #imported and specified parameters
        nn              = 1.0E16*(brnd.rho[:,0]**5) #np.zeros(ntrl.nn[:,0].shape)
        ni              = brnd.ni[:,0]
        ni0             = ni[0]
        Ti_J            = brnd.Ti_J[:,0]
        r               = brnd.r[:,0]
        a               = inp.a
        taue            = 0.1
        fz              = brnd.fracz[0,0]
        Lz              = imp.brnd_emissivity[:,0]
        dLzdT           = imp.brnd_dEmiss_dT[:,0]
        sigv_fus        = brnd.sigv_fus[:,0]
        sigv_ion        = brnd.sigv_ion[:,0]
        dsigv_fus_dT    = brnd.dsigv_fus_dT[:,0]
        dsigv_ion_dT    = brnd.dsigv_ion_dT[:,0] 

        H_aux           = (nbi.pNBe + nbi.pNBi)[:,0]*1.0E6 #converted from MW/m^3 to W/m^3
        H_ohm           = 0.0 #converted from MW/m^3 to W/m^3

        fusion_on = 0
        if fusion_on==1:
            Ualpha_MeV = 2.375 #for DD
            #Ualpha_MeV = 3.5 #for DT
        else:
            Ualpha_MeV = 0.0
        Ualpha_J = Ualpha_MeV * 1.0E6 * 1.6021E-19

        #parameters calculated for thermal instability
        n_av        = J0_weight(ni, r, a)
        
        g           = ni / ni0
        f           = ni0 / n_av
        chi         = a**2.0 / taue
        D           = a**2.0 / (2*taue)
        
        dHdn        = 0.0
        dHdT        = 3.0/2.0 * (1.0 / Ti_J) * ( H_ohm - H_aux )
        dSdn        = nn * sigv_ion
        dSdT        = 0.0 #ni*nn * dsigv_ion_dT
        
        g_av        = J0_weight(g, r, a)
        L_z_av      = 1.1*fz*J0_weight(g * Lz, r, a)
        fus_av      = J0_weight( 0.25 * Ualpha_J * g * sigv_fus, r, a)
        fus_Lz_av   = J0_weight((0.25 * Ualpha_J * dsigv_fus_dT - fz*dLzdT) * g**2, r, a)

        T_J_av      = J0_weight(Ti_J, r, a)
        chi_hat     = J0_weight(ni * chi, r, a) / n_av
        D_hat       = J0_weight(ni * D, r, a) / n_av
        dHdn_av     = J0_weight( dHdn, r, a)
        dHdT_av     = J0_weight(-dHdT, r, a)
        dSdn_av     = J0_weight( dSdn, r, a)
        dSdT_av     = J0_weight( dSdT, r, a)

        # y calculation
        #TODO: fix the ni vs ne sloppiness here. dSdT_av should really take ne, not ni
        y_a = 3.0 * T_J_av * (dSdn_av - (5.5/a)**2 * D_hat) \
                - J0_weight(dHdn + 2.0*ni* (0.25 * Ualpha_J * sigv_fus     - fz * Lz   ),r,a)
        y_b = 3.0 *  n_av  * (dSdn_av - (5.5/a)**2 * D_hat) + 3.0 * T_J_av * dSdT_av \
                - J0_weight(dHdT + ni**2 * (0.25 * Ualpha_J * dsigv_fus_dT - fz * dLzdT),r,a)
        y_c = 3.0 * n_av * dSdT_av
        
        y1 = -y_b*(1.0 + sqrt(1.0 - 4.0 * y_a * y_c / y_b**2)) / (2.0 * y_a)
        y2 = -y_b*(1.0 - sqrt(1.0 - 4.0 * y_a * y_c / y_b**2)) / (2.0 * y_a)
        
        #calculate density limits
        n_limit_y1 = f**-1 * (chi_hat * (5.5/a)**2 * g_av + 2.0 * y1 * (L_z_av - fus_av)) \
                        / (2.0 * fus_Lz_av) * (1.0 + sqrt(1.0 \
                        +  (4.0*(dHdT_av-y1*dHdn_av)*fus_Lz_av) \
                        /    (chi_hat * (5.5/a)**2 * g_av + 2.0 * y1 * (L_z_av - fus_av)**2)))
                        
        n_limit_y2 = f**-1 * (chi_hat * (5.5/a)**2 * g_av + 2.0 * y2 * (L_z_av - fus_av)) \
                        / (2.0 * fus_Lz_av) * (1.0 + sqrt(1.0 + (4.0*(dHdT_av-y2*dHdn_av)*fus_Lz_av) \
                        / (chi_hat*(5.5/a)**2*g_av + 2.0 * y2 * (L_z_av-fus_av)**2)))


        #attempt at simplified version
        ys_a = 3.0 * T_J_av * (dSdn_av - (5.5/a)**2 * D_hat) \
                - J0_weight(dHdn + 2.0*ni* (0.25 * Ualpha_J * sigv_fus     - fz * Lz   ),r,a)
        ys_b = 3.0 *  n_av  * (dSdn_av - (5.5/a)**2 * D_hat) + 3.0 * T_J_av * dSdT_av \
                - J0_weight(dHdT + ni**2 * (0.25 * Ualpha_J * dsigv_fus_dT - fz * dLzdT),r,a)
        ys_c = 3.0 * n_av * dSdT_av
        
        ys1 = -ys_b*(1.0 + sqrt(1.0 - 4.0 * ys_a * ys_c / ys_b**2)) / (2.0 * ys_a)
        ys2 = -ys_b*(1.0 - sqrt(1.0 - 4.0 * ys_a * ys_c / ys_b**2)) / (2.0 * ys_a)
                        
        #n_simp_y1  = (chi_hat * (5.5/a)**2 * g_av + 2.0 * y1 * L_z_av) \
        #                / (2.0 * fus_Lz_av) * (1.0 + sqrt(1.0 \
        #                +  (4.0*(dHdT_av-y1*dHdn_av)*fus_Lz_av) \
        #                /    (chi_hat * (5.5/a)**2 * g_av + 2.0 * y1 * L_z_av**2))) 
        n_simp_y1  = (chi_hat * (5.5/a)**2 * g_av + 2.0 * y1 * L_z_av) \
                        / (2.0 * fus_Lz_av) * (1.0 + sqrt(1.0 \
                        +  (4.0*(dHdT_av-y1*dHdn_av)*fus_Lz_av) \
                        /    (chi_hat * (5.5/a)**2 * g_av + 2.0 * y1 * L_z_av**2))) 


        #equation 9
        tau_9 = 1.0
        a_9 = 0.8
        chi_9 = a_9**2 / tau_9
        n_9 = 1.0E19
        k_9 = chi_9*n_9
        fz_9 = 0.05
        U_a = 0.0*2.375*1.0E6*1.6021E-19
        dLzdT_9 = imp.brnd_dEmiss_dT_eq9
        dfusdT_9 = brnd.dsigv_fus_dT_eq9
        
        #equation 22
        eta_22 = 5.0E-8
        fz_22 = 0.05
        dLzdT_22 = imp.brnd_dEmiss_dT_eq22
        T_22 = 21.0*1.6021E-19
        print '####################################'
        print '{0:15} {1: 5.3E}'.format('dens_limit_22 = ',sqrt(eta_22/(2.0*fz_22*T_22*(-dLzdT_22)))*(5.5E6)*1E-20)
        #print '{0:15} {1: 5.3E}'.format('dens_limit_22 = ',sqrt(eta_22/(2.0*fz_22*5.0E-33))*(5.5E6)*1E-20)
        print '{0:15} {1: 5.3E}'.format('(-dLzdT_22) = ',(-dLzdT_22))
        print '{0:15} {1: 5.3E}'.format('T_22*(-dLzdT_22) = ',T_22*(-dLzdT_22))
        print '{0:15} {1: 5.3E}'.format('in_sqrt = ',eta_22/(2.0**fz_22*T_22*(-dLzdT_22)))
        print
        print '####################################'
        print '{0:15} {1: 5.3E}'.format('chi_9 = ',chi_9)
        print '{0:15} {1: 5.3E}'.format('(5.5/a)**2  = ',(5.5/a)**2)
        print '{0:15} {1: 5.3E}'.format('fz_9 = ',fz_9)
        print '{0:15} {1: 5.3E}'.format('-dLzdT_9 = ',-dLzdT_9)
        print '{0:15} {1: 5.3E}'.format('dfusdT_9 = ',1.0*dfusdT_9)
        print '{0:15} {1: 5.3E}'.format('0.25*U_a*dfusdT_9 - fz_9*dLzdT_9 = ',0.25*U_a*dfusdT_9 - fz_9*dLzdT_9)
        print '{0:15} {1: 5.3E}'.format('dens_limit_9 = ',sqrt(k_9*(5.5/a_9)**2 / (0.25*U_a*dfusdT_9 - fz_9*dLzdT_9)))
        print '####################################'
        #print 'g = ',g
        print '{0:15} {1: 5.3E}'.format('g_av = ',g_av)
        print '{0:15} {1: 5.3E}'.format('n_av = ',n_av)
        print '{0:15} {1: 5.3E}'.format('T_J_av = ',T_J_av)
        print '{0:15} {1: 5.3E}'.format('f = ',f)
        print
        #print '{0:15} {1: 5.3E}'.format('chi_hat = ',chi_hat)
        #print '{0:15} {1: 5.3E}'.format('D_hat = ',D_hat)
        print
        print '{0:15} {1: 5.3E}'.format('L_z_av = ',L_z_av)
        print '{0:15} {1: 5.3E}'.format('fus_av = ',fus_av)
        print '{0:15} {1: 5.3E}'.format('fus_Lz_av = ',fus_Lz_av)
        print '{0:15} {1: 5.3E}'.format('dHdT_av = ',dHdT_av)
        print '{0:15} {1: 5.3E}'.format('dSdn_av = ',dSdn_av)
        print '{0:15} {1: 5.3E}'.format('dSdT_av = ',dSdT_av)
        print
        print '{0:15} {1: 5.3E}'.format('y_a = ',y_a)
        print '{0:15} {1: 5.3E}'.format('y_b = ',y_b)
        #print '{0:15} {1: 5.3E}'.format('y_c = ',y_c)
        #print '{0:15} {1: 5.3E}'.format('y_a*y_c/y_b = ',y_a*y_c/y_b)
        print '{0:15} {1: 5.3E}'.format('y1 = ',y1)
        print '{0:15} {1: 5.3E}'.format('y2 = ',y2)
        print
        print '{0:15} {1: 5.3E}'.format('n_limit_y1 = ',n_limit_y1)
        print '{0:15} {1: 5.3E}'.format('n_limit_y2 = ',n_limit_y2)
        print '{0:15} {1: 5.3E}'.format('n_simp_y1 = ',n_simp_y1)
        print '####################################'
        print
        print '{0:15} {1: 5.3E}'.format('chi_hat * (5.5/a)**2 * g_av = ',chi_hat * (5.5/a)**2 * g_av)
        print '{0:15} {1: 5.3E}'.format('2.0 * y1 * L_z_av = ',2.0 * y1 * L_z_av)
        print '{0:15} {1: 5.3E}'.format('2.0 * fus_Lz_av = ',2.0 * fus_Lz_av)
        print '{0:15} {1: 5.3E}'.format('last term = ',(1.0 + sqrt(1.0 \
                        +  (4.0*(dHdT_av-y1*dHdn_av)*fus_Lz_av) \
                        /    (chi_hat * (5.5/a)**2 * g_av + 2.0 * y1 * L_z_av**2))) )
        sys.exit()
        