#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 08:58:46 2018

@author: max
"""
import numpy as np
from math import pi,cos
import sys
import matplotlib.pyplot as plt
#from math import sqrt
class rad_trans():
    """
    
    Methods:
        calc_xnuc           
        calc_src_rates      
        calc_rad_part_flx   
        calc_y11_y22        
        calc_rotation       
        calc_xnudrag        
        calc_chi_i          
    
    Attributes:
        atnum_D     atomic number of Deuterium (should move somewhere else or use scipy)
        atnum_C     atomic number of Carbon (should move somewhere else or use scipy)
        xnuc_DD     Momentum transfer frequency from Deuterium to Deuterium
        xnuc_DC     Momentum transfer frequency from Deuterium to Carbon
        xnuc_CD     Momentum transfer frequency from Carbon to Deuterium
        xnuc_CC     Momentum transfer frequency from Carbon to Carbon
        snbi        
        snbi_iol    
        zeff        Z_effective (should move to background calculation)
        se          
        delma       
        volm        
        qnbi        
        qnbe        
        anomtorq    
        xmomtor1    
        xmomtor2    
        qie         
        gamma       
        gammahat    
        qheati      
        qhatheati   
        gammaC      
        y11         
        y22         
        vtorDPert   
        intrin_D    
        intrin_D    
        vtorChat    
        vtorD       
        vtorDhat    
        vpolC       
        vpolD       
        xnudrag_D   
        xnudrag_C   
        chi1        
        chi2        
        chi3        
        chi4        
        chi5        
    """
    
    def __init__(self,inp,brnd,tiol,fiol,ntrl,nbi):
        sys.dont_write_bytecode = True 
        self.calc_xnuc(inp,brnd)
        self.calc_src_rates(inp,brnd,nbi,tiol,fiol)
        self.calc_rad_part_flx(inp,brnd,tiol)
        self.calc_y11_y22(inp,brnd)
        self.calc_rotation(inp,brnd,tiol)
        self.calc_xnudrag(inp,brnd)
        self.calc_chi_i(inp,brnd)
        return
    
    #Calculate scattering cross sections
    def calc_xnuc(self,inp,brnd):
        """
        """
        def calc_coul_log(z1,z2,Ti,nC):
            """
            """
            ep0=8.854E-12
            eq=1.6E-19
            # TODO: raise error if carbon density gets below 1E9 anywhere   
            # TODO: also, verify this whole calculation
            return np.log(12*pi*(Ti**1.5)*(ep0/eq)**1.5)/(np.sqrt(nC)*(z2**2)*z1)
        
        self.atnum_D = 1.0
        self.atnum_C = 6.0
        
        xmr_DD = inp.xmas1 * (1 + (inp.xmas1 / inp.xmas1)) * 1.0E3 #xmrArray[0,0]
        xmr_DC = inp.xmas1 * (1 + (inp.xmas1 / inp.xmas2)) * 1.0E3 #xmrArray[0,1]
        xmr_CD = inp.xmas2 * (1 + (inp.xmas2 / inp.xmas1)) * 1.0E3 #xmrArray[1,0]
        xmr_CC = inp.xmas2 * (1 + (inp.xmas2 / inp.xmas2)) * 1.0E3 #xmrArray[1,1]
        
        coul_DD = calc_coul_log(1,1,brnd.Ti_kev,brnd.nC)
        coul_DC = calc_coul_log(1,6,brnd.Ti_kev,brnd.nC)
        coul_CD = calc_coul_log(6,1,brnd.Ti_kev,brnd.nC)
        coul_CC = calc_coul_log(6,6,brnd.Ti_kev,brnd.nC)

        C1 = 1.0/((((4.8E-10)/(1.6E-12))**1.5)*((4.8E-10)**2.5)) 
    
        self.xnuc_DD = 3.34 * (coul_DD * (self.atnum_D**4) * 1.0E-6 * brnd.ni) / \
                    (C1 * np.sqrt(xmr_DD) * (brnd.Ti_ev**1.5))
                                
        self.xnuc_DC = 3.34 * (coul_DC * (self.atnum_D * brnd.zbar2)**2) * 1.0E-6 * brnd.nC / \
                    (C1 * np.sqrt(xmr_DC) * (brnd.Ti_ev**1.5))
                                
        self.xnuc_CD = 3.34 * (coul_CD * (self.atnum_D * brnd.zbar2)**2) * 1.0E-6 * brnd.ni / \
                    (C1 * np.sqrt(xmr_CD) * (brnd.Ti_ev**1.5))
                                
        self.xnuc_CC = 3.34 * (coul_CC * (brnd.zbar2**4)    * 1E-6 * brnd.nC) / \
                    (C1 * np.sqrt(xmr_CC) * (brnd.Ti_ev**1.5))

        #plot, mostly for debugging purposes
        #xnuc_plot = plt.figure(figsize=(6,4))
        ##ax1 = xnuc_plot.add_subplot(1,1,1)
        ##ax1.set_title('WTF is xnuc??')
        ##ax1.set_xlabel(r'$\rho$')
        #ax1.set_ylabel(r'???')
        #ax1.plot(brnd.rho[:,0],self.xnuc_DD[:,0],label='xnuc_DD')
        #ax1.plot(brnd.rho[:,0],self.xnuc_DC[:,0],label='xnuc_DC')
        #ax1.plot(brnd.rho[:,0],self.xnuc_CD[:,0],label='xnuc_CD')
        #ax1.plot(brnd.rho[:,0],self.xnuc_CC[:,0],label='xnuc_CC')
        #ax1.legend()

        return self.xnuc_DD, self.xnuc_DC, self.xnuc_CD, self.xnuc_CC
    
    #calculate particle, momentum, and energy source rates from beams
    def calc_src_rates(self,inp,brnd,nbi,tiol,fiol):
        """
        """
        def calc_snbi(self,inp,brnd,nbi):
            """
            """
            vp = 19.72 * inp.R0_a  * inp.a**2 * 0.5 * (1+brnd.kappa**2) / 2.0 #TODO: both 0.5 and /2.0 ??
            #TODO: verify that ebeam is the correct thing to use in the following. Used to be 'nbi'
            snbi = (0.624E25*inp.pbeam / (vp*inp.ebeam)) * inp.pwrfrac1 * nbi.dep_prof1 + \
                                                       2 * inp.pwrfrac2 * nbi.dep_prof2 + \
                                                       3 * inp.pwrfrac3 * nbi.dep_prof3
            return snbi                                    
        
        #nbi     = 75000.
        self.snbi       = calc_snbi(self,inp,brnd,nbi)
        self.snbi_iol   = self.snbi*(1.0-tiol.F_orb)

        self.zeff       = (brnd.ni * self.atnum_D**2 + brnd.nC * brnd.zbar2**2) / \
                          (brnd.ni * self.atnum_D    + brnd.nC * brnd.zbar2   )

        #TODO: not sure what se is, but I don't think it's meant to be radially dependent
        #TODO: also, this whole volume calculation should probably use the stuff in brnd,
        #      which is much more accurate
        self.se     = np.sqrt(0.5 * (1+brnd.kappa**2))
        self.delma  = 1.0 / (inp.rpts - 1.0)
        self.volm   = 4.0 * pi**2 * inp.R0_a * inp.a * self.delma * self.se
        
        xlam1  = 5.5E17 * (inp.ebeam/1) / (inp.abeam * 0.5 * brnd.ni * self.zeff**0.5)
        atten1 = 1.0 - np.exp(-self.delma / cos(inp.alphain) / xlam1)
        unatten1 = 1-np.flipud(atten1)
    
        xlam2  = 5.5E17 * (inp.ebeam/2) / (inp.abeam * 0.5 * brnd.ni * self.zeff**0.5)
        atten2 = 1.0 - np.exp(-self.delma / cos(inp.alphain) / xlam2)
        unatten2 = 1-np.flipud(atten2)
        
        xlam3  = 5.5E17 * (inp.ebeam/3) / (inp.abeam * 0.5 * brnd.ni * self.zeff**0.5)
        atten3 = 1.0 - np.exp(-self.delma / cos(inp.alphain) / xlam3)
        unatten3 = 1-np.flipud(atten3)   
            
        xpartdot1 = unatten1 * atten1 * inp.pwrfrac1 * inp.pbeam * 1E6 / \
                                          (1.6E-19 * (inp.ebeam/1) * 1E3) / self.volm 
                                          
        xpartdot2 = unatten2 * atten2 * inp.pwrfrac2 * inp.pbeam * 1E6 / \
                                          (1.6E-19 * (inp.ebeam/2) * 1E3) / self.volm
                                          
        xpartdot3 = unatten3 * atten3 * inp.pwrfrac3 * inp.pbeam * 1E6 / \
                                          (1.6E-19 * (inp.ebeam/3) * 1E3) / self.volm
        
        # qnb calculation. Note that NBIeloss = 0, reducing this formula significantly
        
        self.qnbi = (xpartdot1 + xpartdot2/2 + xpartdot3/3) * 1.6E-19 * 1E3 * inp.ebeam  
        self.qnbe = 0 #TODO: This is most definitely not true... MH
        
        NBIspin=1.
        abeam = 2.
        xmbeam= abeam*1.673E-27 
        
        #TODO: verify what fiol.F_orb are actually supposed to be. 
        torque1 = np.sqrt(2 * xmbeam / (1.6E-19 * inp.ebeam * 1E3 / 1)) * inp.rtang * \
                    inp.pwrfrac1 * inp.pbeam * 1E6 * (1-fiol.F_orb)*NBIspin 
        torque2 = np.sqrt(2 * xmbeam / (1.6E-19 * inp.ebeam * 1E3 / 2)) * inp.rtang * \
                    inp.pwrfrac2 * inp.pbeam * 1E6 * (1-fiol.F_orb)*NBIspin
        torque3 = np.sqrt(2 * xmbeam / (1.6E-19 * inp.ebeam * 1E3 / 3)) * inp.rtang * \
                    inp.pwrfrac3 * inp.pbeam * 1E6 * (1-fiol.F_orb)*NBIspin
        
        xmphi1 = unatten1*atten1*torque1/inp.R0_a
        xmphi2 = unatten2*atten2*torque2/inp.R0_a
        xmphi3 = unatten3*atten3*torque3/inp.R0_a
    
        totmomphi = (xmphi1 + xmphi2 + xmphi3) / self.volm
        
        frac = brnd.ni / (brnd.ni + brnd.nC)
        
        # TODO: Anomolous torque
        self.anomtorq = np.zeros(brnd.rho.shape)
        
        self.xmomtor1 = (1-frac) * (totmomphi+self.anomtorq)
        self.xmomtor2 =   frac   * (totmomphi+self.anomtorq)

        #plot, mostly for debugging purposes
        #src_rate_plot = plt.figure(figsize=(6,4))
        #ax1 = src_rate_plot.add_subplot(1,1,1)
        #ax1.set_title('Neutral beam particle source rates (I think)')
        #ax1.set_xlabel(r'$\rho$')
        #ax1.set_ylabel(r'???')
        #ax1.plot(brnd.rho[:,0],self.snbi[:,0],label='snbi')
        #ax1.plot(brnd.rho[:,0],self.snbi_iol[:,0],label='snbi_iol')
        #ax1.legend()
        
        #zeff_plot = plt.figure(figsize=(6,4))
        #ax1 = zeff_plot.add_subplot(1,1,1)
        #ax1.set_title('zeff')
        #ax1.set_xlabel(r'$\rho$')
        #ax1.set_ylabel(r'???')
        #ax1.plot(brnd.rho[:,0],self.zeff[:,0],label='zeff')
        #ax1.legend()
        
        #qnbi_plot = plt.figure(figsize=(6,4))
        #ax1 = qnbi_plot.add_subplot(1,1,1)
        #ax1.set_title('qnbi')
        #ax1.set_xlabel(r'$\rho$')
        #ax1.set_ylabel(r'???')
        #ax1.plot(brnd.rho[:,0],self.qnbi[:,0],label='qnbi')
        #ax1.legend()

        #xmomtor_plot = plt.figure(figsize=(6,4))
        #ax1 = xmomtor_plot.add_subplot(1,1,1)
        #ax1.set_title('xmomtor')
        #ax1.set_xlabel(r'$\rho$')
        #ax1.set_ylabel(r'???')
        #ax1.plot(brnd.rho[:,0],self.xmomtor1[:,0],label='xmomtor1')
        #ax1.plot(brnd.rho[:,0],self.xmomtor2[:,0],label='xmomtor2')
        #ax1.legend()
        
    #calculate gamma, qheati, qheate (with and without IOL)
    def calc_rad_part_flx(self,inp,brnd,tiol):
        """
        """
        eq1=1.6E-19
        ep0 = 8.854E-12 #TODO: WTF is this? Also the 7.9E-42? WTF this whole calculation.        
        couloge = np.log(12.0 * pi * (brnd.Te_kev**1.5) * (ep0/eq1)**1.5  / \
                         np.sqrt(brnd.ni*self.atnum_D + brnd.nC * brnd.zbar2))
        cequil = 7.9e-42 * couloge * self.zeff / inp.xmas1
        self.qie = cequil*brnd.ne*(brnd.Ti_kev-brnd.Te_kev)/(brnd.Te_kev**1.5)
        
        
    #    rhor=[0.+a*1./mesh for a in range(mesh+1)]
    #    for n in range(mesh-2,0,-1):
    #        se=sqrt(0.5*(1.+elong**2))
    #        
    #        rhor[n] = rhor[n+1] - delma/(aminor*se)

        #TODO: Read in xnuioni and xnuati and conver to background mesh
        cxcool          = np.zeros(brnd.rho.shape)
        srprim          = np.zeros(brnd.rho.shape)
        srprimq         = np.zeros(brnd.rho.shape)
        xpon            = np.zeros(brnd.rho.shape)
        xponq           = np.zeros(brnd.rho.shape)
        self.qheati     = np.zeros(brnd.rho.shape)
        self.qhatheati  = np.zeros(brnd.rho.shape)
        self.gamma      = np.zeros(brnd.rho.shape)
        self.gammahat   = np.zeros(brnd.rho.shape)
        self.gammaC     = np.zeros(brnd.rho.shape) #TODO: gammaC
        for i,rhoval in enumerate(brnd.rho[:,0]):
            if i==0:
                cxcool[i,:]         = 0
                srprim[i,:]         = 0
                srprimq[i,:]        = 0
                xpon[i,:]           = 0
                xponq[i,:]          = 0
                self.qheati[i,:]    = 0
                self.qhatheati[i,:] = 0
                self.gamma[i,:]     = 0
                self.gammahat[i,:]  = 0
            else:
                #TODO: figure out what xnuati and xnuioni are. Currently leaving them zeroed out in inp.
                cxcool[i,:]         = 1.5 * (brnd.ni[i,:] + brnd.ni[i-1,:]) / 2 * \
                                        (brnd.Ti_kev[i,:] + brnd.Ti_kev[i-1,:]) / 2 * inp.xk * inp.xnuati
                srprim[i,:]         = self.snbi[i,:] + 0.5*(brnd.ni[i,:]+brnd.ni[i-1,:]) * \
                                        inp.xnuioni * (1. + brnd.fracz[i,:] * brnd.zbar2[i,:])
                                        
                srprimq[i,:]        = self.qnbi[i,:] - cxcool[i,:] - self.qie[i,:]
                xpon[i,:]           = np.exp(-2.0*(tiol.F_orb[i,:]-tiol.F_orb[i-1,:]))
                xponq[i,:]          = np.exp(-1.0*(tiol.E_orb[i,:]-tiol.E_orb[i-1,:]))
                self.gamma[i,:]     = brnd.rho[i-1,:] / brnd.rho[i,:] * self.gamma[i-1,:] \
                                        *    (1.0)   + srprim[i,:] * self.delma #What is delma? Is this correct here?
                self.gammahat[i,:]  = brnd.rho[i-1,:] / brnd.rho[i,:] * self.gamma[i-1,:] \
                                        *  xpon[i,:] + srprim[i,:] * self.delma #What is delma? Is this correct here?
                self.qheati[i,:]    = (brnd.rho[i-1,:]/brnd.rho[i,:]) *self.qheati[i-1,:]    \
                                        *    (1.0)   + srprimq[i,:] * self.delma
                self.qhatheati[i,:] = (brnd.rho[i-1,:]/brnd.rho[i,:]) *self.qhatheati[i-1,:] \
                                        * xponq[i,:] + srprimq[i,:] * self.delma
        #TODO: electron heating
        
        
        
        #gamma_plot = plt.figure(figsize=(6,4))
        #ax1 = gamma_plot.add_subplot(1,1,1)
        #ax1.set_title('gamma')
        #ax1.set_xlabel(r'$\rho$')
        #ax1.set_ylabel(r'???')
        #ax1.plot(brnd.rho[:,0],self.gamma[:,0],label='gamma')
        #ax1.plot(brnd.rho[:,0],self.gammahat[:,0],label='gammahat')
        #ax1.legend()
        
        #qheati_plot = plt.figure(figsize=(6,4))
        #ax1 = qheati_plot.add_subplot(1,1,1)
        #ax1.set_title('qheati')
        #ax1.set_xlabel(r'$\rho$')
        #ax1.set_ylabel(r'???')
        #ax1.plot(brnd.rho[:,0],self.qheati[:,0],label='qheati')
        #ax1.plot(brnd.rho[:,0],self.qhatheati[:,0],label='qhatheati')
        #ax1.legend()
        

                
        return self.gamma,self.gammahat,self.qheati,self.qhatheati
        
    
    #momentum terms (y11,y22)
    def calc_y11_y22(self,inp,brnd):
        """
        """
        self.y11 = self.xmomtor1 + inp.xk * (brnd.ni * inp.ephia + brnd.B_p * self.gammahat)
        self.y22 = self.xmomtor2 + brnd.zbar2 * inp.xk * (brnd.nC * inp.ephia + brnd.B_p * self.gammaC)
        
        #y11y22_plot = plt.figure(figsize=(6,4))
        #ax1 = y11y22_plot.add_subplot(1,1,1)
        #ax1.set_title('y11 and y22')
        #ax1.set_xlabel(r'$\rho$')
        #ax1.set_ylabel(r'???')
        #ax1.plot(brnd.rho[:,0],self.y11[:,0],label='y11')
        #ax1.plot(brnd.rho[:,0],self.y22[:,0],label='y22')
        #ax1.legend()

        
        return self.y11, self.y22

    
    #Calculate intrinsic rotation
    def calc_rotation(self,inp,brnd,tiol):
        """
        """
        def calc_vtordPert(self,inp,brnd):
            """
            """
            torv        = self.vtorChat + self.intrin_C
            delv0       = self.intrin_D - self.intrin_C
            xnudrageff1 = (self.y11+self.y22) / ((brnd.ni * inp.xmas1 + brnd.nC * inp.xmas2) \
                            * torv + brnd.ni * inp.xmas1 * delv0)
            delv1       = (self.y11-brnd.ni * inp.xmas1 * xnudrageff1 * torv) / \
                            (brnd.ni * inp.xmas1 * (self.xnuc_DC + xnudrageff1))
            self.vtorDPert   = torv + delv1
            return self.vtorDPert
        
        self.intrin_D = 2 / np.sqrt(pi) * tiol.M_orb  * np.sqrt(2 * inp.xk * brnd.Ti_kev / inp.xmas1)
        self.intrin_C = 2 / np.sqrt(pi) * tiol.M_orb_C * np.sqrt(2 * inp.xk * brnd.Ti_kev / inp.xmas2)
        self.vtorChat = brnd.vtorC - self.intrin_C

        if hasattr(inp,'vtorD_rho'):
            self.vtorD    = brnd.vtorD #TODO: Need to interpolate this onto brnd.rho
            self.vtorDhat = self.vtorD - self.intrin_D #TODO: Need to interpolate this onto brnd.rho
        else:
            self.vtorD    = calc_vtordPert(self,inp,brnd)
            self.vtorDhat = calc_vtordPert(self,inp,brnd) - self.intrin_D
            
        self.vpolC = brnd.vpolC #TODO: Need to interpolate this onto brnd.rho
        self.vpolD = self.vpolC / 0.4
        
        #intrin_plot = plt.figure(figsize=(6,4))
        #ax1 = intrin_plot.add_subplot(1,1,1)
        #ax1.set_title('intrinsic rotation')
        #ax1.set_xlabel(r'$\rho$')
        #ax1.set_ylabel(r'???')
        #ax1.plot(brnd.rho[:,0],self.intrin_D[:,0],label='intrin_D')
        #ax1.plot(brnd.rho[:,0],self.intrin_C[:,0],label='intrin_C')
        #ax1.legend()
        
        #vtorChat_plot = plt.figure(figsize=(6,4))
        #ax1 = vtorChat_plot.add_subplot(1,1,1)
        #ax1.set_title('vtorChat')
        #ax1.set_xlabel(r'$\rho$')
        #ax1.set_ylabel(r'???')
        #ax1.plot(brnd.rho[:,0],self.vtorChat[:,0],label='vtorChat')
        #ax1.legend()
        
        #vtorDPert_plot = plt.figure(figsize=(6,4))
        #ax1 = vtorDPert_plot.add_subplot(1,1,1)
        #ax1.set_title('vtorDPert')
        #ax1.set_xlabel(r'$\rho$')
        #ax1.set_ylabel(r'???')
        #ax1.plot(brnd.rho[:,0],self.vtorDPert[:,0],label='vtorDPert')
        #ax1.legend()
        
        return self.intrin_D, self.intrin_C, self.vtorChat
        
    def calc_xnudrag(self,inp,brnd):
        """
        """
        # PROBLEM: When nudrag is negative, it can come close to xnuc12 in magnitude
        #          and blow up the pertrubation theory.            
        torv        = self.vtorChat + self.intrin_C
        delv0       = self.intrin_D - self.intrin_C
        xnudrageff1 = (self.y11+self.y22) / ((brnd.ni * inp.xmas1 + brnd.nC * inp.xmas2) \
                        * torv + brnd.ni * inp.xmas1 * delv0)
        delv1       = (self.y11-brnd.ni * inp.xmas1 * xnudrageff1 * torv) / \
                           (brnd.ni * inp.xmas1 * (self.xnuc_DC + xnudrageff1))
        #NOTE: xnudrageff2 isn't actually used anywhere - MH
        xnudrageff2     = (self.y22 + brnd.nC * inp.xmas2 * self.xnuc_CD * delv1) / \
                            (brnd.nC * inp.xmas2 * torv)
        xnudtot1        = (self.y11 + self.y22 - brnd.ni * inp.xmas1 * xnudrageff1 * delv1) / \
                            ((brnd.ni * inp.xmas1 + brnd.nC * inp.xmas2) * torv)
        delv2           = (self.y11 - brnd.ni * inp.xmas1 * xnudtot1 * torv) / \
                            (brnd.ni * inp.xmas1 * (self.xnuc_DC * xnudtot1))
        self.xnudrag_D  = (self.y11 - brnd.ni * inp.xmas1 * self.xnuc_DC * delv2) / \
                            (brnd.ni * inp.xmas1 * (torv + delv2))
        self.xnudrag_C  = (self.y22 + brnd.nC * inp.xmas2 * self.xnuc_CD * delv2) / \
                            (brnd.nC * inp.xmas2 * torv)
                            
        #nudrag_plot = plt.figure(figsize=(6,4))
        #ax1 = nudrag_plot.add_subplot(1,1,1)
        #ax1.set_title('nudrag')
        #ax1.set_xlabel(r'$\rho$')
        #ax1.set_ylabel(r'???')
        #ax1.plot(brnd.rho[:,0],self.xnudrag_D[:,0],label='xnudrag_D')
        #ax1.plot(brnd.rho[:,0],self.xnudrag_C[:,0],label='xnudrag_C')
        #ax1.legend()
            
        return self.xnudrag_D, self.xnudrag_C
        
    #Infer Chi from experimental data
    def calc_chi_i(self,inp,brnd):
        """
        """
        def calc_qcond(self,inp,conv15=False,conv25=False,visc=False,pressure=False):
            """
            """
            def calc_conv15(self,inp,brnd):
                return 1.5 * inp.xk * self.gamma * brnd.Ti_kev

            def calc_conv25(self,inp,brnd):
                return 2.5 * inp.xk * self.gamma * brnd.Ti_kev

            def calc_heatin(self,inp,brnd):
                return self.gamma * 0.5 * inp.xmas1 * (self.vtorD**2 + self.vpolD**2)
            
            def calc_heatvisc(self,inp,brnd):
                se = np.sqrt(0.5 * (1 + brnd.kappa**2))
                ep = inp.a * se / inp.R0_a
                
                fp = brnd.B_p / brnd.B_phi #these vary with theta, not sure about the best way to make 1D - MH
                xnustar11=0.
        
                # QUESTIONABLE CALCULATION COMMENTED OUT
        #        for a,b in zip(data.xnuc[0,1],data.vpolD):
        #            xnustar11=xnustar11+a*abs(data.q95)*data.rmajor/b
        #
        #       Also, why the fuck is it xnuc12 in GTEDGE?! WHAT ARE YOU DOING WITH YOUR LIFE, STACEY?!
        
                xnustar11 = self.xnuc_DD * abs(inp.q95) * inp.R0_a / self.vpolD   
                
                eff = xnustar11/((1+xnustar11)*(ep**1.5+xnustar11))
                vrad1 = self.gammahat / brnd.ni
                  
                eta0 = brnd.ni * inp.xmas1 * self.vpolD * brnd.q * inp.R0_a * eff
                eta4 = brnd.ni * inp.xmas1 * brnd.Ti_kev * inp.xk / (inp.eq1 * abs(brnd.B_phi))
        #       Calculate viscous heating:  a=vtord, b=fp, c = eta0, d=vrad1, f = eta4, g= vpold
        #   TODO: THIs does not match POP17-052504
                asymR = 0.1 / inp.R0_a
                return asymR * self.vtorDhat * (fp * eta0 * vrad1 - 0.5 * eta4 * (4.* self.vtorDhat + self.vpolD)) - \
                        0.5 * self.vpolD * (eta0 * vrad1 + eta4 * (self.vtorDhat + 0.5 * self.vpolD))
                
            qcond = self.qhatheati
            if ((conv15==True) and (conv25==True)):
                raise("Can't use Conv15 and Conv25 simultaneously")
            if conv25==True:
                qcond = qcond - calc_conv25(self,inp,brnd)
            if conv15==True:
                qcond = qcond - calc_conv15(self,inp,brnd)
            if pressure==True:
                qcond = qcond - calc_heatin(self,inp,brnd)
            if visc==True:
                qcond = qcond - calc_heatvisc(self,inp,brnd)
            return qcond

        self.chi1=brnd.exlti/(brnd.ni*brnd.Ti_kev*inp.xk)*calc_qcond(self,inp)
        self.chi2=brnd.exlti/(brnd.ni*brnd.Ti_kev*inp.xk)*calc_qcond(self,inp)
        self.chi3=brnd.exlti/(brnd.ni*brnd.Ti_kev*inp.xk)*calc_qcond(self,inp,conv25=True)
        self.chi4=brnd.exlti/(brnd.ni*brnd.Ti_kev*inp.xk)*calc_qcond(self,inp,conv25=True,pressure=True)
        self.chi5=brnd.exlti/(brnd.ni*brnd.Ti_kev*inp.xk)*calc_qcond(self,inp,conv25=True,pressure=True,visc=True)
        
        #chi_plot = plt.figure(figsize=(6,4))
        #ax1 = chi_plot.add_subplot(1,1,1)
        #ax1.set_title('chi')
        #ax1.set_xlabel(r'$\rho$')
        #ax1.set_ylabel(r'???')
        #ax1.plot(brnd.rho[:,0],self.chi1[:,0],label='chi1')
        #ax1.plot(brnd.rho[:,0],self.chi2[:,0],label='chi2')
        #ax1.plot(brnd.rho[:,0],self.chi2[:,0],label='chi3')
        #ax1.plot(brnd.rho[:,0],self.chi2[:,0],label='chi4')
        #ax1.plot(brnd.rho[:,0],self.chi2[:,0],label='chi5')
        #ax1.legend()
        
        #for i,v in enumerate(self.chi1):
        #    print i,v
            
        #sys.exit()        