#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  5 16:05:08 2017

@author: max
"""
from __future__ import division
from math import pi, sin, acos, ceil
import numpy as np
from scipy.constants import mu_0, elementary_charge, k
from scipy.interpolate import UnivariateSpline, interp1d, interp2d, griddata
from scipy.integrate import quad
import matplotlib.pyplot as plt
import sys

def PolyArea(x, y):
    """finds the area of a polygon"""
    return 0.5*np.abs(np.dot(x, np.roll(y, 1))-np.dot(y, np.roll(x, 1)))

class mil_core_brnd():
    """Calculates various plasma properties using a modified Miller geometry
    
    Methods:
        createbackround
        xmiller
    
    Attributes:
        r
        theta
        rho
        ni
        ne
        Ti_kev
        Ti_K
        Ti_J
        Te_kev
        Te_K
        Te_J
        nC
        E_pot
        pressure
        j_r
        kappa
        tri
        R
        Z
        diff_vol
        IP
        B_phi
        psi
        psi_norm
        B_p
        B_tot
        f_phi
    """
    def __init__(self, inp, ntrl_switch):
        sys.dont_write_bytecode = True 
        self.miller(inp)
        self.xsec(inp)

    def miller(self, p):
        """Create background plasma using the miller model.
        
        Note:
            
        Args:
        """
        ## CREATE MAIN r AND theta MATRICES
        try:
            rho1d = np.concatenate((np.linspace(0, p.edge_rho*p.a, p.rpts_core, endpoint=False),
                                 np.linspace(p.edge_rho*p.a, p.a, p.rpts_edge)), axis=0)
        except AttributeError:
            try:
                rho1d = np.linspace(0, 1, p.rhopts)
            except AttributeError:
                raise AttributeError("You haven't specified the number of radial points.")

        # define thetapts and xpt
        self.thetapts = int(4 * ceil(float(p.thetapts_approx)/4))+1

        theta1d = np.linspace(0, 2*pi, self.thetapts)
        self.theta, self.rho = np.meshgrid(theta1d, rho1d)
        self.r = self.rho * p.a

        ped_loc = 1.0
        ##########################################################################################
        # CREATE DENSITY, TEMPERATURE, PRESSURE, AND CURRENT DENSITY ARRAYS
        ##########################################################################################
        try:
            self.ni = UnivariateSpline(p.ni_data[:, 0], p.ni_data[:, 1], k=3, s=2.0)(self.rho)
        except AttributeError:
            self.ni = np.where(self.r < ped_loc*p.a,
                               (p.ni0-p.ni9)*(1-self.rho**2)**p.nu_ni + p.ni9, 
                               (p.ni_sep-p.ni9)/(0.1*p.a)*(self.r-ped_loc*p.a)+p.ni9)
        #gradient scale length
        # self.dni_dr  = np.gradient(self.ni, self.r[:, 0], axis=0)
        # self.L_ni = -self.dni_dr / self.ni
        #############################################

        try:
            self.ne = UnivariateSpline(p.ne_data[:, 0], p.ne_data[:, 1], k=3, s=2.0)(self.rho)
        except AttributeError:
            self.ne = np.where(self.r < ped_loc*p.a,
                               (p.ne0-p.ne9)*(1-self.rho**2)**p.nu_ne + p.ne9, 
                               (p.ne_sep-p.ne9)/(0.1*p.a)*(self.r-ped_loc*p.a)+p.ne9)

        #gradient scale length
        # TODO: Redo gradient scale length calculation to be in 2d, possibly with scipy.interpolate.SmoothBivariateSpline
        # self.dne_dr  = np.gradient(self.ne, self.r[:, 0], axis=0)
        # self.L_ne = -self.dne_dr / self.ne
        #############################################

        try:
            self.fracz = UnivariateSpline(p.fracz_data[:, 0], p.fracz_data[:, 1], k=5, s=2.0)(self.rho)
        except AttributeError:
            self.fracz = np.zeros(self.rho.shape) + 0.025
            
        self.nC = self.ne * self.fracz
        
        #gradient scale length
        # self.dnC_dr  = np.gradient(self.nC, self.r[:, 0], axis=0)
        # self.L_nC = -self.dnC_dr / self.nC
        
        self.z_eff = (self.ni*1.0**2 + self.nC*6.0**2) / self.ne
        
        # TODO: calculate z_0 over all charge states from imp_rad.
        # Might need to move this calculation there.
        self.z_0 = self.nC*6.0**2 / self.ni
        #############################################

        try:
            self.Ti_kev = UnivariateSpline(p.Ti_data[:, 0], p.Ti_data[:, 1], k=5, s=2.0)(self.rho)
        except AttributeError:
            self.Ti_kev = np.where(self.rho<ped_loc,
                             (p.Ti0-p.Ti9)*(1-self.rho**2)**p.nu_Ti + p.Ti9, 
                             (p.Ti_sep-p.Ti9)/(0.1*p.a)*(self.rho-ped_loc)+p.Ti9)
        self.Ti_K  = self.Ti_kev * 1.159E7
        self.Ti_ev = self.Ti_kev * 1000
        self.Ti_J  = self.Ti_ev  * elementary_charge

        #gradient scale length
        # self.dTi_J_dr  = np.gradient(self.Ti_J, self.r[:, 0], axis=0)
        # self.L_Ti_J = -self.dTi_J_dr / self.Ti_J
        #############################################

        try:
            self.Te_kev = UnivariateSpline(p.Te_data[:, 0], p.Te_data[:, 1], k=5, s=2.0)(self.rho)
        except AttributeError:
            self.Te_kev = np.where(self.rho<ped_loc,
                             (p.Te0-p.Te9)*(1-self.rho**2)**p.nu_Te + p.Te9, 
                             (p.Te_sep-p.Te9)/(0.1*p.a)*(self.rho-ped_loc)+p.Te9)
        self.Te_K  = self.Te_kev * 1.159E7
        self.Te_ev = self.Te_kev * 1000
        self.Te_J  = self.Te_ev  * elementary_charge
        
        #gradient scale length
        # self.dTe_J_dr  = np.gradient(self.Te_J, self.r[:, 0], axis=0)
        # self.L_Te_J = -self.dTe_J_dr / self.Te_J
        #############################################

        try:
            E_r_fit = UnivariateSpline(p.er_data[:, 0], p.er_data[:, 1], k=5, s=2.0)
            self.E_r = E_r_fit(self.rho)
            self.E_pot = np.zeros(self.rho.shape)
            for i, rhoval in enumerate(rho1d):
                self.E_pot[i] = E_r_fit.integral(rhoval, 1.0)
        except AttributeError:
            raise AttributeError("You need E_r data")
            sys.exit()

        #############################################

        try:
            self.j_r = p.j0*(1-(self.rho)**2)**p.nu_j
        except AttributeError:
            raise AttributeError("You haven't specified a current distribution.")  

        #############################################

        try:
            self.fz1 = UnivariateSpline(p.fz1_data[:, 0], p.fz1_data[:, 1], k=5, s=2.0)(self.rho)
        except AttributeError:
            self.fz1 = 0.025*self.ne

        #############################################

        try:
            self.vpolC = UnivariateSpline(p.vpolC_data[:, 0], p.vpolC_data[:, 1], k=5, s=2.0)(self.rho)
        except AttributeError:
            self.vpolC = 0.0

        #############################################

        try:
            self.vtorC = UnivariateSpline(p.vtorC_data[:, 0], p.vtorC_data[:, 1], k=5, s=2.0)(self.rho)
        except AttributeError:
            self.vtorC = 0.0

        #############################################

        try:
            self.vpolD = UnivariateSpline(p.vpolD_data[:, 0], p.vpolD_data[:, 1], k=5, s=2.0)(self.rho)
        except AttributeError:
            self.vpolD = 0.0

        #############################################

        try:
            self.vtorD = UnivariateSpline(p.vtorD_data[:, 0], p.vtorD_data[:, 1], k=5, s=2.0)(self.rho)
        except AttributeError:
            self.vtorD = 0.0
        #############################################

        try:
            self.q = UnivariateSpline(p.q_data[:, 0], p.q_data[:, 1], k=5, s=2.0)(self.rho)
        except AttributeError:
            self.q = np.zeros(self.rho.shape) #will calculated later with the other miller stuff

        #############################################

        try:
            self.zbar2 = UnivariateSpline(p.zbar2_data[:, 0], p.zbar2_data[:, 1], k=5, s=2.0)(self.rho)
        except AttributeError:
            self.zbar2 = np.zeros(self.rho.shape) + 0.025


        self.pressure = self.ni * k * self.Ti_K
        
        ##########################################################################################
        ## CREATE kappa, tri AND RELATED MATRICES
        ##########################################################################################
        upperhalf = (self.theta >= 0)&(self.theta < pi)
        self.kappa = np.where(upperhalf,
                         p.kappa_up / (p.a**p.s_k_up) * self.r**p.s_k_up, 
                         p.kappa_lo / (p.a**p.s_k_lo) * self.r**p.s_k_lo)
        
        
        ## All we're doing with kappa in this next part is making the derivative between upper and lower
        ## elongation continuous by "smoothing out" the "step function"
        ## using f(x) = tanh(B*sin(x)), where be controlls how smooth or squre the function is.
        ## Plot that function and you'll see what we're doing. This is necessary 
        ## to prevent shafranov shift from producing ugly pictures with high poloidal
        ## resolution. It also makes Richard's stuff easier. Just deal with it 
        ## and don't put this in any papers. It's just a bandaid. We do the same 
        ## thing with triangularity. - MH
        
        #B_kappa = 0.0
        #self.kappa  = (((p.kappa_up / (p.a**p.s_k_up) * self.r**p.s_k_up) - (p.kappa_lo / (p.a**p.s_k_lo) * self.r**p.s_k_lo))/2.0 
        #        * np.tanh(B_kappa*np.sin(self.theta))
        #        + ((p.kappa_up / (p.a**p.s_k_up) * self.r**p.s_k_up) + (p.kappa_lo / (p.a**p.s_k_lo) * self.r**p.s_k_lo))/2.0)
         
        if p.xmil==1:
            self.xpt = np.array([p.xpt_R, p.xpt_Z])
            self.kappa = self.xmiller(self.kappa, p)
            tri_lo = sin(3*pi/2 - acos((self.xpt[0]-p.R0_a)/p.a))
            tri_up = p.tri_up
        else:
            tri_lo = p.tri_lo
            tri_up = p.tri_up
            

        tri = np.where(upperhalf,
                         tri_up * (self.rho)**1,
                         tri_lo * (self.rho)**1)

        # s_tri = np.where(upperhalf,
        #                 self.r*p.tri_up/(p.a*np.sqrt(1-tri)),
        #                 self.r*tri_lo/(p.a*np.sqrt(1-tri)))
        
        ## CALCULATE INITIAL R, Z WITH NO SHAFRANOV SHIFT
        ## (NECESSARY TO GET ESTIMATES OF L_r WHEN CALCULATING SHAFRANOV SHIFT)
        R0 = np.ones(self.rho.shape) * p.R0_a
        self.R = R0 + self.r * np.cos(self.theta+np.arcsin(tri*np.sin(self.theta)))
        self.Z = self.kappa*self.r*np.sin(self.theta)
        
        # THIS CALCULATES A MATRIX OF THE LENGTHS OF EACH SECTION OF EACH FLUX
        # SURFACE AND THEN SUMS THEM TO GET THE PERIMETER IN 2D OF EACH FLUX
        # SURFACE (VALUE OF r).
        self.L_seg = np.sqrt((self.Z-np.roll(self.Z, -1, axis=1))**2 + (self.R-np.roll(self.R, -1, axis=1))**2)
        self.L_seg [:, -1] = 0        
        self.L_r = np.tile(np.sum(self.L_seg, axis=1), (self.thetapts, 1)).T
        
        #CALCULATE CROSS-SECTIONAL AREA CORRESPONDING TO EACH r AND ASSOCIATED
        #DIFFERENTIAL AREAS
        area = np.zeros(self.rho.shape)
        for i in range(0, len(self.rho)):
            area[i, :] = PolyArea(self.R[i, :], self.Z[i, :])
    
        diff_area = area - np.roll(area, 1, axis=0)
        diff_area[0, :]=0
        
        self.diff_vol = diff_area * 2*pi*p.R0_a #approx because it uses R0_a instead of shifted R0
        vol = np.cumsum(self.diff_vol, axis=0)
        
        #Calculate each differential I and sum to get cumulative I
        j_r_ave = np.roll((self.j_r + np.roll(self.j_r, -1, axis=0))/2.0, 1, axis=0)
        j_r_ave[0, :]=0
        diff_I = diff_area * j_r_ave
        self.I = np.cumsum(diff_I, axis=0)
        self.IP = self.I[-1, 0]  
        print 'IP = ',self.IP
        #Calculate B_p_bar
        B_p_bar = mu_0 * self.I / self.L_r
        B_p_bar[0, :]=0
        
        #Calculate li
        li = (np.cumsum(B_p_bar**2 * self.diff_vol, axis=0) / vol) / (2*B_p_bar**2)
        li[0, :]=0
        
        #Calculate beta_p
        beta_p = 2*mu_0*(np.cumsum(self.pressure*self.diff_vol, axis=0)/vol-self.pressure) / B_p_bar**2
    
        #Calculate dR0dr
        self.dR0dr = np.zeros(self.rho.shape)
        self.R0 = np.zeros(self.rho.shape)
    
        f = 2*(self.kappa**2+1)/(3*self.kappa**2+1)*(beta_p+li/2)+1/2*(self.kappa**2-1)/(3*self.kappa**2+1)
        f[0, :] = f[1, :] ############ NEED TO REVISIT, SHOULD EXTRAPOLATE SOMEHOW
        
        self.dR0dr[-1, :] = -2.0*p.a*f[-1, :]/p.R0_a
        self.R0[-1, :] = p.R0_a
        
        for i in range(len(self.rho)-2, -1, -1):
            self.R0[i, :] = self.dR0dr[i+1, :] * (self.r[i, :]-self.r[i+1, :]) + R0[i+1, :]
            self.dR0dr[i, :] = -2.0*self.r[i, :]*f[i, :]/R0[i, :]
        
        #NOW USE UPDATED R0 AND dR0dr to get new R, Z.
        self.R = self.R0 + self.r * np.cos(self.theta+np.arcsin(tri*np.sin(self.theta)))
        self.Z = self.kappa*self.r*np.sin(self.theta) + p.Z0

        #RECALCULATE L_seg and L_r
        self.L_seg = np.sqrt((self.Z-np.roll(self.Z, -1, axis=1))**2 + (self.R-np.roll(self.R, -1, axis=1))**2)
        self.L_seg [:, -1] = 0        
        self.L_r = np.tile(np.sum(self.L_seg, axis=1), (self.thetapts, 1)).T
        
        ## RECALCULATE GRAD-r
        dkappa_dtheta = np.gradient(self.kappa, edge_order=1)[1] * self.thetapts/(2*pi)
        dkappa_dr = np.gradient(self.kappa, edge_order=1)[0] * p.rhopts/p.a
    
        dkappa_dtheta[-1] = dkappa_dtheta[-2]
        dkappa_dr[-1] = dkappa_dr[-2]
    
        dZ_dtheta = np.gradient(self.Z, edge_order=2)[1] * self.thetapts/(2*pi)  # self.r*(self.kappa*np.cos(self.theta)+dkappa_dtheta*np.sin(self.theta))
        dZ_dr = np.gradient(self.Z, edge_order=2)[0] * p.rhopts/p.a  # np.sin(self.theta)*(self.r*dkappa_dr + self.kappa)
        dR_dr = np.gradient(self.R, edge_order=2)[0] * p.rhopts/p.a  # dR0dr - np.sin(self.theta + np.sin(self.theta)*np.arcsin(tri))*(np.sin(self.theta)*s_tri) + np.cos(self.theta+np.sin(self.theta)*np.arcsin(tri))
        dR_dtheta = np.gradient(self.R, edge_order=2)[1] * self.thetapts/(2*pi)  # -self.r*np.sin(self.theta+np.sin(self.theta)*np.arcsin(tri))*(1+np.cos(self.theta)*np.arcsin(tri))
    
        abs_grad_r = np.sqrt(dZ_dtheta**2 + dR_dtheta**2) / np.abs(dR_dr*dZ_dtheta - dR_dtheta*dZ_dr)
        
        # WE WANT TO CALCULATE THE POLOIDAL FIELD STRENGTH EVERYWHERE
        # THE PROBLEM IS THAT WE'VE GOT 2 EQUATIONS IN 3 UNKNOWNS. HOWEVER, IF WE ASSUME THAT THE POLOIDAL
        # INTEGRAL OF THE FLUX SURFACE AVERAGE OF THE POLOIDAL MAGNETIC FIELD IS APPROX. THE SAME AS THE
        # POLOIDAL INTEGRAL OF THE ACTUAL POLOIDAL MAGNETIC FIELD, THEN WE CAN CALCULATE THE Q PROFILE
        self.BT = p.BT0 * self.R[0, 0] / self.R
        
        # Calculate initial crappy guess on q
        q_mil = p.BT0*self.R[0, 0] / (2*pi*B_p_bar) * np.tile(np.sum(self.L_seg/self.R**2, axis=1), (self.thetapts, 1)).T #Equation 16 in the miller paper. The last term is how I'm doing a flux surface average
        q_mil[0, :] = q_mil[1, :]
        
        dpsidr = (p.BT0 * self.R[0, 0]) / (2*pi*q_mil)*np.tile(np.sum(self.L_seg/(self.R*abs_grad_r), axis=1), (self.thetapts, 1)).T
        
        self.psi = np.zeros(self.r.shape)
        for index, row in enumerate(self.r):
            if index >= 1:
                self.psi[index] = dpsidr[index]*(self.r[index, 0]-self.r[index-1, 0]) + self.psi[index-1]
        self.psi_norm = self.psi / self.psi[-1, 0]
        
        self.B_p = dpsidr * 1/self.R * abs_grad_r
        self.B_p[0, :] = 0
        
        
        self.B_t = p.BT0 * self.R[0, 0] / self.R
        self.B_tot = np.sqrt(self.B_p**2 + self.B_t**2)
        self.f_phi = self.B_t/self.B_tot
        #######################################################################
        ## CALCULATE ELECTRIC POTENTIAL FROM EXPERIMENTAL RADIAL ELECTRIC FIELD DATA
        

    def xmiller(self, kappa, p):
        """
        """
        ##########################################################################
        ## PART 1: SELECT CONVOLUTION KERNEL. WE'LL USE A STANDARD BUMP FUNCTION.
        ##########################################################################
        def bump(x, epsilon):
            """
            """
            #define eta0
            def eta0(x2):
                """
                """
                eta0 = np.where(np.logical_and((x2>-1), (x2<1)), np.exp(-1.0/(1.0-np.power(np.abs(x2), 2.0))), 0.)
                return eta0
            #calculate normalization coefficient
            C = quad(eta0, -1, 1)[0]
            #calculate eta_eps
            eta_eps = 1/epsilon/C*eta0(x/epsilon)
            return eta_eps  
            
        #############################################################            
        ## PART 2: DEFINE SEPERATRIX FUNCTION
        #############################################################          
        def kappa_sep(x, gamma1, gamma2):
            """
            """
            #Necessary kappa to get the z-value of the x-point
            kappa_x = (p.xpt[1]-p.Z0) / (p.a*sin(3.*pi/2.))
            # Amount of "extra" kappa needed at the x-point, i.e. kappa_tilda
            delta_kappa = kappa_x - p.kappa_lo
    
            kappa_sep = np.piecewise(
                                x, 
                                [x <= pi, #condition 1
                                 np.logical_and((pi < x), (x <= 3.*pi/2.)), #condition 2
                                 np.logical_and(( 3.*pi/2. < x), (x <= 2.*pi)), #condition 3
                                 x>2*pi], #condition 4
                                 [lambda x: 0, #function for condition 1
                                  lambda x: delta_kappa*(1.-abs(2.*x/pi - 3.)**(1.0))**gamma1, #function for condition 2
                                  lambda x: delta_kappa*(1.-abs(2.*x/pi - 3.)**(1.0))**gamma2, #function for condition 3
                                  lambda x: 0] #function for condition 4
                                 )
            return kappa_sep  
        
        ##########################################################################
        ## PART 2: RESCALE THETA (OPTIONAL) (HOW QUICKLY TO DO YOU LIMIT
        ## THETA RANGE OF CONVLUTION, IF AT ALL)
        ##########################################################################
        def rescale_theta(theta, epsilon):
            """
            """
            return (theta-3*pi/2)/(1-epsilon)+3*pi/2
            
        ##########################################################################
        ## PART 3: DEFINE Y-SCALE (HOW QUICKLY DO FLUX SURFACES
        ## APPROACH THE SEPERATRIX AS A FUNCTION OF r)
        ##########################################################################
        def yscale(r, a):
            """
            """
            nu=5.0 #nu=10ish works well
            return np.power(r/a, nu) #* (xpt[1] / (a*sin(3.*pi/2.)) - kappa_lo)
            
        ##########################################################################
        ## PART 4: DEFINE EPSILON AS A RUNCTION OF r (HOW QUICKLY DO FLUX SURFACES
        ## GET 'POINTY' AS THEY APPROACH THE SEPERATRIX)
        ##########################################################################
        def calc_eps(r, a):
            """
            """
            D=2.0 #What is this?
            nu=3.0
            return D*(1-np.power(r/a, nu))
        ##########################################################################
        ## PART 5: POST TRANSFORM (SUBTRACT ENDPOINTS, ETC.)
        ##########################################################################
        def posttransform(x, f): #here x is from pi to 2pi
            """
            """
            f_pt = f - ((f[-1]-f[0])/(2*pi-pi) * (x-pi) + f[0])
            return f_pt
        
        #INITIALIZE KAPPA_TILDA ARRAY
        kappa_tilda = np.zeros(kappa.shape)
    
        #DEFINE POINTS FOR THE CONVOLUTION. THIS IS SEPARATE FROM THETA POINTS.
        xnum = 1001
        
        #DEFINE SEPERATRIX FUNCTION
        gamma1 = 5.0 #3.8
        gamma2 = 2.0 #1.0
        k_sep = kappa_sep(np.linspace(pi, 2*pi, xnum), gamma1, gamma2)    
    
        #For each flux surface (i.e. r value)
        for i, rhoval in enumerate(self.rho.T[0]):
            if rhoval < 1:
                #Calculate epsilon
                epsilon = calc_eps(rhoval, p.a)
                #OPTIONALLY - modify domain to ensure a constant domain of compact support based on current epsilon
                scale_theta=0
                if scale_theta==1:
                    thetamod = rescale_theta(self.theta, epsilon)
                
                #define convolution kernel for the flux surface, eta_epsilon (eta for short)
                eta = bump(np.linspace(-1, 1, xnum), epsilon)
                
                #scale eta. The convolution operation doesn't 
                scaled_eta = eta * yscale(rhoval, p.a)
                
                #convolve seperatrix function and bump function
                kappa_tilda_pre = np.convolve(k_sep, scaled_eta, 'same')/((xnum-1)/2) #Still don't understand why we need to divide by this, but we definitely need to.
        
                #post processing
                kappa_tilda_post = posttransform(np.linspace(pi, 2*pi, xnum), kappa_tilda_pre)
                
                #create a 1D interpolation function for everywhere except for the seperatrix
                kappa_tilda_f = interp1d(np.linspace(pi, 2*pi, xnum), kappa_tilda_post, kind="linear")
                for j in range(0, kappa_tilda.shape[1]):
                    if self.theta[i, j] > pi and self.theta[i, j] <= 2*pi:
                        kappa_tilda[i, j] = kappa_tilda_f(self.theta[i, j])
            else:
                kappa_tilda_f = interp1d(np.linspace(pi, 2*pi, xnum), k_sep, kind="linear")
                for j in range(0, kappa_tilda.shape[1]):
                    if self.theta[i, j] > pi and self.theta[i, j] <= 2*pi:
                        kappa_tilda[i, j] = kappa_tilda_f(self.theta[i, j])
        
        kappa = kappa + kappa_tilda
        return kappa
    
    def xsec(self, inp):
        #Fusion Reactivity calculation  
        def calc_sigv_fus(mode='dd'):   
            def sigv(Ti, mode): #function takes T in kev
                if mode=='dt':
                    B_G = 34.3827    
                    m_rc2 = 1124656
                    
                    C1 = 1.17302E-9
                    C2 = 1.51361E-2
                    C3 = 7.51886E-2
                    C4 = 4.60643E-3
                    C5 = 1.35000E-2
                    C6 = -1.06750E-4
                    C7 = 1.36600E-5
                
                    theta = Ti/(1.0-(Ti*(C2+Ti*(C4+Ti*C6)))/(1.0+Ti*(C3+Ti*(C5+Ti*C7))))
                    xi = (B_G**2.0/(4.0*theta))**(1.0/3.0)
                    sigv = C1 * theta * np.sqrt(xi/(m_rc2 * Ti**3.0)) * np.exp(-3.0*xi)
                    sigv = sigv/1.0E6 #convert from cm^3/s to m^3/s
                    
                elif mode=='dd':
                    
                    B_G = 31.3970 
                    m_rc2 = 937814
                    
                    #first for the D(d, p)T reaction
                    C1_1 = 5.65718E-12
                    C2_1 = 3.41267E-3
                    C3_1 = 1.99167E-3
                    C4_1 = 0.0
                    C5_1 = 1.05060E-5
                    C6_1 = 0.0
                    C7_1 = 0.0
                
                    theta_1 = Ti/(1.0-(Ti*(C2_1+Ti*(C4_1+Ti*C6_1)))/(1.0+Ti*(C3_1+Ti*(C5_1+Ti*C7_1))))
                    xi_1 = (B_G**2.0/(4.0*theta_1))**(1.0/3.0)
                    sigv_1 = C1_1 * theta_1 * np.sqrt(xi_1/(m_rc2 * Ti**3.0)) * np.exp(-3.0*xi_1)
                    
                    #then for the D(d, n)He3 reaction
                    
                    C1_2 = 5.43360E-12
                    C2_2 = 5.85778E-3
                    C3_2 = 7.68222E-3
                    C4_2 = 0.0
                    C5_2 = -2.96400E-6
                    C6_2 = 0.0
                    C7_2 = 0.0
                
                    theta_2 = Ti/(1.0-(Ti*(C2_2+Ti*(C4_2+Ti*C6_2)))/(1.0+Ti*(C3_2+Ti*(C5_2+Ti*C7_2))))
                    xi_2 = (B_G**2.0/(4.0*theta_2))**(1.0/3.0)
                    sigv_2 = C1_2 * theta_2 * np.sqrt(xi_2/(m_rc2 * Ti**3.0)) * np.exp(-3.0*xi_2)                
                    
                    sigv = (0.5*sigv_1 + 0.5*sigv_2) / 1.0E6 #convert from cm^3/s to m^3/s                
                return sigv
            
            #create logspace over the relevant temperature range
            #(bosch hale technically only valid over 0.2 - 100 kev)
            Ti_range            = np.logspace(-1, 2, 1000) #values in kev
            sigv_fus_range      = sigv(Ti_range, mode='dd') #in m^3/s
            sigv_fus_interp     = UnivariateSpline(Ti_range*1.0E3*1.6021E-19, sigv_fus_range, s=0) #converted to Joules
            self.sv_fus       = sigv_fus_interp(self.Ti_J)
            self.dsv_fus_dT   = sigv_fus_interp.derivative()(self.Ti_J)
            self.dsv_fus_dT_eq9   = sigv_fus_interp.derivative()(5.0E2*1.6021E-19)
        
        def calc_sigv_ion():
            #TODO: configure so it can use any of the cross section libraries
            #currently using the Stacey-Thomas cross sections
            T_exps_fit          = np.array([-1, 0, 1, 2, 3, 4, 5])
            sigv_exps_fit       = np.array([-2.8523E+01, -1.7745E+01, -1.3620E+01, 
                                            -1.3097E+01, -1.3301E+01, -1.3301E+01, -1.3301E+01])
            interp1             = UnivariateSpline(T_exps_fit, sigv_exps_fit, s=0)
            
            T_exps_range        = np.linspace(-1, 5, 1000)
            sigv_vals_range     = 10.0**interp1(T_exps_range) #in m^3/s
            
            T_vals_range        = np.logspace(-1, 5, 1000)*1.6021E-19 #in joules
            interp2             = UnivariateSpline(T_vals_range, sigv_vals_range, s=0)
            
            self.sv_ion       = interp2(self.Ti_J)
            self.dsv_ion_dT   = interp2.derivative()(self.Ti_J)

        def calc_svel():
        
            tint = np.array([-1, 0, 1, 2, 3])
            tnnt = np.array([0, 1, 2])

            elast = np.array([[-1.3569E+01, -1.3337E+01, -1.3036E+01, -1.3569E+01, -1.3337E+01], 
                              [-1.3036E+01, -1.3337E+01, -1.3167E+01, -1.3046E+01, -1.3036E+01], 
                              [-1.3046E+01, -1.2796E+01, -1.3036E+01, -1.3046E+01, -1.2796E+01]])
                
            interp1 = interp2d(tint, tnnt, elast)
            
            Ti_exps   = np.linspace(-1, 3, 100)
            Tn_exps   = np.linspace( 0, 2, 100)
            svel_vals = 10.0**(interp1(Ti_exps, Tn_exps)) #in m^3/s
            
            Ti_vals   = np.logspace(-1, 3, 100)*1.6021E-19 #in joules
            Tn_vals   = np.logspace( 0, 2, 100)*1.6021E-19 #in joules

            dsvel_dTi_vals  = np.gradient(svel_vals, Ti_vals, axis=0)

            Ti_vals2d, Tn_vals2d = np.meshgrid(Ti_vals, Tn_vals)
            
            Ti_mod = np.where(self.Ti_ev>1E3, 1E3 * 1.6021E-19, self.Ti_ev*1.6021E-19)
            Tn_mod = np.zeros(Ti_mod.shape) + 2.0*1.6021E-19

            self.sv_el       = griddata(np.column_stack((Ti_vals2d.flatten(), Tn_vals2d.flatten())), 
                                       svel_vals.flatten(), 
                                       (Ti_mod, Tn_mod), 
                                       method='linear', rescale=False)
            self.dsv_el_dT   = griddata(np.column_stack((Ti_vals2d.flatten(), Tn_vals2d.flatten())), 
                                       dsvel_dTi_vals.flatten(), 
                                       (Ti_mod, Tn_mod), 
                                       method='linear', rescale=False)
         
        def calc_svcxi_st():
            
            tint = np.array([-1, 0, 1, 2, 3])
            tnnt = np.array([0, 1, 2])
            
            cx = np.array([[-1.4097E+01, -1.3921E+01, -1.3553E+01, -1.4097E+01, -1.3921E+01], 
                           [-1.3553E+01, -1.3921E+01, -1.3824E+01, -1.3538E+01, -1.3553E+01], 
                           [-1.3538E+01, -1.3432E+01, -1.3553E+01, -1.3538E+01, -1.3432E+01]])
                
            interp1 = interp2d(tint, tnnt, cx)
            
            Ti_exps   = np.linspace(-1, 3, 100)
            Tn_exps   = np.linspace( 0, 2, 100)
            svcx_vals = 10.0**(interp1(Ti_exps, Tn_exps)) #in m^3/s
            
            Ti_vals   = np.logspace(-1, 3, 100)*1.6021E-19 #in joules
            Tn_vals   = np.logspace( 0, 2, 100)*1.6021E-19 #in joules

            dsvcx_dTi_vals  = np.gradient(svcx_vals, Ti_vals, axis=0)

            Ti_vals2d, Tn_vals2d = np.meshgrid(Ti_vals, Tn_vals)
            
            Ti_mod = np.where(self.Ti_ev>1E3, 1E3 * 1.6021E-19, self.Ti_ev*1.6021E-19)
            Tn_mod = np.zeros(Ti_mod.shape) + 2.0*1.6021E-19

            self.sv_cx       = griddata(np.column_stack((Ti_vals2d.flatten(), Tn_vals2d.flatten())), 
                                       svcx_vals.flatten(), 
                                       (Ti_mod, Tn_mod), 
                                       method='linear', rescale=False)
            
            self.dsv_cx_dT   = griddata(np.column_stack((Ti_vals2d.flatten(), Tn_vals2d.flatten())), 
                                       dsvcx_dTi_vals.flatten(), 
                                       (Ti_mod, Tn_mod), 
                                       method='linear', rescale=False)

        def calc_svrec_st(): 
            #TODO: check this calculation. -MH
            znint = np.array([16, 18, 20, 21, 22])
            Tint  = np.array([-1, 0, 1, 2, 3])
        
            rec = np.array([[-1.7523E+01, -1.6745E+01, -1.5155E+01, -1.4222E+01, -1.3301E+01], 
                            [-1.8409E+01, -1.8398E+01, -1.8398E+01, -1.7886E+01, -1.7000E+01], 
                            [-1.9398E+01, -1.9398E+01, -1.9398E+01, -1.9398E+01, -1.9398E+01], 
                            [-2.0155E+01, -2.0155E+01, -2.0155E+01, -2.0155E+01, -2.0155E+01], 
                            [-2.1000E+01, -2.1000E+01, -2.1000E+01, -2.1000E+01, -2.1000E+01]])
            
            interp1 = interp2d(znint, Tint, rec)
            
            zni_exps   = np.linspace(16, 22, 100)
            Ti_exps    = np.linspace(-1, 3, 100)
            svrec_vals = 10.0**(interp1(zni_exps, Ti_exps)) #in m^3/s
            
            zni_vals   = np.logspace(16, 22, 100)
            Ti_vals    = np.logspace(-1, 3, 100)*1.6021E-19 #in joules
            
            dsvrec_dTi_vals  = np.gradient(svrec_vals, Ti_vals, axis=0)
                
            zni_vals2d, Ti_vals2d = np.meshgrid(zni_vals, Ti_vals)
            
            zni_mod = np.where(self.ni>1E22, 1E22, self.ni)
            zni_mod = np.where(self.ni<1E16, 1E16, zni_mod)
            Ti_mod  = np.where(self.Ti_ev>1E3, 1E3 * 1.6021E-19, self.Ti_ev*1.6021E-19)
            Ti_mod  = np.where(self.Ti_ev<1E-1, 1E-1 * 1.6021E-19, Ti_mod)
                
            self.sv_rec     = griddata(np.column_stack((zni_vals2d.flatten(), Ti_vals2d.flatten())), 
                                       svrec_vals.flatten(), 
                                       (zni_mod, Ti_mod), 
                                       method='linear', rescale=False)
            
            self.dsv_rec_dT = griddata(np.column_stack((zni_vals2d.flatten(), Ti_vals2d.flatten())), 
                                       dsvrec_dTi_vals.flatten(), 
                                       (zni_mod, Ti_mod), 
                                       method='linear', rescale=False)

        calc_sigv_fus()
        calc_sigv_ion()
        calc_svel()
        calc_svcxi_st()
        #calc_svrec_st()
        