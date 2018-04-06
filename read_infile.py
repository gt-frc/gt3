# -*- coding: utf-8 -*-
"""
Created on Sat Aug  5 16:10:25 2017

@author: max
"""
from math import ceil
import os
import re
import sys
import numpy as np
from math import pi
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from shapely.geometry import LineString, Polygon
from scipy.interpolate import UnivariateSpline, interp1d

class read_infile():
    """Reads main GT3 input file.
    
    Methods:
        read_vars
        read_exp
        wall_prep
        showparams
    
    Attributes:
        d3d_iter         (int)      
        d3d_input        (int)      
        a                (float)    tokamak minor radius (m)
        B_phi_0          (float)    toroidal field strength at mag. axis (T)
        R0_a             (float)    tokamak major radius (m)
        Z0               (float)    vertical height of the magnetic axis (m)
        kappa_up         (float)    upper elongation at the seperatrix
        kappa_lo         (float)    lower elongation at the seperatrix
        tri_up           (float)
        tri_lo           (float)
        xmil             (int)
        xpt_R            (float)
        xpt_Z            (float)
        xpt
        thetapts_approx  (int)
        thetapts
        rmeshnum_p       (int)
        rpts             (int)
        ni0              (float)
        ni9              (float)
        ni_sep           (float)
        nu_ni            (float)
        ne0              (float)
        ne9              (float)
        ne_sep           (float)
        nu_ne            (float)  
        Ti0              (float)
        Ti9              (float)
        Ti_sep           (float)
        nu_Ti            (float)
        Te0              (float)
        Te9              (float)
        Te_sep           (float)
        nu_Te            (float)
        j0               (float)
        j_sep            (float)
        nu_j             (float)
        s_k_up           (float)
        s_k_lo           (float)
        xtheta1          (float)
        xtheta2          (float)
        xtheta3          (float)
        xtheta4          (float)
        wallfile         (str)         
    """
    
    def __init__(self, infile):
        sys.dont_write_bytecode = True 
        self.read_vars(infile)
        self.read_exp()
        if hasattr(self, 'lim_file'):
            self.wall_prep()

    def read_vars(self, infile):
        """
        """
        
        #some regex commands we'll use when reading stuff in from the input file
        r0di = "r'%s *= *([ ,\d]*) *'%(v)"
        r0df = "r'%s *= *([+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?) *'%(v)"
        r0ds = "r'%s *= *((?:/?\.?\w+\.?)+/?) *'%(v)"
        r1di = "r'%s\( *(\d*) *\) *= *(\d*) *'%(v)"
        r1df = "r'%s\( *(\d*)\) *= *([+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?) *'%(v)"
        r2df = "r'%s\( *(\d*)\) *= *((?:[+\-]?\d*\.?\d*(?:[eE]?[+\-]?\d+)?,?)*) *'%(v)"

        self.invars = {}
        self.invars["d3d_iter"]         = ["int",r0di]
        self.invars["d3d_input"]        = ["int",r0di]
        self.invars["a"]                = ["float",r0df]
        self.invars["B_phi_0"]          = ["float",r0df]
        self.invars["R0_a"]             = ["float",r0df]
        self.invars["Z0"]               = ["float",r0df]
        self.invars["kappa_up"]         = ["float",r0df]
        self.invars["kappa_lo"]         = ["float",r0df]
        self.invars["tri_up"]           = ["float",r0df]
        self.invars["tri_lo"]           = ["float",r0df]
        self.invars["xmil"]             = ["int",r0di]
        self.invars["xpt_R"]            = ["float",r0df]
        self.invars["xpt_Z"]            = ["float",r0df]
        self.invars["thetapts_approx"]  = ["int",r0di]
        self.invars["rmeshnum_p"]       = ["int",r0di]
        self.invars["rpts"]             = ["int",r0di]
        self.invars["ni0"]              = ["float",r0df]
        self.invars["ni9"]              = ["float",r0df]
        self.invars["ni_sep"]           = ["float",r0df]
        self.invars["nu_ni"]            = ["float",r0df]
        self.invars["ne0"]              = ["float",r0df]
        self.invars["ne9"]              = ["float",r0df]
        self.invars["ne_sep"]           = ["float",r0df]
        self.invars["nu_ne"]            = ["float",r0df]  
        self.invars["Ti0"]              = ["float",r0df]
        self.invars["Ti9"]              = ["float",r0df]
        self.invars["Ti_sep"]           = ["float",r0df]
        self.invars["nu_Ti"]            = ["float",r0df]
        self.invars["Te0"]              = ["float",r0df]
        self.invars["Te9"]              = ["float",r0df]
        self.invars["Te_sep"]           = ["float",r0df]
        self.invars["nu_Te"]            = ["float",r0df]
        self.invars["j0"]               = ["float",r0df]
        self.invars["j_sep"]            = ["float",r0df]
        self.invars["nu_j"]             = ["float",r0df]
        self.invars["q95"]              = ["float",r0df]
        self.invars["s_k_up"]           = ["float",r0df]
        self.invars["s_k_lo"]           = ["float",r0df]
        self.invars["xtheta1"]          = ["float",r0df]
        self.invars["xtheta2"]          = ["float",r0df]
        self.invars["xtheta3"]          = ["float",r0df]
        self.invars["xtheta4"]          = ["float",r0df]
     
        #NEW FROM CONSTS
        self.invars["eq1"]              = ["float",r0df]
        self.invars["eq2"]              = ["float",r0df]
        self.invars["xmas1"]            = ["float",r0df]
        self.invars["xmas2"]            = ["float",r0df]
        self.invars["ephia"]            = ["float",r0df]
        self.invars["xk"]               = ["float",r0df]
        self.invars["delma"]            = ["float",r0df]
        self.invars["xnuati"]           = ["float",r0df]
        self.invars["xnuioni"]          = ["float",r0df]
        
        #NEUTRALS CALCULATION
        self.invars["neut_outfile"]     = ["str",r0ds]        
        self.invars["ntrl_rho_start"]   = ["str",r0ds]        

        #NEUTRAL BEAM CALCULATION
        self.invars["nbeams_loc"]       = ["str",r0ds]        
        self.invars["adpak_loc"]        = ["str",r0ds]        
        self.invars["ebeam"]            = ["float",r0df]
        self.invars["abeam"]            = ["float",r0df]
        self.invars["alphain"]          = ["float",r0df]
        self.invars["pbeam"]            = ["float",r0df]
        self.invars["rtang"]            = ["float",r0df]
        self.invars["bknot"]            = ["float",r0df]
        self.invars["pwrfrac1"]         = ["float",r0df]
        self.invars["pwrfrac2"]         = ["float",r0df]
        self.invars["pwrfrac3"]         = ["float",r0df]
        self.invars["epsknot"]          = ["float",r0df]
        self.invars["eps_sep"]          = ["float",r0df]
        self.invars["shftknot"]         = ["float",r0df]
    
        self.in_prof = {}
        self.in_prof["er_file"]         = ["str",r0ds,'er_rho']
        self.in_prof["jr_file"]         = ["str",r0ds,'jr_rho']
        self.in_prof["ne_file"]         = ["str",r0ds,'ne_rho']
        self.in_prof["ni_file"]         = ["str",r0ds,'ni_rho']
        self.in_prof["Te_file"]         = ["str",r0ds,'Te_rho']
        self.in_prof["Ti_file"]         = ["str",r0ds,'Ti_rho']
        self.in_prof["fz1_file"]        = ["str",r0ds,'fz1_rho']
        self.in_prof["fracz_file"]      = ["str",r0ds,'fracz_rho']
        self.in_prof["exlti_file"]      = ["str",r0ds,'exlti_rho']
        self.in_prof["exlte_file"]      = ["str",r0ds,'exlte_rho']
        self.in_prof["exlni_file"]      = ["str",r0ds,'exlni_rho']
        self.in_prof["vpolC_file"]      = ["str",r0ds,'vpolC_rho']
        self.in_prof["vtorC_file"]      = ["str",r0ds,'vtorC_rho']
        self.in_prof["vpolD_file"]      = ["str",r0ds,'vpolD_rho']
        self.in_prof["vtorD_file"]      = ["str",r0ds,'vtorD_rho']
        self.in_prof["q_file"]          = ["str",r0ds,'q_rho']
        self.in_prof["zbar2_file"]      = ["str",r0ds,'zbar2_rho']

        self.in_map2d = {}
        self.in_map2d["psirz_file"]     = ["str",r0ds,'psirz_exp']
        self.in_map2d["bpol_file"]      = ["str",r0ds,'bpol_exp']
        self.in_map2d["btor_file"]      = ["str",r0ds,'btor_exp']
        
        self.in_line2d = {}
        self.in_line2d["lim_file"]      = ["str",r0ds,'lim_exp']
        self.in_line2d["sep_file"]      = ["str",r0ds,'sep_exp']

        with open(os.getcwd() + '/' + infile, 'r') as f:
            for count, line in enumerate(f):
                if not line.startswith("#"):
                    #read in 0d variables
                    for v in self.invars:
                        exec("result = re.match(%s,line)"%(self.invars[v][1]))
                        if result:
                            exec("self.%s = %s(result.group(1))"%(v,self.invars[v][0]))
                            
                    #read in the names of radial profile input files 
                    for v in self.in_prof:
                        exec("result = re.match(%s,line)"%(self.in_prof[v][1]))
                        if result:
                            exec("self.%s = %s(result.group(1))"%(v,self.in_prof[v][0]))
   
                    #read in the names of input files that map a quantity on the R-Z plane
                    for v in self.in_map2d:
                        exec("result = re.match(%s,line)"%(self.in_map2d[v][1]))
                        if result:
                            exec("self.%s = %s(result.group(1))"%(v,self.in_map2d[v][0]))

                    #read in the names of input files that define a line in the R-Z plane 
                    for v in self.in_line2d:
                        exec("result = re.match(%s,line)"%(self.in_line2d[v][1]))
                        if result:
                            exec("self.%s = %s(result.group(1))"%(v,self.in_line2d[v][0])) 

        self.thetapts =  int(4 * ceil(float(self.thetapts_approx)/4))+1 
        self.xpt = np.array([self.xpt_R,self.xpt_Z])

    def read_exp(self):
        
        #read in additional input files 
        for infile in self.in_prof:
            try:
                exec("filename = self.%s"%(infile))
                filepath = os.getcwd()+'/inputs/'+ filename
                exec("self.%s = np.loadtxt('%s')"%(self.in_prof[infile][2],filepath))
            except:
                pass
            
        for infile in self.in_map2d:
            try:

                exec("filename = self.%s"%(infile))
                filepath = os.getcwd()+'/inputs/'+ filename
                exec("self.%s = np.loadtxt('%s')"%(self.in_map2d[infile][2],filepath))
            except:
                pass
            
        for infile in self.in_line2d:
            try:
                exec("filename = self.%s"%(infile))
                filepath = os.getcwd()+'/inputs/'+ filename
                exec("self.%s = np.loadtxt('%s')"%(self.in_line2d[infile][2],filepath))        
            except:
                pass
            
    def wall_prep(self):
        """
        """   
    
        adotb = (self.lim_exp[:,0]-np.roll(self.lim_exp[:,0],1))*(self.lim_exp[:,0]-np.roll(self.lim_exp[:,0],-1)) + \
                (self.lim_exp[:,1]-np.roll(self.lim_exp[:,1],1))*(self.lim_exp[:,1]-np.roll(self.lim_exp[:,1],-1))
        mag_a = np.sqrt((self.lim_exp[:,0]-np.roll(self.lim_exp[:,0], 1))**2+(self.lim_exp[:,1]-np.roll(self.lim_exp[:,1], 1))**2)
        mag_b = np.sqrt((self.lim_exp[:,0]-np.roll(self.lim_exp[:,0],-1))**2+(self.lim_exp[:,1]-np.roll(self.lim_exp[:,1],-1))**2)
        
        lim_angles = np.arccos(adotb/(mag_a*mag_b))/pi
        self.lim_vertex = np.zeros((self.lim_exp.shape[0],2))
        
        for i in range(0,self.lim_exp.shape[0]):
            if lim_angles[i] <= 0.99:
                self.lim_vertex[i,:]=self.lim_exp[i,:]
            else:
                self.lim_vertex[i,:]=0
        self.lim_vertex = self.lim_vertex[np.all(self.lim_vertex != 0, axis=1)] #removing zeros from array
        #need to add in an additional criteria to also remove points that are extremely close to other points, even if they create a sufficiently large angle
        self.lim_line = LineString(np.vstack((self.lim_vertex,self.lim_vertex[0])))
        #self.lim_line = LineString(self.lim_vertex)
 

    def showparams(self):
        """
        """
        print '**PARAMETERS FOR SHOT \'{}\'.'.format(self.shotlabel)
        for key in vars(self).iteritems():
            if key[0][0]!='_' and key[0]!='line' and key[0]!='infile' and key[0]!='variable' and key[0]!='value':
                print ('{} = {}'.format(key[0],key[1]))
        print '**END OF PARAMETERS**'