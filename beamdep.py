#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 18 21:51:19 2018

@author: max
"""
import numpy as np
from scipy.interpolate import interp1d
from subprocess import call
import os
import re
import sys
from scipy.interpolate import interp1d

class beamdep():
    """
    Methods:
    prep_nbi_infile
    read_nbi_outfile
    
    Attributes:
    P_abs_tot
    P_lst_tot
    I_nbi_tot
    I_nbi_eff
    Beta_nbi_tot
    Taus_nbi
    Volp_nbi
    P_abs_1
    P_abs_2
    P_abs_3
    P_lst_1
    P_lst_2
    P_lst_3
    I_nbi_1
    I_nbi_2
    I_nbi_3
    I_nbi_eff_1
    I_nbi_eff_2
    I_nbi_eff_3
    I_nbi_gam_1
    I_nbi_gam_2
    I_nbi_gan_3
    st_en1_1
    st_en1_2
    st_en1_3
    st_en2_1
    st_en2_2
    st_en2_3
    st_en3_1
    st_en3_2
    st_en3_3
    fus_pwr_bt
    cp_pwr_tot
    rate_dt_n
    rate_dd_n
    dep_prof1
    dep_prof2
    dep_prof3
    ptch_angl1
    ptch_angl2
    ptch_angl3
    jnbtot
    pNBe
    pNBi
    nbfast
    pressb
    pfusb
    dvol
    dA
    
    
    """
    
    def __init__(self,inp,brnd):
        sys.dont_write_bytecode = True 
        self.prep_nbi_infile(inp,brnd)
        #call nbeams. Note to those familiar with the old nbeams, I modified
        #the code to take the input file as a commandline argument.
        call([inp.nbeams_loc+'nbeams', os.getcwd()+'/inbeams.dat'])
        self.read_nbi_outfile(inp,brnd)
        pass
    
    def prep_nbi_infile(self,inp,brnd):
        #f1=open(inp.nbeams_loc+"inbeams.dat","w")
        f=open("inbeams.dat","w")
        f.write("$nbin\n")
        f.write("nbeams = 1\n")
        f.write("inbfus = 1\n")
        f.write("amb = 2.0\n")
        f.write("zbeam = 1.0\n")
        f.write("ebeam = "+str(inp.ebeam)+"\n")
        f.write("pbeam = "+str(inp.pbeam)+"\n")
        f.write("rtang = "+str(inp.rtang)+"\n")
        f.write("nbshape = 1\n")
        f.write("bwidth = 0.12\n")                                                   # Is this default?
        f.write("bheigh = 0.48\n")
        f.write("bgaussR = 0.066\n")
        f.write("bgaussZ = 0.18\n")
        f.write("bzpos = 0.0\n")
        f.write("nbptype = 1\n")
        f.write("maxiter = 2\n")
        f.write("pwrfrac(1,1) = "+str(inp.pwrfrac1)+"   "+str(inp.pwrfrac2)+"   "+str(inp.pwrfrac3)+"\n")
        f.write("a = "+str(inp.a)+"\n")
        f.write("r0 = "+str(inp.R0_a)+"\n")
        f.write("b0 = "+str(inp.bknot)+"\n")
        f.write("n = 51\n")
        f.write("e0 = "+str(inp.epsknot)+"\n")
        f.write("ea = "+str(inp.eps_sep)+"\n")
        f.write("shft0 = "+str(inp.shftknot)+"\n")
        f.write("nion = 2\n")
        f.write("aion = 2.0 12.0\n")
        f.write("zion = 1.0 6.0\n")
        
        rho_nbi = np.linspace(0,1,51)
        ni_nbi = interp1d(inp.ni_rho[:,0],inp.ni_rho[:,1])(rho_nbi)
        ne_nbi = interp1d(inp.ne_rho[:,0],inp.ne_rho[:,1])(rho_nbi)
        Ti_nbi = interp1d(inp.Ti_rho[:,0],inp.Ti_rho[:,1])(rho_nbi)
        Te_nbi = interp1d(inp.Te_rho[:,0],inp.Te_rho[:,1])(rho_nbi)
        
        for i,v in enumerate(rho_nbi):
            f.write('ni20('+str(i+1)+',1) = '+str(ni_nbi[i]*1E-20)+'\n')
        for i,v in enumerate(rho_nbi):
            f.write('ne20('+str(i+1)+') = '+str(ne_nbi[i]*1E-20)+'\n')
        for i,v in enumerate(rho_nbi):
            f.write('tikev('+str(i+1)+') = '+str(Ti_nbi[i])+'\n')
        for i,v in enumerate(rho_nbi):
            f.write('tekev('+str(i+1)+') = '+str(Te_nbi[i])+'\n')
        f.write("$end\n")
        f.close()
    
    def read_nbi_outfile(self,inp,brnd):
        print 'reading nbi outfile'
        with open(os.getcwd() + '/outbeams.dat', 'r') as f:
            for count, line in enumerate(f):
                #print line
                if line.startswith(" Total Absorbed Power"):
                    result = re.match(r'.*= *((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*',line).group(1)
                    try:
                        self.P_abs_tot = float(result)
                    except:
                        self.P_abs_tot = np.NaN

                if line.startswith(" Total Lost Power"):
                    result = re.match(r'.*= *((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*',line).group(1)
                    try:
                        self.P_lst_tot = float(result)
                    except:
                        self.P_lst_tot = np.NaN

                if line.startswith(" Total NB Driven Current"):
                    result = re.match(r'.*= *((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*',line).group(1)
                    try:
                        self.I_nbi_tot = float(result)
                    except:
                        self.I_nbi_tot = np.NaN
                        
                if line.startswith(" Total NBCD Efficiency"):
                    result = re.match(r'.*= *((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*',line).group(1)
                    try:
                        self.I_nbi_eff = float(result)
                    except:
                        self.I_nbi_eff = np.NaN
                        

                if line.startswith(" Total Beam Beta"):
                    result = re.match(r'.*= *((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*',line).group(1)
                    try:
                        self.Beta_nbi_tot = float(result)
                    except:
                        self.Beta_nbi_tot = np.NaN
                        
                if line.startswith(" Taus"):
                    result = re.match(r'.*= *((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*',line).group(1)
                    try:
                        self.Taus_nbi = float(result)
                    except:
                        self.Taus_nbi = np.NaN
                        
                if line.startswith(" Volp"):
                    result = re.match(r'.*= *((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*',line).group(1)
                    try:
                        self.Volp_nbi = float(result)
                    except:
                        self.Volp_nbi = np.NaN

                if line.startswith(" Absorbed Power"):
                    #this will need to be modified for multiple beams
                    print 'line = ',line
                    result = re.match(r'.*((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*',line).group(1)
                    try:
                        self.P_abs_1 = float(result)
                    except:
                        self.P_abs_1 = np.NaN
                        
                if line.startswith(" Lost Power"):                   
                    #this will need to be modified for multiple beams
                    
                    result = re.match(r'.*((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*',line).group(1)
                    try:
                        self.P_lst_1 = float(result)
                    except:
                        self.P_lst_1 = np.NaN
                        
                if line.startswith(" NB driven current"):
                    #this will need to be modified for multiple beams
                    result = re.match(r'.*((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*',line).group(1)
                    try:
                        self.I_nbi_1 = float(result)
                    except:
                        self.I_nbi_1 = np.NaN
                        
                if line.startswith(" NBCD efficiency"):
                    #this will need to be modified for multiple beams
                    result = re.match(r'.*((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*',line).group(1)
                    try:
                        self.I_nbi_eff_1 = float(result)
                    except:
                        self.I_nbi_eff_1 = np.NaN
                        
                if line.startswith(" NBCD gamma"):
                    #this will need to be modified for multiple beams
                    result = re.match(r'.*((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*',line).group(1)
                    try:
                        self.I_nbi_gam_1 = float(result)
                    except:
                        self.I_nbi_gam_1 = np.NaN
                        
                if line.startswith("    energy group 1"):
                    #this will need to be modified for multiple beams
                    result = re.match(r'.*((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*',line).group(1)
                    try:
                        self.st_en1_1 = float(result)
                    except:
                        self.st_en1_1 = np.NaN
                        
                if line.startswith("    energy group 2"):
                    #this will need to be modified for multiple beams
                    result = re.match(r'.*((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*',line).group(1)
                    try:
                        self.st_en2_1 = float(result)
                    except:
                        self.st_en2_1 = np.NaN
                        
                if line.startswith("    energy group 3"):
                    #this will need to be modified for multiple beams
                    result = re.match(r'.*((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*',line).group(1)
                    try:
                        self.st_en3_1 = float(result)
                    except:
                        self.st_en3_1 = np.NaN
                        
                if line.startswith(" Total Beam-Target Fusion Power"):
                    #this will need to be modified for multiple beams
                    result = re.match(r'.*= *((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*',line).group(1)
                    try:
                        self.fus_pwr_bt = float(result)
                    except:
                        self.fus_pwr_bt = np.NaN
                        
                if line.startswith(" Total Power to Charged Particles"):
                    #this will need to be modified for multiple beams
                    result = re.match(r'.*= *((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*',line).group(1)
                    try:
                        self.cp_pwr_tot = float(result)
                    except:
                        self.cp_pwr_tot = np.NaN
                        
                if line.startswith(" Total DT Neutron Rate"):
                    #this will need to be modified for multiple beams
                    result = re.match(r'.*= *((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*',line).group(1)
                    try:
                        self.rate_dt_n = float(result)
                    except:
                        self.rate_dt_n = np.NaN
                        
                if line.startswith(" Total DD Neutron Rate"):
                    #this will need to be modified for multiple beams
                    result = re.match(r'.*= *((?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)?)|NaN).*',line).group(1)
                    try:
                        self.rate_dd_n = float(result)
                    except:
                        self.rate_dd_n = np.NaN
                        
        with open(os.getcwd() + '/outbeams.dat', 'r') as f:
            data = f.read().replace('\n',' ')
            #print data
            result = re.match(r'.*hofr_3 *((?:(?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)? +)|(?:NaN +))+) +Pitch.*',data).group(1)
            array = np.reshape(np.asarray(result.split(),dtype=float),(-1,4))
            self.dep_prof1  = interp1d(array[:,0],array[:,1])(brnd.rho)
            self.dep_prof2  = interp1d(array[:,0],array[:,2])(brnd.rho)
            self.dep_prof3  = interp1d(array[:,0],array[:,3])(brnd.rho)

            result = re.match(r'.*zeta_3 *((?:(?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)? +)|(?:NaN +))+) +rho.*',data).group(1)
            array = np.reshape(np.asarray(result.split(),dtype=float),(-1,4))
            self.ptch_angl1 = interp1d(array[:,0],array[:,1])(brnd.rho)
            self.ptch_angl2 = interp1d(array[:,0],array[:,2])(brnd.rho)
            self.ptch_angl3 = interp1d(array[:,0],array[:,3])(brnd.rho)
            
            result = re.match(r'.*dA *((?:(?:[-\+]?\d*(?:.?\d+)?(?:[Ee][-\+]?\d+)? +)|(?:NaN +))+) *',data).group(1)
            array = np.reshape(np.asarray(result.split(),dtype=float),(-1,9))
            self.jnbtot     = interp1d(array[:,0],array[:,1])(brnd.rho)
            self.pNBe_dvol       = interp1d(array[:,0],array[:,2] * array[:,7])(brnd.rho)
            self.pNBi_dvol       = interp1d(array[:,0],array[:,3] * array[:,7])(brnd.rho)
            self.nbfast     = interp1d(array[:,0],array[:,4])(brnd.rho)
            self.pressb     = interp1d(array[:,0],array[:,5])(brnd.rho)
            self.pfusb      = interp1d(array[:,0],array[:,6])(brnd.rho)
            self.dvol       = interp1d(array[:,0],array[:,7])(brnd.rho)
            self.dA         = interp1d(array[:,0],array[:,8])(brnd.rho)