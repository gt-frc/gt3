#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 24 13:25:51 2018

@author: max
"""
from subprocess import call
import numpy as np
from scipy.interpolate import interp2d, interp1d, UnivariateSpline
import os
import re
import sys
import matplotlib.pyplot as plt

class imp_rad():
    """
    Attributes:
        
    Methods:
        
    """
    def __init__(self,inp,brnd):
        sys.dont_write_bytecode = True 
        self.prep_adpak_infile(inp,brnd)
        
        try:
            call(['adpak', os.getcwd()+'/toadpak'])
        except:
            try:
                call([inp.adpak_loc+'adpak', os.getcwd()+'/toadpak'])
            except AttributeError:
                print 'Unable to locate adpak. Stopping.
                sys.exit()

        self.read_adpak_outfile(inp,brnd)
        self.coronal_eq(inp,brnd)
        pass
    
    def prep_adpak_infile(self,inp,brnd):
        
        #TODO: Convert adpack routines to python binary using f2py so we can
        #eliminate the main.f driver program and do what we need to in python.
        
        #TODO: Extend this to do elements other than carbon
        self.inucz = 6
        self.nte = 25
        self.nne = 2
        self.nmin = -3
        self.nmax = 2
        self.tei = np.logspace(self.nmin,self.nmax,self.nte) #actual temperatures
        self.tei_lin = np.linspace(self.nmin,self.nmax,self.nte)
        self.anei = np.array([1.0E13, 1.0E+15])
        f = open('./toadpak','w')
        f.write(' &inp')
        f.write('\n' + '  inucz = ' + str(self.inucz))
        f.write('\n' + '  zte = 1.0') #this parameter doesn't matter because imode = 1
        f.write('\n' + '  zne = 1.0e14') #this parameter doesn't matter because imode = 1
        f.write('\n' + '  laden = 1')
        f.write('\n' + '  ladtip = 1')
        f.write('\n' + '  leci = 1')
        f.write('\n' + '  ldrmlt = 2')
        f.write('\n' + '  ncxb = 0')
        f.write('\n' + '  ncxopt = 1')
        f.write('\n' + '  ivunit = 2')
        f.write('\n' + '  anneut = 1.0e18')
        f.write('\n' + '  vneut = 0.001')
        f.write('\n' + '  imode = 1')
        f.write('\n' + '  nte = ' + str(self.nte))
        f.write('\n' + '  nne = ' + str(self.nne))
        f.write('\n' + '  tei = ' + ' '.join(map(str,self.tei )))
        f.write('\n' + '  anei = '+ ' '.join(map(str,self.anei)))
        f.write('\n' + '  nmin = ' + str(self.nmin))
        f.write('\n' + '  nmax = ' + str(self.nmax))
        f.write('\n' + ' $end')
        f.write('\n' + '')
        f.write('\n')
        f.close()
    
    def read_adpak_outfile(self,inp,brnd):
        
        #some regex commands we'll use when reading stuff in from the input file
        regex1 = "r'.*?data \(+%s.*?\/(.*?)\/.*'%(v)"
        regex2 = "r'.*?data \(+%s\(\d+ +%s.*?\/(.*?)\/.*'%(v,cs+1)"
        self.v1d = {}
        self.v1d["altei"]  = [regex1]
        self.v1d["alnei"]  = [regex1]
        
        self.v3d = {}
        self.v3d["alinzr"] = [regex2]
        self.v3d["alradr"] = [regex2]
        self.v3d["alrecr"] = [regex2]
        
        self.alinzr = np.zeros((self.inucz+1,self.nne,self.nte))
        self.alradr = np.zeros((self.inucz+1,self.nne,self.nte))
        self.alrecr = np.zeros((self.inucz+1,self.nne,self.nte))
        
        with open(os.getcwd() + '/outblk.dat', 'r') as f:
            data = f.read().replace('\n',' ').replace('. ',' ').replace(',',' ')
            for v in self.v1d:
                exec("result = re.match(%s,data).group(1)"%(self.v1d[v][0]))
                exec("self.%s = np.asarray(result.split(),dtype=float)"%(v))                
            for v in self.v3d:
                for cs in np.arange(self.inucz+1): #for each charge state plus the ground state,
                    exec("result = re.match(%s,data).group(1)"%(self.v3d[v][0]))
                    exec("self.%s[cs,:,:] = np.asarray(result.split(),dtype=float).reshape(-1,self.nte)"%(v))

    def coronal_eq(self,inp,brnd):
        def null(M, eps=1e-20):
            u, s, vh = np.linalg.svd(M)
            return vh[np.argmin(np.abs(s))]
        
        #Create nchrgsr x nchrgsr array
        M = np.zeros((self.inucz+1,self.inucz+1))

        frac_abun = np.zeros((self.nne,self.nte,self.inucz+1))
        for i,dens in enumerate([1.0E13, 1.0E+15]):
            for j,temp in enumerate(np.logspace(self.nmin,self.nmax,self.nte)):
    
                for (cs,col),value in np.ndenumerate(M):
                    if cs==0: #do the first row differently
                        M[cs,cs  ] = -(10**self.alinzr[cs  ,i,j] + 10**self.alrecr[cs  ,i,j]) 
                        M[cs,cs+1] =                               10**self.alrecr[cs+1,i,j]
                    elif cs==self.inucz: #do the last row differently
                        M[cs,cs-1] =   10**self.alinzr[cs-1,i,j]
                        M[cs,cs  ] = -(10**self.alinzr[cs  ,i,j] + 10**self.alrecr[cs  ,i,j]) 
                    else:
                        M[cs,cs-1] =   10**self.alinzr[cs-1,i,j]
                        M[cs,cs  ] = -(10**self.alinzr[cs  ,i,j] + 10**self.alrecr[cs  ,i,j]) 
                        M[cs,cs+1] =                               10**self.alrecr[cs+1,i,j]
                
                #the following line may help convergence. Uncomment it if necessary. Everything
                #gets normalized anyway.
                #M = M*1E10
    
                result = np.abs(null(M).T) #solve(M,b)
                frac_abun[i,j,:] = result / np.sum(result) #solve(M,b)

        #CALCULATE IMPURITY DENSITIES BY CHARGE STATE
        self.brnd_nC_cs     = []
        self.brnd_rad_cs    = []
        emiss_tot = np.zeros(self.tei.shape)
        #loop over charge states
        for cs in np.arange(self.inucz+1):
            #create the interpolation function to get fractional abundances of each 
            #charge state as a function of electron density and electron temperature
            frac_abun_interp = interp2d(self.tei_lin,self.anei*1E6,frac_abun[:,:,cs]) 
            
            #for each charge state, get the fractional abundance for the charge state
            # under consideration at each point in the plasma
            brnd_frac_abun  = np.zeros(brnd.nC.shape)
            for (r,theta),nC in np.ndenumerate(brnd.nC):
                brnd_frac_abun[r,theta] = frac_abun_interp( np.log10(brnd.Te_kev[r,theta]) , brnd.ne[r,theta] )
            
            #multiply the fractional abundance by the overal impurity density to get the
            #density of this charge state
            self.brnd_nC_cs.append( brnd_frac_abun * brnd.nC )

            #create the emissivity function over the entire range of temperatures
            emiss_cs_interp = interp2d(self.tei_lin,self.anei*1E6,self.alradr[cs,:,:])
            for i,T in enumerate(self.tei_lin):
                emiss_tot[i] = emiss_tot[i] + \
                                    10**emiss_cs_interp(T,self.anei[0]*1E6) * \
                                    frac_abun_interp(T,self.anei[0]*1E6)
                                    #note that the density doesn't actuall matter here.
        
        #spline fit the logarithm of the emissivity. We'll spline fit the actual values next
        emiss_tot_interp = UnivariateSpline(self.tei,np.log10(emiss_tot)+100.0,s=0)
        
        #the above function gives the right values, but we also need the derivatives
        #we will now do a spline fit of the values themselves, rather than their base-10 logarithm
        new_T_kev               = np.logspace(-3,2,10000)
        new_Em                  = 10.0**(emiss_tot_interp(new_T_kev)-100.0)
        
        new_T_J                 = new_T_kev * 1.0E3 * 1.6021E-19
        self.emiss_tot_interp2  = UnivariateSpline(new_T_J,new_Em,s=0)
        
        self.brnd_emissivity    = self.emiss_tot_interp2(brnd.Te_J)
        self.brnd_dEmiss_dT     = self.emiss_tot_interp2.derivative()(brnd.Te_J)
        self.brnd_dEmiss_dT_eq9     = self.emiss_tot_interp2.derivative()(5.0E2*1.6021E-19)
        self.brnd_dEmiss_dT_eq22     = self.emiss_tot_interp2.derivative()(21.0*1.6021E-19)
        print
        print '#######################'
        print 'Emiss_marfe = ',self.emiss_tot_interp2(250.0*1.6021E-19)
        print 'dEmiss_dT_marfe = ',self.emiss_tot_interp2.derivative()(250.0*1.6021E-19)
        print '#######################'
        print
        emiss_fig = plt.figure(figsize=(6,4))
        ax1 = emiss_fig.add_subplot(1,1,1)
        #ax1.set_xlim(4.5,5.5)
        #ax1.set_yscale('symlog')
        #ax1.set_xscale('log')
        #ax1.set_ylim(-1E-33,1E-33)
        #ax1.loglog(new_T_kev,self.emiss_tot_interp2.derivative()(new_T_J))
        ax1.set_xlabel(r'Electron Temperature ($keV$)')
        ax1.set_ylabel(r'Carbon Radiation Emissivity ($W*m^3$)')
        ax1.loglog(new_T_kev,self.emiss_tot_interp2(new_T_J))
    
        #ax1.semilogy(np.linspace(1.0E-3,1.0E2,10000),10**(emiss_tot_interp(np.linspace(1.0E-3,1.0E2,10000))-100.0),lw=1)

        #ax1.set_xlim(2,10)
        #ax1.set_ylim(0,1E-34)
        #ax1.plot(np.linspace(1.0E-3,1.0E2,10000),10**(emiss_tot_interp(np.linspace(1.0E-3,1.0E2,10000))-100.0),lw=1)

        #emiss_deriv_fig = plt.figure(figsize=(6,6))
        #ax1 = emiss_deriv_fig.add_subplot(1,1,1)
        #ax1.loglog(np.logspace(-3,2,1000),10**(emiss_tot_interp.derivative()(np.logspace(-3,2,1000))-100.0))
