#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 24 13:25:51 2018

@author: max
"""
from subprocess import call
import numpy as np
from scipy.interpolate import interp2d, interp1d
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
        call([inp.adpak_loc+'adpak', os.getcwd()+'/toadpak'])
        self.read_adpak_outfile(inp,brnd)
        self.coronal_eq(inp,brnd)
        pass
    
    def prep_adpak_infile(self,inp,brnd):
        
        #TODO: Convert adpack routines to python binary using f2py so we can
        #eliminate the main.f driver program and do what we need to in python.
        
        #TODO: Extend this to do elements other than carbon
        self.inucz = 6
        self.nte = 21
        self.nne = 2
        self.nmin = -3
        self.nmax = 2
        self.tei = np.logspace(self.nmin,self.nmax,self.nte) #actual temperatures
        self.tei_lin = np.linspace(self.nmin,self.nmax,self.nte)
        print self.tei
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
        f.write('\n' + '  anneut = 1.0e11')
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
        print
        print frac_abun
        print
        print frac_abun.shape
        print
        print frac_abun[0,:,6]
        #CREATE INTERPOLATION FUNCTION FOR EACH CHARGE STATE
        #self.frac_abun_interp[cs] takes two arguments: 
        #    1) log10 of the T_e in kev
        #    2) n_e in m^-3
        self.brnd_nC_cs     = []
        self.brnd_rad_cs    = []
        self.brnd_rad = np.zeros(brnd.nC.shape)
        for i in np.arange(self.inucz+1):
            #create the interpolation function to get fractional abundances of each charge state as a function
            # of electron density and electron temperature
            #calculate those fractional abundances for all points in the plasma
            #calculate impurity densities by charge state
            frac_abun_interp = interp2d(self.tei_lin,self.anei*1E6,frac_abun[:,:,i]) 
            
            brnd_frac_abun  = np.zeros(brnd.nC.shape)
            for (r,theta),nC in np.ndenumerate(brnd.nC):
                brnd_frac_abun[r,theta] = frac_abun_interp( np.log10(brnd.Te_kev[r,theta]) , brnd.ne[r,theta] )
            #if i==6:
            #    print '###########################3'
            #    print brnd.Te_kev[-1,0]
            #    print np.log10(brnd.Te_kev[-1,0])
            #    print brnd.ne[-1,0]
            #    print frac_abun_interp(np.log10(brnd.Te_kev[-1,0]), brnd.ne[-1,0])
            #    print 
            #    print brnd_frac_abun[-1,0]
            #    print '###########################3'
            #    sys.exit()
            self.brnd_nC_cs.append( brnd_frac_abun * brnd.nC )

            #create the interpolation function to get radiation rate(?) as a function of electron density and temperature            
            rad_interp = interp2d(self.tei_lin,self.anei*1E6,self.alradr[i,:,:]) 
            
            #get those radiation rates for all points in the plasma
            brnd_rad_rate = np.zeros(brnd.nC.shape)
            for (r,theta),nC in np.ndenumerate(brnd.nC):
                brnd_rad_rate[r,theta] = 10**rad_interp( np.log10(brnd.Te_kev[r,theta]) , brnd.ne[r,theta] )
            self.brnd_rad_cs.append( brnd_rad_rate * self.brnd_nC_cs[i] )
            self.brnd_rad = self.brnd_rad + self.brnd_rad_cs[i]
            
        radfig = plt.figure(figsize=(6,6))
        ax1 = radfig.add_subplot(1,1,1)
        ax1.axis('equal')
        ax1.contourf(brnd.R,brnd.Z,self.brnd_rad,500,cmap='jet')

        
        sys.exit()
