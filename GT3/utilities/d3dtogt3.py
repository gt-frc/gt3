#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 17 10:52:52 2018

@author: max
"""
import numpy as np
import matplotlib.pyplot as plt
def d3dtogt3(infile,input_type,outfile=None,Rfile=None,Zfile=None):
    if input_type == 'rzv': #R,Z,Value (i.e. btor, bpol, psirz, etc.)
        with open(Rfile) as f:
            R1D = np.asarray([float(i) for i in np.asarray(f.read().split())])
        with open(Zfile) as f:
            Z1D = np.asarray([float(i) for i in np.asarray(f.read().split())])
        with open(infile) as f:
            V1D = np.asarray([float(i) for i in np.asarray(f.read().split())])
        R,Z = np.meshgrid(R1D,Z1D)
        
        array = np.column_stack((R.flatten(),Z.flatten(),V1D))
        
        if outfile<>None:
            np.savetxt(outfile,array)
            
        V = np.resize(V1D,(65,65))
        plt.contourf(R,Z,V,500)
    elif input_type == 'rz': #R,Z (i.e. sep, lim)
        #not much to do but clean up and write a new file
        array = np.loadtxt(infile)
        np.savetxt(outfile,array)
    elif input_type == 'profiles':
        with open(infile) as f:
            first_line = f.readline().split(',')
        array = np.loadtxt(infile,skiprows=1, delimiter=',')

        for i,v in enumerate(first_line):
            if v<>'rhor':
                array2 = np.column_stack((array[:,0],array[:,i]))
                np.savetxt('gt3_118888_1525_'+v+'.dat',array2)
    elif input_type == 'erspl':
        array = np.loadtxt(infile,skiprows=2)
        array2 = np.column_stack((array[:,0],array[:,2]))
        np.savetxt(outfile,array2)
        
    return array

if __name__ == "__main__":
    infile = 'p118888_1525_r88_ersplrhob_sc.fit'
    d3dtogt3(infile,'erspl',outfile='gt3_118888_1525_erspl.dat',Rfile='R.txt',Zfile='Z.txt')