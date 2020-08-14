#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from GT3.Functions.GetNum import getNum
from collections import namedtuple

def getVals(s, t, f):
    """

    This function interactively obtains the 0D parameters for your input file to be created. The commented-out section
    is data for a shot just for debugging purposes

    :param s: shot id
    :param t: time id
    :param f: file to be generated
    :return: Your mom
    """

    data ={}
    """ Function that gets 0D plasma values from user input """

    print """Enter 0D plasma parameters for shot %s.%s

             File will be generated at inputs/%s""" % (str(s), str(t), str(f))

    # data={'eq1' : 1.6E-19,
    #       'eq2' : 9.6E-19,
    #       'xmas1' : 3.35E-27,
    #       'xmas2' : 2.01E-26,
    #       'ephia' : .04,
    #       'xk' : 1.6E-19,
    #       'delma' : 0.005,
    #       'aminor' : 0.598,
    #       'bphi' : -2.04,
    #       'rmajor' : 1.7,
    #       'kappaup': 1.82,
    #       'kappalo' : 1.82,
    #       'triup': .237,
    #       'trilo' : 0.,
    #       'thetapts' : 30,
    #       'rhopts_edge' : 100,
    #       'rhopts_core' : 10,
    #       'xptR' : 1.48,
    #       'xptZ' : -1.24,
    #       'jknot' : 2850000,
    #       'plasmaCur' : 1.38,
    #       'q95' : 3.57,
    #       'ebeam' : 77.48,
    #       'pbeam' : 4.6,
    #       'rtang' : 1.09,
    #       'bknot' : 2.0,
    #       'pwrfrac1' : .76,
    #       'pwrfrac2' : .13,
    #       'pwrfrac3' : .11,
    #       'epsknot' : 1.265,
    #       'epssep' : 1.82,
    #       'shftknot' : 0.033
    #       }

    data['eq1'] = getNum("Enter charge of main ion: ", "f")
    data['eq2'] = getNum("Enter charge of main impurity species: ", "f")
    data['xmas1'] = getNum("Enter main ion mass: ", "f")
    data['xmas2'] = getNum("Enter main impurity mass: ", "f")
    data['ephia'] = getNum("Enter the toroidal electric field: ", "f")
    data['xk'] = 1.6E-19
    data['delma'] = .005
    data['aminor'] = getNum("Enter the plasma radius (a minor): " ,"f")
    data['bphi'] = getNum("Enter the toroidal plasma field strength in T (abs mag): ", "f")
    data['rmajor'] = getNum("Enter R major: ", "f")
    data['kappaup'] = getNum("Enter upper elongation (kappa): ", "f")
    data['kappalo'] = getNum("Enter lower elongation (kappa): ", "f")
    data['triup'] = getNum("Enter upper triangularity (delta): ", "f")
    data['trilo'] = getNum("Enter lower triangularity (delta): ", "f")
    data['thetapts'] = getNum("Enter approximate number of theta points (typically 30): ", 'i')
    data['rhopts_edge'] = getNum("Enter rho points in the edge (typically 100): ", 'i')
    data['rhopts_core'] = getNum("Enter rho points in the core (typically 10): ", 'i')
    data['xptR'] = getNum("Enter the X-point R coordinate: ", 'f')
    data['xptZ'] = getNum("Enter the X-point Z coordinate: ", 'f')
    data['jknot'] = getNum("Enter the r=0 plasma current density (in A/m^3): ", 'f')
    data['plasmaCur'] = getNum("Enter the plasma current (in MA): ", 'f')
    data['q95'] = getNum("Enter q95: ", 'f')
    data['ebeam'] = getNum("Enter beam ion energy in eV: ", 'f')
    data['pbeam'] = getNum("Enter beam power in MWYup: ", 'f')
    data['rtang'] = getNum("Enter radius of tangency in cm: ", 'f')
    data['bknot'] = abs(data['bphi'])
    data['pwrfrac1'] = getNum("Enter fraction of beam power to D1: ", 'f')
    data['pwrfrac2'] = getNum("Enter fraction of beam power to D2: ", 'f')
    data['pwrfrac3'] = getNum("Enter fraction of beam power to D3: ", 'f')
    data['epsknot'] = getNum("Enter epsilon at the plasma center: ", 'f')
    data['epssep'] = getNum("Enter epsilon at separatrix: ", 'f')
    data['shftknot'] = getNum("Enter the Shavranof shift at the plasma center: ", 'f')

    return namedtuple('data', sorted(data))(**data)