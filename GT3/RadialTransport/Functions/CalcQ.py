#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from scipy.interpolate import UnivariateSpline
from math import sqrt, pi, exp
from GT3.RadialTransport.Functions.CalcCoulLog import calc_coul_log
from scipy import constants
import numpy as np
import matplotlib.pyplot as plt

eps_0 = constants.epsilon_0
e = constants.elementary_charge
m_e = constants.electron_mass
m_d = constants.physical_constants['deuteron mass'][0]
m_c = 12 / constants.N_A / 1E3  # in kg


def calc_qie(n, T, ion_species='D'):
    """Calculates collisional energy transfer between an ion species and electrons
    Reference Equation 4.90 in Stacey's Fusion Plasma Physics Book

    :param n:
    :param T:
    :param ion_species:
    :return:
    """

    if ion_species == 'C':
        zi = 1
        Ti = T.i.J  # assumed carbon is at the same temperature as main ion species
        mi = m_c
    else:
        zi = 1
        Ti = T.i.J
        mi = m_d

    coul_log = calc_coul_log(zi, 1, T.i.J, n.i)  # TODO: check these inputs

    qie = n.e * (zi * e**2)**2 * m_e * coul_log * (1 - Ti / T.e.J) / \
          (2*pi * eps_0**2 * np.sqrt(2*pi*m_e*T.e.J) * mi * (1 + 4*sqrt(pi)/3 * (3*m_e*Ti/(2*mi*T.e.J))**1.5))

    return qie * n.i.val  #multiply this shit by density
