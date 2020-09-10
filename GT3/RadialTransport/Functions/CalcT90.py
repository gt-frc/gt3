#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from math import sqrt, pi
from GT3.RadialTransport.Functions.CalcReducedMass import calc_reduced_mass
from GT3.RadialTransport.Functions.CalcCoulLog import calc_coul_log
from scipy.constants import constants

e = constants.elementary_charge
eps_0 = constants.epsilon_0


def calc_t90(m1, m2, z1, z2, n2, T_J):
    """

    :param m1:
    :param m2:
    :param z1:
    :param z2:
    :param n2:
    :param T_J:
    :return:
    """
    m_r = calc_reduced_mass(m1, m2)
    coul_log = calc_coul_log(z1, z2, T_J, n2)
    t_90 = 2*pi*sqrt(m_r) * eps_0**2 * (3*T_J)**1.5 / \
           (n2 * ((z1*e) * (z2 * e))**2 * coul_log) # Piper Changes: Put the denominator in parenthises.
    return t_90