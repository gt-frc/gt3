#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from RadialTransport.Functions.CalcReducedMass import calc_reduced_mass
from RadialTransport.Functions.CalcCoulLogJK import calc_coul_log_j_k

def calc_nu_j_k(m_j, m_k, z_j, z_k, T_j, n_k):
    m_r = calc_reduced_mass(m_j, m_k)
    coul_log_j_k = calc_coul_log_j_k(z_j, z_k, T_j, n_k)
    C1 = 1/((((4.8E-10)/(1.6E-12))**1.5)*((4.8E-10)**2.5))
    nu_j_k = 3.34*(coul_log_j_k * (z_j**2) * (z_k**2) * 1E-6 * n_k) / \
            (C1*np.sqrt(m_r*1E3)*(T_j**1.5))
    # nu_j_k = (1.2734E14) * (z_j**2) * (z_k**2) * n_k * coul_log_j_k / \
    #        (np.sqrt(m_r * 1E3) * (T_j**1.5))
    return nu_j_k
