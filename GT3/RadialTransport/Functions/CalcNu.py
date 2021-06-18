#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from GT3.RadialTransport.Functions.CalcReducedMass import calc_reduced_mass
from GT3.RadialTransport.Functions.CalcCoulLog import calc_coul_log_j_k



def calc_nu_drag(n_j, m_j, v_tor_j, v_tor_k, mbal_rhs, nu_c):
    nu_drag = (mbal_rhs + (n_j * m_j * nu_c * v_tor_k)) / (v_tor_j * n_j * m_j) - nu_c  # Piper Changes: Multiplied V_tor_k in the numerator by n_j*m_j. The algebra was wrong.
    return nu_drag

def calc_nu_j_k(m_j, m_k, z_j, z_k, T_j, n_k):
    m_r = calc_reduced_mass(m_j, m_k)
    coul_log_j_k = calc_coul_log_j_k(z_j, z_k, T_j, n_k)
    C1 = 1/((((4.8E-10)/(1.6E-12))**1.5)*((4.8E-10)**2.5))
    nu_j_k = 3.34*(coul_log_j_k * (z_j**2) * (z_k**2) * 1E-6 * n_k) / \
            (C1*np.sqrt(m_r*1E3)*(T_j**1.5))
    # nu_j_k = (1.2734E14) * (z_j**2) * (z_k**2) * n_k * coul_log_j_k / \
    #        (np.sqrt(m_r * 1E3) * (T_j**1.5))
    return nu_j_k

def calc_nustar(nu90, q, R0_a, vpol):
    nustar = nu90 * abs(q.val) * R0_a / np.abs(vpol)
    return nustar