#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import numpy as np

# Calculates the radial electric field including the J x B force from IOL

def calc_Er_iol(n_i, n_e, m_i, n_n, Bphi, p, e_i, T_i, dn_dr, izn_rate, Jr_iol):
    Jr_visc = 0.0
    Jr_neut = -Jr_iol - Jr_visc
    iol_term = (Jr_neut * Bphi**2.0) / (m_i * izn_rate)  # Both the ion and neutral densities are built into izn_rate
    diamag_term = -1.0 * (p.i.L * T_i.J) / e_i  # Pressure gradient is negative, so multiply by -1
    diamag_term_orig = -1.0 * (p.i.derivative() / (e_i * n_i))  # Original form of diamagnetic term
    try:
        dn_dr
        neut_dens_term = -1.0 * (T_i.J * dn_dr) / (e_i * n_n.t)
    except:
        neut_dens_term = [0.] * len(n_i)
    E_r_iol = iol_term + diamag_term + neut_dens_term

    return E_r_iol, iol_term, diamag_term, diamag_term_orig, neut_dens_term
