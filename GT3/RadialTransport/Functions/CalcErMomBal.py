#!/usr/bin/env python2
# -*- coding: utf-8 -*-


def calc_Er_mom_bal(n, e_z, dp_dr, v_tor, v_pol, B_tor, B_pol):
    pres_term = -1.0 / (n * e_z) * dp_dr
    # pres_term_simp = -1.0 * (Lp * T.J) / e_z # alternative method using grad scale length.
    vxb_term = -1.0 * (v_pol * B_tor - v_tor * B_pol)
    E_r = pres_term + vxb_term
    # E_r_simp = pres_term_simp + vxb_term
    return E_r, pres_term, vxb_term