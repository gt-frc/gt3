#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from GT3.RadialTransport.Functions.CalcT90 import calc_t90
from GT3.RadialTransport.Functions.CalcMbalRHS import calc_mbal_rhs
from scipy import constants

m_d = constants.physical_constants['deuteron mass'][0]
m_c = 12 / constants.N_A / 1E3  # in kg

def calc_vtor_d_pert(vtor_C_total, vtor_C_intrin, vtor_D_intrin, mom_src_ext, z, n, T, B_p, gamma_D, gamma_C):
    """
    """
    # Piper Changes: Changed n to n.i and n.C.
    # Changed z for mbal_rhs_C to 6.
    mbal_rhs_D = calc_mbal_rhs(mom_src_ext, z, n.i, B_p, gamma_D)
    mbal_rhs_C = calc_mbal_rhs(mom_src_ext, z, n.C, B_p, gamma_C) # TODO: This z should be zbar.

    nu_c_DC = 1 / calc_t90(m_d, m_c, 1, 6, n.C, T.i.J)

    #vtor_C_total = vtor_fluid_C + vtor_intrin_C
    del_v0 = vtor_D_intrin - vtor_C_intrin

    nu_drag_approx = (mbal_rhs_D + mbal_rhs_C) / ((n.i * m_d + n.C * m_c) * vtor_C_total + n.i * m_d * del_v0)

    del_v1 = (mbal_rhs_D - n.i.val * m_d * nu_drag_approx * vtor_C_total) / (n.i.val * m_d * (nu_c_DC + nu_drag_approx))

    vtorDPert = vtor_C_total + del_v1
    return vtorDPert
