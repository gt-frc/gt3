#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np


def calc_nu_drag(n_j, m_j, v_tor_j, v_tor_k, mbal_rhs, nu_c):
    nu_drag = (mbal_rhs + (n_j * m_j * nu_c * v_tor_k)) / (v_tor_j * n_j * m_j) - nu_c  # Piper Changes: Multiplied V_tor_k in the numerator by n_j*m_j. The algebra was wrong.
    return nu_drag