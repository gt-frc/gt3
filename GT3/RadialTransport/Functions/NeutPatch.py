#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from GT3.Core.Functions.CalcFSA import calc_fsa
from collections import namedtuple
import numpy as np


def neutPatch(core):
    """

    Re-runs FSA of quantities needed in rad transport that are not available because neutrals update is run after
    these values, which depend on neutral data, are calculated in core
    :param core:
    :return:
    """
    core.izn_rate_fsa_s = calc_fsa(core.izn_rate.s, core.R, core.Z)
    core.izn_rate_fsa_t = calc_fsa(core.izn_rate.t, core.R, core.Z)
    core.izn_rate_fsa = calc_fsa(core.izn_rate.s + core.izn_rate.t, core.R, core.Z)

    # Piper chages: Use carbon Lz for the cool_rate calculation
    core.cool_rate_fsa = calc_fsa(core.n.e * core.n.C * np.nan_to_num(core.Lz_C.s) + \
    +                            core.n.e * core.n.C * np.nan_to_num(core.Lz_C.t),
                                  core.R, core.Z)

    """
    Calculate the neutrals now
    """

    core.dn_dr_fsa = calc_fsa(core.dn_dr, core.R, core.Z)

    """
    Cuts off neutrals calculations at 0.8 to give realistic neutral shit
    """

    core.izn_rate_fsa = np.array(map(lambda x: core.izn_rate_fsa[-1] * x**10, core.r[:,0]/core.a))
    core.cool_rate_fsa = np.array(map(lambda x: core.cool_rate_fsa[-1] * x**10, core.r[:,0]/core.a))
    core.dn_dr_fsa = np.array(map(lambda x: core.dn_dr_fsa[-1] * x**10, core.r[:,0]/core.a))


    core.n_fsa = namedtuple('n', 'i e n C')(
        calc_fsa(core.n.i, core.R, core.Z),
        calc_fsa(core.n.e, core.R, core.Z),
        namedtuple('n', 's t tot')(
            np.array(map(lambda x: core.n_fsa.n.s[-1] * x ** 10, core.r[:, 0] / core.a))/1E3,  # slow
            np.array(map(lambda x: core.n_fsa.n.t[-1] * x ** 10, core.r[:, 0] / core.a))/1E3,  # thermal
            np.array(map(lambda x: core.n_fsa.n.tot[-1] * x ** 10, core.r[:, 0] / core.a))/1E3),  # total
        calc_fsa(core.n.C, core.R, core.Z))