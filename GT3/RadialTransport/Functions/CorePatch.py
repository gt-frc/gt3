#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from collections import namedtuple
from GT3 import Core
from GT3.Core.Functions.CalcFSA import calc_fsa
import numpy as np

def corePatch(core: Core, neutFlag=True):
    """
    Updates deuterium ion density and zeff. D density is invalid if no D density file is given because it will be
    set to 0. This screws up subsequent calculations in non-obvious ways. As an example, z_eff calculation will
    include D density as 0s but still have a carbon density from the fracz input file, giving an invalid result

    This is addressed by recreating the n_fsa namedtuple and z_eff_fsa in full, as namedtuples cannot be changed piecemeal

    :param neutFlag:
    :param core:
    :return:
    """


    if neutFlag:
        core.n = namedtuple('n', 'i e n C')(core.n.e/(1.+.025*6.0), core.n.e, core.n.n, 0.025 * core.n.e/(1.+.025*6.0))   # TODO: Update 0.025 and 6.0 to actual fracz and zbar2
    else:
        core.n.update_neutrals()
        core.n = namedtuple('n', 'i e n C')(core.n.e / (1. + .025 * 6.0), core.n.e,
                                   namedtuple('n', 's t tot')(
                                       np.zeros(core.n.i.shape),  # slow
                                       np.zeros(core.n.i.shape),  # thermal
                                       np.zeros(core.n.i.shape)),  # total
                                   0.025 * core.n.e / (1. + .025 * 6.0))

    core.z_eff_fsa = calc_fsa((core.n.i * (1.**2) + core.n.C * (6.0**2))/(core.n.i * 1.0 + core.n.C * 6.0), core.R, core.Z) #TODO Similar updates (1.0 = atnum, 6.0 = zbar2)
