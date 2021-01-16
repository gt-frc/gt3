#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from deprecation import deprecated

@deprecated(deprecated_in="0.0.3", removed_in="0.0.4", details="NBI-related calculations are now performed in the NBI module")
def calc_en_src_nbi(beam, iol_adjusted=False, E_orb_nbi=None):
    """
    """

    # Piper Changes: verify this calculation. Assumes the same things as the particle source eq.

    try:
        # nbi energy source in Joules/sec. Surprisingly, no need to actually calculate anything here.
        en_src_nbi = sum([beam.dPdr.v1D.W[i] for i in range(len(beam.E.J))])

        # nbi energy source in Joules/(m^3 * sec)

        en_src_nbi_tot = sum([beam.dPdV.v1D.W[i] for i in range(len(beam.E.J))])

    except:
        en_src_nbi = beam.dPdr.v1D.W
        en_src_nbi_tot = beam.dPdV.v1D.W

    if iol_adjusted:
        en_src_nbi = en_src_nbi * (1 - E_orb_nbi)
        en_src_nbi_lost = en_src_nbi_tot * (E_orb_nbi)
        en_src_nbi_kept = en_src_nbi_tot * (1 - E_orb_nbi)
    else:
        en_src_nbi_lost = np.zeros(en_src_nbi_tot.shape)
        en_src_nbi_kept = en_src_nbi_tot

    return en_src_nbi, en_src_nbi_tot, en_src_nbi_lost, en_src_nbi_kept
