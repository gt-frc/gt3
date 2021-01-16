#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from deprecation import deprecated


@deprecated(deprecated_in="0.0.3", removed_in="0.0.4", details="NBI-related calculations are now performed in the NBI module")
def calc_part_src_nbi(beam, iol_adjusted=False, F_orb_nbi=None):
    """
    """
    # # old suspicious calculation
    # snbi = (0.624E25 * beam.dPdV.v1D.W / (vol * beam.E.J)) * beam.num * beam.dp

    # Piper Changes: Changed function to calculate particle source and source density.
    # tot, lost, and kept are all densities. The first return variable is a total source.
    # particle source needs to be split into a total, amount kept, and amount lost to calculate the return current.

    # These equations assume:
    # beam.dPdV.v1D.W = H(rho) * Pbeam * pwrfrac / Volp_nbi
    # beam.dPdr.v1D.W = beam.dPdV.v1D.W * dV

    # The factor of 2 can be applied to F_orb in the continuity equation later.

    # nbi source in particles/sec. Needed to calculate gamma with sources.
    try:

        part_src_nbi = sum([beam.dPdr.v1D.W[i] / beam.E.J[i] for i in range(len(beam.E.J))])

        # nbi source in # of particles/(m^3 * sec). Needed to calculate gamma with source densities.

        part_src_nbi_tot = sum([beam.dPdV.v1D.W[i] / beam.E.J[i] for i in range(len(beam.E.J))])

    except:
        part_src_nbi = beam.dPdr.v1D.W / beam.E.J
        part_src_nbi_tot = beam.dPdV.v1D.W / beam.E.J


    if iol_adjusted:
        part_src_nbi = part_src_nbi * (1 - F_orb_nbi)
        part_src_nbi_lost = part_src_nbi_tot * (F_orb_nbi)
        part_src_nbi_kept = part_src_nbi_tot * (1 - F_orb_nbi)
    else:
        part_src_nbi_lost = np.zeros(part_src_nbi_tot.shape)
        part_src_nbi_kept = part_src_nbi_tot

    return part_src_nbi, part_src_nbi_tot, part_src_nbi_lost, part_src_nbi_kept
