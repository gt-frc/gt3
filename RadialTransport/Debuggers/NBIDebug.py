#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np

import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline



def nbiDebug(rho, a, Qi, Qe, Si, Mi, dVdrho):
    """ """
    dVdr = dVdrho(rho) / a

    fig = plt.figure(figsize=(12, 8))
    fig.tight_layout()
    fig.suptitle(r'Radial Transport Debug Info: NBI')
    ax1 = fig.add_subplot(221)
    ax1.set_title(r'$Q^{nbi}_i$', fontsize=16)
    ax1.set_ylabel(r"""$Q^{nbi}_i$
                    $\left[\frac{W}{m^3}\right]$""", fontsize=16, rotation=0, ha='right')
    ax1.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax1.plot(rho, Qi)

    ax2 = fig.add_subplot(222)
    ax2.set_title(r'$Q^{nbi}_e$', fontsize=16)
    ax2.set_ylabel(r"""$Q^{nbi}_e$
                    $\left[J/{m^3}\right]$""", fontsize=16, rotation=0, ha='right')
    ax2.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax2.plot(rho, Qe)

    ax3 = fig.add_subplot(223)
    ax3.set_title(r'$S^{nbi}_i$', fontsize=16)
    ax3.set_ylabel(r"""$S^{nbi}_i$
                    $\left[\frac{\#}{m^3}\right]$""", fontsize=16, rotation=0, ha='right')
    ax3.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax3.plot(rho, Si)

    ax4 = fig.add_subplot(224)
    ax4.set_title(r'$M^{nbi}_i$', fontsize=16)
    ax4.set_ylabel(r"""$M^{nbi}_i$
                   $\left[\frac{N s}{m^3}\right]$""", fontsize=16, rotation=0, ha='right')
    ax4.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax4.plot(rho, Mi)

    enerTot = UnivariateSpline(rho * a, Qi * dVdr, k=3, s=0).integral(0., a)
    partTot = UnivariateSpline(rho * a, Si * dVdr, k=3, s=0).integral(0., a)
    volTot = UnivariateSpline(rho * a, dVdr, k=3, s=0).integral(0., a)

    print r"""NBI Debug info:

            {:<16}          {}
            {:<16}          {} MW
            {:<16}          {} $m^3$
            """.format("Total # particles", str(partTot),
                       "Total energy", str(enerTot / (1E6)),
                       "Total volume", str(volTot)
                       )

    plt.show(block=False)