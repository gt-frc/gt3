#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

def gammaDebug(rho, a, r2sa, dVdrho, Snbi_d_i, Snbi_kept_i, Sizn, gamma_i_D, gamma_D):

    dVdr = dVdrho(rho)/a

    fig = plt.figure(figsize=(12, 8))
    fig.tight_layout()
    fig.suptitle(r'Radial Transport Debug Info: Deuterium Radial Particle Flux')

    ax1 = fig.add_subplot(231)
    ax1.set_title(r'$S^{nbi}_{diff,i}$', fontsize=16)
    ax1.set_ylabel(r"""$S_i$
                       $\left[\frac{\#}{m^3 s}\right]$""", fontsize=16, rotation=0, ha='right')
    ax1.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax1.plot(rho, Snbi_d_i)

    ax2 = fig.add_subplot(232)
    ax2.set_title(r'$S^{nbi}_{int,i}$', fontsize=16)
    ax2.set_ylabel(r"""$S_i$
                       $\left[\frac{\#}{m^3 s}\right]$""", fontsize=16, rotation=0, ha='right')
    ax2.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax2.plot(rho, Snbi_kept_i)

    ax3 = fig.add_subplot(233)
    ax3.set_title(r'$S_{ion}$', fontsize=16)
    ax3.set_ylabel(r"""$S_i$
                       $\left[\frac{\#}{m^3 s}\right]$""", fontsize=16, rotation=0, ha='right')
    ax3.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax3.plot(rho, Sizn)

    ax4 = fig.add_subplot(234)
    ax4.set_title(r'$\Gamma_{r,int}$', fontsize=16)
    ax4.set_ylabel(r"""$\Gamma_{r}$
                        $\left[\frac{\#}{m^2 s}\right]$""", fontsize=16, rotation=0, ha='right')
    ax4.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax4.plot(rho, gamma_i_D)

    ax5 = fig.add_subplot(235)
    ax5.set_title(r'$\Gamma_{r,diff}$', fontsize=16)
    ax5.set_ylabel(r"""$\Gamma_{r}$
                        $\left[\frac{\#}{m^2 s}\right]$""", fontsize=16, rotation=0, ha='right')
    ax5.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax5.plot(rho, gamma_D)

    plt.show(block=False)

    partIn = UnivariateSpline(rho*a, (Snbi_kept_i + Sizn) * dVdr, k=3, s=0).integral(0., a)
    partOut = gamma_i_D[-1]*r2sa(a)
    totVol = UnivariateSpline(rho*a, dVdr, k=3, s=0).integral(0., a)

    print r"""Radial Particle Flux Debug info:

            {:<16}          {}
            {:<16}          {}
            {:<16}          {}
            {:<16}          {}
            """.format("Total # particles in", str(partIn),
                       "Total # particles out", str(partOut),
                       "Total volume", str(totVol),
                       "Surface area at LCFS", str(r2sa(a))
                       )
    print r"""r          Particles in      Particles out"""
    # for x in range(len(rho)):
    #     strList=balance(gamma_i_D, (Snbi_kept_i + Sizn) * dVdr, rho*a, r2sa, x)
    #     strDiff = 100.*np.abs(strList[1]-strList[2])/np.average(np.array(strList))
    #     print """{:<15}    {:<15}     {:<15}      {:<15}%""".format(strList[0], strList[1], strList[2], str(strDiff))