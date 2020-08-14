#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

def QDebug(rho, a, r2sa, dVdrho, Qie, Qnbi_d_i, Qnbi_kept_i, coolrate, Eorb, Qi, T):
    # TODO: FINISH THIS
    dVdr = dVdrho(rho)/a

    fig = plt.figure(figsize=(12, 8))
    fig.tight_layout()
    fig.suptitle(r'Radial Transport Debug Info: Deuterium Radial Heat Flux')

    ax1 = fig.add_subplot(231)
    ax1.set_title(r'$Q^{nbi}_{diff,i}$', fontsize=16)
    ax1.set_ylabel(r"""$Q_i$
                       $\left[\frac{J}{m^3 s}\right]$""", fontsize=16, rotation=0, ha='right')
    ax1.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax1.plot(rho, Qnbi_d_i)

    ax2 = fig.add_subplot(232)
    ax2.set_title(r'$Q^{nbi}_{int,i}$', fontsize=16)
    ax2.set_ylabel(r"""$Q_i$
                       $\left[\frac{J}{m^3 s}\right]$""", fontsize=16, rotation=0, ha='right')
    ax2.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax2.plot(rho, Qnbi_kept_i)

    ax3 = fig.add_subplot(233)
    ax3.set_title(r'$Q_{ioniz}$', fontsize=16)
    ax3.set_ylabel(r"""$Q_i$
                       $\left[\frac{J}{m^3 s}\right]$""", fontsize=16, rotation=0, ha='right')
    ax3.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax3.plot(rho, coolrate)

    ax4 = fig.add_subplot(234)
    ax4.set_title(r'$Q_{ie}$', fontsize=16)
    ax4.set_ylabel(r"""$Q_{r,ie}$
                        $\left[\frac{J}{m^2 s}\right]$""", fontsize=16, rotation=0, ha='right')
    ax4.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax4.plot(rho, Qie)

    ax5 = fig.add_subplot(235)
    ax5.set_title(r'$dEdr$', fontsize=16)
    ax5.set_ylabel(r"""$dEdr{r}$
                        $\left[\frac{1}{m}\right]$""", fontsize=16, rotation=0, ha='right')
    ax5.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax5.plot(rho, Eorb)

    ax6 = fig.add_subplot(236)
    ax6.set_title(r'Temperatures (red = ion)', fontsize=16)
    ax6.set_ylabel(r"""T{r}$
                        $\left[eV\right]$""", fontsize=16, rotation=0, ha='right')
    ax6.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax6.plot(rho, T.i.ev, color='red')
    ax6.plot(rho, T.e.ev, color='blue')

    plt.show(block=False)

    energyIn = UnivariateSpline(rho*a, (Qnbi_d_i - coolrate) * dVdr, k=3, s=0).integral(0., a)
    energyOut = Qi[-1]*r2sa(a)
    totVol = UnivariateSpline(rho*a, dVdr, k=3, s=0).integral(0., a)

    print r"""Radial Energy Flux Debug info:

            {:<16}          {}  MW
            {:<16}          {}  MW
            {:<16}          {}
            {:<16}          {}
            """.format("Total energy in", str(energyIn / (1E6)),
                       "Total energy out", str(energyOut / (1E6)),
                       "Total volume", str(totVol),
                       "Surface area at LCFS", str(r2sa(a))
                       )
    # print r"""r          Energy in      Energy out"""
    # for x in range(len(rho)):
    #     strList=balance(Qi, (Qnbi_kept_i - coolrate) * dVdr, rho*a, r2sa, x)
    #     strDiff = 100.*np.abs(strList[1]-strList[2])/np.average(np.array(strList))
    #     print """{:<15}    {:<15}     {:<15}      {:<15}%""".format(strList[0], strList[1], strList[2], str(strDiff))