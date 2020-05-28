#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt

def IOLDebug(rhor, F_orb_d, E_orb_d, M_orb_d):
    ticks = [0.8, 0.85, 0.9, 0.95, 1.0]
    fig = plt.figure(figsize=(12, 8))
    fig.tight_layout()
#    fig.suptitle(r'Deuterium loss fractions')
    fig.subplots_adjust(wspace=.25)

    ax1 = fig.add_subplot(131)
#    ax1.set_title(r'', fontsize=16)
    ax1.set_ylabel(r"""$\frac{\partial F}{\partial r}$""", fontsize=16, rotation=0, ha='right')
    ax1.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10, x=-.05)
    ax1.set_xticks(ticks)
    ax1.plot(rhor[-100:], F_orb_d[-100:], color="blue")


    ax2 = fig.add_subplot(132)
#    ax2.set_title(r'$Energy loss fraction$', fontsize=16)
    ax2.set_ylabel(r"""$\frac{\partial E}{\partial r}$""", fontsize=16, rotation=0, ha='right')
    ax2.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10, x=-.05)
    ax2.set_xticks(ticks)
    ax2.plot(rhor[-100:], E_orb_d[-100:], color="red")

    ax3 = fig.add_subplot(133)
#    ax3.set_title(r'$Momentum loss fraction$', fontsize=16)
    ax3.set_ylabel(r"""$\frac{\partial M}{\partial r}$""", fontsize=16, rotation=0, ha='right')
    ax3.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10, y=-10, x=-.05)
    ax3.set_xticks(ticks)
    ax3.plot(rhor[-100:], M_orb_d[-100:], color="green")
