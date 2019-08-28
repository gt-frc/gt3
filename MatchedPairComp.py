#!/usr/bin/python
"""
Created on Thu Dec 28 02:13:04 2017

@author: Jonathan
"""
###############################################################################
#
#    Matched shot pair analysis and comparison module based on GT3
#
###############################################################################

import GTEDGE3_cli
import graphs as graphs
import matplotlib.pyplot as plt
import pickle
import numpy as np
from math import ceil, floor, log10

shotid, timeids = [123301, 123302], [2800, 2810]
runid=''

###############################################################################
#
#   Utilities
#
##############################################################################

def yRangeFind(inList):
    x=np.concatenate(inList)
    flat=x.flatten()
    flat.sort()
    #return [flat[2]*.75,flat[-2]*1.25]
    a=long(flat[2]*.75)
    b=long(flat[-2]*1.25)
    try:
        aval=1*10**int(floor(log10(a)))
        if (ceil(log10(b))-(log10(b))) <=.25:
            return [0 if (aval == 1) else aval, .75*10**int(ceil(log10(b)))]
        elif (ceil(log10(b))-(log10(b))) <=.5:
            return [0 if (aval == 1) else aval, .5*10**int(ceil(log10(b)))]
        elif (ceil(log10(b))-(log10(b))) <=.75:
            return [0 if (aval == 1) else aval, .25*10**int(ceil(log10(b)))]
        else:
            return [0 if (aval == 1) else aval, .1*10**int(ceil(log10(b)))]
    # negative arg in log exception
    except:
        return [flat[2]*.75,flat[-2]*1.25]



###############################################################################
#
#   Chi comparison modules
#
##############################################################################

def chiComp(shotA, shotB):
    #prettyID = [r'$Q^{diff}_{j,r}=Q^{total}_j-Q^{conv}_j-Q^{heatin}_j-Q^{visc}_j$ w/ IOL',
    #            r'$Q^{diff}_{j,r}=Q^{total}_j-Q^{conv}_j-Q^{heatin}_j-Q^{visc}_j$ w/out IOL']
    prettyID = [r'$\chi_{j,r}^{%s}$' % str(shotid[0]),
                r'$\chi_{j,r}^{%s}$' % str(shotid[1])]

    title = r"Ion Heat Conductivity Comparison"

    adjustments = {0: 1.}
    caption = ["Main ion heat conductivity comparison in the edge for DIII-D shot %s.%s, %s.%s w/IOL correction" % (str(shotid[0]), str(timeids[0]), str(shotid[1]), str(timeids[1])),
               "convective heat flow,  work done on plasma and visc. heat"]

    graphs.prettyCompare(('rhor', shotA.rtrans.rhor[:-1]), [(prettyID[0], shotA.rtrans.chi.chi.i.chi4[:-1]), (prettyID[1], shotB.rtrans.chi.chi.i.chi4[:-1])],
                         yrange=yRangeFind([shotA.rtrans.chi.chi.i.chi4[-60:], shotB.rtrans.chi.chi.i.chi4[-60:]]),
                         #yrange=[-5., 5],
                         datalabels=[prettyID[0], prettyID[1]],
                         title=title,
                         ylabel=[r"$\chi_j$", r'$\left[\frac{m^2}{s}\right]$'],
                         caption=caption,
                         adjust=adjustments,
                         xTickremove=None,
                         textSpace=.025,
                         marginBottom=.20,
                         marginLeft=0.095,
                         capAdj=-0.15,
                         yLabelAdj=-.25,
                         xLabelAdj=-.05,
                         size=(16, 12),
                         legend=True)

def fluxComp(shotA, shotB):
    prettyID = [r'$\Gamma_{r,j}^{%s}$' % str(shotid[0]),
                r'$\Gamma_{r,j}^{%s}$' % str(shotid[1])]

    caption = ["Ion radial particle flux for DIII-D shot %s.%s and %s.%s" % (str(shotid[0]), str(timeids[0]), str(shotid[1]), str(timeids[1]))]
    title = r"$\Gamma_{r,j}$ comparison"
    adjustments = {}
    graphs.prettyCompare(('rhor', shotA.rtrans.rhor), [(prettyID[0], shotA.rtrans.gamma_diff_D), (prettyID[1], shotB.rtrans.gamma_diff_D)],
                         #yrange=(yRangeFind([shotA.rtrans.gamma_diff_D[-60:], shotB.rtrans.gamma_diff_D[-60:]])),
                         yrange=[1E17, 1E20],
                         datalabels=[prettyID[0], prettyID[1]],
                         title=title,
                         ylabel=[r"$\Gamma_{r,j}$", r"$\left[\frac{\#}{m^2 s}\right]$"],
                         caption=caption,
                         adjust=adjustments,
                         textSpace=.02,
                         yLabelAdj=-.8,
                         xLabelAdj=-0.025,
                         marginLeft=0.15,
                         marginBottom=0.2,
                         size=(16, 19),
                         legend=True)

###############################################################################
#
#   Heating comparison modules
#
##############################################################################

def heatsComp(shotA, shotB):

    fig = plt.figure(figsize=(12, 8))
    fig.tight_layout()
    fig.suptitle(r'Heating sources')
    ax1 = fig.add_subplot(221)
    ax1.set_title(r'$Q^{conv}_i$', fontsize=16)
    ax1.set_ylabel(r"""$Q^{conv}_i$
                    $\left[\frac{W}{m^3}\right]$""", fontsize=16, rotation=0, ha='right')
    ax1.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax1.plot(shotA.rtrans.rhor[:-1][-60:], shotA.rtrans.chi.conv25[:-1][-60:], label="%s" % str(shotid[0]))
    ax1.plot(shotB.rtrans.rhor[:-1][-60:], shotB.rtrans.chi.conv25[:-1][-60:], label="%s" % str(shotid[1]))
    ax1.legend()

    ax2 = fig.add_subplot(222)
    ax2.set_title(r'$Q^{visc}_i$', fontsize=16)
    ax2.set_ylabel(r"""$Q^{visc}_i$
                    $\left[\frac{W}{m^3}\right]$""", fontsize=16, rotation=0, ha='right')
    ax2.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax2.plot(shotA.rtrans.rhor[:-1][-60:], shotA.rtrans.chi.heatvisc[:-1][-60:], label="%s" % str(shotid[0]))
    ax2.plot(shotB.rtrans.rhor[:-1][-60:], shotB.rtrans.chi.heatvisc[:-1][-60:], label="%s" % str(shotid[1]))
    ax2.legend()

    ax3 = fig.add_subplot(223)
    ax3.set_title(r'$Q^{heatin}_i$', fontsize=16)
    ax3.set_ylabel(r"""$Q^{heatin}_i$
                    $\left[\frac{W}{m^3}\right]$""", fontsize=16, rotation=0, ha='right')
    ax3.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax3.plot(shotA.rtrans.rhor[:-1][-60:], shotA.rtrans.chi.heatin[:-1][-60:], label="%s" % str(shotid[0]))
    ax3.plot(shotB.rtrans.rhor[:-1][-60:], shotB.rtrans.chi.heatin[:-1][-60:], label="%s" % str(shotid[1]))
    ax3.legend()

    ax4 = fig.add_subplot(224)
    ax4.set_title(r'$Q^{total}_i$', fontsize=16)
    ax4.set_ylabel(r"""$Q^{total}_i$
                    $\left[\frac{W}{m^3}\right]$""", fontsize=16, rotation=0, ha='right')
    ax4.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax4.plot(shotA.rtrans.rhor[:-1][-60:], shotA.rtrans.chi.Qi[:-1][-60:], label="%s" % str(shotid[0]))
    ax4.plot(shotB.rtrans.rhor[:-1][-60:], shotB.rtrans.chi.Qi[:-1][-60:], label="%s" % str(shotid[1]))
    ax4.legend()


def convComp(shotA, shotB):

    prettyID = [r'$Q^{conv}_{j,r}^{%s}$' % str(shotid[0]),
                r'$Q^{conv}_{j,r}^{%s}$' % str(shotid[1])]

    title = r"L-H Convective Heating Comparison"

    adjustments = {0: 1.}
    caption = ["Comparison of Convective heating for DIII-D shot %s.%s, %s.%s w/IOL correction" % (str(shotid[0]), str(timeids[0]), str(shotid[1]), str(timeids[1]))]
    graphs.prettyCompare(('rhor', shotA.rtrans.rhor[:-1]), [(prettyID[0], shotA.rtrans.chi.conv25[:-1]), (prettyID[1], shotB.rtrans.chi.conv25[:-1])],
                         yrange=yRangeFind([shotA.rtrans.chi.conv25[-60:], shotB.rtrans.chi.conv25[-60:]]),
                         #yrange=[-5., 5],
                         datalabels=[prettyID[0], prettyID[1]],
                         title=title,
                         ylabel=[r"$Q^{conv}_j$", r'$\left[\frac{W}{m^3}\right]$'],
                         caption=caption,
                         adjust=adjustments,
                         xTickremove=None,
                         textSpace=.025,
                         marginBottom=.20,
                         marginLeft=0.095,
                         capAdj=-0.15,
                         yLabelAdj=-.25,
                         xLabelAdj=-.05,
                         size=(16, 12),
                         legend=True)


##############################################################################
#
#   Plasma profile comparison modules
#
###############################################################################

def plasmaComp(shotA, shotB):

    fig = plt.figure(figsize=(12, 8))
    fig.tight_layout()
    fig.suptitle(r'Plasma profiles')
    ax1 = fig.add_subplot(231)
    ax1.set_title(r'$T_i$', fontsize=16)
    ax1.set_ylabel(r"""$T_i$
                    $[eV]$""", fontsize=16, rotation=0, ha='right')
    ax1.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax1.plot(shotA.rtrans.rhor[:-1][-60:], shotA.rtrans.profiles.T.i.ev[-60:], label="%s" % str(shotid[0]))
    ax1.plot(shotB.rtrans.rhor[:-1][-60:], shotB.rtrans.profiles.T.i.ev[-60:], label="%s" % str(shotid[1]))
    ax1.legend()

    ax2 = fig.add_subplot(232)
    ax2.set_title(r'$T_e$', fontsize=16)
    ax2.set_ylabel(r"""$T_e$
                    $[eV]$""", fontsize=16, rotation=0, ha='right')
    ax2.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax2.plot(shotA.rtrans.rhor[:-1][-60:], shotA.rtrans.profiles.T.e.ev[-60:], label="%s" % str(shotid[0]))
    ax2.plot(shotB.rtrans.rhor[:-1][-60:], shotB.rtrans.profiles.T.e.ev[-60:], label="%s" % str(shotid[1]))
    ax2.legend()

    ax3 = fig.add_subplot(233)
    ax3.set_title(r'$n_i$', fontsize=16)
    ax3.set_ylabel(r"""$n_i$
                    $\left[\frac{1}{m^3}\right]$""", fontsize=16, rotation=0, ha='right')
    ax3.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax3.set_ylim(yRangeFind([shotA.rtrans.profiles.n.i[-60:], shotB.rtrans.profiles.n.i[-60:]]))
    ax3.plot(shotA.rtrans.rhor[:-1][-60:], shotA.rtrans.profiles.n.i[-60:], label="%s" % str(shotid[0]))
    ax3.plot(shotB.rtrans.rhor[:-1][-60:], shotB.rtrans.profiles.n.i[-60:], label="%s" % str(shotid[1]))
    ax3.legend()

    ax4 = fig.add_subplot(234)
    ax4.set_title(r'$n_e$', fontsize=16)
    ax4.set_ylabel(r"""$n_e$
                    $\left[\frac{1}{m^3}\right]$""", fontsize=16, rotation=0, ha='right')
    ax4.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax4.plot(shotA.rtrans.rhor[:-1][-60:], shotA.rtrans.profiles.n.e[-60:], label="%s" % str(shotid[0]))
    ax4.plot(shotB.rtrans.rhor[:-1][-60:], shotB.rtrans.profiles.n.e[-60:], label="%s" % str(shotid[1]))
    ax4.legend()

    ax5 = fig.add_subplot(235)
    ax5.set_title(r'$n_n^{total}$', fontsize=16)
    ax5.set_ylabel(r"""$n_n$
                    $\left[\frac{1}{m^3}\right]$""", fontsize=16, rotation=0, ha='right')
    ax5.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax5.set_ylim(yRangeFind([shotA.rtrans.profiles.nn.tot[-60:], shotB.rtrans.profiles.nn.tot[-60:]]))
    ax5.plot(shotA.rtrans.rhor[:-1][-60:], shotA.rtrans.profiles.nn.tot[-60:], label="%s" % str(shotid[0]))
    ax5.plot(shotB.rtrans.rhor[:-1][-60:], shotB.rtrans.profiles.nn.tot[-60:], label="%s" % str(shotid[1]))
    ax5.legend()

    ax6 = fig.add_subplot(236)
    ax6.set_title(r'$E_r$', fontsize=16)
    ax6.set_ylabel(r"""$E_r$
                    $\left[\frac{V}{m}\right]$""", fontsize=16, rotation=0, ha='right')
    ax6.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax6.set_ylim(yRangeFind([shotA.core.E_r_fsa[-60:], shotB.core.E_r_fsa[-60:]]))
    ax6.plot(shotA.rtrans.rhor[:-1][-60:], shotA.core.E_r_fsa[-60:], label="%s" % str(shotid[0]))
    ax6.plot(shotB.rtrans.rhor[:-1][-60:], shotB.core.E_r_fsa[-60:], label="%s" % str(shotid[1]))
    ax6.legend()

def rotComp(shotA, shotB):

    fig = plt.figure(figsize=(12, 8))
    fig.tight_layout()
    fig.suptitle(r'Plasma rotation profiles')

    ax1 = fig.add_subplot(221)
    ax1.set_title(r'$v^{d}_{\phi}$', fontsize=16)
    ax1.set_ylabel(r"""$v_{\phi}$
                    $[\frac{m}{s}]$""", fontsize=16, rotation=0, ha='right')
    ax1.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax1.plot(shotA.rtrans.rhor[:-1][-60:], shotA.rtrans.vtor_D_total[-60:], label="%s" % str(shotid[0]))
    ax1.plot(shotB.rtrans.rhor[:-1][-60:], shotB.rtrans.vtor_D_total[-60:], label="%s" % str(shotid[1]))
    ax1.legend()

    ax2 = fig.add_subplot(222)
    ax2.set_title(r'$v^{d}_{\theta}$', fontsize=16)
    ax2.set_ylabel(r"""$v_{\theta}$
                    $[\frac{m}{s}]$""", fontsize=16, rotation=0, ha='right')
    ax2.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax2.plot(shotA.rtrans.rhor[:-1][-60:], shotA.rtrans.vpol_D[-60:], label="%s" % str(shotid[0]))
    ax2.plot(shotB.rtrans.rhor[:-1][-60:], shotB.rtrans.vpol_D[-60:], label="%s" % str(shotid[1]))
    ax2.legend()

    ax3 = fig.add_subplot(223)
    ax3.set_title(r'$v^{C}_{\phi}$', fontsize=16)
    ax3.set_ylabel(r"""$v_{\phi}$
                    $[\frac{m}{s}]$""", fontsize=16, rotation=0, ha='right')
    ax3.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax3.plot(shotA.rtrans.rhor[:-1][-60:], shotA.rtrans.vtor_C_total[-60:], label="%s" % str(shotid[0]))
    ax3.plot(shotB.rtrans.rhor[:-1][-60:], shotB.rtrans.vtor_C_total[-60:], label="%s" % str(shotid[1]))
    ax3.legend()

    ax4 = fig.add_subplot(224)
    ax4.set_title(r'$v^{C}_{\theta}$', fontsize=16)
    ax4.set_ylabel(r"""$v_{\theta}$
                    $[\frac{m}{s}]$""", fontsize=16, rotation=0, ha='right')
    ax4.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax4.set_ylim(-20000,20000)
    ax4.plot(shotA.rtrans.rhor[:-1][-60:], shotA.rtrans.vpol_C[-60:], label="%s" % str(shotid[0]))
    ax4.plot(shotB.rtrans.rhor[:-1][-60:], shotB.rtrans.vpol_C[-60:], label="%s" % str(shotid[1]))
    ax4.legend()



def IOLcomp(shotA, shotB):

    fig = plt.figure(figsize=(12, 8))
    fig.tight_layout()
    fig.suptitle(r'Ion orbit loss Profiles')

    ax1 = fig.add_subplot(231)
    ax1.set_title(r'$F_{iol}$', fontsize=16)
    ax1.set_ylabel(r"""$F_{iol}$
                    """, fontsize=16, rotation=0, ha='right')
    ax1.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax1.set_ylim(0,.5)
    ax1.plot(shotA.rtrans.rhor[:-1][-60:], shotA.iol.forb_d_therm_1D[-60:], label="%s" % str(shotid[0]))
    ax1.plot(shotB.rtrans.rhor[:-1][-60:], shotB.iol.forb_d_therm_1D[-60:], label="%s" % str(shotid[1]))
    ax1.legend()

    ax2 = fig.add_subplot(232)
    ax2.set_title(r'$E_{iol}$', fontsize=16)
    ax2.set_ylabel(r"""$E_{iol}$
                    """, fontsize=16, rotation=0, ha='right')
    ax2.set_ylim(0,.5)
    ax2.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax2.plot(shotA.rtrans.rhor[:-1][-60:], shotA.iol.eorb_d_therm_1D[-60:], label="%s" % str(shotid[0]))
    ax2.plot(shotB.rtrans.rhor[:-1][-60:], shotB.iol.eorb_d_therm_1D[-60:], label="%s" % str(shotid[1]))
    ax2.legend()

    ax3 = fig.add_subplot(233)
    ax3.set_title(r'$E_{iol}$', fontsize=16)
    ax3.set_ylabel(r"""$M_{iol}$
                    """, fontsize=16, rotation=0, ha='right')
    ax3.set_ylim(-1,1.)
    ax3.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax3.plot(shotA.rtrans.rhor[:-1][-60:], shotA.iol.morb_d_therm_1D[-60:], label="%s" % str(shotid[0]))
    ax3.plot(shotB.rtrans.rhor[:-1][-60:], shotB.iol.morb_d_therm_1D[-60:], label="%s" % str(shotid[1]))
    ax3.legend()

    ax4 = fig.add_subplot(234)
    ax4.set_title(r'$F^C_{iol}$', fontsize=16)
    ax4.set_ylabel(r"""$F_{iol}$
                    """, fontsize=16, rotation=0, ha='right')
    ax4.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax4.set_ylim(0, .5)
    ax4.plot(shotA.rtrans.rhor[:-1][-60:], shotA.iol.forb_c_therm_1D[-60:], label="%s" % str(shotid[0]))
    ax4.plot(shotB.rtrans.rhor[:-1][-60:], shotB.iol.forb_c_therm_1D[-60:], label="%s" % str(shotid[1]))
    ax4.legend()

    ax5 = fig.add_subplot(235)
    ax5.set_title(r'$E^C_{iol}$', fontsize=16)
    ax5.set_ylabel(r"""$E_{iol}$
                    """, fontsize=16, rotation=0, ha='right')
    ax5.set_ylim(0, .5)
    ax5.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax5.plot(shotA.rtrans.rhor[:-1][-60:], shotA.iol.eorb_c_therm_1D[-60:], label="%s" % str(shotid[0]))
    ax5.plot(shotB.rtrans.rhor[:-1][-60:], shotB.iol.eorb_c_therm_1D[-60:], label="%s" % str(shotid[1]))
    ax5.legend()

    ax6 = fig.add_subplot(236)
    ax6.set_title(r'$M^C_{iol}$', fontsize=16)
    ax6.set_ylabel(r"""$E_{iol}$
                    """, fontsize=16, rotation=0, ha='right')
    ax6.set_ylim(-1, 1.)
    ax6.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax6.plot(shotA.rtrans.rhor[:-1][-60:], shotA.iol.morb_c_therm_1D[-60:], label="%s" % str(shotid[0]))
    ax6.plot(shotB.rtrans.rhor[:-1][-60:], shotB.iol.morb_c_therm_1D[-60:], label="%s" % str(shotid[1]))
    ax6.legend()

if __name__ == "__main__":

    #   Multi-shot compa

    # L mode

    shotargs = {'shotid': shotid[0],
                'timeid': timeids[0],
                'runid': runid,
                'nbRun': True,
                'IOL': True,
                'quiet': True,
                'reNeut': False,
                'neutrals': True,
                'gt3Method': 'radialtrans',
                'debug' : False}
    try:
        shotA = pickle.load(open("outputs/s%s.%s.dat" % (str(shotid[0]), str(timeids[0])), "rb"))
    except:
        try:
            shotA = GTEDGE3_cli.runGT3(shotargs)
        except Exception as e:
            raise Exception(e)
        try:
            pass
#            with open("outputs/s%s.%s.dat"  % (str(shotid[0]), str(timeid[0])), "wb") as f:
#                pickle.dump(shotA, f)
        except Exception as e:
             print e
        try:
            f.close()
        except:
            pass

    # H-mode

    shotargs = {'shotid': shotid[1],
                'timeid': timeids[1],
                'runid': runid,
                'nbRun': True,
                'IOL': True,
                'quiet': True,
                'reNeut': False,
                'neutrals': True,
                'gt3Method': 'radialtrans',
                'debug' : False}
    try:
        shotB = pickle.load(open("outputs/s%s.%s.dat" % (str(shotid[1]), str(timeids[1])), "rb"))
    except:
        try:
            shotB = GTEDGE3_cli.runGT3(shotargs)
        except Exception as e:
            raise Exception(e)
        try:
            pass
#            with open("outputs/s%s.%s.dat"  % (str(shotid[1]), str(timeid[1])), "wb") as f:
#                pickle.dump(shotB, f)
        except Exception as e:
             print e
        try:
            f.close()
        except:
            pass

    #chiComp(shotA, shotB)
    #fluxComp(shotA, shotB)
    #convComp(shotA, shotB)
    #IOLcomp(shotA, shotB)
    heatsComp(shotA, shotB)
    #plasmaComp(shotA, shotB)
    #rotComp(shotA, shotB)
    plt.show(block=True)
