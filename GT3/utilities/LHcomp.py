#!/usr/bin/python
"""
Created on Thu Dec 28 02:13:04 2017

@author: Jonathan
"""
###############################################################################
#
#    L-H transition analysis and comparison module based on GT3
#
###############################################################################

from utilities import GTEDGE3_cli
import graphs as graphs
import matplotlib.pyplot as plt
import pickle
import numpy as np
from math import ceil, floor, log10

shotid, timeids = 164436, [3720, 3740]
#Sshotid, timeids = 118888, [1525, 1570]
runid=None


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

def chiComp(shotL, shotH):
    #prettyID = [r'$Q^{diff}_{j,r}=Q^{total}_j-Q^{conv}_j-Q^{heatin}_j-Q^{visc}_j$ w/ IOL',
    #            r'$Q^{diff}_{j,r}=Q^{total}_j-Q^{conv}_j-Q^{heatin}_j-Q^{visc}_j$ w/out IOL']
    prettyID = [r'$\chi_{j,r,L}$',
                r'$\chi_{j,r,H}$']

    title = r"L-H Mode Ion Heat Conductivity Comparison"

    adjustments = {0: 1.}
    caption = ["Main ion heat conductivity L-H mode comparison in the edge for DIII-D shot %s.%s, %s.%s w/IOL correction" % (str(shotid), str(timeids[0]), str(shotid), str(timeids[1])),
               "convective heat flow,  work done on plasma and visc. heat"]

    graphs.prettyCompare(('rhor', shotL.rtrans.rhor[:-1]), [(prettyID[0], shotL.rtrans.chi.chi.i.chi4[:-1]), (prettyID[1], shotH.rtrans.chi.chi.i.chi4[:-1])],
                         yrange=yRangeFind([shotL.rtrans.chi.chi.i.chi4[-30:], shotH.rtrans.chi.chi.i.chi4[-30:]]),
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

def fluxComp(shotL, shotH):
    prettyID = [r'$\Gamma_{r,j,L}$',
                r'$\Gamma_{r,j,H}$']

    caption = ["Ion radial particle flux for DIII-D shot %s.%s and %s.%s" % (str(shotid), str(timeids[0]), str(shotid), str(timeids[1]))]
    title = r"$\Gamma_{r,j}$ L-H transition"
    adjustments = {}
    graphs.prettyCompare(('rhor', shotL.rtrans.rhor), [(prettyID[0], shotL.rtrans.gamma_diff_D), (prettyID[1], shotH.rtrans.gamma_diff_D)],
                         #yrange=(yRangeFind([shotL.rtrans.gamma_diff_D[-30:], shotH.rtrans.gamma_diff_D[-30:]])),
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

def heatsComp(shotL, shotH):

    fig = plt.figure(figsize=(12, 8))
    fig.tight_layout()
    fig.suptitle(r'Heating sources')
    ax1 = fig.add_subplot(221)
    ax1.set_title(r'$Q^{conv}_i$', fontsize=16)
    ax1.set_ylabel(r"""$Q^{conv}_i$
                    $\left[\frac{W}{m^3}\right]$""", fontsize=16, rotation=0, ha='right')
    ax1.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax1.plot(shotL.rtrans.rhor[:-1][-30:], shotL.rtrans.chi.conv25[:-1][-30:], label="L mode")
    ax1.plot(shotH.rtrans.rhor[:-1][-30:], shotH.rtrans.chi.conv25[:-1][-30:], label="H mode")
    ax1.legend()

    ax2 = fig.add_subplot(222)
    ax2.set_title(r'$Q^{visc}_i$', fontsize=16)
    ax2.set_ylabel(r"""$Q^{visc}_i$
                    $\left[\frac{W}{m^3}\right]$""", fontsize=16, rotation=0, ha='right')
    ax2.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax2.plot(shotL.rtrans.rhor[:-1][-30:], shotL.rtrans.chi.heatvisc[:-1][-30:], label="L mode")
    ax2.plot(shotH.rtrans.rhor[:-1][-30:], shotH.rtrans.chi.heatvisc[:-1][-30:], label="H mode")
    ax2.legend()

    ax3 = fig.add_subplot(223)
    ax3.set_title(r'$Q^{heatin}_i$', fontsize=16)
    ax3.set_ylabel(r"""$Q^{heatin}_i$
                    $\left[\frac{W}{m^3}\right]$""", fontsize=16, rotation=0, ha='right')
    ax3.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax3.plot(shotL.rtrans.rhor[:-1][-30:], shotL.rtrans.chi.heatin[:-1][-30:], label="L mode")
    ax3.plot(shotH.rtrans.rhor[:-1][-30:], shotH.rtrans.chi.heatin[:-1][-30:], label="H mode")
    ax3.legend()

    ax4 = fig.add_subplot(224)
    ax4.set_title(r'$Q^{total}_i$', fontsize=16)
    ax4.set_ylabel(r"""$Q^{total}_i$
                    $\left[\frac{W}{m^3}\right]$""", fontsize=16, rotation=0, ha='right')
    ax4.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax4.plot(shotL.rtrans.rhor[:-1][-30:], shotL.rtrans.chi.Qi[:-1][-30:], label="L mode")
    ax4.plot(shotH.rtrans.rhor[:-1][-30:], shotH.rtrans.chi.Qi[:-1][-30:], label="H mode")
    ax4.legend()
def plasmaComp(shotL, shotH):

    fig = plt.figure(figsize=(12, 8))
    fig.tight_layout()
    fig.suptitle(r'Plasma profiles')
    ax1 = fig.add_subplot(231)
    ax1.set_title(r'$T_i$', fontsize=16)
    ax1.set_ylabel(r"""$T_i$
                    $[eV]$""", fontsize=16, rotation=0, ha='right')
    ax1.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax1.plot(shotL.rtrans.rhor[:-1][-30:], shotL.rtrans.profiles.T.i.ev[-30:], label="L mode")
    ax1.plot(shotH.rtrans.rhor[:-1][-30:], shotH.rtrans.profiles.T.i.ev[-30:], label="H mode")
    ax1.legend()

    ax2 = fig.add_subplot(232)
    ax2.set_title(r'$T_e$', fontsize=16)
    ax2.set_ylabel(r"""$T_e$
                    $[eV]$""", fontsize=16, rotation=0, ha='right')
    ax2.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax2.plot(shotL.rtrans.rhor[:-1][-30:], shotL.rtrans.profiles.T.e.ev[-30:], label="L mode")
    ax2.plot(shotH.rtrans.rhor[:-1][-30:], shotH.rtrans.profiles.T.e.ev[-30:], label="H mode")
    ax2.legend()

    ax3 = fig.add_subplot(233)
    ax3.set_title(r'$n_i$', fontsize=16)
    ax3.set_ylabel(r"""$n_i$
                    $\left[\frac{1}{m^3}\right]$""", fontsize=16, rotation=0, ha='right')
    ax3.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax3.set_ylim(yRangeFind([shotL.rtrans.profiles.n.i[-30:], shotH.rtrans.profiles.n.i[-30:]]))
    ax3.plot(shotL.rtrans.rhor[:-1][-30:], shotL.rtrans.profiles.n.i[-30:], label="L mode")
    ax3.plot(shotH.rtrans.rhor[:-1][-30:], shotH.rtrans.profiles.n.i[-30:], label="H mode")
    ax3.legend()

    ax4 = fig.add_subplot(234)
    ax4.set_title(r'$n_e$', fontsize=16)
    ax4.set_ylabel(r"""$n_e$
                    $\left[\frac{1}{m^3}\right]$""", fontsize=16, rotation=0, ha='right')
    ax4.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax4.plot(shotL.rtrans.rhor[:-1][-30:], shotL.rtrans.profiles.n.e[-30:], label="L mode")
    ax4.plot(shotH.rtrans.rhor[:-1][-30:], shotH.rtrans.profiles.n.e[-30:], label="H mode")
    ax4.legend()

    ax5 = fig.add_subplot(235)
    ax5.set_title(r'$n_n$', fontsize=16)
    ax5.set_ylabel(r"""$n_i$
                    $\left[\frac{1}{m^3}\right]$""", fontsize=16, rotation=0, ha='right')
    ax5.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax5.set_ylim(yRangeFind([shotL.rtrans.profiles.nn.tot[-30:], shotH.rtrans.profiles.nn.tot[-30:]]))
    ax5.plot(shotL.rtrans.rhor[:-1][-30:], shotL.rtrans.profiles.nn.tot[-30:], label="L mode")
    ax5.plot(shotH.rtrans.rhor[:-1][-30:], shotH.rtrans.profiles.nn.tot[-30:], label="H mode")
    ax5.legend()

    ax6 = fig.add_subplot(236)
    ax6.set_title(r'$E_r$', fontsize=16)
    ax6.set_ylabel(r"""$E_r$
                    $\left[\frac{V}{m}\right]$""", fontsize=16, rotation=0, ha='right')
    ax6.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax6.set_ylim(yRangeFind([shotL.core.E_r_fsa[-60:], shotH.core.E_r_fsa[-60:]]))
    ax6.plot(shotL.rtrans.rhor[:-1][-60:], shotL.core.E_r_fsa[-60:], label="L mode")
    ax6.plot(shotH.rtrans.rhor[:-1][-60:], shotH.core.E_r_fsa[-60:], label="H mode")
    ax6.legend()

def convComp(shotL, shotH):





    prettyID = [r'$Q^{conv}_{j,r,L}$',
                r'$Q^{conv}_{j,r,H}$']

    title = r"L-H Convective Heating Comparison"

    adjustments = {0: 1.}
    caption = ["Comparison of Convective heating for DIII-D shot %s.%s, %s.%s w/IOL correction" % (str(shotid), str(timeids[0]), str(shotid), str(timeids[1]))]
    graphs.prettyCompare(('rhor', shotL.rtrans.rhor[:-1]), [(prettyID[0], shotL.rtrans.chi.conv25[:-1]), (prettyID[1], shotH.rtrans.chi.conv25[:-1])],
                         yrange=yRangeFind([shotL.rtrans.chi.conv25[-30:], shotH.rtrans.chi.conv25[-30:]]),
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

def heatViscComp(shotL, shotH):
    prettyID = [r'$Q^{visc}_{j,r,L}$',
                r'$Q^{visc}_{j,r,H}$']

    caption = ["Comparison of Viscous heating for DIII-D shot %s.%s, %s.%s w/IOL correction" % (str(shotid), str(timeids[0]), str(shotid), str(timeids[1]))]
    title = r"$\Gamma_{r,j}$ L-H transition"
    adjustments = {}
    graphs.prettyCompare(('rhor', shotL.rtrans.rhor), [(prettyID[0], shotL.rtrans.chi.heatvisc), (prettyID[1], shotH.rtrans.chi.heatvisc)],
                         yrange=yRangeFind([shotL.rtrans.chi.heatvisc[-30:], shotH.rtrans.chi.heatvisc[-30:]]),
                         #yrange=[1E17, 1E20],
                         datalabels=[prettyID[0], prettyID[1]],
                         title=title,
                         ylabel=[r"$Q^{conv}_j$", r'$\left[\frac{W}{m^3}\right]$'],
                         caption=caption,
                         adjust=adjustments,
                         textSpace=.02,
                         yLabelAdj=-.8,
                         xLabelAdj=-0.025,
                         marginLeft=0.15,
                         marginBottom=0.2,
                         size=(16, 19),
                         legend=True)

def IOLcomp(shotL, shotH):

    fig = plt.figure(figsize=(12, 8))
    fig.tight_layout()
    fig.suptitle(r'Ion orbit loss Profiles')

    ax1 = fig.add_subplot(231)
    ax1.set_title(r'$F_{iol}$', fontsize=16)
    ax1.set_ylabel(r"""$F_{iol}$
                    """, fontsize=16, rotation=0, ha='right')
    ax1.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax1.set_ylim(0,1.)
    ax1.plot(shotL.rtrans.rhor[:-1][-30:], shotL.iol.forb_d_therm_1D[-30:], label="L mode")
    ax1.plot(shotH.rtrans.rhor[:-1][-30:], shotH.iol.forb_d_therm_1D[-30:], label="H mode")
    ax1.legend()

    ax2 = fig.add_subplot(232)
    ax2.set_title(r'$E_{iol}$', fontsize=16)
    ax2.set_ylabel(r"""$E_{iol}$
                    """, fontsize=16, rotation=0, ha='right')
    ax2.set_ylim(0,1.)
    ax2.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax2.plot(shotL.rtrans.rhor[:-1][-30:], shotL.iol.eorb_d_therm_1D[-30:], label="L mode")
    ax2.plot(shotH.rtrans.rhor[:-1][-30:], shotH.iol.eorb_d_therm_1D[-30:], label="H mode")
    ax2.legend()

    ax3 = fig.add_subplot(233)
    ax3.set_title(r'$E_{iol}$', fontsize=16)
    ax3.set_ylabel(r"""$M_{iol}$
                    """, fontsize=16, rotation=0, ha='right')
    ax3.set_ylim(-1,1.)
    ax3.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax3.plot(shotL.rtrans.rhor[:-1][-30:], shotL.iol.morb_d_therm_1D[-30:], label="L mode")
    ax3.plot(shotH.rtrans.rhor[:-1][-30:], shotH.iol.morb_d_therm_1D[-30:], label="H mode")
    ax3.legend()

    ax4 = fig.add_subplot(234)
    ax4.set_title(r'$F^C_{iol}$', fontsize=16)
    ax4.set_ylabel(r"""$F_{iol}$
                    """, fontsize=16, rotation=0, ha='right')
    ax4.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax4.set_ylim(0, 1.)
    ax4.plot(shotL.rtrans.rhor[:-1][-30:], shotL.iol.forb_c_therm_1D[-30:], label="L mode")
    ax4.plot(shotH.rtrans.rhor[:-1][-30:], shotH.iol.forb_c_therm_1D[-30:], label="H mode")
    ax4.legend()

    ax5 = fig.add_subplot(235)
    ax5.set_title(r'$E^C_{iol}$', fontsize=16)
    ax5.set_ylabel(r"""$E_{iol}$
                    """, fontsize=16, rotation=0, ha='right')
    ax5.set_ylim(0, 1.)
    ax5.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax5.plot(shotL.rtrans.rhor[:-1][-30:], shotL.iol.eorb_c_therm_1D[-30:], label="L mode")
    ax5.plot(shotH.rtrans.rhor[:-1][-30:], shotH.iol.eorb_c_therm_1D[-30:], label="H mode")
    ax5.legend()

    ax6 = fig.add_subplot(236)
    ax6.set_title(r'$M^C_{iol}$', fontsize=16)
    ax6.set_ylabel(r"""$E_{iol}$
                    """, fontsize=16, rotation=0, ha='right')
    ax6.set_ylim(-1, 1.)
    ax6.set_xlabel(r'$\rho$', fontsize=16, labelpad=-10)
    ax6.plot(shotL.rtrans.rhor[:-1][-30:], shotL.iol.morb_c_therm_1D[-30:], label="L mode")
    ax6.plot(shotH.rtrans.rhor[:-1][-30:], shotH.iol.morb_c_therm_1D[-30:], label="H mode")
    ax6.legend()

if __name__ == "__main__":

    #   L-H Mode senstivity analysis

    # L mode

    shotargs = {'shotid': shotid,
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
        shotL = pickle.load(open("outputs/s%s.%s.dat" % (str(shotid), str(timeids[0])), "rb"))
    except:
        try:
            shotL = GTEDGE3_cli.runGT3(shotargs)
        except Exception as e:
            raise Exception(e)
        try:
            pass
#            with open("outputs/s%s.%s.dat"  % (str(shotid), str(timeid[0])), "wb") as f:
#                pickle.dump(shotL, f)
        except Exception as e:
             print e
        try:
            f.close()
        except:
            pass

    # H-mode

    shotargs = {'shotid': shotid,
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
        shotH = pickle.load(open("outputs/s%s.%s.dat" % (str(shotid), str(timeids[1])), "rb"))
    except:
        try:
            shotH = GTEDGE3_cli.runGT3(shotargs)
        except Exception as e:
            raise Exception(e)
        try:
            pass
#            with open("outputs/s%s.%s.dat"  % (str(shotid), str(timeid[1])), "wb") as f:
#                pickle.dump(shotH, f)
        except Exception as e:
             print e
        try:
            f.close()
        except:
            pass

    chiComp(shotL, shotH)
    fluxComp(shotL, shotH)
    #convComp(shotL, shotH)
    #heatViscComp(shotL, shotH)
    IOLcomp(shotL, shotH)
    heatsComp(shotL, shotH)
    plasmaComp(shotL, shotH)
    plt.show(block=True)
