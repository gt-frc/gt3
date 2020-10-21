#!/usr/bin/python
"""
Created on Thu Dec 28 02:13:04 2017

@author: Jonathan
"""
###############################################################################
#
#    sensitivity study of IOL and neutrals
#
###############################################################################

from . import graphs as graphs
import matplotlib.pyplot as plt
import pickle
import numpy as np
from math import ceil, floor, log10
#164436_3720 good shot
shotid=164436

timeid=3720

runid='j1900'


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
    a=int(flat[2]*.75)
    b=int(flat[-2]*1.25)
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
#    Chi Comparison
#
###############################################################################

def chiComp1(shot, shotnoIOL):
    prettyID = [r'$Q^{diff}_{j,r}=Q^{total}_j$ w/ IOL',
                r'$Q^{diff}_{j,r}=Q^{total}_j$ w/out IOL']

    title = r"Comparison of Ion Heat Conductivity w/ and w/out IOL"

    adjustments = {}
    caption = ["Ion heat conductivity in the edge for DIII-D shot %s.%s w/ and w/out IOL correction" % (str(shotid), str(timeid)),
               "Neutrals calculation included"]

    graphs.prettyCompare(('rhor', shot.rtrans.rhor), [(prettyID[0], shot.rtrans.chi.chi.i.chi1), (prettyID[1], shotnoIOL.rtrans.chi.chi.i.chi1)],
                         yrange=yRangeFind([shot.rtrans.chi.chi.i.chi1[-30:], shotnoIOL.rtrans.chi.chi.i.chi1[-30:]]),
                         #yrange=[0., 15],
                         datalabels=[prettyID[0], prettyID[1]],
                         title=title,
                         ylabel=[r"$\chi_j$", r'$\left[\frac{m^2}{s}\right]$'],
                         caption=caption,
                         adjust=adjustments,
                         xTickremove=None,
                         textSpace=.065,
                         marginBottom=.2,
                         marginLeft=0.095,
                         capAdj=0.15,
                         yLabelAdj=-.25,
                         xLabelAdj=-.05,
                         size=(16, 12))


def chiComp2(shot, shotnoIOL):
    prettyID = [r'$Q^{diff}_{j,r}=Q^{total}_j-Q^{conv}_j$ w/ IOL',
                r'$Q^{diff}_{j,r}=Q^{total}_j-Q^{conv}_j$ w/out IOL']

    title = r"Comparison of Ion Heat Conductivity w/ and w/out IOL correction"

    adjustments = {0: 0.}
    caption = ["Main ion heat conductivity in the edge for DIII-D shot %s.%s w/ and w/out IOL" % (str(shotid), str(timeid)),
               "and convective heat flow, neutrals calculation included"]

    graphs.prettyCompare(('rhor', shot.rtrans.rhor), [(prettyID[0], shot.rtrans.chi.chi.i.chi2), (prettyID[1], shotnoIOL.rtrans.chi.chi.i.chi2)],
                         yrange=yRangeFind([shot.rtrans.chi.chi.i.chi2[-30:], shotnoIOL.rtrans.chi.chi.i.chi2[-30:]]),
                         #yrange=[0., 15],
                         datalabels=[prettyID[0], prettyID[1]],
                         title=title,
                         ylabel=[r"$\chi_j$", r'$\left[\frac{m^2}{s}\right]$'],
                         caption=caption,
                         adjust=adjustments,
                         xTickremove=None,
                         textSpace=.065,
                         marginBottom=.2,
                         marginLeft=0.095,
                         capAdj=0.15,
                         yLabelAdj=-.25,
                         xLabelAdj=-.05,
                         size=(16, 12))


def chiComp3(shot, shotnoIOL):
    prettyID = [r'$Q^{diff}_{j,r}=Q^{total}_j-Q^{conv}_j-Q^{heatin}_j$ w/ IOL',
                r'$Q^{diff}_{j,r}=Q^{total}_j-Q^{conv}_j-Q^{heatin}_j$ w/out IOL']

    title = r"Comparison of Ion Heat Conductivity w/ and w/out IOL correction"

    adjustments = {0: 1.}
    caption = ["Main ion heat conductivity in the edge for DIII-D shot %s.%s w/ and w/out IOL" % (str(shotid), str(timeid)),
               "convective heat flow, and work done on plasma, neutrals calculation included"]

    graphs.prettyCompare(('rhor', shot.rtrans.rhor), [(prettyID[0], shot.rtrans.chi.chi.i.chi3), (prettyID[1], shotnoIOL.rtrans.chi.chi.i.chi3)],
                         yrange=yRangeFind([shot.rtrans.chi.chi.i.chi3[-30:], shotnoIOL.rtrans.chi.chi.i.chi3[-30:]]),
                         datalabels=[prettyID[0], prettyID[1]],
                         title=title,
                         ylabel=[r"$\chi_j$", r'$\left[\frac{m^2}{s}\right]$'],
                         caption=caption,
                         adjust=adjustments,
                         xTickremove=None,
                         textSpace=.065,
                         marginBottom=.2,
                         marginLeft=0.095,
                         capAdj=-0.15,
                         yLabelAdj=-.25,
                         xLabelAdj=-.05,
                         size=(16, 12))


def chiComp4(shot, shotnoIOL):
    #prettyID = [r'$Q^{diff}_{j,r}=Q^{total}_j-Q^{conv}_j-Q^{heatin}_j-Q^{visc}_j$ w/ IOL',
    #            r'$Q^{diff}_{j,r}=Q^{total}_j-Q^{conv}_j-Q^{heatin}_j-Q^{visc}_j$ w/out IOL']
    prettyID = [r'$\chi_{j,r}$ w/ IOL correction',
                r'$\chi_{j,r}$ w/out IOL correction']

    title = r"Comparison of Ion Heat Conductivity w/ and w/out IOL correction"

    adjustments = {0: 1.}
    caption = ["Main ion heat conductivity in the edge for DIII-D shot %s.%s w/ and w/out IOL correction" % (str(shotid), str(timeid)),
               "convective heat flow,  work done on plasma and visc. heat"]

    graphs.prettyCompare(('rhor', shot.rtrans.rhor[:-1]), [(prettyID[0], shot.rtrans.chi.chi.i.chi4[:-1]), (prettyID[1], shotnoIOL.rtrans.chi.chi.i.chi4[:-1])],
                         yrange=yRangeFind([shot.rtrans.chi.chi.i.chi4[-30:], shotnoIOL.rtrans.chi.chi.i.chi4[-30:]]),
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
                         size=(16, 12))


def chieComp(shot, shotnoIOL):
    prettyID = [r'$\chi_e$ w/ IOL',
                r'$\chi_e$ w/out IOL']

    title = r"Comparison of Electron Heat Conductivity w/ and w/out IOL correction"

    adjustments = {0: 1.}
    caption = ["Main electron heat conductivity in the edge for DIII-D shot %s.%s w/ and w/out IOL correction" % (str(shotid), str(timeid)),
               "neutrals calculation included"]

    graphs.prettyCompare(('rhor', shot.rtrans.rhor), [(prettyID[0], shot.rtrans.chi.chi.e), (prettyID[1], shotnoIOL.rtrans.chi.chi.e)],
                         yrange=yRangeFind([shot.rtrans.chi.chi.e[-30:], shotnoIOL.rtrans.chi.chi.e[-30:]]),
                         #yrange=[0., 50],
                         datalabels=[prettyID[0], prettyID[1]],
                         title=title,
                         ylabel=[r"$\chi_e$", r'$\left[\frac{m^2}{s}\right]$'],
                         caption=caption,
                         adjust=adjustments,
                         xTickremove=None,
                         textSpace=.02,
                         marginBottom=.25,
                         marginLeft=0.095,
                         capAdj=-.15,
                         xLabelAdj=-0.05,
                         yLabelAdj=-0.,
                         size=(16, 12))


###############################################################################
#
#   Deactivate neutrals
#
###############################################################################

def chiComp1neuts(shot, shotnoneuts):
    prettyID = [r'$Q^{diff}_{i}=Q^{total}_j$ w/ neutrals',
                r'$Q^{diff}_{i}=Q^{total}_j$ w/out neutrals']

    title = r"Comparison of Ion Heat Conductivity w/ and w/out neutrals"

    adjustments = {0: -1.}
    caption = ["Ion heat conductivity in the edge for DIII-D shot %s.%s w/ and w/out neutrals, IOL corrected" % (str(shotid), str(timeid))]

    graphs.prettyCompare(('rhor', shot.rtrans.rhor), [(prettyID[0], shot.rtrans.chi.chi.i.chi1), (prettyID[1], shotnoneuts.rtrans.chi.chi.i.chi1)],
                         yrange=yRangeFind([shot.rtrans.chi.chi.i.chi1[-30:], shotnoneuts.rtrans.chi.chi.i.chi1[-30:]]),
                         #yrange=[0, 15.],
                         datalabels=[prettyID[0], prettyID[1]],
                         title=title,
                         ylabel=[r"$\chi_j$", r'$\left[\frac{m^2}{s}\right]$'],
                         caption=caption,
                         adjust=adjustments,
                         xTickremove=None,
                         textSpace=.055,
                         #                    marginBottom=.15,
                         marginLeft=.1,
                         yLabelAdj=-.25,
                         size=(16, 12))


def chiComp2neuts(shot, shotnoneuts):
    prettyID = [r'$Q^{diff}_{i}=Q^{total}_j-Q^{conv}_j$ w/ neutrals',
                r'$Q^{diff}_{i}=Q^{total}_j-Q^{conv}_j$ w/out neutrals']

    title = r"Comparison of Ion Heat Conductivity w/ and w/out neutrals"

    adjustments = {0: -.5}
    caption = ["Main ion heat conductivity in the edge for DIII-D shot %s.%sw/ and w/out neutrals" % (str(shotid), str(timeid)),
               "IOL corrected, convective heat subtracted"]

    graphs.prettyCompare(('rhor', shot.rtrans.rhor), [(prettyID[0], shot.rtrans.chi.chi.i.chi2), (prettyID[1], shotnoneuts.rtrans.chi.chi.i.chi2)],
                         yrange=yRangeFind([shot.rtrans.chi.chi.i.chi2[-30:], shotnoneuts.rtrans.chi.chi.i.chi2[-30:]]),
                         #yrange=[0, 15.],
                         datalabels=[prettyID[0], prettyID[1]],
                         title=title,
                         ylabel=[r"$\chi_j$", r'$\left[\frac{m^2}{s}\right]$'],
                         caption=caption,
                         adjust=adjustments,
                         xTickremove=None,
                         textSpace=.055,
                         marginBottom=.15,
                         marginLeft=0.095,
                         yLabelAdj=-.25,
                         size=(16, 12))


def chiComp3neuts(shot, shotnoneuts):
    prettyID = [r'$Q^{diff}_{j,r}=Q^{total}_j-Q^{conv}_j-Q^{heatin}_j$ w/ neutrals',
                r'$Q^{diff}_{j,r}=Q^{total}_j-Q^{conv}_j-Q^{heatin}_j$ w/out neutrals']

    title = r"Comparison of Ion Heat Conductivity w/ and w/out neutrals"

    adjustments = {0: -1.}
    caption = ["Main ion heat conductivity in the edge for DIII-D shot %s.%s w/ and w/out neutrals" % (str(shotid), str(timeid)),
               "convective heat flow, and work done on plasma, IOL corrected"]

    graphs.prettyCompare(('rhor', shot.rtrans.rhor), [(prettyID[0], shot.rtrans.chi3), (prettyID[1], shotnoneuts.rtrans.chi.chi.chi3)],
                         yrange=yRangeFind([shot.rtrans.chi.chi.i.chi3[-30:], shotnoneuts.rtranschi.chi.i.chi3[-30:]]),
                         datalabels=[prettyID[0], prettyID[1]],
                         title=title,
                         ylabel=[r"$\chi_j$", r'$\left[\frac{m^2}{s}\right]$'],
                         caption=caption,
                         adjust=adjustments,
                         xTickremove=None,
                         textSpace=.065,
                         marginBottom=.15,
                         marginLeft=0.095,
                         yLabelAdj=-.25,
                         size=(16, 12))


def chiComp4neuts(shot, shotnoneuts):
    prettyID = [r'$Q^{diff}_{j,r}=Q^{total}_j-Q^{conv}_j-Q^{heatin}_j-Q^{visc}_j$ w/ neutrals',
                r'$Q^{diff}_{j,r}=Q^{total}_j-Q^{conv}_j-Q^{heatin}_j-Q^{visc}_j$ w/out neutrals']

    title = r"Comparison of Ion Heat Conductivity w/ and w/out neutrals"

    adjustments = {0: -1.}
    caption = ["Main ion heat conductivity in the edge for DIII-D shot %s.%s w/ and w/out neutrals" % (str(shotid), str(timeid)),
               "convective heat flow,  work done on plasma and visc., IOL corrected"]

    graphs.prettyCompare(('rhor', shot.rtrans.rhor), [(prettyID[0], shot.rtrans.chi.chi.i.chi4), (prettyID[1], shotnoneuts.rtrans.chi.chi.i.chi4)],
                         yrange=yRangeFind([shot.rtrans.chi.chi.i.chi4[-30:], shotnoneuts.rtrans.chi.chi.i.chi4[-30:]]),
                         datalabels=[prettyID[0], prettyID[1]],
                         title=title,
                         ylabel=[r"$\chi_j$", r'$\left[\frac{m^2}{s}\right]$'],
                         caption=caption,
                         adjust=adjustments,
                         xTickremove=None,
                         textSpace=.065,
                         marginBottom=.2,
                         marginLeft=0.095,
                         capAdj=0.1,
                         xLabelAdj=-0.05,
                         yLabelAdj=-.25,
                         size=(16, 12))

def paperComp(shotIOL, shot):
    prettyID = [r'$\chi_{j,r}$ no corr',
                r'$\chi_{j,r}$ w/ IOL',
                r'$\chi_{j,r}$ w/ IOL, con',
                r'$\chi_{j,r}$ w/ IOL, conv, int']

    title = r"Comparison of Ion Heat Conductivity with various corrections"

    adjustments = {0: -1.}
    caption = ["Main ion heat conductivity in the edge for DIII-D shot %s.%s with no corrections (red)" % (str(shotid), str(timeid)),
               "corrected for IOL (blue), corrected for convective heat flux (green), and corrected for work done on plasma (purple)"]

    graphs.prettyCompare(('rhor', shot.rtrans.rhor[:-2]), [(prettyID[0], shot.rtrans.chi.chi.i.chi1[:-2]), (prettyID[1], shotIOL.rtrans.chi.chi.i.chi1[:-2]), (prettyID[2], shotIOL.rtrans.chi.chi.i.chi2[:-2]),
        (prettyID[3], shotIOL.rtrans.chi.chi.i.chi3[:-2])],
                         #yrange=yRangeFind([shot.rtrans.chi.chi.i.chi1[-30:], shotIOL.rtrans.chi.chi.i.chi3[-30:]]),
                         yrange=[-1.5, 3.5],
                         datalabels=[prettyID[0], prettyID[1], prettyID[2], prettyID[3]],
                         title=title,
                         ylabel=[r"$\chi_j$", r'$\left[\frac{m^2}{s}\right]$'],
                         caption=caption,
                         adjust=adjustments,
                         xTickremove=None,
                         textSpace=.065,
                         marginBottom=.15,
                         marginLeft=0.095,
                         yLabelAdj=-.25,
                         size=(16, 12),
                         legend=True)



    prettyID = [r'$D_{drag,j}$']

    title = r"Calculation of $D$"

    adjustments = {0: -1.}
    caption = ["Main ion particle diffusion coefficient in the edge for DIII-D shot %s.%s with no corrections (red)" % (str(shotid), str(timeid))]

    graphs.prettyCompare(('rhor', shot.rtrans.nu_drag_D[:-2]), [(prettyID[0], shotIOL.rtrans.nu_drag_D[:-2])],
                         #yrange=yRangeFind([shot.rtrans.chi.chi.i.chi1[-30:], shotIOL.rtrans.chi.chi.i.chi3[-30:]]),
                         yrange=[-15, 35],
                         datalabels=[prettyID[0]],
                         title=title,
                         ylabel=[r"$D_{j}$", r'$\left[\frac{m^2}{s}\right]$'],
                         caption=caption,
                         adjust=adjustments,
                         xTickremove=None,
                         textSpace=.065,
                         marginBottom=.15,
                         marginLeft=0.095,
                         yLabelAdj=-.25,
                         size=(16, 12),
                         legend=True)

def chieCompneuts(shot, shotnoneuts):
    prettyID = [r"$\chi_e$ w/ neutrals",
                r"$\chi_e$ w/out neutrals"]

    title = r"Comparison of Electron Heat Conductivity w/ and w/out neutrals"

    adjustments = {0: -1.}
    caption = ["Main electron heat conductivity in the edge for DIII-D shot %s.%s w/ and w/out neutrals" % (str(shotid), str(timeid)),
               "IOL corrected"]

    graphs.prettyCompare(('rhor', shot.rtrans.rhor), [(prettyID[0], shot.rtrans.chi.chi.e), (prettyID[1], shotnoneuts.rtrans.chi.chi.e)],
                         yrange=yRangeFind([shot.rtrans.chi.chi.e[-30:], shotnoneuts.rtrans.chi.chi.e[-30:]]),
                         #yrange=[0, 50.],
                         datalabels=[prettyID[0], prettyID[1]],
                         title=title,
                         ylabel=[r"$\chi_e$", r'$\left[\frac{m^2}{s}\right]$'],
                         caption=caption,
                         adjust=adjustments,
                         xTickremove=None,
                         textSpace=.02,
                         marginBottom=.2,
                         marginLeft=0.095,
                         xLabelAdj=-0.025,
                         capAdj=0.15,
                         yLabelAdj=-.05,
                         size=(16, 12))

    #
    # prettyID = [r"Term 1",
    #             r"Term 2"]
    #
    # title = r"Comparison of Electron Heat Conductivity Terms"
    #
    # adjustments = {}
    # caption = ["Terms in electron heat conductivity in the edge for DIII-D shot %s.%s" % (str(shotid), str(timeid)),
    #            "IOL corrected. Term1 = 'diffusive' piece, Term2 = 'convective' piece"]
    #
    # term1 = [a / (b * shot.rtrans.xk * c) for a, b, c in zip(shot.rtrans.qHeate, shot.rtrans.xne, shot.rtrans.xte)]
    # term2 = [-1.5 * a / b for a, b in zip(shot.rtrans.gameltemp, shot.rtrans.xne)]
    #
    # graphs.prettyCompare(('rhor', shot.rtrans.rhor), [(prettyID[0], term1), (prettyID[1], term2)],
    #                      yrange=[-10, 50],
    #                      datalabels=[prettyID[0], prettyID[1]],
    #                      title=title,
    #                      ylabel=[r"$\frac{\chi_e}{L_{T_e}}$", r'$\left[\frac{m^2}{s}\right]$'],
    #                      caption=caption,
    #                      adjust=adjustments,
    #                      xTickremove=None,
    #                      textSpace=.02,
    #                      marginBottom=.2,
    #                      marginLeft=0.095,
    #                      xLabelAdj=-0.025,
    #                      capAdj=0.15,
    #                      yLabelAdj=-.05,
    #                      size=(16, 12))


def fluxComp(shot, shotnoIOL, shotnoneuts):
    prettyID = [r'$\Gamma_{r,j}$ w/ neutrals',
                r'$\Gamma_{r,j}$ w/out neutrals']

    caption = ["Ion radial particle flux for DIII-D shot %s.%s w/ and w/out neutrals" % (str(shotid), str(timeid))]
    title = r"$\Gamma_{r,j}$  w/ and w/out neutrals"
    adjustments = {}
    graphs.prettyCompare(('rhor', shot.rtrans.rhor), [(prettyID[0], shot.rtrans.gamma_diff_D), (prettyID[1], shotnoneuts.rtrans.gamma_diff_D)],
                         #yrange=(yRangeFind([shot.rtrans.gamma_diff_D[-30:], shotnoneuts.rtrans.gamma_diff_D[-30:]])),
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
                         size=(16, 19))

    prettyID = [r'$\Gamma_{r,j}$ w/ IOL corr.',
                r'$\Gamma_{r,j}$ w/out IOL corr']

    caption = ["Ion radial particle flux for DIII-D shot %s.%s w/ and w/out IOL Correction" % (str(shotid), str(timeid))]
    title = r"$\Gamma_{r,j}$  w/ and w/out IOL corr"
    adjustments = {}
    graphs.prettyCompare(('rhor', shot.rtrans.rhor), [(prettyID[0], shot.rtrans.gamma_diff_D), (prettyID[1], shotnoIOL.rtrans.gamma_diff_D)],
                         #yrange=yRangeFind([shot.rtrans.gamma_diff_D[-30:], shotnoIOL.rtrans.gamma_diff_D[-30:]]),
                         yrange=[0.,4E20],
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
                         size=(16, 19))


if __name__ == "__main__":

    #   IOL Sensitivity

    shotargs = {'shotid': shotid,
                'timeid': timeid,
                'runid': runid,
                'nbRun': True,
                'IOL': True,
                'quiet': False,
                'reNeut': False,
                'neutrals': True,
                'gt3Method': 'radialtrans',
                'debug' : True}
    try:
        shot = pickle.load(open("outputs/s%s.%s.dat" % (str(shotid), str(timeid)), "rb"))
    except:
        try:
            shot = GTEDGE3_cli.runGT3(shotargs)
        except BaseException as e:
            raise(e)
        try:
            pass
#            with open("outputs/s%s.%s.dat"  % (str(shotid), str(timeid)), "wb") as f:
#                pickle.dump(shot, f)
        except BaseException as e:
             raise(e)
        try:
            f.close()
        except:
            pass


    shotargs = {'shotid': shotid,
                'timeid': timeid,
                'runid': runid,
                'nbRun': True,
                'IOL': False,
                'quiet': True,
                'reNeut': False,
                'neutrals': True,
                'gt3Method': 'radialtrans',
                'debug' : False}

    try:
        shotnoIOL = pickle.load(open("outputs/s%s.%s.noIOL.dat" % (str(shotid), str(timeid)), "rb"))
    except:
        try:
            shotnoIOL = GTEDGE3_cli.runGT3(shotargs)
        except Exception as e:
            raise Exception(e)
        try:
            pass
#            with open("outputs/s%s.%s.noIOL.dat" % (str(shotid), str(timeid)), "wb") as f:
#                pickle.dump(shot, f)
        except Exception as e:
            print(e)
        try:
            f.close()
        except:
            pass

#    chiComp1(shot,shotnoIOL)
#    chiComp2(shot,shotnoIOL)
    # chiComp3(shot,shotnoIOL)
    #chiComp4(shot, shotnoIOL)

    #chieComp(shot, shotnoIOL)
    #    qComp1(shot,shotnoIOL)
    #    qComp2(shot,shotnoIOL)
    #    qComp3(shot,shotnoIOL)

    #   Neutrals sensitivity

    shotargs = {'shotid': shotid,
                'timeid': timeid,
                'runid': runid,
                'nbRun': True,
                'IOL': True,
                'quiet': True,
                'reNeut': False,
                'neutrals': False,
                'gt3Method': 'radialtrans',
                'debug' : False}
    try:
        shotnoneuts = pickle.load(open("outputs/s%s.%s.noneuts.dat"  % (str(shotid), str(timeid)), "rb"))
    except:
        try:
            shotnoneuts = GTEDGE3_cli.runGT3(shotargs)
        except Exception as e:
            raise Exception(e)
        try:
            pass
#            with open("outputs/s%s.%s.noneuts.dat" % (str(shotid), str(timeid)), "wb") as f:
#                pickle.dump(shot, f)
        except Exception as e:
            print(e)
        try:
            f.close()
        except:
            pass

#    chiComp1neuts(shot, shotnoneuts)
#    chiComp2neuts(shot, shotnoneuts)
    #    chiComp3neuts(shot,shotnoneuts)
    #chiComp4neuts(shot, shotnoneuts)



    #    qComp1neuts(shot,shotnoneuts)
    #    qComp2neuts(shot,shotnoneuts)
    #    qComp3neuts(shot,shotnoneuts)
    #    qComp4neuts(shot,shotnoIOL)
    #
    #qieComp1(shot, shotnoIOL)
    paperComp(shot, shotnoIOL)
    fluxComp(shot, shotnoIOL, shotnoneuts)

    ####################################################################################
    #
    #
    #   TODO: Electron chi
    #
    ####################################################################################
#    chieCompneuts(shot, shotnoneuts.rtrans.)
    plt.show(block=True)