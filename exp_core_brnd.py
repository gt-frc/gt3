#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 18 13:22:31 2018

@author: max
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
# import matplotlib._cntr as cntr
from scipy.interpolate import griddata, UnivariateSpline, interp1d, interp2d, Rbf
from scipy.constants import elementary_charge
from shapely.geometry import Point, MultiPoint, LineString, LinearRing
from shapely.ops import polygonize, linemerge
import sys
from math import atan2, pi, ceil, sin, cos
from collections import namedtuple
from contours.quad import QuadContourGenerator
from scipy import constants

e = constants.elementary_charge

def calc_svfus(T, mode='dd'):
    def sigv(T, mode):  # function takes T in kev
        if mode == 'dt':
            B_G = 34.3827
            m_rc2 = 1124656

            C1 = 1.17302E-9
            C2 = 1.51361E-2
            C3 = 7.51886E-2
            C4 = 4.60643E-3
            C5 = 1.35000E-2
            C6 = -1.06750E-4
            C7 = 1.36600E-5

            theta = T / (1.0 - (T * (C2 + T * (C4 + T * C6))) / (1.0 + T * (C3 + T * (C5 + T * C7))))
            xi = (B_G ** 2.0 / (4.0 * theta)) ** (1.0 / 3.0)
            sigv = C1 * theta * np.sqrt(xi / (m_rc2 * T ** 3.0)) * np.exp(-3.0 * xi)
            sigv = sigv / 1.0E6  # convert from cm^3/s to m^3/s

        elif mode == 'dd':

            B_G = 31.3970
            m_rc2 = 937814

            # first for the D(d, p)T reaction
            C1_1 = 5.65718E-12
            C2_1 = 3.41267E-3
            C3_1 = 1.99167E-3
            C4_1 = 0.0
            C5_1 = 1.05060E-5
            C6_1 = 0.0
            C7_1 = 0.0

            theta_1 = T / (
                        1.0 - (T * (C2_1 + T * (C4_1 + T * C6_1))) / (1.0 + T * (C3_1 + T * (C5_1 + T * C7_1))))
            xi_1 = (B_G ** 2.0 / (4.0 * theta_1)) ** (1.0 / 3.0)
            sigv_1 = C1_1 * theta_1 * np.sqrt(xi_1 / (m_rc2 * T ** 3.0)) * np.exp(-3.0 * xi_1)

            # then for the D(d, n)He3 reaction

            C1_2 = 5.43360E-12
            C2_2 = 5.85778E-3
            C3_2 = 7.68222E-3
            C4_2 = 0.0
            C5_2 = -2.96400E-6
            C6_2 = 0.0
            C7_2 = 0.0

            theta_2 = T / (
                        1.0 - (T * (C2_2 + T * (C4_2 + T * C6_2))) / (1.0 + T * (C3_2 + T * (C5_2 + T * C7_2))))
            xi_2 = (B_G ** 2.0 / (4.0 * theta_2)) ** (1.0 / 3.0)
            sigv_2 = C1_2 * theta_2 * np.sqrt(xi_2 / (m_rc2 * T ** 3.0)) * np.exp(-3.0 * xi_2)

            sigv = (0.5 * sigv_1 + 0.5 * sigv_2) / 1.0E6  # convert from cm^3/s to m^3/s
        return sigv

    # create logspace over the relevant temperature range
    # (bosch hale technically only valid over 0.2 - 100 kev)
    Ti_range = np.logspace(-1, 2, 1000)  # values in kev
    sigv_fus_range = sigv(Ti_range, mode=mode)  # in m^3/s
    sigv_fus_interp = UnivariateSpline(Ti_range * 1.0E3 * e, sigv_fus_range, s=0)  # converted to Joules
    sv_fus = sigv_fus_interp(T.i.J)
    dsv_fus_dT = sigv_fus_interp.derivative()(T.i.J)

    return sv_fus, dsv_fus_dT


def calc_svrec_st(n, T):
    # # TODO: check this calculation. -MH
    # znint = np.array([16, 18, 20, 21, 22])
    # Tint = np.array([-1, 0, 1, 2, 3])
    #
    # rec = np.array([[-1.7523E+01, -1.6745E+01, -1.5155E+01, -1.4222E+01, -1.3301E+01],
    #                 [-1.8409E+01, -1.8398E+01, -1.8398E+01, -1.7886E+01, -1.7000E+01],
    #                 [-1.9398E+01, -1.9398E+01, -1.9398E+01, -1.9398E+01, -1.9398E+01],
    #                 [-2.0155E+01, -2.0155E+01, -2.0155E+01, -2.0155E+01, -2.0155E+01],
    #                 [-2.1000E+01, -2.1000E+01, -2.1000E+01, -2.1000E+01, -2.1000E+01]])
    #
    # interp1 = interp2d(znint, Tint, rec)
    #
    # zni_exps = np.linspace(16, 22, 100)
    # Ti_exps = np.linspace(-1, 3, 100)
    # svrec_vals = 10.0 ** (interp1(zni_exps, Ti_exps))  # in m^3/s
    #
    # zni_vals = np.logspace(16, 22, 100)
    # Ti_vals = np.logspace(-1, 3, 100) * e  # in joules
    #
    # dsvrec_dTi_vals = np.gradient(svrec_vals, Ti_vals, axis=0)
    #
    # zni_vals2d, Ti_vals2d = np.meshgrid(zni_vals, Ti_vals)
    #
    # zni_mod = np.where(n.i > 1E22, 1E22, n.i)
    # zni_mod = np.where(n.i < 1E16, 1E16, zni_mod)
    # Ti_mod = np.where(T.i.ev > 1E3, 1E3 * e, T.i.ev * e)
    # Ti_mod = np.where(T.i.ev < 1E-1, 1E-1 * e, Ti_mod)
    #
    # plt.semilogx(zni_vals2d.flatten(), Ti_vals2d.flatten())
    # plt.show()
    # print np.column_stack((zni_vals2d.flatten(), Ti_vals2d.flatten()))
    # sys.exit()
    # sv_rec = griddata(np.column_stack((zni_vals2d.flatten(), Ti_vals2d.flatten())),
    #                   svrec_vals.flatten(),
    #                   (zni_mod, Ti_mod),
    #                   method='linear', rescale=False)
    #
    # dsv_rec_dT = griddata(np.column_stack((zni_vals2d.flatten(), Ti_vals2d.flatten())),
    #                            dsvrec_dTi_vals.flatten(),
    #                            (zni_mod, Ti_mod),
    #                            method='linear', rescale=False)

    # return sv_rec, dsv_rec_dT
    return 0, 0


def calc_svcx_st(T):
    tint = np.array([-1, 0, 1, 2, 3])
    tnnt = np.array([0, 1, 2])

    cx = np.array([[-1.4097E+01, -1.3921E+01, -1.3553E+01, -1.4097E+01, -1.3921E+01],
                   [-1.3553E+01, -1.3921E+01, -1.3824E+01, -1.3538E+01, -1.3553E+01],
                   [-1.3538E+01, -1.3432E+01, -1.3553E+01, -1.3538E+01, -1.3432E+01]])

    interp1 = interp2d(tint, tnnt, cx)

    Ti_exps = np.linspace(-1, 3, 100)
    Tn_exps = np.linspace(0, 2, 100)
    svcx_vals = 10.0 ** (interp1(Ti_exps, Tn_exps))  # in m^3/s

    Ti_vals = np.logspace(-1, 3, 100) * e  # in joules
    Tn_vals = np.logspace(0, 2, 100) * e  # in joules

    dsvcx_dTi_vals = np.gradient(svcx_vals, Ti_vals, axis=0)

    Ti_vals2d, Tn_vals2d = np.meshgrid(Ti_vals, Tn_vals)

    Ti_mod = np.where(T.i.ev > 1E3, 1E3 * e, T.i.ev * e)
    Tn_mod = np.zeros(Ti_mod.shape) + 2.0 * e

    sv_cx = griddata(np.column_stack((Ti_vals2d.flatten(), Tn_vals2d.flatten())),
                          svcx_vals.flatten(),
                          (Ti_mod, Tn_mod),
                          method='linear', rescale=False)

    dsv_cx_dT = griddata(np.column_stack((Ti_vals2d.flatten(), Tn_vals2d.flatten())),
                              dsvcx_dTi_vals.flatten(),
                              (Ti_mod, Tn_mod),
                              method='linear', rescale=False)

    return sv_cx, dsv_cx_dT


def calc_svion_st(T):
    # TODO: configure so it can use any of the cross section libraries
    # currently using the Stacey-Thomas cross sections
    T_exps_fit = np.array([-1, 0, 1, 2, 3, 4, 5])
    sigv_exps_fit = np.array([-2.8523E+01, -1.7745E+01, -1.3620E+01,
                              -1.3097E+01, -1.3301E+01, -1.3301E+01, -1.3301E+01])
    interp1 = UnivariateSpline(T_exps_fit, sigv_exps_fit, s=0)

    T_exps_range = np.linspace(-1, 5, 1000)
    sigv_vals_range = 10.0 ** interp1(T_exps_range)  # in m^3/s

    T_vals_range = np.logspace(-1, 5, 1000) * e  # in joules
    interp2 = UnivariateSpline(T_vals_range, sigv_vals_range, s=0)

    sv_ion = interp2(T.i.J)
    dsv_ion_dT = interp2.derivative()(T.i.J)

    return sv_ion, dsv_ion_dT


def calc_svel_st(T):
    tint = np.array([-1, 0, 1, 2, 3])
    tnnt = np.array([0, 1, 2])

    elast = np.array([[-1.3569E+01, -1.3337E+01, -1.3036E+01, -1.3569E+01, -1.3337E+01],
                      [-1.3036E+01, -1.3337E+01, -1.3167E+01, -1.3046E+01, -1.3036E+01],
                      [-1.3046E+01, -1.2796E+01, -1.3036E+01, -1.3046E+01, -1.2796E+01]])

    interp1 = interp2d(tint, tnnt, elast)

    Ti_exps = np.linspace(-1, 3, 100)
    Tn_exps = np.linspace(0, 2, 100)
    svel_vals = 10.0 ** (interp1(Ti_exps, Tn_exps))  # in m^3/s

    Ti_vals = np.logspace(-1, 3, 100) * e  # in joules
    Tn_vals = np.logspace(0, 2, 100) * e  # in joules

    dsvel_dTi_vals = np.gradient(svel_vals, Ti_vals, axis=0)

    Ti_vals2d, Tn_vals2d = np.meshgrid(Ti_vals, Tn_vals)

    Ti_mod = np.where(T.i.ev > 1E3, 1E3 * e, T.i.ev * e)
    Tn_mod = np.zeros(Ti_mod.shape) + 2.0 * e

    sv_el = griddata(np.column_stack((Ti_vals2d.flatten(), Tn_vals2d.flatten())),
                          svel_vals.flatten(),
                          (Ti_mod, Tn_mod),
                          method='linear', rescale=False)
    dsv_el_dT = griddata(np.column_stack((Ti_vals2d.flatten(), Tn_vals2d.flatten())),
                              dsvel_dTi_vals.flatten(),
                              (Ti_mod, Tn_mod),
                              method='linear', rescale=False)

    return sv_el, dsv_el_dT


########################################

# isclose is included in python3.5+, so you can delete this if the code ever gets ported into python3.5+
def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def calc_fsa(x, R, Z):
    R1 = R[:, :-1]
    R2 = np.roll(R[:, :-1], -1, axis=1)
    Z1 = Z[:, :-1]
    Z2 = np.roll(Z[:, :-1], -1, axis=1)
    x1 = x[:, :-1]
    x2 = np.roll(x[:, :-1], -1, axis=1)

    dl = np.sqrt((R2 - R1) ** 2 + (Z2 - Z1) ** 2)

    R_av = (R1 + R2)/2

    dA = dl * (2 * pi * R_av)

    x_av = (x1 + x2)/2

    fsa = np.sum(x_av * dA, axis=1) / np.sum(dA, axis=1)
    fsa[0] = x[0,0]
    return fsa


def draw_contour_line(R, Z, array, val, pathnum):
    c = QuadContourGenerator.from_rectilinear(R[0], Z[:, 0], array)

    res = c.contour(val)[pathnum]
    x = res[:, 0]
    y = res[:, 1]
    return x, y


def cut(line, distance):
    # Cuts a line in two at a distance from its starting point
    if distance <= 0.0 or distance >= 1.0:
        return [LineString(line)]
    coords = list(line.coords)
    for i, p in enumerate(coords):
        pd = line.project(Point(p), normalized=True)
        if pd == distance:
            return [
                LineString(coords[:i+1]), 
                LineString(coords[i:])]
        if pd > distance:
            cp = line.interpolate(distance, normalized=True)
            return [
                LineString(coords[:i] + [(cp.x, cp.y)]), 
                LineString([(cp.x, cp.y)] + coords[i:])]


def draw_core_line(R, Z, psi, psi_val, sep_pts):
    # create contour generator
    c = QuadContourGenerator.from_rectilinear(R[0], Z[:, 0], psi)

    # draw contours with psi_val
    contours = c.contour(psi_val)

    if len(contours) == 0 and isclose(psi_val, 0, abs_tol=1E-9):
        # This is probably the magnetic axis
        line = None
        fs_axis = None
    elif len(contours) == 0 and not isclose(psi_val, 0, abs_tol=1E-9):
        # This either means that either:
        #    A) psi is a value not present in the psi data or
        #    B) psi is very close to the magnetic axis, but is too small for the contours
        #       package to give us a contour. This can happen if you you have an very fine
        #       radial mesh in the vicinity of the magnetic axis.
        # Passing back None for now, but should probably raise an exception in the future  # TODO
        line = None
        fs_axis = None
    else:
        if len(contours) == 1:
            # then we're definitely dealing with a surface inside the seperatrix
            x, y = draw_contour_line(R, Z, psi, psi_val, 0)
        else:
            # we need to find which of the surfaces is inside the seperatrix
            for j, line in enumerate(contours):
                x, y = draw_contour_line(R, Z, psi, psi_val, j)

                if (np.amax(x) < np.amax(sep_pts[:, 0]) and
                    np.amin(x) > np.amin(sep_pts[:, 0]) and
                    np.amax(y) < np.amax(sep_pts[:, 1]) and
                    np.amin(y) > np.amin(sep_pts[:, 1])):
                    # then it's an internal flux surface
                    break

        pts = np.column_stack((x, y))
        line = LineString(pts)
        out_pt = pts[np.argmax(pts, axis=0)[0]]
        in_pt = pts[np.argmin(pts, axis=0)[0]]
        top_pt = pts[np.argmax(pts, axis=0)[1]]
        bot_pt = pts[np.argmin(pts, axis=0)[1]]
        fs_axis = [(out_pt[0]+in_pt[0])/2, (out_pt[1]+in_pt[1])/2]
    return line, fs_axis


def calc_rho1d(edge_rho=None, rhopts_core=None, rhopts_edge=None, rhopts=None):
    # define rho points
    try:
        rho1d = np.concatenate((np.linspace(0, edge_rho, rhopts_core, endpoint=False),
                                np.linspace(edge_rho, 1, rhopts_edge, endpoint=True)))
    except:
        try:
            rho1d = np.linspace(0, 1, rhopts)
        except:
            print 'rho parameters not defined. Using 100 evenly spaced rho values'
            rho1d = np.linspace(0, 1, 100)
    return rho1d


def find_xpt_mag_axis(R, Z, psi):
    # calculate gradients of psi in R and Z directions
    dpsidR = np.gradient(psi, R[0, :], axis=1)
    dpsidZ = np.gradient(psi, Z[:, 0], axis=0)

    # find line(s) where dpsidR=0
    c_dpsidR = QuadContourGenerator.from_rectilinear(R[0], Z[:,0], dpsidR)
    c_dpsidZ = QuadContourGenerator.from_rectilinear(R[0], Z[:,0], dpsidZ)
    dpsidR_0 = c_dpsidR.contour(0.0)
    dpsidZ_0 = c_dpsidZ.contour(0.0)

    # populate list all intersection points (as Shapely points)
    int_pts = []
    for i, path1 in enumerate(dpsidR_0):
        for j, path2 in enumerate(dpsidZ_0):
            try:
                # find intersection points between curves for dpsidR=0 and dpsidZ=0
                # these correspond to local minima, local maxima, or saddle-points in psi
                ints = LineString(path1).intersection(LineString(path2))
                if ints.type == 'Point':
                    int_pts.append(ints)
                elif ints.type == 'MultiPoint':
                    for pt in ints:
                        int_pts.append(pt)
            except:
                pass

    # magnetic axis is the intersection point with the minimum value of psi
    psi_int_pts = griddata(np.column_stack((R.flatten(), Z.flatten())),
                           psi.flatten(),
                           np.asarray(MultiPoint(int_pts)),
                           method='cubic')
    mag_axis = np.asarray(MultiPoint(int_pts))[np.argmin(psi_int_pts)]

    # find the dpsidR_0 = 0 contour that contains the magnetic axis.
    # Of the other intersection points along the same dpsidR_0 = 0 contour, the x-point
    # will be the highest (in Z) of those points that are below the magnetic axis.

    # create list of potential x-points (i.e. points that also fall on the same dpsidR_0 = 0 contour
    # as the magnetic axis and whose y-values are lower than that of the magnetic axis.
    pot_xpts = []
    for i, path in enumerate(dpsidR_0):
        line = LineString(path)
        if line.distance(Point(mag_axis)) < 1E-8:
            # then we've found the contour that includes the magnetic axis
            for pt in int_pts:
                if line.distance(pt) < 1E-8 and pt.y < Point(mag_axis).y:
                    pot_xpts.append(pt)

    # of the potential x-points, take the one with the largest y value (there might be only one potential x-point)
    pts_xpts_array = np.asarray(MultiPoint(pot_xpts))
    xpt = pts_xpts_array[pts_xpts_array[:, 1].argsort()][-1]

    return xpt, mag_axis


def calc_psi_norm(R, Z, psi, xpt, axis_mag):
    # normalize psi
    # psi_interp = Rbf(R, Z, psi)
    # psi_min = psi_interp(axis_mag[0], axis_mag[1])
    #
    # psi_shifted = psi - psi_min  # set center to zero
    # psi_shifted_interp = Rbf(R, Z, psi_shifted)
    # psi_shifted_xpt = psi_shifted_interp(xpt[0], xpt[1])
    #
    # psi_norm = psi_shifted / psi_shifted_xpt

    psi_min = griddata(np.column_stack((R.flatten(), Z.flatten())),
                       psi.flatten(),
                       [axis_mag[0], axis_mag[1]],
                       method='cubic')

    psi_shifted = psi - psi_min  # set center to zero

    psi_shifted_xpt = griddata(np.column_stack((R.flatten(), Z.flatten())),
                               psi_shifted.flatten(),
                               [xpt[0], xpt[1]],
                               method='cubic')

    psi_norm = psi_shifted / psi_shifted_xpt

    return psi_norm


def calc_pts_lines(psi_data, xpt, wall, mag_axis):
    # create lines for seperatrix and divertor legs of seperatrix

    c = QuadContourGenerator.from_rectilinear(psi_data.R[0], psi_data.Z[:, 0], psi_data.psi_norm)
    contours = c.contour(1.0)

    # Replace points in the vic. of the xpt with the xpt
    contour_new = []
    for contour in contours:
        dist = np.sqrt((contour[:, 0] - xpt[0])**2 + (contour[:, 1] - xpt[1])**2)
        contour[dist < 0.05] = xpt
        contour_new.append(contour[np.where((contour != np.roll(contour, -1, axis=0)).all(axis=1))])

    # if there is only one contour, then it's drawing one continuous line from ib strike point to ob strike point
    # if there are two contours, then it's drawing separate seperatrix and divertor legs
    if len(contour_new) == 1:
        contour = contour_new[0]

        # make sure the line goes from inboard to outboard
        if contour[-1, 0] < contour[0, 0]:
            contour = np.flipud(contour_new)

        xpt_loc = np.where((contour == xpt).all(axis=1))[0]

        ib = np.flipud(contour[:xpt_loc[0]+1])
        lcfs = contour[xpt_loc[0]:xpt_loc[1]]
        ob = contour[xpt_loc[1]:]

    elif len(contour_new) == 2:
        for contour in contour_new:
            # determine if contour is seperatrix or divertor legs by checking if
            # the largest y value is larger than the y value of the magnetic axis
            if np.amax(contour[:, 1]) > mag_axis[1]:  # we have the main seperatrix
                lcfs = contour

            else:  # we have the divertor legs
                # make sure the line goes from inboard to outboard
                if contour[-1, 0] < contour[0, 0]:
                    contour = np.flipud(contour)

                # create ib and ob divertor legs, each starting at the x-point
                xpt_loc = np.where((contour == xpt).all(axis=1))[0][0]
                ob = contour[xpt_loc:]
                ib = np.flipud(contour[:xpt_loc+1])

    # create lcfs lines
    lcfs_line = LineString(lcfs)
    lcfs_line_closed = LinearRing(lcfs)

    # create ib and ob linestrings and truncate at the wall
    ib_line = LineString(ib)
    ob_line = LineString(ob)

    ib_strike = ib_line.intersection(wall)
    ob_strike = ob_line.intersection(wall)

    ib_line_cut = cut(ib_line, ib_line.project(ib_strike, normalized=True))[0]
    ob_line_cut = cut(ob_line, ib_line.project(ob_strike, normalized=True))[0]

    entire_sep = np.vstack((np.flipud(ib), lcfs[1:], ob))
    entire_sep_line = LineString(entire_sep)

    # create points object
    obmp_pt = lcfs[np.argmax(lcfs, axis=0)[0]]
    ibmp_pt = lcfs[np.argmin(lcfs, axis=0)[0]]
    top_pt = lcfs[np.argmax(lcfs, axis=0)[1]]
    bot_pt = lcfs[np.argmin(lcfs, axis=0)[1]]
    axis_geo = [(obmp_pt[0] + ibmp_pt[0]) / 2, (obmp_pt[1] + ibmp_pt[1]) / 2]

    axis_nt = namedtuple('axis', 'geo mag')
    axis = axis_nt(axis_geo, mag_axis)

    strike_nt = namedtuple('strike', 'ib ob')
    strike = strike_nt(ib_strike, ob_strike)

    pts_nt = namedtuple('pts', 'ibmp obmp top bottom xpt axis strike')
    pts = pts_nt(ibmp_pt, obmp_pt, top_pt, bot_pt, xpt, axis, strike)

    # create lines object
    div_lines_nt = namedtuple('div_lines', 'ib ob')
    div_lines = div_lines_nt(ib_line_cut, ob_line_cut)

    lines_nt = namedtuple('lines', 'sep sep_closed div ib2ob')
    lines = lines_nt(lcfs_line, lcfs_line_closed, div_lines, entire_sep_line)

    # plt.contourf(psi_data.R, psi_data.Z, psi_data.psi_norm, 500)
    # plt.colorbar()
    # plt.plot(ib[:,0],ib[:,1],color='blue')
    # plt.plot(ob[:,0],ob[:,1],color='red')
    # plt.plot(lcfs[:,0],lcfs[:,1],color='yellow')
    # plt.scatter(xpt[0],xpt[1],color='black')
    # plt.show()

    return pts, lines


def calc_theta1d(pts, thetapts_approx):
    # define theta points
    def atan3(y, x):
        result = atan2(y, x)
        if result < 0:
            result = result + 2 * pi
        return result

    # these theta markers correspond to the obmp, top, ibmp, bot, and obmp+2pi, respectively
    theta_markers = np.zeros(5)
    theta_markers[0] = atan2((pts.obmp[1] - pts.axis.geo[1]), (pts.obmp[0] - pts.axis.geo[0]))
    theta_markers[1] = atan3((pts.top[1] - pts.axis.geo[1]), (pts.top[0] - pts.axis.geo[0]))
    theta_markers[2] = atan3((pts.ibmp[1] - pts.axis.geo[1]), (pts.ibmp[0] - pts.axis.geo[0]))
    theta_markers[3] = atan3((pts.xpt[1] - pts.axis.geo[1]), (pts.xpt[0] - pts.axis.geo[0]))
    theta_markers[4] = theta_markers[0] + 2 * pi

    try:
        min_delta_theta = 2 * pi / thetapts_approx
    except:
        print 'thetapts_approx not defined. Setting to 30'
        min_delta_theta = 2 * pi / 30

    theta1d = np.zeros(0)
    for i in range(4):
        quad_pts = int(ceil((theta_markers[i + 1] - theta_markers[i]) / min_delta_theta))
        if i == 3:
            quad_theta = np.linspace(theta_markers[i], theta_markers[i + 1], quad_pts+1, endpoint=True)
        else:
            quad_theta = np.linspace(theta_markers[i], theta_markers[i + 1], quad_pts, endpoint=False)
        theta1d = np.concatenate((theta1d, quad_theta))
    return theta1d, theta_markers


def calc_psi(rho, pts, psi_data):
    # move along line between m_axis and obmp and define psi values corresponding to rho values
    rho_line = LineString([Point(pts.axis.mag), Point(pts.obmp)])
    init_coords = np.zeros((0, 2))
    for i, rhoval in enumerate(rho[:,0]):
        pt_coords = np.asarray(rho_line.interpolate(rhoval, normalized=True).coords)[0]
        init_coords = np.vstack((init_coords, pt_coords))

    psi_vals = griddata(np.column_stack((psi_data.R.flatten(), psi_data.Z.flatten())),
                             psi_data.psi.flatten(),
                             init_coords,
                             method='cubic')
    psi_norm_vals = griddata(np.column_stack((psi_data.R.flatten(), psi_data.Z.flatten())),
                             psi_data.psi_norm.flatten(),
                             init_coords,
                             method='cubic')

    psi = interp1d(rho[:,0], psi_vals)(rho)
    psi_norm = interp1d(rho[:, 0], psi_norm_vals)(rho)

    return psi, psi_norm


def calc_RZ(rho, theta, theta_xpt, pts, psi_data, psi_norm, lines):
    # get parameters that depend on both rho and theta

    sep_pts = np.asarray(lines.sep_closed.coords)

    R = np.zeros(rho.shape)
    Z = np.zeros(rho.shape)
    for i, psi_norm_val in enumerate(psi_norm[:, 0]):
        if i == 0:
            R[i] = pts.axis.mag[0]
            Z[i] = pts.axis.mag[1]
        else:
            # attempt to draw flux surface line through that point
            # (may not work for flux surfaces very close to the magnetic axis)
            fs_line, fs_axis = draw_core_line(psi_data.R, psi_data.Z, psi_data.psi_norm, psi_norm_val, sep_pts)
            if fs_line == None and fs_axis == None:
                # then draw_core_line didn't find any contours, which probably means it was trying
                # to draw a surface closer to psi_norm=0 than the contours package would cooperate with.
                # When this happens, the best thing to do for now is decrease the radial resolution in
                # the vicinity of the magnetic axis. Ideas for the fixing this in the future:
                #   1) resample the raw psi data (maybe an Rbf interpolation) onto a finer mesh. May or may not work.
                #   2) some kind of interpolation based on the flux surfaces it can draw.  # TODO
                print '\nGT3 had trouble drawing a contour line when getting the R and Z points. This ' \
                      'is most likely due to an overly fine radial mesh in the vicnity of the magnetic ' \
                      'axis. Try reducing your number of radial meshes in the core and try again. This ' \
                      'will hopefully be fixed in a future update. Stopping.'
                sys.exit()

            for j, thetaval in enumerate(theta[0]):
                if psi_norm_val < 1.0:
                    thetaline = LineString([Point(fs_axis),
                                            Point([3.0 * cos(thetaval) + fs_axis[0],
                                                   3.0 * sin(thetaval) + fs_axis[1]])])
                    int_pt = fs_line.intersection(thetaline)
                else:
                    if thetaval == theta_xpt:
                        int_pt = Point(pts.xpt)

                    else:
                        thetaline = LineString([Point(pts.axis.geo),
                                                Point([3.0 * cos(thetaval) + pts.axis.geo[0],
                                                       3.0 * sin(thetaval) + pts.axis.geo[1]])])
                        int_pt = LineString(sep_pts).intersection(thetaline)

                R[i, j] = int_pt.x
                Z[i, j] = int_pt.y

    return R, Z


def calc_nT_grad(rho, quant, psi, R, Z, psi_data):
    rhovals = rho[:,0]
    psivals = psi[:,0]
    vals = quant[:,0]

    # calculate values as function of psi and get psi derivative function
    psi_fit = UnivariateSpline(psivals, vals, k=3, s=2.0)
    d_dpsi_fit = psi_fit.derivative()

    # get value of dval_dpsi on the main computational grid
    dval_dpsi = d_dpsi_fit(psi)

    # calculate dpsi_norm_dr everywhere on the main computational grid

    #plt.contourf(psi_data.R, psi_data.Z, psi_data.dpsidr, 500)

    #plt.show()
    #sys.exit()

    dpsi_dr = griddata(np.column_stack((psi_data.R.flatten(),
                                        psi_data.Z.flatten())),
                       psi_data.dpsidr.flatten(),
                       (R,Z),
                       method='linear')

    dval_dr = dval_dpsi * dpsi_dr

    #plt.plot(rhovals, vals,label='rho')
    #plt.plot(psivals, vals,label='psi')
    #plt.legend()
    #plt.show()
    #sys.exit()

    return dval_dr


def calc_rho2psi_interp(pts, psi_data):
    rhovals = np.linspace(0, 1, 100)
    ptsRZ = np.zeros((len(rhovals), 2))

    obmp_line = LineString([Point(pts.axis.mag), Point(pts.obmp)])
    for i, rho in enumerate(rhovals):
        ptsRZ[i] = np.asarray(obmp_line.interpolate(rho, normalized=True).coords)

    psivals = griddata(np.column_stack((psi_data.R.flatten(), psi_data.Z.flatten())),
                        psi_data.psi_norm.flatten(),
                        ptsRZ,
                        method='linear')

    rho2psi = interp1d(rhovals, psivals)
    psi2rho = interp1d(psivals, rhovals)

    return rho2psi, psi2rho


class ExpCoreBrnd():
    def __init__(self, inp):

        raw_psi_R = inp.psirz_exp[:, 0].reshape(-1, 65)
        raw_psi_Z = inp.psirz_exp[:, 1].reshape(-1, 65)
        raw_psi = inp.psirz_exp[:, 2].reshape(-1, 65)

        xpt, mag_axis = find_xpt_mag_axis(raw_psi_R, raw_psi_Z, raw_psi)
        raw_psi_norm = calc_psi_norm(raw_psi_R, raw_psi_Z, raw_psi, xpt, mag_axis)

        raw_dpsidR = np.abs(np.gradient(raw_psi_norm, raw_psi_R[0, :], axis=1))
        raw_dpsidZ = np.abs(np.gradient(raw_psi_norm, raw_psi_Z[:, 0], axis=0))
        raw_dpsidr = raw_dpsidR + raw_dpsidZ

        psi_data_nt = namedtuple('psi_data', 'R Z psi psi_norm dpsidR dpsidZ dpsidr')
        self.psi_data = psi_data_nt(
            raw_psi_R,
            raw_psi_Z,
            raw_psi,
            raw_psi_norm,
            raw_dpsidR,
            raw_dpsidZ,
            raw_dpsidr
        )

        # calculate some important points and lines
        self.pts, self.lines = calc_pts_lines(self.psi_data, xpt, inp.wall_line, mag_axis)
        self.R0_a = self.pts.axis.geo[0]

        # TODO: Figure out which of these is the definition of 'a'
        # self.a = self.pts.obmp[0] - self.pts.axis.geo
        self.a = self.pts.obmp[0] - self.pts.axis.mag[0]
        self.vol = (pi*self.a**2) * 2 * pi * self.R0_a  # TODO: This is only a rough estimate. Replace later.

        # calculate rho and theta values and number of thetapts

        # specify rho values
        try:
            rho1d = np.concatenate((np.linspace(0, inp.edge_rho, inp.rhopts_core, endpoint=False),
                                    np.linspace(inp.edge_rho, 1, inp.rhopts_edge)), axis=0)
        except AttributeError:
            try:
                rho1d = np.linspace(0, 1, inp.rhopts)
            except AttributeError:
                raise AttributeError("You haven't specified the number of radial points.")
        self.rhopts = len(rho1d)

        theta1d, theta_markers = calc_theta1d(self.pts, inp.thetapts_approx)
        theta_xpt = theta_markers[3]
        self.thetapts = len(theta1d)

        # create rho and theta arrays (initializing the main computational grid)
        self.theta, self.rho = np.meshgrid(theta1d, rho1d)
        self.r = self.rho * self.a
        self.psi, self.psi_norm = calc_psi(self.rho, self.pts, self.psi_data)

        self.R, self.Z = calc_RZ(self.rho, self.theta, theta_xpt, self.pts, self.psi_data, self.psi_norm, self.lines)

        # create interpolation functions to convert rho to psi and vice versa
        self.rho2psi, self.psi2rho = calc_rho2psi_interp(self.pts, self.psi_data)

        # populate q-profile from input data
        # TODO: in the future, get this from psi data
        self.q = UnivariateSpline(inp.q_data[:, 0], inp.q_data[:, 1], k=5, s=2.0)(self.rho)

        # initialize ionization rate arrays with zero
        self.izn_rate = namedtuple('izn_rate', 's t tot')(
            np.zeros(self.rho.shape),  # slow
            np.zeros(self.rho.shape),  # thermal
            np.zeros(self.rho.shape)   # total
        )

        # initialize cooling rate array with zero
        self.cool_rate = np.zeros(self.rho.shape)

        # initialize main density object
        fracz = UnivariateSpline(inp.fracz_data[:, 0], inp.fracz_data[:, 1], k=5, s=2.0)(self.rho)
        ni = UnivariateSpline(inp.ni_data[:, 0], inp.ni_data[:, 1], k=3, s=2.0)(self.rho)
        ne = UnivariateSpline(inp.ne_data[:, 0], inp.ne_data[:, 1], k=3, s=2.0)(self.rho)
        nC = ne * fracz
        nn = namedtuple('nn', 's t tot')(
            np.zeros(self.rho.shape),  # slow
            np.zeros(self.rho.shape),  # thermal
            np.zeros(self.rho.shape)   # total
        )
        self.n = namedtuple('n', 'i e n C')(ni, ne, nn, nC)
        self.z_0 = nC*6.0**2 / ni
        self.z_eff = (ni*1.0**2 + nC*6.0**2) / ne

        # populate temperature namedtuples, including the main T object
        Ti_kev = UnivariateSpline(inp.Ti_data[:, 0], inp.Ti_data[:, 1], k=5, s=2.0)(self.rho)
        Te_kev = UnivariateSpline(inp.Te_data[:, 0], inp.Te_data[:, 1], k=5, s=2.0)(self.rho)
        Tns_kev = np.full(self.rho.shape, 0.002)
        Tnt_kev = Ti_kev
        TC_kev = Ti_kev

        self.T = namedtuple('T', 'i e n C')(
            namedtuple('Ti', 'kev ev J')(
                Ti_kev,
                Ti_kev * 1E3,
                Ti_kev * 1E3 * e),
            namedtuple('eT', 'kev ev J')(
                Te_kev,
                Te_kev * 1E3,
                Te_kev * 1E3 * e),
            namedtuple('Tn', 's t')(
                namedtuple('Tn_s', 'kev ev J')(
                    Tns_kev,
                    Tns_kev * 1E3,
                    Tns_kev * 1E3 * e),
                namedtuple('Tn_t', 'kev ev J')(
                    Tnt_kev,
                    Tnt_kev * 1E3,
                    Tnt_kev * 1E3 * e)
            ),
            namedtuple('TC', 'kev ev J')(
                TC_kev,
                TC_kev * 1E3,
                TC_kev * 1E3 * e))

        # initialize spatial gradients and gradient scale lengths
        dni_dr = calc_nT_grad(self.rho, self.n.i, self.psi_norm, self.R, self.Z, self.psi_data)
        dne_dr = calc_nT_grad(self.rho, self.n.e, self.psi_norm, self.R, self.Z, self.psi_data)
        dnC_dr = calc_nT_grad(self.rho, self.n.C, self.psi_norm, self.R, self.Z, self.psi_data)
        dTi_kev_dr = calc_nT_grad(self.rho, self.T.i.kev, self.psi_norm, self.R, self.Z, self.psi_data)
        dTi_ev_dr = dTi_kev_dr * 1E3
        dTi_J_dr = dTi_kev_dr * 1E3 * e
        dTe_kev_dr = calc_nT_grad(self.rho, self.T.i.kev, self.psi_norm, self.R, self.Z, self.psi_data)
        dTe_ev_dr = dTe_kev_dr * 1E3
        dTe_J_dr = dTe_kev_dr * 1E3 * e

        L_ni = -dni_dr / self.n.i
        L_ne = -dne_dr / self.n.e
        L_Ti = -dTi_kev_dr / self.T.i.kev  # note: independent of units
        L_Te = -dTe_kev_dr / self.T.e.kev  # note: independent of units

        Ln_nt = namedtuple('Ln', 'i e')
        LT_nt = namedtuple('LT', 'i e')
        L_nt = namedtuple('L', 'n T')

        self.L = L_nt(Ln_nt(L_ni, L_ne),
                 LT_nt(L_Ti, L_Te))

        # initialize E_r and the corresponding electric potential
        E_r_fit = UnivariateSpline(inp.er_data[:, 0], inp.er_data[:, 1], k=5, s=2.0)
        self.E_r = E_r_fit(self.rho)
        self.E_pot = np.zeros(self.rho.shape)
        for i, rhoval in enumerate(rho1d):
            self.E_pot[i] = E_r_fit.integral(rhoval, 1.0)

        # initialize rotation velocities from data
        try:
            vpolD = UnivariateSpline(inp.vpolD_data[:, 0], inp.vpolD_data[:, 1], k=5, s=2.0)(self.rho)
        except AttributeError:
            vpolD = np.zeros(self.rho.shape)
        try:
            vpolC = UnivariateSpline(inp.vpolC_data[:, 0], inp.vpolC_data[:, 1], k=5, s=2.0)(self.rho)
        except AttributeError:
            vpolC = np.zeros(self.rho.shape)
        try:
            vtorD = UnivariateSpline(inp.vtorD_data[:, 0], inp.vtorD_data[:, 1], k=5, s=2.0)(self.rho)
        except AttributeError:
            vtorD = np.zeros(self.rho.shape)
        try:
            vtorC = UnivariateSpline(inp.vtorC_data[:, 0], inp.vtorC_data[:, 1], k=5, s=2.0)(self.rho)
        except AttributeError:
            vtorC = np.zeros(self.rho.shape)

        self.v = namedtuple('v_comp', 'pol tor')(
            namedtuple('v_spec', 'D C')(vpolD, vpolC),
            namedtuple('v_spec', 'D C')(vtorD, vtorC),
        )

        self.v_1D = namedtuple('v_comp', 'pol tor')(
            namedtuple('v_spec', 'D C')(vpolD[:, 0], vpolC[:, 0]),
            namedtuple('v_spec', 'D C')(vtorD[:, 0], vtorC[:, 0]),
        )

        # initialize magnetic field-related quantities
        B_pol_raw = np.sqrt((np.gradient(self.psi_data.psi, axis=1) / self.psi_data.R) ** 2 +
                            (-np.gradient(self.psi_data.psi, axis=0) / self.psi_data.R) ** 2)

        self.B_p = griddata(np.column_stack((raw_psi_R.flatten(), raw_psi_Z.flatten())),
                            B_pol_raw.flatten(),
                            (self.R, self.Z),
                            method='linear')

        self.B_t = inp.BT0 * self.pts.axis.mag[0] / self.R
        self.B_tot = np.sqrt(self.B_p**2 + self.B_t**2)
        self.f_phi = self.B_t/self.B_tot

        # create Lz-related variables. These will remain zero unless set by the ImpRad module
        Lz_nt = namedtuple('Lz', 's t ddT')
        Lz_ddT_nt = namedtuple('Lz_ddT', 's t')
        Lz_ddT = Lz_ddT_nt(np.zeros(self.rho.shape), np.zeros(self.rho.shape))
        self.Lz = Lz_nt(np.zeros(self.rho.shape), np.zeros(self.rho.shape), Lz_ddT)

        # calculate cross sections on the main computational grid
        svfus_dd, svfus_dd_ddT = calc_svfus(self.T, mode='dd')
        svfus_dt, svfus_dt_ddT = calc_svfus(self.T, mode='dt')
        svrec_st, svrec_st_ddT = calc_svrec_st(self.n, self.T)
        svcx_st, svcx_st_ddT = calc_svcx_st(self.T)
        svion_st, svion_st_ddT = calc_svion_st(self.T)
        svel_st, svel_st_ddT = calc_svel_st(self.T)

        sv_nt = namedtuple('sv', 'fus rec cx ion el')

        sv_fus_nt = namedtuple('sv_fus', 'dd dt d_dT')
        sv_fus_d_dT_nt = namedtuple('sv_fus_d_dT', 'dd dt')

        sv_rec_nt = namedtuple('sv_rec', 'st d_dT')
        sv_rec_d_dT_nt = namedtuple('sv_rec_d_dT', 'st')

        sv_cx_nt = namedtuple('sv_cx', 'st d_dT')
        sv_cx_d_dT_nt = namedtuple('sv_cx_d_dT', 'st')

        sv_ion_nt = namedtuple('sv_ion', 'st d_dT')
        sv_ion_d_dT_nt = namedtuple('sv_ion_d_dT', 'st')

        sv_el_nt = namedtuple('sv_el', 'st d_dT')
        sv_el_d_dT_nt = namedtuple('sv_el_d_dT', 'st')

        sv_fus_d_dT = sv_fus_d_dT_nt(svfus_dd_ddT, svfus_dt_ddT)
        sv_fus = sv_fus_nt(svfus_dd, svfus_dt, sv_fus_d_dT)

        sv_rec_d_dT = sv_rec_d_dT_nt(svrec_st_ddT)
        sv_rec = sv_rec_nt(svrec_st, sv_rec_d_dT)

        sv_cx_d_dT = sv_cx_d_dT_nt(svcx_st_ddT)
        sv_cx = sv_cx_nt(svcx_st, sv_cx_d_dT)

        sv_ion_d_dT = sv_ion_d_dT_nt(svion_st_ddT)
        sv_ion = sv_ion_nt(svion_st, sv_ion_d_dT)

        sv_el_d_dT = sv_el_d_dT_nt(svel_st_ddT)
        sv_el = sv_el_nt(svel_st, sv_el_d_dT)

        self.sv = sv_nt(sv_fus, sv_rec, sv_cx, sv_ion, sv_el)

        # Calculate some 1D Flux-surface averaged quantities
        self.izn_rate_fsa_s = calc_fsa(self.izn_rate.s, self.R, self.Z)
        self.izn_rate_fsa_t = calc_fsa(self.izn_rate.t, self.R, self.Z)
        self.izn_rate_fsa = calc_fsa(self.izn_rate.s + self.izn_rate.t, self.R, self.Z)
        self.cool_rate_fsa = calc_fsa(self.cool_rate, self.R, self.Z)

        self.z_eff_fsa = calc_fsa(self.z_eff, self.R, self.Z)
        self.B_p_fsa = calc_fsa(self.B_p, self.R, self.Z)

        self.n_fsa = namedtuple('n', 'i e n C')(
            calc_fsa(self.n.i, self.R, self.Z),
            calc_fsa(self.n.e, self.R, self.Z),
            namedtuple('nn', 's t tot')(
                calc_fsa(self.n.n.s, self.R, self.Z),  # slow
                calc_fsa(self.n.n.t, self.R, self.Z),  # thermal
                calc_fsa(self.n.n.tot, self.R, self.Z)  # total
            ),
            calc_fsa(self.n.C, self.R, self.Z))

        Ti_kev_fsa = calc_fsa(self.T.i.kev, self.R, self.Z)
        Te_kev_fsa = calc_fsa(self.T.e.kev, self.R, self.Z)
        Tns_kev_fsa = calc_fsa(self.T.n.s.kev, self.R, self.Z)
        Tnt_kev_fsa = calc_fsa(self.T.n.t.kev, self.R, self.Z)
        TC_kev_fsa = calc_fsa(self.T.C.kev, self.R, self.Z)

        self.T_fsa = namedtuple('T', 'i e n C')(
            namedtuple('Ti', 'kev ev J')(
                Ti_kev_fsa,
                Ti_kev_fsa * 1E3,
                Ti_kev_fsa * 1E3 * e),
            namedtuple('eT', 'kev ev J')(
                Te_kev_fsa,
                Te_kev_fsa * 1E3,
                Te_kev_fsa * 1E3 * e),
            namedtuple('Tn', 's t')(
                namedtuple('Tn_s', 'kev ev J')(
                    Tns_kev_fsa,
                    Tns_kev_fsa * 1E3,
                    Tns_kev_fsa * 1E3 * e),
                namedtuple('Tn_t', 'kev ev J')(
                    Tnt_kev_fsa,
                    Tnt_kev_fsa * 1E3,
                    Tnt_kev_fsa * 1E3 * e)
            ),
            namedtuple('TC', 'kev ev J')(
                TC_kev_fsa,
                TC_kev_fsa * 1E3,
                TC_kev_fsa * 1E3 * e))

        # initialize chi_r. This might get updated later
        self.chi_r = np.full(self.rho.shape, 2.0)
    
    def update_ntrl_data(self, data):
        try:
            n_n_s = griddata(np.column_stack((data.R, data.Z)),
                             data.n_n_slow,
                             (self.R, self.Z),
                             method='linear')
        except:
            n_n_s = self.n.n.s

        try:
            n_n_t = griddata(np.column_stack((data.R, data.Z)),
                             data.n_n_thermal,
                             (self.R, self.Z),
                             method='linear')
        except:
            n_n_t = self.n.n.t

        try:
            izn_rate_s = griddata(np.column_stack((data.R, data.Z)),
                                  data.izn_rate_slow,
                                  (self.R, self.Z),
                                  method='linear')
        except:
            izn_rate_s = self.izn_rate.s

        try:
            izn_rate_t = griddata(np.column_stack((data.R, data.Z)),
                                             data.izn_rate_thermal,
                                             (self.R, self.Z),
                                             method='linear')
        except:
            izn_rate_t = self.izn_rate.t

        nn = namedtuple('nn', 's t tot')(
            n_n_s,  # slow
            n_n_t,  # thermal
            n_n_s + n_n_t  # total
        )

        self.n = namedtuple('n', 'i e n C')(self.n.i, self.n.e, nn, self.n.C)

        self.izn_rate = namedtuple('izn_rate', 's t tot')(
            izn_rate_s,  # slow
            izn_rate_t,  # thermal
            izn_rate_s + izn_rate_t  # total
        )

    def update_Lz_data(self, z, Lz, dLzdT):

        self.Lz_slow = Lz(np.log10(self.T.n.s),
                          np.log10(self.n.n.s / self.n.e),
                          np.log10(self.T.e.kev))

        self.dLzdT_slow = dLzdT(np.log10(self.T.n.s),
                          np.log10(self.n.n.s / self.n.e),
                          np.log10(self.T.e.kev))

        self.Lz_thermal = Lz(np.log10(self.T.n.t),
                          np.log10(self.n.n.t / self.n.e),
                          np.log10(self.T.e.kev))

        self.dLzdT_thermal = dLzdT(np.log10(self.T.n.t),
                          np.log10(self.n.n.t / self.n.e),
                          np.log10(self.T.e.kev))



