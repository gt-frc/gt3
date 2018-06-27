#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 18 13:22:31 2018

@author: max
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib._cntr as cntr
from scipy.interpolate import griddata, UnivariateSpline, interp1d, interp2d, Rbf
from scipy.constants import elementary_charge
from shapely.geometry import Point, LineString
from shapely.ops import polygonize, linemerge
import sys
from math import atan2, pi, ceil, sin, cos
from collections import namedtuple


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
    sigv_fus_interp = UnivariateSpline(Ti_range * 1.0E3 * 1.6021E-19, sigv_fus_range, s=0)  # converted to Joules
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
    # Ti_vals = np.logspace(-1, 3, 100) * 1.6021E-19  # in joules
    #
    # dsvrec_dTi_vals = np.gradient(svrec_vals, Ti_vals, axis=0)
    #
    # zni_vals2d, Ti_vals2d = np.meshgrid(zni_vals, Ti_vals)
    #
    # zni_mod = np.where(n.i > 1E22, 1E22, n.i)
    # zni_mod = np.where(n.i < 1E16, 1E16, zni_mod)
    # Ti_mod = np.where(T.i.ev > 1E3, 1E3 * 1.6021E-19, T.i.ev * 1.6021E-19)
    # Ti_mod = np.where(T.i.ev < 1E-1, 1E-1 * 1.6021E-19, Ti_mod)
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

    Ti_vals = np.logspace(-1, 3, 100) * 1.6021E-19  # in joules
    Tn_vals = np.logspace(0, 2, 100) * 1.6021E-19  # in joules

    dsvcx_dTi_vals = np.gradient(svcx_vals, Ti_vals, axis=0)

    Ti_vals2d, Tn_vals2d = np.meshgrid(Ti_vals, Tn_vals)

    Ti_mod = np.where(T.i.ev > 1E3, 1E3 * 1.6021E-19, T.i.ev * 1.6021E-19)
    Tn_mod = np.zeros(Ti_mod.shape) + 2.0 * 1.6021E-19

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

    T_vals_range = np.logspace(-1, 5, 1000) * 1.6021E-19  # in joules
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

    Ti_vals = np.logspace(-1, 3, 100) * 1.6021E-19  # in joules
    Tn_vals = np.logspace(0, 2, 100) * 1.6021E-19  # in joules

    dsvel_dTi_vals = np.gradient(svel_vals, Ti_vals, axis=0)

    Ti_vals2d, Tn_vals2d = np.meshgrid(Ti_vals, Tn_vals)

    Ti_mod = np.where(T.i.ev > 1E3, 1E3 * 1.6021E-19, T.i.ev * 1.6021E-19)
    Tn_mod = np.zeros(Ti_mod.shape) + 2.0 * 1.6021E-19

    sv_el = griddata(np.column_stack((Ti_vals2d.flatten(), Tn_vals2d.flatten())),
                          svel_vals.flatten(),
                          (Ti_mod, Tn_mod),
                          method='linear', rescale=False)
    dsv_el_dT = griddata(np.column_stack((Ti_vals2d.flatten(), Tn_vals2d.flatten())),
                              dsvel_dTi_vals.flatten(),
                              (Ti_mod, Tn_mod),
                              method='linear', rescale=False)

    return sv_el, dsv_el_dT


def draw_contour_line(R, Z, array, val, pathnum):
    res = cntr.Cntr(R, Z, array).trace(val)[pathnum]
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
    num_lines = int(len(cntr.Cntr(R, Z, psi).trace(psi_val))/2)
    if num_lines == 1:
        # then we're definitely dealing with a surface inside the seperatrix
        x, y = draw_contour_line(R, Z, psi, psi_val, 0)
    else:
        # we need to find which of the surfaces is inside the seperatrix
        for j, line in enumerate(cntr.Cntr(R, Z, psi).trace(psi_val)[:num_lines]):
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
    # find x-point location
    dpsidR = np.gradient(psi, R[0, :], axis=1)
    dpsidZ = np.gradient(psi, Z[:, 0], axis=0)
    d2psidR2 = np.gradient(dpsidR, R[0, :], axis=1)
    d2psidZ2 = np.gradient(dpsidZ, Z[:, 0], axis=0)

    # find line(s) where dpsidR=0
    dpsidR_0 = cntr.Cntr(R, Z, dpsidR).trace(0.0)
    # find line(s) where dpsidZ=0
    dpsidZ_0 = cntr.Cntr(R, Z, dpsidZ).trace(0.0)

    for i, path1 in enumerate(dpsidR_0):
        for j, path2 in enumerate(dpsidZ_0):
            try:
                # find intersection points between curves for dpsidR=0 and dpsidZ=0
                ints = LineString(path1).intersection(LineString(path2))
                # if there is only one intersection ('Point'), then we're probably not
                # dealing with irrelevant noise in psi
                if ints.type == 'Point':
                    # check if local maximum or minimum
                    d2psidR2_pt = griddata(np.column_stack((R.flatten(), Z.flatten())),
                                           d2psidR2.flatten(),
                                           [ints.x, ints.y],
                                           method='cubic')
                    d2psidZ2_pt = griddata(np.column_stack((R.flatten(), Z.flatten())),
                                           d2psidZ2.flatten(),
                                           [ints.x, ints.y],
                                           method='cubic')

                    if d2psidR2_pt > 0 and d2psidZ2_pt > 0:
                        # we've found the magnetic axis
                        mag_axis = np.array([ints.x, ints.y])
                    elif d2psidR2_pt < 0 and d2psidZ2_pt < 0:
                        # we've found a magnet. Do nothing.
                        pass
                    elif ints.y < 0:
                        # we've probably found our x-point, although this isn't super robust
                        # and obviously only applies to a single-diverted, lower-null configuration
                        # TODO: make this more robust, I could easily see this failing on some shots
                        xpt = np.array([ints.x, ints.y])

                    # uncomment this line when debugging
                    # print list(ints.coords), d2psidR2(ints.x, ints.y), d2psidZ2(ints.x, ints.y)
            except:
                pass
    return xpt, mag_axis


def calc_psi_norm(R, Z, psi, xpt, axis_mag):
    # normalize psi
    psi_interp = Rbf(R, Z, psi)
    psi_min = psi_interp(axis_mag[0], axis_mag[1])

    psi_shifted = psi - psi_min  # set center to zero
    psi_shifted_interp = Rbf(R, Z, psi_shifted)
    psi_shifted_xpt = psi_shifted_interp(xpt[0], xpt[1])

    psi_norm = psi_shifted / psi_shifted_xpt

    return psi_norm


def calc_pts_lines(psi_data, xpt, wall, mag_axis):
    # create lines for seperatrix and divertor legs of seperatrix

    R = psi_data.R
    Z = psi_data.Z
    psi_norm = psi_data.psi_norm
    num_lines = int(len(cntr.Cntr(R, Z, psi_norm).trace(1.0)) / 2)
    if num_lines == 1:
        # in this case, the contour points that matplotlib returned constitute
        # a single line from inboard divertor to outboard divertor. We need to
        # add in the x-point in at the appropriate locations and split into a
        # main and a lower seperatrix line, each of which will include the x-point.
        x_psi, y_psi = draw_contour_line(R, Z, psi_norm, 1.0, 0)

        loc1 = np.argmax(y_psi > xpt[1])
        loc2 = len(y_psi) - np.argmin(y_psi[::-1] < xpt[1])

        x_psi = np.insert(x_psi, (loc1, loc2), xpt[0])
        y_psi = np.insert(y_psi, (loc1, loc2), xpt[1])

        psi_1_pts = np.column_stack((x_psi, y_psi))
        main_sep_pts = psi_1_pts[loc1:loc2 + 1, :]
        main_sep_line = LineString(main_sep_pts[:-1])
        main_sep_line_closed = LineString(main_sep_pts)

        # get the inboard and outboard divertor legs seperately. This is so that
        # everything that includes the x-point can start with the x-point, which
        # elliminates the risk of tiny triangles in the vicinity of the x-point
        inboard_div_sep = np.flipud(psi_1_pts[:loc1 + 1])
        outboard_div_sep = psi_1_pts[loc2 + 1:]

        # cut inboard line at the wall and add intersection point to wall_line
        line = LineString(inboard_div_sep)
        int_pt = line.intersection(wall)
        ib_div_line = line
        ib_div_line_cut = cut(line, line.project(int_pt, normalized=True))[0]

        # cut outboard line at the wall and add intersection point to wall_line
        line = LineString(outboard_div_sep)
        int_pt = line.intersection(wall)
        ob_div_line = line
        ob_div_line_cut = cut(line, line.project(int_pt, normalized=True))[0]

    elif num_lines == 2:
        # in this case, we have a lower seperatrix trace (line 0), and a main
        # seperatrix trace (line 1).

        # first do lower seperatrix line
        x_psi, y_psi = draw_contour_line(R, Z, psi_norm, 1.0, 0)
        loc = np.argmax(x_psi > xpt[0])

        x_psi = np.insert(x_psi, loc, xpt[0])
        y_psi = np.insert(y_psi, loc, xpt[1])
        psi_1_pts = np.column_stack((x_psi, y_psi))

        inboard_div_sep = np.flipud(psi_1_pts[:loc + 1])
        outboard_div_sep = psi_1_pts[loc + 1:]

        # cut inboard line at the wall and add intersection point to wall_line
        line = LineString(inboard_div_sep)
        int_pt = line.intersection(wall)
        ib_div_line = line
        ib_div_line_cut = cut(line, line.project(int_pt, normalized=True))[0]

        # cut inboard line at the wall and add intersection point to wall_line
        line = LineString(outboard_div_sep)
        int_pt = line.intersection(wall)
        ob_div_line = line
        ob_div_line_cut = cut(line, line.project(int_pt, normalized=True))[0]
        # TODO: add point to wall line

        # now to main seperatrix line
        x_psi, y_psi = draw_contour_line(R, Z, psi_norm, 1.0, 1)
        main_sep_pts = np.insert(np.column_stack((x_psi, y_psi)), 0, xpt, axis=0)
        main_sep_line = LineString(main_sep_pts[:-1])
        main_sep_line_closed = LineString(main_sep_pts)

        # entire_sep_pts = np.vstack((ib_div_pts, sep_pts[1:, :], ob_div_pts))
        # entire_sep_line = LineString(entire_sep_pts)

    ib_div_pts = np.flipud(np.asarray(ib_div_line_cut.xy).T)
    sep_pts = np.asarray(main_sep_line.xy).T
    ob_div_pts = np.asarray(ob_div_line_cut.xy).T

    entire_sep_pts = np.vstack((ib_div_pts, sep_pts[1:, :], ob_div_pts))
    entire_sep_line = LineString(entire_sep_pts)

    # create points object
    obmp_pt = main_sep_pts[np.argmax(main_sep_pts, axis=0)[0]]
    ibmp_pt = main_sep_pts[np.argmin(main_sep_pts, axis=0)[0]]
    top_pt = main_sep_pts[np.argmax(main_sep_pts, axis=0)[1]]
    bot_pt = main_sep_pts[np.argmin(main_sep_pts, axis=0)[1]]
    axis_geo = [(obmp_pt[0] + ibmp_pt[0]) / 2, (obmp_pt[1] + ibmp_pt[1]) / 2]

    axis_nt = namedtuple('axis', 'geo mag')
    axis = axis_nt(axis_geo, mag_axis)

    pts_nt = namedtuple('pts', 'ibmp obmp top bottom xpt axis')
    pts = pts_nt(ibmp_pt, obmp_pt, top_pt, bot_pt, xpt, axis)

    # create lines object
    div_lines_nt = namedtuple('div_lines', 'ib ob')
    div_lines = div_lines_nt(ib_div_line_cut, ob_div_line_cut)

    lines_nt = namedtuple('lines', 'sep sep_closed div ib2ob')
    lines = lines_nt(main_sep_line, main_sep_line_closed, div_lines, entire_sep_line)

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
        quad_pts = ceil((theta_markers[i + 1] - theta_markers[i]) / min_delta_theta)
        quad_theta = np.linspace(theta_markers[i], theta_markers[i + 1], quad_pts)
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
                             method='linear')
    psi_norm_vals = griddata(np.column_stack((psi_data.R.flatten(), psi_data.Z.flatten())),
                             psi_data.psi_norm.flatten(),
                             init_coords,
                             method='linear')

    psi = interp1d(rho[:,0], psi_vals)(rho)
    psi_norm = interp1d(rho[:, 0], psi_norm_vals)(rho)

    return psi, psi_norm


def calc_RZ(rho, theta, theta_xpt, pts, psi_data, psi_norm, lines):
    # get parameters that depend on both rho and theta

    sep_pts = np.asarray(lines.sep_closed.coords)

    R = np.zeros(rho.shape)
    Z = np.zeros(rho.shape)
    for i, psi_norm_val in enumerate(psi_norm[:, 0]):
        # draw flux surface line through that point
        fs_line, fs_axis = draw_core_line(psi_data.R, psi_data.Z, psi_data.psi_norm, psi_norm_val, sep_pts)

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

    R[0] = pts.axis.mag[0]
    Z[0] = pts.axis.mag[1]
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


def rho2psi(pts, rhovals, R, Z, psi):

    psi_vals_R = np.zeros(rhovals.shape)
    psi_vals_Z = np.zeros(rhovals.shape)
    psi_vals = np.zeros(rhovals.shape)

    obmp_line = LineString([Point(pts.axis.mag), Point(pts.obmp)])
    for i, rho in enumerate(rhovals):
        pt = np.asarray(obmp_line.interpolate(rho, normalized=True).coords)[0]
        psi_vals_R[i] = pt[0]
        psi_vals_Z[i] = pt[1]

    psivals = griddata(np.column_stack((R.flatten(), Z.flatten())),
                        psi.flatten(),
                        (psi_vals_R, psi_vals_Z),
                        method='linear')
    # self.rho2psi = interp1d(inp.ni_data[:, 0], psi_vals)
    # self.psi2rho = interp1d(psi_vals, inp.ni_data[:, 0])
    return psivals


class exp_core_brnd():
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
        psi_data = psi_data_nt(
            raw_psi_R,
            raw_psi_Z,
            raw_psi,
            raw_psi_norm,
            raw_dpsidR,
            raw_dpsidZ,
            raw_dpsidr
        )

        # calculate some important points and lines
        self.pts, self.lines = calc_pts_lines(psi_data, xpt, inp.wall_line, mag_axis)

        self.R0_a = self.pts.axis.geo[0]

        # TODO: Figure out which of these is the definition of 'a'
        # self.a = self.pts.obmp[0] - self.pts.axis.geo
        self.a = self.pts.obmp[0] - self.pts.axis.mag[0]

        # calculate rho and theta values and number of thetapts
        rho1d = calc_rho1d(rhopts=inp.rhopts)
        theta1d, theta_markers = calc_theta1d(self.pts, inp.thetapts_approx)
        theta_xpt = theta_markers[3]
        thetapts = len(theta1d)

        # create rho and theta arrays (initializing the main computational grid)
        self.theta, self.rho = np.meshgrid(theta1d, rho1d)
        self.r = self.rho * self.a
        self.psi, self.psi_norm = calc_psi(self.rho, self.pts, psi_data)

        print 'calc_RZ'
        self.R, self.Z = calc_RZ(self.rho, self.theta, theta_xpt, self.pts, psi_data, self.psi_norm, self.lines)

        # populate q-profile from input data
        # TODO: in the future, get this from psi data
        self.q = UnivariateSpline(inp.q_data[:, 0], inp.q_data[:, 1], k=5, s=2.0)(self.rho)

        # define some namedtuples
        nn_nt = namedtuple('nn', 's t tot')
        n_nt = namedtuple('n', 'i e n C')
        T_units_nt = namedtuple('T', 'kev ev J')
        Tn_nt = namedtuple('Tn', 's t')
        T_nt = namedtuple('T', 'i e n C')
        v_spec_nt = namedtuple('v_spec', 'D C')
        v_comp_nt = namedtuple('v_comp', 'pol tor')

        # initialize neutral density arrays with zero
        nn = nn_nt(
            np.zeros(self.rho.shape),  # slow
            np.zeros(self.rho.shape),  # thermal
            np.zeros(self.rho.shape)   # total
        )

        # initialize main density object
        fracz = UnivariateSpline(inp.fracz_data[:, 0], inp.fracz_data[:, 1], k=5, s=2.0)(self.rho)
        ni = UnivariateSpline(inp.ni_data[:, 0], inp.ni_data[:, 1], k=3, s=2.0)(self.rho)
        ne = UnivariateSpline(inp.ne_data[:, 0], inp.ne_data[:, 1], k=3, s=2.0)(self.rho)
        nC = ne * fracz
        self.n = n_nt(ni, ne, nn, nC)
        self.z_0 = nC*6.0**2 / ni
        self.z_eff = (ni*1.0**2 + nC*6.0**2) / ne

        # populate temperature namedtuples, including the main T object
        Ti_kev = UnivariateSpline(inp.Ti_data[:, 0], inp.Ti_data[:, 1], k=5, s=2.0)(self.rho)
        Ti_ev = Ti_kev * 1E3
        Ti_J = Ti_kev * 1E3 * 1.6021E-19

        Te_kev = UnivariateSpline(inp.Te_data[:, 0], inp.Te_data[:, 1], k=5, s=2.0)(self.rho)
        Te_ev = Te_kev * 1E3
        Te_J = Te_kev * 1E3 * 1.6021E-19

        Tns_kev = np.full(self.rho.shape, 0.002)
        Tns_ev = Tns_kev * 1E3
        Tns_J = Tns_kev * 1E3 * 1.6021E-19

        Tnt_kev = Ti_kev
        Tnt_ev = Tnt_kev * 1E3
        Tnt_J = Tnt_kev * 1E3 * 1.6021E-19

        TC_kev = Ti_kev
        TC_ev = Ti_ev
        TC_J = Ti_J

        Ti = T_units_nt(Ti_kev, Ti_ev, Ti_J)
        Te = T_units_nt(Te_kev, Te_ev, Te_J)
        Tn = Tn_nt(
            T_units_nt(Tns_kev, Tns_ev, Tns_J),
            T_units_nt(Tnt_kev, Tnt_ev, Tnt_J)
        )
        TC = T_units_nt(TC_kev, TC_ev, TC_J)
        self.T = T_nt(Ti, Te, Tn, TC)

        # initialize spatial gradients and gradient scale lengths
        dni_dr = calc_nT_grad(self.rho, self.n.i, self.psi_norm, self.R, self.Z, psi_data)
        dne_dr = calc_nT_grad(self.rho, self.n.e, self.psi_norm, self.R, self.Z, psi_data)
        dnC_dr = calc_nT_grad(self.rho, self.n.C, self.psi_norm, self.R, self.Z, psi_data)
        dTi_kev_dr = calc_nT_grad(self.rho, self.T.i.kev, self.psi_norm, self.R, self.Z, psi_data)
        dTi_ev_dr = dTi_kev_dr * 1E3
        dTi_J_dr = dTi_kev_dr * 1E3 * 1.6021E-19
        dTe_kev_dr = calc_nT_grad(self.rho, self.T.i.kev, self.psi_norm, self.R, self.Z, psi_data)
        dTe_ev_dr = dTe_kev_dr * 1E3
        dTe_J_dr = dTe_kev_dr * 1E3 * 1.6021E-19

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

        self.v = v_comp_nt(
            v_spec_nt(vpolD, vpolC),
            v_spec_nt(vtorD, vtorC),
        )

        # initialize magnetic field-related quantities
        B_pol_raw = np.sqrt((np.gradient(psi_data.psi, axis=1) / psi_data.R) ** 2 +
                            (-np.gradient(psi_data.psi, axis=0) / psi_data.R) ** 2)

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

        # initialize chi_r. This might get updated later
        self.chi_r = np.full(self.rho.shape, 2.0)
    
    def update_ntrl_data(self, data):
        try:
            self.n_n_slow = griddata(np.column_stack((data.R, data.Z)),
                                     data.n_n_slow,
                                     (self.R, self.Z),
                                     method='linear')
        except:
            pass

        try:
            self.n_n_thermal = griddata(np.column_stack((data.R, data.Z)),
                                        data.n_n_thermal,
                                        (self.R, self.Z),
                                        method='linear')
        except:
            pass

        try:
            self.izn_rate_slow = griddata(np.column_stack((data.R, data.Z)),
                                          data.izn_rate_slow,
                                          (self.R, self.Z),
                                          method='linear')
        except:
            pass

        try:
            self.izn_rate_thermal = griddata(np.column_stack((data.R, data.Z)),
                                             data.izn_rate_thermal,
                                             (self.R, self.Z),
                                             method='linear')
        except:
            pass

        self.n.n.tot = self.n.n.s + self.n_n_t
        self.izn_rate.tot = self.izn_rate.s + self.izn_rate.t

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


    def core_lines_main(self, inp, R, Z, psi):                    
        # define lines to be used for the main computational grid
        self.core_main_lines = []
        psi_pts_main = np.concatenate((np.linspace(0, 0.8, 20, endpoint=False), np.linspace(0.8, 1.0, 20, endpoint=False)))
        for i, v in enumerate(psi_pts_main):
            num_lines = int(len(cntr.Cntr(R, Z, self.psi_norm_raw).trace(v))/2)
            if num_lines==1:
                # then we're definitely dealing with a surface inside the seperatrix
                x, y = draw_contour_line(R, Z, self.psi_norm_raw, v, 0)
                self.core_main_lines.append(LineString(np.column_stack((x[:-1], y[:-1]))))
            else:
                # we need to find which of the surfaces is inside the seperatrix
                for j, line in enumerate(cntr.Cntr(R, Z, self.psi_norm_raw).trace(v)[:num_lines]):
                # for j, line in enumerate(cntr.Cntr(R, Z, self.psi_norm_raw).trace(v)):
                    x, y = draw_contour_line(R, Z, self.psi_norm_raw, v, j)
                    if (np.amax(x) < np.amax(self.main_sep_pts[:, 0]) and
                        np.amin(x) > np.amin(self.main_sep_pts[:, 0]) and
                        np.amax(y) < np.amax(self.main_sep_pts[:, 1]) and
                        np.amin(y) > np.amin(self.main_sep_pts[:, 1])):
                        # then it's an internal flux surface
                        self.core_main_lines.append(LineString(np.column_stack((x[:-1], y[:-1]))))
                        break

    def core_main(self, inp, R_psi, Z_psi, psi, B_pol_raw):







        #plt.axis('equal')
        #for i, (Rvals, Zvals) in enumerate(zip(self.R, self.Z)):
        #    plt.plot(Rvals, Zvals, lw=0.5)

        #calculate gradients and gradient scale lengths
        #get quantities on a fairly fine R, Z grid for the purpose of taking gradients, etc.
        # R_temp, Z_temp = np.meshgrid(np.linspace(0.98*self.ibmp_pt[0], 1.02*self.obmp_pt[0], 500),
        #                             np.linspace(1.02*self.top_pt[1], 1.02*self.bot_pt[1], 500))
        #
        # ni_grid = griddata(np.column_stack((self.R.flatten(), self.Z.flatten())),
        #                    self.ni.flatten(),
        #                    (R_temp, Z_temp),
        #                    method='cubic',
        #                    fill_value = self.ni[-1, 0])
        # ne_grid = griddata(np.column_stack((self.R.flatten(), self.Z.flatten())),
        #                    self.ne.flatten(),
        #                    (R_temp, Z_temp),
        #                    method='cubic',
        #                    fill_value = self.ne[-1, 0])
        # Ti_kev_grid = griddata(np.column_stack((self.R.flatten(), self.Z.flatten())),
        #                    self.Ti_kev.flatten(),
        #                    (R_temp, Z_temp),
        #                    method='cubic',
        #                    fill_value = self.Ti_kev[-1, 0])
        # Te_kev_grid = griddata(np.column_stack((self.R.flatten(), self.Z.flatten())),
        #                    self.Te_kev.flatten(),
        #                    (R_temp, Z_temp),
        #                    method='cubic',
        #                    fill_value = self.Te_kev[-1, 0])
        
        # dnidr_temp = -1.0*(np.abs(np.gradient(ni_grid, Z_temp[:, 0], axis=1)) + np.abs(np.gradient(ni_grid, R_temp[0, :], axis=0)))
        # dnedr_temp = -1.0*(np.abs(np.gradient(ne_grid, Z_temp[:, 0], axis=1)) + np.abs(np.gradient(ne_grid, R_temp[0, :], axis=0)))
        # dTidr_temp = -1.0*(np.abs(np.gradient(Ti_kev_grid, Z_temp[:, 0], axis=1)) + np.abs(np.gradient(Ti_kev_grid, R_temp[0, :], axis=0)))
        # dTedr_temp = -1.0*(np.abs(np.gradient(Te_kev_grid, Z_temp[:, 0], axis=1)) + np.abs(np.gradient(Te_kev_grid, R_temp[0, :], axis=0)))
        #
        # self.dni_dr = griddata(np.column_stack((R_temp.flatten(), Z_temp.flatten())),
        #                       dnidr_temp.flatten(),
        #                       (self.R, self.Z),
        #                       method='cubic')
        # self.dne_dr = griddata(np.column_stack((R_temp.flatten(), Z_temp.flatten())),
        #                       dnedr_temp.flatten(),
        #                       (self.R, self.Z),
        #                       method='cubic')
        # self.dTi_kev_dr = griddata(np.column_stack((R_temp.flatten(), Z_temp.flatten())),
        #                       dTidr_temp.flatten(),
        #                       (self.R, self.Z),
        #                       method='cubic')
        # self.dTe_kev_dr = griddata(np.column_stack((R_temp.flatten(), Z_temp.flatten())),
        #                       dTedr_temp.flatten(),
        #                       (self.R, self.Z),
        #                       method='cubic')


        self.dTi_J_dr = self.dTi_kev_dr * 1.6021E-16
        self.dTi_ev_dr = self.dTi_kev_dr * 1E3
        self.dTe_J_dr = self.dTe_kev_dr * 1.6021E-16
        self.dTe_ev_dr = self.dTe_kev_dr * 1E3
        self.L_ni = -self.dni_dr / self.ni
        self.L_ne = -self.dne_dr / self.ne
        self.L_Ti_J = -self.dTi_J_dr / self.Ti_J
        self.L_Te_J = -self.dTe_J_dr / self.Te_J
        #plt.axis('equal')
        #plt.contourf(self.R, self.Z, self.dnidr, 500)
        #plt.colorbar()
        #sys.exit()
            
        #create neutrals-related variables. These will remain zero unless set by exp_neutpy_prep or read_ntrl_data modules
        # self.n_n_slow = np.zeros(self.rho.shape)
        # self.n_n_thermal = np.zeros(self.rho.shape)
        # self.n_n_total = np.zeros(self.rho.shape)
        #
        # self.izn_rate_slow = np.zeros(self.rho.shape)
        # self.izn_rate_thermal = np.zeros(self.rho.shape)
        # self.izn_rate_total = np.zeros(self.rho.shape)
        #


    def core_nT_ntrl(self, inp, R, Z, psi):
        #CREATE ARRAYS OF POINTS, DENSITIES AND TEMPERATURES FOR THE NEUTRALS CALCULATION
        
        #Master arrays that will contain all the points we'll use to get n, T
        #throughout the plasma chamber via 2-D interpolation
        self.ni_pts = np.zeros((0, 3), dtype='float')
        self.ne_pts = np.zeros((0, 3), dtype='float')
        self.Ti_kev_pts = np.zeros((0, 3), dtype='float')
        self.Te_kev_pts = np.zeros((0, 3), dtype='float')
        
        ##########################################
        #Calculate n, T throughout the core plasma using radial profile input files, uniform on flux surface
        ni = UnivariateSpline(inp.ni_data[:, 0], inp.ni_data[:, 1], k=5, s=2.0)
        ne = UnivariateSpline(inp.ne_data[:, 0], inp.ne_data[:, 1], k=5, s=2.0)
        Ti_kev = UnivariateSpline(inp.Ti_data[:, 0], inp.Ti_data[:, 1], k=5, s=2.0)
        Te_kev = UnivariateSpline(inp.Te_data[:, 0], inp.Te_data[:, 1], k=5, s=2.0)
        
        #get approximate rho values associated with the psi values we're using
        #draw line between magnetic axis and the seperatrix at the outboard midplane
        rho_line = LineString([Point(self.m_axis), Point(self.obmp_pt)])

        rho_pts = np.linspace(0.7, 1, 5, endpoint=False)
        thetapts = np.linspace(0, 1, 100, endpoint=False)

        for i, rho in enumerate(rho_pts): 
            #get n, T information at the point by interpolating the rho-based input file data
            ni_val = ni(rho)
            ne_val = ne(rho)
            Ti_kev_val = Ti_kev(rho)
            Te_kev_val = Te_kev(rho)
            #get R, Z coordinates of each point along the rho_line
            pt_coords = np.asarray(rho_line.interpolate(rho, normalized=True).coords)[0]

            #get psi value at that point
            psi_val = griddata(np.column_stack((R.flatten(), Z.flatten())), 
                                 self.psi_norm_raw.flatten(), 
                                 pt_coords, 
                                 method='linear')
            #map this n, T data to every point on the corresponding flux surface
            num_lines = int(len(cntr.Cntr(R, Z, self.psi_norm_raw).trace(psi_val))/2)

            if num_lines==1:
                #then we're definitely dealing with a surface inside the seperatrix
                x, y = draw_contour_line(R, Z, self.psi_norm_raw, psi_val, 0)
                surf = LineString(np.column_stack((x, y)))
            else:
                #we need to find which of the surfaces is inside the seperatrix
                for j, line in enumerate(cntr.Cntr(R, Z, self.psi_norm_raw).trace(psi_val)[:num_lines]):
                    #for j, line in enumerate(cntr.Cntr(R, Z, self.psi_norm_raw).trace(v)):
                    x, y = draw_contour_line(R, Z, self.psi_norm_raw, psi_val, j)
                    if (np.amax(x) < np.amax(self.main_sep_pts[:, 0]) and \
                        np.amin(x) > np.amin(self.main_sep_pts[:, 0]) and \
                        np.amax(y) < np.amax(self.main_sep_pts[:, 1]) and \
                        np.amin(y) > np.amin(self.main_sep_pts[:, 1])):
                        #then it's an internal flux surface
                        surf = LineString(np.column_stack((x, y)))
                        break
            
            for j, theta_norm in enumerate(thetapts):
                pt = np.asarray(surf.interpolate(theta_norm, normalized=True).coords).T
                self.ni_pts = np.vstack((self.ni_pts, np.append(pt, ni_val)))
                self.ne_pts = np.vstack((self.ne_pts, np.append(pt, ne_val)))
                self.Ti_kev_pts = np.vstack((self.Ti_kev_pts, np.append(pt, Ti_kev_val)))
                self.Te_kev_pts = np.vstack((self.Te_kev_pts, np.append(pt, Te_kev_val)))

        #Do seperatrix separately so we don't accidentally assign the input n, T data to the divertor legs
        self.ni_sep_val = ni(1.0)
        self.ne_sep_val = ne(1.0)
        self.Ti_kev_sep_val = Ti_kev(1.0)
        self.Te_kev_sep_val = Te_kev(1.0)
        self.Ti_J_sep_val = self.Ti_kev_sep_val * 1.0E3 * 1.6021E-19
        self.Te_J_sep_val = self.Te_kev_sep_val * 1.0E3 * 1.6021E-19
        for j, theta_norm in enumerate(thetapts): 
            pt = np.asarray(self.main_sep_line.interpolate(theta_norm, normalized=False).coords, dtype='float').T
            self.ni_pts = np.vstack((self.ni_pts, np.append(pt, self.ni_sep_val)))
            self.ne_pts = np.vstack((self.ne_pts, np.append(pt, self.ne_sep_val)))
            self.Ti_kev_pts = np.vstack((self.Ti_kev_pts, np.append(pt, self.Ti_kev_sep_val)))
            self.Te_kev_pts = np.vstack((self.Te_kev_pts, np.append(pt, self.Te_kev_sep_val)))
