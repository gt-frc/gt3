#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from shapely.geometry import Point, LineString, LinearRing, MultiPoint
from shapely.ops import nearest_points, split, snap
from GT3.Core.Functions.Cut import cut
from GT3.utilities.GT3LineString import GT3LineString
from collections import namedtuple
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.path import Path
from GT3.utilities.PlotBase import PlotBase


def calc_pts_lines(psi_data, xpt, wall, mag_axis, sep_val, norm_xpt, debug=False):

    """
    Calculates points and lines for the flux surfaces

    :param debug: Activate debugging function and plots
    :type debug: bool
    :param psi_data:
    :type psi_data: psi
    :param xpt:
    :param wall:
    :type wall: LineString
    :param mag_axis:
    :return:
    :rtype: object
    """
    Points = namedtuple('Points', 'ibmp obmp top bottom xpt axis strike')
    PlasmaAxis = namedtuple('PlasmaAxis', 'geo mag')
    StrikePoints = namedtuple('StrikePoints', 'ib ob')
    DivertorLines = namedtuple('DivertorLines', 'ib ob ib_long ob_long')
    PlasmaLines = namedtuple('PlasmaLines', 'sep sep_closed div ib2ob')

    if sep_val < 1.0:
        # We have manually overridden the separatrix to exclude the X-points and SOL. We're only able to get a few
        # parameters. We will assume only 1 contour will have this value
        contours = plt.contour(psi_data.R[0], psi_data.Z[:, 0], psi_data.psi_norm, [sep_val]).collections[0].get_paths()
        if debug:
            _plot_contours(contours)
        for contour in contours:
            ls = GT3LineString(contour.vertices)
            if wall.convex_hull.contains(ls):
                ibmp = [ls.bounds[0], mag_axis[1]]
                obmp = [ls.bounds[2], mag_axis[1]]
                geo = [ls.centroid.xy[0][0], ls.centroid.xy[1][0]]
                top = [ls.centroid.xy[0][0], ls.bounds[3]]
                bottom = [ls.centroid.xy[0][0], ls.bounds[1]]
                mag = mag_axis
                pts = Points(ibmp, obmp, top, bottom, xpt,
                             PlasmaAxis(geo, mag),
                             StrikePoints(None, None))
                lines = PlasmaLines(ls, ls, None, None)
                return pts, lines
        raise Exception("Could not find an LCFS with overwritten sep_val")
    #We find the potential sepatartrix lines first
    contours_paths = plt.contour(psi_data.R[0], psi_data.Z[:, 0], psi_data.psi_norm, [sep_val]).collections[0].get_paths()

    if xpt[0] is not None and xpt[1] is not None:
        # We have 2 X-points. We'll basically ignore one, and SOL will be disabled in the rest of GT3. Any neutrals
        # calculations will have to be omitted. We also know that this line (possibly by definition) will contain
        # an x-point and have IB/OB strike points. We'll overwrite contours_paths to say that this is the only
        # contour worth caring about in the rest of the calculation.

        temp_contours = plt.contour(psi_data.R[0], psi_data.Z[:, 0], psi_data.psi_norm, [sep_val]).collections[0].get_paths()
        if debug:
            _plot_contours(temp_contours, title="X-point contours")
        if len(temp_contours) == 1:
            ls = GT3LineString(temp_contours.vertices)
            if ls.convex_hull.contains(Point(xpt[0])):
                # The upper X-point will define our separatrix
                xpt_sep = xpt[0]
            elif ls.convex_hull.contains(Point(xpt[1])):
                # The lower X-point will define our separatrix
                xpt_sep = xpt[1]
            else:
                raise Exception("There are 2 x-points, but neither are on the separatrix as found in CalcPtsLines.")
        else:
            for path in temp_contours:
                ls = GT3LineString(path.vertices, wall=wall)
                if debug:
                    _plot_contours(ls)
                if ls.convex_hull.contains(Point(xpt[0])):
                    # The lower X-point will define our separatrix
                    xpt_sep = xpt[0]
                    xpt_temp = xpt_sep
                    if ls.contains(xpt_sep):
                        contours_paths = path
                    else:
                        print("FUck")
                    break
                elif ls.convex_hull.contains(Point(xpt[1])):
                    # The upper X-point will define our separatrix
                    xpt_sep = xpt[1]
                    xpt_temp = xpt_sep
                    contours_paths = path
                    break
    else:
        # Grab the actual x-point
        if debug:
            _plot_contours(contours_paths)
        if xpt[0] is not None:
            xpt_temp = xpt[0]
        else:
            xpt_temp = xpt[1]

    # This line might render code above unnecessary?

    xpt_temp = norm_xpt

    # create lines for seperatrix and divertor legs of seperatrix

    #Convert back to arrays since I don't want to deal with changing everything to shapely since this already works
    if isinstance(contours_paths, Path):
        contours = [contours_paths.vertices]
    else:
        contours = [a.vertices for a in contours_paths]
    # Replace points in the vic. of the xpt with the xpt
    contour_new = []
    for contour in contours:
        dist = np.sqrt((contour[:, 0] - xpt_temp[0])**2 + (contour[:, 1] - xpt_temp[1])**2)
        # TODO: Put this distance in the input file or make distance selection smarter
        contour[dist < 0.03] = xpt_temp
        contour_new.append(contour[np.where((contour != np.roll(contour, -1, axis=0)).all(axis=1))])

    # if there is only one contour, then it's drawing one continuous line from ib strike point to ob strike point
    # if there are two contours, then it's drawing separate seperatrix and divertor legs

    if len(contour_new) == 1:
        contour = contour_new[0]

        # make sure the line goes from inboard to outboard

        if contour[-1, 0] < contour[0, 0]:
            contour = np.flipud(contour)

        xpt_loc = np.where((contour == xpt_temp).all(axis=1))[0]
        _plot_contours(GT3LineString(contour), title="LCFS Contour")
        try:
            ib = np.flipud(contour[:xpt_loc[0]+1])
        except IndexError:
            # The separatrix did not converge to the x-point. We should cut and project 2 times to create the LCFS.
            raise Exception("The separatrix does not include the x-point.")
            ls = GT3LineString(contour)
            pt = Point(xpt_temp)
            nearest = nearest_points(ls, pt)[0]
            cuts = split(ls, nearest)
            temp_ls = GT3LineString(cuts[0].coords[:]+pt.coords[:]+cuts[1].coords[:])
            _plot_contours(temp_ls, title="Temporary LS from cuts")

        if len(xpt_loc) == 1:
            second_point = _find_ob(contour, xpt_sep)
            xpt_loc = np.append(xpt_loc, second_point)
            # raise Exception("Only 1 xpoint index found. Should be 2. Fix me")

        lcfs = contour[xpt_loc[0]:xpt_loc[1]]
        ob = contour[xpt_loc[1]:]

    elif len(contour_new) > 1:
        for contour in contour_new:
            # determine if contour is seperatrix, divertor legs, a main ib2ob line, or something else
            contour = np.array(contour)
            if not GT3LineString(contour).intersects(wall) and wall.convex_hull.contains(GT3LineString(contour)):
                lcfs = contour
                continue
            # Meaningless noise lines will not pass near the x-point, so check that first to elliminate them
            if LineString(contour).distance(Point(xpt_temp)) < 0.01:
                # if the largest y value is larger than the y value of the magnetic axis AND the line
                # intersects the wall twice, then it's an ib2ob line and there are more than one psi_norm=1
                # contours because the others are noise. Treat this one the same way we would if there were
                # only one contour.

                # count number of intersections with the wall
                wall_ints = len(GT3LineString(contour).intersection(wall))

                if wall_ints >= 2 and np.amax(contour[:, 1]) > mag_axis[1]:
                    # then we probably have an ib2ob line. Treat the same way as above

                    # make sure the line goes from inboard to outboard
                    if contour[-1, 0] < contour[0, 0]:
                        contour = np.flipud(contour_new)

                    xpt_loc = np.where((contour == xpt_temp).all(axis=1))[0]

                    ib = np.flipud(contour[:xpt_loc[0] + 1])
                    lcfs = contour[xpt_loc[0]:xpt_loc[1]]
                    ob = contour[xpt_loc[1]:]

                # if the largest y value is larger than the y value of the magnetic axis
                elif np.amax(contour[:, 1]) > mag_axis[1]:
                    # we have the main seperatrix
                    lcfs = contour

                else:  # we have the divertor legs
                    # make sure the line goes from inboard to outboard
                    if contour[-1, 0] < contour[0, 0]:
                        contour = np.flipud(contour)

                    # create ib and ob divertor legs, each starting at the x-point
                    xpt_loc = np.where((contour == xpt_temp).all(axis=1))[0][0]
                    ob = contour[xpt_loc:]
                    ib = np.flipud(contour[:xpt_loc+1])

    # create lcfs lines
    lcfs_line = GT3LineString(lcfs)
    lcfs_line_closed = LinearRing(lcfs)

    # create ib and ob linestrings and truncate at the wall
    ib_line = GT3LineString(ib)
    ob_line = GT3LineString(ob)

    ib_strike = ib_line.intersection(wall)
    ob_strike = ob_line.intersection(wall)

    ib_line_cut = cut(ib_line, ib_line.project(ib_strike, normalized=True))[0]
    ob_line_cut = cut(ob_line, ob_line.project(ob_strike, normalized=True))[0]

    entire_sep = np.vstack((np.flipud(ib), lcfs[1:], ob))
    entire_sep_line = GT3LineString(entire_sep)

    # create points object
    obmp_pt = lcfs[np.argmax(lcfs, axis=0)[0]]
    ibmp_pt = lcfs[np.argmin(lcfs, axis=0)[0]]
    top_pt = lcfs[np.argmax(lcfs, axis=0)[1]]
    bot_pt = lcfs[np.argmin(lcfs, axis=0)[1]]
    axis_geo = [(obmp_pt[0] + ibmp_pt[0]) / 2, (obmp_pt[1] + ibmp_pt[1]) / 2]

    pts = Points(
        ibmp_pt,
        obmp_pt,
        top_pt,
        bot_pt,
        xpt,
        PlasmaAxis(
            axis_geo,
            mag_axis),
        StrikePoints(
            ib_strike,
            ob_strike))

    # create lines object
    div_lines = DivertorLines(
        ib_line_cut,
        ob_line_cut,
        ib_line,
        ob_line)

    lines = PlasmaLines(
        lcfs_line,
        lcfs_line_closed,
        div_lines,
        entire_sep_line)

    return pts, lines

def _plot_contours(contours, title="No Title"):
    colors = ['r', 'g','b','darkgreen']
    fig = plt.figure()
    plt.ion()

    try:
        contours.__iter__
        num = len(contours)
        if len(contours) == 1:
            fig, ax = plt.subplots(ncols=num, sharex=True, sharey=True)
            for n, c in enumerate(contours):
                contourPlotter = PlotBase()
                contourPlotter.set_default_color(colors[n])
                ax.set_title(title)
                contourPlotter._unknown_data_plot_helper(c, ax)
        else :

            fig, ax = plt.subplots(ncols=num, sharex=True, sharey=True)
            for n,c in enumerate(contours):
                contourPlotter = PlotBase()
                contourPlotter.set_default_color(colors[n])
                ax[n].set_title(title)
                contourPlotter._unknown_data_plot_helper(c, ax[n])
    except AttributeError:
        ax = fig.add_subplot(111)
        ax.set_title(title)
        contourPlotter = PlotBase()
        contourPlotter._unknown_data_plot_helper(contours, ax)
    plt.show()

def _find_ob(contour, xpt):
    xpt_Pt = Point(xpt)
    dist = 10
    res_pt = None
    for coord in contour:
        if coord[0] > xpt[0]:
            if Point(coord).distance(xpt_Pt) < dist:
                res_pt = coord
                dist = Point(coord).distance(xpt_Pt)

    if not res_pt.all():
        raise Exception("Failed to find the second splice point for separatrix splicing.")
    else:
        return np.where((contour==res_pt).all(axis=1))[0][0]