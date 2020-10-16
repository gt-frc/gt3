#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from shapely.geometry import Point, LineString, LinearRing
from GT3.Core.Functions.Cut import cut
from collections import namedtuple
import matplotlib.pyplot as plt
from matplotlib.path import Path


def calc_pts_lines(psi_data, xpt, wall, mag_axis, sep_val, core=None):

    """
    Calculates points and lines for the flux surfaces

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
        for contour in contours:
            ls = LineString(contour.vertices)
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
        if len(temp_contours) == 1:
            ls = LineString(temp_contours.vertices)
            if ls.convex_hull.contains(Point(xpt[0])):
                # The upper X-point will define our separatrix
                xpt = xpt[0]
            elif ls.convex_hull.contains(Point(xpt[1])):
                # The lower X-point will define our separatrix
                xpt = xpt[1]
            else:
                raise Exception("There are 2 x-points, but neither are on the separatrix as found in CalcPtsLines.")
        else:
            for path in temp_contours:
                ls = LineString(path.vertices)
                if ls.convex_hull.contains(Point(xpt[0])):
                    # The upper X-point will define our separatrix
                    xpt = xpt[0]
                    contours_paths = path
                    break
                elif ls.convex_hull.contains(Point(xpt[1])):
                    # The lower X-point will define our separatrix
                    xpt = xpt[1]
                    contours_paths = path
                    break
    else:
        # Grab the actual x-point
        if xpt[0] is not None:
            xpt_temp = xpt[0]
        else:
            xpt_temp = xpt[1]

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
        contour[dist < 0.1] = xpt_temp
        contour_new.append(contour[np.where((contour != np.roll(contour, -1, axis=0)).all(axis=1))])

    # if there is only one contour, then it's drawing one continuous line from ib strike point to ob strike point
    # if there are two contours, then it's drawing separate seperatrix and divertor legs

    if len(contour_new) == 1:
        contour = contour_new[0]

        # make sure the line goes from inboard to outboard
        if contour[-1, 0] < contour[0, 0]:
            contour = np.flipud(contour_new)[0]

        xpt_loc = np.where((contour == xpt_temp).all(axis=1))[0]
        ib = np.flipud(contour[:xpt_loc[0]+1])
        lcfs = contour[xpt_loc[0]:xpt_loc[1]]
        ob = contour[xpt_loc[1]:]

    elif len(contour_new) > 1:
        for contour in contour_new:
            # determine if contour is seperatrix, divertor legs, a main ib2ob line, or something else
            contour = np.array(contour)
            if not LineString(contour).intersects(wall) and wall.convex_hull.contains(LineString(contour)):
                lcfs = contour
                continue
            # Meaningless noise lines will not pass near the x-point, so check that first to elliminate them
            if LineString(contour).distance(Point(xpt_temp)) < 0.01:
                # if the largest y value is larger than the y value of the magnetic axis AND the line
                # intersects the wall twice, then it's an ib2ob line and there are more than one psi_norm=1
                # contours because the others are noise. Treat this one the same way we would if there were
                # only one contour.

                # count number of intersections with the wall
                wall_ints = len(LineString(contour).intersection(wall))

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
    lcfs_line = LineString(lcfs)
    lcfs_line_closed = LinearRing(lcfs)

    # create ib and ob linestrings and truncate at the wall
    ib_line = LineString(ib)
    ob_line = LineString(ob)

    ib_strike = ib_line.intersection(wall)
    ob_strike = ob_line.intersection(wall)

    ib_line_cut = cut(ib_line, ib_line.project(ib_strike, normalized=True))[0]
    ob_line_cut = cut(ob_line, ob_line.project(ob_strike, normalized=True))[0]

    entire_sep = np.vstack((np.flipud(ib), lcfs[1:], ob))
    entire_sep_line = LineString(entire_sep)

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