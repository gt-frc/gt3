#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from contours.quad import QuadContourGenerator
from shapely.geometry import Point, LineString, LinearRing
from Core.Functions.Cut import cut
from collections import namedtuple


def calc_pts_lines(psi_data, xpt, wall, mag_axis):

    """
    Calculates points and lines for the flux surfaces

    :param psi_data:
    :type psi_data: psi
    :param xpt:
    :param wall:
    :param mag_axis:
    :return:
    :rtype: object
    """
    # create lines for seperatrix and divertor legs of seperatrix

    c = QuadContourGenerator.from_rectilinear(psi_data.R[0], psi_data.Z[:, 0], psi_data.psi_norm)
    contours = c.contour(1.0)

    # Replace points in the vic. of the xpt with the xpt
    contour_new = []
    for contour in contours:
        dist = np.sqrt((contour[:, 0] - xpt[0])**2 + (contour[:, 1] - xpt[1])**2)
        # TODO: Put this distance in the input file or make distance selection smarter
        contour[dist < 0.1] = xpt
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

    elif len(contour_new) > 1:
        for contour in contour_new:
            # determine if contour is seperatrix, divertor legs, a main ib2ob line, or something else

            # Meaningless noise lines will not pass near the x-point, so check that first to elliminate them
            if LineString(contour).distance(Point(xpt)) < 0.01:
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

                    xpt_loc = np.where((contour == xpt).all(axis=1))[0]

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
    ob_line_cut = cut(ob_line, ob_line.project(ob_strike, normalized=True))[0]

    entire_sep = np.vstack((np.flipud(ib), lcfs[1:], ob))
    entire_sep_line = LineString(entire_sep)

    # create points object
    obmp_pt = lcfs[np.argmax(lcfs, axis=0)[0]]
    ibmp_pt = lcfs[np.argmin(lcfs, axis=0)[0]]
    top_pt = lcfs[np.argmax(lcfs, axis=0)[1]]
    bot_pt = lcfs[np.argmin(lcfs, axis=0)[1]]
    axis_geo = [(obmp_pt[0] + ibmp_pt[0]) / 2, (obmp_pt[1] + ibmp_pt[1]) / 2]

    pts = namedtuple('pts', 'ibmp obmp top bottom xpt axis strike')(
        ibmp_pt,
        obmp_pt,
        top_pt,
        bot_pt,
        xpt,
        namedtuple('axis', 'geo mag')(
            axis_geo,
            mag_axis),
        namedtuple('strike', 'ib ob')(
            ib_strike,
            ob_strike))

    # create lines object
    div_lines = namedtuple('div_lines', 'ib ob ib_long ob_long')(
        ib_line_cut,
        ob_line_cut,
        ib_line,
        ob_line)

    lines = namedtuple('lines', 'sep sep_closed div ib2ob')(
        lcfs_line,
        lcfs_line_closed,
        div_lines,
        entire_sep_line)

    return pts, lines