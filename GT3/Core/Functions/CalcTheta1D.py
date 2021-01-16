#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from math import pi, ceil, atan2


def calc_theta1d(pts, thetapts_approx):
    # define theta points
    def atan3(y, x):
        result = atan2(y, x)
        if result < 0:
            result = result + 2 * pi
        return result

    if pts.xpt[1] is None or pts.xpt[0] is None:
        # Grab the actual x-point
        if pts.xpt[0] is not None:
            xpt_temp = pts.xpt[0]
        else:
            xpt_temp = pts.xpt[1]
    else:
        # THere are 2 xpoints. Use the bottom one.
        xpt_temp = pts.xpt[0]

    if xpt_temp[1] > 0.:
    # these theta markers correspond to the obmp, xpt, ibmp, bot, and obmp+2pi, respectively
        theta_markers = np.zeros(5)
        theta_markers[0] = atan2((pts.obmp[1] - pts.axis.geo[1]), (pts.obmp[0] - pts.axis.geo[0]))
        theta_markers[1] = atan3((xpt_temp[1] - pts.axis.geo[1]), (xpt_temp[0] - pts.axis.geo[0]))
        theta_markers[2] = atan3((pts.ibmp[1] - pts.axis.geo[1]), (pts.ibmp[0] - pts.axis.geo[0]))
        theta_markers[3] = atan3((pts.bottom[1] - pts.axis.geo[1]), (pts.bottom[0] - pts.axis.geo[0]))
        theta_markers[4] = theta_markers[0] + 2 * pi

    else:
        # these theta markers correspond to the obmp, top, ibmp, xpt, and obmp+2pi, respectively
        theta_markers = np.zeros(5)
        theta_markers[0] = atan2((pts.obmp[1] - pts.axis.geo[1]), (pts.obmp[0] - pts.axis.geo[0]))
        theta_markers[1] = atan3((pts.top[1] - pts.axis.geo[1]), (pts.top[0] - pts.axis.geo[0]))
        theta_markers[2] = atan3((pts.ibmp[1] - pts.axis.geo[1]), (pts.ibmp[0] - pts.axis.geo[0]))
        theta_markers[3] = atan3((xpt_temp[1] - pts.axis.geo[1]), (xpt_temp[0] - pts.axis.geo[0]))
        theta_markers[4] = theta_markers[0] + 2 * pi

    try:
        min_delta_theta = 2 * pi / thetapts_approx
    except:
        print('thetapts_approx not defined. Setting to 30')
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