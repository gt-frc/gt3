#!/usr/bin/env python2
# -*- coding: utf-8 -*-

def create_triangle_opts(inp):
    """
    Create input options for Triangle.

    Refer to https://www.cs.cmu.edu/~quake/triangle.html
    :param inp:
    :return:
    """
    tri_options = '-p'
    try:
        tri_options = tri_options + 'q' + str(inp.tri_min_angle)
    except:
        pass

    try:
        tri_options = tri_options + 'a' + str(inp.tri_min_area)
    except:
        pass

    tri_options = tri_options + 'nz'
    return tri_options