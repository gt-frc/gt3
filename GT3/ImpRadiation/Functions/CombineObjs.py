#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from collections import namedtuple

def combine_objs(obj1, obj2):
    """Combines dictionaries of two instantiated classes or namedtuples into a new namedtuple"""

    #this function might need some work, just FYI.
    try:
        dict1 = obj1._asdict()
    except AttributeError:
        dict1 = obj1.__dict__

    try:
        dict2 = obj2._asdict()
    except AttributeError:
        dict2 = obj2.__dict__

    temp = dict2.copy()
    new_dict = dict1.update(dict2)
    new_nt = namedtuple('new', list(new_dict.keys()))(*list(new_dict.values()))
    return new_nt
