# Copyright (C) 2015 Institute for Theoretical Physics, ETH Zurich
#               1994-2009 by Bela Bauer <bauerb@phys.ethz.ch>
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

import numpy as np

def dict_intersect(dicts):
    """ computes the intersection of a list of dicts
    
        this function takes a list of dicts as input and returns a dict containing all those key-value pairs that appear with identical values in all dicts 
    """
    sets = [set(q.keys()) for q in dicts]
    intersection = sets[0]
    for iset in sets:
        intersection &= iset
    ret = {}
    for key in intersection:
        take = True
        val0 = dicts[0][key]
        for idict in dicts:
            try:
                if val0 != idict[key]:
                    take = False
            except:
                if np.all(val0 != idict[key]):
                    take = False
        if take:
            ret[key] = dicts[0][key]
    return ret

def dict_difference(dicts):
    sets = [set(q.keys()) for q in dicts]
    intersection = sets[0]
    for iset in sets:
        intersection &= iset
    ret = []
    for key in intersection:
        take = True
        val0 = dicts[0][key]
        for idict in dicts:
            if val0 != idict[key]:
                take = False
        if not take:
            ret.append(key)
    return ret
