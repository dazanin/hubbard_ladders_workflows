# Copyright (C) 2015 Institute for Theoretical Physics, ETH Zurich
#               2015 by Michele Dolfi <dolfim@phys.ethz.ch>
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

import numpy as np
import pyalps_dset

class ObservableNotFound(Exception):
    pass

def select_obs(sets, obs):
    sel = pyalps_dset.select(pyalps_dset.flatten(sets), lambda d: d.props['observable'] == obs)
    if len(sel) == 0:
        raise ObservableNotFound()
    return sel[0]

def index_map(L, W):
    a = np.arange(L*W, dtype=int)
    return a.reshape((L, W))

def select_to_1d(sets, obs, idx):
    d = select_obs(sets, obs)
    a = np.zeros((idx.size,))
    for i, x in enumerate(d.x):
        a[idx[x]] = d.y[0][i]
    return a

def select_to_2d(sets, obs, idx):
    d = select_obs(sets, obs)
    a = np.zeros((idx.size,idx.size))
    for i, x in enumerate(d.x):
        a[idx[tuple(x[0])], idx[tuple(x[1])]] = d.y[0][i]
        a[idx[tuple(x[1])], idx[tuple(x[0])]] = d.y[0][i]
    return a

def select_to_nd(sets, obs, idx, dims):
    print '...','loading ',obs,'...'
    d = select_obs(sets, obs)
    a = np.zeros([idx.size]*dims)
    for i, x in enumerate(d.x):
        key = tuple([ idx[tuple(xi)] for xi in x.reshape(x.size/2,2) ])
        a[key] = d.y[0][i]
    return a

def average_around_middle(corr, L, shifts):
    l = range(1, int(L-1-2*max(shifts)-2), 1)
    res = np.zeros_like(l, dtype=float)
    for ii,li in enumerate(l):
        tmp = []
        for ss in shifts:
            i1 = int(L/2. - li/2. + ss)
            i2 = int(i1 + li)
            tmp.append( corr[i1,i2] )
        res[ii] = np.mean(tmp)
    return l, res

def to_dataset(d, props):
    dd = pyalps_dset.DataSet()
    dd.props = props
    if d is not None and d.size != 0:
        if len(d.shape) == 2:
            dd.x = d[:,0]
            dd.y = d[:,1]
        else:
            dd.y = np.array([d[1]])
    return dd