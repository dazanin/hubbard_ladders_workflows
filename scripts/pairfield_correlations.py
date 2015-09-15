# Copyright (C) 2015 Institute for Theoretical Physics, ETH Zurich
#               2015 by Michele Dolfi <dolfim@phys.ethz.ch>
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

import numpy as np
from copy import deepcopy

import pyalps_dset

import utils, load_raw_data
from corr_helpers import *
import extrapolate_local

def compute(fname):
    print 'loading', fname
    data = load_raw_data.loadEigenstateMeasurements([fname],
                    what=[
                          'pair field 1',
                          'pair field 2',
                          'pair field 3',
                          'pair field 4',
                          ])
    flat_data = pyalps_dset.flatten(data)
    flat_data = filter(lambda d: type(d) != list or (len(d) > 0 and type(d[0]) != list), flat_data)
    if len(flat_data) < 4:
        print 'WARNING:', 'Not enough datasets loaded in', fname
        return None
        
    common_props = pyalps_dset.dict_intersect([d.props for d in pyalps_dset.flatten(data)])
    
    L = int(common_props['L'])
    W = int(common_props['W']) if 'W' in common_props else 2
    idx = index_map(L, W)
    
    try:
        pairfield_1 = select_to_nd(data, 'pair field 1', idx, dims=4)
        pairfield_2 = select_to_nd(data, 'pair field 2', idx, dims=4)
        pairfield_3 = select_to_nd(data, 'pair field 3', idx, dims=4)
        pairfield_4 = select_to_nd(data, 'pair field 4', idx, dims=4)
    except ObservableNotFound:
        print 'WARNING:', 'Measurement not found in', fname
        return None
    
    ix_lower = idx[:,0].reshape(L,1)
    ix_upper = idx[:,1].reshape(L,1)
    jx_lower = idx[:,0].reshape(1,L)
    jx_upper = idx[:,1].reshape(1,L)
    
    corr = np.zeros((L,L))
    
    corr += +pairfield_1[ ix_lower, ix_upper, jx_lower, jx_upper ]
    corr += -pairfield_2[ ix_lower, ix_upper, jx_lower, jx_upper ]
    corr += -pairfield_3[ ix_lower, ix_upper, jx_lower, jx_upper ]
    corr += +pairfield_4[ ix_lower, ix_upper, jx_lower, jx_upper ]


    d = pyalps_dset.DataSet()
    d.props = deepcopy(common_props)
    d.props['observable'] = 'Pairfield Correlation'
    d.y = corr
    d.idx = idx

    return d

def compute_extrap(L, filling, bond_dim, correlation_type, at_x=None):
    ## Load all bond dimensions
    sets = map(lambda parms: to_dataset(*utils.load_or_evaluate('pairfield', evaluate, **parms)),
               utils.iter_bond_dim(L=L, filling=filling, correlation_type=correlation_type))
    sets = filter(lambda d: d.props is not None, sets)
    
    ## Extrapolate
    return extrapolate_local.extrapolate(sets, 'Pairfield Correlation', foreach=['L', 't\'', 'Nup_total', 'Ndown_total'], extrap_type=bond_dim.type, deg=bond_dim.deg, num_points=bond_dim.num_points, full_output_at=[at_x])


def evaluate_single(L, filling, bond_dim, correlation_type):
    ## Get filenames
    fname = utils.find_resfile(L, filling, bond_dim)
    if fname is None:
        print 'No data available for L={}, n={}, M={}, {}'.format(L, filling, bond_dim, correlation_type)
        return None, None
    
    ## Get correlation matrix
    corr = compute(fname)
    if corr is None:
        print 'No data available for L={}, n={}, M={}, {}'.format(L, filling, bond_dim, correlation_type)
        return None, None
    
    if isinstance(correlation_type, utils.FixedStart):
        start_site = correlation_type.start
        y = corr.y[start_site, start_site+1:]
        x = np.arange(1, len(y)+1, dtype=float)
        d = np.column_stack([x, y])
        props = deepcopy(corr.props)
    elif isinstance(correlation_type, utils.Averaged):
        shifts = range(-5, 6)
        x,y = average_around_middle(corr.y, corr.props['L'], shifts)
        d = np.column_stack([x, y])
        props = deepcopy(corr.props)
    
    props['correlation_type'] = str(correlation_type)
    return d, props


def evaluate_extrap_at(L, filling, bond_dim, correlation_type, at_x):
    extrap, obs_vs_extrap, fits = compute_extrap(L, filling, bond_dim, correlation_type, at_x)
    if len(obs_vs_extrap) == 0: return None, None
    d = np.column_stack([obs_vs_extrap[0].x, obs_vs_extrap[0].y])
    props = obs_vs_extrap[0].props
    props['extrap_type'] = str(bond_dim)
    return d, props

def evaluate_extrap(L, filling, bond_dim, correlation_type):
    extrap, obs_vs_extrap, fits = compute_extrap(L, filling, bond_dim, correlation_type)
    if len(extrap) == 0: return None, None
    d = np.column_stack([extrap[0].x, extrap[0].y])
    props = extrap[0].props
    props['extrap_type'] = str(bond_dim)
    return d, props


def evaluate(**kwargs):
    if isinstance(kwargs['bond_dim'], utils.Extrapolation):
        return evaluate_extrap(**kwargs)
    return evaluate_single(**kwargs)

