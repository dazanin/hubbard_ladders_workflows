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
    data = load_raw_data.loadEigenstateMeasurements([fname],
                    what=[
                          'dens corr up-up',
                          'dens corr up-down',
                          'dens corr down-up',
                          'dens corr down-down',
                          'Local density up',
                          'Local density down',
                          ])
    flat_data = pyalps_dset.flatten(data)
    flat_data = filter(lambda d: type(d) != list or (len(d) > 0 and type(d[0]) != list), flat_data)
    if len(flat_data) < 6:
        print 'WARNING:', 'Not enough datasets loaded in', fname
        return None
    
    common_props = pyalps_dset.dict_intersect([d.props for d in pyalps_dset.flatten(data)])
    
    L = int(common_props['L'])
    W = int(common_props['W']) if 'W' in common_props else 2
    idx = index_map(L, W)
    
    try:
        dcor_up_up     = select_to_2d(data, 'dens corr up-up'    , idx)
        dcor_up_down   = select_to_2d(data, 'dens corr up-down'  , idx)
        dcor_down_up   = select_to_2d(data, 'dens corr down-up'  , idx)
        dcor_down_down = select_to_2d(data, 'dens corr down-down', idx)
        dens_up        = select_to_1d(data, 'Local density up'   , idx)
        dens_down      = select_to_1d(data, 'Local density down' , idx)
    except ObservableNotFound:
        print 'WARNING:', 'Measurement not found in', fname
        return None
    
    ix_lower_chain = idx[:,0]
    ix_upper_chain = idx[:,1]
    
    total_dens = dens_up + dens_down
    dens_on_rungs = total_dens[ix_lower_chain] + total_dens[ix_upper_chain]
    
    
    dcor = np.zeros((L,L))
    ## combine density correlators for <N(i)*N(j)> between rungs
    dcor += dcor_up_up    [np.ix_(ix_lower_chain, ix_lower_chain)]
    dcor += dcor_up_down  [np.ix_(ix_lower_chain, ix_lower_chain)]
    dcor += dcor_down_up  [np.ix_(ix_lower_chain, ix_lower_chain)]
    dcor += dcor_down_down[np.ix_(ix_lower_chain, ix_lower_chain)]
    
    dcor += dcor_up_up    [np.ix_(ix_lower_chain, ix_upper_chain)]
    dcor += dcor_up_down  [np.ix_(ix_lower_chain, ix_upper_chain)]
    dcor += dcor_down_up  [np.ix_(ix_lower_chain, ix_upper_chain)]
    dcor += dcor_down_down[np.ix_(ix_lower_chain, ix_upper_chain)]
    
    dcor += dcor_up_up    [np.ix_(ix_upper_chain, ix_lower_chain)]
    dcor += dcor_up_down  [np.ix_(ix_upper_chain, ix_lower_chain)]
    dcor += dcor_down_up  [np.ix_(ix_upper_chain, ix_lower_chain)]
    dcor += dcor_down_down[np.ix_(ix_upper_chain, ix_lower_chain)]
    
    dcor += dcor_up_up    [np.ix_(ix_upper_chain, ix_upper_chain)]
    dcor += dcor_up_down  [np.ix_(ix_upper_chain, ix_upper_chain)]
    dcor += dcor_down_up  [np.ix_(ix_upper_chain, ix_upper_chain)]
    dcor += dcor_down_down[np.ix_(ix_upper_chain, ix_upper_chain)]
    
    ## connected correlator
    dcor += -np.outer(dens_on_rungs, dens_on_rungs)
    
    
    d = pyalps_dset.DataSet()
    d.props = deepcopy(common_props)
    d.props['observable'] = 'Density Correlation'
    d.y = dcor
    d.idx = idx
    
    return d



def compute_extrap(L, filling, bond_dim, correlation_type, at_x=None):
    ## Load all bond dimensions
    sets = map(lambda parms: to_dataset(*utils.load_or_evaluate('densdens', evaluate, **parms)),
               utils.iter_bond_dim(L=L, filling=filling, correlation_type=correlation_type))
    sets = filter(lambda d: d.props is not None, sets)
    
    ## Extrapolate
    return extrapolate_local.extrapolate(sets, 'Density Correlation', foreach=['L', 't\'', 'Nup_total', 'Ndown_total'], extrap_type=bond_dim.type, deg=bond_dim.deg, num_points=bond_dim.num_points, full_output_at=[at_x])


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
    
    props['correlation_type'] = correlation_type
    return d, props


def evaluate_extrap_at(L, filling, bond_dim, correlation_type, at_x):
    extrap, obs_vs_extrap, fits = compute_extrap(L, filling, bond_dim, correlation_type, at_x)
    if len(obs_vs_extrap) == 0: return None, None
    d = np.column_stack([obs_vs_extrap[0].x, obs_vs_extrap[0].y])
    props = obs_vs_extrap[0].props
    return d, props

def evaluate_extrap(L, filling, bond_dim, correlation_type):
    extrap, obs_vs_extrap, fits = compute_extrap(L, filling, bond_dim, correlation_type)
    if len(extrap) == 0: return None, None
    d = np.column_stack([extrap[0].x, extrap[0].y])
    props = extrap[0].props
    return d, props


def evaluate(**kwargs):
    if isinstance(kwargs['bond_dim'], utils.Extrapolation):
        return evaluate_extrap(**kwargs)
    return evaluate_single(**kwargs)

