# Copyright (C) 2015 Institute for Theoretical Physics, ETH Zurich
#               2015 by Michele Dolfi <dolfim@phys.ethz.ch>
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

import numpy as np
import scipy
from copy import deepcopy
import pyalps_dset

import pyalps_dset.fit_wrapper as fw

import utils, load_raw_data
from corr_helpers import *
import extrapolate



def compute(fname):
    data = load_raw_data.loadEigenstateMeasurements([fname], what='Energy')
    return data[0][0]


def compute_extrap(L, filling, bond_dim):
    ## Load all bond dimensions
    sets = map(lambda parms: to_dataset(*utils.load_or_evaluate('energy', evaluate, **parms)),
               utils.iter_bond_dim(L=L, filling=filling))
    sets = filter(lambda d: d.props is not None, sets)
    
    ## Extrapolate
    return extrapolate.extrapolate(sets, foreach=['L', 't\'', 'Nup_total', 'Ndown_total'], extrap_type=bond_dim.type, deg=bond_dim.deg, num_points=bond_dim.num_points)


def evaluate_single(L, filling, bond_dim):
    ## Get filenames
    fname = utils.find_resfile(L, filling, bond_dim)
    if fname is None:
        print 'No data available for L={}, n={}, M={}'.format(L, filling, bond_dim)
        return None, None
    
    en = compute(fname)
    
    d = np.column_stack([en.x,en.y])
    props = en.props
    return d, props


def evaluate_extrap(L, filling, bond_dim):
    extrap, obs_vs_extrap, fits = compute_extrap(L, filling, bond_dim)
    if len(obs_vs_extrap) == 0: return None, None
    d = np.column_stack([obs_vs_extrap[0].x, obs_vs_extrap[0].y])
    props = obs_vs_extrap[0].props
    props['extrap_type'] = str(bond_dim)
    return d, props


def evaluate(**kwargs):
    if isinstance(kwargs['bond_dim'], utils.Extrapolation):
        return evaluate_extrap(**kwargs)
    return evaluate_single(**kwargs)

