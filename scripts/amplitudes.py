# Copyright (C) 2015 Institute for Theoretical Physics, ETH Zurich
#               2015 by Michele Dolfi <dolfim@phys.ethz.ch>
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

import sys
import numpy as np
import scipy
from copy import deepcopy
import pyalps_dset

import utils
from corr_helpers import *
import density


def rsquared(x, y, coeffs):
    # r-squared
    p = np.poly1d(coeffs)
    # fit values, and mean
    yhat = p(x)                         # or [p(z) for z in x]
    ybar = np.sum(y)/len(y)          # or sum(y)/len(y)
    ssreg = np.sum((yhat-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
    return ssreg / sstot


def compute_amplitude(d):
    ret = pyalps_dset.DataSet()
    ret.props = deepcopy(d.props)
    ret.props['observable'] = 'Amplitude'
    ret.y = np.array([ abs(d.props['density_fitted_at_middle'] - d.props['density_fitted_n0']) ])
    return ret


def compute(filling, bond_dim, odd_sizes=False, amplitude_points=None):
    ## Load all system sizes
    if odd_sizes:
        densities = map(lambda parms: to_dataset(*utils.load_or_evaluate('density_fit', density.evaluate_fit, **parms)),
                   utils.iter_odd_system_size(bond_dim=bond_dim, filling=filling))
    else:
        densities = map(lambda parms: to_dataset(*utils.load_or_evaluate('density_fit', density.evaluate_fit, **parms)),
                   utils.iter_system_size(bond_dim=bond_dim, filling=filling))
    densities = filter(lambda d: d.props is not None, densities)
    if len(densities) < 3:
        print 'WARNING:', 'Only got {} valid densities, you need at least 3.'.format(len(densities))
        d = pyalps_dset.DataSet()
        d.props = None
        return d
    
    # check that particle density in the middle is symmetric, otherwise exclude dataset
    # for L=192 it is better to use only M > 2000
    # def symm_middle_filter(d):
    #     x = int(d.props['L'] / 2)
    #     delta = abs(d.y[x] - d.y[x-1]) if d.props['L']%2==0 else abs(d.y[x+1] - d.y[x-1])
    #     if delta > 5e-5:
    #         print '# Discard L={L:.0f}, n={filling}, tperp={t\'}, M={max_bond_dimension:.0f} : {delta}'.format(delta=delta, **d.props)
    #         return False
    #     return True
    # densities = filter(symm_middle_filter, densities)
    
    
    ## Only even number of pairs
    densities = filter(lambda d: d.props['L']%2 != 0 or d.props['Nup_total'] % 2 == 0, densities)

    amplitudes = map(compute_amplitude, densities)
    for d in amplitudes:
        if d.props['density_fit_error'] > abs(d.y[0]):
            print 'WARNING:', 'Error in density fit is higher than amplitude magnitude.'
            print 'Error per site is', d.props['density_fit_error'], 'Amplidue is', d.y[0]
            print 'Dataset was:', 'L=%s filling=%s t_perp=%s M=%s' % (d.props['L'], d.props['filling'], d.props['t\''], d.props['max_bond_dimension'])
    
    common_props = pyalps_dset.dict_intersect([d.props for d in amplitudes])
    amp_vs_L = pyalps_dset.collectXY(amplitudes, 'L', 'Amplitude')
    d = amp_vs_L[0]
    
    ## fit finite size analysis
    sel = np.ones(len(d.x), dtype=bool)
    if amplitude_points is not None and amplitude_points < len(d.x): sel[:-amplitude_points] = False
    
    cov = np.ones((2,2))*np.nan
    if sum(sel) > 2:
        coeff, cov = np.polyfit(np.log(d.x[sel]), np.log(d.y[sel]), deg=1, cov=True)
    else:
        coeff = np.polyfit(np.log(d.x[sel]), np.log(d.y[sel]), deg=1)
    r2 = rsquared(np.log(d.x), np.log(d.y), coeff)
    fit_range = [d.x[0], d.x[-1]]
    if amplitude_points is not None and amplitude_points < len(d.x): fit_range[0] = d.x[-amplitude_points]
    ## compute new Krho
    slope = coeff[0]
    errslope = np.sqrt(cov[0,0])
    Krho = -2. * slope
    errKrho = 2. * errslope
    print 'M={} t_perp={} filling={:4.4f}  ->  Krho={:.3f} +/- {:.3f} :  R^2 {:.4f}'.format(common_props['max_bond_dimension'], common_props['t\''], common_props['filling'], Krho, errKrho, r2)
    
    d.props['fit_range']         = fit_range
    d.props['fit_numpoints']     = amplitude_points
    d.props['fitted_coeff']      = coeff
    d.props['fitted_Krho']       = Krho
    d.props['fitted_Krho_error'] = errKrho
    d.props['fitted_r2']         = r2
    
    return d



def evaluate(**kwargs):
    data = compute(**kwargs)
    if data.props is None:
        return None, None
    d = np.column_stack([data.x,data.y])
    props = data.props
    return d, props
