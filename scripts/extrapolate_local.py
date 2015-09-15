# Copyright (C) 2015 Institute for Theoretical Physics, ETH Zurich
#               2015 by Michele Dolfi <dolfim@phys.ethz.ch>
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

import pyalps_dset
import numpy as np
from copy import deepcopy
import warnings


def rsquared(x, y, coeffs):
    # r-squared
    p = np.poly1d(coeffs)
    # fit values, and mean
    yhat = p(x)                         # or [p(z) for z in x]
    ybar = np.sum(y)/len(y)          # or sum(y)/len(y)
    ssreg = np.sum((yhat-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
    return ssreg / sstot

def extrapolate_with_error(x, y, deg, num_points):
    sel = np.ones(len(x), dtype=bool)
    if num_points is not None and num_points < len(x): sel[num_points:] = False
    coeff = np.polyfit(x[sel], y[sel], deg=deg)
    val = coeff[-1]
    err = abs(val - y[0]) * .5 # error defined as 50% of the distance from the last point
    r2 = rsquared(x, y, coeff)
    
    xfit = np.linspace(0, max(x))
    yfit = np.polyval(coeff, xfit)
    
    return val, err, r2, coeff, xfit, yfit


def do_extrapolate(data, obs, xgetter, res_props, foreach, deg, num_points, full_output_at=[]):
    extrap = []
    fits   = []
    obs_vs_extrap = []
    
    groups = pyalps_dset.groupSets(data, for_each=foreach)
    for gg in groups:
        if (num_points is not None and num_points < len(gg) and num_points < deg+1) or len(gg) < deg+1:
            print 'WARNING:', 'Extrapolation not possible.', 'len() < deg+1, len={}, deg={}'.format(len(gg), deg)
            continue
        
        common_props = pyalps_dset.dict_intersect([d.props for d in gg])
        
        extrap_x    = []
        bond_dims   = []
        observables = []
        xval = None
        for d in gg:
            # get x values
            if   xval is None: xval = d.x
            elif np.all(abs(d.x - xval) > 1e-10):  raise Exception('`x` values do not match between the extrapolation group.')
            # get y values
            observables.append( d.y )
            # get variance
            extrap_x.append( xgetter(d) )
            bond_dims.append(d.props['max_bond_dimension'])
        extrap_x    = np.array(extrap_x)
        bond_dims   = np.array(bond_dims)
        observables = np.array(observables)
        
        order = np.argsort(extrap_x)
        extrap_x    = extrap_x[order]
        bond_dims   = bond_dims[order]
        observables = observables[order]
        
        dd = pyalps_dset.DataSet()
        dd.props = deepcopy(common_props)
        dd.props['fit_deg'] = deg
        dd.props['fit_numpoints'] = num_points
        dd.x = deepcopy(xval)
        dd.y = [0.]*len(xval)
        for i, xi in enumerate(xval):
            # warnings.filterwarnings('error', category=np.RankWarning)
                # extrapolate x_i
            val, err, r2, coeff, xfit, yfit = extrapolate_with_error(extrap_x, observables[:,i], deg, num_points)
            dd.y[i] = val
            # except np.RankWarning:
            #     print 'RankWarning: num_points={}, deg={}, xi={}'.format(num_points, deg, xi)
            #     print '  ', ' '.join(['{}={}'.format(k,common_props[k]) for k in foreach])
            #     dd.y[i] = np.nan
            #     continue
            
            if xi in full_output_at:
                fit_cut = extrap_x[num_points-1] if num_points is not None and num_points < len(extrap_x) else extrap_x[-1]
                
                dfit = pyalps_dset.DataSet()
                dfit.props = deepcopy(common_props)
                dfit.props['fitted_x']  = xi
                dfit.props['fitted_r2'] = r2
                dfit.props['fit_deg'] = deg
                dfit.props['fit_numpoints'] = num_points
                dfit.props['fit_cut']       = fit_cut
                dfit.props['fitted_coeff'] = coeff
                dfit.x     = xfit
                dfit.y     = yfit
                fits.append(dfit)
            
                dvals = pyalps_dset.DataSet()
                dvals.props = deepcopy(common_props)
                dvals.props['fitted_x'] = xi
                dvals.props['line'] = 'scatter'
                dvals.props['bond_dims'] = deepcopy(bond_dims)
                dvals.props['fit_deg'] = deg
                dvals.props['fit_numpoints'] = num_points
                dvals.props['fit_cut']       = fit_cut
                dvals.props['fitted_coeff'] = coeff
                dvals.props['fitted_r2'] = r2
                dvals.x     = deepcopy(extrap_x)
                dvals.y     = observables[:,i]
                obs_vs_extrap.append(dvals)
        
        dd.props.update(res_props)
        extrap.append(dd)
    
    return extrap, obs_vs_extrap, fits



def bond_dimension(data, obs, foreach=[], deg=2, num_points=None, full_output_at=[]):
    get_inv_bond_dim = lambda ss: 1./ss.props['max_bond_dimension']
    props = {
        'max_bond_dimension' : 'inf',
    }
    return do_extrapolate(data=data, obs=obs, xgetter=get_inv_bond_dim, res_props=props, foreach=foreach, deg=deg, num_points=num_points, full_output_at=full_output_at)


def variance(data, obs, foreach=[], deg=2, num_points=None, full_output_at=[]):
    props = {
        'max_bond_dimension' : 'inf',
        'EnergyVariance'     : 0.,
    }
    return do_extrapolate(data=data, obs=obs, xgetter=lambda ss: ss.props['EnergyVariance'], res_props=props, foreach=foreach, deg=deg, num_points=num_points, full_output_at=full_output_at)


def truncation(data, obs, foreach=[], deg=2, num_points=None, full_output_at=[]):
    props = {
        'max_bond_dimension' : 'inf',
        'TruncatedWeight'    : 0.,
    }
    return do_extrapolate(data=data, obs=obs, xgetter=lambda ss: ss.props['TruncatedWeight'], res_props=props, foreach=foreach, deg=deg, num_points=num_points, full_output_at=full_output_at)


def extrapolate(data, obs, foreach=[], extrap_type='variance', deg=2, num_points=None, full_output_at=[]):
    
    if extrap_type == 'variance':
        return variance(data, obs, foreach, deg=deg, num_points=num_points, full_output_at=full_output_at)
    elif extrap_type == 'bonddim':
        return bond_dimension(data, obs, foreach, deg=deg, num_points=num_points, full_output_at=full_output_at)
    elif extrap_type == 'truncation':
        return truncation(data, obs, foreach, deg=deg, num_points=num_points, full_output_at=full_output_at)
    else:
        raise Exception('Wrong extrapolation type `%s`' % extrap_type)
