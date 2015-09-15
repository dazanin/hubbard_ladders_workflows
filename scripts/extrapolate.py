# Copyright (C) 2015 Institute for Theoretical Physics, ETH Zurich
#               2015 by Michele Dolfi <dolfim@phys.ethz.ch>
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

import pyalps_dset
import numpy as np
from copy import deepcopy

def rsquared(x, y, coeffs):
    # r-squared
    p = np.poly1d(coeffs)
    # fit values, and mean
    yhat = p(x)                         # or [p(z) for z in x]
    ybar = np.sum(y)/len(y)          # or sum(y)/len(y)
    ssreg = np.sum((yhat-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
    return ssreg / sstot

def bond_dimension(data, foreach=[], deg=1, num_points=None):
    en_vs_bond = pyalps_dset.collectXY(pyalps_dset.flatten(data), 'max_bond_dimension', 'Energy', foreach=foreach)
    extrap = []
    fits   = []
    for d in en_vs_bond:
        bond_dims = deepcopy(d.x)
        d.x = 1./d.x
        d.x = d.x[::-1]
        d.y = d.y[::-1]
        
        sel = np.ones(len(d.x), dtype=bool)
        if num_points is not None and num_points < len(d.x): sel[num_points:] = False
        coeff = np.polyfit(d.x[sel], d.y[sel], deg=deg)
        r2 = rsquared(d.x[sel], d.y[sel], coeff)
        fit_cut = d.x[num_points-1] if num_points is not None and num_points < len(d.x) else d.x[-1]
        
        d.props['fitted_r2'] = r2
        d.props['fit_deg'] = deg
        d.props['fit_numpoints'] = num_points
        d.props['fit_cut']       = fit_cut
        d.props['fitted_coeff']  = coeff
        d.props['fitted_energy'] = coeff[-1]
        
        dd = pyalps_dset.DataSet()
        dd.props = deepcopy(d.props)
        dd.props['fitted_r2'] = r2
        dd.props['fit_deg'] = deg
        dd.props['fitted_coeff'] = coeff
        dd.x = np.linspace(0, max(d.x))
        dd.y = np.polyval(coeff, dd.x)
        fits.append(dd)
        
        dd = pyalps_dset.DataSet()
        dd.props = deepcopy(d.props)
        dd.props['max_bond_dimension'] = 'inf'
        dd.x = np.array([])
        dd.y = np.array([ coeff[-1] ])
        extrap.append(dd)
        
        d.props['line'] = 'scatter'
        d.props['bond_dims'] = bond_dims
        
    
    return extrap, en_vs_bond, fits


def variance(data, foreach=[], deg=1, num_points=None):
    en_vs_variance = pyalps_dset.collectXY(data, 'EnergyVariance', 'Energy', foreach=foreach)
    extrap = []
    fits   = []
    for d in en_vs_variance:
        sel = np.ones(len(d.x), dtype=bool)
        if num_points is not None and num_points < len(d.x): sel[num_points:] = False
        coeff = np.polyfit(d.x[sel], d.y[sel], deg=deg)
        r2 = rsquared(d.x[sel], d.y[sel], coeff)
        fit_cut = d.x[num_points-1] if num_points is not None and num_points < len(d.x) else d.x[-1]
        
        d.props['fitted_r2'] = r2
        d.props['fit_deg'] = deg
        d.props['fit_numpoints'] = num_points
        d.props['fit_cut']       = fit_cut
        d.props['fitted_coeff']  = coeff
        d.props['fitted_energy'] = coeff[-1]
        
        dd = pyalps_dset.DataSet()
        dd.props = deepcopy(d.props)
        dd.props['fitted_r2'] = r2
        dd.props['fit_deg'] = deg
        dd.props['fitted_coeff'] = coeff
        dd.x = np.linspace(0, max(d.x))
        dd.y = np.polyval(coeff, dd.x)
        fits.append(dd)
        
        dd = pyalps_dset.DataSet()
        dd.props = deepcopy(d.props)
        dd.props['max_bond_dimension'] = 'inf'
        dd.props['EnergyVariance'] = 0.
        dd.x = np.array([])
        dd.y = np.array([ coeff[-1] ])
        extrap.append(dd)
        
        d.props['line'] = 'scatter'
    
    return extrap, en_vs_variance, fits


def truncation(data, foreach=[], deg=1, num_points=None):
    en_vs_truncation = pyalps_dset.collectXY(pyalps_dset.flatten(data), 'TruncatedWeight', 'Energy', foreach=foreach)
    extrap = []
    fits   = []
    for d in en_vs_truncation:
        sel = np.ones(len(d.x), dtype=bool)
        if num_points is not None and num_points < len(d.x): sel[num_points:] = False
        coeff = np.polyfit(d.x[sel], d.y[sel], deg=deg)
        r2 = rsquared(d.x[sel], d.y[sel], coeff)
        fit_cut = d.x[num_points-1] if num_points is not None and num_points < len(d.x) else d.x[-1]
        
        d.props['fitted_r2'] = r2
        d.props['fit_deg'] = deg
        d.props['fit_numpoints'] = num_points
        d.props['fit_cut']       = fit_cut
        d.props['fitted_coeff']  = coeff
        d.props['fitted_energy'] = coeff[-1]
        
        dd = pyalps_dset.DataSet()
        dd.props = deepcopy(d.props)
        dd.props['fitted_r2'] = r2
        dd.props['fit_deg'] = deg
        dd.props['fit_numpoints'] = num_points
        dd.props['fit_cut']       = fit_cut
        dd.props['fitted_coeff']  = coeff
        dd.x = np.linspace(0, max(d.x))
        dd.y = np.polyval(coeff, dd.x)
        fits.append(dd)
        
        dd = pyalps_dset.DataSet()
        dd.props = deepcopy(d.props)
        dd.props['max_bond_dimension'] = 'inf'
        dd.x = np.array([])
        dd.y = np.array([ coeff[-1] ])
        extrap.append(dd)
        
        d.props['line'] = 'scatter'
    
    return extrap, en_vs_truncation, fits


def extrapolate(data, foreach=[], extrap_type='variance', deg=2, num_points=None):
    
    if extrap_type == 'variance':
        return variance(data, foreach, deg=deg, num_points=num_points)
    elif extrap_type == 'bonddim':
        return bond_dimension(data, foreach, deg=deg, num_points=num_points)
    elif extrap_type == 'truncation':
        return truncation(data, foreach, deg=deg, num_points=num_points)
    else:
        raise Exception('Wrong extrapolation type `%s`' % extrap_type)

