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
import extrapolate_local


def get_density(data):
    def mean_density_at(j):
        amp = 0.
        all_values = np.array([])
        for d in pyalps_dset.flatten(data):
            all_values = np.concatenate([ all_values, d.y[0][( d.x==j )[:,0]] ])
        amp = np.sum(all_values)
        std = np.std(all_values)
        return amp, std
    
    ret = pyalps_dset.DataSet()
    ret.props = pyalps_dset.dict_intersect([d.props for d in pyalps_dset.flatten(data)])
    ret.props['observable'] = 'Rung density'
    
    bond_dim = int(ret.props['max_bond_dimension'])
    L = int(ret.props['L'])
    Leff = L if L%2 == 0 else L-1
    nup = ret.props['Nup_total']
    nholes = L - nup
    if nup % 2 == 0: ## for the case with two additional full sites
        filling = float(nup) / Leff
    else:
        filling = float(nup-1) / Leff
    ret.props['nholes']  = nholes
    ret.props['filling'] = filling
    
    density = np.empty((L,2))
    for j in xrange(L):
        density[j] = mean_density_at(j)
    x = np.arange(0.5, L, 1)
    
    ret.x = x
    ret.y = density[:,0]
    
    return ret
    

def dens_fit(d, ff):
    L = d.props['L']
    bond_dim = d.props['max_bond_dimension']
    filling = d.props['filling']
    nholes = d.props['nholes']
    
    xgrid=np.linspace(0.5,L-0.5,1000)
    
    doloop = True
    delta_sel = L/4.
    while doloop:
        sel = np.array([ x > L/2. - delta_sel and x < L/2. + delta_sel for x in d.x])
        fw.fit(None, lambda self, x, args: ff.func(x, args), ff.pars, d.y[sel], d.x[sel])
        
        amp = abs(ff.func(L/2., ff.pars) - ff.n0())
        if amp < 1e-8:
            delta_sel /= 1.5
            ff.reset()
            print 'WARNING:', 'Reducing fit range.', 'Amplitude too small:', amp, 'Very likely the fitting routine has problems converging the result.'
            if delta_sel < L/12.:
                doloop=False
                print 'PROBLEM:', 'Failed to fit.'
        else:
            doloop = False
    
    sel = np.array([ x > L/2. - delta_sel and x < L/2. + delta_sel for x in d.x])
    chi2 = sum( (d.y[sel] - ff.func(d.x[sel], ff.pars))**2 )
    error_per_point = np.sqrt(chi2 / sum(sel))
    
    n0    = ff.n0()
    K_rho = 2*ff.alpha()
    
    # print '============='
    print '## FIT `{}` ##'.format(ff.name)
    print '#', 'L      :', L
    print '#', 'fill   :', filling
    print '#', 'nholes :', nholes
    # print 'n0    :', n0
    print 'K_rho :', K_rho
    
    for n,p in zip(ff.parname, ff.pars):
        print '--', '{} : {}'.format(n, p())

    dd = pyalps_dset.DataSet()
    dd.props = deepcopy(d.props)
    dd.x = xgrid
    dd.y = ff.func(xgrid, ff.pars)
    
    dd.props['density_fit_function'] = ff.name
    dd.props['density_fit_sel']      = np.array([L/2.-delta_sel, L/2.+delta_sel])
    dd.props['density_fit_chi2']     = chi2
    dd.props['density_fit_error']    = error_per_point
    dd.props['density_fitted_n0']    = n0
    dd.props['density_fitted_Krho']  = K_rho
    dd.props['density_fitted_at_middle'] = ff.func(L/2., ff.pars)
    
    return dd

class fit_func:
    name = 'fit_func'
    def __init__(self, L, nholes):
        self.L = float(L)
        self.nholes = float(nholes)
        self.reset()
    def reset(self):
        self.pars = [fw.Parameter(self.nholes/self.L),fw.Parameter(1.),fw.Parameter(0.5)]
        self.parname = ['n0', 'A', 'alpha']
    
    def n0(self):
        return self.pars[0]()
    def alpha(self):
        return self.pars[2]()
    
    def func(self, x, args):
        n0 = args[0]()
        A = args[1]()
        alpha = args[2]()
        L_eff = self.L - 2.
        kk = self.nholes / L_eff
        shift = - np.pi*self.nholes * self.L / L_eff
        if self.L % 2 != 0:
            shift += np.pi
        shift2 = np.pi/2 * (1. - self.L / L_eff)
        return A * np.cos(2*np.pi*kk*x + shift) / ((2*L_eff/np.pi) * np.sin(np.pi*x/L_eff + shift2))**alpha + n0


def compute(fname, nup_name='Local density up', ndown_name='Local density down'):
    data = load_raw_data.loadEigenstateMeasurements([fname], what=[nup_name, ndown_name])
    d = get_density(data)
    return d

def compute_fit(L, filling, bond_dim):
    ## Get data
    data = to_dataset(*utils.load_or_evaluate('density', evaluate, L=L, filling=filling, bond_dim=bond_dim))
    if data.props is None:
        print 'Data not found for L={}, n={}, M={}'.format(L, filling, bond_dim)
        return data
    ## Compute fit
    d = dens_fit(data, fit_func(data.props['L'], data.props['nholes']))
    return d




def compute_extrap(L, filling, bond_dim, at_x=None):
    ## Load all bond dimensions
    sets = map(lambda parms: to_dataset(*utils.load_or_evaluate('density', evaluate, **parms)),
               utils.iter_bond_dim(L=L, filling=filling))
    sets = filter(lambda d: d.props is not None, sets)
    
    # check that particle density in the middle is symmetric, otherwise exclude dataset
    # for L=192 it is better to use only M > 2000
    def symm_middle_filter(d):
        x = int(d.props['L'] / 2)
        delta = abs(d.y[x] - d.y[x-1]) if d.props['L']%2==0 else abs(d.y[x+1] - d.y[x-1])
        if delta > 5e-5:
            print '# Discard L={L:.0f}, n={filling}, tperp={t\'}, M={max_bond_dimension:.0f} : middle density not symmetry, delta={delta}'.format(delta=delta, **d.props)
            return False
        return True
    sets = filter(symm_middle_filter, sets)
    
    
    ## Extrapolate
    return extrapolate_local.extrapolate(sets, 'Rung density', foreach=['L', 't\'', 'Nup_total', 'Ndown_total'], extrap_type=bond_dim.type, deg=bond_dim.deg, num_points=bond_dim.num_points, full_output_at=[at_x])


def evaluate_single(L, filling, bond_dim):
    ## Get filenames
    fname = utils.find_resfile(L, filling, bond_dim)
    if fname is None:
        print 'No data available for L={}, n={}, M={}'.format(L, filling, bond_dim)
        return None, None
    
    dens = compute(fname)
    
    d = np.column_stack([dens.x,dens.y])
    props = dens.props
    return d, props


def evaluate_extrap_at(L, filling, bond_dim, at_x):
    extrap, obs_vs_extrap, fits = compute_extrap(L, filling, bond_dim, at_x+0.5)
    if len(obs_vs_extrap) == 0: return None, None
    d = np.column_stack([obs_vs_extrap[0].x, obs_vs_extrap[0].y])
    props = obs_vs_extrap[0].props
    props['extrap_type'] = str(bond_dim)
    return d, props

def evaluate_extrap(L, filling, bond_dim):
    extrap, obs_vs_extrap, fits = compute_extrap(L, filling, bond_dim)
    if len(extrap) == 0: return None, None
    d = np.column_stack([extrap[0].x, extrap[0].y])
    props = extrap[0].props
    props['extrap_type'] = str(bond_dim)
    return d, props


def evaluate(**kwargs):
    if isinstance(kwargs['bond_dim'], utils.Extrapolation):
        return evaluate_extrap(**kwargs)
    return evaluate_single(**kwargs)

def evaluate_fit(**kwargs):
    data = compute_fit(**kwargs)
    if data.props is None:
        return None, None
    d = np.column_stack([data.x,data.y])
    props = data.props
    return d, props

