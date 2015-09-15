# Copyright (C) 2015 Institute for Theoretical Physics, ETH Zurich
#               2015 by Michele Dolfi <dolfim@phys.ethz.ch>
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

import numpy as np
from copy import deepcopy
from os import path


PROJECT_ROOT = path.dirname(path.dirname(path.abspath(__file__)))
ENABLE_EXTRAPOLATION_CACHE = False
VERBOSE_LOADING = False
VERBOSE_EMPTY_DATASET = False

def filename(what, **kwargs):
    fname = what
    if 'correlation_type' in kwargs:
        fname += '_{correlation_type}'
    if 'odd_sizes' in kwargs:
        fname += '_Lodd'
    if 'L' in kwargs:
        fname += '_L{L}'
    if 'filling' in kwargs:
        fname += '_n{filling}'
    if 'bond_dim' in kwargs:
        fname += '_M{bond_dim}'
    if 'at_x' in kwargs:
        fname += '_at{at_x:.0f}'
    if 'amplitude_points' in kwargs:
        fitnum = 'All' if kwargs['amplitude_points'] is None else kwargs['amplitude_points']
        fname += '_ampl_fit{fitnum}'.format(fitnum=fitnum)
    
    fname += '.txt'
    return path.join(path.join(PROJECT_ROOT,'data_extracted'),fname.format(**kwargs))

import re
reFloat = r'([+-]?\d+(?:\.\d+)?(?:[eE][+-]\d+)?)'
fit_float_pattern = re.compile(r'# ([^=]+)\s+= %s' % reFloat)
fit_string_pattern = re.compile(r'# ([^=]+)\s+= "([^"]+)"')
fit_array_pattern = re.compile(r'# ([^=]+)\s+= \[([^\]]+)\]')
def load_props(fname):
    res = {}
    for line in open(fname):
        match = fit_float_pattern.search(line)
        if match:
            name = match.group(1).strip()
            val  = float(match.group(2))
            res[name] = val
            continue
        match = fit_string_pattern.search(line)
        if match:
            name = match.group(1).strip()
            val  = match.group(2).strip().strip('"')
            res[name] = val
            continue
        match = fit_array_pattern.search(line)
        if match:
            name = match.group(1).strip()
            val  = np.array([float(vi) for vi in match.group(2).split(',')])
            res[name] = val
            continue
    if len(res) == 0:
        return None
    return res

def load(fname):
    props = load_props(fname)
    d = None
    if props is not None:
        d = np.loadtxt(fname)
    elif VERBOSE_EMPTY_DATASET:
        print 'Warning: file {} does not contain any dataset.'.format(fname)
    return d, props

def save(fname, d, props):
    import sys, datetime
    with open(fname, 'w') as ff:
        ff.write('# Generated on {}\n'.format(datetime.datetime.now()))
        ff.write('# Command line : {}\n'.format(sys.argv))
        if d is not None:
            for k,v in sorted(props.items()):
                if isinstance(v,str):
                    ff.write('# {} = "{}"\n'.format(k,v))
                elif isinstance(v, np.ndarray):
                    ff.write('# {} = [{}]\n'.format( k, ', '.join([str(vi) for vi in v]) ))
                else:
                    ff.write('# {} = {}\n'.format(k,v))
        if d is not None:
            np.savetxt(ff, d)
    return d, props

def load_or_evaluate(what, evaluator, **kwargs):
    fname = filename(what, **kwargs)
    if path.exists(fname):
        d,props = load(fname)
        if VERBOSE_LOADING:
            print 'Loaded', fname
        return load(fname)
    elif ENABLE_EXTRAPOLATION_CACHE or not 'bond_dim' in kwargs or not isinstance(kwargs['bond_dim'], Extrapolation):
        return save(fname, *evaluator(**kwargs))
    return evaluator(**kwargs)

def find_resfile(L, filling, bond_dim, t_perp=1.0):
    from glob import glob
    
    if L % 2 == 0:
        N = int(filling * L)
        if N / filling != L: return None ## L is not a nice multiples
    else:
        N = int(filling * (L-1)) + 1
        if (N-1) / filling != L-1: return None ## L is not a nice multiples
    
    tp = str(t_perp).replace('.','')
    
    flist = glob(PROJECT_ROOT + '/data_raw/L{L}Nu{N}Nd{N}/t{tp}U8/*M{M}.out.res.h5'.format(L=L, N=N, M=bond_dim, tp=tp))
    if len(flist) == 0: return None
    return sorted(flist)[-1] # TODO: temporary fix to get sim with more nsweeps
    # if len(flist) != 1: return None
    # return flist[0]


## Extrapolation types
class Extrapolation(object):
    def __init__(self, extrap_type, deg, num_points):
        self.type = str(extrap_type)
        self.deg = int(deg)
        self.num_points = num_points
    
    def __str__(self):
        nump = self.num_points or 'All'
        return 'extrap_{}_deg{}_num{}'.format(self.type, self.deg, nump)


## Correlation types
class Averaged(object):
    def __str__(self):
        return 'avg'

class FixedStart(object):
    def __init__(self, start):
        self.start = int(start)

    
    def __str__(self):
        return 'start{}'.format(self.start)

## Iterate all parameters for analysis
def iter_bond_dim(**kwargs):
    all_bond_dim = [1200, 1600, 2000, 2800, 3200, 3600, 4000, 4800]
    ll = [dict(kwargs.items() + [('bond_dim',m)]) for m in all_bond_dim]
    return iter(ll)
def iter_system_size(**kwargs):
    all_L = [32, 48, 64, 80, 96, 128, 160, 192]
    ll = [dict(kwargs.items() + [('L',li)]) for li in all_L]
    return iter(ll)
def iter_odd_system_size(**kwargs):
    all_L = [33, 49, 65, 81, 97, 129]
    ll = [dict(kwargs.items() + [('L',li)]) for li in all_L]
    return iter(ll)

