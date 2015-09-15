# Copyright (C) 2015 Institute for Theoretical Physics, ETH Zurich
#               2015 by Michele Dolfi <dolfim@phys.ethz.ch>
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

import numpy as np

def load_variance_for_dset(ss):
    try:
        import pyalps
    except ImportError, e:
        print 'ERROR: To extract new observbales from the raw data you need the ALPS.Python library.'
        raise e

    variance = pyalps.loadEigenstateMeasurements([ss.props['filename']], what=['EnergyVariance'])
    if len(variance) < 1 or len(variance[0]) < 1:
        raise Exception('EnergyVariance not found in', ss.props['filename'])
    return variance[0][0].y[0]

def load_truncated_weight_for_dset(ss):
    try:
        import pyalps
    except ImportError, e:
        print 'ERROR: To extract new observbales from the raw data you need the ALPS.Python library.'
        raise e

    ar = pyalps.hdf5.archive(ss.props['filename'])
    try:
        if 'simulation' in ar.list_children('/'):
            iteration_path = '/simulation/iteration'
        else:
            iteration_path = '/spectrum/iteration'
        sweeps = ar.list_children(iteration_path)
        sweeps = [int(s) for s in sweeps]
        max_sweep = max(sweeps)
        
        truncated_weight = ar[iteration_path+'/'+str(max_sweep)+'/results/TruncatedWeight/mean/value']
        ss.props['TruncatedWeight'] = max(truncated_weight)
        return ss.props['TruncatedWeight']
    except Exception as e:
        print 'Warning:', 'no TruncatedWeight found in', ss.props['filename']
        print e

def loadEigenstateMeasurements(*args, **kwargs):
    try:
        import pyalps
    except ImportError, e:
        print 'ERROR: To extract new observbales from the raw data you need the ALPS.Python library.'
        raise e
    
    data = pyalps.loadEigenstateMeasurements(*args, **kwargs)
    for d in pyalps.flatten(data):
        if isinstance(d, pyalps.DataSet):
            d.props['TruncatedWeight'] = load_truncated_weight_for_dset(d)
            d.props['EnergyVariance']  = load_variance_for_dset(d)
    
    return data
