# Copyright (C) 2015 Institute for Theoretical Physics, ETH Zurich
#               2015 by Michele Dolfi <dolfim@phys.ethz.ch>
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

from os import path
import utils


def result(what, **kwargs):
    fname = utils.filename(what, **kwargs)
    if path.exists(fname):
        return utils.load(fname)
    else:
        try:
            import pyalps_dset
            import pairfield_correlations, density_correlations, density, amplitudes
            if   what == 'pairfield':
                return utils.load_or_evaluate(what, pairfield_correlations.evaluate, **kwargs)
            elif what == 'densdens':
                return utils.load_or_evaluate(what, density_correlations.evaluate, **kwargs)
            elif what == 'density':
                return utils.load_or_evaluate(what, density.evaluate, **kwargs)
            elif what == 'density_fit':
                return utils.load_or_evaluate(what, density.evaluate_fit, **kwargs)
            elif what == 'density_amplitudes':
                return utils.load_or_evaluate(what, amplitudes.evaluate, **kwargs)
            else:
                raise Exception('`%s` is not a valid measurement.' % what)
            
        except ImportError, e:
            print 'To execute new evaluations you need the ALPS.Python library.'
            print e

def extrapolation(what, **kwargs):
    what = 'extrap_'+what
    fname = utils.filename(what, **kwargs)
    if path.exists(fname):
        return utils.load(fname)
    else:
        import pairfield_correlations, density_correlations, density, energy
        if   what == 'extrap_pairfield':
            return utils.load_or_evaluate(what, pairfield_correlations.evaluate_extrap_at, **kwargs)
        elif what == 'extrap_densdens':
            return utils.load_or_evaluate(what, density_correlations.evaluate_extrap_at, **kwargs)
        elif what == 'extrap_density':
            return utils.load_or_evaluate(what, density.evaluate_extrap_at, **kwargs)
        elif what == 'extrap_energy':
            return utils.load_or_evaluate(what, energy.evaluate_extrap, **kwargs)
        else:
            raise Exception('`%s` is not a valid measurement.' % what)
            
        # except ImportError, e:
        #     print 'To execute new evaluations you need the ALPS.Python library.'
        #     print e

def krho(**kwargs):
    _,props = result('density_amplitudes', **kwargs)
    return props['fitted_Krho']

