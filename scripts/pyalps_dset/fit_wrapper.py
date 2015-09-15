# Copyright (C) 2015 Institute for Theoretical Physics, ETH Zurich
#               1994-2009 by Bela Bauer <bauerb@phys.ethz.ch>
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

import numpy as np
from scipy import optimize


__all__ = ['fit', 'Parameter']

class Parameter:    
    def __init__(self, value):
        self.value = value

    def set(self, value):
        self.value = value

    def get(self):
        return self.value

    def __call__(self):
        return self.value

def fit(self,function, parameters, y, x = None):
    def f(params):
        i = 0
        for p in parameters:
            p.set(params[i])
            i += 1
        return y - function(self,x,parameters)

    if x == None: x = np.arange(y.shape[0])
    p = [param() for param in parameters]
    optimize.leastsq(f, p)
