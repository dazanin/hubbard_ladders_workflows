# Copyright (C) 2015 Institute for Theoretical Physics, ETH Zurich
#               1994-2009 by Bela Bauer <bauerb@phys.ethz.ch>
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

import copy
import numpy as np

class ResultProperties:
    def __init__(self):
        self.props = {}

class DataSet(ResultProperties):
    """
    The DataSet class stores a set of data, usually in XY format, along with all the properties
    describing the data, such as input parameters to the simulation etc.
    
    Members are:
     * x, y - These contain the data and are expected to come as lists of Numpy arrays
              by many functions operating on DataSets. However, for user-supplied functions,
              other ways of representing data may be used.
     * props - This is a dictionary of properties describing the dataset.
    """
    def __init__(self,x=None,y=None,props=None):
        ResultProperties.__init__(self)
        if x == None:   self.x = np.array([])
        else:           self.x = x
        if y == None:   self.y = np.array([])
        else:           self.y = y
        if props != None:   self.props = props
    
    def __repr__(self):
        return "x=%s\ny=%s\nprops=%s" % (self.x, self.y, self.props)
        
class ResultFile(ResultProperties):
    def __init__(self,fn=None):
        ResultProperties.__init__(self)
        if fn != None:
            self.props['filename'] = fn


