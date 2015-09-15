# Copyright (C) 2015 Institute for Theoretical Physics, ETH Zurich
#               1994-2009 by Bela Bauer <bauerb@phys.ethz.ch>
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

import os.path
import sys
import glob
import math
import copy
import numpy as np

from hlist import deep_flatten, flatten, depth
from dict_intersect import dict_intersect
from dataset import DataSet

def make_list(infiles):
    if type(infiles) == list:
      return infiles
    else:
      return [infiles]

def size(lst):
    try:
      return len(lst)
    except:
      return 1

def recursiveGlob(dirname,pattern):
    ret = glob.glob(os.path.join(dirname, pattern))
    for d in os.listdir(dirname):
        d = os.path.join(dirname, d)
        if os.path.isdir(d):
            ret += recursiveGlob(d, pattern)
    return ret

def collectXY(sets,x,y,foreach=[],ignoreProperties=False):
      """ collects specified data from a list of DataSet objects
         
          this function is used to collect data from a list of DataSet objects, to prepare plots or evaluation. The parameters are:
    
            sets:    the list of datasets
            x:       the name of the property or measurement to be used as x-value of the collected results 
            y:       the name of the property or measurement to be used as y-value of the collected results 
            foreach: an optional list of properties used for grouping the results. A separate DataSet object is created for each unique set of values of the specified parameers.
            ignoreProperties: setting ignoreProperties=True prevents collectXY() from collecting properties.
            
          The function returns a list of DataSet objects.
      """
      foreach_sets = {}
      for iset in flatten(sets):
          if iset.props['observable'] != y and not y in iset.props:
              continue
          
          fe_par_set = tuple((iset.props[m] for m in foreach))
          if fe_par_set in foreach_sets:
              foreach_sets[fe_par_set].append(iset)
          else:
              foreach_sets[fe_par_set] = [iset]
      for k,v in foreach_sets.items():
          common_props = dict_intersect([q.props for q in v])
          res = DataSet()
          res.props = common_props
          for im in range(0,len(foreach)):
              m = foreach[im]
              res.props[m] = k[im]
          res.props['xlabel'] = x
          res.props['ylabel'] = y
          
          for data in v:
              if data.props['observable'] == y:
                  if len(data.y)>1:
                      res.props['line'] = '.'
                  xvalue = np.array([data.props[x] for i in range(len(data.y))])
                  if len(res.x) > 0 and len(res.y) > 0:
                      res.x = np.concatenate((res.x, xvalue ))
                      res.y = np.concatenate((res.y, data.y))
                  else:
                      res.x = xvalue
                      res.y = data.y
              elif not ignoreProperties:
                  res.props['line'] = '.'
                  xvalue = np.array([ data.props[x] ])
                  if len(res.x) > 0 and len(res.y) > 0:
                      res.x = np.concatenate((res.x, xvalue ))
                      res.y = np.concatenate((res.y, np.array([ data.props[y] ])))
                  else:
                      res.x = xvalue
                      res.y = np.array([ data.props[y] ])
          
          order = np.argsort(res.x, kind = 'mergesort')
          res.x = res.x[order]
          res.y = res.y[order]
          res.props['label'] = ''
          for im in range(0,len(foreach)):
              res.props['label'] += '%s = %s ' % (foreach[im], k[im])
          
          foreach_sets[k] = res
      return foreach_sets.values()

def ResultsToXY(sets,x,y,foreach=[]):
    """ combines observable x and y to build a list of DataSet with y vs x
 
    this function is used to collect data from a hierarchy of DataSet objects, to prepare plots or evaluation.
    the inner-most list has to contain one DataSet with props['observable'] = x and one props['observable'] = y,
    this will be the pair x-y used in the collection.

    The parameters are:
      sets:    hierarchy of datasets where the inner-most list must contain to pair x-y
      x:       the name of the observable to be used as x-value of the collected results 
      y:       the name of the observable to be used as y-value of the collected results 
      foreach: an optional list of properties used for grouping the results. A separate DataSet object is created for each unique set of values of the specified parameers.

    The function returns a list of DataSet objects.
    """
    
    dd = depth(sets)
    if dd < 2:
        raise Exception('The input hierarchy does not provide a unique pair x-y. The input structure has to be a list of lists as minimum. pyalps.groupSets might help you.')
    
    hgroups = flatten(sets, fdepth=-1)
    
    foreach_sets = {}
    for gg in hgroups:
        xset = None
        yset = None
        for d in gg:
            if d.props['observable'] == x:
                xset = d
            if d.props['observable'] == y:
                yset = d
        if xset is None or yset is None:
            continue
        
        common_props = dict_intersect([d.props for d in gg])
        fe_par_set = tuple((common_props[m] for m in foreach))
        
        if not fe_par_set in foreach_sets:
            foreach_sets[fe_par_set] = DataSet()
            foreach_sets[fe_par_set].props = common_props
            foreach_sets[fe_par_set].props['xlabel'] = x
            foreach_sets[fe_par_set].props['ylabel'] = y
        
        if len(xset.y) == len(yset.y):
            foreach_sets[fe_par_set].x = np.concatenate((foreach_sets[fe_par_set].x, xset.y))
            foreach_sets[fe_par_set].y = np.concatenate((foreach_sets[fe_par_set].y, yset.y))
        elif len(xset.y) == 1:
            foreach_sets[fe_par_set].x = np.concatenate((foreach_sets[fe_par_set].x, np.array( [xset.y[0]]*len(yset.y) )))
            foreach_sets[fe_par_set].y = np.concatenate((foreach_sets[fe_par_set].y, yset.y))
    
    for k, res in foreach_sets.items():
        order = np.argsort(res.x, kind = 'mergesort')
        res.x = res.x[order]
        res.y = res.y[order]
        res.props['label'] = ''
        for p in foreach:
            res.props['label'] += '%s = %s ' % (p, res.props[p])
        
    return foreach_sets.values()

def groupSets(groups, for_each = []):
    """ groups a list of DataSet objects into a list of lists
        
        this function groups a list of DataSet objects into a list of lists, according to the values of the properties given in the for_ech argument. DataSet objects with the same values of the properties given in for_each are grouped together.
        The parameters are:
          data: the data to be grouped
          for_each: the properties according to which the data is grouped
    """
    dd = depth(groups)

    if dd > 1:
        hgroups = flatten(groups, -1)
        hgroups_idcs = hgroups.indices()
    else:
        hgroups = [groups]
        hgroups_idcs = [0]

    for idx in hgroups_idcs:
        sets = hgroups[idx]

        for_each_sets = {}
        for iset in sets:
            fe_par_set = tuple((iset.props[m] for m in for_each))

            if fe_par_set in for_each_sets:
                for_each_sets[fe_par_set].append(iset)
            else:
                for_each_sets[fe_par_set] = [iset]

        hgroups[idx] = for_each_sets.values()

    if dd > 1:
        return groups
    else:
        return hgroups[0]

def select(inp,condition):
    data_ = []
    for ds in flatten(inp):
        if condition(ds):
            data_.append(ds)
    return data_

def select_by_property(data,proplist):
    for k,v in proplist.items():
        data = select(data, lambda ds: ds.props[k]==v)
    return data

def values(data, key):
    vals = []
    if type(key) == list:
        for ds in flatten(data):
            keyv = tuple([ds.props[kk] for kk in key])
            if keyv not in vals:
                vals.append(keyv)
        return vals
    else:
        for ds in flatten(data):
            if ds.props[key] not in vals:
                vals.append(ds.props[key])
        return np.sort(vals)

def mergeDataSets(dsets):
    props = dict_intersect([d.props for d in dsets])
    merged = copy.deepcopy(dsets.pop())
    for d in dsets:
        if (d.x != merged.x).any():   raise ValueError('cannot merge datasets: x values mismatch')
        if len(d.y) != len(merged.y): raise ValueError('cannot merge datasets: y shape mismatch')
        if isinstance(merged.y,alea.MCVectorData):
            merged.y.merge(d.y)
        else:
            for i in range(len(merged.y)):
                merged.y[i].merge(d.y[i])
    merged.props = props
    return merged

def mergeMeasurements(measurements):
    byname = {}
    for mset in measurements:
        for m in mset:
            key = m.props['observable']
            if key not in byname:   byname[key] = [m]
            else:                   byname[key].append(m)
    merged = [mergeDataSets(v) for v in byname.values()]
    return merged

def SetLabels (data, proplist):
    """
    Set labels according to the properties given in 'proplist'.
    """
    if type(proplist) == str:
        proplist = [proplist]
    for d in flatten(data):
        if 'label' in d.props:
            del d.props['label']
        label_items = []
        for propname in proplist:
            try:
                label_items.append(propname+' = %s' % (d.props[propname]))
            except:
                pass
        d.props['label'] = ', '.join(label_items)
    return data

def CycleColors (data, foreach,
                 colors=['k','b','g','m','c','y']):
    """
    Cyclically assign colors to the lines/markers that will be used to
    display the DataSets, based on the properties in 'foreach'. This means
    that DataSet instances that have the same values for the properties
    in 'foreach' will receive the same color.
    """
    all = {}
    for q in flatten(data):
        key = tuple([q.props[k] for k in foreach])
        all[key] = ''

    icolor = 0
    for k in all.keys():
        all[k] = colors[icolor]
        icolor = (icolor+1)%len(colors)

    for q in flatten(data):
        key = tuple([q.props[k] for k in foreach])
        q.props['color'] = all[key]

    return data

def CycleMarkers (data, foreach,
                  markers=['s', 'o', '^', '>', 'v', '<', 'd', 'p', 'h', '+', 'x']):
    """
    Cyclically assign markers to the lines/markers that will be used to
    display the DataSets, based on the properties in 'foreach'. This means
    that DataSet instances that have the same values for the properties
    in 'foreach' will receive the same marker.
    """
    all = {}
    for q in flatten(data):
        key = tuple([q.props[k] for k in foreach])
        all[key] = ''

    imarker = 0
    for k in all.keys():
        all[k] = markers[imarker]
        imarker = (imarker+1)%len(markers)

    for q in flatten(data):
        key = tuple([q.props[k] for k in foreach])
        q.props['marker'] = all[key]
        if 'line' in q.props:
            ll = list(q.props['line'])
            ll[0] = all[key]
            q.props['line'] = ''.join(ll)
        else:
            q.props['line'] = all[key] + '-'

    return data

