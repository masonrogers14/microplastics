#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 08 2023

@author: Mason Rogers

read_MITgcm.py uses xmitgcm to read MITgcm output (or analogous) into xarray
gridded datasets.
"""

#imports
import xarray as xr
import numpy as np
import xgcm as xg
from mds import rdmds
from gs_param import *

#data locations
#dirs = ['/pool001/masonr/eddy/', '/pool001/masonr/eddy3/']#,'../mitgcm/run_c/']
#fnames = [['230713_fluid_4_diff'], ['230713_parti_4_diff']]#, ['fluid']]
# dirs = ['/pool001/masonr/eddy/', '/pool001/masonr/eddyD/']
# fnames = [['230713_fluid_4_diff'], ['231205_parti_4_diff']]
dirs = ['/pool001/masonr/eddyF/', '/pool001/masonr/eddyD/']
fnames = [['231205_fluid_4_diff'], ['231205_parti_4_diff']]
#iters = [np.arange(0,23000+1,144)]*2
iters = [np.arange(0, 93311+1, 864)]*2

#grid info
gm = {'RAC': ['XC', 'YC'], 'RAW': ['XG', 'YC'], 'RAS': ['XC', 'YG'], 'RAZ': ['XG', 'YG'], #keys: metrics
      'DXG': ['XC', 'YG'], 'DYG': ['XG', 'YC'], 'DXC': ['XG', 'YC'], 'DYC': ['XC', 'YG'], #values: dims
      'DRC': ['Zp1'], 'DRF': ['Z'], 'Depth': ['XC', 'YC'],
      'hFacC': ['XC', 'YC', 'Z'], 'hFacS': ['XC', 'YG', 'Z'], 'hFacW': ['XG', 'YC', 'Z']}

#load data
ds = {}
gr = {}
for d, fs, it in zip(dirs, fnames, iters):
    nSoFar = 0
    tmp = {}
    ds[d] = xr.Dataset()
    tmp['XC'] = xr.DataArray(rdmds(d+'XC')[0,:], dims='XC'); #tmp['XC'] = tmp['XC'].assign_coords(XC=tmp['XC'])
    tmp['YC'] = xr.DataArray(rdmds(d+'YC')[:,0], dims='YC'); #tmp['YC'] = tmp['YC'].assign_coords(YC=tmp['YC'])
    tmp['XG'] = xr.DataArray(rdmds(d+'XG')[0,:], dims='XG'); #tmp['XG'] = tmp['XG'].assign_coords(XG=tmp['XG'])
    tmp['YG'] = xr.DataArray(rdmds(d+'YG')[:,0], dims='YG'); #tmp['YG'] = tmp['YG'].assign_coords(YG=tmp['YG'])
    tmp['Z'] = xr.DataArray(rdmds(d+'RC')[:,0,0], dims='Z'); #tmp['Z'] = tmp['Z'].assign_coords(Z=tmp['Z'])
    tmp['Zp1'] = xr.DataArray(rdmds(d+'RF')[:,0,0], dims='Zp1'); #tmp['Zp1'] = tmp['Zp1'].assign_coords(Zp1=tmp['Zp1'])
    tmp['Zu'] = tmp['Zp1'][1:].rename(Zp1='Zu')
    tmp['Zl'] = tmp['Zp1'][:-1].rename(Zp1='Zl')
    for k, v in gm.items():
        ds[d][k] = xr.DataArray(np.squeeze(rdmds(d+k)), dims=v[::-1], coords={v_j: tmp[v_j] for v_j in v[::-1]}, name=k)
    for f in fs[:1]: 
        data = rdmds(d+f, list(it))
        if it.size == 1: 
            ds[d]['TRAC01'] = xr.DataArray(data, dims=['Z','YC','XC'], coords={'Z': tmp['Z'], 'YC': tmp['YC'], 'XC': tmp['XC']}, name=f)
        else:
            ds[d]['TRAC01'] = xr.DataArray(data, dims=['time', 'Z','YC','XC'], coords={'time': it*dt, 'Z': tmp['Z'], 'YC': tmp['YC'], 'XC': tmp['XC']}, name=f)
        #n_j = len(tmp[f].data_vars)
        #for j in range(1, n_j+1):
        #    tmp[f] = tmp[f].rename({'TRAC{0:02d}'.format(j): 'TRAC{0:02d}'.format(j+nSoFar)})
        #nSoFar += n_j
        
    gr[d] = xg.Grid(ds[d], periodic=None,
                    coords={'Z': {'center': 'Z', 'outer':'Zp1'},
                            'X': {'center':'XC', 'left':'XG'},
                            'Y': {'center':'YC', 'left':'YG'}},
                    metrics={('X',): ['DXC', 'DXG'], ('Y',): ['DYC', 'DYG'], ('Z',): ['DRC', 'DRF'],
                             ('X', 'Y'): ['RAC', 'RAW', 'RAS', 'RAZ']})

