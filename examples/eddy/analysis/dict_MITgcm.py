#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 2021

@author: Mason Rogers

read_MITgcm.py uses xmitgcm to read MITgcm output (or analogous) into xarray
gridded datasets.
"""

#imports
import xmitgcm as xm
import xarray as xr
import numpy as np
import xgcm as xg

#data locations
dirs = ['/pool001/masonr/eddyF/', '/pool001/masonr/eddyD/']
fnames = [['231205_fluid_4_diff'], ['231205_parti_4_diff']]
# dirs = ['/pool001/masonr/eddyD/', '/pool001/masonr/eddyRevisions/']
# fnames = [['231205_parti_4_diff'], ['B0.999996_d5e-03']]

#iters = [np.arange(0, 93311+1, 864)]*2
iters = [np.arange(0, 93311+1, 864)]*2

#load
ds = {}
gr = {}
for d, fs, i in zip(dirs, fnames, iters):
    nSoFar = 0
    tmp = {}
    for f in fs: 
        tmp[f] = xm.open_mdsdataset(d, iters=i, prefix=f, geometry='cartesian') 
        tmp[f]['time'] = tmp[f]['time'] / np.timedelta64(1,'s') / 1000
        tmp[f].coords['drCl'] = xr.DataArray(data=tmp[f]['drC'].values[:-1],
                                             coords={'Zl': tmp[f]['Zl'].values},
                                             dims='Zl')
        tmp[f] = tmp[f].drop_dims(['Zp1','Zu'])
        
        n_j = len(tmp[f].data_vars)
        for j in range(1, n_j+1):
            tmp[f] = tmp[f].rename(
                {'TRAC{0:02d}'.format(j): 'TRAC{0:02d}'.format(j+nSoFar)}
            )
        nSoFar += n_j

    ds[d] = xr.merge(tmp.values())
        
    gr[d] = xg.Grid(ds[d],
                    coords={'Z': {'center': 'Z', 'left':'Zl'},
                            'X': {'center':'XC', 'left':'XG'},
                            'Y': {'center':'YC', 'left':'YG'}})


