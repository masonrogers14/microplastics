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
dirs = ['../output/exp1/', '../output/exp2/']
fnames = [['julia4dx', 'julia2dx', 'dedalus2d', 'mitgcm2d'], ['julia4dx', 'julia2dx', 'dedalus2d', 'mitgcm2d']]
iters = [np.arange(0,100000,1000), np.arange(0,100000,1000)]

#dirs = ['../mitgcm/run/']; fnames=[['mitgcm2d_quad']]; iters = [np.arange(0, 30001, 1000)]
#dirs = ['../output/for_movies/']; fnames = [['julia4d', 'mitgcm2d_movie']]; iters = [np.arange(0, 30001, 250)];

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
            tmp[f] = tmp[f].rename({'TRAC{0:02d}'.format(j): 'TRAC{0:02d}'.format(j+nSoFar)})
        nSoFar += n_j

    ds[d] = xr.combine_by_coords(tmp.values())
        
    gr[d] = xg.Grid(ds[d],
                    coords={'Z': {'center': 'Z', 'left':'Zl'},
                            'X': {'center':'XC', 'left':'XG'},
                            'Y': {'center':'YC', 'left':'YG'}})


#hasty garbage code
x = ds[dirs[0]]
#x['TRAC02'] = x['TRAC01']
#x['TRAC03'] = x['TRAC01']
#x['TRAC04'] = x['TRAC01']
