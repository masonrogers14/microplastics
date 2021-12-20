#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 2021

@author: Mason Rogers
"""

#imports
import xmitgcm as xm
import xarray as xr
import numpy as np
import xgcm as xg

#data locations
data_dirs = ["../traj/bin/kv/", "../spec/bin/kv/"] #["../run_alt/", "../traj/bin/"]
file_names = [['4dprob'], ['2dspec']] #[['tfast'], ['twide']]

#load
tmp = {}
nSoFar = 0
for dir, fname in zip(data_dirs, file_names):
    print(dir)
    nTracs = 0
    iters = np.arange(0,1000000,20000)
    tmp[dir] = xm.open_mdsdataset(dir,iters=iters,
                                  prefix=[*fname],#'t1bud','t2bud','t3bud','t4bud'],
                                  geometry='cartesian') 
    nTracs = len(tmp[dir].data_vars)
    for j in range(1, nTracs+1): 
        tmp[dir] = tmp[dir].rename({'TRAC{0:02d}'.format(j): 'TRAC{0:02d}'.format(j+nSoFar)})
    nSoFar += nTracs

#merge
ds = xr.combine_by_coords(tmp.values(), compat='override', combine_attrs='drop')

#edit grid
ds.coords['drCl'] = xr.DataArray(data=ds['drC'].values[:-1],
                                 coords={'Zl': ds['Zl'].values},
                                 dims='Zl')
ds = ds.drop_dims(['Zp1','Zu'])
gr = xg.Grid(ds,
             coords={'Z': {'center': 'Z', 'left':'Zl'},
                     'X': {'center':'XC', 'left':'XG'},
                     'Y': {'center':'YC', 'left':'YG'}})
