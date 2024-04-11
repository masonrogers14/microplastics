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
#data_dirs = ['../output/', '../output/']
#file_names = [['2dspec'], ['4dtraj']] #[['small4d'], ['small2d']] #[['noad4d'], ['d32d']]
#data_dirs = ['../big_output/']; file_names = [['wayoffcen']]
#data_dirs = ['../output/']; file_names = [['trapped']]
data_dirs = ['../big_output/', '../output/']; file_names = [['offcen'], ['trapped']]

data_dirs = ['../mitgcm/run/']#, '../output/exp2/']
file_names = [['mitgcm2d_quad']]#, ['dedalus2d']]
#load
tmp = {}
nSoFar = 0
for dir, fname in zip(data_dirs, file_names):
    print(dir) 
    iters = np.hstack([np.arange(0,30001,1000)])#, np.arange(500000, 10500000, 10000)])
    tmp[nSoFar] = xm.open_mdsdataset(dir, iters=iters, prefix=fname, geometry='cartesian') 
    nTracs = 1 #len(tmp[nSoFar].data_vars)
    for j in range(1, nTracs+1):
        tmp[nSoFar] = tmp[nSoFar].rename({'TRAC{0:02d}'.format(j): 'TRAC{0:02d}'.format(j+nSoFar)})
    nSoFar += nTracs

#merge
ds = xr.combine_by_coords(tmp.values(), compat='override', combine_attrs='drop')

#edit grid
ds['time'] = ds['time'] / np.timedelta64(1,'s') / 1000
ds.coords['drCl'] = xr.DataArray(data=ds['drC'].values[:-1],
                                 coords={'Zl': ds['Zl'].values},
                                 dims='Zl')
ds = ds.drop_dims(['Zp1','Zu'])
gr = xg.Grid(ds,
             coords={'Z': {'center': 'Z', 'left':'Zl'},
                     'X': {'center':'XC', 'left':'XG'},
                     'Y': {'center':'YC', 'left':'YG'}})
