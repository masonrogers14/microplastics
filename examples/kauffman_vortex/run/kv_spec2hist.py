#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 12 2021

@author: Mason Rogers

kv_io.py: handle output from the Kauffman vortex experiment.
"""

#tinker
in_dir = "kv_snaps/"
in_file = in_dir+"kv_snaps_s1.h5"
out_dir = "bin/kv/"
h_prefix = out_dir+"2dspec"

#imports
import h5py
import numpy as np
import xarray as xr
import xmitgcm as xm
import matplotlib.pyplot as plt
from scipy.special import eval_chebyt
from matplotlib.colors import LogNorm

#read data from HDF5 file
with h5py.File(in_file, mode='r') as file:
    t = file['scales']['sim_time']
    Tr = file['scales']['Tr']
    ks = file['scales']['ks']
    p_c = xr.DataArray(data = file['tasks']['p_c'],
                       dims = ['t', 'ks', 'Tr'],
                       coords = {'t': t, 'ks': ks, 'Tr': Tr},
                       name = 'TRAC01')

#round to correct times
time = np.arange(0,1001,20)
p_c = p_c.interp(t=time)

#Kaufmann vortex parameters
R = .5

#make grid
dx = .01
dy = .01
dt = 1e-3
nx = round(2*R/dx) + 4
ny = round(2*R/dy) + 4
XC = np.arange((1-nx)/2*dx, nx/2*dx, dx)
YC = np.arange((1-ny)/2*dy, ny/2*dy, dy)
XC = xr.DataArray(data=XC, dims=['XC'], coords={'XC': XC}) 
YC = xr.DataArray(data=YC, dims=['YC'], coords={'YC': YC})
s = np.arctan2(YC,XC)
rn = 2*np.sqrt(XC**2+YC**2)/R - 1
ks_p = p_c['ks'].sel(ks=slice(1,None))
Tr = p_c['Tr']

#pack data files (can get grid/meta from ../traj/ode_io.jl)
for t_j in time: 
    p_c_j = p_c.sel(t=t_j)

    p_cp_j = p_c_j.sel(ks=slice(1,None))
    p_gp_j = (p_cp_j * eval_chebyt(Tr,rn) * np.exp(1j*ks_p*s)).sum(dim=['Tr','ks'])
    p_g0_j = (p_c_j.sel(ks=0) * eval_chebyt(Tr,rn)).sum(dim='Tr')
    p_g_j = np.real(2*p_gp_j + p_g0_j)

    p_g_j = p_g_j.where(XC**2 + YC**2 < R**2, 0)
    
    xm.utils.write_to_binary(p_g_j.transpose('YC','XC').values.flatten(),
                             h_prefix+".{0:010d}".format(round(t_j/dt))+".data") 
