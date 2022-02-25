#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 18 2021

@author: Mason Rogers

kv_spec2hist.py takes grid-space output from Dedalus in HDF5 format
and generates gridded output analogous to MITgcm output.
"""

#tinker
in_dir = "kv_snaps/"
in_file = in_dir+"kv_snaps_s1.h5"
out_dir = "../output/"
h_prefix = out_dir+"d32d"

#imports
import h5py
import numpy as np
import xarray as xr
import xmitgcm as xm
from ../input/kv_param import *
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

#read data from HDF5 file
with h5py.File(in_file, mode='r') as file:
    t = np.array(file['scales']['sim_time'])
    r = np.array(file['tasks']['p'].dims[2][0][:])
    s = np.array(file['tasks']['p'].dims[1][0][:])
    p_g = xr.DataArray(data = np.array(file['tasks']['p'][:,:,:]),
                       dims = ['t', 's', 'r'],
                       coords = {'t': t, 's': s, 'r': r},
                       name = 'TRAC01')

#round to correct times
time = np.arange(0, tStop+wFreq/2, wFreq)
p_g = p_g.interp(t=time)

#wrap at 2π
p_w = p_g.sel(s=0).expand_dims('s').assign_coords(s=np.array([2*np.pi]))
p_g = xr.concat([p_g, p_w], dim='s')

#make grid
nx = round(2*R/dx) + 4
ny = round(2*R/dy) + 4
XC = np.arange((1-nx)/2*dx, nx/2*dx, dx)
YC = np.arange((1-ny)/2*dy, ny/2*dy, dy)
XC = xr.DataArray(data=XC, dims=['XC'], coords={'XC': XC})
YC = xr.DataArray(data=YC, dims=['YC'], coords={'YC': YC})

#interpolated polar coordinates
s = np.mod(np.arctan2(YC,XC), 2*np.pi)
r = np.sqrt(XC**2+YC**2)

#pack data files (can get grid/meta from traj2hist)
for t_j in time:
    p_g_j = p_g.sel(t=t_j).interp(s=s,r=r)
    p_g_j = p_g_j.where(XC**2 + YC**2 < R**2, 0)

    xm.utils.write_to_binary(p_g_j.transpose('YC','XC').values.flatten(),
                             h_prefix+".{0:010d}".format(round(t_j/dt))+".data")

    print(t_j)
