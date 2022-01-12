#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 18 2021

@author: Mason Rogers

kv_spec2hist.py takes coefficient-space output from Dedalus in HDF5 format
and generates gridded output analogous to MITgcm output.
"""

#tinker
in_dir = "kv_snaps/"
in_file = in_dir+"kv_snaps_s1.h5"
out_dir = "../output/"
h_prefix = out_dir+"fine2d"

#imports
import h5py
import numpy as np
import xarray as xr
import xmitgcm as xm
from kv_param import *
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
time = np.arange(0, tStop+wFreq/2, wFreq)
p_c = p_c.interp(t=time)

#make grid
nx = round(2*R/dx) + 4
ny = round(2*R/dy) + 4
XC = np.arange((1-nx)/2*dx, nx/2*dx, dx)
YC = np.arange((1-ny)/2*dy, ny/2*dy, dy)
XC = xr.DataArray(data=XC, dims=['XC'], coords={'XC': XC})
YC = xr.DataArray(data=YC, dims=['YC'], coords={'YC': YC})

#normalized s (0, 2π) and r (-1, 1) for dedalus
s = np.arctan2(YC,XC)
rn = 2*np.sqrt(XC**2+YC**2)/R - 1
ks = p_c['ks']
Tr = p_c['Tr']

#pack data files (can get grid/meta from traj2hist)
for t_j in time:
    p_g = xr.zeros_like(rn) #initialize time step histogram
    for Tr_j in Tr:
        p_g += np.real(p_c.sel(t=t_j, ks=0, Tr=Tr_j) * eval_chebyt(Tr_j, rn)) #compute k=0 component
        for ks_j in ks[1:]: #compute k>0 components
            p_g += 2*np.real(p_c.sel(t=t_j, ks=ks_j, Tr=Tr_j) * eval_chebyt(Tr_j, rn) * np.exp(1j*ks_j*s))
    p_g = p_g.where(XC**2 + YC**2 < R**2, 0)

    xm.utils.write_to_binary(p_g.transpose('YC','XC').values.flatten(),
                             h_prefix+".{0:010d}".format(round(t_j/dt))+".data")

    print(t_j)
