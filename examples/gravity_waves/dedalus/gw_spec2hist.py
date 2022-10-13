#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 18 2021

@author: Mason Rogers

kv_spec2hist.py takes coefficient-space output from Dedalus in HDF5 format
and generates gridded output analogous to MITgcm output.

This file is untested because other examples revealed that grid-space
interpolation is sufficiently accurate and much faster. Consider using grid2hist
instead.
"""

#tinker
in_dir = "gw_snaps/"
in_file = in_dir+"gw_snaps_s1.h5"
out_dir = "../output/"
h_prefix = out_dir+"2d"

#imports
import h5py
import numpy as np
import xarray as xr
import xmitgcm as xm
from ../input/gw_param import *
import matplotlib.pyplot as plt
from scipy.special import eval_chebyt
from matplotlib.colors import LogNorm

#read data from HDF5 file
with h5py.File(in_file, mode='r') as file:
    t = file['scales']['sim_time']
    Tζ = file['scales']['Tζ']
    kx = file['scales']['kx']
    p_c = xr.DataArray(data = file['tasks']['p_c'],
                       dims = ['t', 'kx', 'Tz'],
                       coords = {'t': t, 'ks': kx, 'Tr': Tz},
                       name = 'TRAC01')

#round to correct times
time = np.arange(0, tStop+wFreq/2, wFreq)
p_c = p_c.interp(t=time)

#make grid
nx = round(Lx/dx)
nz = round(Lz/dz) + 2
XC = np.arange(dx/2, nx*dx, dx)
Z = np.arange(dz/2, -nz*dz, -dz)
XC = xr.DataArray(data=XC, dims=['XC'], coords={'XC': XC})
Z = xr.DataArray(data=Z, dims=['Z'], coords={'Z': Z})

#normalized x (0, 2π) and r (-1, 1) for dedalus
x = XC/Lx * 2*np.pi
ζn = 2*(Z-η0*np.cos(k*XC-ω*time))/Lz + 1
kx = p_c['kx']
Tz = p_c['Tz']

#pack data files (can get grid/meta from traj2hist)
for t_j in time:
    p_g = xr.zeros_like(zn) #initialize time step histogram
    for Tz_j in Tz:
        p_g += np.real(p_c.sel(t=t_j, kx=0, Tz=Tz_j) * eval_chebyt(Tz_j, zn)) #compute k=0 component
        for kx_j in kx[1:]: #compute k>0 components
            p_g += 2*np.real(p_c.sel(t=t_j, kx=kx_j, Tz=Tz_j) * eval_chebyt(Tz_j, zn) * np.exp(1j*kx_j*s))
    p_g = p_g.where(-Lz + η0*np.cos(k*XC-ω*t_j) < Z < η0*np.cos(k*XC-ω*t_j), 0)

    xm.utils.write_to_binary(p_g.transpose('YC','XC').values.flatten(),
                             h_prefix+".{0:010d}".format(round(t_j/dt))+".data")

    print(t_j)
