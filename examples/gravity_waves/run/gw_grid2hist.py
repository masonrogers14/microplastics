#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 18 2021

@author: Mason Rogers

gw_grid2hist.py takes grid-space output from Dedalus in HDF5 format
and generates gridded output analogous to MITgcm output.
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
from gw_param import *
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

#read data from HDF5 file
with h5py.File(in_file, mode='r') as file:
    t = file['scales']['sim_time']
    x = file['scales']['x']['1.0']
    ζ = file['scales']['ζ']['1.0']
    p_g = xr.DataArray(data = file['tasks']['p'],
                       dims = ['t', 'x', 'ζ'],
                       coords = {'t': t, 'x': x, 'ζ': ζ},
                       name = 'TRAC01')

#round to correct times
time = np.arange(0, tStop+wFreq/2, wFreq)
p_g = p_g.interp(t=time)

#wrap at 2π
p_w = p_g.sel(x=0).expand_dims('x').assign_coords(x=np.array([Lx]))
p_g = xr.concat([p_g, p_w], dim='x')

#make grid
nx = round(Lx/dx)
nz = round(Lz/dz) + 2
XC = np.arange(dx/2, nx*dx, dx)
Z = np.arange(dz/2, (1-nz)*dz, -dz)
XC = xr.DataArray(data=XC, dims=['XC'], coords={'XC': XC})
Z = xr.DataArray(data=Z, dims=['Z'], coords={'Z': Z})

#wavefield parameters
k = nk*2*np.pi/Lx
ω = np.sqrt(g*k)

#pack data files (can get grid/meta from traj2hist)
for t_j in time:
    ζ = Z - η0*np.cos(k*XC-ω*t_j)
    p_g_j = p_g.sel(t=t_j).interp(x=XC, ζ=ζ)
    p_g_j = p_g_j.where((-Lz < ζ) & (ζ < 0), 0)

    xm.utils.write_to_binary(p_g_j.transpose('Z','XC').values.flatten(),
                             h_prefix+".{0:010d}".format(round(t_j/dt))+".data")

    print(t_j)
