#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 2021

@author: Mason Rogers
"""

#tinker
save_figures = True
nTracs = 3
out_fname = "../figures/1216_contract.png"
θ_fname = "eps00/att_theta.bin"
z_fname = "eps00/att_z.bin"
r_fname = "eps00/att_r.bin"

#imports
import numpy as np
import xarray as xr
import dask.array as da
import matplotlib.pyplot as plt
from read_MITgcm import ds
from scipy.interpolate import CubicSpline

#extract dataset
ds = list(ds.values())[0]
XC = ds['XC']
YC = ds['YC']
Z = ds['Z']

#obtain attractor as callable
θ_att = np.fromfile(θ_fname)
z_att = np.fromfile(z_fname)
r_att = np.fromfile(r_fname)
θ_att = np.append(θ_att, θ_att[0]+2*np.pi)
z_att = np.append(z_att, z_att[0])
r_att = np.append(r_att, r_att[0])
z_of_θ = CubicSpline(θ_att, z_att, bc_type='periodic')
r_of_θ = CubicSpline(θ_att, r_att, bc_type='periodic')

#compute distribution spread around attractor
d_rms = [None]*nTracs
for i in range(1, nTracs+1):
    p = ds['TRAC0'+str(i)]
    XC = xr.DataArray(data=XC, dims=['XC'])
    YC = xr.DataArray(data=YC, dims=['YC'])
    Z = xr.DataArray(data=Z, dims=['Z'])
    θ_obs = da.arctan2(YC, XC)
    z_obs = Z
    r_obs = da.sqrt(da.square(XC) + da.square(YC))
    z_exp = xr.DataArray(data=z_of_θ(θ_obs), dims=θ_obs.dims)
    r_exp = xr.DataArray(data=r_of_θ(θ_obs), dims=θ_obs.dims)
    z_err = z_obs - z_exp
    r_err = r_obs - r_exp
    d_rms[i-1] = da.sqrt(((z_err**2+r_err**2)*p).sum(dim=['XC', 'YC', 'Z']) / p.sum(dim=['XC', 'YC', 'Z'])).compute()



'''-----------------------------------------------------------------------------
--------PLOTTING FUNCTIONS------------------------------------------------------
-----------------------------------------------------------------------------'''
#tinker
bfs = 14
blw = 1.5

def initialize_plots():
    #declare variables
    global f_dis, a_dis, p_dis

    #declare plots
    f_dis, a_dis = plt.subplots()

    #label axes
    a_dis.set_xlabel(r"$t$")
    a_dis.set_ylabel(r"$d_{RMS}$")

    #titles
    a_dis.set_title("RMS distance from attractor over time")

    #prepare to store plots for legends
    p_dis = None

def tidy_up_plots():
    #legends
    plt.legend(*a_dis.get_legend_handles_labels(), title=r"$\kappa$")

    #save stuff
    if save_figures:
        plt.figure(f_dis.number)
        plt.savefig(out_fname)

if __name__ == "__main__":
    try:
        initialize_plots()
        for i in range(nTracs):
            a_dis.plot(d_rms[i]['time'], d_rms[i], label="{0:.1e}".format(i*1e-6))
        tidy_up_plots()

        plt.show()
    finally:
        plt.close('all')
