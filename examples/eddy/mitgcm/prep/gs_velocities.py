#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 2021

@author: Mason Rogers

gs_uf_to_u.py is designed to read LLC4320 velocity and grid data in netcdf format
and create meta/data files for plastic particle advection.

we use the form:
u^P = u^F + εL/U * 2(1-B)/(1+2*B) * (u^F_t + ([u^F - gεL/Ue_3] . ∇)u^F) 
and lump the buoyant term directly into the vertical advective velocity
"""

#imports
import dask as da
import xgcm as xg
import numpy as np
import xarray as xr
from gs_param import *

#tinker
dates = np.arange('2012-06-07', '2012-07-31', dtype='datetime64[D]')
dfiles = ['../gulf_stream/gs_{0:02d}{1:02d}.nc'.format(d.item().month, d.item().day) for d in dates]
gfile = '../gulf_stream/gs_grid_fixed.nc'
outdir = '/pool001/masonr/eddyRevisions/off/'
prefix = outdir + 'B{0:.3f}_d{1:.3e}_'.format(B, d)
cutInd = {'i': [None, None], 'j': [None, None], 'k': [None, None]} 
writeOutput = True

#init
ε = (1+2*B)*d**2*Us / (36*ν*Ls)
print(ε)
D = -5681.69

def uf_to_u(u, v, w, u_t, v_t, w_t, gr):
    #derivatives of velocity field
    #x---
    u_x = gr.derivative(u, 'X', boundary='extend')
    v_x = gr.derivative(v, 'X', boundary='extend')
    w_x = gr.derivative(w, 'X', boundary='extend')
    #y---
    u_y = gr.derivative(u, 'Y', boundary='extend')
    v_y = gr.derivative(v, 'Y', boundary='extend')
    w_y = gr.derivative(w, 'Y', boundary='extend')
    #z---
    u_z = gr.derivative(u, 'Z', boundary='extend')
    v_z = gr.derivative(v, 'Z', boundary='extend')
    w_z = gr.derivative(w, 'Z', boundary='extend')
    
    #interpolations
    #x---
    u_xu = gr.interp(u_x, 'X', boundary='extend')
    v_ζ3 = gr.interp(v, 'X', boundary='extend')
    w_ζ2 = gr.interp(w, 'X', boundary='extend')
    #y---
    u_ζ3 = gr.interp(u, 'Y', boundary='extend')
    v_yv = gr.interp(v_y, 'Y', boundary='extend')
    w_ζ1 = gr.interp(w, 'Y', boundary='extend')
    #z---
    u_ζ2 = gr.interp(u, 'Z', boundary='extend')
    v_ζ1 = gr.interp(v, 'Z', boundary='extend')
    w_zw = gr.interp(w_z, 'Z', boundary='extend')
    
    #advective fluxes
    #x---
    Fux = u * u_xu
    Fvx = u_ζ3 * v_x
    Fwx = u_ζ2 * w_x
    #y---
    Fvy = v * v_yv
    Fuy = v_ζ3 * u_y
    Fwy = v_ζ1 * w_y
    #z---
    Fuz = (w_ζ2 - g*ε*Ls/Us) * u_z
    Fvz = (w_ζ1 - g*ε*Ls/Us) * v_z
    Fwz = (w - g*ε*Ls/Us) * w_zw
    
    #interpolations
    #x---
    Fwx_w = gr.interp(Fwx, 'X', boundary='extend')
    Fvx_v = gr.interp(Fvx, 'X', boundary='extend')
    #y---
    Fuy_u = gr.interp(Fuy, 'Y', boundary='extend')
    Fwy_w = gr.interp(Fwy, 'Y', boundary='extend')
    #z---
    Fuz_u = gr.interp(Fuz, 'Z', boundary='extend')
    Fvz_v = gr.interp(Fvz, 'Z', boundary='extend')
    
    #correction velocities
    uC = ε*Ls/Us * 2*(1-B)/(1+2*B) * (u_t + Fux + Fuy_u + Fuz_u)
    vC = ε*Ls/Us * 2*(1-B)/(1+2*B) * (v_t + Fvx_v + Fvy + Fvz_v)
    wC = ε*Ls/Us * 2*(1-B)/(1+2*B) * (w_t + Fwx_w + Fwy_w + Fwz + g)
    
    #particle velocities
    uP = (u + uC).where(~np.isnan(ds['U']), 0)
    vP = (v + vC).where(~np.isnan(ds['V']), 0)
    wP = (w + wC).where(~np.isnan(ds['W']), 0)
    wP.loc[{'k_l': 0}] = 0 * xr.ones_like(wP.isel(k_l=0)) #np.minimum(wP.isel(k_l=0), 0.) #correct for surface

    return uP, vP, wP

#write output
#print('printing {0:d} timesteps:'.format(nt))
for i in range(len(dates)): #nt 
    #read data
    ds = xr.open_dataset(dfiles[i])
    ds = xr.combine_by_coords([ds, xr.open_dataset(gfile)])
    #ds = ds.chunk(time=6, i=36, i_g=36, j=64, j_g=64)
    ds['time'] = dates[i]
    ds = ds.isel(i=slice(cutInd['i'][0], cutInd['i'][1]),
                 i_g=slice(cutInd['i'][0], cutInd['i'][1]),
                 j=slice(cutInd['j'][0], cutInd['j'][1]),
                 j_g=slice(cutInd['j'][0], cutInd['j'][1]),
                 k=slice(cutInd['k'][0], cutInd['k'][1]),
                 k_l=slice(cutInd['k'][0], cutInd['k'][1]))
    
    #generate xgcm grid with metrics
    gr = xg.Grid(ds, periodic=[],
                 coords={'X': {'left': 'i_g', 'center': 'i'},
                         'Y': {'left': 'j_g', 'center': 'j'},
                         'Z': {'left': 'k_l', 'center': 'k'}})
    ds['DRF'] = gr.diff(ds['Zl'], 'Z', boundary='fill', fill_value=D)
    ds['DRC'] = gr.diff(ds['Z'], 'Z', boundary='fill') 
    gr = xg.Grid(ds, periodic=[],
                 coords={'X': {'left': 'i_g', 'center': 'i'},
                         'Y': {'left': 'j_g', 'center': 'j'},
                         'Z': {'left': 'k_l', 'center': 'k'}},
                 metrics = {('X',): ['DXC', 'DXG', 'DXV'],
                            ('Y',): ['DYC', 'DYG', 'DYU'],
                            ('Z',): ['DRC', 'DRF']})

    #calculate velocities
    u = ds['U'].load()
    v = ds['V'].load()
    w = ds['W'].load()
    if i == 0:
        dsP = xr.open_dataset(dfiles[1])
        u_t = (dsP['U'].load() - u) / ((dates[1] - dates[0]) / np.timedelta64(1, 's'))
        v_t = (dsP['V'].load() - v) / ((dates[1] - dates[0]) / np.timedelta64(1, 's'))
        w_t = (dsP['W'].load() - w) / ((dates[1] - dates[0]) / np.timedelta64(1, 's'))
    elif i < (len(dates) - 1):
        dsP = xr.open_dataset(dfiles[i+1])
        dsM = xr.open_dataset(dfiles[i-1])
        u_t = (dsP['U'].load() - dsM['U'].load()) / ((dates[i+1] - dates[i-1]) / np.timedelta64(1, 's'))
        v_t = (dsP['V'].load() - dsM['V'].load()) / ((dates[i+1] - dates[i-1]) / np.timedelta64(1, 's'))
        w_t = (dsP['W'].load() - dsM['W'].load()) / ((dates[i+1] - dates[i-1]) / np.timedelta64(1, 's'))
    else:
        dsM = xr.open_dataset(dfiles[-2])
        u_t = (u - dsM['U'].load()) / ((dates[-1] - dates[-2]) / np.timedelta64(1, 's'))
        v_t = (v - dsM['V'].load()) / ((dates[-1] - dates[-2]) / np.timedelta64(1, 's'))
        w_t = (w - dsM['W'].load()) / ((dates[-1] - dates[-2]) / np.timedelta64(1, 's'))
    uP, vP, wP = uf_to_u(u, v, w, u_t, v_t, w_t, gr)
    uP = uP.fillna(0).transpose('k','j','i_g',...)
    vP = vP.fillna(0).transpose('k','j_g','i',...)
    wP = wP.fillna(0).transpose('k_l','j','i',...)

    #write output
    n = int(np.round((dates[i] - dates[0]) / np.timedelta64(1, 's') / dt))
    n += int(np.round((dates[1] - dates[0]) / np.timedelta64(1, 's') / dt))
    with open(prefix+'U.{0:010d}'.format(n), "w") as ud:
        uP.values.astype('>f4').tofile(ud)
    with open(prefix+'V.{0:010d}'.format(n), "w") as vd:
        vP.values.astype('>f4').tofile(vd)
    with open(prefix+'W.{0:010d}'.format(n), "w") as wd:
        wP.values.astype('>f4').tofile(wd)
    print('{:8d}'.format(i))
