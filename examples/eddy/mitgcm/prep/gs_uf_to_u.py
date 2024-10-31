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
prefix = outdir + 'B{0:.6f}_d{1:.0e}_'.format(B, d)
cutInd = {'i': [None, None], 'j': [None, None], 'k': [None, None]} 
writeOutput = True

#init
ε = (1+2*B)*d**2*Us / (36*ν*Ls)
print(ε)
D = -5681.69

#read data
ds = xr.open_mfdataset(dfiles, combine='nested', concat_dim='time')
ds = xr.combine_by_coords([ds, xr.open_dataset(gfile)])
ds = ds.chunk(time=6, i=36, i_g=36, j=64, j_g=64)
ds['time'] = dates
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

##generate masks
#Zu = ds['Zl'].shift(k_l=-1, fill_value=D).rename(k_l='k')
#ds['hFacC'] = np.maximum(0, np.minimum(1, (Zu+ds['Depth'])/ds['DRF'])) 
#ds['maskC'] = ds['hFacC'] > 0

#fill velocities at land points
u = ds['U']
v = ds['V']
w = ds['W']

#derivatives of velocity field
#t---
u_t = u.differentiate('time', datetime_unit='s')
v_t = v.differentiate('time', datetime_unit='s')
w_t = w.differentiate('time', datetime_unit='s')
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

#CFLs
dtMax = {}
dtMax['uF'] = np.abs(ds['DXC'] / u).min(dim=['i_g','j','k'])
dtMax['uC'] = np.abs(ds['DXC'] / uC).min(dim=['i_g','j','k'])
dtMax['uP'] = np.abs(ds['DXC'] / uP).min(dim=['i_g','j','k'])
dtMax['vF'] = np.abs(ds['DYC'] / v).min(dim=['i','j_g','k'])
dtMax['vC'] = np.abs(ds['DYC'] / vC).min(dim=['i','j_g','k'])
dtMax['vP'] = np.abs(ds['DYC'] / vP).min(dim=['i','j_g','k'])
dtMax['wF'] = np.abs(ds['DRC'] / w).min(dim=['i','j','k_l'])
dtMax['wC'] = np.abs(ds['DRC'] / wC).min(dim=['i','j','k_l'])
dtMax['wP'] = np.abs(ds['DRC'] / wP).min(dim=['i','j','k_l'])
for k, v in dtMax.items():
    print('{0}: {1:.4f}'.format(k, v.min().values))

#cut to size and transpose
uP = uP.fillna(0).transpose('k','j','i_g',...)
vP = vP.fillna(0).transpose('k','j_g','i',...)
wP = wP.fillna(0).transpose('k_l','j','i',...)

#write output
if writeOutput:
    nt = uP.time.size
    nx = uP.i_g.size
    ny = uP.j.size
    nz = uP.k.size
    print('printing {0:d} timesteps:'.format(nt))
    for i in range(nt):
        n = int(np.round((ds['time'][i] - ds['time'][0]) / np.timedelta64(1, 's') / dt))
        n += int(np.round((ds['time'][1] - ds['time'][0]) / np.timedelta64(1, 's') / dt))
        with open(prefix+'U.{0:010d}'.format(n), "w") as ud:
            uP.isel(time=i).values.astype('>f4').tofile(ud)
        print('U', end=', ')
        with open(prefix+'V.{0:010d}'.format(n), "w") as vd:
            vP.isel(time=i).values.astype('>f4').tofile(vd)
        print('V', end=', ')
        with open(prefix+'W.{0:010d}'.format(n), "w") as wd:
            wP.isel(time=i).values.astype('>f4').tofile(wd)
        print('W')
        #    with open(prefix+'U.{0:010d}.meta'.format(n), "w") as um:
        #        um.write(" nDims = [3];\n")
        #        um.write(" dimList = [\n {0:d}, 1, {0:d}, \
        #                              \n {1:d}, 1, {1:d}, \
        #                              \n {2:d}, 1, {2:d} \
        #                              \n ];\n".format(nx, ny, nz))
        #        um.write(" dataprec = ['float32'];\n")
        #        um.write(" nrecords = [1];\n")
        #        um.write(" timestepnumber = [1];\n")
        #    with open(prefix+'V.{0:010d}.meta'.format(n), "w") as vm:
        #        vm.write(" nDims = [3];\n")
        #        vm.write(" dimList = [\n {0:d}, 1, {0:d}, \
        #                              \n {1:d}, 1, {1:d}, \
        #                              \n {2:d}, 1, {2:d} \
        #                              \n ];\n".format(nx, ny, nz))
        #        vm.write(" dataprec = ['float32'];\n")
        #        vm.write(" nrecords = [1];\n")
        #        vm.write(" timestepnumber = [1];\n")
        #    with open(prefix+'W.{0:010d}.meta'.format(n), "w") as wm:
        #        wm.write(" nDims = [3];\n")
        #        wm.write(" dimList = [\n {0:d}, 1, {0:d}, \
        #                              \n {1:d}, 1, {1:d}, \
        #                              \n {2:d}, 1, {2:d} \
        #                              \n ];\n".format(nx, ny, nz))
        #        wm.write(" dataprec = ['float32'];\n")
        #        wm.write(" nrecords = [1];\n")
        #        wm.write(" timestepnumber = [1];\n") 
        #print('    {0:d}: '.format(i) + ' '.join(['{} {:.0f}'.format(k, v[i]) for k, v in dtMax.items()]))
        print('{:8d}'.format(i))
