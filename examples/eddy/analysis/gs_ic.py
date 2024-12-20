#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gs_ic.py plots the initial conditions

Created on Thu Jul 13 2023

@author: Mason Rogers
"""
'''-----------------------------------------------------------------------------
----------INIT------------------------------------------------------------------
-----------------------------------------------------------------------------'''

#tinker
saveFigures = True
uFile = '../mitgcm/gulf_stream/gs_0607.nc'
gFile = '../mitgcm/gulf_stream/gs_grid.nc'
pFile = '/pool001/masonr/eddyF/p_230712_blob.bin'
AFile = '/pool001/masonr/eddyF/RAC.data'
dzFile = '/pool001/masonr/eddyF/DRF.data'

#imports
import xgcm as xg
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba, ListedColormap

#style
def gen_op_cmap(c, α0=0):
    if isinstance(c, str):
        rgb = to_rgba(c)[:-1]
    else:
        rgb = c[0:3]
    arr = np.hstack([np.outer(np.ones(256), rgb), np.outer(np.linspace(α0,1,256), np.ones(1))])
    return ListedColormap(arr)
cTopo = gen_op_cmap('grey')

'''-----------------------------------------------------------------------------
----------CODE------------------------------------------------------------------
-----------------------------------------------------------------------------'''
#read data
ds = xr.open_dataset(uFile)
me = xr.open_dataset(gFile)
gr = xg.Grid(me, periodic=None,
             coords={'Z': {'center': 'k', 'left':'k_l'},
                     'X': {'center':'i', 'left':'i_g'},
                     'Y': {'center':'j', 'left':'j_g'}},
             metrics={('X',): ['DXC', 'DXG', 'DXV'],
                      ('Y',): ['DYC', 'DYG', 'DYU']})
p0 = np.fromfile(pFile, '>f4').reshape(87, 384, 216)
A = np.fromfile(AFile, '>f4').reshape(384, 216)
dz = np.fromfile(dzFile, '>f4').reshape(87)

#take derivatives
u = ds['U']
v = ds['V']
v_x = gr.derivative(v, 'X', boundary='fill')
u_y = gr.derivative(u, 'Y', boundary='fill')
ζ = (v_x - u_y).isel(k=0, i_g=slice(1, -1), j_g=slice(1,-1))
#ζ0 = (v_x - u_y).isel(i_g=171, j_g=242).values

#integrate initial release
p0_i = (p0*A).sum(axis=(-1,-2))
p0_i = p0_i / (p0_i*dz).sum()

#initial distribution width
p0_iz = (p0*dz.reshape((87,1,1))).sum(axis=0)
p0_iz = p0_iz / p0_iz.max()
p0_iz = p0_iz[1:-1,1:-1]

#initial distribution thickness
p0_sy = p0[:, 243, :]
p0_sy = p0_sy / p0_sy.max()

#for plots
x = me['XG'].isel(i_g=slice(1,-1), j_g=slice(1,-1))
y = me['YG'].isel(i_g=slice(1,-1), j_g=slice(1,-1))
z = me['Z']
top = (me['Depth'] == 0).isel(i=slice(1,-1), j=slice(1,-1))




'''-----------------------------------------------------------------------------
----------PLOT------------------------------------------------------------------
-----------------------------------------------------------------------------'''
def initialize_plots():
    #declare variables
    global f, a, p

    #declare plots
    f, a = plt.subplots(figsize=(10, 3), nrows=1, ncols=4,
                        layout='constrained',
                        gridspec_kw={'width_ratios': [1, .1, .3, 1]})

    #hide dummy plot for spacing
    a[1].set_visible(False)

    #label axes
    # f.suptitle('Initial Release', fontsize='large')
    # a[0].set_title('vorticity')
    a[0].set_xlabel('longitude [deg]', fontsize='small')
    a[0].set_ylabel('latitude [deg]', fontsize='small')
    # a[3].set_title('initial release')
    a[2].set_xlabel(r'$p_0(z) \ [{\sf m^{-1}}]$', fontsize='small')
    a[2].set_ylabel('depth [m]', fontsize='small')
    a[3].set_xlabel('longitude [deg]', fontsize='small')

    # #column labels
    # a[0].set_title('(a)', fontsize='small', loc='left')
    # a[2].set_title('(b)', fontsize='small', loc='left')
    # a[3].set_title('(c)', fontsize='small', loc='left')

    #tick formatting
    for aa in a: aa.tick_params(labelsize='x-small', which='both')

    #spines
    a[2].spines.right.set_visible(False)
    a[2].spines.top.set_visible(False)
        
    #prepare to store plots for legends
    p = [None, None]

def tidy_up_plots():
    #limits
    a[2].set_ylim([-500,0])
    a[3].set_ylim([-500,0])

    #legends
    #a[0].legend(p[1], ['release'])
    #a[3].legend(p[2:], [r'$p_0$', r'$\zeta_3$'], loc=4)

    #tick labels
    a[3].set_yticks(a[2].get_yticks())
    a[3].set_yticklabels([None]*len(a[2].get_yticks()))

    #colorbars
    plt.colorbar(p[0], ax=a[0], label=r'$\zeta_3 \ [{\sf 10^{-4} \ s^{-1}}]$')
    plt.colorbar(p[1], ax=a[3], label=r'$v \ [{\sf m~s^{-1}}]$')

    #save
    if saveFigures:
        today = np.datetime64('today').item()
        todayStr = '{0:02d}{1:02d}'.format(today.month, today.day)
        plt.figure(f.number) 
        plt.savefig('../figures/gsIC.png', dpi=200, transparent=False)

if __name__ == "__main__":
    try:
        plt.style.use('mason')
        initialize_plots()

        #transpose everything
        x = x.transpose('i_g', ...)
        y = y.transpose('i_g', ...)
        ζ = ζ.transpose('i_g', ...)
        p0_iz = np.transpose(p0_iz)
        top = top.transpose('i', ...)
        u = u.transpose('k', ...)
        v = v.transpose('k', ...)

        #colormaps
        vMax = np.max(np.abs(v.isel(j_g=243)))

        #vorticity plot
        p[0] = a[0].pcolormesh(x, y, ζ*1e4,
                               cmap='RdBu_r', vmin=-1, vmax=1)
        a[0].pcolormesh(x, y, top, cmap=cTopo)
        a[0].contour(x, y, p0_iz,
                     levels=np.array([np.exp(-4)]), colors='black', linewidths=1)

        #depth distribution/velocity plot
        p[1] = a[3].pcolormesh(x.isel(j_g=242), z, v.isel(i=slice(1,-1), j_g=243),
                               cmap='RdBu_r', vmin=-vMax, vmax=vMax)
        a[3].pcolormesh(x.isel(j_g=242), z,
                        z*xr.ones_like(me['i'][1:-1]) < -me['Depth'].isel(i=slice(1,-1), j=243),
                        cmap=cTopo)
        a[3].contour(x.isel(j_g=242), z, p0_sy[:, 1:-1],
                     levels=np.array([np.exp(-4)]), colors='black', linewidths=1)
        a[2].plot(p0_i, me['Z'], color='black')

        tidy_up_plots() 
        plt.show()
    finally:
        plt.close('all')

