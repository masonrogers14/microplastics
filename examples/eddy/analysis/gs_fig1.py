#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
kv_fig2.py plots the variance growth for two configurations of the Kauffman vortex
experiment.

Created on Fri Jul 22 2022

@author: Mason Rogers
"""
'''-----------------------------------------------------------------------------
----------INIT------------------------------------------------------------------
-----------------------------------------------------------------------------'''

#imports
from matplotlib.colors import LogNorm, to_rgba, ListedColormap

#tinker
saveFigures = True
pFiles = ['p_small.py', 'p_large.py']

#imports
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from read_MITgcm import ds, gr, dirs

#relevant variables for plot parameters
nConfs = len(ds.keys())
zSnaps = np.array([0, 15, 27])
tSnap = -1
xSlice = -1 # < 0 --> integrate
ySlice = -1
zSlice = -1
x0 = -76.489586  
y0 = 32.71517

#init
nSnaps = zSnaps.size

def gen_op_cmap(c, α0=0):
    if isinstance(c, str):
        rgb = to_rgba(c)[:-1]
    else:
        rgb = c[0:3]
    arr = np.hstack([np.outer(np.ones(256), rgb), np.outer(np.linspace(α0,1,256), np.ones(1))])
    return ListedColormap(arr)
def gen_white_cmap(c, α0=0):
    if isinstance(c, str):
        rgb = to_rgba(c)[:-1]
    else:
        rgb = c[0:3]
    rgb = np.array(rgb)
    arr = np.hstack([np.linspace((1-α0)*(1-rgb)+rgb, rgb, 256), np.ones((256,1))])
    return ListedColormap(arr)
cmaps = gen_white_cmap(plt.cm.get_cmap('tab10')(0))
cTopo = gen_op_cmap('black')



'''-----------------------------------------------------------------------------
----------CODE------------------------------------------------------------------
-----------------------------------------------------------------------------'''
p = {}
pMax = 0
for k in ds.keys():
    pTot = (ds[k]['TRAC01']*ds[k]['DRF']*ds[k]['RAC']).isel(time=0).sum()
    p[k] = (ds[k]['TRAC01']*ds[k]['DRF']).isel(time=tSnap, Z=zSnaps) / pTot
    pMax = np.maximum(pMax, p[k].max())


'''-----------------------------------------------------------------------------
----------PLOT------------------------------------------------------------------
-----------------------------------------------------------------------------'''
#tinker
bfs = 14
blw = 2

def initialize_plots():
    #declare variables
    global fC, aC, pC, cC

    #declare plots
    fC = plt.figure(figsize=(5,5), layout='constrained')
    gC = fC.add_gridspec(nSnaps, nConfs+1, width_ratios=[9, 9, 1])
    aC = np.array([np.array([fC.add_subplot(gC[i,j])
                   for j in range(0, nConfs)])
                   for i in range(nSnaps)])
    cC = fC.add_subplot(gC[:,-1])

    #label axes
    for a in aC[-1]: a.set_xlabel(r'longitude', fontsize=bfs)
    for j in range(nSnaps): 
        aC[j,0].set_ylabel('$z = {0:.0f}$ m\nlatitude'.format(ds[k]['Z'][zSnaps[j]]),
                            fontsize=bfs)
    aC[0,0].set_title('fluid parcels', fontsize=bfs)
    aC[0,1].set_title('microplastics\n$(B=.99, d=.1 \ {\sfm m})$', fontsize=bfs)

    #tick formatting
    for a in aC.flatten():
        a.set_xticks(np.unique(np.round(ds[k]['XC'] + x0)))
        a.set_yticks(np.unique(np.round(ds[k]['YC'] + y0)))
        a.xaxis.set_major_formatter(lambda x, pos: '{0:.0f}'.format(x))
        a.yaxis.set_major_formatter(lambda x, pos: '{0:.0f}'.format(x))
        
    #prepare to store plots for legends
    pC = [[None for j in range(nConfs)] for i in range(nSnaps)]

def tidy_up_plots():
    #colorbars
    cbar = plt.colorbar(pC[0][0], cax=cC)
    cbar.set_label('$p(x, y) \ [{\sf 10^{-12} m^{-2}}]$', fontsize=bfs)

    #share axes
    for a in aC.flatten():
        a.set_xlim(aC[0,0].get_xlim())
        a.set_ylim(aC[0,0].get_ylim())
    for a in aC:
        for j in range(1,nConfs):
            a[j].set_yticklabels(['']*len(a[j].get_yticks()))
    for a in aC[:-1,:].flatten():
        a.set_xticklabels(['']*len(a.get_xticks()))
    

    #save
    if saveFigures:
        today = np.datetime64('today').item()
        todayStr = '{0:02d}{1:02d}'.format(today.month, today.day)
        plt.figure(fC.number) 
        plt.savefig('../figures/'+todayStr+'_fig1.png')

if __name__ == "__main__":
    try:
        initialize_plots()

        for r in range(nSnaps):
            for k, c in zip(ds.keys(), range(nConfs)):
                pC[r][c] = aC[r,c].pcolormesh(ds[k]['XC'] + x0,
                                              ds[k]['YC'] + y0,
                                              p[k].isel(Z=r) * 1e12,
                                              cmap=cmaps,
                                              vmin=0,
                                              vmax=pMax * 1e12)
                aC[r,c].pcolormesh(ds[k]['XC'] + x0,
                                   ds[k]['YC'] + y0,
                                   ds[k]['hFacC'].isel(Z=zSnaps[r]) == 0,
                                   cmap=cTopo)

        tidy_up_plots() 
        plt.show()
    finally:
        plt.close('all')

