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

#tinker
saveFigures = True
conf = '../output/exp2/'
pFile = 'p_large.py'
tracs = [1, 4]
times = [0, 25, 50, 75]

#imports
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from dict_MITgcm import ds, gr
from matplotlib.colors import to_rgba, ListedColormap

#read files
with open(pFile, 'r') as f:
    exec(f.read())
    ϵ = ((1+2*B)*d**2*Us)/(36*ν*Ls) #small parameter

#initialize
labels = ['full Monte Carlo', 'reduced Monte Carlo',
          'reduced Dedalus', 'reduced MITgcm']
labels = [labels[j-1] for j in tracs]
nTracs = len(tracs)
nTimes = len(times)
p = [None] * nTracs
for j in range(nTracs):
    p[j] = ds[conf]['TRAC0{0:d}'.format(tracs[j])]
    if 'Z' in p[j].dims:
        p[j] = p[j].sum(dim='Z') * dz




'''-----------------------------------------------------------------------------
----------CODE------------------------------------------------------------------
-----------------------------------------------------------------------------'''
#interpolate to polar coordinates
θ = xr.DataArray(data=np.linspace(0,2*np.pi,endpoint=True),
                 dims=['θ'],
                 coords={'θ': np.linspace(0,2*np.pi,endpoint=True)})
r = xr.DataArray(data=np.linspace(0,2*R/3,endpoint=False),
                 dims=['r'],
                 coords={'r': np.linspace(0,2*R/3,endpoint=False)})
for j in range(nTracs):
    p[j] = p[j].interp(XC=r*np.cos(θ), YC=r*np.sin(θ)).transpose('r','θ',...)
v = np.max(np.array([p[j].max().values for j in range(nTracs)]))



'''-----------------------------------------------------------------------------
----------PLOT------------------------------------------------------------------
-----------------------------------------------------------------------------'''
#tinker
bfs = 14
blw = 2

def initialize_plots():
    #declare variables
    global fSep, aSep, pSep, cSep
    global cmapList

    #declare plots
    fSep = plt.figure(figsize=(10,6), constrained_layout=True)
    gSep = fSep.add_gridspec(nTracs, nTimes+1, width_ratios=(1, 1, 1, 1, .15))
    aSep = np.array([[None] * nTimes for j in range(nTracs)])
    cSep = [None] * nTracs
    for j1 in range(nTracs):
        cSep[j1] = fSep.add_subplot(gSep[j1,-1])
        for j2 in range(nTimes):
            aSep[j1,j2] = fSep.add_subplot(gSep[j1,j2], projection='polar')

    #label axes
    for a in aSep[-1]: a.set_xlabel(r'$x$', fontsize=bfs)
    for a in aSep: a[0].set_ylabel(r'$y$', fontsize=bfs)
    for j in range(nTimes): aSep[0,j].set_title(r'$t = {0:d}$'.format(times[j]))

    #turn labels off
    for a in aSep.flatten():
        a.set_xticklabels([])
        a.set_yticklabels([])

    #colors
    plotBl = plt.cm.get_cmap('tab10')(0.)[:-1]
    plotOr = plt.cm.get_cmap('tab10')(.1)[:-1]
    arrayBl = np.hstack([np.outer(np.ones(256), plotBl), np.outer(np.linspace(0,1,256), np.ones(1))])
    arrayOr = np.hstack([np.outer(np.ones(256), plotOr), np.outer(np.linspace(0,1,256), np.ones(1))])
    cmapBl = ListedColormap(arrayBl)
    cmapOr = ListedColormap(arrayOr)
    cmapList = [cmapBl, cmapOr]

    #prepare to store plots for legends
    pSep = [None] * nTracs

def tidy_up_plots():
    #colorbar
    for j in range(nTracs):
        cbar = plt.colorbar(pSep[j], cax=cSep[j])

    #grids
    for a in aSep.flatten(): a.grid()

    #save
    if saveFigures:
        today = np.datetime64('today').item()
        todayStr = '{0:02d}{1:02d}'.format(today.month, today.day)
        plt.figure(fSep.number)
        #plt.savefig(fname)
        plt.savefig('../figures/'+todayStr+'_fig1.png')

if __name__ == "__main__":
    try:
        initialize_plots()

        for j1 in range(nTracs):
            for j2 in range(nTimes):
                p_j = p[j1].sel(time=times[j2], method='nearest')
                pSep[j1] = aSep[j1,j2].pcolormesh(θ, r, p_j, vmin=0, vmax=v,
                                                  shading='gouraud', cmap=cmapList[j1])

        tidy_up_plots() 
        plt.show()
    finally:
        plt.close('all')

