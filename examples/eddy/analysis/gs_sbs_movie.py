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
saveStill = True
saveMovie = True
pFiles = ['p_small.py', 'p_large.py']

#imports
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.animation as movie
from dict_MITgcm import ds, gr, dirs

#relevant variables for plot parameters
nConfs = len(ds.keys())
zSnaps = np.array([0, 19, 25])
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



'''-----------------------------------------------------------------------------
----------CODE------------------------------------------------------------------
-----------------------------------------------------------------------------'''
p = {}
pMax = 0
for k in ds.keys():
    pTot = (ds[k]['TRAC01']*ds[k]['drF']*ds[k]['rA']).isel(time=0).sum()
    p[k] = (ds[k]['TRAC01']).isel(Z=zSnaps) / pTot
    pMax = np.maximum(pMax, p[k].max())


'''-----------------------------------------------------------------------------
----------PLOT------------------------------------------------------------------
-----------------------------------------------------------------------------'''
#tinker
blw = 2

def initialize_plots():
    #declare plots
    fC = plt.figure(figsize=(8, 6), layout='constrained')
    gC = fC.add_gridspec(nSnaps+1, nConfs+1,
                         height_ratios=[1]+[9]*nSnaps, width_ratios=[9]*nConfs+[1])
    aC = np.array([np.array([fC.add_subplot(gC[i,j])
                   for j in range(nConfs)])
                   for i in range(1, nSnaps+1)])
    cC = fC.add_subplot(gC[1:,-1])
    tC = fC.add_subplot(gC[0, :-1])

    #label axes
    for a in aC[-1]: a.set_xlabel('longitude [deg]')
    for j in range(nSnaps): 
        aC[j,0].set_ylabel(
            '$z = {0:.0f}$ m\nlatitude [deg]'.format(ds[k]['Z'][zSnaps[j]])
        )
    aC[0,0].set_title('fluid parcels')
    aC[0,1].set_title('inertial particles')
    tC.set_title('time since release: {0:02d} days'.format(0))

    #tick formatting
    for a in aC.flatten():
        a.set_xticks(np.unique(np.round(ds[k]['XC'] + x0)))
        a.set_yticks(np.unique(np.round(ds[k]['YC'] + y0)))
        a.xaxis.set_major_formatter(lambda x, pos: '{0:.0f}'.format(x))
        a.yaxis.set_major_formatter(lambda x, pos: '{0:.0f}'.format(x))
    tC.set_xticks([])
    tC.set_yticks([])

    #spines
    tC.spines[['top', 'right']].set_visible(True)
    
    #prepare to store plots for legends
    pC = [None for _ in range(nConfs*nSnaps + 1)]

    return fC, aC, pC, cC, tC

def tidy_up_plots(fC, aC, pC, cC):
    #colorbars
    cbar = plt.colorbar(pC[0], cax=cC)
    cbar.set_label('$p(x, y) \ [{\sf 10^{-12} m^{-3}}]$')
    yTicks = np.union1d(cC.get_yticks(), cC.get_ylim())
    yTicks = yTicks[(yTicks >= cC.get_ylim()[0]) & (yTicks <= cC.get_ylim()[1])]
    yExp = np.floor(np.log10(yTicks))
    yMul = np.round(yTicks / 10**yExp)
    yLabels = [r'${0:.0f} \times '.format(m) if m > 1 else
                '$' for m in yMul]
    yLabels = [(l + '10^{' + '{0:.0f}'.format(e) + '}$') for e, l in zip(yExp, yLabels)]
    cC.set_yticks(yTicks)
    cC.set_yticklabels(yLabels)

    #share axes
    for a in aC.flatten():
        a.set_xlim(aC[0,0].get_xlim())
        a.set_ylim(aC[0,0].get_ylim())
    for a in aC:
        for j in range(1,nConfs):
            a[j].set_yticklabels(['']*len(a[j].get_yticks()))
    for a in aC[:-1,:].flatten():
        a.set_xticklabels(['']*len(a.get_xticks()))

def make_movie(i):
    for r in range(nSnaps):
        for k, c in zip(ds.keys(), range(nConfs)):
            j = r*nConfs + c
            pC[j].set_array(p[k].isel(Z=r, time=i) * 1e12)
    pC[-1].set_array(np.concatenate(
        (np.ones((2, i)), np.zeros((2, p[k]['time'].size-1-i))),
        axis=1
    ))
    tC.set_title('time since release: {0:02d} days'.format((i+1)//2))
    return pC

if __name__ == "__main__":
    try:
        fC, aC, pC, cC, tC = initialize_plots()

        #styling
        cmaps = gen_white_cmap(plt.cm.get_cmap('tab10')(0))
        cTopo = gen_op_cmap('black')
        cnorm = LogNorm(vmin=0.04, vmax=40)
        
        for r in range(nSnaps):
            for k, c in zip(ds.keys(), range(nConfs)):
                j = r*nConfs + c
                pC[j] = aC[r,c].pcolormesh(ds[k]['XC'] + x0,
                                           ds[k]['YC'] + y0,
                                           p[k].isel(Z=r, time=0) * 1e12,
                                           cmap=cmaps,
                                           norm=cnorm)
                aC[r,c].pcolormesh(ds[k]['XC'] + x0,
                                   ds[k]['YC'] + y0,
                                   ds[k]['hFacC'].isel(Z=zSnaps[r]) == 0,
                                   cmap=cTopo)
        pC[-1] = tC.pcolormesh(p[k]['time'].values[:-1],
                               np.arange(2),
                               np.zeros((2, p[k]['time'].size-1)),
                               cmap=cTopo,
                               vmin=0,
                               vmax=1)

        tidy_up_plots(fC, aC, pC, cC) 

        today = np.datetime64('today').item()
        if saveMovie:
            m = movie.FuncAnimation(fC, make_movie, frames=p[k].time.size, blit=True)
            Writer = movie.writers['ffmpeg_file']
            writer = Writer(fps=15, metadata=dict(artist='Mason'), bitrate=1500)
            mname = '../figures/gs_sbs_movie_x.mp4'
            m.save(mname, writer=writer, dpi=200)
        if saveStill:
            fC.delaxes(tC)
            fname = '../figures/gs_sbs_still_x.png'
            plt.savefig(fname, dpi=200, transparent=False)

        plt.show()
    finally:
        plt.close('all')


'''
if e != 0 else
                (l + '$') if l != '$' else
                    '1'
'''
