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
import h5py
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.animation as movie
from matplotlib.colors import Normalize, LogNorm, to_rgba, ListedColormap
from functools import partial

#tinker
saveStill = True
saveMovie = True
filenames = ['../output/trac2d_top.jld2', '../output/trac2d_mid.jld2']

#relevant variables for plot parameters
ntr = 4
ndp = len(filenames)

#init

#style
def gen_op_cmap(c, α0=0):
    if isinstance(c, str):
        rgb = to_rgba(c)[:-1]
    else:
        rgb = c[0:3]
    arr = np.hstack([
        np.outer(np.ones(256), rgb),
        np.outer(np.linspace(α0,1,256), np.ones(1))
    ])
    return ListedColormap(arr)

def gen_white_cmap(c, α0=0):
    if isinstance(c, str):
        rgb = to_rgba(c)[:-1]
    else:
        rgb = c[0:3]
    rgb = np.array(rgb)
    arr = np.hstack([
        np.linspace((1-α0)*(1-rgb)+rgb, rgb, 256),
        np.ones((256,1))
    ])
    return ListedColormap(arr)



'''-----------------------------------------------------------------------------
----------CODE------------------------------------------------------------------
-----------------------------------------------------------------------------'''
#read data
csep = []
for i, filename in enumerate(filenames):
    with h5py.File(filename, 'r') as f:
        if i == 0:
            #time
            t_keys = list(f['timeseries/t'].keys())
            t_keys = [t_keys[i] for i in np.argsort([float(ti) for ti in t_keys])]
            t_keys = t_keys[::10]
            t = np.array([f['timeseries/t/{0}'.format(t_key)][()] for t_key in t_keys])

            #space
            Hx = f['grid/Hx'][()]
            Hy = f['grid/Hy'][()]
            x = f['grid/xᶜᵃᵃ'][Hx:-Hx]
            y = f['grid/yᵃᶜᵃ'][Hy:-Hy]

        #tracers
        csep.append(np.stack([np.stack([
            np.squeeze(f['timeseries/c{0}/{1}'.format(j, t_key)])
            for t_key in t_keys])
            for j in range(ntr)
        ]))

ccom = np.stack(csep)
cmax = np.nanmax(ccom)
cmin = 0 
    


'''-----------------------------------------------------------------------------
----------PLOT------------------------------------------------------------------
-----------------------------------------------------------------------------'''
def initialize_plots():
    #declare plots
    plt.style.use('mason')
    fC = plt.figure(figsize=(10, 8), layout='constrained')
    gC = fC.add_gridspec(ndp+1, ntr+1,
                         height_ratios=[1]+[9]*ndp, width_ratios=[9]*ntr+[1])
    aC = np.array([np.array([fC.add_subplot(gC[i, j])
                   for j in range(ntr)])
                   for i in range(1, ndp+1)])
    cC = fC.add_subplot(gC[1:, -1])
    tC = fC.add_subplot(gC[0, :-1])

    #share axes
    for i in range(ndp):
        for j in range(ntr):
            if i < ndp - 1:
                aC[i, j].sharex(aC[-1, j])
                aC[i, j].tick_params(labelbottom=False)
            if j > 0:
                # aC[i, j].sharey(aC[i, 0])
                aC[i, j].tick_params(labelleft=False)

    #label axes
    for a in aC[-1]:
        a.set_xlabel(r'$x$ [m]')
    for j in range(ndp): 
        aC[j, 0].set_ylabel(
            r'$z =$' + ' {0:.0f}\n'.format(j) + r'$y$ [deg]'
        )
    tC.set_title('time since release: {0:02d} hours'.format(0))

    #tick formatting
    tC.set_xticks([])
    tC.set_yticks([])
    
    #prepare to store plots for legends
    pC = [None for _ in range(ntr*ndp + 1)]

    return fC, aC, pC, cC, tC

def tidy_up_plots(fC, aC, pC, cC):
    #colorbars
    cbar = plt.colorbar(pC[0], cax=cC)
    cbar.set_label('$p(x, y)$')

def make_movie_plots(i, plots):
    for r in range(ndp):
        for c in range(ntr):
            j = r*ntr + c
            plots[j].set_array(ccom[r, c, i, ...])
    plots[-1].set_array(np.concatenate(
        (np.ones((2, i)), np.zeros((2, t.size - 1 - i))),
        axis=1
    ))
    # tC.set_title('time since release: {0:02d} hours'.format((i+1)//2))
    return pC

fC, aC, pC, cC, tC = initialize_plots()
make_movie = partial(make_movie_plots, plots=pC)

#styling
cmap = gen_white_cmap(plt.cm.get_cmap('tab10')(0))
ctime = gen_op_cmap('black')
cnorm = Normalize(vmin=cmin, vmax=cmax)

for r in range(ndp):
    for c in range(ntr):
        j = r*ntr + c
        pC[j] = aC[r, c].pcolormesh(x, y, ccom[r, c, 0, ...], norm=cnorm, cmap=cmap)#,
                                    # cmap=cmaps,
                                    # norm=cnorm)
pC[-1] = tC.pcolormesh(t[:-1], np.arange(2), np.zeros((2, t.size-1)),
                       cmap=ctime, vmin=0, vmax=1)

tidy_up_plots(fC, aC, pC, cC) 

if saveMovie:
    m = movie.FuncAnimation(fC, make_movie, frames=t.size, blit=False)
    Writer = movie.writers['ffmpeg_file']
    writer = Writer(fps=15, metadata=dict(artist='Mason'), bitrate=1500)
    mname = 'ht_sbs_movie.mp4'
    m.save(mname, writer=writer)
if saveStill:
    fname = 'ht_sbs_still.png'
    plt.savefig(fname, dpi=200, transparent=False)

# plt.show()


'''
if e != 0 else
                (l + '$') if l != '$' else
                    '1'
'''
