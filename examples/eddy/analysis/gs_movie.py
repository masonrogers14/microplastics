#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 08 2023

@author: Mason Rogers
"""

#tinker
writeMovie = True

#imports
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.animation as movie
from read_MITgcm import ds, gr, dirs
from matplotlib.colors import LogNorm, to_rgba, ListedColormap

#define plot variable
nTracs = 1
ds = ds[dirs[0]]
gr = gr[dirs[0]]
p = ds
for i in range(nTracs):
    tName = 'TRAC{0:02d}'.format(i+1)
    p[tName] = gr.integrate(p[tName], 'Z') 
v = np.max(np.array([p['TRAC0'+str(i)].isel(time=0).max().values for i in range(1, nTracs+1)])) 
print(v)



'''-----------------------------------------------------------------------------
--------PLOTTING FUNCTIONS------------------------------------------------------
-----------------------------------------------------------------------------'''
#tinker
bfs = 14
blw = 1.5

def gen_op_cmap(c):
    if isinstance(c, str):
        rgb = to_rgba(c)[:-1]
    else:
        rgb = c[0:3]
    arr = np.hstack([np.outer(np.ones(256), rgb), np.outer(np.linspace(0,1,256), np.ones(1))])
    return ListedColormap(arr)

def initialize_plots():
    #declare variables
    global f_p, a_p, p_p
    global new_pinks, new_oranges, new_limes, new_aquas
    global cmaps

    #declare plots
    f_p, a_p = plt.subplots(figsize=(10,7), constrained_layout=True)

    #titles
    a_p.set_title("fluid parcels in Gulf Stream", fontsize=bfs+2)

    #prepare to store plots for legends
    p_p = [None] * nTracs

    #define opacity colormaps
    #    deeppink = to_rgba('deeppink')[:-1]
    #    orange = to_rgba('orange')[:-1]
    #    lime = to_rgba('lime')[:-1]
    #    aqua = to_rgba('aqua')[:-1]
    #    array_pinks = np.hstack([np.outer(np.ones(256),deeppink), np.outer(np.linspace(0,1,256), np.ones(1))])
    #    array_oranges = np.hstack([np.outer(np.ones(256),orange), np.outer(np.linspace(0,1,256), np.ones(1))])
    #    array_limes = np.hstack([np.outer(np.ones(256),lime), np.outer(np.linspace(0,1,256), np.ones(1))])
    #    array_aquas = np.hstack([np.outer(np.ones(256),aqua), np.outer(np.linspace(0,1,256), np.ones(1))])
    #    new_pinks = ListedColormap(array_pinks)
    #    new_oranges = ListedColormap(array_oranges)
    #    new_limes = ListedColormap(array_limes)
    #    new_aquas = ListedColormap(array_aquas)
    #    cmaps = [new_pinks, new_oranges, new_limes, new_aquas]
    cmaps = [gen_op_cmap(plt.cm.get_cmap('tab10')(i)) for i in range(nTracs)]
    cmaps.append(gen_op_cmap('black'))

def tidy_up_plots():
    #colorbars
    plt.colorbar(p_p[0], ax=a_p, label=r"$p$ [1/unit$^2$]")    

    #save a still
    today = np.datetime64('today').item()
    todayStr = '{0:02d}{1:02d}'.format(today.month, today.day)
    plt.savefig('../figures/'+todayStr+'_still.png')

def makemovie(i):
    for j in range(nTracs):
        p_p[j].set_array(p['TRAC0'+str(j+1)].isel(time=i)) 
    return p_p

def startmovie():
    return p_p

if __name__ == "__main__":
    try:
        initialize_plots()
        for i in range(nTracs-1,-1,-1):
            a_p.pcolormesh(p['XC'], p['YC'], p['Depth'] == 0, cmap=cmaps[-1])
            p_p[i] = a_p.pcolormesh(p['XC'], p['YC'], p['TRAC0'+str(i+1)].isel(time=0),
                                    cmap=cmaps[i], shading='gouraud', vmin=0, vmax=v,  animated=True)
        tidy_up_plots()

        if writeMovie:
            m = movie.FuncAnimation(f_p, makemovie, init_func=startmovie, frames=p.time.size, blit=True)
            Writer = movie.writers['ffmpeg_file']
            writer = Writer(fps=15, metadata=dict(artist='Mason'), bitrate=1500)
            today = np.datetime64('today').item()
            mname = '../figures/{0:02d}{1:02d}_movie.mp4'.format(today.month, today.day)
            m.save(mname, writer=writer)

        plt.show()
    finally:
        plt.close('all')
