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
y_j = 175

ds = ds[dirs[0]]
gr = gr[dirs[0]]
p = ds
v = np.max(np.array([p['TRAC0'+str(i)].isel(time=0).max().values for i in range(1, nTracs+1)])) 
for i in range(nTracs):
    tName = 'TRAC{0:02d}'.format(i+1)
    p[tName] = p[tName].isel(YC=y_j)
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
    global cmaps

    #declare plots
    f_p, a_p = plt.subplots(figsize=(10,7), constrained_layout=True)

    #titles
    a_p.set_title("Microplastics (vertical slice)", fontsize=bfs+2)
    a_p.set_xlabel("longitude", fontsize=bfs)
    a_p.set_ylabel("depth [m]", fontsize=bfs)

    #limits
    a_p.set_ylim([-100,0])

    #prepare to store plots for legends
    p_p = [None] * nTracs

    #define opacity colormaps
    cmaps = [gen_op_cmap(plt.cm.get_cmap('tab10')(i+2)) for i in range(nTracs)]
    cmaps.append(gen_op_cmap('black'))

def tidy_up_plots():
    #colorbars
    plt.colorbar(p_p[0], ax=a_p, label=r"$p$ [1/unit$^2$]")    

    #make a map
    f_m, a_m = plt.subplots(figsize=(10,7), layout='constrained')
    a_m.pcolormesh(p['XC'], p['YC'], p['Depth'] == 0, cmap=cmaps[-1])
    a_m.axhline(p['YC'][y_j])
    plt.figure(f_m)
    today = np.datetime64('today').item()
    plt.savefig('../figures/{0:02d}{1:02d}_map.png'.format(today.month, today.day))

def makemovie(i):
    for j in range(nTracs):
        p_p[j].set_array(p['TRAC{0:02d}'.format(j+1)].isel(time=i)) 
    return p_p

def startmovie():
    return p_p

if __name__ == "__main__":
    try:
        initialize_plots()
        for i in range(nTracs-1,-1,-1):
            p_p[i] = a_p.pcolormesh(p['XC'], p['Z'], p['TRAC{0:02d}'.format(i+1)].isel(time=0),
                                    cmap=cmaps[i], shading='gouraud', vmin=0, vmax=v,  animated=True)
        tidy_up_plots()

        if writeMovie:
            m = movie.FuncAnimation(f_p, makemovie, init_func=startmovie, frames=p.time.size, blit=True)
            Writer = movie.writers['ffmpeg_file']
            writer = Writer(fps=15, metadata=dict(artist='Mason'), bitrate=1500)
            today = np.datetime64('today').item()
            mname = '../figures/{0:02d}{1:02d}_slice.mp4'.format(today.month, today.day)
            m.save(mname, writer=writer)

        plt.show()
    finally:
        plt.close('all')
