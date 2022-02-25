#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 30 2021

@author: Mason Rogers
"""

#tinker
write_movie = True
mname = "../figures/0208_short.mp4"

#imports
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.animation as movie
from read_MITgcm import ds
from matplotlib.colors import LogNorm, to_rgba, ListedColormap

#polar coordinates
θ = xr.DataArray(data=np.linspace(0,2*np.pi,endpoint=True),
                 dims=['θ'],
                 coords={'θ': np.linspace(0,2*np.pi,endpoint=True)})
r = xr.DataArray(data=np.linspace(0,0.5,endpoint=False),
                 dims=['r'],
                 coords={'r': np.linspace(0,0.5,endpoint=False)})
p = ds.sum(dim='Z').interp(XC=r*np.cos(θ), YC=r*np.sin(θ))
v = np.maximum(p['TRAC01'].max().values, p['TRAC02'].max().values)
#cnorm = LogNorm(v*1e-4,v)
print(v)

#transpose for pyplot
p = p.transpose('r','θ',...)

#break out plotting variables
p1 = p['TRAC01']
p2 = p['TRAC02']
p3 = p1-p2


'''-----------------------------------------------------------------------------
--------PLOTTING FUNCTIONS------------------------------------------------------
-----------------------------------------------------------------------------'''
#tinker
bfs = 14
blw = 1.5

def initialize_plots():
    #declare variables
    global f_p, a_p, p_p
    global new_pinks, new_oranges, new_limes, new_blues, new_orbu

    #declare plots
    f_p = plt.figure(figsize=(10,4), constrained_layout=True)
    g_p = f_p.add_gridspec(1,3)
    a_p = np.array([f_p.add_subplot(g_p[0,i], projection='polar') for i in range(3)])

    #label axes
    for a in a_p: a.set_rlabel_position(0)
    for a in a_p: a.text(0,0.25,r"$r$", fontsize=bfs, va='top', ha='center')

    #titles
    f_p.suptitle("Microplastics in Kaufmann vortex", fontsize=bfs+2)
    a_p[0].set_title("Monte Carlo", fontsize=bfs)
    a_p[1].set_title("Simplified", fontsize=bfs)
    a_p[2].set_title("Difference", fontsize=bfs)

    #ticks
    for a in a_p: a.tick_params(labelsize=bfs-6)

    #prepare to store plots for legends
    p_p = [None] * 3

    #define opacity colormaps
    deeppink = to_rgba('deeppink')[:-1]
    orange = to_rgba('orange')[:-1]
    lime = to_rgba('lime')[:-1]
    blue = to_rgba('deepskyblue')[:-1]
    array_pinks = np.hstack([np.outer(np.ones(256),deeppink), np.outer(np.linspace(0,1,256), np.ones(1))])
    array_oranges = np.hstack([np.outer(np.ones(256),orange), np.outer(np.linspace(0,1,256), np.ones(1))])
    array_limes = np.hstack([np.outer(np.ones(256),lime), np.outer(np.linspace(0,1,256), np.ones(1))])
    array_blues = np.hstack([np.outer(np.ones(256),blue), np.outer(np.linspace(0,1,256), np.ones(1))])
    array_orbu = np.vstack([array_blues[::-1,:], array_oranges])
    new_pinks = ListedColormap(array_pinks)
    new_oranges = ListedColormap(array_oranges)
    new_limes = ListedColormap(array_limes)
    new_blues = ListedColormap(array_blues)
    new_orbu = ListedColormap(array_orbu)

def tidy_up_plots():
    #colorbars
    plt.colorbar(p_p[0], ax=a_p[0], label=r"$p$ [1/unit$^2$]", location='bottom')    
    plt.colorbar(p_p[1], ax=a_p[1], label=r"$p$ [1/unit$^2$]", location='bottom')    
    plt.colorbar(p_p[2], ax=a_p[2], label=r"$p$ [1/unit$^2$]", location='bottom')    
    
    #grid
    for a in a_p: a.grid(axis='y', alpha=0.4)

def makemovie(i):
    p_p[0].set_array(p1.isel(time=i)) 
    p_p[1].set_array(p2.isel(time=i)) 
    p_p[2].set_array(p3.isel(time=i)) 
    return p_p

def startmovie():
    return p_p

if __name__ == "__main__":
    try:
        initialize_plots()
        p_p[0] = a_p[0].pcolormesh(θ, r, p1.isel(time=0), cmap=new_oranges, shading='gouraud', vmin=0, vmax=v, animated=True)
        p_p[1] = a_p[1].pcolormesh(θ, r, p2.isel(time=0), cmap=new_blues, shading='gouraud', vmin=0, vmax=v, animated=True)
        p_p[2] = a_p[2].pcolormesh(θ, r, p3.isel(time=0), cmap=new_orbu, shading='gouraud', vmin=-v/10, vmax=v/10, animated=True)
        tidy_up_plots()

        if write_movie:
            m = movie.FuncAnimation(f_p, makemovie, init_func=startmovie, frames=p.time.size, blit=True)
            Writer = movie.writers['ffmpeg_file']
            writer = Writer(fps=15, metadata=dict(artist='Mason'), bitrate=1500)
            m.save(mname, writer=writer)

        plt.show()
    finally:
        plt.close('all')
