#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 30 2021

@author: Mason Rogers
"""

#tinker
write_movie = True
mname = "../figures/1121_still.mp4"

#imports
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.animation as movie
from rc_param import *
from read_MITgcm import ds
from matplotlib.colors import LogNorm, to_rgba, ListedColormap

#define plot variable
nTracs = 1

#select plane
p = list(ds.values())[0].interp(YC=0) * dy
v = np.max(np.array([p['TRAC0'+str(i)].max().values for i in range(1, nTracs+1)])) 
#cnorm = LogNorm(v*1e-4,v)
print(v)

#fix Z-coordinate
p['Z'] = np.arange(H-dz/2,0,-dz)
p['Zl'] = np.arange(H,0,-dz)

#transpose for pyplot
p = p.transpose('Z','XC',...)



'''-----------------------------------------------------------------------------
--------PLOTTING FUNCTIONS------------------------------------------------------
-----------------------------------------------------------------------------'''
#tinker
bfs = 14
blw = 1.5

def initialize_plots():
    #declare variables
    global f_p, a_p, p_p
    global new_pinks, new_oranges, new_limes, new_aquas
    global cmaps

    #declare plots
    f_p = plt.figure(figsize=(10,7), constrained_layout=True)
    g_p = f_p.add_gridspec(1,1)
    a_p = f_p.add_subplot(g_p[0,0], projection='rectilinear')

    #label axes
    a_p.set_xlabel(r'$x$', fontsize=bfs)
    a_p.set_ylabel(r'$z$', fontsize=bfs) 

    #titles
    a_p.set_title("Microplastics in rotating can", fontsize=bfs+2)

    #prepare to store plots for legends
    p_p = [None] * nTracs

    #define opacity colormaps
    deeppink = to_rgba('deeppink')[:-1]
    orange = to_rgba('orange')[:-1]
    lime = to_rgba('lime')[:-1]
    aqua = to_rgba('aqua')[:-1]
    array_pinks = np.hstack([np.outer(np.ones(256),deeppink), np.outer(np.linspace(0,1,256), np.ones(1))])
    array_oranges = np.hstack([np.outer(np.ones(256),orange), np.outer(np.linspace(0,1,256), np.ones(1))])
    array_limes = np.hstack([np.outer(np.ones(256),lime), np.outer(np.linspace(0,1,256), np.ones(1))])
    array_aquas = np.hstack([np.outer(np.ones(256),aqua), np.outer(np.linspace(0,1,256), np.ones(1))])
    new_pinks = ListedColormap(array_pinks)
    new_oranges = ListedColormap(array_oranges)
    new_limes = ListedColormap(array_limes)
    new_aquas = ListedColormap(array_aquas)
    cmaps = [new_pinks, new_oranges, new_limes, new_aquas]

def tidy_up_plots():
    #colorbars
    plt.colorbar(p_p[0], ax=a_p, label=r"$p$ [1/unit$^2$]")

def makemovie(i):
    for j in range(nTracs):
        p_p[j].set_array(p['TRAC0'+str(j+1)].isel(time=i)) 
    return p_p

def startmovie():
    return p_p

if __name__ == "__main__":
    try:
        initialize_plots()
        dname = "eps"+"{0:02d}".format(int(100*δ))+"/"
        for i in range(20):
            x_fname = dname+"poincare_x_"+str(i)+".bin"
            z_fname = dname+"poincare_z_"+str(i)+".bin"
            with open(x_fname, "r") as file:
                x_pts = np.fromfile(file)
            with open(z_fname, "r") as file:
                z_pts = np.fromfile(file)
            a_p.scatter(x_pts, z_pts, c='black', s=4) 
        for i in range(nTracs-1,-1,-1):
            p_p[i] = a_p.pcolormesh(p.XC, p.Z, p['TRAC0'+str(i+1)].isel(time=0),
                                    cmap=cmaps[i], shading='gouraud', vmin=0, vmax=v, animated=True)
        plt.scatter(1/3,.5,c='green',s=10)
        tidy_up_plots()

        if write_movie:
            m = movie.FuncAnimation(f_p, makemovie, init_func=startmovie, frames=p.time.size, blit=True)
            Writer = movie.writers['ffmpeg_file']
            writer = Writer(fps=15, metadata=dict(artist='Mason'), bitrate=1500)
            m.save(mname, writer=writer)

        plt.show()
    finally:
        plt.close('all')
