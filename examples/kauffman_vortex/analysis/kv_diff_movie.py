#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 30 2021

@author: Mason Rogers
"""

#tinker
write_movie = True
mname = "../figures/1222_kv_diff.mp4"

#imports
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.animation as movie
from read_MITgcm import ds
from matplotlib.colors import LogNorm, to_rgba, ListedColormap

#define plot variable
nTracs = 1

#polar coordinates
θ = xr.DataArray(data=np.linspace(0,2*np.pi,endpoint=True),
                 dims=['θ'],
                 coords={'θ': np.linspace(0,2*np.pi,endpoint=True)})
r = xr.DataArray(data=np.linspace(0,0.5,endpoint=False),
                 dims=['r'],
                 coords={'r': np.linspace(0,0.5,endpoint=False)})
p = (ds['TRAC01']-ds['TRAC02']).sum(dim='Z').interp(XC=r*np.cos(θ), YC=r*np.sin(θ))
v = np.max(np.abs(p)).values
#cnorm = LogNorm(v*1e-4,v)
print(v)

#transpose for pyplot
p = p.transpose('r','θ',...)



'''-----------------------------------------------------------------------------
--------PLOTTING FUNCTIONS------------------------------------------------------
-----------------------------------------------------------------------------'''
#tinker
bfs = 14
blw = 1.5

def initialize_plots():
    #declare variables
    global f_p, a_p, p_p
    global cnorm

    #declare plots
    f_p = plt.figure(figsize=(10,7), constrained_layout=True)
    g_p = f_p.add_gridspec(1,1)
    a_p = f_p.add_subplot(g_p[0,0], projection='polar')

    #label axes
    a_p.set_rlabel_position(0)
    a_p.text(0,0.25,r"$r$", fontsize=bfs, va='top', ha='center')

    #titles
    a_p.set_title("Microplastics in Kaufmann vortex", fontsize=bfs+2)

    #prepare to store plots for legends
    p_p = [None] * nTracs

    #diverging colormap norm
    cnorm = plt.cm.colors.TwoSlopeNorm(vmin=-v, vcenter=0, vmax=v)

def tidy_up_plots():
    #colorbars
    plt.colorbar(p_p[0], ax=a_p, label=r"$p$ [1/unit$^2$]")
    #save stuff

def makemovie(i):
    for j in range(nTracs):
        p_p[j].set_array(p.isel(time=i))
        #p_p[j].set_clim(v[i]/50, v[i])
    return p_p

def startmovie():
    return p_p

if __name__ == "__main__":
    try:
        initialize_plots()
        for i in range(nTracs-1,-1,-1):
            p_p[i] = a_p.pcolormesh(θ, r, p.isel(time=0),
                                    cmap='RdBu_r', shading='gouraud', norm=cnorm,  animated=True)
        tidy_up_plots()

        if write_movie:
            m = movie.FuncAnimation(f_p, makemovie, init_func=startmovie, frames=p.time.size, blit=True)
            Writer = movie.writers['ffmpeg_file']
            writer = Writer(fps=15, metadata=dict(artist='Mason'), bitrate=1500)
            m.save(mname, writer=writer)

        plt.show()
    finally:
        plt.close('all')
