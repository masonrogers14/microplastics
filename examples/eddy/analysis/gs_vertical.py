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
    p[tName] = gr.integrate(p[tName].where(p['hFacC'] > 0), ['X','Y']) 
v = np.max(np.array([p['TRAC0'+str(i)].isel(time=0).max().values for i in range(1, nTracs+1)])) 
print(v)



'''-----------------------------------------------------------------------------
--------PLOTTING FUNCTIONS------------------------------------------------------
-----------------------------------------------------------------------------'''
#tinker
bfs = 14
blw = 1.5

def initialize_plots():
    #declare variables
    global f_p, a_p, p_p

    #declare plots
    f_p, a_p = plt.subplots(figsize=(10,7), constrained_layout=True)

    #titles
    a_p.set_title("Microplastics (horizontally integrated)", fontsize=bfs+2)
    a_p.set_xlabel("$p$ [1/unit^2]", fontsize=bfs)
    a_p.set_ylabel("$z$ [m]", fontsize=bfs)

    #limits
    a_p.set_ylim([-100,0])

    #plots
    p_p = [None]*nTracs
    for j in range(nTracs):
        p_p[j], = a_p.plot(p['TRAC{0:02d}'.format(j+1)].isel(time=0), p['Z'])

def start_movie():
    return p_p

def make_movie(i):
    for j in range(nTracs):
        p_p[j].set_xdata(p['TRAC{0:02d}'.format(j+1)].isel(time=i)) 
    return p_p

if __name__ == "__main__":
    try:
        initialize_plots()
        if writeMovie:
            m = movie.FuncAnimation(f_p, make_movie, init_func=start_movie, frames=p.time.size, blit=True)
            Writer = movie.writers['ffmpeg_file']
            writer = Writer(fps=15, metadata=dict(artist='Mason'), bitrate=1500)
            today = np.datetime64('today').item()
            mname = '../figures/{0:02d}{1:02d}_vertical.mp4'.format(today.month, today.day)
            m.save(mname, writer=writer)
        plt.show()
    finally:
        plt.close('all')
