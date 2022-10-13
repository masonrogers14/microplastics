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
fname = '../figures/kv_fig1.png'

#imports
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
#from read_MITgcm import ds, gr

#relevant variables for plot parameters
nTracs = 2



'''-----------------------------------------------------------------------------
----------CODE------------------------------------------------------------------
-----------------------------------------------------------------------------'''



'''-----------------------------------------------------------------------------
----------PLOT------------------------------------------------------------------
-----------------------------------------------------------------------------'''
#tinker
bfs = 14
blw = 2

def initialize_plots():
    #declare variables
    global f_sep, a_sep, p_sep

    #declare plots
    f_sep, a_sep = plt.subplots(figsize=(10,6), ncols=4, nrows=2,
                                sharey=True, sharex=True,  constrained_layout=True)

    #label axes
    for a in a_sep[-1]: a.set_xlabel(r'$x$', fontsize=bfs)
    for a in a_sep: a[0].set_ylabel(r'$y$', fontsize=bfs)
    for j in range(4): a_sep[0][j].set_title(r'$t = $ {0:d}'.format(j*500))

    #prepare to store plots for legends
    p_sep = [[None] * 4] * 2

def tidy_up_plots():
    #save
    if saveFigures:
        today = np.datetime64('today').item()
        todayStr = '{0:02d}{1:02d}'.format(today.month, today.day)
        plt.figure(f_sep.number)
        plt.savefig(fname)
        #plt.savefig('../figures/'+todayStr+'_fig2.png')

if __name__ == "__main__":
    try:
        initialize_plots()


        tidy_up_plots() 
        plt.show()
    finally:
        plt.close('all')

