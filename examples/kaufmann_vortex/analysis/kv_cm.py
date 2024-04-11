#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
kv_cm.py

Created on Fri Jul 01 2022

@author: Mason Rogers
"""
'''-----------------------------------------------------------------------------
----------INIT------------------------------------------------------------------
-----------------------------------------------------------------------------'''

#tinker
saveFigures = True

#imports
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from kv_param import *
from moments import moms
from read_MITgcm import ds

#relevant variables for plot parameters
nTracs = 1
t = moms['TRAC01']['P']['time']



'''-----------------------------------------------------------------------------
----------CODE------------------------------------------------------------------
-----------------------------------------------------------------------------'''

x = moms['TRAC01']['μ_XC']
y = moms['TRAC01']['μ_YC']
R = np.sqrt(ds['XC']**2 + ds['YC']**2)
ER = (R*ds['TRAC01']*ds['rA']*ds['drF']).sum(dim=['XC','YC','Z'])


'''-----------------------------------------------------------------------------
----------PLOT------------------------------------------------------------------
-----------------------------------------------------------------------------'''
#tinker
bfs = 14
blw = 1.5

def initialize_plots():
    #declare variables
    global f_eor, a_eor, p_eor
    global f_com, a_com, p_com

    #declare plots
    f_eor, a_eor = plt.subplots(figsize=(7,7), constrained_layout=True)
    f_com, a_com = plt.subplots(figsize=(7,7), constrained_layout=True)

    #label axes
    a_eor.set_title(r'Growth in $\mathbb{E}[R]$', fontsize=bfs+2)
    a_eor.set_xlabel('time', fontsize=bfs)
    a_eor.set_ylabel('$\mathbb{E}[R]$', fontsize=bfs)
    a_com.set_title(r'center of mass', fontsize=bfs+2)
    a_com.set_xlabel('$x$', fontsize=bfs)
    a_com.set_ylabel('$y$', fontsize=bfs)

    #prepare to store plots for legends
    p_eor = None
    p_com = None

def tidy_up_plots():
    #save
    if saveFigures:
        today = np.datetime64('today').item()
        todayStr = '{0:02d}{1:02d}'.format(today.month, today.day)
        plt.figure(f_eor.number)
        plt.savefig('../figures/'+todayStr+'_er.png')
        plt.figure(f_com.number)
        plt.savefig('../figures/'+todayStr+'_cm.png')

if __name__ == "__main__":
    try:
        initialize_plots() 
        
        p_eor = a_eor.plot(t, ER)
        p_com = a_com.plot(x, y)

        tidy_up_plots()
        plt.show()
    finally:
        plt.close('all')

