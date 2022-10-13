#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
kv_var.py computes E[R^2] for the Kauffman vortex and compares it to a
theoretical prediction.

For P(R^2 << a^2) ~ 1, E[R^2] -> 2κ/c where c = lim(r -> 0) -v_r/r is the
small-r linear coefficient of v_r(r). The convergence is exponential with an
e-folding time of 1/(2c).

For P(R^2 >> a^2) ~ 1, d/dt E[R^2] -> 4κ.

This script measures relevant parameters:
* a^2 (units L^2)
* 2κ/c (L^2)
* dx^2 (L^2) -- need to have adequate resultion for gridded output!
* 1/(2c) (T)
* 4κ (L^2/T)

Created on Tue May 17 2022

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

#compute parameters
ϵ = ((1+2*B)*d**2*Us)/(36*ν*Ls) #small parameter
c = Ls*ϵ/Us * 2*(1-B)/(1+2*B) * (Γ/(2*np.pi))**2 * a**-4
Σc = 2*κ/c
s = 4*κ
τ = 1/(2*c)

#compute variance timeseries
Σo = xr.DataArray(np.zeros((t.size, nTracs)), dims=['time', 'n'],
                  coords={'n': np.arange(1, nTracs+1), 'time': t})
for j in range(1, nTracs+1):
    EX2 = moms['TRAC{0:02d}'.format(j)]['Σ_XC_XC'] + moms['TRAC{0:02d}'.format(j)]['μ_XC']**2 
    EY2 = moms['TRAC{0:02d}'.format(j)]['Σ_YC_YC'] + moms['TRAC{0:02d}'.format(j)]['μ_YC']**2
    Σo.loc[{'n': j}] = EX2 + EY2

#predicted timeseries for large and small Σ limits
def Σp_large(Σ0):
    return s*t + Σ0
def Σp_small(Σ0):
    return (Σ0 - Σc)*np.exp(-2*c*t) + Σc

#compute p-in-target timeseries
P = xr.DataArray(np.zeros((t.size, nTracs, 5)), dims=['time', 'n', 'm'],
                 coords={'n': np.arange(1, nTracs+1), 'time': t,
                         'm': np.array([16, 64, 144, 256, 400])})
R2 = ds['XC']**2 + ds['YC']**2
mask = R2 < np.minimum(R**2, P['m']*a**2) 
vol = ds['rA'] * ds['drF'] * ds['maskC']
for j in range(1, nTracs+1):
    P.loc[{'n': j}] = np.squeeze((ds['TRAC{0:02d}'.format(j)]*vol).where(mask).sum(dim=['XC','YC'])) 



'''-----------------------------------------------------------------------------
----------PLOT------------------------------------------------------------------
-----------------------------------------------------------------------------'''
#tinker
bfs = 14
blw = 1.5

def initialize_plots():
    #declare variables
    global f_var, a_var, p_var
    global f_pin, a_pin, v_pin

    #declare plots
    f_var, a_var = plt.subplots(figsize=(10,7), constrained_layout=True)
    f_pin, a_pin = plt.subplots(figsize=(10,7), constrained_layout=True)

    #label axes
    a_var.set_title(r'Growth in $\mathbb{E}[R^2]$', fontsize=bfs+2)
    a_var.set_xlabel('time', fontsize=bfs)
    a_var.set_ylabel('$\mathbb{E}[R^2]$', fontsize=bfs)
    a_pin.set_title(r'$\mathbb{P}[R^2 < ca^2]$', fontsize=bfs+2)
    a_pin.set_xlabel('time', fontsize=bfs)
    a_pin.set_ylabel('P', fontsize=bfs)

    #prepare to store plots for legends
    p_var = [None] * (nTracs + 4)
    p_pin = [None] * 5

def Σp_lines_small():
    ll, ul = a_var.get_ylim()
    #3 cases: Σc < ll, ll < Σc < ul, ul < Σc
    lΣ = ll#np.minimum((ll-Σc)*np.exp(2*c*t.max().values) + Σc, ll)
    uΣ = ul#np.maximum((ul-Σc)*np.exp(2*c*t.max().values) + Σc, ul)
    Σ0 = np.linspace(lΣ, uΣ, 20)
    for Σ0_j in Σ0:
        p_var[nTracs] = a_var.plot(t, Σp_small(Σ0_j),
                                   color='grey', alpha=0.5, lw=blw/2)
    a_var.set_ylim(ll, ul)

def Σp_lines_large():
    ll, ul = a_var.get_ylim()
    lΣ = ll - s*t.max().values
    uΣ = ul
    Σ0 = np.linspace(lΣ, uΣ, 20)
    for Σ0_j in Σ0:
        p_var[nTracs+1] = a_var.plot(t, Σp_large(Σ0_j),
                                     color='grey', alpha=0.5, lw=blw/2,
                                     linestyle='dashed')
    a_var.set_ylim(ll, ul)

def tidy_up_plots():
    #save
    if saveFigures:
        today = np.datetime64('today').item()
        todayStr = '{0:02d}{1:02d}'.format(today.month, today.day)
        plt.figure(f_var.number)
        plt.savefig('../figures/'+todayStr+'_v.png')
        plt.figure(f_pin.number)
        plt.savefig('../figures/'+todayStr+'_p.png')

if __name__ == "__main__":
    try:
        initialize_plots()

        p_var[:nTracs] = a_var.plot(t, Σo, lw=blw)
        Σp_lines_large()
        Σp_lines_small()
        p_var[-2] = a_var.axhline(Σc, color='red', lw=blw/2)

        p_pin = a_pin.plot(t, P.isel(n=0))    

        tidy_up_plots()
        plt.show()
    finally:
        plt.close('all')

