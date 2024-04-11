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
pFiles = ['p_small.py', 'p_large.py']

#imports
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from dict_moments import calc_moms
from dict_MITgcm import ds, gr, dirs

#relevant variables for plot parameters
nConfs = 2
nPerCf = 4

#read files
ϵ = np.zeros(nConfs)
c = np.zeros(nConfs)
Σc = np.zeros(nConfs)
s = np.zeros(nConfs)
τ = np.zeros(nConfs)
BB = np.zeros(nConfs)
for j in range(nConfs):
    with open(pFiles[j], 'r') as f:
        exec(f.read())
        BB[j] = B
        ϵ[j] = ((1+2*B)*d**2*Us)/(36*ν*Ls) #small parameter
        c[j] = Ls*ϵ[j]/Us * 2*(1-B)/(1+2*B) * (Γ/(2*np.pi))**2 * a**-4
        Σc[j] = 2*κ/c[j]
        s[j] = 4*κ
        τ[j] = 1/(2*c[j])
B = BB



'''-----------------------------------------------------------------------------
----------CODE------------------------------------------------------------------
-----------------------------------------------------------------------------'''

#compute variance timeseries
for ds_j, gr_j in zip(ds.values(), gr.values()):
    moms = calc_moms(ds_j, gr_j)
    for j in range(1, nPerCf+1):
        EX2 = moms['TRAC{0:02d}'.format(j)]['Σ_XC_XC'] + moms['TRAC{0:02d}'.format(j)]['μ_XC']**2 
        EY2 = moms['TRAC{0:02d}'.format(j)]['Σ_YC_YC'] + moms['TRAC{0:02d}'.format(j)]['μ_YC']**2
        ds_j['Σo{0:02d}'.format(j)] = EX2 + EY2

#predicted timeseries for large and small Σ limits
def Σp_large(t, j, Σ0):
    return s[j]*t + Σ0
def Σp_small(t, j, Σ0):
    return (Σ0 - Σc[j])*np.exp(-2*c[j]*t) + Σc[j]



'''-----------------------------------------------------------------------------
----------PLOT------------------------------------------------------------------
-----------------------------------------------------------------------------'''
#tinker
bfs = 14
blw = 2

def initialize_plots():
    #declare variables
    global f_var, a_var, p_var

    #declare plots
    f_var, a_var = plt.subplots(figsize=(10,6), ncols=2, constrained_layout=True)

    #label axes
    for a in a_var:
        a.set_xlabel('time', fontsize=bfs)
        a.set_ylabel('$\mathbb{E}[R^2]$', fontsize=bfs)

    #prepare to store plots for legends
    p_var = [[None] * (nPerCf+1) for j in range(2)]

def Σp_lines_small(t, j):
    ll, ul = a_var[j].get_ylim()
    #3 cases: Σc < ll, ll < Σc < ul, ul < Σc
    lΣ = np.minimum((ll-Σc[j])*np.exp(2*c[j]*np.max(t)) + Σc[j], ll)
    uΣ = np.maximum((ul-Σc[j])*np.exp(2*c[j]*np.max(t)) + Σc[j], ul)
    Σ0 = np.sinh(np.linspace(np.arcsinh(lΣ), np.arcsinh(uΣ), 20)) 
    for Σ0_j in Σ0:
        p_var[j][nPerCf] = a_var[j].plot(t, Σp_small(t, j, Σ0_j), 
                                         color='grey', alpha=0.5, lw=blw/2)[0]
    a_var[j].set_ylim(ll, ul)

def Σp_lines_large(t, j):
    ll, ul = a_var[j].get_ylim()
    lΣ = ll - s[j]*np.max(t)
    uΣ = ul
    Σ0 = np.linspace(lΣ, uΣ, 20)
    for Σ0_j in Σ0:
        p_var[j][nPerCf] = a_var[j].plot(t, Σp_large(t, j, Σ0_j),
                                         color='grey', alpha=0.5, lw=blw/2,
                                         linestyle='dashed')[0]
    a_var[j].set_ylim(ll, ul)

def tidy_up_plots():
    #legends
    labels0 = ['full Monte Carlo', 'reduced Monte Carlo', 'reduced Dedalus', 'reduced MITgcm', r'$R \ll a$ prediction']
    labels1 = ['full Monte Carlo', 'reduced Monte Carlo', 'reduced Dedalus', 'reduced MITgcm', r'$R \gg a$ prediction']
    a_var[0].legend(p_var[0], labels0, title='Simulation', fontsize=bfs-2, title_fontsize=bfs-2)
    a_var[1].legend(p_var[1], labels1, title='Simulation', fontsize=bfs-2, title_fontsize=bfs-2)

    #save
    if saveFigures:
        today = np.datetime64('today').item()
        todayStr = '{0:02d}{1:02d}'.format(today.month, today.day)
        plt.figure(f_var.number) 
        plt.savefig('../figures/'+todayStr+'_fig2.png')

if __name__ == "__main__":
    try:
        initialize_plots()

        t0 = ds[dirs[0]]['time']
        t1 = ds[dirs[1]]['time']

        a_var[0].set_xlim([np.min(t0), np.max(t0)])
        a_var[1].set_xlim([np.min(t1), np.max(t1)])

        for j in range(1, nPerCf+1): 
            p_var[0][j-1], = a_var[0].plot(t0, ds[dirs[0]]['Σo{0:02d}'.format(j)], lw=blw)
            p_var[1][j-1], = a_var[1].plot(t1, ds[dirs[1]]['Σo{0:02d}'.format(j)], lw=blw)

        Σp_lines_small(ds[dirs[0]]['time'].values, 0)
        Σp_lines_large(ds[dirs[1]]['time'].values, 1)
        
        #p_var[-2] = a_var.axhline(Σc, color='red', lw=blw/2)

        tidy_up_plots() 
        plt.show()
    finally:
        plt.close('all')

