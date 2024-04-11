#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
kv_figc.py plots the variance growth and snapshots for two configurations of the
Kauffman vortex experiment.

Created on Mon Oct 31 2022

@author: Mason Rogers
"""
'''-----------------------------------------------------------------------------
----------INIT------------------------------------------------------------------
-----------------------------------------------------------------------------'''

#tinker
saveFigures = True
pFiles = ['p_small.py', 'p_large.py']
logNorm = False

#imports
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from dict_moments import calc_moms
from dict_MITgcm import ds, gr, dirs
from matplotlib.colors import to_rgba, ListedColormap, LogNorm, Normalize

#relevant variables for plot parameters
nConfs = 2
nPerCf = 4
nSnaps = 1

#read files
ϵ = np.zeros(nConfs)
c = np.zeros(nConfs)
Σc = np.zeros(nConfs)
s = np.zeros(nConfs)
τ = np.zeros(nConfs)
BB = np.zeros(nConfs)
RR = np.zeros(nConfs)
aa = np.zeros(nConfs)
for j in range(nConfs):
    with open(pFiles[j], 'r') as f:
        exec(f.read())
        BB[j] = B
        ϵ[j] = ((1+2*B)*d**2*Us)/(36*ν*Ls) #small parameter
        c[j] = Ls*ϵ[j]/Us * 2*(1-B)/(1+2*B) * (Γ/(2*np.pi))**2 * a**-4
        Σc[j] = 2*κ/c[j]
        s[j] = 4*κ
        τ[j] = 1/(2*c[j])
        RR[j] = R
        aa[j] = a
B = BB
R = RR
a = aa
R[0] *= 1/5
R[1] *= 2/3

#prepare to store polar arrays
p = [[[None]*2 for j2 in range(nSnaps)] for j1 in range(nConfs)]
v = np.zeros(nConfs)



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

#interpolate to polar coordinates
for ds_j, j1 in zip(ds.values(), range(nConfs)):
    θ = xr.DataArray(data=np.linspace(0,2*np.pi,endpoint=True),
                     dims=['θ'],
                     coords={'θ': np.linspace(0,2*np.pi,endpoint=True)})
    r = xr.DataArray(data=np.linspace(0,R[j1],endpoint=False),
                     dims=['r'],
                     coords={'r': np.linspace(0,R[j1],endpoint=False)})
    ds_j = ds_j.isel(Z=0)
    for j2 in range(nSnaps):
        tmp = ds_j['TRAC{0:02d}'.format(3*j1+1)].interp(XC=r*np.cos(θ), YC=r*np.sin(θ))
        p[j1][j2][0] = tmp.isel(time=0)
        p[j1][j2][1] = tmp.isel(time=-1)
        v[j1] = np.maximum(v[j1], np.maximum(p[j1][j2][0].max().values, p[j1][j2][1].max().values))



'''-----------------------------------------------------------------------------
----------PLOT------------------------------------------------------------------
-----------------------------------------------------------------------------'''
#tinker
bfs = 14
blw = 2

def initialize_plots():
    #declare variables
    global fTot, aVar, pVar, aSnp, cSnp, pSnp
    global cmapList, cnorms

    #declare plots
    fTot = plt.figure(figsize=(13,8), constrained_layout=True)
    gTot = fTot.add_gridspec(nSnaps+1, 3*nConfs,
                             height_ratios=[1.]+[.5]*nSnaps, width_ratios=[.4,.4,.1]*nConfs)
    aVar = np.array([fTot.add_subplot(gTot[0,3*j:3*j+3]) for j in range(nConfs)])
    aSnp = np.array([np.array([fTot.add_subplot(gTot[j2,j1+(j1+1)//3], projection='polar') \
                     for j1 in range(2*nConfs)]) for j2 in range(1, nSnaps+1)])

    #label axes
    for a in aVar:
        a.set_xlabel('time', fontsize=bfs)
        a.set_ylabel('$\mathbb{E}[R^2]$', fontsize=bfs)
    #aSnp[-1,0].set_xlabel('$x$', fontsize=bfs)
    #aSnp[-1,0].set_ylabel('$y$', fontsize=bfs)

    #turn labels off
    for a in aSnp.flatten():
        a.set_xticklabels([])
        a.set_yticklabels([])

    #titles
    aVar[0].set_title('(a) centered release', fontsize=bfs)
    aVar[1].set_title('(b) far-field release', fontsize=bfs)
    for j in range(nConfs):
        aSnp[0,2*j].set_title('initial', fontsize=bfs-2)
        aSnp[0,2*j+1].set_title('final', fontsize=bfs-2)

    #colors
    plotRGBs = [plt.cm.get_cmap('tab10')(3*j/10)[:-1] for j in range(nSnaps)]
    cmapArrs = [np.hstack([np.outer(np.ones(256), pRGB_j), np.outer(np.linspace(0,1,256), np.ones(1))]) \
                for pRGB_j in plotRGBs]
    cmapList = [ListedColormap(cArr_j) for cArr_j in cmapArrs]

    #colormap norms
    cnorms = [LogNorm(v_j**1e-2, v_j) for v_j in v] if logNorm else [Normalize(0, v_j) for v_j in v]

    #colorbars
    cSnp = [[None for j2 in range(nSnaps)] for j1 in range(nConfs)]
    for j1 in range(nConfs):
        for j2 in range(nSnaps):
            cSnp[j1][j2] = fTot.add_subplot(gTot[j2+1, 3*j1+2])

    #prepare to store plots for legends
    pVar = [[None] * (nPerCf+1) for j in range(2)]
    pSnp = [[[None] * 2 for j2 in range(nSnaps)] for j1 in range(nConfs)]

def Σp_lines_small(t, j):
    ll, ul = aVar[j].get_ylim()
    #3 cases: Σc < ll, ll < Σc < ul, ul < Σc
    lΣ = np.minimum((ll-Σc[j])*np.exp(2*c[j]*np.max(t)) + Σc[j], ll)
    uΣ = np.maximum((ul-Σc[j])*np.exp(2*c[j]*np.max(t)) + Σc[j], ul)
    Σ0 = np.sinh(np.linspace(np.arcsinh(lΣ), np.arcsinh(uΣ), 20)) 
    for Σ0_j in Σ0:
        pVar[j][nPerCf] = aVar[j].plot(t, Σp_small(t, j, Σ0_j), 
                                         color='grey', alpha=0.5, lw=blw/2)[0]
    aVar[j].set_ylim(ll, ul)

def Σp_lines_large(t, j):
    ll, ul = aVar[j].get_ylim()
    lΣ = ll - s[j]*np.max(t)
    uΣ = ul
    Σ0 = np.linspace(lΣ, uΣ, 20)
    for Σ0_j in Σ0:
        pVar[j][nPerCf] = aVar[j].plot(t, Σp_large(t, j, Σ0_j),
                                         color='grey', alpha=0.5, lw=blw/2,
                                         linestyle='dashed')[0]
    aVar[j].set_ylim(ll, ul)

def tidy_up_plots():
    #legends
    labels0 = ['full Monte Carlo', 'reduced Monte Carlo', 'reduced Dedalus', 'reduced MITgcm', r'$R \ll a$ prediction']
    labels1 = ['full Monte Carlo', 'reduced Monte Carlo', 'reduced Dedalus', 'reduced MITgcm', r'$R \gg a$ prediction']
    aVar[0].legend(pVar[0], labels0, title='Simulation', fontsize=bfs-2, title_fontsize=bfs-2)
    aVar[1].legend(pVar[1], labels1, title='Simulation', fontsize=bfs-2, title_fontsize=bfs-2)

    #colorbars
    for j1 in range(nConfs):
        for j2 in range(nSnaps):
            cbar = plt.colorbar(pSnp[j1][j2][0], cax=cSnp[j1][j2])
            cbar.ax.tick_params(labelsize=bfs-4)
            cbar.ax.set_ylabel(r'$p$', fontsize=bfs)

    #grid lines
    for j1 in range(nConfs):
        for j2 in range(nSnaps):
            for a in aSnp[j2,2*j1:2*j1+2]:
                a.set_xticks(np.pi/2 * np.arange(4))
                a.set_yticks([R[j1]/4, R[j1]/2, 3*R[j1]/4])
                a.grid(alpha=.3)
        aSnp[-1,2*j1].set_yticklabels(['{0:.1f}'.format(r) for r in [R[j1]/4, R[j1]/2, 3*R[j1]/4]],
                                      position=(np.pi,0), fontsize=bfs-4, ha='center') 

    #save
    if saveFigures:
        today = np.datetime64('today').item()
        todayStr = '{0:02d}{1:02d}'.format(today.month, today.day)
        plt.figure(fTot.number) 
        plt.savefig('../figures/'+todayStr+'_fig2.png')

if __name__ == "__main__":
    try:
        initialize_plots()

        #variance plots
        t0 = ds[dirs[0]]['time']
        t1 = ds[dirs[1]]['time']
        aVar[0].set_xlim([np.min(t0), np.max(t0)])
        aVar[1].set_xlim([np.min(t1), np.max(t1)])
        for j in range(1, nPerCf+1): 
            pVar[0][j-1], = aVar[0].plot(t0, ds[dirs[0]]['Σo{0:02d}'.format(j)], lw=blw)
            pVar[1][j-1], = aVar[1].plot(t1, ds[dirs[1]]['Σo{0:02d}'.format(j)], lw=blw)
        Σp_lines_small(ds[dirs[0]]['time'].values, 0)
        Σp_lines_large(ds[dirs[1]]['time'].values, 1)
        
        #snapshots
        for j1 in range(nConfs):
            for j2 in range(nSnaps):
                tmp = p[j1][j2]
                pSnp[j1][j2][0] = aSnp[j2,2*j1].pcolormesh(tmp[0].θ, tmp[0].r, tmp[0],
                                  shading='gouraud', cmap=cmapList[j2], norm=cnorms[j1])
                pSnp[j1][j2][1] = aSnp[j2,2*j1+1].pcolormesh(tmp[1].θ, tmp[1].r, tmp[1],
                                  shading='gouraud', cmap=cmapList[j2], norm=cnorms[j1])

        tidy_up_plots() 
        plt.show()
    finally:
        plt.close('all')


