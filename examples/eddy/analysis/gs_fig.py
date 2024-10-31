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

#imports
from matplotlib.colors import LogNorm, to_rgba, ListedColormap

#tinker
saveFigures = True
pFiles = ['p_small.py', 'p_large.py']

#imports
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from read_MITgcm import ds, gr, dirs

#relevant variables for plot parameters
nConfs = len(ds.keys())
tSnaps = np.array([0]) * 86400.
xSlice = -1 # < 0 --> integrate
ySlice = -1
zSlice = -1

#init
nSnaps = tSnaps.size

# #read files
# ϵ = np.zeros(nConfs)
# c = np.zeros(nConfs)
# Σc = np.zeros(nConfs)
# s = np.zeros(nConfs)
# τ = np.zeros(nConfs)
# BB = np.zeros(nConfs)
# for j in range(nConfs):
#     with open(pFiles[j], 'r') as f:
#         exec(f.read())
#         BB[j] = B
#         ϵ[j] = ((1+2*B)*d**2*Us)/(36*ν*Ls) #small parameter
#         c[j] = Ls*ϵ[j]/Us * 2*(1-B)/(1+2*B) * (Γ/(2*np.pi))**2 * a**-4
#         Σc[j] = 2*κ/c[j]
#         s[j] = 4*κ
#         τ[j] = 1/(2*c[j])
# B = BB

def gen_op_cmap(c, α0=0):
    if isinstance(c, str):
        rgb = to_rgba(c)[:-1]
    else:
        rgb = c[0:3]
    arr = np.hstack([np.outer(np.ones(256), rgb), np.outer(np.linspace(α0,1,256), np.ones(1))])
    return ListedColormap(arr)
def gen_white_cmap(c, α0=0):
    if isinstance(c, str):
        rgb = to_rgba(c)[:-1]
    else:
        rgb = c[0:3]
    rgb = np.array(rgb)
    arr = np.hstack([np.linspace((1-α0)*(1-rgb)+rgb, rgb, 256), np.ones((256,1))])
    return ListedColormap(arr)
cHorz = gen_white_cmap('tab:blue')
cVert = gen_white_cmap('tab:green')
cTopo = gen_op_cmap('black')



'''-----------------------------------------------------------------------------
----------CODE------------------------------------------------------------------
-----------------------------------------------------------------------------'''
for k in ds.keys():
    ds[k] = ds[k].sel(Z=slice(0,-300), Zp1=slice(0, -300))

xData = [{} for j in range(nSnaps)]
yData = [{} for j in range(nSnaps)]
zData = [{} for j in range(nSnaps)]
for k in ds.keys():
    for j in range(nSnaps):
        xData[j][k] = gr[k].integrate(ds[k]['TRAC01'].interp(time=tSnaps[j]), 'X')
        yData[j][k] = gr[k].integrate(ds[k]['TRAC01'].interp(time=tSnaps[j]), 'Y')
        zData[j][k] = gr[k].integrate(ds[k]['TRAC01'].interp(time=tSnaps[j]), 'Z')

mHorz = np.maximum(xData[0][k].max().values, yData[0][k].max().values)
mVert = zData[0][k].max().values
lHorz = np.geomspace(mHorz*1e-3, mHorz, 20)
lVert = np.geomspace(mVert*1e-3, mVert, 20)
nHorz = LogNorm(mHorz*1e-3, mHorz)
nVert = LogNorm(mVert*1e-3, mVert)

for k in ds.keys():
    for j in range(nSnaps):
        xData[j][k] = np.maximum(xData[j][k], mHorz*1e-4)
        yData[j][k] = np.maximum(yData[j][k], mHorz*1e-4)
        zData[j][k] = np.maximum(zData[j][k], mVert*1e-4)

xLim = [ds[k]['XC'].min().values, ds[k]['XC'].max().values]
yLim = [ds[k]['YC'].min().values, ds[k]['YC'].max().values]
zLim = [ds[k]['Z'].min().values, ds[k]['Z'].max().values]



'''-----------------------------------------------------------------------------
----------PLOT------------------------------------------------------------------
-----------------------------------------------------------------------------'''
#tinker
bfs = 14
blw = 2

def initialize_plots():
    #declare variables
    global fC, aC, pC

    #declare plots
    fC = plt.figure(figsize=(13,8), layout='constrained')
    gC = fC.add_gridspec(nSnaps, nConfs)
    aC = np.array([np.array([fC.add_subplot(gC[i,j], projection='3d')
                    for j in range(nConfs)])
                    for i in range(nSnaps)])

    #label axes
    for a in aC[-1]:
        a.set_xlabel(r'longitude', fontsize=bfs)
        a.set_ylabel(r'latitude', fontsize=bfs)
    # for a in aC: a[0].set_zlabel(r'depth [m]', fontsize=bfs)

    #TODO: flip where vertical axis is labeled
        
    #prepare to store plots for legends
    pC = [[None for j in range(nConfs)] for i in range(nSnaps)]

def tidy_up_plots():
    #limits:
    for a in aC.flatten():
        a.set_xlim(xLim)
        a.set_ylim(yLim)
        a.set_zlim(zLim)

    #legends

    #save
    if saveFigures:
        today = np.datetime64('today').item()
        todayStr = '{0:02d}{1:02d}'.format(today.month, today.day)
        plt.figure(fC.number) 
        plt.savefig('../figures/'+todayStr+'_fig2.png')

if __name__ == "__main__":
    try:
        initialize_plots()

        for r in range(nSnaps):
            for k, c in zip(ds.keys(), range(nConfs)):
                vx = xData[r][k]
                vy = yData[r][k]
                vz = zData[r][k]
                aC[r,c].contourf(vx.transpose('Z',...), vx['YC'], vx['Z'],
                                 zdir='x', offset=xLim[1],
                                 levels=lHorz, norm=nHorz, extend='both', cmap=cHorz)
                aC[r,c].contourf(vy['XC'], vy.transpose('XC',...), vy['Z'],
                                 zdir='y', offset=yLim[0],
                                 levels=lHorz, norm=nHorz, extend='both', cmap=cHorz)
                aC[r,c].contourf(vz['XC'], vz['YC'], vz.transpose('YC',...),
                                 zdir='z', offset=zLim[1],
                                 levels=lVert, norm=nVert, extend='both', cmap=cVert)

                #stuff that will be moved to tidy_up_plots() tomorrow
                aC[r,c].plot([xLim[0], xLim[1]], [yLim[0], yLim[0]], [zLim[1], zLim[1]],
                             color='black', alpha=0.5, zorder=1e3)
                aC[r,c].plot([xLim[1], xLim[1]], [yLim[0], yLim[1]], [zLim[1], zLim[1]],
                             color='black', alpha=0.5, zorder=1e3)
                aC[r,c].plot([xLim[1], xLim[1]], [yLim[0], yLim[0]], [zLim[0], zLim[1]],
                             color='black', alpha=0.5, zorder=1e3)
                aC[r,c].plot_surface(ds[k]['XC']*xr.ones_like(ds[k]['YC']), xr.ones_like(ds[k]['XC'])*ds[k]['YC'], xr.zeros_like(ds[k]['RAC']).transpose('XC',...),
                                 facecolors=cTopo(xr.where(ds[k]['Depth'] > 0, 0., 1.).transpose('XC',...)))

        tidy_up_plots() 
        plt.show()
    finally:
        plt.close('all')

