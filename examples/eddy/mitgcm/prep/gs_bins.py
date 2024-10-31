#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 2021

@author: Mason Rogers
"""

#imports
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy.stats import binom, norm
from gs_param import *

#parameters
ϵ = ((1+2*B)*d**2*Us)/(36*ν*Ls)

#assemble grids
gr = xr.open_dataset('../gulf_stream/gs_grid.nc')

#depth and walls
d = -gr['Depth']
isWet = (gr['Zl'] > d).any(dim='k_l')
ob = {}
ob['W'] = isWet.isel(i=2)
ob['E'] = isWet.isel(i=-3)
ob['N'] = isWet.isel(j=-3)
ob['S'] = isWet.isel(j=2)
card = ob.keys()

#get rid of Chesapeake creeks that frustrate OBCS
for j in range(ob['W'].size):
    if not(ob['W'][j]):
        d[j,:3] = 0

#make OBCS strings
l = {}
for k in card:
    o = ob[k].values[0]
    l[k] = [0]
    for x in ob[k]:
        if (x and o) or (not(x) and not(o)):
            l[k][-1] += 1
        else:
            l[k].append(1)
            o = not(o)
    o = ob[k].values[0]
    v = 3 if k in ['W','S'] else -3 
    for n in range(len(l[k])):
        l[k][n] = str(l[k][n]) + '*' + str(v if o else 0)
        o = not(o)
    l[k] = ','.join(l[k]) + ','    
print(l)

#choose initial conditions
isDry = d > gr['Zl'][5]
hasDryNbr = isDry.shift(i=1, j=1, fill_value=False) + \
            isDry.shift(i=1, j=0, fill_value=False) + \
            isDry.shift(i=1, j=-1, fill_value=False) + \
            isDry.shift(i=0, j=-1, fill_value=False) + \
            isDry.shift(i=-1, j=-1, fill_value=False) + \
            isDry.shift(i=-1, j=0, fill_value=False) + \
            isDry.shift(i=-1, j=1, fill_value=False) + \
            isDry.shift(i=0, j=1, fill_value=False)
isCoastal = np.logical_not(isDry) & hasDryNbr 
p = np.zeros((87,384,216))
#p[5][isCoastal] = 1.
#hockey puck p[5,200:220,80:100] = 1.
k, j, i = np.meshgrid(np.arange(87), np.arange(384), np.arange(216), indexing='ij')
p = np.exp(-((i-172)**2 + (j-243)**2)/25 - (k-25)**2/9)  #eddy center

#write walls to file
with open("../run/top.bin", "w") as top:
    d.values.astype('>f4').tofile(top)

#write initial conditions
with open("/pool001/masonr/eddy/p_230712_blob.bin", "w") as p_init:
    #p = np.zeros((87,384,216))
    #p[5,190:200,120:130] = 100
    p.astype('>f4').tofile(p_init)
