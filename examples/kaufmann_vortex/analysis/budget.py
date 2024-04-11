#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 31 2021

@author: Mason Rogers
"""

#tinker
nTracs = 1

#imports
import numpy as np
import xarray as xr
from read_MITgcm import ds, gr

#compute budget terms
hb_vars = {}
vol = ds['rA']*ds['drF']

for n in range(1,nTracs+1):
    hb_vars['pt'+str(n)] = ds['Tp_gTr0'+str(n)]
    hb_vars['ax'+str(n)] = gr.diff(ds['ADVxTr0'+str(n)], 'X', boundary='extend') / vol
    hb_vars['ay'+str(n)] = gr.diff(ds['ADVyTr0'+str(n)], 'Y', boundary='extend') / vol
    #hb_vars['az'+str(n)] = -gr.diff(ds['ADVrTr0'+str(n)], 'Z', boundary='extend') / vol
    hb_vars['dx'+str(n)] = gr.diff(ds['DFxETr0'+str(n)], 'X', boundary='extend') / vol
    hb_vars['dy'+str(n)] = gr.diff(ds['DFyETr0'+str(n)], 'Y', boundary='extend') / vol
    #hb_vars['de'+str(n)] = -gr.diff(ds['DFrETr0'+str(n)], 'Z', boundary='extend') / vol
    #hb_vars['di'+str(n)] = -gr.diff(ds['DFrITr0'+str(n)], 'Z', boundary='extend') / vol 
    hb_vars['ad'+str(n)] = hb_vars['ax'+str(n)] + hb_vars['ay'+str(n)]# + hb_vars['az'+str(n)]
    hb_vars['df'+str(n)] = hb_vars['dx'+str(n)] + hb_vars['dy'+str(n)]# + hb_vars['de'+str(n)] + hb_vars['di'+str(n)]
    hb_vars['rs'+str(n)] = hb_vars['pt'+str(n)] + hb_vars['ad'+str(n)] + hb_vars['df'+str(n)]
    
    abs_pt = np.abs(hb_vars['pt'+str(n)])
    abs_ad = np.abs(hb_vars['ad'+str(n)])
    abs_df = np.abs(hb_vars['df'+str(n)])
    norm = abs_pt.where(abs_pt > abs_ad, abs_ad)
    norm = norm.where(norm > abs_df, abs_df)
    norm = norm.where(norm > 0, np.inf)
    hb_vars['ne'+str(n)] = hb_vars['rs'+str(n)] / norm
    
hb = xr.Dataset(data_vars = hb_vars)

