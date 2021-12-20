#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 17 2021

@author: Mason Rogers

moments.py computes spatial moments of probability distributions.
"""

#tinker

#imports
import numpy as np
import xarray as xr
import dask.array as da
import matplotlib.pyplot as plt
import matplotlib.animation as movie
from read_MITgcm import ds
from matplotlib.colors import LogNorm, to_rgba, ListedColormap

#get volume and variable keys
vol = ds['rA'] * ds['drF'] * ds['maskC']
keys = list(ds.data_vars.keys())
sDims = ['XC','YC','Z']

#compute moments
moms = {}
for key in keys:
    p = ds[key]
    P = {'P': (p*vol).sum(dim=sDims)} 
    μ = {'μ_'+d: (p*vol*ds[d]).sum(dim=sDims) for d in sDims}
    Σ = {'Σ_'+d1+'_'+d2: (p*vol*(ds[d1]-μ['μ_'+d1])*(ds[d2]-μ['μ_'+d2])).sum(dim=sDims) for d1 in sDims for d2 in sDims}
    moms[key] = xr.Dataset(data_vars={**P, **μ, **Σ})

