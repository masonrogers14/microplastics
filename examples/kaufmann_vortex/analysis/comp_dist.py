#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 30 2021

@author: Mason Rogers
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

#select relevant variables
p1 = ds['TRAC01']
p2 = ds['TRAC02']

#get volume
vol = ds['rA'] * ds['drF'] * ds['maskC']

#total variation distance
δ = (vol*(p1-p2)).where(p1>p2).sum(dim=['XC','YC','Z'])

#Hellinger distance
H = da.sqrt(1-(vol*da.sqrt(p1*p2)).sum(dim=['XC','YC','Z']))

#average percentage error
η = (da.fabs(p1-p2)/da.maximum(da.maximum(p1,p2), 1e-4/vol)).mean(dim=['XC','YC','Z'])
