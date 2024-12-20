#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
kv_figc.py plots the variance growth and snapshots for two configurations of the
Kaufmann vortex experiment.

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
labels = ['full MC', 'reduced MC', 'Dedalus', 'MITgcm']

#imports
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from dict_MITgcm import ds, gr, dirs



'''-----------------------------------------------------------------------------
----------CODE------------------------------------------------------------------
-----------------------------------------------------------------------------'''
def tvd(p1, p2, gr):
    return gr.integrate((p1 - p2).where(p1 > p2), ['X', 'Y', 'Z'])

M = [[[None for _ in range(4)] for _ in range(4)] for _ in range(2)]
for i in range(4):
    for j in range(i):
        for k, di in enumerate(dirs):
            d = ds[di]
            g = gr[di]
            k1 = 'TRAC{0:02d}'.format(i+1)
            k2 = 'TRAC{0:02d}'.format(j+1)
            M[k][j][i] = tvd(d[k1], d[k2], g)

#read files
ϵ = np.zeros(2)
for j in range(2):
    with open(pFiles[j], 'r') as f:
        exec(f.read())
        ϵ[j] = ((1+2*B)*d**2*Us)/(36*ν*Ls) #small parameter




'''-----------------------------------------------------------------------------
----------PLOT------------------------------------------------------------------
-----------------------------------------------------------------------------'''
plt.style.use('mason')

fig, ax = plt.subplots(ncols=2, sharey=True, sharex=True,
                       figsize=(10, 4), layout='constrained')
ax[0].set_xlabel('time')
ax[1].set_xlabel('time')
ax[0].set_ylabel(r'$\delta(p_a, p_e)$')

colors = [plt.get_cmap('tab10')(i) for i in range(4)]

for k, di in enumerate(dirs):
    t = ds[di]['time']
    for i in range(1, 4):
        ax[k].plot(t, M[k][0][i],
                   label=labels[i],
                   color=colors[i])
    # ax[k].plot(t, ϵ[k]*t, '--', color='grey', lw=1, label=r'$\varepsilon t$')

ax[1].legend(title=r'$p_a$')

plt.savefig('../figures/kv_tvd.png', dpi=200, transparent=False)
