#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 2021

@author: Mason Rogers
"""

#imports
import numpy as np
from scipy.stats import binom, norm
import matplotlib.pyplot as plt
from kv_param import *

#parameters
ϵ = ((1+2*B)*d**2*Us)/(36*ν*Ls)

#assemble grids
nx = int(np.round(2*R/dx)) + 4
ny = int(np.round(2*R/dy)) + 4
nz = 1
XG = np.arange(-nx*dx/2,(nx-1)*dx/2,dx); XC = XG + dx/2
YG = np.arange(-ny*dy/2,(ny-1)*dy/2,dy); YC = YG + dx/2
Zl = np.arange(dz*nz,0,-dz); Z = Zl - dz/2
uz, uy, ux = np.meshgrid(Z,YC,XG,indexing='ij')
vz, vy, vx = np.meshgrid(Z,YG,XC,indexing='ij')
wz, wy, wx = np.meshgrid(Zl,YC,XC,indexing='ij')
ur = np.sqrt(ux**2 + uy**2)
vr = np.sqrt(vx**2 + vy**2)
wr = np.sqrt(wx**2 + wy**2)

#initial probability parameters
vol = dx*dy*dz
mx = np.round(4*Σ/dx**2)
my = np.round(4*Σ/dy**2)
mz = np.round(4*Σ/dz**2)

#make velocity binaries
#fluid velocities
Uf = -Γ/2/np.pi * uy/(a**2 + ur**2)
Vf = Γ/2/np.pi * vx/(a**2 + vr**2)
#x derivatives
Ux = 0
Vx = Γ/2/np.pi/(a**2 + vr**2)
#y derivatives
Uy = -Γ/2/np.pi/(a**2 + ur**2)
Vy = 0
#correction velocities
Uc = (ϵ*Ls/Us) * (3/(1+2*B) - 1)*(Uf*Ux + Vf*Uy)
Vc = (ϵ*Ls/Us) * (3/(1+2*B) - 1)*(Uf*Vx + Vf*Vy)
#Maxey-Riley velocities
U = Uf + Uc
V = Vf + Vc
W = 0 * wr

#walls
land = wr > R
U[land] = 0
V[land] = 0
W[land] = 0

#write flow to files
with open("../run/off/B{0:.3f}_U.{1:010d}.data".format(B,1), "w") as ud:
    U.astype('>f4').tofile(ud)
with open("../run/off/B{0:.3f}_V.{1:010d}.data".format(B,1), "w") as vd:
    V.astype('>f4').tofile(vd)
with open("../run/off/B{0:.3f}_W.{1:010d}.data".format(B,1), "w") as wd:
    W.astype('>f4').tofile(wd)

with open("../run/off/B{0:.3f}_U.{1:010d}.meta".format(B,1), "w") as um:
    um.write(" nDims = [3];\n")
    um.write(" dimList = [\n {0:d}, 1, {0:d}, \
                          \n {1:d}, 1, {1:d}, \
                          \n {2:d}, 1, {2:d} \
                          \n ];\n".format(nx, ny, nz))
    um.write(" dataprec = ['float32'];\n")
    um.write(" nrecords = [1];\n")
    um.write(" timestepnumber = [1];\n")
with open("../run/off/B{0:.3f}_V.{1:010d}.meta".format(B,1), "w") as vm:
    vm.write(" nDims = [3];\n")
    vm.write(" dimList = [\n {0:d}, 1, {0:d}, \
                          \n {1:d}, 1, {1:d}, \
                          \n {2:d}, 1, {2:d} \
                          \n ];\n".format(nx, ny, nz))
    vm.write(" dataprec = ['float32'];\n")
    vm.write(" nrecords = [1];\n")
    vm.write(" timestepnumber = [1];\n")
with open("../run/off/B{0:.3f}_W.{1:010d}.meta".format(B,1), "w") as wm:
    wm.write(" nDims = [3];\n")
    wm.write(" dimList = [\n {0:d}, 1, {0:d}, \
                          \n {1:d}, 1, {1:d}, \
                          \n {2:d}, 1, {2:d} \
                          \n ];\n".format(nx, ny, nz))
    wm.write(" dataprec = ['float32'];\n")
    wm.write(" nrecords = [1];\n")
    wm.write(" timestepnumber = [1];\n")

#write walls to file
with open("../run/top.bin", "w") as top:
    d = (land - 1)[0,:,:] * dz*nz
    d.astype('>f4').tofile(top)

#write initial conditions (now supplied in kv_param.py)
with open("../run/p_init.bin", "w") as p_init:
    #px = binom.pmf(np.round((XC-x0)/dx),mx,0.5,loc=-mx//2)
    #py = binom.pmf(np.round((YC-y0)/dy),my,0.5,loc=-my//2)
    px = norm.pdf(XC, loc=x0, scale=np.sqrt(Σ)) * dx
    py = norm.pdf(YC, loc=y0, scale=np.sqrt(Σ)) * dy 
    pz = np.ones(1)
    p = np.prod(np.stack(np.meshgrid(pz,py,px,indexing='ij')),0) / vol
    #p = (np.sqrt((wx-x0)**2 + (wy-y0)**2 + (uz-z0)**2) < .05)# * (np.abs(wy-y0) < (dy/4))
    #p = p/np.sum(p)/vol
    p[land] = 0
    p.astype('>f4').tofile(p_init)
