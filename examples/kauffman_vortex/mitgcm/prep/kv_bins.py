#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 2021

@author: Mason Rogers
"""

#imports
import numpy as np
from scipy.stats import binom
import matplotlib.pyplot as plt
from ../../input/kv_param import *

#assemble grids
dx = .01; nx = 104
dy = .01; ny = 104
dz = .02; nz = 50
XG = np.arange(-nx*dx/2,nx*dx/2,dx); XC = XG + dx/2
YG = np.arange(-ny*dy/2,ny*dy/2,dy); YC = YG + dx/2
Zl = np.arange(dz*nz,0,-dz); Z = Zl - dz/2
uz, uy, ux = np.meshgrid(Z,YC,XG,indexing='ij')
vz, vy, vx = np.meshgrid(Z,YG,XC,indexing='ij')
wz, wy, wx = np.meshgrid(Zl,YC,XC,indexing='ij')
ur = np.sqrt(ux**2 + uy**2)
vr = np.sqrt(vx**2 + vy**2)
wr = np.sqrt(wx**2 + wy**2)

#rotating can parameters
R = (nx-4)*dx/2
a = .62
b = 7.5
c = .7
β = 1
γ = .2
δ = .2
ys = -.2
σ = 0
nSnaps = 10 if σ != 0 else 1

#Maxey-Riley Parameters
B = 0.98 #density ratio ρ/ρᶠ
ν = 1e-6 #fluid viscosity
d = 1e-3 #particle diameter
g = 9.8 #gravity
Us = 1 #characteristic velocity scale of flow
Ls = 1 #characteristic length scale of flow
ϵ = ((1+2*B)*d**2*Us)/(36*ν*Ls) #should be small
C = (2*g*(B-1)*Ls*ϵ)/((1+2*B)*Us**2) #should be O(1)

#initial probability parameters
vol = dx*dy*dz
Σ = (21/4) * dx**2
mx = np.round(4*Σ/dx**2)
my = np.round(4*Σ/dy**2)
mz = np.round(4*Σ/dz**2)

#make velocity binaries
T = 2*np.pi/σ if σ!= 0 else 1
for t, i in zip(np.linspace(0, T, nSnaps, endpoint=False), range(1,nSnaps+1)):
    #fluid velocities
    Uf = -b*ux*(1-2*uz)*(R-ur)/3 - a*uy*(c+uz**2) + δ*(1-β*uz)*(uy*(uy-ys+γ*np.cos(σ*t)) - (R**2-ur**2)/2)
    Vf = -b*vy*(1-2*vz)*(R-vr)/3 + a*vx*(c+vz**2) - δ*(1-β*vz)*vx*(vy-ys+γ*np.cos(σ*t))
    Wf = b*wz*(1-wz)*(2*R/3-wr)
    #time derivatives
    Ut = -δ*(1-β*uz)*uy*σ*γ*np.sin(σ*t)
    Vt = δ*(1-β*vz)*vx*σ*γ*np.sin(σ*t)
    Wt = np.zeros(Wf.shape)
    #x derivatives
    Ux = b*(1-2*uz)*ux**2/(3*ur) - b*(1-2*uz)*(R-ur)/3 + δ*(1-β*uz)*ux
    Vx = b*(1-2*vz)*vx*vy/(3*vr) + a*(c+vz**2) - δ*(1-β*vz)*(vy-ys+γ*np.cos(σ*t))
    Wx = -b*(1-wz)*wz/wr * wx
    #y derivatives
    Uy = b*(1-2*uz)*ux*uy/(3*ur) - a*(c+uz**2) + δ*(1-β*uz)*(3*uy-ys+γ*np.cos(σ*t))
    Vy = b*(1-2*vz)*vy**2/(3*vr) - b*(1-2*vz)*(R-vr)/3 - δ*(1-β*vz)*vx
    Wy = -b*(1-wz)*wz/wr * wy
    #z derivatives
    Uz = 2*b*ux*(R-ur)/3 - 2*a*uy*uz - δ*β*(uy*(uy-ys+γ*np.cos(σ*t)) - (R**2-ur**2)/2)
    Vz = 2*b*vy*(R-vr)/3 + 2*a*vx*vz + δ*β*vx*(vy-ys+γ*np.cos(σ*t))
    Wz = b*(1-2*wz)*(2*R/3-wr)
    #correction velocities
    M = Wf - C*Us #mean vertical velocity
    Uc = (ϵ*Ls/Us) * ((3/(1+2*B) - 1)*(Ut + Uf*Ux + Vf*Uy) + (3/(1+2*B)*Wf - M)*Uz)
    Vc = (ϵ*Ls/Us) * ((3/(1+2*B) - 1)*(Vt + Uf*Vx + Vf*Vy) + (3/(1+2*B)*Wf - M)*Vz)
    Wc = (ϵ*Ls/Us) * ((3/(1+2*B) - 1)*(Wt + Uf*Wx + Vf*Wy) + (3/(1+2*B)*Wf - M)*Wz)

    #Maxey-Riley velocities
    U = Uf + Uc
    V = Vf + Vc
    W = M + Wc

    #walls
    land = wr > R
    U[land] = 0
    V[land] = 0
    W[land] = 0

    #write flow to files
    with open("../input/input_off/B{0:.3f}_U.{1:010d}.data".format(B,i), "w") as ud:
        U.astype('>f4').tofile(ud)
    with open("../input/input_off/B{0:.3f}_V.{1:010d}.data".format(B,i), "w") as vd:
        V.astype('>f4').tofile(vd)
    with open("../input/input_off/B{0:.3f}_W.{1:010d}.data".format(B,i), "w") as wd:
        W.astype('>f4').tofile(wd)

    with open("../input/input_off/B{0:.3f}_U.{1:010d}.meta".format(B,i), "w") as um:
        um.write(" nDims = [3];\n")
        um.write(" dimList = [\n {0:d}, 1, {0:d}, \
                              \n {1:d}, 1, {1:d}, \
                              \n {2:d}, 1, {2:d} \
                              \n ];\n".format(nx, ny, nz))
        um.write(" dataprec = ['float32'];\n")
        um.write(" nrecords = [1];\n")
        um.write(" timestepnumber = [1];\n")
    with open("../input/input_off/B{0:.3f}_V.{1:010d}.meta".format(B,i), "w") as vm:
        vm.write(" nDims = [3];\n")
        vm.write(" dimList = [\n {0:d}, 1, {0:d}, \
                              \n {1:d}, 1, {1:d}, \
                              \n {2:d}, 1, {2:d} \
                              \n ];\n".format(nx, ny, nz))
        vm.write(" dataprec = ['float32'];\n")
        vm.write(" nrecords = [1];\n")
        vm.write(" timestepnumber = [1];\n")
    with open("../input/input_off/B{0:.3f}_W.{1:010d}.meta".format(B,i), "w") as wm:
        wm.write(" nDims = [3];\n")
        wm.write(" dimList = [\n {0:d}, 1, {0:d}, \
                              \n {1:d}, 1, {1:d}, \
                              \n {2:d}, 1, {2:d} \
                              \n ];\n".format(nx, ny, nz))
        wm.write(" dataprec = ['float32'];\n")
        wm.write(" nrecords = [1];\n")
        wm.write(" timestepnumber = [1];\n")

#write walls to file
with open("../input/top.bin", "w") as top:
    d = (land - 1)[0,:,:] * dz*nz
    d.astype('>f4').tofile(top)

#write initial conditions
with open("../input/w_p_init1.bin", "w") as p_init:
    x0 = 0.33
    y0 = 0.00
    z0 = 0.56
    px = binom.pmf(np.floor((XC-x0)/dx),mx,0.5,loc=-mx//2)
    py = binom.pmf(np.floor((YC-y0)/dy),my,0.5,loc=-my//2)
    pz = binom.pmf(np.floor((Z-z0)/dz),mz,0.5,loc=-mz//2)
    p = np.prod(np.stack(np.meshgrid(pz,py,px,indexing='ij')),0) / vol
    #p = (np.sqrt((wx-x0)**2 + (wy-y0)**2 + (uz-z0)**2) < .05)# * (np.abs(wy-y0) < (dy/4))
    #p = p/np.sum(p)/vol 
    p.astype('>f4').tofile(p_init)
with open("../input/w_p_init2.bin", "w") as p_init:
    x0 = 0.41
    y0 = 0.00
    z0 = 0.17
    px = binom.pmf(np.floor((XC-x0)/dx),mx,0.5,loc=-mx//2)
    py = binom.pmf(np.floor((YC-y0)/dy),my,0.5,loc=-my//2)
    pz = binom.pmf(np.floor((Z-z0)/dz),mz,0.5,loc=-mz//2)
    p = np.prod(np.stack(np.meshgrid(pz,py,px,indexing='ij')),0) / vol
    #p = (np.sqrt((wx-x0)**2 + (wy-y0)**2 + (uz-z0)**2) < .05)# * (np.abs(wy-y0) < (dy/4))
    #p = p/np.sum(p)/vol 
    p.astype('>f4').tofile(p_init)
with open("../input/w_p_init3.bin", "w") as p_init:
    x0 = 0.29
    y0 = 0.00
    z0 = 0.89
    px = binom.pmf(np.floor((XC-x0)/dx),mx,0.5,loc=-mx//2)
    py = binom.pmf(np.floor((YC-y0)/dy),my,0.5,loc=-my//2)
    pz = binom.pmf(np.floor((Z-z0)/dz),mz,0.5,loc=-mz//2)
    p = np.prod(np.stack(np.meshgrid(pz,py,px,indexing='ij')),0) / vol
    #p = (np.sqrt((wx-x0)**2 + (wy-y0)**2 + (uz-z0)**2) < .05)# * (np.abs(wy-y0) < (dy/4))
    #p = p/np.sum(p)/vol 
    p.astype('>f4').tofile(p_init)
with open("../input/w_p_init4.bin", "w") as p_init: 
    x0 = 0.07
    y0 = 0.00
    z0 = 0.12
    px = binom.pmf(np.floor((XC-x0)/dx),mx,0.5,loc=-mx//2)
    py = binom.pmf(np.floor((YC-y0)/dy),my,0.5,loc=-my//2)
    pz = binom.pmf(np.floor((Z-z0)/dz),mz,0.5,loc=-mz//2)
    p = np.prod(np.stack(np.meshgrid(pz,py,px,indexing='ij')),0) / vol
    #p = (np.sqrt((wx-x0)**2 + (wy-y0)**2 + (uz-z0)**2) < .05)# * (np.abs(wy-y0) < (dy/4))
    #p = p/np.sum(p)/vol 
    p.astype('>f4').tofile(p_init)
