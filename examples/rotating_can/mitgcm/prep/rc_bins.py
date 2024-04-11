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
from rc_param import *

#assemble grids
nx = int(np.round(2*R/dx)) + 4
ny = int(np.round(2*R/dy)) + 4
nz = int(np.round(H/dz))
XG = np.arange(-nx*dx/2,(nx-1)*dx/2,dx); XC = XG + dx/2
YG = np.arange(-ny*dy/2,(ny-1)*dy/2,dy); YC = YG + dx/2
Zl = np.arange(dz*nz,0,-dz); Z = Zl - dz/2
uz, uy, ux = np.meshgrid(Z,YC,XG,indexing='ij')
vz, vy, vx = np.meshgrid(Z,YG,XC,indexing='ij')
wz, wy, wx = np.meshgrid(Zl,YC,XC,indexing='ij')
ur = np.sqrt(ux**2 + uy**2)
vr = np.sqrt(vx**2 + vy**2)
wr = np.sqrt(wx**2 + wy**2)
vol = dx*dy*dz

#make velocity binaries
nSnaps = 10 if σ > 0 else 1
T = 2*np.pi/σ if σ != 0 else 1
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
    with open("../run/off/B{0:.3f}_U.{1:010d}.data".format(B,i), "w") as ud:
        U.astype('>f4').tofile(ud)
    with open("../run/off/B{0:.3f}_V.{1:010d}.data".format(B,i), "w") as vd:
        V.astype('>f4').tofile(vd)
    with open("../run/off/B{0:.3f}_W.{1:010d}.data".format(B,i), "w") as wd:
        W.astype('>f4').tofile(wd)

    with open("../run/off/B{0:.3f}_U.{1:010d}.meta".format(B,i), "w") as um:
        um.write(" nDims = [3];\n")
        um.write(" dimList = [\n {0:d}, 1, {0:d}, \
                              \n {1:d}, 1, {1:d}, \
                              \n {2:d}, 1, {2:d} \
                              \n ];\n".format(nx, ny, nz))
        um.write(" dataprec = ['float32'];\n")
        um.write(" nrecords = [1];\n")
        um.write(" timestepnumber = [1];\n")
    with open("../run/off/B{0:.3f}_V.{1:010d}.meta".format(B,i), "w") as vm:
        vm.write(" nDims = [3];\n")
        vm.write(" dimList = [\n {0:d}, 1, {0:d}, \
                              \n {1:d}, 1, {1:d}, \
                              \n {2:d}, 1, {2:d} \
                              \n ];\n".format(nx, ny, nz))
        vm.write(" dataprec = ['float32'];\n")
        vm.write(" nrecords = [1];\n")
        vm.write(" timestepnumber = [1];\n")
    with open("../run/off/B{0:.3f}_W.{1:010d}.meta".format(B,i), "w") as wm:
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

#write initial conditions
with open("../run/p_init1.bin", "w") as p_init:
    px = norm.pdf(XC, loc=x0, scale=np.cbrt(Σ)) * dx
    py = norm.pdf(YC, loc=y0, scale=np.cbrt(Σ)) * dz
    pz = norm.pdf(Z, loc=z0, scale=np.cbrt(Σ)) * dy 
    p = np.prod(np.stack(np.meshgrid(pz,py,px,indexing='ij')),0) / vol
    p.astype('>f4').tofile(p_init)
with open("../run/p_init2.bin", "w") as p_init:
    x0=.41; y0=0.0; z0=.17;
    px = norm.pdf(XC, loc=x0, scale=np.cbrt(Σ)) * dx
    py = norm.pdf(YC, loc=y0, scale=np.cbrt(Σ)) * dz
    pz = norm.pdf(Z, loc=z0, scale=np.cbrt(Σ)) * dy 
    p = np.prod(np.stack(np.meshgrid(pz,py,px,indexing='ij')),0) / vol
    p.astype('>f4').tofile(p_init)

