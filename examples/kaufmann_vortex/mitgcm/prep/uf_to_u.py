#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 2021

@author: Mason Rogers

uf_to_u.py is designed to read MITgcm fluid velocity files and output
Maxey-Riley advection velocity files.
Code limitations:
1.  Grid parameters are not inferred from model output and must be set manually.
2.  Data are assumed to be periodic in time and space, but the code attempts to
    infer boundaries from topographic data.
3.  Only constant-spacing grids with flat bottoms are supported for now.
4.  Only nonrotating, unstratified fluids are supported for now.
"""

#tinker
fourth_order = False

#imports
import numpy as np

#grid parameters
nt = 1
dt = 2*np.pi/(3*nt)
dx = .01
dy = .01
dz = .02
nx = 104
ny = 104
nz = 50

#Maxey-Riley parameters
B = 1.04 #density ratio ρ/ρᶠ
ν = 1e-6 #fluid viscosity
d = 1e-3 #particle diameter
g = 9.8 #gravity
Us = 1 #characteristic velocity scale of flow
Ls = 1 #characteristic length scale of flow
ϵ = ((1+2*B)*d**2*Us)/(36*ν*Ls) #should be small
C = (2*g*(B-1)*Ls*ϵ)/((1+2*B)*Us**2) #should be O(1)

#filenames
top_file = "../input/top.bin"
Uf_files = ["../input/input_off/B{0:.3f}_U.{1:010d}.data".format(1,i) for i in range(1, nt+1)]
Vf_files = ["../input/input_off/B{0:.3f}_V.{1:010d}.data".format(1,i) for i in range(1, nt+1)]
Wf_files = ["../input/input_off/B{0:.3f}_W.{1:010d}.data".format(1,i) for i in range(1, nt+1)]
Ua_files = ["../input/input_off/B{0:.3f}_U.{1:010d}.data".format(B,i) for i in range(1, nt+1)]
Va_files = ["../input/input_off/B{0:.3f}_V.{1:010d}.data".format(B,i) for i in range(1, nt+1)]
Wa_files = ["../input/input_off/B{0:.3f}_W.{1:010d}.data".format(B,i) for i in range(1, nt+1)]

#read topography data
with open(top_file, "r") as f:
    mask = np.fromfile(f, '>f4').reshape((ny,nx)) != 0

#identify x-axis distances to boundary for points near boundary
wb = np.zeros((ny,nx))
eb = np.zeros((ny,nx))
wb += (2*(1-np.roll(mask, 2, axis=1))*np.roll(mask, 1, axis=1)*mask)
wb += ((1-np.roll(mask, 1, axis=1))*mask)
eb += (2*(1-np.roll(mask, -2, axis=1))*np.roll(mask, -1, axis=1)*mask)
eb += ((1-np.roll(mask, -1, axis=1))*mask)

#identify y-axis distances to boundary for points near boundary
sb = np.zeros((ny,nx))
nb = np.zeros((ny,nx))
sb += (2*(1-np.roll(mask, 2, axis=0))*np.roll(mask, 1, axis=0)*mask)
sb += ((1-np.roll(mask, 1, axis=0))*mask)
nb += (2*(1-np.roll(mask, -2, axis=0))*np.roll(mask, -1, axis=0)*mask)
nb += ((1-np.roll(mask, -1, axis=0))*mask)

#prevent nonsensical derivatives
if ((eb==1)*(wb==1) + (nb==1)*(sb==1)).any():
    raise ValueError("Spatial derivative calculation failed (probably due to wet point bounded by coasts on both sides)")
xc4 = mask * (wb==0) * (eb==0)
xc2 = mask * (wb!=1) * (eb!=1) * (1-xc4)
xf2 = mask * (wb==1) * (eb==0)
xf1 = mask * (wb==1) * (eb==2)
xb2 = mask * (wb==0) * (eb==1)
xb1 = mask * (wb==2) * (eb==1)
yc4 = mask * (sb==0) * (nb==0)
yc2 = mask * (sb!=1) * (nb!=1) * (1-yc4)
yf2 = mask * (sb==1) * (nb==0)
yf1 = mask * (sb==1) * (nb==2)
yb2 = mask * (sb==0) * (nb==1)
yb1 = mask * (sb==2) * (nb==1)
zc4 = np.hstack([0,0,np.ones(nz-4),0,0])
zc2 = np.hstack([0,1,np.zeros(nz-4),1,0])
zf2 = np.hstack([np.zeros(nz-1),1])
zb2 = np.hstack([1,np.zeros(nz-1)])
if ~fourth_order:
    xc2 += xc4; xc4 ^= xc4
    yc2 += yc4; yc4 ^= yc4
    zc2 += zc4; zc4 -= zc4
if (((1-mask) + xc4 + xc2 + xf2 + xf1 + xb2 + xb1) != 1).any():
    raise ValueError("Could not uniquely determine x-differencing method for each point.")
if (((1-mask) + yc4 + yc2 + yf2 + yf1 + yb2 + yb1) != 1).any():
    raise ValueError("Could not uniquely determine y-differencing method for each point.")
if ((zc4 + zc2 + zf2 + zb2) != 1).any():
    raise ValueError("Could not uniquely determine z-differencing method for each point.")

#read 4D data
Uf = np.zeros((nz,ny,nx,nt))
Vf = np.zeros((nz,ny,nx,nt))
Wf = np.zeros((nz,ny,nx,nt))
for i in range(nt):
    with open(Uf_files[i], "r") as Uf_file:
        Uf[:,:,:,i] = np.fromfile(Uf_file, '>f4').reshape((nz,ny,nx))
    with open(Vf_files[i], "r") as Vf_file:
        Vf[:,:,:,i] = np.fromfile(Vf_file, '>f4').reshape((nz,ny,nx))
    with open(Wf_files[i], "r") as Wf_file:
        Wf[:,:,:,i] = np.fromfile(Wf_file, '>f4').reshape((nz,ny,nx))

#time derivatives
Ut = np.zeros((nz,ny,nx,nt))
Vt = np.zeros((nz,ny,nx,nt))
Wt = np.zeros((nz,ny,nx,nt))
for i in range(nt):
    imm = np.mod(i-2, nt)
    im = np.mod(i-1, nt)
    ip = np.mod(i+1, nt)
    ipp = np.mod(i+2, nt)
    #centered fourth-order first derivative
    Ut[:,:,:,i] = (Uf[:,:,:,imm]/12 - 2*Uf[:,:,:,im]/3 + 2*Uf[:,:,:,ip]/3 - Uf[:,:,:,ipp]/12) / dt
    Vt[:,:,:,i] = (Vf[:,:,:,imm]/12 - 2*Vf[:,:,:,im]/3 + 2*Vf[:,:,:,ip]/3 - Vf[:,:,:,ipp]/12) / dt
    Wt[:,:,:,i] = (Wf[:,:,:,imm]/12 - 2*Wf[:,:,:,im]/3 + 2*Wf[:,:,:,ip]/3 - Wf[:,:,:,ipp]/12) / dt

#horizontal derivatives (assumes periodic domain and tries to handle walls like MITgcm)
Ux = np.zeros((nz,ny,nx,nt))
Vx = np.zeros((nz,ny,nx,nt))
Wx = np.zeros((nz,ny,nx,nt))
Uy = np.zeros((nz,ny,nx,nt))
Vy = np.zeros((nz,ny,nx,nt))
Wy = np.zeros((nz,ny,nx,nt))
for i in range(nx):
    for j in range(ny):
        imm = np.mod(i-2,nx); jmm = np.mod(j-2,ny)
        im = np.mod(i-1,nx); jm = np.mod(j-1,ny)
        ip = np.mod(i+1,nx); jp = np.mod(j+1,ny)
        ipp = np.mod(i+2,nx); jpp = np.mod(j+2,ny)
        if xc4[j,i]: #centered fourth-order first derivative
            Ux[:,j,i,:] = (Uf[:,j,imm,:]/12 - 2*Uf[:,j,im,:]/3 + 2*Uf[:,j,ip,:]/3 - Uf[:,j,ipp,:]/12) / dx
            Vx[:,j,i,:] = (Vf[:,j,imm,:]/12 - 2*Vf[:,j,im,:]/3 + 2*Vf[:,j,ip,:]/3 - Vf[:,j,ipp,:]/12) / dx
            Wx[:,j,i,:] = (Wf[:,j,imm,:]/12 - 2*Wf[:,j,im,:]/3 + 2*Wf[:,j,ip,:]/3 - Wf[:,j,ipp,:]/12) / dx
        elif xc2[j,i]: #centered second-order first derivative
            Ux[:,j,i,:] = (-Uf[:,j,im,:]/2 + Uf[:,j,ip,:]/2) / dx
            Vx[:,j,i,:] = (-Vf[:,j,im,:]/2 + Vf[:,j,ip,:]/2) / dx
            Wx[:,j,i,:] = (-Wf[:,j,im,:]/2 + Wf[:,j,ip,:]/2) / dx
        elif xf2[j,i]: #forward second-order first derivative
            Ux[:,j,i,:] = (-3*Uf[:,j,i,:]/2 + 2*Uf[:,j,ip,:] - Uf[:,j,ipp,:]/2) / dx
            Vx[:,j,i,:] = (-3*Vf[:,j,i,:]/2 + 2*Vf[:,j,ip,:] - Vf[:,j,ipp,:]/2) / dx
            Wx[:,j,i,:] = (-3*Wf[:,j,i,:]/2 + 2*Wf[:,j,ip,:] - Wf[:,j,ipp,:]/2) / dx
        elif xf1[j,i]: #forward first-order first derivative
            Ux[:,j,i,:] = (-Uf[:,j,i,:] + Uf[:,j,ip,:]) / dx
            Vx[:,j,i,:] = (-Vf[:,j,i,:] + Vf[:,j,ip,:]) / dx
            Wx[:,j,i,:] = (-Wf[:,j,i,:] + Wf[:,j,ip,:]) / dx
        elif xb2[j,i]: #backward second-order first derivative
            Ux[:,j,i,:] = (3*Uf[:,j,i,:]/2 - 2*Uf[:,j,im,:] + Uf[:,j,imm,:]/2) / dx
            Vx[:,j,i,:] = (3*Vf[:,j,i,:]/2 - 2*Vf[:,j,im,:] + Vf[:,j,imm,:]/2) / dx
            Wx[:,j,i,:] = (3*Wf[:,j,i,:]/2 - 2*Wf[:,j,im,:] + Wf[:,j,imm,:]/2) / dx
        elif xb1[j,i]: #backward first-order first derivative
            Ux[:,j,i,:] = (Uf[:,j,i,:] - Uf[:,j,im,:]) / dx
            Vx[:,j,i,:] = (Vf[:,j,i,:] - Vf[:,j,im,:]) / dx
            Wx[:,j,i,:] = (Wf[:,j,i,:] - Wf[:,j,im,:]) / dx
        else:
            assert mask[j,i] == 0
        if yc4[j,i]: #centered fourth-order first derivative
            Uy[:,j,i,:] = (Uf[:,jmm,i,:]/12 - 2*Uf[:,jm,i,:]/3 + 2*Uf[:,jp,i,:]/3 - Uf[:,jpp,i,:]/12) / dy
            Vy[:,j,i,:] = (Vf[:,jmm,i,:]/12 - 2*Vf[:,jm,i,:]/3 + 2*Vf[:,jp,i,:]/3 - Vf[:,jpp,i,:]/12) / dy
            Wy[:,j,i,:] = (Wf[:,jmm,i,:]/12 - 2*Wf[:,jm,i,:]/3 + 2*Wf[:,jp,i,:]/3 - Wf[:,jpp,i,:]/12) / dy
        elif yc2[j,i]: #centered second-order first derivative
            Uy[:,j,i,:] = (-Uf[:,jm,i,:]/2 + Uf[:,jp,i,:]/2) / dy
            Vy[:,j,i,:] = (-Vf[:,jm,i,:]/2 + Vf[:,jp,i,:]/2) / dy
            Wy[:,j,i,:] = (-Wf[:,jm,i,:]/2 + Wf[:,jp,i,:]/2) / dy
        elif yf2[j,i]: #forward second-order first derivative
            Uy[:,j,i,:] = (-3*Uf[:,j,i,:]/2 + 2*Uf[:,jp,i,:] - Uf[:,jpp,i,:]/2) / dy
            Vy[:,j,i,:] = (-3*Vf[:,j,i,:]/2 + 2*Vf[:,jp,i,:] - Vf[:,jpp,i,:]/2) / dy
            Wy[:,j,i,:] = (-3*Wf[:,j,i,:]/2 + 2*Wf[:,jp,i,:] - Wf[:,jpp,i,:]/2) / dy
        elif yf1[j,i]: #forward first-order first derivative
            Uy[:,j,i,:] = (-Uf[:,j,i,:] + Uf[:,jp,i,:]) / dy
            Vy[:,j,i,:] = (-Vf[:,j,i,:] + Vf[:,jp,i,:]) / dy
            Wy[:,j,i,:] = (-Wf[:,j,i,:] + Wf[:,jp,i,:]) / dy
        elif yb2[j,i]: #backward second-order first derivative
            Uy[:,j,i,:] = (3*Uf[:,j,i,:]/2 - 2*Uf[:,jm,i,:] + Uf[:,jmm,i,:]/2) / dy
            Vy[:,j,i,:] = (3*Vf[:,j,i,:]/2 - 2*Vf[:,jm,i,:] + Vf[:,jmm,i,:]/2) / dy
            Wy[:,j,i,:] = (3*Wf[:,j,i,:]/2 - 2*Wf[:,jm,i,:] + Wf[:,jmm,i,:]/2) / dy
        elif yb1[j,i]: #backward first-order first derivative
            Uy[:,j,i,:] = (Uf[:,j,i,:] - Uf[:,jm,i,:]) / dy
            Vy[:,j,i,:] = (Vf[:,j,i,:] - Vf[:,jm,i,:]) / dy
            Wy[:,j,i,:] = (Wf[:,j,i,:] - Wf[:,jm,i,:]) / dy
        else:
            assert mask[j,i] == 0

#vertical derivatives
Uz = np.zeros((nz,ny,nx,nt))
Vz = np.zeros((nz,ny,nx,nt))
Wz = np.zeros((nz,ny,nx,nt))
for i in range(nz):
    imm = np.mod(i+2, nz) #z axis is flipped in MITgcm, so this preserves FDC. Maybe easier to have dz<0???
    im = np.mod(i+1, nz)
    ip = np.mod(i-1, nz)
    ipp = np.mod(i-2, nz)
    if zc4[i]:
        Uz[i,:,:,:] = (Uf[imm,:,:,:]/12 - 2*Uf[im,:,:,:]/3 + 2*Uf[ip,:,:,:]/3 - Uf[ipp,:,:,:]/12) / dz
        Vz[i,:,:,:] = (Vf[imm,:,:,:]/12 - 2*Vf[im,:,:,:]/3 + 2*Vf[ip,:,:,:]/3 - Vf[ipp,:,:,:]/12) / dz
        Wz[i,:,:,:] = (Wf[imm,:,:,:]/12 - 2*Wf[im,:,:,:]/3 + 2*Wf[ip,:,:,:]/3 - Wf[ipp,:,:,:]/12) / dz
    elif zc2[i]:
        Uz[i,:,:,:] = (-Uf[im,:,:,:]/2 + Uf[ip,:,:,:]/2) / dz
        Vz[i,:,:,:] = (-Vf[im,:,:,:]/2 + Vf[ip,:,:,:]/2) / dz
        Wz[i,:,:,:] = (-Wf[im,:,:,:]/2 + Wf[ip,:,:,:]/2) / dz
    elif zf2[i]:
        Uz[i,:,:,:] = (-3*Uf[i,:,:,:]/2 + 2*Uf[ip,:,:,:] - Uf[ipp,:,:,:]/2) / dz
        Vz[i,:,:,:] = (-3*Vf[i,:,:,:]/2 + 2*Vf[ip,:,:,:] - Vf[ipp,:,:,:]/2) / dz
        Wz[i,:,:,:] = (-3*Wf[i,:,:,:]/2 + 2*Wf[ip,:,:,:] - Wf[ipp,:,:,:]/2) / dz
    elif zb2[i]:
        Uz[i,:,:,:] = (3*Uf[i,:,:,:]/2 - 2*Uf[im,:,:,:] + Uf[imm,:,:,:]/2) / dz
        Vz[i,:,:,:] = (3*Vf[i,:,:,:]/2 - 2*Vf[im,:,:,:] + Vf[imm,:,:,:]/2) / dz
        Wz[i,:,:,:] = (3*Wf[i,:,:,:]/2 - 2*Wf[im,:,:,:] + Wf[imm,:,:,:]/2) / dz

#correction velocities
M = (Wf - C*Us) #mean vertical velocity
Uc = (ϵ*Ls/Us) * ((3/(1+2*B) - 1)*(Ut + Uf*Ux + Vf*Uy) + (3/(1+2*B)*Wf - M)*Uz)
Vc = (ϵ*Ls/Us) * ((3/(1+2*B) - 1)*(Vt + Uf*Vx + Vf*Vy) + (3/(1+2*B)*Wf - M)*Vz)
Wc = (ϵ*Ls/Us) * ((3/(1+2*B) - 1)*(Wt + Uf*Wx + Vf*Wy) + (3/(1+2*B)*Wf - M)*Wz)

#Maxey-Riley velocities
Ua = (Uf + Uc) * np.expand_dims(mask,2)
Va = (Vf + Vc) * np.expand_dims(mask,2)
Wa = (M + Wc) * np.expand_dims(mask,2)
for i in range(nt):
    with open(Ua_files[i], "w") as Ua_file:
        Ua[:,:,:,i].astype('>f4').tofile(Ua_file)
    with open(Va_files[i], "w") as Va_file:
        Va[:,:,:,i].astype('>f4').tofile(Va_file)
    with open(Wa_files[i], "w") as Wa_file:
        Wa[:,:,:,i].astype('>f4').tofile(Wa_file)
