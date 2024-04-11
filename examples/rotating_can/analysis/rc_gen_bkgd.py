#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 2021

@author: Mason Rogers
"""

#imports
import numpy as np
import matplotlib.pyplot as plt
from rc_param import *
from scipy.integrate import solve_ivp

#directory
dname = "eps"+"{0:02d}".format(int(100*δ))+"/"

#velocities
def U(t,r):
    ux = r[0]; uy = r[1]; uz = r[2]
    ur = np.sqrt(ux**2 + uy**2)
    return -b*ux*(1-2*uz)*(R-ur)/3 - a*uy*(c+uz**2) + δ*(1-β*uz)*(uy*(uy-ys+γ*np.cos(σ*t)) - (R**2-ur**2)/2)

def V(t,r):
    vx = r[0]; vy = r[1]; vz = r[2]
    vr = np.sqrt(vx**2 + vy**2)
    return -b*vy*(1-2*vz)*(R-vr)/3 + a*vx*(c+vz**2) - δ*(1-β*vz)*(vx*(vy-ys+γ*np.cos(σ*t)))

def W(t,r):
    wx = r[0]; wy = r[1]; wz = r[2]
    wr = np.sqrt(wx**2 + wy**2)
    return b*wz*(1-wz)*(2*R/3-wr)

def vel(t,r):
    return np.array([U(t,r), V(t,r), W(t,r)])

#events for Poincaré section
def y(t,r):
    return r[1]

# solver
for i in range(20):
    ρ0 = np.sqrt(np.random.uniform(0,R**2))
    θ0 = np.random.uniform(0,2*np.pi)
    z0 = np.random.uniform(0,H)
    r0 = np.array([ρ0*np.cos(θ0), ρ0*np.sin(θ0), z0])
    print(r0)
    sol = solve_ivp(vel, [0, 1000], r0, events=y, first_step=1e-2, max_step=1e-1)
    px = sol.y_events[0][:,0]
    pz = sol.y_events[0][:,2] 
    with open(dname+"poincare_x_"+str(i)+".bin", "w") as file:
        px.tofile(file)
    with open(dname+"poincare_z_"+str(i)+".bin", "w") as file:
        pz.tofile(file)
    print(i)

# #plotter
# for i in range(20):
#     x_fname = dname+"poincare_x_"+str(i)+".bin"
#     z_fname = dname+"poincare_z_"+str(i)+".bin"
#     with open(x_fname, "r") as file:
#         x_pts = np.fromfile(file)
#     with open(z_fname, "r") as file:
#         z_pts = np.fromfile(file)
#     plt.scatter(x_pts, z_pts, c='black', s=4)
# plt.show()
