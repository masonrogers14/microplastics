#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 7 2022

@author: Mason Rogers

gw_fpe_2d.py solves the Fokker-Planck equation for a distribution of particles
in a field of deep water gravity waves using Dedalus. Multiprocessing is enabled
natively; use
> mpiexec -n <# of procs> gw_fpe_2d.py
"""

#imports
import time
import pathlib
import logging
import numpy as np
from mpi4py import MPI
from gw_param import *
from dedalus import public as de
from dedalus.extras import flow_tools
logger = logging.getLogger(__name__)

#nondimensional parameters
ϵ = ((1+2*B)*d**2*Us)/(36*ν*Ls) #should be small
C = (2*g*(B-1)*Ls*ϵ)/((1+2*B)*Us**2)

#dimensional parameters
k = nk*2*np.pi/Lx
ω = np.sqrt(g*k)

# Create bases and domain
x_basis = de.Fourier('x', 64, interval=(0, Lx), dealias=3/2)
ζ1_basis = de.Chebyshev('ζ1', 32, interval=(-Lz, -Lz/4), dealias=3/2)
ζ2_basis = de.Chebyshev('ζ2', 32, interval=(-Lz/4, 0), dealias=3/2)
ζ_basis = de.Compound('ζ', (ζ1_basis, ζ2_basis))
# ζ_basis = de.Laguerre('ζ', 1024, edge=0, stretch=-1, dealias=3/2)
domain = de.Domain([x_basis, ζ_basis], grid_dtype=np.float64)
x, ζ = domain.grids(scales=1)

# 2D equations
problem = de.IVP(domain, variables=['p','px','pz'])
problem.parameters['Us'] = Us
problem.parameters['Ls'] = Ls
problem.parameters['B'] = B
problem.parameters['C'] = C
problem.parameters['ε'] = ϵ
problem.parameters['κ'] = κ
problem.parameters['g'] = g
problem.parameters['η0'] = η0
problem.parameters['k'] = k
problem.parameters['ω'] = ω
problem.substitutions['η(t,x)'] = 'η0 * cos(k*x-ω*t)'
problem.substitutions['ηt(t,x)'] = 'ω*η0 * sin(k*x-ω*t)'
problem.substitutions['ηx(t,x)'] = '-k*η0 * sin(k*x-ω*t)'
problem.substitutions['uF(t,x,ζ)'] = 'ω*η0 * cos(k*x-ω*t) * exp(k*(ζ+η(t,x)))'
problem.substitutions['wF(t,x,ζ)'] = 'ω*η0 * sin(k*x-ω*t) * exp(k*(ζ+η(t,x)))'
problem.substitutions['uA(t,x,ζ)'] = 'ω*wF(t,x,ζ)'
problem.substitutions['wA(t,x,ζ)'] = '-ω*uF(t,x,ζ)'
problem.substitutions['ad_x(t,x,ζ)'] = 'p*(uF(t,x,ζ) + 2*ϵ*Ls/Us*(1-B)/(1+2*B)*uA(t,x,ζ))'
problem.substitutions['ad_z(t,x,ζ)'] = 'p*(wF(t,x,ζ) + 2*ϵ*Ls/Us*(1-B)/(1+2*B)*wA(t,x,ζ))' #missing -C because of envelope issues in Laguerre basis
problem.substitutions['adv(t,x,ζ)'] = 'dx(ad_x(t,x,ζ)) - ηx(t,x)*dζ(ad_x(t,x,ζ)) + dζ(ad_z(t,x,ζ)) - C*pz'
problem.substitutions['dif(t,x,ζ)'] = 'κ*(dx(px) - ηx(t,x)*dζ(px) + dζ(pz))'
problem.substitutions['vbf(t,x,ζ)'] = 'p*ηt(t,x)'
# problem.meta['η']['ζ']['envelope'] = False
# problem.meta['ηt']['ζ']['envelope'] = False
# problem.meta['ηx']['ζ']['envelope'] = False
problem.add_equation("dt(p) = -adv(t,x,ζ) + dif(t,x,ζ) + ηt(t,x)*dζ(p)")
problem.add_equation("px = dx(p) - ηx(t,x)*dζ(p)") #true x
problem.add_equation("pz - dζ(p) = 0", tau=2) #true z

#BCs
problem.add_bc("right(-C*p - κ*pz) = right(-ad_z(t,x,ζ) + vbf(t,x,ζ))")
problem.add_bc("left(pz) = 0")

# Build solver
solver = problem.build_solver(de.timesteppers.RK443)
logger.info('Solver built')

# Initial conditions
x, ζ = domain.all_grids()
p = solver.state['p']
px = solver.state['px']
pz = solver.state['pz']
M = Lx**2/1000
p['g'] = (2*np.pi*M)**-1 * np.exp(-((x-Lx/2)**2+(ζ+Lz/4)**2)/(2*M))
p.differentiate('x', out=px)
p.differentiate('ζ', out=pz)

# Timestepping and output
dt = 1e-3
stop_sim_time = 1.01
fh_mode = 'overwrite'

# Integration parameters
solver.stop_sim_time = stop_sim_time

# Analysis
snapshots = solver.evaluator.add_file_handler('gw_snaps', sim_dt=.1, max_writes=101, mode=fh_mode)
snapshots.add_system(solver.state)
snapshots.add_task("integ(p)", name="E[1]")
snapshots.add_task("integ(p*x)", name="E[x]")

# CFL
CFL = flow_tools.CFL(solver, initial_dt=dt, safety=0.1, max_change=1.05)
CFL.add_velocities(('uF(t,x,ζ)', 'wF(t,x,ζ)-C'))

# Main loop
try:
    logger.info('Starting loop')
    start_time = time.time()
    while solver.proceed:
        dt = CFL.compute_dt()
        dt = solver.step(dt)
        if (solver.iteration-1) % 10 == 0:
            if np.max(np.abs(solver.state['p']['g'])) > 1000: raise
            logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    end_time = time.time()
    logger.info('Iterations: %i' %solver.iteration)
    logger.info('Sim end time: %f' %solver.sim_time)
    logger.info('Run time: %.2f sec' %(end_time-start_time))
    logger.info('Run time: %f cpu-hr' %((end_time-start_time)/60/60*domain.dist.comm_cart.size))
