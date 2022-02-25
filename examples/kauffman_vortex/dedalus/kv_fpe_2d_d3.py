#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 2022

@author: Mason Rogers

kv_fpe_2d.py solves the Fokker-Planck equation for a distribution of particles
in a Kauffman vortex using Dedalus. Multiprocessing is enabled natively; use
> mpiexec -n <# of procs> kv_fpe_2d.py
"""

#imports
import time
import pathlib
import logging
import numpy as np
from mpi4py import MPI
from ../input/kv_param import *
from dedalus import public as d3
from dedalus.extras import flow_tools
logger = logging.getLogger(__name__)

#nondimensional parameters
ϵ = ((1+2*B)*d**2*Us)/(36*ν*Ls) #should be small

# Parameters
Ns, Nr = 32, 128
dealias = 3/2
stop_sim_time = tStop + wFreq/100
timestepper = d3.SBDF2
timestep = 1e-4
dtype = np.float64

# Bases
coords = d3.PolarCoordinates('s', 'r')
dist = d3.Distributor(coords, dtype=dtype)
basis = d3.DiskBasis(coords, shape=(Ns, Nr), radius=R, dealias=dealias, dtype=dtype, azimuth_library='matrix')
s, r = dist.local_grids(basis)
S1_basis = basis.S1_basis(radius=R)

# Fields
p = dist.Field(name='p', bases=basis)
pr = dist.Field(name='pr', bases=basis)
tau_p = dist.Field(name='tau_p', bases=S1_basis)
rHat = dist.VectorField(coords, name='rHat', bases=basis)
rHat['g'][0] = 0
rHat['g'][1] = 1
ad = dist.VectorField(coords, bases=basis)
ad['g'][0] = Γ/2/np.pi/(r**2+a**2) * r
ad['g'][1] = -(Γ/2/np.pi/(r**2+a**2))**2 * r * (ϵ*Ls/Us) * 2*(1-B) / (1+2*B)

# Substitutions
integ = lambda A: d3.Integrate(A, coords)
lift_basis = basis.clone_with(k=2) # Natural output basis
lift = lambda A, n: d3.LiftTau(A, lift_basis, n)

# Problem
del(dt)
problem = d3.IVP([p, tau_p, pr], namespace=locals())
problem.add_equation("dt(p) - κ*lap(p) + lift(tau_p,-1) = -div(p*ad)")
problem.add_equation("pr = dot(rHat, grad(p))")
problem.add_equation("p(r=R) + κ*pr(r=R) = p(r=R) + (dot(rHat, ad)*p)(r=R)")

# Solver
solver = problem.build_solver(timestepper)
solver.stop_sim_time = stop_sim_time

# Initial conditions
Σ = (401/16)*dx**2
x0 = 0.2
y0 = 0.0
p['g'] = 1/(2*np.pi*Σ) * np.exp(-(r**2 - 2*r*(x0*np.cos(s)+y0*np.sin(s)) + x0**2+y0**2)/(2*Σ))

# Analysis
snapshots = solver.evaluator.add_file_handler('kv_snaps', sim_dt=wFreq, max_writes=1000)
snapshots.add_tasks(solver.state)

# Flow properties
flow = d3.GlobalFlowProperty(solver, cadence=100)
flow.add_property(p, name='p')

# Main loop
try:
    logger.info('Starting main loop')
    while solver.proceed:
        solver.step(timestep)
        if (solver.iteration-1) % 1000 == 0:
            max_p = flow.max('p')
            logger.info("Iteration=%i, Time=%e, dt=%e, max(p)=%e" %(solver.iteration, solver.sim_time, timestep, max_p))
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    solver.log_stats()

# Post-processing
if dist.comm.rank == 0:
    snapshots.process_virtual_file()

