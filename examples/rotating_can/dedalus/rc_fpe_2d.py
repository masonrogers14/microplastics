#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 2022

@author: Mason Rogers

rc_fpe_2d.py solves the Fokker-Planck equation for a distribution of particles
in the rotating can flow using Dedalus. Multiprocessing is enabled natively; use
> mpiexec -n <# of procs> rc_fpe_2d.py
"""

#imports
import time
import pathlib
import logging
import numpy as np
from mpi4py import MPI
from rc_param import *
from dedalus import public as d3
from dedalus.extras import flow_tools
logger = logging.getLogger(__name__)

#nondimensional parameters
ϵ = ((1+2*B)*d**2*Us)/(36*ν*Ls) #should be small

# Parameters
Ns, Nr, Nz = 32, 32, 32 
dealias = 3/2
stop_sim_time = tStop + wFreq/100
timestepper = d3.RK443
timestep = 1e-3
dtype = np.float64

# Bases
coords_sr = d3.PolarCoordinates('s', 'r')
coords_z = d3.CartesianCoordinates('z')
dist = d3.Distributor(coords_sr, dtype=dtype)
disk = d3.DiskBasis(coords_sr, shape=(Ns, Nr), radius=R, dealias=dealias, dtype=dtype, azimuth_library='matrix')
edge = disk.S1_basis()
#interval = d3.ChebyshevT(coords_z, size=Nz, bounds=(0.,1.), dealias=3/2) #TODO: change 1 to z_max
s, r = dist.local_grids(disk)
#z = dist.local_grids(interval)

# Fields
p = dist.Field(name='p', bases=disk)
tau_p = dist.Field(name='tau_p', bases=edge)

# Substitutions
#integ = lambda A: d3.Integrate(A, (coords_sr, coords_z))
lift = lambda A, n: d3.Lift(A, disk, n)

# Problem
del(dt)
problem = d3.IVP([p, tau_p], namespace=locals())
problem.add_equation("dt(p) - κ*lap(p) + lift(tau_p,-1) = 0")
problem.add_equation("radial(κ*grad(p)(r=R)) = 0")

# Solver
solver = problem.build_solver(timestepper)
solver.stop_sim_time = stop_sim_time

# Initial conditions
p['g'] = 1/(2*np.pi*Σ) * np.exp(-(r**2 - 2*r*(x0*np.cos(s)+y0*np.sin(s)) + x0**2+y0**2)/(2*Σ))

# Analysis
snapshots = solver.evaluator.add_file_handler('rc_snaps', sim_dt=wFreq, max_writes=1000)
snapshots.add_tasks(solver.state)

# Flow properties
flow = d3.GlobalFlowProperty(solver, cadence=1000)
flow.add_property(p, name='p')

# Main loop
try:
    logger.info('Starting main loop')
    while solver.proceed:
        solver.step(timestep)
        if solver.iteration % 10 == 0:
            one = flow.volume_integral('p')
            logger.info("Iteration=%i, Time=%e, dt=%e, int(p dA)=%e" %(solver.iteration, solver.sim_time, timestep, one))
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    solver.log_stats()

# Post-processing
if dist.comm.rank == 0:
    snapshots.process_virtual_file()

