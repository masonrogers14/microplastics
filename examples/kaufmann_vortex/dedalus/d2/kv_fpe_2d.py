#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 2021

@author: Mason Rogers

kv_fpe_2d.py solves the Fokker-Planck equation for a distribution of particles
in a Kaufmann vortex using Dedalus. Multiprocessing is enabled natively; use
> mpiexec -n <# of procs> kv_fpe_2d.py
"""

#imports
import time
import pathlib
import logging
import numpy as np
from mpi4py import MPI
from ../input/kv_param import *
from dedalus import public as de
from dedalus.extras import flow_tools
logger = logging.getLogger(__name__)

#nondimensional parameters
ϵ = ((1+2*B)*d**2*Us)/(36*ν*Ls) #should be small

# Create bases and domain
s_basis = de.Fourier('s', 64, interval=(0, 2*np.pi), dealias=3/2) #s <--> θ
r_basis = de.Chebyshev('r', 64, interval=(0.0, R), dealias=3/2)
domain = de.Domain([s_basis, r_basis], grid_dtype=np.float64)
s, r = domain.grids(scales=1)

#NCCs
ad_s = domain.new_field(name='ad_s')
ad_s['g'] = Γ/2/np.pi/(r**2+a**2) * r
ad_s.meta['s']['constant'] = True
ad_r = domain.new_field(name='ad_r')
ad_r['g'] = -(Γ/2/np.pi/(r**2+a**2))**2 * r * (ϵ*Ls/Us) * 2*(1-B) / (1+2*B)
ad_r.meta['s']['constant'] = True

# 2D equations
problem = de.IVP(domain, variables=['p','ps','pr'])
problem.parameters['Us'] = Us
problem.parameters['Ls'] = Ls
problem.parameters['eps'] = ϵ
problem.parameters['kap'] = κ
problem.parameters['ad_r'] = ad_r
problem.parameters['ad_s'] = ad_s
problem.add_equation("dt(p)*r**2 - kap*(dr(r*pr)*r + ds(ps)) = -(ad_r*pr*r**2 + ad_s*ps*r + p*dr(r*ad_r)*r + p*ds(ad_s)*r)")
problem.add_equation("ps - ds(p) = 0")
problem.add_equation("pr - dr(p) = 0")

#BCs
problem.add_bc("left(p) = 0", condition="ns != 0")
problem.add_bc("left(pr) = 0", condition="ns == 0")
problem.add_bc("right(kap*pr - ad_r*p) = 0")

# Build solver
solver = problem.build_solver(de.timesteppers.RK443)
logger.info('Solver built')

# Initial conditions
s, r = domain.all_grids()
p = solver.state['p']
pr = solver.state['pr']
ps = solver.state['ps']
Σ = (401/16)*dx**2
x0 = 0.2
y0 = 0.0
p['g'] = 1/(2*np.pi*Σ) * np.exp(-(r**2 - 2*r*(x0*np.cos(s)+y0*np.sin(s)) + x0**2+y0**2)/(2*Σ))
p.differentiate('r', out=pr)
p.differentiate('s', out=ps)

# p['c']
# from dedalus.extras.plot_tools import plot_bot_2d
# import matplotlib.pyplot as plt
# log_mag = lambda xmesh, ymesh, data: (xmesh, ymesh, np.log10(np.abs(data)))
# plot_bot_2d(p, func=log_mag, clim=(-20, 0), cmap='viridis', title="log10(abs(f['c'])");
# plt.show()

# Timestepping and output
stop_sim_time = tStop + wFreq/100
fh_mode = 'overwrite'

# Integration parameters
solver.stop_sim_time = stop_sim_time

# Analysis
snapshots = solver.evaluator.add_file_handler('kv_snaps', sim_dt=wFreq, max_writes=1001, mode=fh_mode)
snapshots.add_system(solver.state)
snapshots.add_task("integ(p*r)", name="E[1]")
snapshots.add_task("integ(p*r**2)", name="E[r]")
snapshots.add_task("integ(p*r*s)", name="E[s]")
snapshots.add_task("integ(p*r**3)", name="E[r^2]")
snapshots.add_task("integ(p*r*s**2)", name="E[s^2]")
snapshots.add_task("integ(p*r**2*s)", name="E[rs]")
snapshots.add_task("p", name="p_c", layout="c")

# CFL
CFL = flow_tools.CFL(solver, initial_dt=dt, safety=0.5, max_change=1.05)
CFL.add_velocities(('ad_s/r', 'ad_r'))

# Main loop
try:
    logger.info('Starting loop')
    start_time = time.time()
    while solver.proceed:
        dt = CFL.compute_dt()
        dt = solver.step(dt)
        if (solver.iteration-1) % 100 == 0:
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
