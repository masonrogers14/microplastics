"""
Dedalus script for 2D advection-diffusion
To run, merge, and plot using 4 processes, for instance, you could use:
    $ mpiexec -n 2 python3 dww.py
    $ mpiexec -n 2 python3 -m dedalus merge_procs test_rv
    $ mpiexec -n 2 python3 plot_slices.py test_rv/*.h5
The simulations should take a few process-minutes to run.
"""

import numpy as np
from mpi4py import MPI
import time
import pathlib

from dedalus import public as de
from dedalus.extras import flow_tools

import logging
logger = logging.getLogger(__name__)

from scipy.stats import norm

#Maxey-Riley parameters
B = 0.98 #density ratio ρ/ρᶠ
ν = 1e-6 #fluid viscosity
d = 1e-3 #particle diameter
g = 9.8 #gravity
Us = 1.0 #characteristic velocity scale of flow
Ls = 1.0 #characteristic length scale of flow

#dimensional 2d diffusivity
κ = 4e-5

#nondimensional parameters
ϵ = ((1+2*B)*d**2*Us)/(36*ν*Ls) #should be small

#Rankine vortex parameters
R = .5
a = .2
Γ = 1.

# Create bases and domain
s_basis = de.Fourier('s', 256, interval=(0, 2*np.pi), dealias=3/2) #s <--> θ
r_basis = de.Chebyshev('r', 256, interval=(0.0, R), dealias=3/2)
domain = de.Domain([s_basis, r_basis], grid_dtype=np.float64)
s, r = domain.grids(scales=1)

#NCCs
ad_s = domain.new_field(name='ad_s')
ad_s['g'] = Γ/2/np.pi/(r**2+a**2) * r
ad_s.meta['s']['constant'] = True
ad_r = domain.new_field(name='ad_r')
ad_r['g'] = -(Γ/2/np.pi/(r**2+a**2))**2 * r * (ϵ*Ls/Us) * (1-B) / (1+2*B)
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
Σ = (401/4)*(.01)**2
p['g'] = norm.pdf(r, scale=np.sqrt(Σ)) / np.sqrt(2*np.pi*Σ)
p.differentiate('r', out=pr)
p.differentiate('s', out=ps)

# p['c']
# from dedalus.extras.plot_tools import plot_bot_2d
# import matplotlib.pyplot as plt
# log_mag = lambda xmesh, ymesh, data: (xmesh, ymesh, np.log10(np.abs(data)))
# plot_bot_2d(p, func=log_mag, clim=(-20, 0), cmap='viridis', title="log10(abs(f['c'])");
# plt.show()

# Timestepping and output
dt = 1e-3
stop_sim_time = 1000.1
fh_mode = 'overwrite'

# Integration parameters
solver.stop_sim_time = stop_sim_time

# Analysis
snapshots = solver.evaluator.add_file_handler('kv_snaps', sim_dt=20, max_writes=101, mode=fh_mode)
snapshots.add_system(solver.state)
#snapshots.add_task("integ(p*r)", name="E[1]")
#snapshots.add_task("integ(p*r**2)", name="E[r]")
#snapshots.add_task("integ(p*r*s)", name="E[s]")
#snapshots.add_task("integ(p*r**3)", name="E[r^2]")
#snapshots.add_task("integ(p*r*s**2)", name="E[s^2]")
#snapshots.add_task("integ(p*r**2*s)", name="E[rs]")
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
        if (solver.iteration-1) % 10 == 0:
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

