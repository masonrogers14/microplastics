#!/usr/bin/env julia
# -*- coding utf-8 -*-
#=
Created on Wed Apr 17 2024

@author Mason Rogers

new_kv_sde_4d.jl is the parent script for a simulation of an ensemble of
particles in a Kaufmann vortex. Running via SLURM automatically enables
multiprocessing with all available processors.
=#

#imports
using Pkg
Pkg.activate(".")
using Distributed
using ClusterManagers

#multiprocessor setup
addprocs_slurm(parse(Int, ENV["SLURM_NTASKS"]) - 1,
               nodes=parse(Int, ENV["SLURM_NNODES"]),
               exename="/home/software/julia/1.8.5/bin/julia")

#tinker
@everywhere nDim = 4
@everywhere nTraj = 10000000
@everywhere saveTraj = true
@everywhere saveHist = false
@everywhere packGrid = true
@everywhere dir = "/pool001/masonr/kv4d/"
@everywhere t_prefix = dir*"parallel"
@everywhere h_prefix = dir*"parallel"
@everywhere initTime = 0.

#initialize
@everywhere include("kv_sde_param.jl")
@everywhere include("kv_sde_init.jl")
@everywhere include("kv_traj2hist.jl")

#save grid if necessary
if packGrid
    pack_grid()
end

#run
@everywhere run_sde()
