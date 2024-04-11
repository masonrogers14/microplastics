#!/usr/bin/env julia
# -*- coding utf-8 -*-
#=
Created on Tue Dec 7 2021

@author Mason Rogers

kv_sde_4d.jl is the parent script for a simulation of an ensemble of particles
in a Kaufmann vortex. Running via SLURM automatically enables multiprocessing
with all available processors.
=#

#tinker
nTraj = 10000000
saveTraj = true
saveHist = true
packGrid = true
dir = "/pool001/masonr/kv4d/"
t_prefix = dir*"julia4d7"
h_prefix = dir*"julia4d7"
initTime = 0.

#activate environment
using Pkg
Pkg.activate(".")
using Printf

#add processors
using Distributed
if "SLURM_NTASKS" in keys(ENV)
    @printf "slurm detected, adding %s processors\n" ENV["SLURM_NTASKS"]
    using ClusterManagers
    addprocs(SlurmManager(parse(Int,ENV["SLURM_NTASKS"])-1))
end

#imports
@everywhere include("kv_sde_init.jl")
include("kv_traj2hist.jl")

#initialize storage arrays
temp_arr = NaN * zeros(4,nTraj)
full_arr = NaN * zeros(4, nTraj, nOuts)

#save grid if necessary
if packGrid
    pack_grid()
end

#run solver (either in memory or in chunks of wFreq)
println("running")
@everywhere x₀ = [x0,y0] #supplied in kv_param.py
prob = SDEProblem(mre_det_4d!, mre_sto_4d!, zeros(4), (0, wFreq*(nOuts-1)), saveat=0:wFreq:wFreq*(nOuts-1))
ense = EnsembleProblem(prob, prob_func=rand_ic_4d!)
solu = solve(ense, SOSRI(), EnsembleDistributed(), trajectories=nTraj, callback=cb_set, dt=5e-4, adaptive=false)
for i in 1:nTraj
    full_arr[:,i,:] = solu[i][:,:]
end
for j in 1:nOuts
    temp_arr[:,:] = full_arr[:,:,j] 
    #save trajectories
    if saveTraj || j == nOuts
        save_trajectories(j)
    end
    #compute and save histogram data in MITgcm format
    if saveHist 
        save_histogram(j)
    end
end
