#!/usr/bin/env julia
# -*- coding utf-8 -*-
#=
Created on Tue Dec 7 2021

@author Mason Rogers

kv_sde_4d.jl is the parent script for a simulation of an ensemble of particles
in a Kauffman vortex. Running via SLURM automatically enables multiprocessing
with all available processors.
=#

#tinker
nTraj = 100000
saveTraj = false
saveHist = true
packGrid = true
dir = "../output/"
t_prefix = dir*"4d"
h_prefix = dir*"4d"

#activate environment
using Pkg
Pkg.activate(".")

#add processors
using Distributed
if "SLURM_NTASKS" in keys(ENV)
    using ClusterManagers
    addprocs(SlurmManager(parse(Int,ENV["SLURM_NTASKS"])-1))
end

#imports
@everywhere include("gw_sde_init.jl")
include("gw_traj2hist.jl")

#initialize storage arrays
temp_arr = NaN*zeros(4,nTraj)
inMemory = nTraj * nOuts < 1e5
if inMemory
    println("in memory")
    full_arr = NaN * zeros(4, nTraj, nOuts)
else
    println("not in memory")
end

#initial conditions
@everywhere x₀ = [Lx/2, -Lz/4]
u₀ = fluid_vel(0, x₀...)
ξ₀ = vcat(x₀, u₀)

#save grid if necessary
if packGrid
    pack_grid()
end

#add randomness to initial conditions
init_prob = SDEProblem(mre_det_4d!, mre_sto_4d!, ξ₀, (0.0,wFreq), save_everystep=false, save_end=false)
init_ense = EnsembleProblem(init_prob, prob_func=rand_ic!)
init_solu = solve(init_ense, SOSRI(), EnsembleDistributed(), trajectories=nTraj, dt=1e-3, adaptive=false)
for i in 1:nTraj
    temp_arr[:,i] = init_solu[i][1]
    if inMemory
        full_arr[:,i,1] = init_solu[i][1]
    end
    #save trajectories
    if saveTraj
        save_trajectories(0)
    end
    #compute and save histogram data in MITgcm format
    if saveHist
        save_histogram(0)
    end
end

#run solver in chunks of wFreq
for j in 1:nOuts-1
    prob = SDEProblem(mre_det_4d!, mre_sto_4d!, ξ₀, (wFreq*(j-1),wFreq*j), save_everystep=false, save_end=true)
    ense = EnsembleProblem(prob, prob_func=renew!)
    solu = solve(ense, SOSRI(), EnsembleDistributed(), trajectories=nTraj, callback=cb_set, dt=1e-3, adaptive=false)
    for i in 1:nTraj
        temp_arr[:,i] = solu[i][end]
        if inMemory
            full_arr[:,i,j+1] = solu[i][end]
        end
    end
    #save trajectories
    if saveTraj
        save_trajectories(j)
    end
    #compute and save histogram data in MITgcm format
    if saveHist
        save_histogram(j)
    end
end
