#!/usr/bin/env julia
# -*- coding utf-8 -*-
#=
Created on Tue Dec 7 2021

@author Mason Rogers
=#

#tinker
nTraj = 1
saveTraj = false
saveHist = false
packGrid = false
dir = "../output/"
t_fname = dir*"4dtraj.bin"
h_prefix = dir*"4dtraj"

#activate environment
using Pkg
Pkg.activate(".")

#add processors
using DifferentialEquations, StatsBase, Distributed
if "SLURM_NTASKS" in keys(ENV)
    using ClusterManagers
    addprocs(SlurmManager(parse(Int,ENV["SLURM_NTASKS"])-1))
end

#imports
@everywhere include("kv_sde_init.jl")
include("kv_traj2hist.jl")

#parameters
const nOuts = Int(floor(tStop/wFreq) + 1)

#---------------------------------
@everywhere x₀ = [0.,0.]
u₀ = fluid_vel(0, x₀...)
ξ₀ = vcat(x₀, u₀)

prob = SDEProblem(mre_det_4d!, mre_sto_4d!, ξ₀, (0.0,tStop+wFreq/100), save_everystep=false, save_end=false)
ense = EnsembleProblem(prob, prob_func=change_ic!)
solu = solve(ense, SOSRI(), EnsembleDistributed(), saveat=0:wFreq:tStop, save_idxs=[1,2], trajectories=nTraj, callback=cb_set)

#package ensemble in array
solu_arr = NaN * zeros(nTraj, nOuts, 2)
for i in 1:nTraj
    for j in 1:size(solu[i])[2]
        solu_arr[i,j,:] = solu[i][j]
    end
end

#save trajectories
if saveTraj
    save_trajectories()
end

#compute and save histogram data in MITgcm format
if saveHist
    save_histogram()
end

#save grid if necessary
if packGrid
    pack_grid()
end
