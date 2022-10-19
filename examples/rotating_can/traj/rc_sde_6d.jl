#!/usr/bin/env julia
# -*- coding utf-8 -*-
#=
Created on Mon Oct 17 2022 

@author Mason Rogers

rc_sde_6d.jl is the parent script for a simulation of an ensemble of particles
in the rotating can flow. Running via SLURM automatically enables multiprocessing
with all available processors.
=#

#tinker
nTraj = 500000
saveTraj = false
saveHist = true
packGrid = true
dir = "../output/exp1/"
t_prefix = dir*"julia6dx"
h_prefix = dir*"julia6dx"
initTime = 0
nDims = 6

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
@everywhere include("rc_sde_init.jl")
include("rc_traj2hist.jl")

#initialize storage arrays
temp_arr = NaN*zeros(nDims, nTraj)
inMemory = nTraj * nOuts < 1e9
if inMemory
    println("in memory")
    full_arr = NaN * zeros(nDims, nTraj, nOuts)
else
    println("not in memory")
end

#save grid if necessary
if packGrid
    pack_grid()
end

#choose initial conditions: either random or from file
println("setting initial conditions")
if initTime == 0
    #initial conditions
    @everywhere x₀ = [x0,y0,z0] #supplied in rc_param.py
    u₀ = fluid_vel(0, x₀...)
    ξ₀ = vcat(x₀, u₀)
    
    init_prob = SDEProblem(mre_det_6d!, mre_sto_6d!, ξ₀, (0.0, eps()), save_everystep=false, save_end=false)
    init_ense = EnsembleProblem(init_prob, prob_func=rand_ic_6d!)
    init_solu = solve(init_ense, SOSRI(), EnsembleDistributed(), trajectories=nTraj, dt=5e-4, adaptive=false)
    for i in 1:nTraj
	temp_arr[:,i] = init_solu[i][1]
	#save trajectories
	if saveTraj
	    save_trajectories(0)
	end
	#compute and save histogram data in MITgcm format
	if saveHist
	    save_histogram(0)
	end
    end
else
    t_suffix = @sprintf ".%010d.bin" Int(round(initTime/dt))
    initFile = t_prefix*t_suffix
    read!(initFile, temp_arr)
    #compute and save histogram data in MITgcm format
    if saveHist
	save_histogram(0)
    end
end
if inMemory
    full_arr[:,:,1] = temp_arr
end

#run solver (either in memory or in chunks of wFreq)
println("running")
if inMemory
    prob = SDEProblem(mre_det_6d!, mre_sto_6d!, zeros(nDims), (0, wFreq*(nOuts-1)), saveat=0:wFreq:wFreq*(nOuts-1))
    ense = EnsembleProblem(prob, prob_func=renew!)
    solu = solve(ense, SOSRI(), EnsembleDistributed(), trajectories=nTraj, callback=cb_set, dt=5e-4, adaptive=false)
    for i in 1:nTraj
        full_arr[:,i,:] = solu[i][:,:]
    end
    for j in 1:nOuts-1
        temp_arr[:,:] = full_arr[:,:,j+1] 
        #save trajectories
        if saveTraj || j == nOuts-1
            save_trajectories(j)
        end
        #compute and save histogram data in MITgcm format
        if saveHist 
            save_histogram(j)
        end
    end
else 
    for j in 1:nOuts-1
        prob = SDEProblem(mre_det_6d!, mre_sto_6d!, zeros(nDims), (wFreq*(j-1),wFreq*j), save_everystep=false, save_end=true)
        ense = EnsembleProblem(prob, prob_func=renew!)
        solu = solve(ense, SOSRI(), EnsembleDistributed(), trajectories=nTraj, callback=cb_set, dt=5e-4, adaptive=false)
        for i in 1:nTraj
            temp_arr[:,i] = solu[i][end]
        end
        #save trajectories
        if saveTraj || j == nOuts-1
            save_trajectories(j)
        end
        #compute and save histogram data in MITgcm format
        if saveHist 
            save_histogram(j)
        end
    end
end
