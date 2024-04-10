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
nTraj = 100
saveTraj = true
saveHist = true
packGrid = true
dir = "/pool001/masonr/kv4d/"
t_prefix = dir*"julia4d"
h_prefix = dir*"julia4d"
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
temp_arr = NaN*zeros(4,nTraj)
inMemory = nTraj * nOuts < 1e9
if inMemory
    println("in memory")
    full_arr = NaN * zeros(4, nTraj, nOuts)
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
    @everywhere x₀ = [x0,y0] #supplied in kv_param.py
    u₀ = fluid_vel(0, x₀...)
    ξ₀ = vcat(x₀, u₀)
    
    init_prob = SDEProblem(mre_det_4d!, mre_sto_4d!, ξ₀, (0.0,wFreq), save_everystep=false, save_end=false)
    init_ense = EnsembleProblem(init_prob, prob_func=rand_ic_4d!)
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
    prob = SDEProblem(mre_det_4d!, mre_sto_4d!, zeros(4), (0, wFreq*(nOuts-1)), saveat=0:wFreq:wFreq*(nOuts-1))
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
        prob = SDEProblem(mre_det_4d!, mre_sto_4d!, zeros(4), (wFreq*(j-1),wFreq*j), save_everystep=false, save_end=true)
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
