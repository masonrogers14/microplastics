#!/usr/bin/env julia
# -*- coding utf-8 -*-
#=
Created on Wed Apr 17 2024

@author Mason Rogers

new_kv_sde_init.jl defines the constants and functions required on every
processor used in a multiprocessor ensemble solution of SDEs.
=#

#imports
using Pkg
Pkg.activate(".")
using StochasticDiffEq, DiffEqCallbacks, Printf, Distributed

##########   FLUID FIELDS   ####################################################
#fluid velocities
function fluid_vel(t, x, y)
    ρ = sqrt(x^2 + y^2)
    θ = atan(y, x)
    ω = Γ/2/π/(ρ^2+a^2)
    u = -ρ*ω*sin(θ)
    v = ρ*ω*cos(θ)
    return [u, v]
end

#material fluid acceleration
function fluid_acc(t, x, y)
    ρ = sqrt(x^2 + y^2)
    θ = atan(y, x)
    ω = Γ/2/π/(ρ^2+a^2)
    au = -ρ*ω^2*cos(θ)
    av = -ρ*ω^2*sin(θ)
    return [au, av]
end

#2d particle velocity
function parti_vel(t, x, y)
    uᶠ, vᶠ = fluid_vel(t, x, y)
    au, av = fluid_acc(t, x, y)
    u = uᶠ + 2*(ϵ*Ls/Us)*(1-B)/(1+2*B)*au
    v = vᶠ + 2*(ϵ*Ls/Us)*(1-B)/(1+2*B)*av
    return [u, v]
end

##########   DIFFERENTIAL EQUATIONS   ##########################################
#4d deterministic equations
function mre_det_4d!(ξ̇, ξ, q, t)
    x, y, u, v = ξ
    uᶠ, vᶠ = fluid_vel(t, x, y)
    au, av = fluid_acc(t, x, y)
    ξ̇[1] = u
    ξ̇[2] = v
    ξ̇[3] = 3/(1+2*B)*au + (Us/(Ls*ϵ))*(uᶠ-u)
    ξ̇[4] = 3/(1+2*B)*av + (Us/(Ls*ϵ))*(vᶠ-v)
end

#4d stochastic terms
function mre_sto_4d!(ξ̇, ξ, q, t)
    ξ̇[1] = 0
    ξ̇[2] = 0
    ξ̇[3] = α
    ξ̇[4] = α
end

#2d deterministic equations
function mre_det_2d!(ξ̇, ξ, q, t)
    x, y = ξ
    u, v = parti_vel(t, x, y)
    ξ̇[1] = u
    ξ̇[2] = v
end

#2d stochastic terms
function mre_sto_2d!(ξ̇, ξ, q, t)
    ξ̇[1] = sqrt(2*κ)
    ξ̇[2] = sqrt(2*κ)
end

##########   INITIAL AND BOUNDARY CONDITIONS   #################################
#4d random ensemble initial conditions
function rand_ic_4d!(p, i, r)
    x₁ = x₀ + sqrt(Σ)*randn(2)
    u₁ = fluid_vel(0, x₁...)
    p.u0 .= vcat(x₁, u₁)
    return p
end

function rand_ic_2d!(p, i, r)
    x₁ = x₀ + sqrt(Σ)*randn(2)
    p.u0 .= x₁
    return p
end

#continue ensemble
function renew!(p, i, r)
    p.u0 .= step_arr[:,i]
    return p
end

#particle exits domain event
function out_of_domain(ξ, t, integrator)
    ρ = sqrt(ξ[1]^2+ξ[2]^2)
    return ρ-R
end

#particle reflects off boundary
function reflect!(integrator)
    #x, y, u, v = integrator.u
    #sn = (u*x + v*y)/sqrt(x^2 + y^2)
    #integrator.u[3] -= 2*sn*x/sqrt(x^2+y^2)
    #integrator.u[4] -= 2*sn*y/sqrt(x^2+y^2)
end
cb_out = ContinuousCallback(out_of_domain,
                            reflect!,
                            save_positions=(false,false))

#package callbacks
cb_set = CallbackSet(cb_out)

##########   ENSEMBLE PROBLEMS   ###############################################
#initialize storage arrays
step_arr = NaN * zeros(nDim, nPerProc)

#initial conditions
x₀ = [x0, y0] #supplied in kv_param.py
u₀ = fluid_vel(0, x₀...)

function run_sde()
    #decide between 2D or 4D
    if nDim == 4
        ic_func! = rand_ic_4d!
        det_func! = mre_det_4d!
        sto_func! = mre_sto_4d!
        ξ₀ = vcat(x₀, u₀)
    elseif nDim == 2
        ic_func! = rand_ic_2d!
        det_func! = mre_det_2d!
        sto_func! = mre_sto_2d!
        ξ₀ = x₀
    end

    #choose initial conditions: either random or from file
    println("setting initial conditions")
    if initTime == 0
        #initial conditions        
        init_prob = SDEProblem(det_func!, sto_func!, ξ₀, (0.0, 5e-4),
                               save_everystep=false, save_end=false)
        init_ense = EnsembleProblem(init_prob, prob_func=ic_func!)
        init_solu = solve(init_ense, SOSRI(), EnsembleThreads(), 
                          trajectories=nPerProc, dt=5e-4, adaptive=false)
        for i in 1:nPerProc
            step_arr[:,i] = init_solu[i][1]
        end
        #save trajectories
        if saveTraj
            save_trajectories(0)
        end
        #compute and save histogram data in MITgcm format
        if saveHist
            save_histogram(0)
        end 
    else
        t_suffix = @sprintf ".%010d_%04d.bin" Int(round(initTime/dt)) myid()
        initFile = t_prefix*t_suffix
        read!(initFile, step_arr)
        #compute and save histogram data in MITgcm format
        if saveHist
            save_histogram(0)
        end
    end

    #run solver (either in memory or in chunks of wFreq)
    println("running")
    for j in 1:nOuts-1
        prob = SDEProblem(det_func!,
                          sto_func!, 
                          zeros(nDim), 
                          (wFreq*(j-1), wFreq*j),
                          save_everystep=false, save_end=true)
        ense = EnsembleProblem(prob, prob_func=renew!)
        solu = solve(ense, SOSRI(), EnsembleThreads(),
                    trajectories=nPerProc, callback=cb_set, dt=5e-4,
                    adaptive=false)
        for i in 1:nPerProc
            step_arr[:,i] = solu[i][end]
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
