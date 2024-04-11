#!/usr/bin/env julia
# -*- coding utf-8 -*-
#=
Created on Tue Dec 7 2021

@author Mason Rogers

kv_sde_init.jl defines the constants and functions required on every processor
used in a multiprocessor ensemble solution of SDEs.
=#

#imports
using Pkg
Pkg.activate(".")
using DifferentialEquations, Printf, Distributed

##########   PARAMETERS   ############################################
#read parameters
io = open("kv_param.py", "r")
param_str_list = readlines(io)
close(io)
for param_str in param_str_list
    param_expr = Meta.parse(param_str)
    if !(typeof(param_expr) == Nothing)
        if param_expr.head == :(=)
            if isinteractive()
                eval(param_expr)
            else
                eval(Expr(:const, param_expr))
            end
        end
    end
end

#nondimensional parameters
const ϵ = ((1+2*B)*d^2*Us)/(36*ν*Ls) #should be small
const C = (2*g*(B-1)*Ls*ϵ)/((1+2*B)*Us^2) #should be O(1)
const A = κ/ϵ/Us/Ls
@printf "A: %.6f\n" A
@printf "C: %.6f\n" C
@printf "ϵ: %.6f\n" ϵ
@printf "κ: %.6f\n" κ

#4d noise magnitude
const α = sqrt(2*Us^3*A/Ls/ϵ)

##########   FLUID DYNAMICS   ########################################
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

##########   DIFFERENTIAL EQUATIONS   ################################
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

##########   INITIAL AND BOUNDARY CONDITIONS   ########################
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

#retrieve initial conditions from process 1
function get_ic(i::Int64)
    return temp_arr[:,i]
end

#continue ensemble
function renew!(p, i, r)
    @printf "proc: %d \n" myid()
    p.u0 .= remotecall_fetch(get_ic, 1, i)
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
cb_out = ContinuousCallback(out_of_domain, reflect!, save_positions=(false,false))

#package callbacks
cb_set = CallbackSet(cb_out)
