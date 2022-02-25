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

#read parameters
io = open("../input/gw_param.py", "r")
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

#dimensional parameters
const k = nk*2π/Lx
const ω = sqrt(g*k)

#4d noise magnitude
const α = sqrt(2*Us^3*A/Ls/ϵ)

#fluid velocities
function fluid_vel(t, x, z)
    u = η0*ω * exp(k*z) * cos(k*x-ω*t)
    w = η0*ω * exp(k*z) * sin(k*x-ω*t)
    return [u, w]
end

#material fluid acceleration
function fluid_acc(t, x, z)
    u, w = fluid_vel(t, x, z)
    au = ω*w
    aw = -ω*u
    return [au, aw]
end

#4d deterministic equations
function mre_det_4d!(ξ̇, ξ, q, t)
    x, z, u, w = ξ
    uᶠ, wᶠ = fluid_vel(t, x, z)
    au, aw = fluid_acc(t, x, z)
    ξ̇[1] = u
    ξ̇[2] = w
    ξ̇[3] = 3/(1+2*B)*au + (36*ν)/((1+2*B)*d^2)*(uᶠ-u)
    ξ̇[4] = 3/(1+2*B)*aw + (36*ν)/((1+2*B)*d^2)*(wᶠ-w) - (2*g*(B-1))/(1+2*B)
end

#4d stochastic terms
function mre_sto_4d!(ξ̇, ξ, q, t)
    ξ̇[1] = 0
    ξ̇[2] = 0
    ξ̇[3] = α
    ξ̇[4] = α
end

#modify ensemble initial conditions
function rand_ic!(p, i, r)
    Σ = Lx^2/1000
    x₁ = x₀ + sqrt(Σ)*randn(2)
    u₁ = fluid_vel(0, x₁...)
    p.u0 .= vcat(x₁, u₁)
    return p
end

#retrieve initial conditions from process 1
function get_ic(i::Int64)
    return temp_arr[:,i]
end

#continue ensemble
function renew!(p, i, r)
    p.u0 .= remotecall_fetch(get_ic, 1, i)
    return p
end

#particle exits domain event
function out_of_domain(ξ,t,integrator)
    x, z, u, w = ξ
    return z - η0*cos(k*x-ω*t)
end
#particle reflects off boundary
# function reflect!(integrator)
#     x, z, u, w = integrator.u
#     sn = (u*x + v*y)/sqrt(x^2 + y^2)
#     integrator.u[3] -= 2*sn*x/sqrt(x^2+y^2)
#     integrator.u[4] -= 2*sn*y/sqrt(x^2+y^2)
# end
cb_out = ContinuousCallback(out_of_domain, terminate!, save_positions=(false,false))

#package callbacks
cb_set = CallbackSet(cb_out)

#= XXXX
#2d particle velocity
#To do: implement function parti_vel(t, x, y)

#2d deterministic equations
function mre_det_2d!(ξ̇, ξ, q, t)
    x, z = ξ
    u, w = parti_vel(t, x, z)
    ξ̇[1] = u
    ξ̇[2] = w
end

#2d stochastic terms
function mre_sto_2d!(ξ̇, ξ, q, t)
    ξ̇[1] = sqrt(2*κ)
    ξ̇[2] = sqrt(2*κ)
end
=#
