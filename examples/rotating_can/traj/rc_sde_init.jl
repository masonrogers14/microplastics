#!/usr/bin/env julia
# -*- coding utf-8 -*-
#=
Created on Mon Oct 17 2022 

@author Mason Rogers

rc_sde_init.jl defines the constants and functions required on every processor
used in a multiprocessor ensemble solution of SDEs.
=#

#imports
using Pkg
Pkg.activate(".")
using DifferentialEquations, Printf, Distributed

##########   PARAMETERS   ############################################
#read parameters
io = open("rc_param.py", "r")
param_str_list = readlines(io)
close(io)
for param_str in param_str_list
    param_expr = Meta.parse(replace(param_str,"**"=>"^"))
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
# const ϵ = ((1+2*B)*d^2*Us)/(36*ν*Ls) #should be small
# const C = (2*g*(B-1)*Ls*ϵ)/((1+2*B)*Us^2) #should be O(1)
const A = κ/ϵ/Us/Ls
@printf "A: %.6f\n" A
@printf "C: %.6f\n" C
@printf "ϵ: %.6f\n" ϵ
@printf "κ: %.6f\n" κ

#6d noise magnitude
const α = sqrt(2*Us^3*A/Ls/ϵ)

##########   FLUID DYNAMICS   ########################################
#fluid velocities
function fluid_vel(t, x, y, z)
    ρ = sqrt(x^2 + y^2)
    u = -b*x*(1-2*z)*(R-ρ)/3 - a*y*(c+z^2) + δ*(1-β*z)*(y*(y-ys+γ*cos(σ*t)) - (R^2-ρ^2)/2)
    v = -b*y*(1-2*z)*(R-ρ)/3 + a*x*(c+z^2) - δ*(1-β*z)*x*(y-ys+γ*cos(σ*t))
    w = b*z*(1-z)*(2*R/3-ρ)    
    return [u, v, w]
end

#t derivatives
function fluid_vel_t(t, x, y, z)
    ρ = sqrt(x^2 + y^2)
    ut = -δ*(1-β*z)*y*σ*γ*sin(σ*t)
    vt = δ*(1-β*z)*x*σ*γ*sin(σ*t)
    wt = 0
    return [ut, vt, wt]
end

#x derivatives
function fluid_vel_x(t, x, y, z)
    ρ = sqrt(x^2 + y^2)
    ux = b*(1-2*z)*x^2/(3*ρ) - b*(1-2*z)*(R-ρ)/3 + δ*(1-β*z)*x
    vx = b*(1-2*z)*x*y/(3*ρ) + a*(c+z^2) - δ*(1-β*z)*(y-ys+γ*cos(σ*t))
    wx = -b*(1-z)*z/ρ * x
    return [ux, vx, wx]
end

#y derivatives
function fluid_vel_y(t, x, y, z)
    ρ = sqrt(x^2 + y^2)
    uy = b*(1-2*z)*x*y/(3*ρ) - a*(c+z^2) + δ*(1-β*z)*(3*y-ys+γ*cos(σ*t))
    vy = b*(1-2*z)*y^2/(3*ρ) - b*(1-2*z)*(R-ρ)/3 - δ*(1-β*z)*x
    wy = -b*(1-z)*z/ρ * y
    return [uy, vy, wy]
end

#z derivatives
function fluid_vel_z(t, x, y, z)
    ρ = sqrt(x^2 + y^2)
    uz = 2*b*x*(R-ρ)/3 - 2*a*y*z - δ*β*(y*(y-ys+γ*cos(σ*t)) - (R^2-ρ^2)/2)
    vz = 2*b*y*(R-ρ)/3 + 2*a*x*z + δ*β*x*(y-ys+γ*cos(σ*t))
    wz = b*(1-2*z)*(2*R/3-ρ)
    return [uz, vz, wz]
end

#material fluid acceleration
function fluid_acc(t, x, y, z)
    u, v, w = fluid_vel(t, x, y, z)
    ut, vt, wt = fluid_vel_t(t, x, y, z)
    ux, vx, wx = fluid_vel_x(t, x, y, z)
    uy, vy, wy = fluid_vel_y(t, x, y, z)
    uz, vz, wz = fluid_vel_z(t, x, y, z)
    au = ut + u*ux + v*uy + w*uz
    av = vt + u*vx + v*vy + w*vz
    aw = wt + u*wx + v*wy + w*wz
    return [au, av, aw]
end

#3d particle velocity
function parti_vel(t, x, y, z)
    uᶠ, vᶠ, wᶠ = fluid_vel(t, x, y, z)
    au, av, aw = fluid_acc(t, x, y, z)
    uz, vz, wz = fluid_vel_z(t, x, y, z)
    u = uᶠ + (ϵ*Ls/Us) * (2*(1-B)/(1+2*B)*au + Us*C*uz)
    v = vᶠ + (ϵ*Ls/Us) * (2*(1-B)/(1+2*B)*av + Us*C*vz)
    w = wᶠ + (ϵ*Ls/Us) * (2*(1-B)/(1+2*B)*aw + Us*C*wz)
    return [u, v, w]
end

##########   DIFFERENTIAL EQUATIONS   ################################
#6d deterministic equations
function mre_det_6d!(ξ̇, ξ, q, t)
    x, y, z, u, v, w = ξ
    uᶠ, vᶠ, wᶠ = fluid_vel(t, x, y, z)
    au, av, aw = fluid_acc(t, x, y, z)
    ξ̇[1] = u
    ξ̇[2] = v
    ξ̇[3] = w
    ξ̇[4] = 3/(1+2*B)*au + (Ls/(Us*ϵ))*(uᶠ-u)
    ξ̇[5] = 3/(1+2*B)*av + (Ls/(Us*ϵ))*(vᶠ-v)
    ξ̇[6] = 3/(1+2*B)*aw + (Ls/(Us*ϵ))*(wᶠ-w-C*Us)
end

#4d stochastic terms
function mre_sto_6d!(ξ̇, ξ, q, t)
    ξ̇[1] = 0
    ξ̇[2] = 0
    ξ̇[3] = 0
    ξ̇[4] = α
    ξ̇[5] = α
    ξ̇[6] = α
end

#3d deterministic equations
function mre_det_3d!(ξ̇, ξ, q, t)
    x, y, z = ξ
    u, v, w = parti_vel(t, x, y, z)
    ξ̇[1] = u
    ξ̇[2] = v
    ξ̇[3] = w
end

#3d stochastic terms
function mre_sto_3d!(ξ̇, ξ, q, t)
    ξ̇[1] = sqrt(2*κ)
    ξ̇[2] = sqrt(2*κ)
    ξ̇[3] = sqrt(2*κ)
end

##########   INITIAL AND BOUNDARY CONDITIONS   ########################
#6d random ensemble initial conditions
function rand_ic_6d!(p, i, r)
    x₁ = x₀ + cbrt(Σ)*randn(3)
    u₁ = fluid_vel(0, x₁...)
    p.u0 .= vcat(x₁, u₁)
    return p
end

function rand_ic_3d!(p, i, r)
    x₁ = x₀ + cbrt(Σ)*randn(3)
    p.u0 .= x₁
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

#particle exits lateral wall
function out_of_domain_r(ξ, t, integrator)
    ρ = sqrt(ξ[1]^2+ξ[2]^2)
    return ρ-R
end

#particle exits top or bottom
function out_of_domain_z(ξ, t, integrator)
    z = ξ[3]
    return z*(z-H)
end

#particle reflects off lateral boundary
function reflect_r!(integrator)
    #x, y, z, u, v, w = integrator.u
    #sn = (u*x + v*y)/sqrt(x^2 + y^2)
    #integrator.u[4] -= 2*sn*x/sqrt(x^2+y^2)
    #integrator.u[5] -= 2*sn*y/sqrt(x^2+y^2)
end
cb_out_r = ContinuousCallback(out_of_domain_r, reflect_r!, save_positions=(false,false))

#particle reflects off top or bottom boundary
function reflect_z!(integrator)
    #x, y, z, u, v, w = integrator.u
    #integrator.u[6] *= -1 
end
cb_out_z = ContinuousCallback(out_of_domain_z, reflect_z!, save_positions=(false,false))

#package callbacks
cb_set = CallbackSet(cb_out_r, cb_out_z)
