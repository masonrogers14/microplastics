#!/usr/bin/env julia
# -*- coding utf-8 -*-
#=
Created on Thu Oct 21 2021

@author Mason Rogers
=#

#imports
using Pkg
Pkg.activate(".")
using DifferentialEquations, Printf
include("../traj/rc_sde_init.jl")

#tinker
tMax = 10000
dir = @sprintf "eps%02d/" Int(round(δ*100))
θ_fname = dir*"att_theta.bin"
z_fname = dir*"att_z.bin"
r_fname = dir*"att_r.bin"

x₀ = [.33, 0., .56] #should be very close to the attractor!
u₀ = fluid_vel(0., x₀...)
ξ₀ = vcat(x₀, u₀)

prob = ODEProblem(mre_det_6d!, ξ₀, (0.0,tMax), save_everystep=false, save_end=false)
solu = solve(prob, Tsit5(), saveat=99*tMax/100:.1:tMax)

x = [solu.u[i][1] for i in 1:length(solu)]
y = [solu.u[i][2] for i in 1:length(solu)]
z = [solu.u[i][3] for i in 1:length(solu)]
θ = mod.(atan.(y,x), 2π)
r = sqrt.(x.^2+y.^2)

i = findall(diff(θ) .< 0)[end-1:end]
i[1] += 1
z = z[i[1]:i[2]]
θ = θ[i[1]:i[2]]
r = r[i[1]:i[2]]

io = open(θ_fname, "w")
write(io, θ)
close(io)
io = open(z_fname, "w")
write(io, z)
close(io)
io = open(r_fname, "w")
write(io, r)
close(io)

