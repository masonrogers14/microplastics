#!/usr/bin/env julia
# -*- coding utf-8 -*-
#=
Created on Thu Apr 18 2024

@author Mason Rogers

new_kv_sde_param.jl reads in parameters for solving trajectories.
=#

#imports
using Pkg
Pkg.activate(".")
using Printf

##########   PARAMETERS   ######################################################
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

#grid parameters
const nx = Int(round(2*R/dx)) + 4
const ny = Int(round(2*R/dy)) + 4
const nz = 1
const vol = dx*dy*dz
const XG = dx .* collect(-nx/2:nx/2-1)
const XC = XG .+ dx/2
const YG = dy .* collect(-ny/2:ny/2-1)
const YC = YG .+ dy/2
const RF = dz .* collect(nz:-1:0) #length nz+1
const RC = (RF[1:end-1] .+ RF[2:end]) ./ 2

#histogram bin edges and flips
const edges = (vcat(XG, maximum(XG)+dx), vcat(YG, maximum(YG)+dy))
const flip_last = false

#simulation parameters
const nPerProc = Int(nTraj / nprocs())
const nOuts = Int(floor(tStop/wFreq) + 1)