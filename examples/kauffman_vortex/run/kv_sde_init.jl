using Pkg
Pkg.activate(".")
using DifferentialEquations

#Maxey-Riley parameters
const B = 0.98 #density ratio ρ/ρᶠ
const ν = 1e-6 #fluid viscosity
const d = 1e-3 #particle diameter
const g = 9.8 #gravity
const Us = 1.0 #characteristic velocity scale of flow
const Ls = 1.0 #characteristic length scale of flow

#dimensional 2d diffusivity
const κ = 4e-5

#nondimensional parameters
const ϵ = ((1+2*B)*d^2*Us)/(36*ν*Ls) #should be small
const C = (2*g*(B-1)*Ls*ϵ)/((1+2*B)*Us^2) #should be O(1)
const A = κ/ϵ/Us/Ls

#4d noise magnitude
const α = sqrt(2*Us^3*A/Ls/ϵ)

using Printf
@printf "A: %.6f\n" A
@printf "C: %.6f\n" C
@printf "ϵ: %.6f\n" ϵ
@printf "κ: %.6f\n" κ

#Rankine vortex parameters
const R = .5
const a = .2
const Γ = 1.

#output grid parameters
const dx = .01
const dy = .01
const dz = 1.0
const nx = Int(round(2*R/dx)) + 4
const ny = Int(round(2*R/dy)) + 4
const nz = 1
const dt = 1e-3
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
#To do: implement function parti_vel(t, x, y)

#4d deterministic equations
function mre_det_4d!(ξ̇, ξ, q, t)
    x, y, u, v = ξ
    uᶠ, vᶠ = fluid_vel(t, x, y)
    au, av = fluid_acc(t, x, y)
    ξ̇[1] = u
    ξ̇[2] = v
    ξ̇[3] = 3/(1+2*B)*au + (36*ν)/((1+2*B)*d^2)*(uᶠ-u)
    ξ̇[4] = 3/(1+2*B)*av + (36*ν)/((1+2*B)*d^2)*(vᶠ-v)
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

#modify ensemble initial conditions
function change_ic!(p, i, r)
    Σ = (401/4)*dx^2
    x₁ = x₀ + sqrt(Σ)*randn(2)
    u₁ = fluid_vel(0, x₁...)
    p.u0 .= vcat(x₁, u₁)
    return p
end

#y=0 event
# function y_axis(ξ,t,integrator)
#     x, y, u, v = ξ
#     return y
# end
# cb_y_axis = ContinuousCallback(y_axis, (integrator -> integrator), save_positions=(true,false))

#particle exits domain event
function out(ξ,t,integrator)
    x, y, u, v = ξ 
    ρ = sqrt(x^2+y^2)
    return ρ-R
end
cb_out = ContinuousCallback(out, terminate!, save_positions=(false,false))

#package callbacks
cb_set = CallbackSet(cb_out)
