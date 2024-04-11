using Pkg
Pkg.activate(".")
using Random, StatsBase, Printf
include("kv_param.py")

function fluid_vel(t, x, y)
    ρ = sqrt(x^2 + y^2)
    θ = atan(y, x)
    ω = Γ/2/π/(ρ^2+a^2)
    u = -ρ*ω*sin(θ)
    v = ρ*ω*cos(θ)
    return [u, v]
end

function find_rng_seed(n)
    tol = 0.000012
    
    pred = 2*Σ + x0^2 + y0^2
    rng = Xoshiro()
    nx = Int(round(2*R/dx)) + 4
    ny = Int(round(2*R/dy)) + 4
    XG = dx .* collect(-nx/2:nx/2-1)
    YG = dy .* collect(-ny/2:ny/2-1)
    XC = ones(size(YG)) * (XG .+ dx/2)'
    YC = (YG .+ dy/2) * ones(size(XG))'
    edges = (vcat(XG, maximum(XG)+dx), vcat(YG, maximum(YG)+dy))
    
    rv = zeros(2, n)
    seed = 0
    err = tol+1
    while err >= tol 
        seed += 1
        Random.seed!(rng, seed)
        rv = (sqrt(Σ) * randn(rng, (2, n))) .+ [x0, y0]
        p = fit(Histogram, (rv[1,:], rv[2,:]), edges).weights / n
        v = sum(p .* (XC.^2 + YC.^2))
        
        err = sum(p) == 1 ? abs(v - pred) : tol+1
        @printf "seed: %d, error: %.5f \n" seed err
    end

    return rv
end

function write_files(rv)
    dir = "../output/exp1/"
    fn2d = dir*"julia2dx.0000000000.bin"
    fn4d = dir*"julia4dx.0000000000.bin"

    vd = reduce(hcat, fluid_vel.(zeros(n), rv[1,:], rv[2,:]))
    ic = vcat(rv, vd)

    io = open(fn2d, "w")
    write(io, rv)
    close(io) 
    io = open(fn4d, "w")
    write(io, ic)
    close(io)

    println(ic[:, 1:10])
    println()
    println(rv[:, 1:10])
end 

n = 500000
rv = find_rng_seed(n)
write_files(rv)
