using Pkg
Pkg.activate(".")
include("rc_sde_init.jl")
   
function run_test()
    #choose random points uniformly distributed through can
    n = 100
    r = R*sqrt.(rand(n))
    θ = 2π*rand(n)
    z = H*rand(n)
    x = r.*cos.(θ)
    y = r.*sin.(θ)
    t = zeros(n)
    
    #preallocate
    uᶠ_o = zeros(n,3)
    u̇ᶠ_o = zeros(n,3)
    uᶠ_n = zeros(n,3)
    ∂tuᶠ_n = zeros(n,3)
    ∇uᶠ_n = zeros(n,3,3)
    u̇ᶠ_n = zeros(n,3)
    
    #fluid velocity and material derivatives
    @time for j in 1:n
        uᶠ_o[j,:] = fluid_vel(t[j], x[j], y[j], z[j])
        u̇ᶠ_o[j,:] = fluid_acc(t[j], x[j], y[j], z[j])
    end
    
    @time for j in 1:n
        uᶠ_n[j,:], ∂tuᶠ_n[j,:], ∇uᶠ_n[j,:,:] = uf_and_derivatives(t[j], x[j], y[j], z[j])
        u̇ᶠ_n[j,:] = ∂tuᶠ_n[j,:] + ∇uᶠ_n[j,:,:]*uᶠ_n[j,:]
    end

    Δuᶠ = uᶠ_o - uᶠ_n
    Δu̇ᶠ = u̇ᶠ_o - u̇ᶠ_n

    println(maximum(abs.(Δuᶠ)))
    println(maximum(abs.(Δu̇ᶠ)))
end

