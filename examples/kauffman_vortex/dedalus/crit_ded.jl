using Printf
include("kv_param.py")

C = 3 #safety factor
ϵ = ((1+2*B)*d^2*Us)/(36*ν*Ls)

nr = C*R/sqrt(Σ)
nθ = 2*π*C*R/sqrt(Σ)
nρ = max(nr, nθ)

vθ = Γ/2/π * R/(nρ*a^2) 
vr = (Γ/4/π/a)^2 / a * (ϵ*Ls/Us) * 2*(1-B) / (1+2*B)

Δt = zeros(4)
Δt[1] = R/(C*nρ*vr)
Δt[2] = 2*π*R/(C*vθ*nθ*nρ)
Δt[3] = (R/nρ)^2/(C*κ)
Δt[4] = (2*π*R/(nρ*nθ))^2/(C*κ) 

@printf "nr > %d\n" nr
@printf "nθ > %d\n" nθ
@printf "Δt < %.6f, %.6f, %.6f, %.6f\n" Δt...
