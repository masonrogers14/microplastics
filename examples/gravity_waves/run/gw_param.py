#Maxey-Riley parameters
B = 0.98 #density ratio ρ/ρᶠ
ν = 1e-6 #fluid viscosity
d = 1e-3 #particle diameter
g = 9.8 #gravity
Us = 1.0 #characteristic velocity scale of flow
Ls = 1.0 #characteristic length scale of flow

#stochastic parameters
κ = 4e-4

#wavefield parameters
Lx, Lz = (10., 10.)
nk = 1
η0 = .1

#solver parameters
tStop = 100.0
wFreq = 0.25

#output grid parameters
dx = .01
dy = 1.0
dz = 0.1
dt = 1e-3
