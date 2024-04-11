#Maxey-Riley parameters
B = 0.8 #density ratio ρ/ρᶠ
ν = 1e-6 #fluid viscosity
d = 1e-4 #particle diameter
g = 9.8 #gravity
Us = 1.0 #characteristic velocity scale of flow
Ls = 1.0 #characteristic length scale of flow

#stochastic parameters
κ = 4e-5

#Kaufmann vortex parameters
R = 2.
a = 1.
Γ = 50.

#solver parameters
tStop = 100.0
wFreq = 1.0

#output grid parameters
dx = .01
dy = .01
dz = 1.0
dt = 1e-3

#initial conditions
Σ = .0005
x0 = 0.
y0 = 0.
