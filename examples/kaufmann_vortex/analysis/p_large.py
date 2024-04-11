#Maxey-Riley parameters
B = 0.98 #density ratio ρ/ρᶠ
ν = 1e-6 #fluid viscosity
d = 1e-4 #particle diameter
g = 9.8 #gravity
Us = 1.0 #characteristic velocity scale of flow
Ls = 1.0 #characteristic length scale of flow

#stochastic parameters
κ = 4e-5

#Kaufmann vortex parameters
R = 5.
a = .2
Γ = 5.

#solver parameters
#tStop = 200.0
#wFreq = 0.25
tStop = 1000.0
wFreq = 10.0

#output grid parameters
dx = .01
dy = .01
dz = 1.0
dt = 1e-3
