#Maxey-Riley parameters
B = 0.98 #density ratio ρ/ρᶠ
ν = 1e-6 #fluid viscosity
d = 1e-4 #particle diameter
g = 9.8 #gravity
Us = 1.0 #characteristic velocity scale of flow
Ls = 1.0 #characteristic length scale of flow

#stochastic parameters
κ = 4e-5

#Kauffman vortex parameters
#R = .5
R = 4.
a = .5
Γ = 8.

#solver parameters
#tStop = 200.0
#wFreq = 0.25
tStop = 4000.0
wFreq = 10.0

#output grid parameters
dx = .01
dy = .01
dz = 1.0
dt = 1e-3
