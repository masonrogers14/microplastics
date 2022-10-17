#Maxey-Riley parameters
B = 0.98 #density ratio ρ/ρᶠ
ν = 1e-6 #fluid viscosity
d = 1e-3 #particle diameter
g = 9.8 #gravity
Us = 1 #characteristic velocity scale of flow
Ls = 1 #characteristic length scale of flow
ϵ = ((1+2*B)*d**2*Us)/(36*ν*Ls) #should be small
C = (2*g*(B-1)*Ls*ϵ)/((1+2*B)*Us**2) #should be O(1)

#stochastic parameters
κ = 4e-5

#rotating can parameters
R = 0.5 
H = 1.0
a = .62
b = 7.5
c = .7
β = 1
γ = .2
δ = .2
ys = -.2
σ = 0

#solver parameters
tStop = 100.0
wFreq = 1.0

#output grid parameters
dx = .01
dy = .01
dz = .02
dt = 1e-3

#initial conditions
Σ = 1e-5
x0 = 0.33
y0 = 0.
z0 = 0.56
