import numpy as np
from scipy.optimize import root

def tau_w0(d, B=0.999999, g=9.81, nu=1e-6):
    y = d**3 * (1-B) * g / 18 / nu**2
    Re = root(lambda Re: Re + 0.15*Re**1.687 - y, 800).x[0]
    if Re < 800:
        f = 1 + 0.15*Re**0.687
        tau = (1 + 2*B) * d**2 / 36 / nu / f
        w0 = d**2 * (1-B) * g / 18 / nu / f
        return tau, w0
    else:
        return np.nan, np.nan

for d in np.geomspace(1e-5, 1e-1, 100):
    tau, w0 = tau_w0(d)
    print('{0:.4e}, {1:.4e}'.format(tau, w0))