import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.animation as movie

dates = np.arange('2012-06-07', '2012-08-01', dtype='datetime64[D]')
files = ['gs_{0:02d}{1:02d}.nc'.format(d.item().month, d.item().day) for d in dates]
gr = xr.open_dataset('gs_grid.nc')

f, a = plt.subplots(figsize=(8,8), constrained_layout=True)
global p
x = gr['XG'].isel(i_g=slice(None,-1), j_g=slice(None, -1))
y = gr['YG'].isel(i_g=slice(None,-1), j_g=slice(None,-1))
dx = gr['DXV']
dy = gr['DYU']

p = plt.pcolormesh(x, y, np.zeros(x.shape), vmin=-2e-5, vmax=2e-5, cmap='PiYG')
plt.scatter(-72.9, 36.55)
title = plt.title(dates[0])

def make_movie(i):
    print(i)
    ds = xr.open_dataset(i)
    u = ds['U'].isel(k=0)
    v = ds['V'].isel(k=0)
    v_x = v.diff(dim='i').values / dy.isel(i=slice(None,-1)).values
    u_y = u.diff(dim='j').values / dx.isel(j=slice(None,-1)).values
    ζ = v_x[:-1,:] - u_y[:,:-1]
    
    p.set_array(ζ)

    title.set_text(str(i))
    return [p,title]

m = movie.FuncAnimation(f, make_movie, frames=files, blit=True)
Writer = movie.writers['ffmpeg_file']
writer = Writer(fps=3, metadata=dict(artist='Mason'), bitrate=1500)
m.save('zeta_movie.mp4', writer=writer)
