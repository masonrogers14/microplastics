'''
subset_gs_grid.py finds the subset of LLC4320 points for which we have
velocity data.

assumptions:
1. the 
'''


#imports
import numpy as np
import xarray as xr

#init
check = False
tile = 10 #0-indexed
writeBin = True
outFile = '../run/grid.bin'
cutInd = {'i': [60, 100], 'j': [60, 100]}
outDict = {}
keyPairs = {'XC': 'XC', 'YC': 'YC', 'XG': 'XG', 'YG': 'YG',
            'DXC': 'DYC', 'DYC': 'DXC', 'DXG': 'DYG', 'DYG': 'DXG',
            'DXF': 'DYF', 'DYF': 'DXF', 'DXV': 'DYU', 'DYU': 'DXV',
            'RAC': 'RAC', 'RAW': 'RAS', 'RAS': 'RAW', 'RAZ': 'RAZ'}
keyShift = {'XC': False, 'YC': False, 'XG': True, 'YG': True,
            'DXC': False, 'DYC': True, 'DXG': True, 'DYG': False,
            'DXF': False, 'DYF': False, 'DXV': True, 'DYU': True,
            'RAC': False, 'RAW': False, 'RAS': True, 'RAZ': True}
if tile < 7:
    for k in keyShift.keys():
        keyShift[k] = False
ds = xr.open_dataset('gs_grid.nc')

#read data
def read_data(k):
    x = np.fromfile(k+'.data', '>f4')
    s = x.size
    ny = int(np.sqrt(s*13))
    nx = int(ny / 13)
    x = x.reshape(ny, nx)
    if tile < 7:
        x = x[tile*nx:(tile+1)*nx,:]
    else:
        if tile < 10:
            x = x[7*nx+(tile-7)%3:10*nx:3]
        else:
            x = x[10*nx+(tile-7)%3::3]
        x = np.rot90(x)
    return x 

#get matching index pattern
def get_match(k1, k2, x, y):
    if k1 in ds.data_vars and k2 in ds.data_vars:
        gx = ds[k1].values
        gy = ds[k2].values
        try:
            m0 = (gx[0,0] == x) & (gy[0,0] == y)
            m1 = (gx[-1,-1] == x) & (gy[-1,-1] == y)
            assert np.sum(m0) == 1
            assert np.sum(m1) == 1
            i0 = m0.max(axis=1).argmax()
            j0 = m0.max(axis=0).argmax()
            i1 = m1.max(axis=1).argmax() + 1
            j1 = m1.max(axis=0).argmax() + 1
        except AssertionError:
            print('match failed')
            i0 = np.nan
            j0 = np.nan
            i1 = np.nan
            j1 = np.nan
        finally:
            return i0, i1, j0, j1
    else:
        raise ValueError('variables used to get match don\'t exist in both datasets')

#check the validity of the match
def check_match(k):
    if k in ds.data_vars:
        return (ds[k].values == outDict[k]).all()
    else:
        return True

#get the matching indices from XC, YC
print('reading XC and YC')
outDict['XC'] = read_data('XC')
outDict['YC'] = read_data('YC')
print('getting match indices')
i0, i1, j0, j1 = get_match('XC','YC',outDict['XC'],outDict['YC'])
print('checking match')
outDict['XC'] = outDict['XC'][i0:i1,j0:j1]
outDict['YC'] = outDict['YC'][i0:i1,j0:j1]
print('XC: {0:d}'.format(check_match('XC')))
print('YC: {0:d}'.format(check_match('YC')))

#verify for everything else except the angles
for k, v in keyPairs.items():
    if not writeBin:
        print('match failed; halting without writing')
        break
    outDict[k] = read_data(v)
    if keyShift[k]:
        outDict[k] = outDict[k][i0-1:i1-1,j0:j1]
    else:
        outDict[k] = outDict[k][i0:i1,j0:j1]
    if check:
        writeBin = check_match(k)
    print('{0}: {1:d}'.format(k, writeBin))

#make angles
outDict['AngleCS'] = ds['AngleCS'].values
outDict['AngleSN'] = ds['AngleSN'].values

#write output
outOrder = ['XC', 'YC', 'DXF', 'DYF', 'RAC', 'XG', 'YG', 'DXV', 'DYU', 'RAZ', 'DXC', 'DYC', 'RAW', 'RAS', 'DXG', 'DYG', 'AngleCS', 'AngleSN']
if writeBin:
    outArray = np.stack([outDict[k][cutInd['j'][0]:cutInd['j'][1], cutInd['i'][0]:cutInd['i'][1]] for k in outOrder]) 
    outArray.astype('>f8').tofile(outFile)
