import os
import numpy as np
import xarray as xr

for d in np.arange('2012-07-01', '2012-09-01', dtype='datetime64[D]'):
    print(d, end=': ')
    
    #define url
    mo = d.item().month
    da = d.item().day
    url = 'https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-protected/MITgcm_LLC4320_Pre-SWOT_JPL_L4_WestAtlantic_v1.0/LLC4320_pre-SWOT_WestAtlantic_2012{0:02d}{1:02d}.nc'.format(mo,da)
    fnOld = 'LLC4320_pre-SWOT_WestAtlantic_2012{0:02d}{1:02d}.nc'.format(mo,da)
    fnNew = 'gs_{0:02d}{1:02d}.nc'.format(mo,da)

    #get contents
    print('getting data', end=' ')
    os.system('wget {0}'.format(url))

    #open dataset and extract contents
    print('-> extracting data', end=' ')
    dsOld = xr.open_dataset(fnOld)
    dsNew = xr.Dataset()
    for k in ['U', 'V', 'W', 'Eta']:
        dsNew[k] = dsOld[k].mean(dim='time')

    #write new dataset, delete old
    print('-> writing data')
    dsNew.to_netcdf(fnNew)
    os.system('rm {0}'.format(fnOld))
