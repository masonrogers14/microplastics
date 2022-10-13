import os
import numpy as np

for j in np.arange(0, 1000000, 10000):
    no = '{0:010d}'.format(j)
    nn = '{0:010d}'.format((j+1)//10)
    fo = 'mitgcm2d.{0}.data'.format(no)
    fn = 'mitgcm2d.{0}.data'.format(nn)
    os.system('mv {0} {1}'.format(fo, fn))
    fo = 'mitgcm2d.{0}.meta'.format(no)
    fn = 'mitgcm2d.{0}.meta'.format(nn)
    os.system('mv {0} {1}'.format(fo, fn))

#for j in np.arange(0, 100001, 1000):
#    fo = 'mitgcm2dl.{0:010d}.data'.format(j)
#    fn = 'mitgcm2d.{0:010d}.data'.format(j)
#    os.system('mv {0} {1}'.format(fo, fn))
#    fo = 'mitgcm2dl.{0:010d}.meta'.format(j)
#    fn = 'mitgcm2d.{0:010d}.meta'.format(j)
#    os.system('mv {0} {1}'.format(fo, fn))
# 
#
