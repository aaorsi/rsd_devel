# A script to plot the correlation between sigmav and the host halo mass

import numpy as np
import pylab as pl

datadir = '/home/CEFCA/aaorsi/data/GalCat/MXXL_extended/'

Gtype = ['Mstellar','Starforming']
ndensarr = ['1e-3','5e-4','1e-4']
sz = '0.0'

ng = len(Gtype)
ndens = len(ndensarr)




for igt in range(ng):
  for nd in range(ndens):
    
    fname = '%s%s_sz_%s_q02.8e7_g0-1.3_ndens%s.sigma_M' % (datadir, Gtype[igt], sz, ndens[nd])
    data = np.loadtxt(fname)

    

