from os.path import expanduser
import numpy as np

import matplotlib
matplotlib.use('Agg')
import pylab as pl
import matplotlib.gridspec as gsc
import matplotlib.colors as colors
import matplotlib.cm as cmx

from scipy.optimize import curve_fit

plottype = 'vhalo'  # vhalo or sigma

home = expanduser("~")



#datadir = '/home/CEFCA/aaorsi/data/GalCat/MXXL_extended/'

if home == '/Users/aaorsi': #That means laptop copy
  datadir = '/Users/aaorsi/work/rsd_devel/data/cat/'
else:
  datadir = '/home/CEFCA/aaorsi/data/GalCat/MXXL_extended/'

Galtype = ['Starforming','Mstellar']#,'Halpha','OII_3727']
ngtype = len(Galtype)
sz = ['0.0']
ndens = ['1e-4','5e-4','1e-3']
nnd   = len(ndens)
zspace = [1]#,0]
zspace_type = ['None']
q0 = ['2.8e7']
g0 = ['-1.3']
sfog = ['20.0']
ngmax = 3e7
outdir = datadir

pl.rc('font', family='STIXGeneral')
gs = gsc.GridSpec(1,nnd)
gs.update(left=0.1,wspace=0.0,hspace=0.0,top=0.5)

axarr = []

for i in range(nnd):
  axarr.append(pl.subplot(gs[i]))

idens = 0

nmass = 5
nrbins = 10

marr = np.zeros(nmass)
rarr = np.zeros(nrbins)

svarr = np.zeros([nmass,nrbins])

col = pl.cm.Paired(np.linspace(0.0,0.75,nmass))

for igt in range(ngtype):
  gt = Galtype[igt]
      
  for _sz in sz:
    for idens in range(nnd):
      nd = ndens[idens]
      ax = axarr[idens]
      if idens == 0:
        ax.set_ylabel(r'$\log(\sigma_{\rm std})$',fontsize=15)
      else:
        ax.get_yaxis().set_ticks([])

      for zsp in zspace:
        if zsp == 1:
          for ztype in zspace_type:
            for _q0 in q0:
              for _g0 in g0:
                for s2 in sfog:
                  
                  outf = "%s%s_sz_%s_q0%s_g0%s_ndens%s.sigma_RM" % (datadir,gt,_sz,_q0,_g0,nd)
                   
                  data = np.loadtxt(outf)
                  j = 0
                  for imm in range(nmass):
                    marr[imm] = data[j,0]
                    for irr in range(nrbins):
                      if imm == 0:
                        rarr[irr] = data[j,1]
                      svarr[imm,irr] = data[j,2]
                      j += 1
                  

                  for imm in range(nmass):
                    ax.plot(rarr,svarr[imm,:],color=col[imm],
                    linewidth=2,label=gt+'M=%.2f' % (marr[imm]))
                    if igt == 0:
                      mhc = np.where(data[:,0] > 13.0)[0]
                      f = lambda x, *p: p[0] * x**(1./3)
                      coef,pcov = curve_fit(f,10**data[mhc,0],data[mhc,1],[100.0])
                      fx = f(10**data[:,0],coef)
                      ax.plot(data[:,0],fx,'--',color='black')
                    #ax.fill_between(data[:,0],np.log10(data[:,2]),np.log10(data[:,3]),color=col,alpha=0.25)
                    ax.set_xlabel(r'$\log(M_{\rm halo})$',fontsize=15)
                    ax.set_title('n = '+nd,fontsize=15)
                    if idens == 0:
                      ax.legend(loc='upper left')
                    
                    ax.set_title('n = '+nd,fontsize=15)
                    if idens == 0:
                      ax.legend(loc='upper left')
                      ax.set_ylabel(r'$\Delta V$',fontsize=15)



pl.savefig('../plots/sigmarm.pdf')

