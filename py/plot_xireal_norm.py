# Plot real space correlation function and bias as a function of scale


import sys
sys.path.append('/home/CEFCA/aaorsi/')
import matplotlib
matplotlib.use('Agg')

import pylab as pl
import matplotlib.gridspec as gsc
import matplotlib.colors as colors
import matplotlib.cm as cmx

import cosmocodes.modelling.rsd.model as model
from cosmocodes.modelling.rsd.get_pars import *
import cosmocodes.io as cio
from scipy.interpolate import interp1d


sz = '0.0'
pl.rc('font', family='STIXGeneral')

datatype = ['Real','Monopole','Quadrupole']
ndtype = len(datatype)

ndensarr = ['1e-3','5e-4','1e-4']
nstrings = [r'10^{-3}',r'5\times 10^{-4}',r'10^{-4}']
ndens = len(ndensarr)
galtypearr = ['Starforming','Mstellar']#,'OII_3727']

#datadir_sfr = '/home/CEFCA/aaorsi/cosmocodes/data/datatest/'
datadir = '/home/CEFCA/aaorsi/data/'
mcmap = pl.get_cmap('bone')
scmap = pl.get_cmap('cool')

cNorm = colors.Normalize(vmin = -1, vmax=ndens+1)

mscalmap = cmx.ScalarMappable(norm=cNorm,cmap=mcmap)
sscalmap = cmx.ScalarMappable(norm=cNorm,cmap=scmap)

#pl.figure('xireal')
#gs = gsc.GridSpec(1,6)

f, axarr = pl.subplots(3,6)
f.subplots_adjust(wspace=0.0,hspace=0.0,right=0.85)

xlin_lim  = [20,99]
xlog_lim  = [0.1,20]

xilim = [[-5,23],[-5,24],[-15,44]]
blim  = [1,3.9] 

sfcol = plt.cm.Blues(np.linspace(0.3, 1, 4))
mscol = plt.cm.Oranges(np.linspace(0.3, 1, 4))
o2col = plt.cm.Greens(np.linspace(0.3,1,4))

#gs.update(wspace=0,hspace=0,top=0.5)

lstyle = ['-','--',':']
msize = 5

for nd in range(ndens):
  xidata = {}
  for gtype in galtypearr:
    col = sfcol if gtype == 'Starforming' else mscol
    if gtype == 'OII_3727':
      col = o2col
    for idt in range(ndtype):
      ctype = datatype[idt]
      #fname = gtype + 'galdata_sz_'+sz+'_ndens'+ndensarr[nd]+'_bin0_rspacedm_' if ctype == 'real' else  gtype + 'galdata_sz_'+sz+'_zspace_ndens'+ndensarr[nd]+'_bin0_rspacedm_'
      fname = gtype + 'galdata_sz_'+sz+'_ndens'+ndensarr[nd]+'_bin0_rspacedm_' if ctype == 'Real' else gtype+'galdata_sz_'+sz+'_zspace_typeNone_q02.8e7_g0-1.3sfog20.0_ndens'+ndensarr[nd]+'_bin0_rspacedm_'


      xilog = cio.load_correlation_functions(datadir,fname,outdata='xi',logbins=True)
      xilin = cio.load_correlation_functions(datadir,fname,outdata='xi',logbins=False)
      rlog = xilog['r']
      rlin = xilin['r']
      
      if idt > 0:
        xilin.update(cio.load_correlation_functions(datadir,fname,outdata='xi2d'))
      
        clim = xlin_lim[0]
        clog = np.where(rlog < clim)
        clin = np.where(rlin >= clim)
        
        rr = np.concatenate((rlog[clog],rlin[clin]))
        xilin.update(xiutils.build_multipoles_sigpi(rr,xilin,r_int=xilin['r']))
        xileg = xilin['xileg']
        ss    = xilin['s']

      _str = ['1.0', '40', '80',sz,ndensarr[nd],gtype,'Monopole']
      _cosmo, _model, _data = get_pars(_str)

      bias = _cosmo['bias_1']
      fg  =  _cosmo['growthrate']
      beta = fg/bias

      if idt == 0:
        den = bias**2
      if idt == 1:
        den = bias**2 * (1 + (2./3)*beta + (1./5)*beta**2)
      if idt == 2:
        den = bias**2 * (4./3*beta + 4./7*beta**2)
      
  #    xidata.update(gtype:{'xilin':xilin,'xilog':xilog,'r':r,'bias':bias})
     
      sig = 1 if idt == 1 else -1
      ax = axarr[idt,2 * nd]
      if idt == 0:
        ax.semilogx(rlog,rlog**2*xilog['xi']/den,'-',color=col[nd],linewidth=2,label=gtype,markersize=msize)
      else:
        ax.semilogx(ss,sig*ss**2*xileg[idt-1,:]/den,'-',color=col[nd],linewidth=2,label=gtype,markersize=msize)
        
      ax.set_ylim(xilim[idt])
      ax.set_xlim(xlog_lim)
      ax.fill_between(xlog_lim,[xilim[idt][0],xilim[idt][0]],[xilim[idt][1],xilim[idt][1]],color='gray',alpha=0.2)
      ax.set_xticks([0.1,1.0,10.0])
      if nd == 0:
        ax.yaxis.set_ticks_position('left')
      else:
        ax.set_yticks([])

      if nd == 0 and idt == 0:
        ax.text(xlog_lim[0]*0.02,xilim[idt][1]*0.15,r'$r^2\xi(r)/\mathcal{K}(\beta) [({\rm Mpc}/h)^2]$',fontsize=15,rotation=90)

      if idt < ndtype-1:
        ax.set_xticks([])

      ax.spines['right'].set_visible(False)

      cc = np.where(rlin > xlin_lim[0])
      ax = axarr[idt,2 * nd+1]
      if idt ==0:
        ax.plot(rlin[cc],rlin[cc]**2*xilin['xi'][cc]/den,
        '-',color=col[nd],linewidth=2,markersize=msize)
      else:
        ax.plot(ss,sig*ss**2*xileg[idt-1,:]/den,
        '-',color=col[nd],linewidth=2,markersize=msize)
        
      ax.set_ylim(xilim[idt])
      ax.set_xlim(xlin_lim)
      ax.spines['left'].set_visible(False)
      ax.yaxis.set_ticks_position('right')
      ax.set_yticks([])
      ax.set_xticks([50,70])
      if idt == 2 and nd == 1:
        ax.text(-0.1,-0.2,r'$r[{\rm Mpc/}h]$',fontsize=15,transform=ax.transAxes)
      
      if idt == 0:
        ax.text(xlin_lim[0]-30,xilim[idt][1]+2,r'$n = '+nstrings[nd]+'$',fontsize=15)
    
      if nd == 2:
        ax.text(.5,.15,datatype[idt],fontsize=7,transform=ax.transAxes)

#      if nd == 1 and idt == 0:
#        ax.legend(prop={'size':10},loc='upper center',bbox_to_anchor=(0.5,1.05),
#        ncol = 2,fancybox=True,shadow=True)
#        print 'leg here'
    
  
pl.savefig('../plots/xinorm.png',bbox_inches='tight',dpi=300)
#pl.show()
#pl.close("all")


