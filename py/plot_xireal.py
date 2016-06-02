# Plot real space correlation function and bias as a function of scale


import sys
sys.path.append('/home/CEFCA/aaorsi/')
import pylab as pl
import matplotlib.gridspec as gsc
import matplotlib.colors as colors
import matplotlib.cm as cmx

import cosmocodes.modelling.rsd.model as model
from cosmocodes.modelling.rsd.get_pars import *
import cosmocodes.io as cio
from scipy.interpolate import interp1d


sz = '0.0'
pl.rc('font', family='serif')

ndensarr = ['1e-3','5e-4','1e-4']
nstrings = [r'10^{-3}',r'5\times 10^{-4}',r'10^{-4}']
ndens = len(ndensarr)
galtypearr = ['Starforming','Mstellar']

datadir_sfr = '/home/CEFCA/aaorsi/cosmocodes/data/datatest/'
datadir_mstar = '/home/CEFCA/aaorsi/data/'


args = ['','1.0','1','150','0.0','Monopole','Quadrupole','Hexadecapole','PseudoQuadrupole']


cosmo_pars,model_pars,data_pars = get_pars(args[1:])
model_pars['inputpkfile'] = '/home/CEFCA/aaorsi/cosmocodes/data/wmp1_matterpower.dat'



modelfn = model.set_model_function('xilinear', model_pars)


mcmap = pl.get_cmap('bone')
scmap = pl.get_cmap('cool')

cNorm = colors.Normalize(vmin = -1, vmax=ndens+1)

mscalmap = cmx.ScalarMappable(norm=cNorm,cmap=mcmap)
sscalmap = cmx.ScalarMappable(norm=cNorm,cmap=scmap)

pl.figure('xireal')

gs = gsc.GridSpec(3,1)

xlim  = [-1,2.3]
xilim = [-0.99,150]
blim  = [1,3.9] 

sfcol = plt.cm.Blues(np.linspace(0.3, 1, 4))
mscol = plt.cm.Oranges(np.linspace(0.3, 1, 4))

gs.update(wspace=0,hspace=0,right=0.75)

pl.subplot(gs[0:2,0])

lstyle = ['-','--',':']

for nd in range(ndens):
  
  sfrfname = 'Starforminggaldata_sz_'+sz+'_ndens'+ndensarr[nd]+'_bin0_rspacedm_'
  mstellarfname = 'Mstellargaldata_sz_'+sz+'_ndens'+ndensarr[nd]+'_bin0_rspacedm_'
  ximstellar = cio.load_correlation_functions(datadir_mstar,mstellarfname,outdata='xi',logbins=True)
  xisfr      = cio.load_correlation_functions(datadir_sfr,sfrfname,outdata='xi',logbins=True)
  r = ximstellar['r']
  
  if nd == 0:
    ximodel = modelfn(xisfr['r'], cosmo_pars,data_pars=data_pars)
    xidm = interp1d(ximodel['LinearTheory']['r'],ximodel['LinearTheory']['xi_lt'])  

    
  pl.plot(np.log10(r),r**2*ximstellar['xi'],'o',color=mscol[nd],linewidth=2,label=r'$M_{\star},n='+nstrings[nd]+'$')
  pl.plot(np.log10(r),r**2*xisfr['xi'],'o',color=sfcol[nd],linewidth=2,label=r'${\rm SF}, n='+nstrings[nd]+'$')

  pl.xlim(xlim)
  pl.ylim(xilim)
  pl.ylabel(r'$r^2\xi(r) [({\rm Mpc}/h)^2]$',fontsize=15)


pl.plot(np.log10(r),r**2*xidm(r),linewidth=3,color='gray')

#pl.text((xlim[1]-xlim[0])/10 * 7,(xilim[1]-xilim[0])/10 * 8,'Stellar mass',color='red',fontsize=15)
#pl.text((xlim[1]-xlim[0])/10 * 7,(xilim[1]-xilim[0])/10 * 7,'Star-forming',color='blue',fontsize=15)

pl.legend(prop={'size':10},loc='center left',bbox_to_anchor=(1.0,0.5),
          ncol = 1,fancybox=True,shadow=True)

pl.subplot(gs[2,0])


for nd in range(ndens):
  sfrfname = 'Starforminggaldata_sz_'+sz+'_ndens'+ndensarr[nd]+'_bin0_rspacedm_'
  mstellarfname = 'Mstellargaldata_sz_'+sz+'_ndens'+ndensarr[nd]+'_bin0_rspacedm_'
  ximstellar = cio.load_correlation_functions(datadir_mstar,mstellarfname,outdata='xi',logbins=True)
  xisfr      = cio.load_correlation_functions(datadir_sfr,sfrfname,outdata='xi',logbins=True)


  pl.plot(np.log10(r),np.sqrt(ximstellar['xi']/xidm(r)),'-',color=mscol[nd],linewidth=2)
  pl.plot(np.log10(r),np.sqrt(xisfr['xi']/xidm(r)),'-',color=sfcol[nd],linewidth=2)

  pl.xlim(xlim)
  pl.ylim(blim) 

  pl.ylabel(r'$b(r)$',fontsize=15)
  pl.xlabel(r'$\log(r [{\rm Mpc}/h])$',fontsize=15)

pl.draw()
  
pl.savefig('bias_comp.pdf',bbox_inches='tight')
#pl.show()
#pl.close("all")


fig2 = pl.figure('bias2')

brange = [30,60]
for nd in range(ndens):
  
  sfrfname = 'Starforminggaldata_sz_'+sz+'_ndens'+ndensarr[nd]+'_bin0_rspacedm_'
  mstellarfname = 'Mstellargaldata_sz_'+sz+'_ndens'+ndensarr[nd]+'_bin0_rspacedm_'
  ximstellar = cio.load_correlation_functions(datadir_mstar,mstellarfname,outdata='xi',logbins=True)
  xisfr      = cio.load_correlation_functions(datadir_sfr,sfrfname,outdata='xi',logbins=True)
  r = ximstellar['r']
  
  if nd == 0:
    ximodel = modelfn(xisfr['r'], cosmo_pars,data_pars=data_pars)
    xidm = interp1d(ximodel['LinearTheory']['r'],ximodel['LinearTheory']['xi_lt'])  

  xim = ximstellar['xi']
  xis = xisfr['xi']
  
  rrange = np.where((r > brange[0]) & (r < brange[1]))
  
  rsel = r[rrange[0]]

  bm = np.sqrt(np.mean(xim[rrange[0]]/xidm(rsel)))
  bs = np.sqrt(np.mean(xis[rrange[0]]/xidm(rsel)))

  pl.plot(np.log10(r),(np.sqrt(xim/xidm(r)) - bm)/bm,'o',color=mscol[nd])
  pl.plot(np.log10(r),(np.sqrt(xis/xidm(r)) - bs)/bs,'o',color=sfcol[nd])

  pl.xlabel(r'$\log(r[{\rm Mpc}/h])$',fontsize=20)
  pl.ylabel(r'$b_2(r)$',fontsize=20)

  pl.xlim([-2,2.0])
  pl.ylim([-2,30])

  pl.plot([-2,2],[0,0],'-',color='grey',linewidth=4)


pl.savefig('b2.pdf',bbox_inches='tight')



