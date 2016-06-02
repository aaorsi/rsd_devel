# Plot multipoles


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

cosmocode_dir = '/home/CEFCA/aaorsi/cosmocodes/modelling/rsd/'
sz = '0.0'
pl.rc('font', family='serif')

ndensarr = ['1e-3','5e-4','1e-4']
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


sfcol = plt.cm.Blues(np.linspace(0.3, 1, 4))
mscol = plt.cm.Oranges(np.linspace(0.3, 1, 4))


pl.figure()

gs = gsc.GridSpec(2,2)

xlim  = [9,89]
xilim = [20,150]
xilim2 = [0,69]
blim  = [1,2.9] 


lstyle = ['-','--']
gs.update(right=0.75,hspace=0,wspace=0)


ax0 = pl.subplot(gs[0,0])
ax2 = pl.subplot(gs[0,1])
ax4 = pl.subplot(gs[1,0])
axp2 = pl.subplot(gs[1,1])

for nd in range(ndens):
  
  sfrfname = 'Starforminggaldata_sz_'+sz+'_ndens'+ndensarr[nd]+'_bin0_rspacedm_'
  mstellarfname = 'Mstellargaldata_sz_'+sz+'_ndens'+ndensarr[nd]+'_bin0_rspacedm_'

  sfrfname_zspace = 'Starforminggaldata_sz_'+sz+'_zspace_ndens'+ndensarr[nd]+'_bin0_rspacedm_'
  mstellarfname_zspace = 'Mstellargaldata_sz_'+sz+'_zspace_ndens'+ndensarr[nd]+'_bin0_rspacedm_'

  # Read data first

  ximstellar = cio.load_correlation_functions(datadir_mstar,mstellarfname,outdata='xi')
  ximstellar.update(cio.load_correlation_functions(datadir_mstar,mstellarfname_zspace,outdata='xi2d'))
  rr = np.logspace(-3,2.5,num=200)
  ximstellar.update(xiutils.build_multipoles_sigpi(rr,ximstellar,r_int=ximstellar['r']))

  xisfr      = cio.load_correlation_functions(datadir_sfr,sfrfname,outdata='xi')
  xisfr.update(cio.load_correlation_functions(datadir_sfr,sfrfname_zspace,outdata='xiangle2d'))
  xisfr.update(cio.load_correlation_functions(datadir_sfr,sfrfname_zspace,outdata='xi2d'))
  xisfr.update(xiutils.build_multipoles_sigpi(rr,xisfr,r_int=xisfr['r']))

  # Compute models

  sf_str = ['1.0', '40', '80',sz,ndensarr[nd],'Starforming','Monopole']
  sf_cosmo, sf_model, sf_data = get_pars(sf_str) 

  sfile = cosmocode_dir + 'sigmas_' + sf_data['covname']
  sdata = np.loadtxt(sfile,skiprows=1)
  sf_cosmo['s12'] = sdata[0]
  sf_cosmo['s2fog'] = sdata[3]

  ms_str = ['1.0', '40', '80',sz,ndensarr[nd],'Mstellar','Monopole']
  ms_cosmo, ms_model, ms_data = get_pars(ms_str) 

  sfile = cosmocode_dir + 'sigmas_' + ms_data['covname']
  sdata = np.loadtxt(sfile,skiprows=1)
  ms_cosmo['s12'] = sdata[0]
  ms_cosmo['s2fog'] = sdata[3]
  
  modelf_sf = model.set_model_function('xilinear',sf_model)
  modelf_ms = model.set_model_function('xilinear',ms_model)
  
  ximodel_sf = modelf_sf(xisfr['r'],sf_cosmo,data_pars=sf_data)
  ximodel_ms = modelf_ms(ximstellar['r'],ms_cosmo,data_pars=ms_data)


  mdisp_sf_r   = ximodel_sf['Dispersion model']['s']
  mdisp_sf_xi0 = ximodel_sf['Dispersion model']['xileg'][0,:]
  mdisp_sf_xi2 = ximodel_sf['Dispersion model']['xileg'][1,:]
  mdisp_sf_xi4 = ximodel_sf['Dispersion model']['xileg'][2,:]
  mdisp_sf_pxi2 = ximodel_sf['Dispersion model']['xileg'][3,:]

  mdisp_ms_r   = ximodel_ms['Dispersion model']['s']
  mdisp_ms_xi0 = ximodel_ms['Dispersion model']['xileg'][0,:]
  mdisp_ms_xi2 = ximodel_ms['Dispersion model']['xileg'][1,:]
  mdisp_ms_xi4 = ximodel_ms['Dispersion model']['xileg'][2,:]
  mdisp_ms_pxi2 = ximodel_ms['Dispersion model']['xileg'][3,:]

  mgsm_sf_r   = ximodel_sf['GSM']['s']
  mgsm_sf_xi0 = ximodel_sf['GSM']['xileg'][0,:]
  mgsm_sf_xi2 = ximodel_sf['GSM']['xileg'][1,:]
  mgsm_sf_xi4 = ximodel_sf['GSM']['xileg'][2,:]
  mgsm_sf_pxi2 = ximodel_sf['GSM']['xileg'][3,:]

  mgsm_ms_r   = ximodel_ms['GSM']['s']
  mgsm_ms_xi0 = ximodel_ms['GSM']['xileg'][0,:]
  mgsm_ms_xi2 = ximodel_ms['GSM']['xileg'][1,:]
  mgsm_ms_xi4 = ximodel_ms['GSM']['xileg'][2,:]
  mgsm_ms_pxi2 = ximodel_ms['GSM']['xileg'][3,:]



  r = ximstellar['s']
  
  xi0_m =ximstellar['xileg'][0,:] 
  xi2_m =ximstellar['xileg'][1,:] 
  xi4_m =ximstellar['xileg'][2,:] 
  pxi2_m =ximstellar['xileg'][3,:]

  xi0_s =xisfr['xileg'][0,:] 
  xi2_s =xisfr['xileg'][1,:] 
  xi4_s =xisfr['xileg'][2,:] 
  pxi2_s =xisfr['xileg'][3,:]

  ax0.plot(r,r**2*xi0_m,'o',color=mscol[nd],markersize=4)
  ax0.plot(r,r**2*xi0_s,'o',color=sfcol[nd],markersize=4)
  ax0.plot(mgsm_ms_r,mgsm_ms_r**2*mgsm_ms_xi0,'--',color='gray')
  ax0.plot(mgsm_sf_r,mgsm_sf_r**2*mgsm_sf_xi0,'--',color='gray')
  ax0.plot(mdisp_ms_r,mdisp_ms_r**2*mdisp_ms_xi0,'-',color='gray')
  ax0.plot(mdisp_sf_r,mdisp_sf_r**2*mdisp_sf_xi0,'-',color='gray')

  ax0.text(60,120,r'$s^2\xi_0(s)$',fontsize=17)
  #ax0.set_xlabel(r'$s^2[{\rm Mpc}/h]$')
  ax0.set_xlim(xlim)
  ax0.set_ylim(xilim)

  ax0.text(xlim[0]-25,xilim[0]+20,r'${\rm (Mpc/}h)^2}$',fontsize=17,rotation=90) 

  pl.setp(ax0.get_xticklabels(),visible=False)
#  pl.setp(ax0.get_yticklabels(),visible=False)


  ax2.plot(r,-r**2*xi2_m,'o',color=mscol[nd],markersize=4)
  ax2.plot(r,-r**2*xi2_s,'o',color=sfcol[nd],markersize=4)
  ax2.plot(mgsm_ms_r,-mgsm_ms_r**2*mgsm_ms_xi2,'--',color='gray')
  ax2.plot(mgsm_sf_r,-mgsm_sf_r**2*mgsm_sf_xi2,'--',color='gray')
  ax2.plot(mdisp_ms_r,-mdisp_ms_r**2*mdisp_ms_xi2,'-',color='gray')
  ax2.plot(mdisp_sf_r,-mdisp_sf_r**2*mdisp_sf_xi2,'-',color='gray')

  ax2.text(60,120,r'$-s^2\xi_2(s)$',fontsize=17)
  #ax2.set_xlabel(r'$s^2[{\rm Mpc}/h]$')
  ax2.set_xlim(xlim)
  ax2.set_ylim(xilim)
  
  pl.setp(ax2.get_xticklabels(),visible=False)
  pl.setp(ax2.get_yticklabels(),visible=False)

  ldisp = 'Dispersion model' if nd == 0 else ''
  lgsm  = 'GSM' if nd == 0 else '' 

  ax4.plot(r,r**2*xi4_m,'o',color=mscol[nd],markersize=4)
  ax4.plot(r,r**2*xi4_s,'o',color=sfcol[nd],markersize=4)
  ax4.plot(mgsm_ms_r,mgsm_ms_r**2*mgsm_ms_xi4,'--',color='gray',label=lgsm)
  ax4.plot(mgsm_sf_r,mgsm_sf_r**2*mgsm_sf_xi4,'--',color='gray')
  ax4.plot(mdisp_ms_r,mdisp_ms_r**2*mdisp_ms_xi4,'-',color='gray',label=ldisp)
  ax4.plot(mdisp_sf_r,mdisp_sf_r**2*mdisp_sf_xi4,'-',color='gray')

  if nd == 0:
    ax4.legend(loc='upper left',fancybox=True,shadow=True,ncol=1,fontsize=8)

  ax4.text(60,50,r'$s^2\xi_4(s)$',fontsize=17)
#  ax4.set_xlabel(r'$s^2[{\rm Mpc}/h]$')
  ax4.set_xlim(xlim)
  ax4.set_ylim(xilim2)

  ax4.text(xlim[1]-20,xilim[0]-35,r'$s [{\rm Mpc/}h]$',fontsize=17) 
#  pl.setp(ax4.get_xticklabels(),visible=False)
 # pl.setp(ax4.get_yticklabels(),visible=False)
  
  axp2.plot(r,r**2*pxi2_m,'o',color=mscol[nd],markersize=4)
  axp2.plot(r,r**2*pxi2_s,'o',color=sfcol[nd],markersize=4)
  axp2.plot(mgsm_ms_r,mgsm_ms_r**2*mgsm_ms_pxi2,'--',color='gray')
  axp2.plot(mgsm_sf_r,mgsm_sf_r**2*mgsm_sf_pxi2,'--',color='gray')
  axp2.plot(mdisp_ms_r,mdisp_ms_r**2*mdisp_ms_pxi2,'-',color='gray')
  axp2.plot(mdisp_sf_r,mdisp_sf_r**2*mdisp_sf_pxi2,'-',color='gray')

  axp2.text(60,50,r'$s^2\xi_2^{\dagger}(s)$',fontsize=17)
#  axp2.set_xlabel(r'$s^2[{\rm Mpc}/h]$')
  axp2.set_xlim(xlim)
  axp2.set_ylim(xilim2)

#  pl.setp(ax2.get_xticklabels(),visible=False)
  pl.setp(axp2.get_yticklabels(),visible=False)

#  ax0.ylabel(r'$r^2\xi(r) [({\rm Mpc}/h)^2]$',fontsize=15)

  

pl.savefig('multipoles_comp_upd.pdf',bbox_inches='tight')




