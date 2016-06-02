import sys
sys.path.append('/home/CEFCA/aaorsi/')
import matplotlib
#matplotlib.use('Agg')
import pylab as pl
import matplotlib.gridspec as gsc
import matplotlib.colors as colors
import matplotlib.cm as cmx

import cosmocodes.modelling.rsd.model as model
from cosmocodes.modelling.rsd.get_pars import *
import cosmocodes.io as cio
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.optimize import minimize

sz = '0.0'
pl.rc('font', family='STIXGeneral')

showlog = False   # Add log-scale for small scales

ell = 1   # Quadrupole ell = 1

func = model.vel_gauss
sfunc = 'gauss'

sminarr = [30,20,10,1]
nsmin = len(sminarr)

ndensarr = ['1e-3']
nstrings = [r'1\times 10^{-3}']#,r'10^{-4}']
ndens = len(ndensarr)
galtypearr = ['Mstellar','Starforming']
ntype = len(galtypearr)
#zspace_type = ['None', 'halo','halo_disp','halo_shuffle','halo_shuffle_sat','halo_disp_sat','halo_stddev']
zspace_type = ['None', 'halo']
nztype = len(zspace_type)

#datadir_sfr = '/home/CEFCA/aaorsi/cosmocodes/data/datatest/'
datadir = '/home/CEFCA/aaorsi/data/'

nmulti = 1
mult_name = [r'$\xi_2$',r'$\xi_4$',r'$\xi_2^{\dagger}$']

#pl.figure('xireal')

nscales = 2 if showlog else 1

gs1 = gsc.GridSpec(3*nmulti,nscales)
gs1.update(left=0.05,right=0.5,wspace=0.0,hspace=0.0)

gs2 = gsc.GridSpec(3*nmulti,nscales)
gs2.update(left=0.5,right=0.95,wspace=0.0,hspace=0.0)

mscol = plt.cm.Paired(np.linspace(0.0, .75, 5))
fitcol = plt.cm.gist_rainbow(np.linspace(0.2, .95, nsmin))


nd = 0
for igt in range(ntype):
  gtype = galtypearr[igt]

  grid = gs1 if igt == 0 else gs2
  fname = gtype + 'galdata_sz_'+sz+'_ndens'+ndensarr[0]+'_bin0_rspacedm_'
  #fname_zspace = gtype+'galdata_sz_'+sz+'_zspace_ndens'+ndensarr[nd]+'_bin0_rspacedm_'
  fname_zspace = gtype+'galdata_sz_'+sz+'_zspace_typeNone_q02.8e7_g0-1.3sfog20.0_ndens'+ndensarr[nd]+'_bin0_rspacedm_'
  fname_zspace_halo = gtype+'galdata_sz_'+sz+'_zspace_typehalo_q02.8e7_g0-1.3sfog20.0_ndens'+ndensarr[nd]+'_bin0_rspacedm_'

  rr = np.logspace(-3,2.5,num=200)
  
  xidata = cio.load_correlation_functions(datadir,fname_zspace,outdata='xi',logbins=False) 
  xidata.update(cio.load_correlation_functions(datadir,fname_zspace,outdata='xi2d',logbins=False))
  xidata.update(xiutils.build_multipoles_sigpi(rr,xidata,r_int=xidata['r']))

  xidata_halo = cio.load_correlation_functions(datadir,fname_zspace_halo,outdata='xi',logbins=False) 
  xidata_halo.update(cio.load_correlation_functions(datadir,fname_zspace_halo,outdata='xi2d',logbins=False))
  xidata_halo.update(xiutils.build_multipoles_sigpi(rr,xidata_halo,r_int=xidata_halo['r']))
  
  xileg = xidata['xileg']
  s = xidata['s']

  sv_arr = np.linspace(1.0,3000,10)
  ax1 = pl.subplot(grid[0:2,0])
  ax2 = pl.subplot(grid[2,0])

  ax1.plot(s,-s**2*xileg[1,:],color='blue',label='zspace') 
  ax1.plot(xidata_halo['s'],-xidata_halo['s']**2*xidata_halo['xileg'][1,:],color='red',label='no FoG')
  

  args = ['','1.0','0','300',sz,ndensarr[0],gtype,
       'Monopole','Quadrupole']#,'Hexadecapole','PseudoQuadrupole']
  cosmo_pars,model_pars,data_pars = get_pars(args[1:])
  x_arr, data, nkeys, keynames, nrarr = cov.get_subsample(cosmo_pars, data_pars)

  data_cov, data_corr = cov.build(data,nkeys, nrarr)
  data_cov *= 1./data_pars['covnbox']
  options = ['rspace', 'Monopole', 'Quadrupole', 
  'Hexadecapole', 'PseudoQuadrupole']
  scales = {}
 
  

  for key in options: 
    xscale, index = model.get_index(key, data_pars,xidata['r'])
    scales.update({key:xscale})

  xidata.update({'xscale':scales})  
  xidata.update(cov.get_errors(data_cov,xidata,keynames,nrarr))

  errname = 'xi0_error' if ell == 0 else 'xi2_error'
  leg     = 'Monopole' if ell == 0 else 'Quadrupole'

  diag = np.sqrt(xidata[errname])  
  
  def fitmodel(x, p):
    if p < 0:
      return np.zeros(x)
    conv = model.add_disp(xidata_halo,1.0,p,vfunc=func)
    conv.update(xiutils.build_multipoles_sigpi(rr,conv,r_int=xidata_halo['r']))

    return conv['xileg'][1,x]
  
  ax1.set_title(gtype)
  ax2.plot([-100,300],[0.,0],'--',color='black') 
  ax2.plot(xidata_halo['s'],np.abs((xidata_halo['xileg'][1,:] - xileg[1,:])/diag),linewidth=.25,color='red') 
  
  for ism in range(nsmin):
    
    cc = np.where(xidata_halo['s'] > sminarr[ism])
    xdata = xidata_halo['s'][cc[0]]
 
    print 'running curve_fit with smin = %d...' % (sminarr[ism])

    coeff, err = curve_fit(fitmodel,cc[0],xidata_halo['xileg'][1,cc[0]],p0=30.0,
    sigma=diag[cc[0]],bounds=(0.,100.),method='dogbox')
     
     #
    #@# res = minimize(fitmodel,cc[0],method='BFGS',bounds=(0,40.0)
     #
    
    xileg_disp = fitmodel(cc[0],coeff)
    ax1.plot(xdata,-xdata**2*xileg_disp,'-',linewidth=0.25,color=fitcol[ism],
    label=r'$s_{min} = %d Mpc/h, \sigma_v = %.1f km/s$' % (sminarr[ism], coeff)) 
    #ax1.text(0.6,0.9,r'$\sigma_v = %.1f {\rm km/s}$' % (coeff),transform=ax1.transAxes) 
    ax2.plot(xdata,np.abs((xileg_disp - xileg[1,cc[0]])/diag[cc[0]]),linewidth=0.25,color=fitcol[ism]) 

  ax2.set_ylim([-0.05,0.74])
  ax1.set_ylim([-19,120])
 
  ax1.set_xlim([0,150])
  ax2.set_xlim([0,150])
  if igt == 0:
    ax2.set_ylabel(r'$\Delta \xi_{\ell}/\sigma_{\ell}$')
  
  ax1.legend(loc='lower right',fontsize=6)
  
  if igt > 0:
    ax1.get_yaxis().set_ticks([])
    ax2.get_yaxis().set_ticks([])



pl.savefig('../plots/convolution_'+sfunc+'.pdf')

