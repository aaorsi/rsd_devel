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


NoData = False   # Plot the box, not the data
xirspacearr = [True,False]
#xirspace = True

ndensarr = ['1e-3']#,'5e-4','1e-4']
strndens = [r'$10^{-3}$',r'$5\times 10^{-4}$',r'$10^{-4}$']
ndens = len(ndensarr)
galtypearr = ['Mstellar','Starforming']
ngaltype = len(galtypearr)
datadir = '/home/CEFCA/aaorsi/cosmocodes/data/datatest/'

cNorm = colors.Normalize(vmin = -1, vmax=ndens+1)
sfcol = plt.cm.Blues(np.linspace(0.3, 1, 4))
mscol = plt.cm.Oranges(np.linspace(0.3, 1, 4))

f, axarr = pl.subplots(4,1)
f.subplots_adjust(wspace=0.0,hspace=0.0,right=0.75)

xlim  = [0,80]
ylim  = [-0.1,0.09]
lstyle = ['-','--']
#gs.update(right=0.75,hspace=0,wspace=0)

mult_name = [r'$\xi_0$',r'$\xi_2$',r'$\xi_4$',r'$\xi_2^{\dagger}$']
nmulti = 4

plotmodel = 'dispersion'


for nd in range(ndens):
  for iN in range(ngaltype):
    col = sfcol if iN == 1 else mscol
    fname = galtypearr[iN] + 'galdata_sz_'+sz+'_ndens'+ndensarr[nd]+'_bin0_rspacedm_'
    fname_zspace = galtypearr[iN]+'galdata_sz_'+sz+'_zspace_ndens'+ndensarr[nd]+'_bin0_rspacedm_'

# Read data first
    if NoData is False:

      args = ['','1.0','1','200',sz,ndensarr[nd],galtypearr[iN],
             'Monopole','Quadrupole','Hexadecapole','PseudoQuadrupole']
      cosmo_pars,model_pars,data_pars = get_pars(args[1:])
      model_pars['inputpkfile'] = '/home/CEFCA/aaorsi/cosmocodes/data/wmp1_matterpower.dat'

      xidata = cio.load_correlation_functions(datadir,fname,outdata='xi')
      xidata.update(cio.load_correlation_functions(datadir,fname_zspace,outdata='xi2d'))
      rr = np.logspace(-3,2.5,num=200)
      xidata.update(xiutils.build_multipoles_sigpi(rr,xidata,r_int=xidata['r']))

      mult_data = xidata['xileg']
      r = xidata['s']
# Get covariance matrices
      
      x_arr, data, nkeys, keynames, nrarr = cov.get_subsample(cosmo_pars, data_pars)
      data_cov, data_corr = cov.build(data,nkeys, nrarr)
      data_cov *= 1./data_pars['covnbox']
      xidata.update(cov.get_errors(data_cov,xidata,keynames,nrarr))
      
      options = ['rspace', 'Monopole', 'Quadrupole', 
      'Hexadecapole', 'PseudoQuadrupole']
      scales = {}
      for key in options: 
        xscale, index = model.get_index(key, data_pars,xidata['r'])
        scales.update({key:xscale})

      xidata.update({'xscale':scales})  

# Compute models

      _str = ['1.0', '40', '80',sz,ndensarr[nd],galtypearr[iN],'Monopole']
      _cosmo, _model, _data = get_pars(_str) 

      for xirspace in xirspacearr:
            
        xx = 'xirspace' if xirspace else 'xilin'
        legxx = r'$\xi_r^{\rm sim}$' if xirspace else r'$\xi_r^{\rm lin}$'
        linexx = '-' if xirspace else '--'


        sfile = cosmocode_dir + xx + '.sigmas_' + _data['covname']
        sdata = np.loadtxt(sfile,skiprows=1)
        
        nsigma = 3
        s12_var = np.zeros(nsigma)
        s2fog_var = np.zeros(nsigma)
        _cosmo['s12'] = sdata[0]

        s12_var[0] = sdata[0]
        s12_var[1] = sdata[0] + sdata[1]
        s12_var[2] = sdata[0] - sdata[2]
        
        s2fog_var[0] = sdata[3]
        s2fog_var[1] = sdata[3] + sdata[4]
        s2fog_var[2] = sdata[3] - sdata[5]

        _cosmo['s2fog'] = sdata[3]

  #      mdisp_int = interp1d(mdisp_r, mdisp_multi,kind='linear',bounds_error=False,fill_value=0.0)
  #      mgsm_int = interp1d(mgsm_r, mgsm_multi,kind='linear',bounds_error=False,fill_value=0.0)
         
        for im in range(nmulti):
          ax = axarr[im]#,nd]
          for isig in range(1):  
          #_im = im if im < 2 else im + 1
            ldisp = 'Dispersion model' if nd == 0 and iN == 0 else ''
            lgsm  = 'GSM' if nd == 0 and iN == 0 else '' 

            mdata = interp1d(r,r**2*mult_data[im,:],kind='linear',bounds_error=False,fill_value=0.0)

            _cosmo['s12'] = s12_var[isig]
            _cosmo['s2fog'] = s2fog_var[isig]
          

            modelf = model.set_model_function('xilinear',_model)
            ximodel = modelf(xidata['r'],_cosmo,data_pars=_data,use_rspace=xirspace)

            rm   = ximodel['Dispersion model']['s']
            mdisp_multi = ximodel['Dispersion model']['xileg']

            mgsm_r   = ximodel['GSM']['s']
            mgsm_multi = ximodel['GSM']['xileg']

            dm_disp = (mdisp_multi[im,:] - mdata(rm)/rm**2)/(mdata(rm)/rm**2)
            dm_gsm  = (mgsm_multi[im,:] - mdata(rm)/rm**2)/(mdata(rm)/rm**2)

            lwidth = 1 if isig == 0 else 0.5
            sval = _cosmo['s12'] if plotmodel == 'dispersion' else _cosmo['s2fog']
            sval = '%.1f' % (sval)
            short = '12' if plotmodel == 'dispersion' else 'FoG'
            lgt = galtypearr[iN]+r'$, \sigma_{'+short+'} = '+ sval + '$' if im == 0 else ''
            
            if plotmodel == 'dispersion' or plotmodel == 'both': 
              nnan = np.isfinite(dm_disp)
              ax.plot(rm[nnan],dm_disp[nnan] ,linexx,color=col[1],label=legxx,linewidth=lwidth)

            if plotmodel == 'GSM' or plotmodel == 'both':
              nnan = np.isfinite(dm_gsm)
              ax.plot(rm[nnan],dm_gsm[nnan],linexx,color=col[1],label=legxx,linewidth=lwidth)

            if im == 0:
              errname = 'xi0_error'
              leg     = 'Monopole'
            elif im == 1:
              errname = 'xi2_error'
              leg     = 'Quadrupole'
            elif im == 2:
              errname = 'xi4_error'
              leg     = 'Hexadecapole'
            elif im == 3:
              errname = 'pxi2_error'
              leg     = 'PseudoQuadrupole'
            
            xsc   = xidata['xscale'][leg]
            yerr  = np.sqrt(xidata[errname])
            dy1   = yerr/mdata(xsc) * xsc**2
            dy0   = -yerr/mdata(xsc) * xsc**2
          
            ax.fill_between(xsc,dy1,dy0,color=col[0],alpha=0.9)

            if im == 0 and iN == 0:
              ax.legend(loc='center left',bbox_to_anchor=(1,0.4),
              ncol=1,fancybox=True,shadow=True,fontsize=10)

      for im in range(nmulti):
        if im == 2:
          ylim2 = [-1,0.9]
        else:
          ylim2 = ylim
        ax = axarr[im]#,nd]
        ax.plot(xlim,[0,0],'-',color='gray')
#        ax.fill_between(xlim,[-0.02,-0.02],[0.02,0.02],facecolor='pink',alpha = 0.5)
#        ax.fill_between(xlim,[-0.01,-0.01],[0.01,0.01],facecolor='hotpink',alpha =0.5)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim2)
        #if nd == 2 and iN == 0:
        if nd == 0 and iN == 0:
          ax.text(xlim[1] - 10,ylim2[1]*.6,mult_name[im],fontsize=15)
        if im == 1 and iN == 0 and nd == 0:
          #ax.text(xlim[0]-30,ylim[0],r'$\frac{\xi_{\ell}^{\rm mod}-\xi_{\ell}^{\rm dat}}{\xi_{\ell}^{\rm dat}}$',
          ax.text(xlim[0]-30,ylim[0],r'$\frac{\xi_{\ell}^{\rm mod}-\xi_{\ell}^{\rm dat}}{\xi_{\ell}^{\rm dat}}$',
          fontsize=17,rotation=90) 
        if im < 3:
          pl.setp(ax.get_xticklabels(),visible=False)

        #if im ==2 and nd ==1:
        if im ==3 and nd ==0:
          ax.set_xlabel(r'$s[{\rm Mpc/}h]$',fontsize=17) 
 
        if nd > 0:
          pl.setp(ax.get_yticklabels(),visible=False)

        if im == 0:
          ax.set_title(plotmodel + r'$n=$'+strndens[nd] + r'$({\rm Mpc/}h)^{-3}$',fontsize=11)
        
#pl.show()
pl.savefig(plotmodel + '.xirspace.pdf',bbox_inches='tight')

