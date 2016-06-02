# Plot multipoles


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

cosmocode_dir = '/home/CEFCA/aaorsi/cosmocodes/modelling/rsd/'
sz = '0.0'
pl.rc('font', family='STIXGeneral')


LogBins = False
NoData = False   # Plot the box, not the data
xirspacearr = [True]
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

volcol = plt.cm.Greys(np.linspace(0.1,0.4,3))

mult_name = [r'$\xi_0$',r'$\xi_2$',r'$\xi_4$',r'$\xi_2^{\dagger}$']
nmulti = 2

gs1 = gsc.GridSpec(3,1)
gs1.update(left=0.05,right=0.48,wspace=0.0,hspace=0.0)

gs2 = gsc.GridSpec(3,1)
gs2.update(left=0.52,right=0.95,wspace=0.0,hspace=0.0)


#f, axarr = pl.subplots(nmulti*2,2)
#f.subplots_adjust(wspace=0.0,hspace=0.0,right=0.75)

xlim  = [1,120]
ylim  = [-5.,4.9]
lstyle = ['-','--']
#gs.update(right=0.75,hspace=0,wspace=0)

ax_x0 = pl.subplot(gs1[0:2,0])
ax_ex0 = pl.subplot(gs1[2,0])

ax_x2 = pl.subplot(gs2[0:2,0])
ax_ex2 = pl.subplot(gs2[2,0])

for nd in range(ndens):
  for iN in range(ngaltype):
    for xirspace in xirspacearr:

      xx = 'xirspace' if xirspace else ''
      col = sfcol if iN == 1 else mscol
      fname = galtypearr[iN] + 'galdata_sz_'+sz+'_ndens'+ndensarr[nd]+'_bin0_rspacedm_'
      fname_zspace = galtypearr[iN]+'galdata_sz_'+sz+'_zspace_typehalo_q02.8e7_g0-1.3sfog20.0_ndens'+ndensarr[nd]+'_bin0_rspacedm_'

# Read data first
      if NoData is False:

        args = ['','1.0',str(xlim[0]),str(xlim[1]),'0.0',ndensarr[nd],galtypearr[iN],
               'Monopole','Quadrupole']
        cosmo_pars,model_pars,data_pars = get_pars(args[1:])
        model_pars['inputpkfile'] = '/home/CEFCA/aaorsi/cosmocodes/data/wmp1_matterpower.dat'

        xidata = cio.load_correlation_functions(datadir,fname,outdata='xi',logbins=LogBins)
        xidata.update(cio.load_correlation_functions(datadir,fname_zspace,outdata='xi2d'),logbins=LogBins)
        rr = np.logspace(-3,2.5,num=200)
        xidata.update(xiutils.build_multipoles_sigpi(rr,xidata,r_int=xidata['r']))

        mult_data = xidata['xileg']
        r = xidata['s']
# Get covariance matrices
        
        x_arr, data, nkeys, keynames, nrarr = cov.get_subsample(cosmo_pars, data_pars,logbins=LogBins)
        data_cov, data_corr = cov.build(data,nkeys, nrarr)
        data_cov *= 1./data_pars['covnbox']
        xidata.update(cov.get_errors(data_cov,xidata,keynames,nrarr))
        
        options = ['rspace', 'Monopole', 'Quadrupole', 'Hexadecapole', 'PseudoQuadrupole']
        scales = {}
        for key in options: 
          xscale, index = model.get_index(key, data_pars,xidata['r'])
          scales.update({key:xscale})

        xidata.update({'xscale':scales})  

    # Compute models

        _str = ['1.0', str(xlim[0]), str(xlim[1]),sz,ndensarr[nd],galtypearr[iN],'Monopole']
        _cosmo, _model, _data = get_pars(_str) 

        s12file = cosmocode_dir + '/out/xilin.s12_'+galtypearr[iN]+'galdata_sz_0.0_zspace_typeNone_q02.8e7_g0-1.3sfog20.0_ndens'+ndensarr[nd]+'_bin0_rspacedm_' 
        sfogfile = cosmocode_dir + '/out/xilin.s2fog_'+galtypearr[iN]+'galdata_sz_0.0_zspace_ndens'+ndensarr[nd]
        s12data = np.loadtxt(s12file,skiprows=1)
        _cosmo['s12'] = s12data[0]
        sfogdata = np.loadtxt(sfogfile,skiprows=1)
        _cosmo['s2fog'] = sfogdata[0]


        _cosmo['s12'] = 0.0
        _cosmo['s2fog'] = 0.0

        modelf = model.set_model_function('xilinear',_model)
        
        ximodel = modelf(xidata['r'],_cosmo,data_pars=_data,use_rspace=xirspace)

        rm   = ximodel['Dispersion model']['s']
        mdisp_multi = ximodel['Dispersion model']['xileg']

        mgsm_r   = ximodel['GSM']['s']
        mgsm_multi = ximodel['GSM']['xileg']

  #      mdisp_int = interp1d(mdisp_r, mdisp_multi,kind='linear',bounds_error=False,fill_value=0.0)
  #      mgsm_int = interp1d(mgsm_r, mgsm_multi,kind='linear',bounds_error=False,fill_value=0.0)
        

        for im in range(nmulti):
         
          ylim = [0,100] if im == 0 else [-10,100]
          xlim = [0,99]
          sign = 1 if im == 0 else -1
          axl = ax_x0 if im == 0 else ax_x2 
          axs = ax_ex0 if im == 0 else ax_ex2

          #_im = im if im < 2 else im + 1
          ldisp = 'Dispersion model' if nd == 0 and iN == 0 else ''
          lgsm  = 'GSM' if nd == 0 and iN == 0 else '' 
          _im = im 
          #ax1 = axarr[im]#,nd]
          mdata = interp1d(r,r**2*mult_data[_im,:],kind='linear',
          bounds_error=False,fill_value=0.0)

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
          
          sdiag = interp1d(xsc,yerr,fill_value=0.0,bounds_error=False)

          dm_disp = (mdisp_multi[_im,:] - mdata(rm)/rm**2)/sdiag(rm)
          dm_gsm  = (mgsm_multi[_im,:] - mdata(rm)/rm**2)/sdiag(rm)
          
          

          if xirspace:
            axs.plot(rm,dm_disp ,'-',color='gray',linewidth=3,label=ldisp)
            axs.plot(rm,dm_gsm,'--',color='gray',linewidth=3,label=lgsm)

          lgt = galtypearr[iN] if im == 1 else ''
       
          axl.plot(r,sign*r**2*mult_data[_im,:],'o',color=col[1],label=lgt,markersize=3)
          axl.plot(rm,sign*rm**2*mdisp_multi[_im,:],'-',color='gray',label=ldisp)
          axl.plot(rm,sign*rm**2*mgsm_multi[_im,:],'--',color='gray',label=lgsm)
          
          axl.fill_between(r,sign*r**2*(mult_data[_im,:] - sdiag(r)),sign*r**2*(mult_data[_im,:] + sdiag(r)),color=col[1],
          alpha=0.5)

          axs.plot(rm,dm_disp ,'-',color=col[1],linewidth=2)
          axs.plot(rm,dm_gsm,'--',color=col[1],linewidth=2)
  
          if im == 0: 
            axl.set_ylabel(r'$\i^{\ell}s^2\xi_{\ell}[({\rm Mpc/}h)^2]$',fontsize=12)
            axs.set_ylabel(r'$\Delta \xi_{\ell}/\sigma_{\ell}^{\rm diag}$',fontsize=12)

          axl.set_xticklabels([])

          axl.set_ylim(ylim)
          axl.set_xlim(xlim)
          axs.set_xlim(xlim)
          axs.set_ylim([-5,4.9]) 
          #ax.fill_between(xsc,dy1,dy0,color=col[0],alpha=0.9)
          
          if im == 1:
            axl.legend(loc='lower right',ncol=1,fancybox=True,shadow=True,fontsize=10)

          axs.plot(xlim,[0,0],'--',color='gray')
          vols = [27.0]
          nvol = len(vols)

          for ii in range(nvol):
            siglim = np.sqrt(27./vols[ii])  
            axs.fill_between(xlim,[-siglim,-siglim],[siglim,siglim],facecolor=volcol[-(ii+1)],
            alpha = 0.7)
#            if im == 1 and iN == 0: 
#              if vols[ii] == 27.0:
#                axs.text(60,siglim-0.6,r"$%d (Gpc/h)^3 (MXXL)$" % (vols[ii]),fontsize=6) 
#              else:
#                axs.text(60,siglim-0.6,r"$%d (Gpc/h)^3$" % (vols[ii]),fontsize=6) 

  #        ax.fill_between(xlim,[-0.01,-0.01],[0.01,0.01],facecolor='hotpink',alpha =0.5)
          #if nd == 2 and iN == 0:
  #        if im ==1:
  #          axl.text(0.75,.6,mult_name[im],fontsize=15,transform=ax.transAxes)
  #        if im == 1 and iN == 0 and nd == 0:
            #ax.text(xlim[0]-30,ylim[0],r'$\frac{\xi_{\ell}^{\rm mod}-\xi_{\ell}^{\rm dat}}{\xi_{\ell}^{\rm dat}}$',
  #          ax.text(-0.25,0.0,r'$\frac{\xi_{\ell}^{\rm mod}-\xi_{\ell}^{\rm dat}}{\sigma_{\ell}}$',
  #          fontsize=17,rotation=90,transform=ax.transAxes) 
  #        if im < 3:
  #          pl.setp(ax.get_xticklabels(),visible=False)

          #if im ==2 and nd ==1:
          axs.set_xlabel(r'$s[{\rm Mpc/}h]$') 
   
  #        if nd > 0:
  #          pl.setp(ax.get_yticklabels(),visible=False)

          
          axl.text(0.5,0.85,r'$n=$'+strndens[nd] + r'$({\rm Mpc/}h)^{-3}$',fontsize=12,transform=axl.transAxes)
        
          
#pl.show()
pl.savefig('dmulti_halo_smallscales.png',bbox_inches='tight',dpi=300)




