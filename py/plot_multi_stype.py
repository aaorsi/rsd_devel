# Plot real space correlation function and bias as a function of scale

import os.path
import sys
from struct import *
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

showlog = False   # Add log-scale for small scales

ndensarr = ['1e-4']
nstrings = [r'1\times 10^{-4}']#,r'10^{-4}']
ndens = len(ndensarr)
galtypearr = ['Starforming','Mstellar']
#galtypearr = ['Mstellar']
ntype = len(galtypearr)
#zspace_type = ['None', 'halo','halo_disp','halo_shuffle','halo_shuffle_sat','halo_disp_sat','halo_stddev']
#zspace_type = ['None', 'halo','halo_shuffle'],'halo_disp','halo_stddev']#,#'halo_shuffle_sat','halo_disp_sat','halo_stddev']
zspace_type = ['None', 'halo','halo_rvir']#,'halo_stddev']#,#'halo_shuffle_sat','halo_disp_sat','halo_stddev']
nztype = len(zspace_type)

#datadir_sfr = '/home/CEFCA/aaorsi/cosmocodes/data/datatest/'
datadir = '/home/CEFCA/aaorsi/data/'
catdir = '/home/CEFCA/aaorsi/data/GalCat/MXXL_extended/'
mcmap = pl.get_cmap('bone')
scmap = pl.get_cmap('cool')

cNorm = colors.Normalize(vmin = -1, vmax=ndens+1)

mscalmap = cmx.ScalarMappable(norm=cNorm,cmap=mcmap)
sscalmap = cmx.ScalarMappable(norm=cNorm,cmap=scmap)


nmulti = 2
mult_name = [r'$\xi_0$',r'$\xi_2$',r'$\xi_4$',r'$\xi_2^{\dagger}$']

#pl.figure('xireal')

nscales = 2 if showlog else 1

#gs1 = gsc.GridSpec(3*nmulti,nscales)
gs1 = gsc.GridSpec(nmulti,nscales)
gs1.update(left=0.05,right=0.48,wspace=0.0,hspace=0.0)

#gs2 = gsc.GridSpec(3*nmulti,nscales)
gs2 = gsc.GridSpec(nmulti,nscales)
gs2.update(left=0.5,right=0.95,wspace=0.0,hspace=0.0)

#f, axarr = pl.subplots(6,ntype*2)
#f.subplots_adjust(right=0.75)

xlim  = [1,49]

xilim = [[-10,99],
         [-30,99],
         [0,59],
         [0,79]]
         
blim  = [1,3.9] 

xlog = [1e-1,4.99]
dlim = [-24.9,24.9]

sfcol = plt.cm.Paired(np.linspace(0.3, 1, 7))
mscol = plt.cm.Paired(np.linspace(0.0, .75, 7))
o2col = plt.cm.Greens(np.linspace(0.3,1,7))

#gs.update(wspace=0,hspace=0,top=0.5)

lstyle = ['-','-','-','-','-','-','-']
nd = 0
for igt in range(ntype):
  gtype = galtypearr[igt]
  col = sfcol if gtype == 'Starforming' else mscol
  col = mscol
  if gtype == 'OII_3727':
    col = o2col

  grid = gs1 if igt == 0 else gs2

  fname = gtype + 'galdata_sz_'+sz+'_ndens'+ndensarr[0]+'_bin0_rspacedm_'
  for ell in range(nmulti):
    sig = -1 if ell == 1 else 1
     
    if showlog:
      ax1 = pl.subplot(grid[ell*3:ell*3 + 2,0])
      ax2 = pl.subplot(grid[ell*3:ell*3+2,1])
      ax3 = pl.subplot(grid[ell*3+2,0])
      ax4 = pl.subplot(grid[ell*3+2,1])
    else:
      #ax2 = pl.subplot(grid[ell*3:ell*3+2,0])
      #ax4 = pl.subplot(grid[ell*3+2,0])
      ax4 = pl.subplot(grid[ell,0])


    for izt in range(nztype):
      zt = zspace_type[izt]
      if zt is None:
        fname_zspace = gtype+'galdata_sz_'+sz+'_zspace_ndens'+ndensarr[nd]+'_bin0_rspacedm_'
      else:
        fname_zspace = gtype+'galdata_sz_'+sz+'_zspace_type'+zt+'_q02.8e7_g0-1.3sfog20.0_ndens'+ndensarr[nd]+'_bin0_rspacedm_'

      if zt is not 'None' and zt is not 'halo' and zt is not 'halo_rvir':
        fname_zspace = gtype+'galdata_sz_'+sz+'_zspace_typehalo_q02.8e7_g0-1.3sfog20.0_ndens'+ndensarr[nd]+'_shufflev_'+zt+'_bin0_rspacedm_'
        
      if zt is None and gtype == 'OII_3727':
        fname_zspace = gtype+'galdata_sz_'+sz+'_zspace_typeNone_q02.8e7_g0-1.3sfog20.0_ndens'+ndensarr[nd]+'_bin0_rspacedm_'
     
      if zt is 'halo_sigmam':
        fname_zspace = gtype+'galdata_sz_'+sz+'_zspace_typehalo_q02.8e7_g0-1.3sfog20.0_ndens'+ndensarr[nd]+'_shufflevMhalo_bin0_rspacedm_'
      if zt is 'halo_sigmamsat':
        fname_zspace = gtype+'galdata_sz_'+sz+'_zspace_typehalo_q02.8e7_g0-1.3sfog20.0_ndens'+ndensarr[nd]+'_shufflevMhalo_sats_bin0_rspacedm_'


      vcat =  "%sgaldata_sz_%s_zspace_typehalo_q02.8e7_g0-1.3sfog20.0_ndens%s" % (gtype,sz,ndensarr[nd])
      

      


      if zt is 'None':
        label = r'$s(v_z)$' 
        ngmax = 3e7 
        vz = np.zeros(ngmax)
        dvh = np.zeros(ngmax)
        mhalo = np.zeros(ngmax)
        k = 0
        for i in range(20):
          filev = "%s%s.%d.velocities" % (catdir,vcat, i)
          
          f = open(filev, 'rb')
          ngal = (unpack('i',f.read(4)))[0]
          #print 'Ngal %d' % (ngal)
          
          for _v in range(ngal):
            vz[k + _v] = (unpack('f',f.read(4)))[0]
          for _v in range(ngal):
            dvh[k+ _v] = (unpack('f',f.read(4)))[0]
          for _v in range(ngal):
            mhalo[k+ _v] = (unpack('f',f.read(4)))[0]
          
          f.close()
          k += ngal
          
        fsat = len(np.where(dvh != 0)[0])/(k + 0.0)
        print 'satellite fraction: ',fsat

      
      
      
      
      elif zt is 'halo':
        label = r'$s(v_h)$'
      elif zt is 'halo_disp':
        label = r'$s(v_h + \mathcal{G}(0,\sigma))$'
      elif zt is 'halo_disp_s12':
        label = r'$s(v_h + \mathcal{G}(0,\sigma_{12}))$'
      elif zt is 'halo_disp_sat':
        label = r'$s(v_h + \mathcal{G}(0,\sigma)_{sat}}))$'
      elif zt is 'halo_shuffle':
        label = r'$s(v_h + \mathcal{P}(v))$'
      elif zt is 'halo_shuffle_sat':
        label = r'$s(v_h + \mathcal{P}(v)_{sat})$'
      elif zt is 'halo_stddev':
        label = r'$s(v_h + \sigma_{std})$'
      elif zt is 'halo_sigmam':
        label = r'$s(v_h + \mathcal{P}(v,M_{halo}))$'
      elif zt is 'halo_sigmamsat':
        label = r'$s(v_h + \mathcal{P}(v,M_{halo})_{sats})$'
      else:
        label = zt


      if not os.path.isfile(datadir+fname_zspace+'xi.512.0.tree_DD'):
        print datadir + fname_zspace
        continue

      rr = np.logspace(-3,2.5,num=200)
      
      xidata = cio.load_correlation_functions(datadir,fname,outdata='xi',logbins=False) 
      xidata.update(cio.load_correlation_functions(datadir,fname_zspace,outdata='xi2d',logbins=False))
      xidata.update(xiutils.build_multipoles_sigpi(rr,xidata,r_int=xidata['r']))

      xileg = xidata['xileg']
      s = xidata['s']
      if showlog:      
        xidata_log = cio.load_correlation_functions(datadir,fname,outdata='xi',logbins=True) 
        xidata_log.update(cio.load_correlation_functions(datadir,fname_zspace,outdata='xi2d',logbins=True))
        xidata_log.update(xiutils.build_multipoles_sigpi(rr,xidata_log,r_int=xidata_log['r']))
        
        xileg_log = xidata_log['xileg']
        s_log = xidata_log['s']
     
      #import pdb ; pdb.set_trace()

      if zt is 'None':
        xizspace = xileg
        if showlog:
          xizspace_log = xileg_log

        if ell == 0:
          args = ['','1.0','0','300',sz,ndensarr[0],gtype,
               'Monopole','Quadrupole','Hexadecapole','PseudoQuadrupole']
          cosmo_pars,model_pars,data_pars = get_pars(args[1:])

          if showlog:
            datalog, x_arr, data, nkeys, keynames, nrarr = cov.get_subsample(cosmo_pars, data_pars,logbins=True)
            data_cov_log, data_corr_log = cov.build(datalog['data'],nkeys, datalog['nrarr'])
            data_cov_log *= 1./data_pars['covnbox']
            scales_log = {}
          else:
            x_arr, data, nkeys, keynames, nrarr = cov.get_subsample(cosmo_pars, data_pars)

          data_cov, data_corr = cov.build(data,nkeys, nrarr)
          data_cov *= 1./data_pars['covnbox']
          options = ['rspace', 'Monopole', 'Quadrupole', 
          'Hexadecapole', 'PseudoQuadrupole']
          scales = {}
          
          for key in options: 
            xscale, index = model.get_index(key, data_pars,xidata['r'])
            scales.update({key:xscale})
            if showlog:
              xscale_log, index_log = model.get_index(key, data_pars,xidata_log['r'])
              scales_log.update({key:xscale_log})
            

      xidata.update({'xscale':scales})  
      xidata.update(cov.get_errors(data_cov,xidata,keynames,nrarr))
      if showlog:
        xidata_log.update({'xscale':scales_log})  
        xidata_log.update(cov.get_errors(data_cov_log,xidata_log,keynames,datalog['nrarr']))

      errname = 'xi0_error' if ell == 0 else 'xi2_error'
      leg     = 'Monopole' if ell == 0 else 'Quadrupole'
      
      diag = np.sqrt(xidata[errname])
      xsc = xidata['xscale'][leg]
      if showlog:
        diag_log = np.sqrt(xidata_log[errname])
        xsc_log = xidata_log['xscale'][leg]
        ax1.semilogx(s_log,sig*s_log**2*xileg_log[ell,:],lstyle[izt],color=col[izt],linewidth=2,label=label,alpha=0.75)
#      if izt == 0:
#        ax1.fill_between(s,sig*s**2*(xileg[ell,:] - diag),sig*s**2*(xileg[ell,:] + diag),color=col[izt],alpha=0.5)
        ax1.set_xlim(xlog)
        ax1.set_ylim(xilim[ell])
        ax1.set_xticklabels([])
        ax1.spines['right'].set_visible(False)
        ax1.yaxis.set_ticks_position('left')

        ax1.fill_between(xlog,[xilim[ell][0],xilim[ell][0]],[xilim[ell][1],xilim[ell][1]],color='lightgray',alpha=.2)

      
      #ax2.plot(s,sig*s**2*xileg[ell,:],lstyle[izt],color=col[izt],linewidth=1,label=label,alpha=1)
#      ax2.fill_between(s,sig*s**2*(xileg[ell,:] - diag),sig*s**2*(xileg[ell,:] + diag),color=col[izt],alpha=0.5)
      #if showlog:
      #  ax2.set_xlim([5,119])
      #else:
      #  ax2.set_xlim(xlim)
      
      #ax2.set_ylim(xilim[ell])
      #if igt > 0:
      #  ax2.set_yticklabels([])
      
      #ax2.set_xticklabels([])
      
      if showlog:
        ax2.spines['left'].set_visible(False) 
        ax2.yaxis.set_ticks_position('right')
        ax3.semilogx(s_log,(xileg_log[ell,:] - xizspace_log[ell,:])/diag_log,
        color=col[izt],linewidth=2,alpha=1)
        ax3.set_xlim(xlog)
        ax3.set_ylim(dlim)
        if ell < nmulti-1:
          ax3.set_xticklabels([])
        if igt > 0:
          ax3.set_yticklabels([])
          ax1.set_yticklabels([])

        ax3.spines['right'].set_visible(False)
        ax3.yaxis.set_ticks_position('left')
        ax3.fill_between(xlog,[xilim[ell][0],xilim[ell][0]],[xilim[ell][1],xilim[ell][1]],color='lightgray',alpha=.2)
      
#      if zt is 'halo_sigmam':
      lwidth= 3 if zt is 'halo_sigmam' else 1.5
      
      fracerr = [0.01,0.05,0.1]
      fstyles = ':'
      ax4.plot(s,(xileg[ell,:] - xizspace[ell,:])/diag,color=col[izt],linewidth=lwidth,alpha=1,label=label)
      if izt == 0:
        for _ferr in range(len(fracerr)):
          ax4.plot(s,(fracerr[_ferr] * xizspace[ell,:])/diag,fstyles,color='black',linewidth=0.5)
          ax4.plot(s,(-fracerr[_ferr] * xizspace[ell,:])/diag,fstyles,color='black',linewidth=0.5)

      ax4.fill_between([s[0],s[-1]],[-1,-1],[1,1],color='lightgray',alpha=0.5)

      ax4.set_xlim(xlim)
      
      ax4.set_ylim(dlim)
      if igt > 0:
        ax4.set_yticklabels([])
      
      if showlog:
        ax4.spines['left'].set_visible(False) 
        ax4.yaxis.set_ticks_position('right')
        if ell < nmulti-1:
          ax4.set_xticklabels([])
      
#      if ell == 0:
 #       ax.plot(xidata['r'],xidata['r']**2*xidata['xi'],'-',color='gray',linewidth=0.5) 
      
#      ax.spines['left'].set_visible(False)
#      ax.yaxis.set_ticks_position('right')
      
      if ell == 0 and igt == 0:
        if showlog:
          ax1.legend(fontsize=7,loc='upper left',fancybox=True,shadow=True)
        else:
          #ax2.legend(fontsize=7,loc='upper right',fancybox=True,shadow=True)
          ax4.legend(fontsize=7,loc='upper right',fancybox=True,shadow=True)

      if ell == 0 and igt == 1 and izt == 0:
        #ax2.text(.7,.9,r'$n = '+nstrings[0]+'$',fontsize=10,transform=ax2.transAxes)
        ax4.text(.1,.9,r'$n = '+nstrings[0]+'$',fontsize=10,transform=ax4.transAxes)
        
      if zt is 'None' and ell == 0:
        ax4.text(.7,0.1,r'$f_{sat} = %.2f$' % (fsat),fontsize=10,transform=ax4.transAxes)
     
      #if ell == 0 and izt == 0:
      #  if showlog:
      #    ax1.text(.75,1.1,gtype,fontsize=12,transform=ax1.transAxes)
      #  else:
      #     #ax2.text(.5,1.1,gtype,fontsize=12,transform=ax2.transAxes)
      #     ax4.text(.5,1.1,gtype,fontsize=12,transform=ax4.transAxes)
           
      

#      if igt > 0:
#        ax.set_yticks([])
#      if ell != nmulti-1:
#        ax.set_xticks([])

      if igt == 0 and izt == 0:
#        ax1.text(.7,.7,xilim[ell][1]*0.75,mult_name[ell],fontsize=15,transform=ax1.transAxes)
      #if igt == 0 and ell == 0:
        if showlog:
          ax3.set_ylabel(r'$\Delta \xi_{\ell}/\sigma^{\rm diag}_{\ell}$',fontsize=15)
        else:
          subscr = '0' if ell == 0 else '2'
          ax4.set_ylabel(r'$\Delta \xi_{'+subscr+r'}/\sigma^{\rm diag}_{'+subscr+r'}$',fontsize=15)

      #
      #  if ell == 0:
      #    if showlog:
      #      ax1.set_ylabel(r'$s^2\xi_{0}(s)$',fontsize=15)
      #    else:
      #      ax2.set_ylabel(r'$s^2\xi_{0}(s)$',fontsize=15)
      #  if ell == 1: 
      #    if showlog:
      #      ax1.set_ylabel(r'$-s^2\xi_{2}(s)$',fontsize=15)
      #    else:
      #      ax2.set_ylabel(r'$-s^2\xi_{2}(s)$',fontsize=15)
      #"""

#      if igt == 1 and ell == nmulti-1:
      #  ax.text(xlim[0]*.3,xilim[ell][0]-20,r'$s[{\rm Mpc/}h]$',fontsize=15)
      if izt == 0 and ell == nmulti-1: 
        if showlog:
          ax4.text(-0.2,-0.5,r'$s[{\rm Mpc/}h]$',fontsize=15,transform=ax4.transAxes)
        else:  
          ax4.text(.35,-0.15,r'$s[{\rm Mpc/}h]$',fontsize=15,transform=ax4.transAxes)

#      if igt == 1 and ell == 0:
        #ax.text(xlim[1]*.1,xilim[ell][1]+10,r'$n = '+nstrings[0]+'$',fontsize=10)
        
#      ax.set_xticks([10,50,100])
#      ax.text(xlin_lim[0]-30,xilim[0]-4,r'$r[{\rm Mpc/}h]$',fontsize=15)
#      ax.text(xlin_lim[0]-30,xilim[1]+2,r'$n = '+nstrings[nd]+'$',fontsize=15)
     

#    if nd == 0:
#      pl.legend(prop={'size':10},loc='center left',bbox_to_anchor=(1.0,0.5),
#      ncol = 1,fancybox=True,shadow=True)

  
pl.savefig('../plots/short_xileg_zcomp'+ndensarr[0]+'.pdf',bbox_inches='tight')
#pl.show()
#pl.close("all")


