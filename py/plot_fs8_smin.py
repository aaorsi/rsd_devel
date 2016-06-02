# Plot fs8 vs. smin as a function of a number of things: sigma_z, n_dens, galaxy type, RSD model, etc...


import sys
import pickle
import os.path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab as pl
sys.path.append('/home/CEFCA/aaorsi/')
import matplotlib as mpl
import matplotlib.gridspec as gsc
import matplotlib.colors as colors
import matplotlib.cm as cmx

mpl.rc('font', family = 'STIXGeneral')



datadir = '/home/CEFCA/aaorsi/cosmocodes/modelling/rsd/out/'

galtype     = ['Mstellar','Starforming']#,'OII_3727']
ngtype = len(galtype)

ndens_arr   = ['1e-3']#,'1e-4']
ndens_str   = [r'10^{-3}',r'5\times 10^{-4}',r'10^{-4}']
n_ndens = len(ndens_arr)

sigmaz_arr  = ['0.0']#,'0.001','0.003']
nsigma = len(sigmaz_arr)

rmin_arr = [10,20,30,40,50,60]
nrmin = len(rmin_arr)

fs8_arr = np.zeros([ngtype,n_ndens,nsigma,nrmin,3])

#f, axarr = pl.subplots(n_ndens,ngtype,sharex=True,sharey=True)
f, axarr = pl.subplots(1,1,sharex=True,sharey=True)
f.subplots_adjust(hspace=0,wspace=0)


# colormaps for separate galaxy types / colors for different sigmaz
sfcol = pl.cm.Blues(np.linspace(0.3, 1, nsigma)) 
mscol = pl.cm.Oranges(np.linspace(0.3, 1, nsigma)) 
o2col = pl.cm.Greens(np.linspace(0.3,1,nsigma))

#mscol = pl.cm.Spectral(np.linspace(0.2, 0.9, nsigma)) 
for ig in range(ngtype):
  gtype = galtype[ig]
  col = sfcol if gtype == 'Starforming' else mscol
  for ind in range(n_ndens):
    ndens = ndens_arr[ind]
    if ngtype > 1 and n_ndens >1:
      #ax = axarr[ind,ig]
      ax = axarr
    elif ngtype > 1 and n_ndens == 1:
#      ax = axarr[ig]
      ax = axarr
    elif n_ndens > 1 and ngtype == 1:
      ax = axarr[ind]
    else:
      ax = axarr
      

    for isz in range(nsigma):
      sz = sigmaz_arr[isz]

      for ir in range(nrmin):
        rmin = rmin_arr[ir]
        if sz == '0.0':
          fname_s12 = datadir + 'xilin.s12_fs8_'+gtype+'galdata_sz_'+sz+'_zspace_ndens'+ndens+'_bin0_rspacedm__rmin'+str(rmin) 
          fname = datadir + 'xilin.fs8_'+gtype+'galdata_sz_'+sz+'_zspace_ndens'+ndens+'_bin0_rspacedm__rmin'+str(rmin) 
        else:
          fname_s12 = datadir + 'xilin.s12_fs8_'+gtype+'galdata_sz_'+sz+'_zspace_typeNone_q02.8e7_g0-1.3sfog20.0_ndens'+ndens+'_bin0_rspacedm__rmin'+str(rmin) 
          fname = datadir + 'xilin.fs8_'+gtype+'galdata_sz_'+sz+'_zspace_typeNone_q02.8e7_g0-1.3sfog20.0_ndens'+ndens+'_bin0_rspacedm__rmin'+str(rmin) 
          
        if os.path.isfile(fname) is False:
          print '%s not found.' % (fname)
          continue

        f = open(fname,'r')
        raw = f.readlines()
        data = raw[1].split()
        f.close()

        fs8_arr[ig,ind,isz,ir,0] = float(data[0])
        fs8_arr[ig,ind,isz,ir,1] = float(data[1]) + float(data[0])
        fs8_arr[ig,ind,isz,ir,2] = float(data[0]) - float(data[2])
    
        if isz == 0 and ir == 0:
          data = pickle.load(open(fname_s12+'mcmcdata', 'rb'))
          fs8_true = data['mcmc_truths'][1]
          fs8_true_arr = np.zeros(nrmin)
          for j in range(nrmin):
            fs8_true_arr[j] = fs8_true

      zz  = '' if isz == 0 else '(1+z)'
      #lab = r'$\sigma_z = '+sz+zz+'$' if ind == n_ndens-1 else ''
      #lab = r'$\sigma_z = '+sz+zz+'$' if ind == n_ndens-1 else ''
      lab = gtype  if ind == n_ndens-1 else ''
      fsigma = np.abs(fs8_arr[ig,ind,isz,:,1] - fs8_arr[ig,ind,isz,:,2])
 #     ax.plot(rmin_arr,fs8_arr[ig,ind,isz,:,0],'o',color=col[isz],linewidth=3,label=lab)
      #ax.plot(rmin_arr,(fs8_arr[ig,ind,isz,:,0] - fs8_true_arr)/fsigma ,'o',color=col[isz],linewidth=3,label=lab,markersize=5)
      #ax.plot(rmin_arr,(fs8_arr[ig,ind,isz,:,0] - fs8_true_arr)/fs8_true_arr ,'o',color=col[isz],linewidth=3,label=lab,markersize=5)
      ax.plot(rmin_arr,(fs8_arr[ig,ind,isz,:,0] - fs8_true_arr)/fs8_true_arr ,'-o',color=col[ind],linewidth=3,label=lab,markersize=5)
      ax.fill_between(rmin_arr,(fs8_arr[ig,ind,isz,:,2] - fs8_true_arr)/fs8_true_arr,(fs8_arr[ig,ind,isz,:,1] - fs8_true_arr)/fs8_true_arr,color=col[isz],alpha=0.75)

      ax.set_ylim([-0.05,0.29])
      ax.set_xlim([rmin_arr[0]-5,rmin_arr[-1]+5])
     

      #ax.plot([-1,100],[fs8_true,fs8_true],'--',color='gray',linewidth=3)
      ax.plot([-1,100],[0,0],'--',color='gray',linewidth=3)
    if ind == n_ndens-1 and ig == ngtype-1 : 
      #ax.legend(loc='upper right',fancybox=True,shadow=True,fontsize=8)
      ax.legend(loc='upper right',fancybox=True,shadow=True,fontsize=15)
    """
    if ind == 0 and ig == 0:
      ax.text(-5,0.2,r'$f\sigma_8$',rotation=90,fontsize=20)
    """

#    if ind == 1 and ig == 0:
    if ind == 0 and ig == 0:
      ax.set_ylabel(r'$\Delta f\sigma_8 / f\sigma_8^{\rm true} $',rotation=90,fontsize=20)
      #ax.set_ylabel(r'$\Delta f\sigma_8 / \sigma_{f\sigma_8}$',rotation=90,fontsize=20)

  
    if ind == n_ndens-1:
      ax.set_xlabel(r'$s_{\rm min}[{\rm Mpc/}h]$',fontsize=20)


#    if ind == 0:
 #     ax.set_title(galtype[ig],fontsize=15)

    ax.text(.05,.95,r'$n = '+ndens_str[ind]+'$',transform=ax.transAxes, 
    #fontsize=10,verticalalignment='top')
    fontsize=15,verticalalignment='top')


pl.savefig('../plots/fs8_smin.pdf',bbox_inches='tight',dpi=300)

