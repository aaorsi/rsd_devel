from struct import *
from scipy.optimize import curve_fit
import numpy as np
import pylab as pl

import sys
sys.path.append('/home/CEFCA/aaorsi/el-power/python/')

import readutils as read
from scipy.stats import norm

def gaussf(x, mu,sigma):
  return 1./(np.sqrt(2*np.pi)*sigma)* np.exp(-(x-mu)**2/(2.*sigma**2))

""" 
This function reads the velocity files associated
to different galaxy catalogues to construct
a histogram and measure (mu, sigma) fitting a
Gaussian.
"""

datadir = '/home/CEFCA/aaorsi/data/GalCat/MXXL_extended/'
Galtype = ['Starforming','Mstellar']#,'Halpha','OII_3727']
sz = ['0.0']
ndens = ['1e-4','5e-4','1e-3']
zspace = [1]#,0]
zspace_type = ['halo_rvir']
q0 = ['2.8e7']
g0 = ['-1.3']
sfog = ['20.0']
ngmax = 3e7
outdir = datadir

for gt in Galtype:
  for _sz in sz:
    for nd in ndens:
      for zsp in zspace:
        if zsp == 1:
          for ztype in zspace_type:
            for _q0 in q0:
              for _g0 in g0:
                for s2 in sfog:
                  
                  vroot = "%sgaldata_sz_%s_zspace_type%s_q0%s_g0%ssfog%s_ndens%s" % (gt, 
                  _sz,ztype,_q0,_g0,s2,nd)
                  print vroot

                  vz = np.zeros(ngmax)
                  dvh = np.zeros(ngmax)
                  mhalo = np.zeros(ngmax)
                  rdist = np.zeros(ngmax)

                  filepos = "%s%s" % (datadir,vroot)
                  pos = read.readgadget_pos(filepos)
                  
                  k = 0
                  for i in range(20):
                    filev = "%s%s.%d.velocities" % (datadir,vroot, i)

                    f = open(filev, 'rb')
                    ngal = (unpack('i',f.read(4)))[0]
                    print 'Ngal %d' % (ngal)
                    
                    for _v in range(ngal):
                      vz[k + _v] = (unpack('f',f.read(4)))[0]
                    for _v in range(ngal):
                      dvh[k+ _v] = (unpack('f',f.read(4)))[0]
                    for _v in range(ngal):
                      mhalo[k+ _v] = (unpack('f',f.read(4)))[0]
                    for _v in range(ngal):
                      rdist[k+ _v] = (unpack('f',f.read(4)))[0]
                    
                    f.close()
                    k += ngal

                  lmhalo = np.log10(mhalo) + 10.0
                  rvir_im = read.get_rvir(mhalo,1.0)

#                  mhmin = np.min(lmhalo[np.where(mhalo > 0)])
#                  mhmax = np.max(lmhalo)

                  mhmin = 12.5
                  mhmax = 15.5

                  nmass = 5
                  mharr = np.linspace(mhmin,mhmax,num=nmass)
                  mhbin = mharr[1] - mharr[0]

                  outf = "%s%s_sz_%s_q0%s_g0%s_ndens%s.sigma_RM" % (datadir,gt,_sz,_q0,_g0,nd)
                  fout = file(outf,"w")
                  fout.write("#log(M_halo)\tradius\tsigmav\tnhaloes\n")
                  
                  rmin = 0.
                  rmax = 1.2
                  nrbins = 10 
                  rarr = np.linspace(rmin,rmax,num=nrbins)
                  rbin = rarr[1] - rarr[0]

                  for imh in range(nmass):
                    
                    cc = np.where( (lmhalo > mharr[imh]-mhbin/2.) & 
                    (lmhalo < mharr[imh] + mhbin/2.) & (dvh != 0))
                   

                    ncc = len(cc[0])
                    if len(cc[0]) == 0:
                      continue  
                    
                    for irr in range(nrbins):
                      crad = np.where((rdist[cc]/rvir_im[cc] > rarr[irr] - rbin/2.) & (
                      rdist[cc]/rvir_im[cc] < rarr[irr] + rbin/2.))
                    
                      nrad = len(crad[0])

                      if nrad == 0:
                        fout.write("%f %f %f %f\n" % (mharr[imh],rarr[irr],0.0,0.0))
                        continue

                      sigma_m = np.std(dvh[cc[0][crad[0]]])
                      nhtot = nrad
                      """
                      nhtot = 0
                      sig_h = 0.0
                      halocc = lmhalo[crad]
                      dvcc    = dvh[crad]
                      
                      ihh = 0.0
                      h_x = halocc[ihh]
                      v_halo = []
                      sig_h = []
                      ngperh = 0
                      sigma_m = 0.0   
                      



                      for imm in range(nrad):
                        if halocc[imm] == h_x:
                          v_halo.append(dvcc[imm])
                          ngperh += 1
                        else:
                          nhtot += 1
                          
                          ngperh = 1
                          h_x = halocc[imm]

                      if nhtot > 1:
                        sigma_m = np.std(v_halo)
                        medvhalo = np.mean(v_halo)
                        p10 = np.percentile(v_halo,10.0)
                        p90 = np.percentile(v_halo,90.0)
                      """

                      fout.write("%f %f %f %d\n" % (mharr[imh],rarr[irr],sigma_m,nhtot))

                  fout.close()                  

