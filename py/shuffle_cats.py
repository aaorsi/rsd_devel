from os.path import expanduser
from struct import *
from scipy.optimize import curve_fit
import numpy as np
import pylab as pl

import sys
sys.path.append('/home/CEFCA/aaorsi/el-power/python/')

import readutils as read
from scipy.stats import norm

from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0 = 73.0, Om0 = 0.25)
redshift = 1.0

aexp = 1./(1. + redshift)

home = expanduser("~")

def gaussf(x, mu,sigma):
  return 1./(np.sqrt(2*np.pi)*sigma)* np.exp(-(x-mu)**2/(2.*sigma**2))


def get_shuffle():
  """ 
  This function reads the velocity files associated
  to different galaxy catalogues to construct
  a histogram and measure (mu, sigma) fitting a
  Gaussian.
  """

  if home == '/Users/aaorsi': #That means laptop copy
    datadir = '/Users/aaorsi/work/rsd_devel/data/cat/'
  else:
    datadir = '/home/CEFCA/aaorsi/data/GalCat/MXXL_extended/'

  Galtype = ['Starforming','Mstellar']#,'Halpha','OII_3727']
  sz = ['0.0']
  ndens = ['1e-3']
  zspaceall = [1]#,0]
  zspace_type = ['halo_rvir']
  q0 = ['2.8e7']
  g0 = ['-1.3']
  sfog = ['20.0']
  ngmax = 3e7
  outdir = datadir

  write_params = True

  UseRvir      = True     # Only select satellites that are within R_vir.
  splitsats    = False
  Shuffletype  = 'mhalo'  # All, mhalo, radius


  for gt in Galtype:
    for _sz in sz:
      for nd in ndens:
        for zsp in zspaceall:
          if zsp == 1:
            for ztype in zspace_type:
              for _q0 in q0:
                for _g0 in g0:
                  for s2 in sfog:
                   
                    nfiles = 20
                    sats = '_sats' if splitsats else ''
                    rvirstring = '_rvir' if UseRvir else ''
                    vroot = "%sgaldata_sz_%s_zspace_type%s_q0%s_g0%ssfog%s_ndens%s" % (gt,
                    _sz,ztype,_q0,_g0,s2,nd)
                    outfileroot  = "%s_shufflev%s_%s%s" % (vroot,sats,Shuffletype,rvirstring)
                    outf = "%s%s" % (datadir, outfileroot)

                    read.write_paramfile(numfiles = nfiles,
                    outpowerdir='/home/CEFCA/aaorsi/data/',
                    outfile=outfileroot, outdir = datadir,
                    paramdir = '/home/CEFCA/aaorsi/el-power/params/',ranfile=outfileroot,randir = datadir)
                    print vroot

                    vz = np.zeros(ngmax)
                    dvh = np.zeros(ngmax)
                    mhalo = np.zeros(ngmax)
                    rdist = np.zeros(ngmax)

                    k = 0
                    filepos = "%s%s" % (datadir,vroot)
                    pos = read.readgadget_pos(filepos)

                    for i in range(20):
                      filev = "%s%s.%d.velocities" % (datadir,vroot, i)
                      
                      f = open(filev, 'rb')
                      ngal = (unpack('i',f.read(4)))[0]
                      #print 'Ngal %d' % (ngal)
                      
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

                    zspace = np.zeros(k)
                    xypos = np.zeros([k,2])
                    rvir  = read.get_rvir(mhalo,1.0)

                    # conditions to be applied
  
                    czero = mhalo != 0.0
                    ccsat = dvh != 0.0
                    crvir = rdist < rvir

                    if Shuffletype == 'mhalo': 
                      lmhalo = np.log10(mhalo) + 10.0

                      nmass = 10
                      nzz   = np.where(mhalo > 0)
                      mhmin = lmhalo[nzz].min() 
                      mhmax = lmhalo[nzz].max()
                      
                      mhbin = (mhmax - mhmin)/(nmass+0.0)
                      mharr = np.linspace(mhmin + mhbin/2.,mhmax - mhbin/2.,num=nmass)

  #                    outf = "%s%s_sz_%s_q0%s_g0%s_ndens%s.sigma_M" % (datadir,gt,_sz,_q0,_g0,nd)
  #                    fout = file(outf,"w")
  #                    fout.write("#log(M_halo)\tsigma_v\tp10\tp90\n")
                     
                      satc = ccsat if splitsats else 1
                      rvirc = crvir if UseRVir else 1
                      
                      igg = 0
                      
                      for imh in range(nmass):
                        
                        cc = (lmhalo > mharr[imh]-mhbin/2.) & 
                        (lmhalo <= mharr[imh] + mhbin/2.) 

                        
                        allcond = np.where(czero & cc & satc & rvirc)
                        nall = len(allcond[0])
                        if splitsats:
                          
                          centcond = np.where(czero & cc & ~satc & rvirc)
                          ncent = len(centcond[0])
                          if ncent > 0:
                            for ii in range(ncent):
                              ccid = centcond[0][ii]
                              zspace[igg] = pos[ccid,2]
                              xypos[igg,:] = pos[ccid,0:2]
                              igg += 1
                         
                        if nall > 0:
                          shuf = np.random.permutation(allcond[0])
                          print 'introducing shuffled velocities' 
                          for ii in range(ncc):
                            ccid = allcond[0][ii]
                            zspace[igg] = pos[ccid,2] - dvh[shuf[ii]]/(aexp * cosmo.H(redshift).value/cosmo.h)
                            xypos[igg,:] = pos[ccid,0:2]
                            if zspace[igg] < 0:
                              zspace[igg] += 3000.0
                            elif zspace[igg] > 3000.0:
                              zspace[igg] -= 3000.0
                            igg += 1
                        
                          
                    elif Shuffletype == 'All':
                      igg = 0
                     
                      cc = np.where(czero & satc & rvirc)
                      shuf = np.random.permutation(cc[0])
                      ncc = len(cc[0])
                      if splitsats:
                          centcond = np.where(czero & ~satc & rvirc)
                          ncent = len(centcond[0])
                          if ncent > 0:
                            for ii in range(ncent):
                              ccid = centcond[0][ii]
                              zspace[igg] = pos[ccid,2]
                              xypos[igg,:] = pos[ccid,0:2]
                              igg += 1
                      
                      if ncc > 0:
                        for ii in range(ncc):
                          ccid = cc[ii]
                          zspace[igg] = pos[ccid,2] - dvh[shuf[ii]]/(aexp * 
                          cosmo.H(redshift).value/cosmo.h)
                          xypos[igg,:] = pos[ccid,0:2]

                          if zspace[igg] < 0:
                            zspace[igg] += 3000.0
                          elif zspace[igg] > 3000.0:
                            zspace[igg] -= 3000.0
                          igg += 1

                    nperfile = round(k / (nfiles+0.0))
                    for i in range(nfiles):
                      file_i = "%s.%d" % (outf, i)
                      id0 = i * nperfile
                      id1 = (i + 1) * nperfile if i < (nfiles -1) else k
                      
                      ngalfile = id1 - id0
                      
                      GalData = np.empty(ngalfile,dtype=("(3,)f4"))
                      GalData[:,0:2] = xypos[id0:id1,:]
                      GalData[:,2]   = zspace[id0:id1]
#                      print 'writing %s \n' % (file_i)
                      read.make_gadget_output(file_i,ngalfile,GalData,nfiles=nfiles)


                                          

