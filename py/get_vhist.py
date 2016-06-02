from struct import *
from scipy.optimize import curve_fit
import numpy as np
import pylab as pl

from scipy.stats import norm

def gaussf(x, mu,sigma):
  return 1./(np.sqrt(2*np.pi)*sigma)* np.exp(-(x-mu)**2/(2.*sigma**2))

def get_vhist():
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
  zspace_type = ['None']
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
                      
                      f.close()
                      k += ngal

                    lmhalo = np.log10(mhalo) + 10.0

                    mhmin = 9.0
                    mhmax = 16.0
                    
                    nmass = 10
                    mharr = np.linspace(mhmin,mhmax,num=nmass)
                    mhbin = mharr[1] - mharr[0]

                    outf = "%s%s_sz_%s_q0%s_g0%s_ndens%s.sigma_M" % (datadir,gt,_sz,_q0,_g0,nd)
                    fout = file(outf,"w")
                    fout.write("#log(M_halo)\tsigma_v\tp10\tp90\n")
                    
                    for imh in range(nmass):
                      cc = np.where( (lmhalo > mharr[imh]-mhbin/2.) & 
                      (lmhalo < mharr[imh] + mhbin/2.))
                      
                      ncc = len(cc[0])
                      if len(cc[0]) == 0:
                        continue  
                      
                      nhtot = 0
                      sig_h = 0.0
                      halocc = lmhalo[cc]
                      dvcc    = dvh[cc]
                      
                      ihh = 0.0
                      h_x = halocc[ihh]
                      v_halo = []
                      sig_h = []
                      ngperh = 0
                       
                      for imm in range(ncc):
                        if halocc[imm] == h_x:
                          v_halo.append(dvcc[imm])
                          ngperh += 1
                        else:
                          if ngperh > 1:
                            #sig_h.append(np.std(v_halo))
                            nhtot += 1
                          
                          #v_halo = [dvh[imm]]
                          ngperh = 1
                          h_x = halocc[imm]

                      if nhtot > 1:
                        #sigma_m = np.median(sig_h)
                        #p10 = np.percentile(sig_h,16.0)
                        #p90 = np.percentile(sig_h,84.0)
                        sigma_m = np.std(v_halo)
                        medvhalo = np.mean(v_halo)
                        p10 = np.percentile(v_halo,10.0)
                        p90 = np.percentile(v_halo,90.0)

#                     pl.clf()
#                     n, bins, patches = pl.hist(dvh[cc], 100,normed=True)
#                      pl.setp(patches, 'facecolor','g','alpha',0.25)
                   
#                      hist, bin_edges = np.histogram(dvh[cc], bins=100)
#                      bin_centres = (bin_edges[:-1] + bin_edges[1:])/2.
                      
#                      p0 = [0.,300.]
#                      coeff, var_matrix = curve_fit(gaussf, bin_centres,hist,p0=p0)
#                      hist_fit = gaussf(bin_centres, coeff[0],coeff[1])
                      
#                      int_hist = np.sum(hist)
                      
  #                    n2, bins2, patches2 = pl.hist(dvh[cc], 100)
  #                    pl.setp(patches2, 'facecolor','r','alpha',0.75)
#                      pl.plot(bin_centres, hist_fit,'-',color='gray',linewidth=3)
  #                    pl.plot(bin_centres, hist,'--',color='blue',linewidth=3)

#                      mu, sigma = norm.fit(dvh[cc])
#                      sigma_m[imh] = np.std(dvh[cc])
                        fout.write("%f %f %f %f %f\n" % (mharr[imh],sigma_m,medvhalo,p10,p90))
##                      hist2 = gaussf(bin_centres,mu,sigma)  
#                      pl.plot(bin_centres,hist2 ,'-',color='red',linewidth=3)

                    fout.close()                  
#                    out.write("%f %f" % (sigma, sigma_all))
#                    fout.close()
                    
#                    print "sigma %f stddev %f" % (sigma, sigma_all)


#                    pl.draw()
#                    pl.show()

