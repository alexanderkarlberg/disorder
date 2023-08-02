c -*- Fortran -*-
      integer nhist,maxbins,maxmulti
      parameter (nhist=400,maxbins=500,maxmulti=10)
      character * 50 stringhist(nhist)
      real * 8 xhistarr(maxbins+1,nhist),
     1     yhistarr(maxmulti,0:maxbins+1,nhist)
      integer nhits(0:maxbins+1,nhist),nbins(nhist),jhist,
     1     ient1(nhist),nmulti
      real * 8 yhistarr1(maxmulti,0:maxbins+1,nhist),
     1       errhistarr1(maxmulti,0:maxbins+1,nhist),
     1       yhistarr2(maxmulti,0:maxbins+1,nhist),
     3       errhistarr2(maxmulti,0:maxbins+1,nhist)
      common/histnew/xhistarr,yhistarr,
     1       yhistarr1,errhistarr1,
     3       yhistarr2,errhistarr2,
     5       nhits,nbins,jhist,ient1,nmulti,
     6       stringhist
      save /histnew/

