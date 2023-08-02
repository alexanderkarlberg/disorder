!  The next subroutines, open some histograms and prepare them 
!      to receive data 
!  You can substitute these  with your favourite ones
!  init   :  opens the histograms
!  topout :  closes them
!  pwhgfill  :  fills the histograms with data

      subroutine init_histo
      implicit none
      include 'pwhg_bookhist-multi.h'
      integer xnbins, Q2nbins
      
      real * 8 Q2binsize, Q2min, Q2max, logQ2min, logQ2max
      real * 8 xbinsize, logxmin, logxmax, sbeams
      call inihists

            ! Bins in Q^2 and as a function of x
      xnbins = 50
      logxmin = -4d0
      logxmax = 3d0
      xbinsize = (logxmax - logxmin) / dble(xnbins)
      call bookupeqbins('ptj1_Q2_5.0',    xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_Q2_7.0',    xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_Q2_9.0',    xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_Q2_11.0',   xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_Q2_13.0',   xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_Q2_16.0',   xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_Q2_20.0',   xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_Q2_32.0',   xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_Q2_40.0',   xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_Q2_50.0',   xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_Q2_65.0',   xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_Q2_85.0',   xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_Q2_110.0',  xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_Q2_140.0',  xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_Q2_185.0',  xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_Q2_240.0',  xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_Q2_310.0',  xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_Q2_410.0',  xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_Q2_530.0',  xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_Q2_710.0',  xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_Q2_900.0',  xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_Q2_1300.0', xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_Q2_1800.0', xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_Q2_2500.0', xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_Q2_3500.0', xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_Q2_15000.0',xbinsize,logxmin,logxmax)

      call bookupeqbins('ptj1_x_0.00032',xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_x_0.0005', xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_x_0.0008', xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_x_0.0013', xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_x_0.002',  xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_x_0.0032', xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_x_0.005',  xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_x_0.008',  xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_x_0.013',  xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_x_0.02',   xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_x_0.032',  xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_x_0.05',   xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_x_0.08',   xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_x_0.13',   xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_x_0.2',    xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_x_0.32',   xbinsize,logxmin,logxmax)
      call bookupeqbins('ptj1_x_0.5',    xbinsize,logxmin,logxmax)
      
      xbinsize = 10d0 / dble(xnbins)

      call bookupeqbins('etaj1_Q2_5.0',    xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_Q2_7.0',    xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_Q2_9.0',    xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_Q2_11.0',   xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_Q2_13.0',   xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_Q2_16.0',   xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_Q2_20.0',   xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_Q2_32.0',   xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_Q2_40.0',   xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_Q2_50.0',   xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_Q2_65.0',   xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_Q2_85.0',   xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_Q2_110.0',  xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_Q2_140.0',  xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_Q2_185.0',  xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_Q2_240.0',  xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_Q2_310.0',  xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_Q2_410.0',  xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_Q2_530.0',  xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_Q2_710.0',  xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_Q2_900.0',  xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_Q2_1300.0', xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_Q2_1800.0', xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_Q2_2500.0', xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_Q2_3500.0', xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_Q2_15000.0',xbinsize,-5d0,5d0)

      call bookupeqbins('etaj1_x_0.00032',xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_x_0.0005', xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_x_0.0008', xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_x_0.0013', xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_x_0.002',  xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_x_0.0032', xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_x_0.005',  xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_x_0.008',  xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_x_0.013',  xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_x_0.02',   xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_x_0.032',  xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_x_0.05',   xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_x_0.08',   xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_x_0.13',   xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_x_0.2',    xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_x_0.32',   xbinsize,-5d0,5d0)
      call bookupeqbins('etaj1_x_0.5',    xbinsize,-5d0,5d0)
       
      end
      
      subroutine user_analysis(n,dsig,x,y,Q2)
      use parameters
      implicit none
      real * 8 dsig
      
      integer n
      double precision x, y, Q2, plab(0:3,n)
      integer maxjets
      parameter (maxjets=3)
      double precision ppartons(0:3,maxjets),pj(0:3,maxjets)
      double precision ptj1, etaj1
      real * 8 sqrtQ2,sbeams,R, palg, kt, eta, delx

      integer i,npartons,njets

      if(dsig.eq.0d0) return
      
      sqrtQ2 = sqrt(Q2)
      plab = 0d0

      ! Now fill the momenta
      if(n.eq.4) then           ! Born kinematics
        plab = pbornlab
      elseif(n.eq.5) then
        plab = preallab
      elseif(n.eq.6) then
        plab = prreallab
      else
        stop 'Wrong n in analysis'
      endif
      npartons = n - 3

      ppartons(:,1:npartons) = plab(:,4:n)

      R = 0.8d0
      palg = -1d0
      call buildjets(npartons,ppartons,R,palg,pj,njets)

      ptj1 = kt(pj(:,1))
      etaj1 = eta(pj(:,1))
      ptj1 = log(ptj1/sqrt(Q2))
     
      !     Now we fill histograms first in Q2 binning
      if(Q2.lt.4d0) then
         return
      elseif(Q2.lt.5d0) then
         call filld('ptj1_Q2_5.0',    ptj1, dsig / 1d0)
         call filld('etaj1_Q2_5.0',    etaj1, dsig / 1d0)
      elseif(Q2.lt.7d0) then
         call filld('ptj1_Q2_7.0',    ptj1, dsig / 2d0)
         call filld('etaj1_Q2_7.0',    etaj1, dsig / 2d0)
      elseif(Q2.lt.9d0) then
         call filld('ptj1_Q2_9.0',    ptj1, dsig / 2d0)
         call filld('etaj1_Q2_9.0',    etaj1, dsig / 2d0)
      elseif(Q2.lt.11d0) then
         call filld('ptj1_Q2_11.0',    ptj1, dsig / 2d0)
         call filld('etaj1_Q2_11.0',    etaj1, dsig / 2d0)
      elseif(Q2.lt.13d0) then
         call filld('ptj1_Q2_13.0',    ptj1, dsig / 2d0)
         call filld('etaj1_Q2_13.0',    etaj1, dsig / 2d0)
      elseif(Q2.lt.16d0) then
         call filld('ptj1_Q2_16.0',    ptj1, dsig / 3d0)
         call filld('etaj1_Q2_16.0',    etaj1, dsig / 3d0)
      elseif(Q2.lt.20d0) then
         call filld('ptj1_Q2_20.0',    ptj1, dsig / 4d0)
         call filld('etaj1_Q2_20.0',    etaj1, dsig / 4d0)
      elseif(Q2.lt.32d0) then
         call filld('ptj1_Q2_32.0',    ptj1, dsig / 12d0)
         call filld('etaj1_Q2_32.0',    etaj1, dsig / 12d0)
      elseif(Q2.lt.40d0) then
         call filld('ptj1_Q2_40.0',    ptj1, dsig / 8d0)
         call filld('etaj1_Q2_40.0',    etaj1, dsig / 8d0)
      elseif(Q2.lt.50d0) then
         call filld('ptj1_Q2_50.0',    ptj1, dsig / 10d0)
         call filld('etaj1_Q2_50.0',    etaj1, dsig / 10d0)
      elseif(Q2.lt.65d0) then
         call filld('ptj1_Q2_65.0',    ptj1, dsig / 15d0)
         call filld('etaj1_Q2_65.0',    etaj1, dsig / 15d0)
      elseif(Q2.lt.85d0) then
         call filld('ptj1_Q2_85.0',    ptj1, dsig / 20d0)
         call filld('etaj1_Q2_85.0',    etaj1, dsig / 20d0)
      elseif(Q2.lt.110d0) then
         call filld('ptj1_Q2_110.0',    ptj1, dsig / 25d0)
         call filld('etaj1_Q2_110.0',    etaj1, dsig / 25d0)
      elseif(Q2.lt.140d0) then
         call filld('ptj1_Q2_140.0',    ptj1, dsig / 30d0)
         call filld('etaj1_Q2_140.0',    etaj1, dsig / 30d0)
      elseif(Q2.lt.185d0) then
         call filld('ptj1_Q2_185.0',    ptj1, dsig / 45d0)
         call filld('etaj1_Q2_185.0',    etaj1, dsig / 45d0)
      elseif(Q2.lt.240d0) then
         call filld('ptj1_Q2_240.0',    ptj1, dsig / 55d0)
         call filld('etaj1_Q2_240.0',    etaj1, dsig / 55d0)
      elseif(Q2.lt.310d0) then
         call filld('ptj1_Q2_310.0',    ptj1, dsig / 70d0)
         call filld('etaj1_Q2_310.0',    etaj1, dsig / 70d0)
      elseif(Q2.lt.410d0) then
         call filld('ptj1_Q2_410.0',    ptj1, dsig / 100d0)
         call filld('etaj1_Q2_410.0',    etaj1, dsig / 100d0)
      elseif(Q2.lt.530d0) then
         call filld('ptj1_Q2_530.0',    ptj1, dsig / 120d0)
         call filld('etaj1_Q2_530.0',    etaj1, dsig / 120d0)
      elseif(Q2.lt.710d0) then
         call filld('ptj1_Q2_710.0',    ptj1, dsig / 180d0)
         call filld('etaj1_Q2_710.0',    etaj1, dsig / 180d0)
      elseif(Q2.lt.900d0) then
         call filld('ptj1_Q2_900.0',    ptj1, dsig / 190d0)
         call filld('etaj1_Q2_900.0',    etaj1, dsig / 190d0)
      elseif(Q2.lt.1300d0) then
         call filld('ptj1_Q2_1300.0',    ptj1, dsig / 400d0)
         call filld('etaj1_Q2_1300.0',    etaj1, dsig / 400d0)
      elseif(Q2.lt.1800d0) then
         call filld('ptj1_Q2_1800.0',    ptj1, dsig / 500d0)
         call filld('etaj1_Q2_1800.0',    etaj1, dsig / 500d0)
      elseif(Q2.lt.2500d0) then
         call filld('ptj1_Q2_2500.0',    ptj1, dsig / 700d0)
         call filld('etaj1_Q2_2500.0',    etaj1, dsig / 700d0)
      elseif(Q2.lt.3500d0) then
         call filld('ptj1_Q2_3500.0',    ptj1, dsig / 1000d0)
         call filld('etaj1_Q2_3500.0',    etaj1, dsig / 1000d0)
      elseif(Q2.lt.15000d0) then
         call filld('ptj1_Q2_15000.0',    ptj1, dsig / 11500d0)
         call filld('etaj1_Q2_15000.0',    etaj1, dsig / 11500d0)
      endif

      if(x.lt.0.0002d0) then
         return
      elseif(x.lt.0.00032d0) then
         delx = 0.00032d0 - 0.0002d0
         call filld('ptj1_x_0.00032',ptj1, dsig / delx)
         call filld('etaj1_x_0.00032',etaj1, dsig / delx)
      elseif(x.lt.0.0005d0) then
         delx = 0.0005d0 - 0.00032d0
         call filld('ptj1_x_0.0005',ptj1, dsig / delx)
         call filld('etaj1_x_0.0005',etaj1, dsig / delx)
      elseif(x.lt.0.0008d0) then
         delx = 0.0008d0 - 0.0005d0
         call filld('ptj1_x_0.0008',ptj1, dsig / delx)
         call filld('etaj1_x_0.0008',etaj1, dsig / delx)
      elseif(x.lt.0.0013d0) then
         delx = 0.0013d0 - 0.0008d0 
         call filld('ptj1_x_0.0013',ptj1, dsig / delx)
         call filld('etaj1_x_0.0013',etaj1, dsig / delx)
      elseif(x.lt.0.002d0) then
         delx = 0.0020d0 - 0.0013d0
         call filld('ptj1_x_0.002',ptj1, dsig / delx)
         call filld('etaj1_x_0.002',etaj1, dsig / delx)
      elseif(x.lt.0.0032d0) then
         delx = 0.0032d0 - 0.0020d0
         call filld('ptj1_x_0.0032',ptj1, dsig / delx)
         call filld('etaj1_x_0.0032',etaj1, dsig / delx)
      elseif(x.lt.0.005d0) then
         delx = 0.005d0 - 0.0032d0
         call filld('ptj1_x_0.005',ptj1, dsig / delx)
         call filld('etaj1_x_0.005',etaj1, dsig / delx)
      elseif(x.lt.0.008d0) then
         delx = 0.008d0 - 0.005d0
         call filld('ptj1_x_0.008',ptj1, dsig / delx)
         call filld('etaj1_x_0.008',etaj1, dsig / delx)
      elseif(x.lt.0.013d0) then
         delx = 0.013d0 - 0.008d0
         call filld('ptj1_x_0.013',ptj1, dsig / delx)
         call filld('etaj1_x_0.013',etaj1, dsig / delx)
      elseif(x.lt.0.02d0) then 
         delx = 0.02d0 - 0.013d0
        call filld('ptj1_x_0.02',ptj1, dsig / delx)
        call filld('etaj1_x_0.02',etaj1, dsig / delx)
      elseif(x.lt.0.032d0) then
         delx = 0.032d0 - 0.02d0
         call filld('ptj1_x_0.032',ptj1, dsig / delx)
         call filld('etaj1_x_0.032',etaj1, dsig / delx)
      elseif(x.lt.0.05d0) then
         delx = 0.05d0 - 0.032d0
         call filld('ptj1_x_0.05',ptj1, dsig / delx)
         call filld('etaj1_x_0.05',etaj1, dsig / delx)
      elseif(x.lt.0.08d0) then
         delx = 0.08d0 - 0.05d0
         call filld('ptj1_x_0.08',ptj1, dsig / delx)
         call filld('etaj1_x_0.08',etaj1, dsig / delx)
      elseif(x.lt.0.13d0) then
         delx = 0.13d0 - 0.08d0
         call filld('ptj1_x_0.13',ptj1, dsig / delx)
         call filld('etaj1_x_0.13',etaj1, dsig / delx)
      elseif(x.lt.0.2d0) then
         delx = 0.2d0 - 0.13d0
         call filld('ptj1_x_0.2',ptj1, dsig / delx)
         call filld('etaj1_x_0.2',etaj1, dsig / delx)
      elseif(x.lt.0.32d0) then
         delx = 0.32d0 - 0.2d0
         call filld('ptj1_x_0.32',ptj1, dsig / delx)
         call filld('etaj1_x_0.32',etaj1, dsig / delx)
      elseif(x.lt.0.5d0) then
         delx = 0.5d0 - 0.32d0
         call filld('ptj1_x_0.5',ptj1, dsig / delx)
         call filld('etaj1_x_0.5',etaj1, dsig / delx)
      endif
      end

      subroutine buildjets(n, pin, R, palg, pj, njets)
      implicit none
      integer n
      double precision pin(0:3,n), R, palg
      integer maxtrack,maxjet
      parameter (maxtrack=3,maxjet=3)

!     Output
      double precision pj(0:3,maxjet)
!     Internal
      integer mu, njets, ntracks, ijet, j, jetvec(maxtrack)
      double precision pjet(4,maxjet), ptmin
      double precision ptrack(4,maxtrack)

      ptrack = 0d0
      pjet = 0d0
      njets=0
      ntracks = n
      pj = 0d0
      ptmin = 0d0
      jetvec = 0

!     Fast jet needs 0 indexed tracks
      ptrack(4,1:n)=pin(0,1:n)
      do mu=1,3
         ptrack(mu,1:n)=pin(mu,1:n)
      enddo
      call fastjetppgenkt(ptrack,ntracks,r,palg,pjet,njets)

      if(njets.gt.3.or.njets.lt.1) then
         print*, 'njets out of bounds!!', njets
         stop
      endif

!     Back to 1:4 index
      do ijet=1,njets
         do mu=1,3
            pj(mu,ijet)=pjet(mu,ijet)
         enddo
         pj(0,ijet)=pjet(4,ijet)
      enddo
      end

      double precision function kt(p)
      implicit none
      double precision p(0:3)
      
      kt = sqrt(max(p(1)**2 + p(2)**2,0d0))
      end

      function eta(p)
      implicit none
      real * 8 eta, p(0:3), normp, norm

      normp = norm(p)
      if(normp.gt.p(3)) then
         eta = 0.5d0 * log((normp + p(3)) / (normp - p(3)))
      else
         eta = sign(1d100,p(3)) 
      endif
      end

      function norm(p)
      implicit none
      real * 8 norm, p(0:3)

      norm = p(1)**2 + p(2)**2 + p(3)**2
      norm = sqrt(max(0d0,norm))
      end

      function getrapidity0(p)
      implicit none
      real * 8 p(0:3),getrapidity0
      getrapidity0=0.5d0*log((p(0)+p(3))/(p(0)-p(3)))
      end

      subroutine getrapidity(p,y)
      implicit none
      real * 8 p(4),y
      y=0.5d0*log((p(4)+p(3))/(p(4)-p(3)))
      end

      function azi(p)
      implicit none
      real * 8 pi
      parameter(pi = 3.141592653589793D0)
      real * 8 azi,p(0:3)
      azi = atan(p(2)/p(1))
      if (p(1).lt.0d0) then
         if (azi.gt.0d0) then               
            azi = azi - pi
         else
            azi = azi + pi
         endif
      endif    
      end
