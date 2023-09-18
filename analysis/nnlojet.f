!  The next subroutines, open some histograms and prepare them 
!      to receive data 
!  You can substitute these  with your favourite ones
!  init   :  opens the histograms
!  topout :  closes them
!  pwhgfill  :  fills the histograms with data

      subroutine define_histograms
      implicit none
      include 'pwhg_bookhist-multi.h'
      integer xnbins, Q2nbins
      double precision, parameter :: etbins(9) = (/6d0, 8d0, 10d0, 14d0,
     $     21d0, 29d0, 47d0, 71d0, 127d0 /)
      
      real * 8 Qmax

      Qmax = sqrt(1000d0)

!      call bookupeqbins('Etj1',2d0,6d0,106d0)
      call bookup('Etj1',8,etbins)
      call bookupeqbins('etaj1',0.5d0,-1d0,3d0)
       
      end
      
      subroutine user_analysis(n,dsig,x,y,Q2)
      use mod_parameters
      use mod_analysis
      implicit none
      real * 8 dsig(maxscales)
      
      integer n
      double precision x, y, Q2
      integer maxjets
      parameter (maxjets=3)
      double precision ppartons(0:3,maxjets),pj(0:3,maxjets)
      double precision ptj1, etaj1, Etj1, Elep, Et
      real * 8 sqrtQ2,sbeams,R, palg, kt, eta

      integer i,npartons,njets

      if (y < 0.04d0 .or. Q2 < 25d0 )
     $     return
      
      Elep = plab(0,3)

      if(Elep .lt. 10d0) return
      
      sqrtQ2 = sqrt(Q2)

      npartons = n - 3

      ppartons(:,1:npartons) = plab(:,4:n)

      R = 0.7d0
      palg = 1d0
      call buildjets(npartons,ppartons,R,palg,pj,njets)

      do i=1,njets
         Etj1 = Et(pj(:,i))
         etaj1 = eta(pj(:,i))

         if(Etj1.lt.6d0) cycle
         if(etaj1.gt.3d0.or.etaj1.lt.-1d0) cycle

         call filld('Etj1',Etj1,dsig)
         
         call filld('etaj1',etaj1,dsig)
      enddo

     
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
!      call fastjetppgenkt(ptrack,ntracks,r,palg,pjet,njets)
      call fastjetppktetscheme(ptrack,ntracks,r,pjet,njets)

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

      function m2(p)
      implicit none
      double precision m2, p(0:3), norm

      m2 = p(0)**2 - norm(p)**2

      end

      function pt2(p)
      implicit none
      double precision pt2, p(0:3)

      pt2 = p(1)**2 + p(2)**2
      end

      function Et(p)
      implicit none
      double precision Et, p(0:3), m2, pt2

      Et = sqrt(m2(p) + pt2(p))
      end
      
