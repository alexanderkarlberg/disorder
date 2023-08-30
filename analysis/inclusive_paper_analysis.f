      subroutine define_histograms
      implicit none
      include 'pwhg_bookhist-multi.h'
      double precision xmin, logxmin,binsize

      xmin = 9.8456206679269054D-4 ! For Q = 10 and setup in paper
      logxmin = log(xmin)
!     We want 25 bins
      binsize = -logxmin/25d0
      call bookupeqbins('dsigma/dlogx',binsize,logxmin,0d0)
      call bookupeqbins('dsigma/dx',1d0/25d0,0d0,1d0)

      end subroutine
      
      subroutine user_analysis(n,dsig,x,y,Q2)
      use mod_parameters
      use mod_analysis
      implicit none
      double precision dsig(maxscales), x, y, Q2
      integer n

      call filld('dsigma/dlogx',log(x),dsig)
      call filld('dsigma/dx',x,dsig)

      end subroutine
      
