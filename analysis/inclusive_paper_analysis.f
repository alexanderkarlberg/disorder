      subroutine define_histograms
      implicit none
      include 'pwhg_bookhist-multi.h'
      double precision bins
      double precision xmin, logxmin, xbinsize
      double precision Qmax, logQmax, Qbinsize

      bins = 25d0
      xmin = 9.8456206679269054D-4 ! For Q = 10 and setup in paper
      logxmin = log(xmin)
      xbinsize = -logxmin/bins ! xmax = 1d0

      call bookupeqbins('dsigma/dlogx',xbinsize,logxmin,0d0)
      call bookupeqbins('dsigma/dx',1d0/bins,0d0,1d0)

      Qmax = sqrt(1015.68D0)    ! For x = 0.01 and setup in paper
      logQmax = log(Qmax) 
      Qbinsize = logQmax/bins ! Qmin = 1d0

      call bookupeqbins('dsigma/dlogQ^2',Qbinsize,0d0,logQmax)
      call bookupeqbins('dsigma/dQ^2',Qmax/bins,0d0,Qmax)

      end subroutine
      
      subroutine user_analysis(n,dsig,x,y,Q2)
      use mod_parameters
      use mod_analysis
      implicit none
      double precision dsig(maxscales), x, y, Q2
      integer n

      call filld('dsigma/dlogx',log(x),dsig)
      call filld('dsigma/dx',x,dsig)

      call filld('dsigma/dlogQ^2',log(sqrt(Q2)),dsig)
      call filld('dsigma/dQ^2',sqrt(Q2),dsig)

      end subroutine
      
