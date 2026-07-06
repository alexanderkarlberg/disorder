      subroutine define_histograms
      use mod_parameters
      implicit none
      include 'pwhg_bookhist-multi.h'
      double precision bins
      double precision logxmin, xbinsize
      double precision logQmax, Qbinsize, logQmin
      double precision logymin, logymax, ybinsize

      bins = 50d0


      logQmax = log(Q2max)
      logQmin = log(Q2min)
      Qbinsize = (logQmax-logQmin)/bins 

      call bookupeqbins('diff_hist:dsigma_dlogQ2',Qbinsize,logQmin,logQmax)

      logxmin = log(xmin)
      xbinsize = -logxmin/bins ! xmax = 1d0

      call bookupeqbins('diff_hist:dsigma_dlogx',xbinsize,logxmin,0d0)

      logymax = 0d0
      logymin = log(Q2min/S/xmax)
      ybinsize = (logymax-logymin)/bins

      !print*, 'logymin = ', logymin
      !print*, 'logymax = ', logymax
      !stop

      call bookupeqbins('diff_hist:dsigma_dlogy',ybinsize,logymin,logymax)

      call bookupeqbins('diff_hist:sigma',1d0,0d0,1d0)
      end subroutine
      
      subroutine user_analysis(n,dsig,x,y,Q2)
      use mod_parameters
      use mod_analysis
      implicit none
      double precision dsig(maxscales), x, y, Q2
      integer n

      call filld('diff_hist:dsigma_dlogQ2',log(Q2),dsig)

      call filld('diff_hist:dsigma_dlogx',log(x),dsig)

      call filld('diff_hist:dsigma_dlogy',log(y),dsig)

      call filld('diff_hist:sigma',0.5d0,dsig)

      end subroutine
      
