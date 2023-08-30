      subroutine define_histograms
      implicit none
      include 'pwhg_bookhist-multi.h'
      
      call bookupeqbins('sigma',1d0,0d0,1d0)

      end subroutine
      
      subroutine user_analysis(n,dsig,x,y,Q2)
      use mod_parameters
      use mod_analysis
      implicit none
      double precision dsig(maxscales), x, y, Q2
      integer n

      call filld('sigma',0.5d0,dsig)

      end subroutine
      
