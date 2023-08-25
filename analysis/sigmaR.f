      subroutine define_histograms
      use mod_parameters
      implicit none
      include 'pwhg_bookhist-multi.h'
      
      integer xnbins, Q2nbins
      
      real * 8 Q2binsize, logQ2min, logQ2max
      real * 8 xbinsize, logxmin, logxmax
      
      ! Bins in Q^2 and as a function of x
      xnbins = 50
      xbinsize = 0.5d0 / dble(xnbins)
      call bookupeqbins('sigma_r_Q2_5.0',   xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_7.0',   xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_9.0',   xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_11.0',   xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_13.0',   xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_16.0',   xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_20.0',   xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_32.0',   xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_40.0',   xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_50.0',   xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_65.0',   xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_85.0',   xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_110.0',  xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_140.0',  xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_185.0',  xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_240.0',  xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_310.0',  xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_410.0',  xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_530.0',  xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_710.0',  xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_900.0',  xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_1300.0', xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_1800.0', xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_2500.0', xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_3500.0', xbinsize,0d0,0.5d0)
      call bookupeqbins('sigma_r_Q2_15000.0',xbinsize,0d0,0.5d0)
      
      ! Bins in Q^2 and as a function of log(x)
      logxmax = log(0.5d0)!0d0
      logxmin = -12d0
      xbinsize = (logxmax - logxmin) / dble(xnbins)
      call bookupeqbins('log_sigma_r_Q2_5.0',   xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_7.0',   xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_9.0',   xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_11.0',   xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_13.0',   xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_16.0',   xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_20.0',   xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_32.0',   xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_40.0',   xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_50.0',   xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_65.0',   xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_85.0',   xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_110.0',  xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_140.0',  xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_185.0',  xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_240.0',  xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_310.0',  xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_410.0',  xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_530.0',  xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_710.0',  xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_900.0',  xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_1300.0', xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_1800.0', xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_2500.0', xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_3500.0', xbinsize,logxmin,logxmax)
      call bookupeqbins('log_sigma_r_Q2_15000.0',xbinsize,logxmin,logxmax)

      Q2nbins = 50
!     Bins in x and as a function of log(Q^2)
      logQ2min = log(3d0)/log(10d0)
      logQ2max = log(s)/log(10d0)
      Q2binsize = (logQ2max - logQ2min) / dble(Q2nbins)
      
      call bookupeqbins('log_sigma_r_x_0.00032',Q2binsize,logQ2min,logQ2max)
      call bookupeqbins('log_sigma_r_x_0.0005', Q2binsize,logQ2min,logQ2max)
      call bookupeqbins('log_sigma_r_x_0.0008', Q2binsize,logQ2min,logQ2max)
      call bookupeqbins('log_sigma_r_x_0.0013', Q2binsize,logQ2min,logQ2max)
      call bookupeqbins('log_sigma_r_x_0.002',  Q2binsize,logQ2min,logQ2max)
      call bookupeqbins('log_sigma_r_x_0.0032', Q2binsize,logQ2min,logQ2max)
      call bookupeqbins('log_sigma_r_x_0.005',  Q2binsize,logQ2min,logQ2max)
      call bookupeqbins('log_sigma_r_x_0.008',  Q2binsize,logQ2min,logQ2max)
      call bookupeqbins('log_sigma_r_x_0.013',  Q2binsize,logQ2min,logQ2max)
      call bookupeqbins('log_sigma_r_x_0.02',   Q2binsize,logQ2min,logQ2max)
      call bookupeqbins('log_sigma_r_x_0.032',  Q2binsize,logQ2min,logQ2max)
      call bookupeqbins('log_sigma_r_x_0.05',   Q2binsize,logQ2min,logQ2max)
      call bookupeqbins('log_sigma_r_x_0.08',   Q2binsize,logQ2min,logQ2max)
      call bookupeqbins('log_sigma_r_x_0.13',   Q2binsize,logQ2min,logQ2max)
      call bookupeqbins('log_sigma_r_x_0.2',    Q2binsize,logQ2min,logQ2max)
      call bookupeqbins('log_sigma_r_x_0.32',   Q2binsize,logQ2min,logQ2max)
      call bookupeqbins('log_sigma_r_x_0.5',    Q2binsize,logQ2min,logQ2max)

      end subroutine
      
      subroutine user_analysis(n,dsigma,x,y,Q2)
      use mod_analysis
      implicit none
      integer n
      double precision dsigma(maxscales), x, y, Q2
      double precision sigma_r, delx
      double precision, parameter :: pi = 4d0*atan(1d0)

    !     Reduced cross section
      sigma_r =  x * Q2**2
     $     / (2d0 * pi * alpha_em**2 * (1d0 + (1d0 - y)**2))
      sigma_r = sigma_r / 1d6 ! convert to mb

!     Now we fill histograms first in Q2 binning
      if(Q2.lt.4d0) then
         return
      elseif(Q2.lt.5d0) then
         call filld('sigma_r_Q2_5.0',    x, sigma_r * dsigma / 1d0)
         call filld('log_sigma_r_Q2_5.0',    log(x), sigma_r / x * dsigma /  1d0)
      elseif(Q2.lt.7d0) then
         call filld('sigma_r_Q2_7.0',    x, sigma_r * dsigma / 2d0)
         call filld('log_sigma_r_Q2_7.0',    log(x), sigma_r / x * dsigma /  2d0)
      elseif(Q2.lt.9d0) then
         call filld('sigma_r_Q2_9.0',    x, sigma_r * dsigma / 2d0)
         call filld('log_sigma_r_Q2_9.0',    log(x), sigma_r / x * dsigma /  2d0)
      elseif(Q2.lt.11d0) then
         call filld('sigma_r_Q2_11.0',    x, sigma_r * dsigma / 2d0)
         call filld('log_sigma_r_Q2_11.0',    log(x), sigma_r / x * dsigma /  2d0)
      elseif(Q2.lt.13d0) then
         call filld('sigma_r_Q2_13.0',    x, sigma_r * dsigma / 2d0)
         call filld('log_sigma_r_Q2_13.0',    log(x), sigma_r / x * dsigma /  2d0)
      elseif(Q2.lt.16d0) then
         call filld('sigma_r_Q2_16.0',    x, sigma_r * dsigma / 3d0)
         call filld('log_sigma_r_Q2_16.0',    log(x), sigma_r / x * dsigma /  3d0)
      elseif(Q2.lt.20d0) then
         call filld('sigma_r_Q2_20.0',    x, sigma_r * dsigma / 4d0)
         call filld('log_sigma_r_Q2_20.0',    log(x), sigma_r / x * dsigma /  4d0)
      elseif(Q2.lt.32d0) then
         call filld('sigma_r_Q2_32.0',    x, sigma_r * dsigma / 12d0)
         call filld('log_sigma_r_Q2_32.0',    log(x), sigma_r / x * dsigma /  12d0)
      elseif(Q2.lt.40d0) then
         call filld('sigma_r_Q2_40.0',    x, sigma_r * dsigma / 8d0)
         call filld('log_sigma_r_Q2_40.0',    log(x), sigma_r / x * dsigma /  8d0)
      elseif(Q2.lt.50d0) then
         call filld('sigma_r_Q2_50.0',    x, sigma_r * dsigma / 10d0)
         call filld('log_sigma_r_Q2_50.0',    log(x), sigma_r / x * dsigma /  10d0)
      elseif(Q2.lt.65d0) then
         call filld('sigma_r_Q2_65.0',    x, sigma_r * dsigma / 15d0)
         call filld('log_sigma_r_Q2_65.0',    log(x), sigma_r / x * dsigma /  15d0)
      elseif(Q2.lt.85d0) then
         call filld('sigma_r_Q2_85.0',    x, sigma_r * dsigma / 20d0)
         call filld('log_sigma_r_Q2_85.0',    log(x), sigma_r / x * dsigma /  20d0)
      elseif(Q2.lt.110d0) then
         call filld('sigma_r_Q2_110.0',    x, sigma_r * dsigma / 25d0)
         call filld('log_sigma_r_Q2_110.0',    log(x), sigma_r / x * dsigma /  25d0)
      elseif(Q2.lt.140d0) then
         call filld('sigma_r_Q2_140.0',    x, sigma_r * dsigma / 30d0)
         call filld('log_sigma_r_Q2_140.0',    log(x), sigma_r / x * dsigma /  30d0)
      elseif(Q2.lt.185d0) then
         call filld('sigma_r_Q2_185.0',    x, sigma_r * dsigma / 45d0)
         call filld('log_sigma_r_Q2_185.0',    log(x), sigma_r / x * dsigma /  45d0)
      elseif(Q2.lt.240d0) then
         call filld('sigma_r_Q2_240.0',    x, sigma_r * dsigma / 55d0)
         call filld('log_sigma_r_Q2_240.0',    log(x), sigma_r / x * dsigma /  55d0)
      elseif(Q2.lt.310d0) then
         call filld('sigma_r_Q2_310.0',    x, sigma_r * dsigma / 70d0)
         call filld('log_sigma_r_Q2_310.0',    log(x), sigma_r / x * dsigma /  70d0)
      elseif(Q2.lt.410d0) then
         call filld('sigma_r_Q2_410.0',    x, sigma_r * dsigma / 100d0)
         call filld('log_sigma_r_Q2_410.0',    log(x), sigma_r / x * dsigma /  100d0)
      elseif(Q2.lt.530d0) then
         call filld('sigma_r_Q2_530.0',    x, sigma_r * dsigma / 120d0)
         call filld('log_sigma_r_Q2_530.0',    log(x), sigma_r / x * dsigma /  120d0)
      elseif(Q2.lt.710d0) then
         call filld('sigma_r_Q2_710.0',    x, sigma_r * dsigma / 180d0)
         call filld('log_sigma_r_Q2_710.0',    log(x), sigma_r / x * dsigma /  180d0)
      elseif(Q2.lt.900d0) then
         call filld('sigma_r_Q2_900.0',    x, sigma_r * dsigma / 190d0)
         call filld('log_sigma_r_Q2_900.0',    log(x), sigma_r / x * dsigma /  190d0)
      elseif(Q2.lt.1300d0) then
         call filld('sigma_r_Q2_1300.0',    x, sigma_r * dsigma / 400d0)
         call filld('log_sigma_r_Q2_1300.0',    log(x), sigma_r / x * dsigma /  400d0)
      elseif(Q2.lt.1800d0) then
         call filld('sigma_r_Q2_1800.0',    x, sigma_r * dsigma / 500d0)
         call filld('log_sigma_r_Q2_1800.0',    log(x), sigma_r / x * dsigma /  500d0)
      elseif(Q2.lt.2500d0) then
         call filld('sigma_r_Q2_2500.0',    x, sigma_r * dsigma / 700d0)
         call filld('log_sigma_r_Q2_2500.0',    log(x), sigma_r / x * dsigma /  700d0)
      elseif(Q2.lt.3500d0) then
         call filld('sigma_r_Q2_3500.0',    x, sigma_r * dsigma / 1000d0)
         call filld('log_sigma_r_Q2_3500.0',    log(x), sigma_r / x * dsigma /  1000d0)
      elseif(Q2.lt.15000d0) then
         call filld('sigma_r_Q2_15000.0',    x, sigma_r * dsigma / 11500d0)
         call filld('log_sigma_r_Q2_15000.0',    log(x), sigma_r / x * dsigma /  11500d0)
      endif

      if(x.lt.0.0002d0) then
         return
      elseif(x.lt.0.00032d0) then
         delx = 0.00032d0 - 0.0002d0
         call filld('log_sigma_r_x_0.00032',log(Q2)/log(10d0), sigma_r * dsigma / Q2 / delx)
      elseif(x.lt.0.0005d0) then
         delx = 0.0005d0 - 0.00032d0
         call filld('log_sigma_r_x_0.0005',log(Q2)/log(10d0), sigma_r * dsigma / Q2 / delx)
      elseif(x.lt.0.0008d0) then
         delx = 0.0008d0 - 0.0005d0
         call filld('log_sigma_r_x_0.0008',log(Q2)/log(10d0), sigma_r * dsigma / Q2 / delx)
      elseif(x.lt.0.0013d0) then
         delx = 0.0013d0 - 0.0008d0 
         call filld('log_sigma_r_x_0.0013',log(Q2)/log(10d0), sigma_r * dsigma / Q2 / delx)
      elseif(x.lt.0.002d0) then
         delx = 0.0020d0 - 0.0013d0
         call filld('log_sigma_r_x_0.002',log(Q2)/log(10d0), sigma_r * dsigma / Q2 / delx)
      elseif(x.lt.0.0032d0) then
         delx = 0.0032d0 - 0.0020d0
         call filld('log_sigma_r_x_0.0032',log(Q2)/log(10d0), sigma_r * dsigma / Q2 / delx)
      elseif(x.lt.0.005d0) then
         delx = 0.005d0 - 0.0032d0
         call filld('log_sigma_r_x_0.005',log(Q2)/log(10d0), sigma_r * dsigma / Q2 / delx)
      elseif(x.lt.0.008d0) then
         delx = 0.008d0 - 0.005d0
         call filld('log_sigma_r_x_0.008',log(Q2)/log(10d0), sigma_r * dsigma / Q2 / delx)
      elseif(x.lt.0.013d0) then
         delx = 0.013d0 - 0.008d0
         call filld('log_sigma_r_x_0.013',log(Q2)/log(10d0), sigma_r * dsigma / Q2 / delx)
      elseif(x.lt.0.02d0) then 
         delx = 0.02d0 - 0.013d0
        call filld('log_sigma_r_x_0.02',log(Q2)/log(10d0), sigma_r * dsigma / Q2 / delx)
      elseif(x.lt.0.032d0) then
         delx = 0.032d0 - 0.02d0
         call filld('log_sigma_r_x_0.032',log(Q2)/log(10d0), sigma_r * dsigma / Q2 / delx)
      elseif(x.lt.0.05d0) then
         delx = 0.05d0 - 0.032d0
         call filld('log_sigma_r_x_0.05',log(Q2)/log(10d0), sigma_r * dsigma / Q2 / delx)
      elseif(x.lt.0.08d0) then
         delx = 0.08d0 - 0.05d0
         call filld('log_sigma_r_x_0.08',log(Q2)/log(10d0), sigma_r * dsigma / Q2 / delx)
      elseif(x.lt.0.13d0) then
         delx = 0.13d0 - 0.08d0
         call filld('log_sigma_r_x_0.13',log(Q2)/log(10d0), sigma_r * dsigma / Q2 / delx)
      elseif(x.lt.0.2d0) then
         delx = 0.2d0 - 0.13d0
         call filld('log_sigma_r_x_0.2',log(Q2)/log(10d0), sigma_r * dsigma / Q2 / delx)
      elseif(x.lt.0.32d0) then
         delx = 0.32d0 - 0.2d0
         call filld('log_sigma_r_x_0.32',log(Q2)/log(10d0), sigma_r * dsigma / Q2 / delx)
      elseif(x.lt.0.5d0) then
         delx = 0.5d0 - 0.32d0
         call filld('log_sigma_r_x_0.5',log(Q2)/log(10d0), sigma_r * dsigma / Q2 / delx)
      endif


      
      end subroutine
      
