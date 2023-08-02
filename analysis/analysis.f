      subroutine init_histo
      implicit none
      include 'pwhg_bookhist-multi.h'
      
      call inihists
      call bookupeqbins('sigtot', 1d0, 0d0,1d0) 
      call bookupeqbins('sig(log(ptj1) > -5)', 1d0, 0d0, 1d0)
      call bookupeqbins('Q2', 5d2, 0d0,1d4)  
      call bookupeqbins('y', 5d-2, 0d0,1d0) 
      call bookupeqbins('x', 5d-2, 0d0,1d0)       
      call bookupeqbins('logptj1',0.1875d0, -10d0, 2d0)
!      call bookupeqbins('log(broadening)',0.5d0,-10d0,0d0)           
      end subroutine
      
      subroutine user_analysis(n,dsigma,x,y,Qsq)
      use parameters
      implicit none
      integer n
      double precision dsigma, x, y, Qsq, plab(0:3,n), pbreit(0:3,n)
      double precision broad, kt
      integer part

      if(dsigma.eq.0d0) return

      call filld('sigtot',0.5d0, dsigma)
!       First fill the DIS variables which are passed
      call filld('x',x,dsigma)
      call filld('y',y,dsigma)
      call filld('Q2',Qsq,dsigma)
      
      ! Now fill the momenta
      if(n.eq.4) then ! Born kinematics
        pbreit = pbornbreit
        plab = pbornlab
      elseif(n.eq.5) then
        pbreit = prealbreit
        plab = preallab
      elseif(n.eq.6) then
        pbreit = prrealbreit
        plab = prreallab
      else
        stop 'Wrong n in analysis'
      endif      

      broad = 0d0
      do part = 4,n
        broad = broad + kt(pbreit(:,part))
      enddo
      broad = broad /sqrt(Qsq)

      ! For plots with Silvia
      broad = log(max(broad * 0.5d0, 1d-100)) 
!      if (n.eq.5) print*, broad, log(max(broad,1d-100))
      if(broad.gt.-5d0) call filld('sig(log(ptj1) > -5)',0.5d0,dsigma)
      call filld('logptj1',broad,dsigma)
!      call filld('log(broadening)',log(max(broad,1d-100)),dsigma)
      end subroutine
      
      double precision function kt(p)
        implicit none
        double precision p(0:3)

        kt = sqrt(max(p(1)**2 + p(2)**2,0d0))
      end
