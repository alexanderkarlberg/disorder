      subroutine define_histograms
      implicit none
      include 'pwhg_bookhist-multi.h'
      
      call bookupeqbins('sig_Ecurr_cut',1d0,0d0,1d0)

      end subroutine
      
      subroutine user_analysis(n,dsig,x,y,Q2)
      use mod_analysis
      implicit none
      double precision dsig(maxscales), x, y, Q2
      integer n

      double precision Ecut, q(0:3)
      integer npartons
      logical pass_Ecur, CutDIS_Ecur

      q = pbreit(:,1) - pbreit(:,3)

      npartons = n-3

      Ecut = sqrt(Q2)/10d0

      pass_Ecur = CutDIS_Ecur(npartons, pbreit(:,4:), q, Ecut)

      if (pass_Ecur) return

      call filld('sig_Ecurr_cut',0.5d0,dsig)

      end subroutine
      

!     For Ecut < Q*(sqrt(2)-1) this selects only events with no
!     particles in the current hemisphere For Ecut < 0.5*Q, only real
!     radiation is involved, hence the cut selects DIS two-jet events
      logical function CutDIS_Ecur(n,p,q,Ecut) 
      implicit none
      integer n
      double precision p(0:3,n),q(0:3), Ecut
!----------------------------------------
      integer :: i
      double precision :: Ecur
      
      Ecur = 0d0
      do i=1,n
          if (dot_product(p(1:3,i),q(1:3))>0d0) then
            Ecur = Ecur + p(0,i)
         end if
      end do
      cutDIS_Ecur = Ecur<Ecut    
      
      end function CutDIS_Ecur
