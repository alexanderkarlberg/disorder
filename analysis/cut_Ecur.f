      subroutine define_histograms
      implicit none
      include 'pwhg_bookhist-multi.h'
      
      call bookupeqbins('sig_total',1d0,0d0,1d0)
      call bookupeqbins('sig_no_Ecurr',1d0,0d0,1d0)
      call bookupeqbins('sig_Ecurr_cut',1d0,0d0,1d0)
      call bookupeqbins('Ecur',0.1d0,-0.1d0,5d0)

      end subroutine
      
      subroutine user_analysis(n,dsig,x,y,Q2)
      use mod_analysis
      implicit none
      double precision dsig(maxscales), x, y, Q2
      integer n

      double precision Ecut, Ecur, q(0:3)
      integer npartons
      logical Ecur_lt_Ecut, CutDIS_Ecur

!     Cut on electron energy

      if(plab(0,3).lt.11d0) return
!     Always fill total
      call filld('sig_total',0.5d0,dsig)

      q = pbreit(:,1) - pbreit(:,3)

      npartons = n-3

      Ecut = sqrt(Q2)/10d0

      Ecur_lt_Ecut = CutDIS_Ecur(npartons, pbreit(:,4:), q, Ecut, Ecur)

      if (.not.Ecur_lt_Ecut) call filld('sig_Ecurr_cut',0.5d0,dsig)

      if (Ecur .eq. 0d0) call filld('sig_no_Ecurr',0.5d0,dsig)

      call filld('Ecur', Ecur, dsig)

      end subroutine
      

!     For Ecut < Q*(sqrt(2)-1) this selects only events with no
!     particles in the current hemisphere For Ecut < 0.5*Q, only real
!     radiation is involved, hence the cut selects DIS two-jet events
      logical function CutDIS_Ecur(n,p,q,Ecut,Ecur) 
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
