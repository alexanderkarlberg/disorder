      subroutine define_histograms
      implicit none
      include 'pwhg_bookhist-multi.h'
      
      call bookupeqbins('lntauzQ',0.25d0,-30d0,0d0)
      call bookupeqbins('lnBzQ',0.25d0,-20d0,0d0)

      end subroutine
      
      subroutine user_analysis(n,dsig,x,y,Q2)
      use mod_parameters
      use mod_analysis
      implicit none
      double precision dsig(maxscales), x, y, Q2
      integer n
      double precision tauzQ, BzQ, Qvec(4), p_evs(4,n-3), safelog
      double precision lnBzQ, lntauzQ
      integer i
      double precision, parameter :: cut = -2

      interface
         subroutine evs_tauzq_DIS(pcur, Q, var)
         integer, parameter  :: dp = selected_real_kind(15)
         real(dp), intent(in)  :: pcur(:,:)
         real(dp), intent(in)  :: Q(1:4)
         real(dp), intent(out) :: var
         end       
      end interface
      interface
         subroutine evs_brdzq_DIS(pcur, Q, var)
         integer, parameter  :: dp = selected_real_kind(15)
         real(dp), intent(in)  :: pcur(:,:)
         real(dp), intent(in)  :: Q(1:4)
         real(dp), intent(out) :: var
         end
      end interface

      if(n.eq.4) return         ! All event shapes are zero

!     Fill final state particles. Convention is (1,2,3,4) = (x,y,z,E)
!     of the EvsLib
      do i = 4,n
         p_evs(1:3,i-3) = pbreit(1:3,i)
         p_evs(4,i-3)   = pbreit(0,i)
      enddo

      Qvec(1:3) = pbreit(1:3,1) - pbreit(1:3,3)
      Qvec(4)   = pbreit(0,1) - pbreit(0,3)
      
      call evs_tauzQ_DIS(p_evs,Qvec,tauzQ)
      call evs_brdzQ_DIS(p_evs,Qvec,BzQ)

      lnBzQ   = safelog(BzQ)
      lntauzQ = safelog(tauzQ)

!      if((npow1.lt.4.and.safelog(BzQ).gt.cut) ! We fill above the cut
!     $     .or.
!     $     (npow1.ge.4.and.safelog(BzQ).lt.cut) ! We fill below the cut
!     $     ) then                    
         if(tauzQ.gt.0d0) call filld('lntauzQ',log(tauzQ),dsig)
         if(BzQ.gt.0d0)   call filld('lnBzQ',log(BzQ),dsig)
!      end if

      end subroutine
      
      subroutine evs_tauzq_DIS(parr,Qvec,tau_zQ)
      implicit none 
      integer, parameter  :: dp = selected_real_kind(15)
      real(dp), parameter :: zero = 0.0_dp, one = 1.0_dp
      real(dp) :: half = 0.5_dp 
      real(dp), intent(in)  :: parr(:,:)
      real(dp), intent(in) :: Qvec(4)
      real(dp)          :: pcur(size(parr,dim=1), size(parr,dim=2))
      real(dp)          :: Ecur, Q
      real(dp) :: thr, a(3), EoverhalfQ, brd, rho, ptot(4)
      real(dp) :: cpr, costheta
      real(dp), intent(out) :: tau_zQ
      integer           :: n, i, j
      real(dp) :: temp

      Q = sqrt(dot_product(Qvec(:3),Qvec(:3))-Qvec(4)*Qvec(4))
      
      n = 0; Ecur = zero
      do i = 1, size(parr,dim=2)
         if (parr(3,i)*Qvec(3) > zero) then
            !print*, i, ' is in the current', parr(:,i)
             n = n+1
             pcur(:,n) = parr(:,i)
          else
             !print*, i, ' is in the remnant', parr(:,i)
          end if
       end do
!------ tau_z ----
       if (n == 0) then
          thr = one
       else
          thr = sum(abs(pcur(3,:n))) / (half * Q)
       end if
       
       tau_zQ = one - thr

       end subroutine evs_tauzq_DIS      

      subroutine evs_brdzQ_DIS(parr, Qvec, brd_zE)
      implicit none
      integer, parameter  :: dp = selected_real_kind(15)
      real(dp), parameter :: zero = 0.0_dp, one = 1.0_dp
      real(dp) :: half = 0.5_dp
      real(dp), intent(in)  :: parr(:,:)
      real(dp), intent(in)  :: Qvec(4)
      real(dp)          :: pcur(size(parr,dim=1), size(parr,dim=2))
      real(dp)          :: Ecur
      real(dp) :: thr, a(3), EoverhalfQ, brd, rho, ptot(4), Q
      real(dp) :: cpr, costheta
      real(dp), intent(out) :: brd_zE
      integer           :: n, i, j
      real(dp) :: temp

      Q = sqrt(dot_product(Qvec(:3),Qvec(:3))-Qvec(4)*Qvec(4))

      n = 0; Ecur = zero
      do i = 1, size(parr,dim=2)
          if (parr(3,i)*Qvec(3) > zero) then
             n = n+1
             pcur(:,n) = parr(:,i)
!             Ecur = Ecur + parr(4,i)
          end if
       end do
!------ brd_zQ ---
       select case(n)
      case(0)
         brd_zE = zero             !es_undefined_value
      case default
         brd = zero
         do i = 1, n
!brd = brd + sqrt(mod3sq(pcur(:3,i) - a*dot_product(a,pcur(:3,i))))
!--   assume z axis is along 3 direction!
            brd = brd + sqrt(pcur(1,i)**2 + pcur(2,i)**2)
         end do
         brd = brd / Q
         brd_zE = brd
      end select    
      
      
      end subroutine evs_brdzQ_DIS

      double precision function safelog(x)
      implicit none
      double precision x

      if(x.gt.1d-100) then
         safelog = log(x)
      else
         safelog = -1d100
      endif
      end function
