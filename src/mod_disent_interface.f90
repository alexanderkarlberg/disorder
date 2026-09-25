!----------------------------------------------------------------------
! Callbacks handed to DISENTFULL (src/libdisent.f) in p2b mode, plus
! the small helpers they use to translate between DISENT's conventions
! and ours. DISENT stores momenta as P(1:4,i) = (px,py,pz,E) in the
! Breit frame with
!   i = 1: incoming parton, 2-4: outgoing partons,
!   i = 5: exchanged momentum q, 6/7: incoming/outgoing lepton,
! whereas we use p(0:3,i) = (E,px,py,pz) ordered as incoming lepton,
! incoming parton, outgoing lepton, outgoing partons.
module mod_disent_interface
  use hoppet, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF
  use mod_parameters
  use mod_matrix_element
  use mod_phase_space
  use mod_analysis
  use types
  implicit none

  private
  public :: user, dis_cuts, disent_muf, DOT
  public :: disent_kinematics, disent_to_momenta
  public :: disent_imuf, disent_imur

  ! DISENT returns its scale-variation weights for three values of
  ! muF (1, 2 and 1/2 times the central scale, see MUF in
  ! KPFUNS_SCL_VAR in libdisent.f). These arrays map each of our
  ! maxscales points (scales_mur/scales_muf in mod_parameters) onto
  ! the corresponding DISENT muF point and the muR point (1, 2, 1/2)
  ! at which alphas is evaluated.
  integer, parameter :: disent_imuf(maxscales) = (/1, 2, 3, 2, 3, 1, 1/)
  integer, parameter :: disent_imur(maxscales) = (/1, 2, 3, 1, 1, 2, 3/)

contains

 subroutine dis_cuts(s,xminl,xmaxl,Q2minl,Q2maxl,yminl,ymaxl)
   real(dp), intent(in) :: s
   real(dp), intent(out) :: xminl,xmaxl,Q2minl,Q2maxl,yminl,ymaxl

   xminl = xmin
   xmaxl = xmax

   Q2minl = Q2min
   Q2maxl = Q2max

   ! Need to be careful since disent complains if all three variables
   ! are constrained (even if they are compatible).
   if(xmin.eq.xmax.and.Q2min.eq.Q2max) then
      yminl = 0d0
      ymaxl = 1d0
   else
      yminl = ymin
      ymaxl = ymax
   endif

 end subroutine dis_cuts

  subroutine disent_muf(p,s,muF2)
   real(dp), intent(in) :: p(4,7),s
   real(dp), intent(out) :: muF2
   real(dp) :: dummy,x,y,Q2,Qval

   Q2=ABS(DOT(P,5,5))
   Qval=sqrt(Q2)
   Y=DOT(P,1,5)/DOT(P,1,6)
   X=Q2/(y*S)

   call muR_muF(x,y,Qval,dummy,muF2)
   ! muF2 at this point is muf * mu

   muF2 = muF2**2 / Q2

 end subroutine disent_muf

 ! Born variables x, y, Q2 of a DISENT event, given the momentum
 ! fraction eta = 2 p1.p6/S of the incoming parton.
 subroutine disent_kinematics(p,eta,x,y,Q2)
   real(dp), intent(in)  :: p(4,7), eta
   real(dp), intent(out) :: x, y, Q2

   Q2=ABS(DOT(P,5,5))
   X=ETA*Q2/(2*DOT(P,1,5))
   Y=DOT(P,1,5)/DOT(P,1,6)
 end subroutine disent_kinematics

 ! Reorder an n-parton DISENT momentum set into our convention,
 ! pout(0:3,1:n+2), still in the Breit frame.
 subroutine disent_to_momenta(n,p,pout)
   integer,  intent(in)  :: n
   real(dp), intent(in)  :: p(4,7)
   real(dp), intent(out) :: pout(0:3,n+2)
   integer :: i

   pout(:,1) = cshift(p(:,6),-1) ! Incoming lepton
   pout(:,2) = cshift(p(:,1),-1) ! Incoming quark
   pout(:,3) = cshift(p(:,7),-1) ! Outgoing lepton
   do i = 2, n
      pout(:,i+2) = cshift(p(:,i),-1) ! Outgoing partons
   enddo
 end subroutine disent_to_momenta

 ! user-defined event analysis for disent
 ! scale is muF, hence it needs to be divided by xmuf if used as a renormalisation scale
 subroutine user(N,NA,NT,P,S,WEIGHT,SCALE2)
   implicit none
   integer, intent(in) :: N, NA, NT
   real(dp), intent(in) :: s, p(4,7), weight(-6:6),scale2
   real(dp) :: scale
   LOGICAL SCALE_VAR
   DOUBLE PRECISION SCL_WEIGHT(3,-6:6)
   COMMON/cSCALE_VAR/SCL_WEIGHT, SCALE_VAR

   integer isc

   double precision dsig(maxscales), totwgt, wgt_array(maxscales,-6:6)
   double precision, save ::  pdfs(maxscales,-6:6), eta, Q2, Qval, x, y, as2pi(maxscales)

   double precision, save :: etasave = -100d0
   ! For p2b
   logical, save :: recompute = .true.
   double precision, save ::  p2blab(0:3,2+2), p2bbreit(0:3,2+2), Qlab(0:3)

   if (n.eq.0) then ! Disent is done with one event cycle
      call pwhgaccumup
      recompute = .true. ! Signals that next time we have a new event cycle
      return
   endif

   if(p2b.and.n.eq.2) return ! If we do p2b we get the Born and
                             ! virtuals from the structure functions

   ! The following lines are invoked if the user specify only the
   ! nlocoeff or nnlocoeff on the command line
   if(order_max.le.2.and.NA.ge.2) return ! Disregard O(αS**2) if we are doing NLO
   if(order_max.le.1.and.NA.ge.1) return ! Disregard O(αS) if we are doing LO

   if(order_min.gt.2.and.NA.lt.2) return ! Disregard O(αS) if we are doing NNLO coeff
   if(order_min.gt.1.and.NA.lt.1) return ! Disregard O(1) if we are doing NLO coeff

   ! It looks like, in a given set of calls (ie born + real + ...) eta
   ! can change, but x,y,Q2 stay the same. Eta however only changes a
   ! few times, so it is worth recomputing eta (which is cheap) and
   ! saving the pdfs (since they are more expensive to recompute).
   ETA=2*DOT(P,1,6)/S
   SCALE = SQRT(SCALE2)
   if(scaleuncert) then
      ! First we transfer the 3 weights from DISENT (muF variations)
      ! into the full array. We also correct alpha_em which is
      ! hardcoded to 1/137 in DISENT.
      do isc = 1,maxscales
         wgt_array(isc,:) = scl_weight(disent_imuf(isc),:) * weight(:) * (137.0d0 * alpha_em)**2
      enddo

      if(recompute) then
         ! get x and Q2
         call disent_kinematics(p,eta,x,y,Q2)
         Qval=sqrt(Q2)
         do isc = 1,3
            as2pi(isc) = alphasLocal(xmur*scales_mur(isc)*scale/xmuf*Qval)/pi * 0.5d0
         enddo
         ! Then copy the as2pi into the full array
         do isc = 4,maxscales
            as2pi(isc) = as2pi(disent_imur(isc))
         enddo
         call p2bmomenta(x,y,Q2,p2bbreit,p2blab)
         Qlab(:) = p2blab(:,1) - p2blab(:,3)
         recompute = .false.
      endif

      if(eta.ne.etasave) then
         etasave = eta
         do isc = 1,3
            call hoppetEval(eta,scales_muf(isc)*scale*Qval,pdfs(isc,:))
         enddo
         ! Then copy the PDFs into the full array
         do isc = 4,maxscales
            pdfs(isc,:) = pdfs(disent_imuf(isc),:)
         enddo
      endif

      ! First we dress the weights with the PDFs
      do isc = 1, maxscales
         dsig(isc) = dot_product(wgt_array(isc,:),pdfs(isc,:))
      enddo
      ! Now we dress them with αS and scale compensation
      if(order_max.eq.3.and.NA.eq.1) then ! Doing NLO in DISENT and this is the LO term. Include scale compensation.
         dsig(:) = dsig(:) * (as2pi(:) + two * as2pi(:)**2 * b0 * log(xmur*scales_mur(:)*scale/xmuf))
      else
         dsig(:) = dsig(:) * as2pi(:)**NA
      endif
      dsig = dsig * ncall2
   else
      if(recompute) then
         ! get x and Q2
         call disent_kinematics(p,eta,x,y,Q2)
         Qval=sqrt(Q2)

         as2pi(1) = alphasLocal(xmur*scale/xmuf*Qval)/pi * 0.5d0
         call p2bmomenta(x,y,Q2,p2bbreit,p2blab)
         Qlab(:) = p2blab(:,1) - p2blab(:,3)
         recompute = .false.
      endif

      if(eta.ne.etasave) then
         etasave = eta
         call hoppetEval(eta,scale*Qval,pdfs(1,:))
      endif


      if(order_max.eq.3.and.NA.eq.1) then ! Doing NLO in DISENT and this is the LO term. Include scale compensation.
         totwgt = dot_product(weight,pdfs(1,:))*(as2pi(1) + two * as2pi(1)**2 * b0 * log(xmur*scale/xmuf))
      else
         totwgt = dot_product(weight,pdfs(1,:))*as2pi(1)**na
      endif
      ! We correct alpha_em which is hardcoded to 1/137 in DISENT.
      dsig(1) = totwgt * ncall2 * (137.0d0 * alpha_em)**2
   endif

   if(xmin.eq.xmax) dsig = dsig / x ! Because of convention in DISENT where it currently returns x dσ/dx when x is fixed.

   ! First we transfer the DISENT momenta to our convention
   pbornbreit = 0 ! pborn(0:3,2+2)
   prealbreit = 0 ! preal(0:3,3+2)
   prrealbreit = 0 ! prreal(0:3,4+2)

   if(n.eq.2) then ! Born kinematics
      call disent_to_momenta(n,p,pbornbreit)
      call mbreit2lab(n+2,Qlab,pbornbreit,pbornlab,.true.)
   elseif(n.eq.3) then
      call disent_to_momenta(n,p,prealbreit)
      call mbreit2lab(n+2,Qlab,prealbreit,preallab,.true.)
   elseif(n.eq.4) then
      call disent_to_momenta(n,p,prrealbreit)
      call mbreit2lab(n+2,Qlab,prrealbreit,prreallab,.true.)
   else
      print*, 'n = ', n
      stop 'Wrong n in user routine of DISENT'
   endif
   call analysis(n+2, dsig, x, y, Q2)
   ! projection-to-Born analysis call
   pbornbreit = p2bbreit
   pbornlab   = p2blab
   call analysis(2+2, -dsig, x, y, Q2)
 END subroutine user

 ! Taken directly from DISENT
 FUNCTION DOT(P,I,J)
   IMPLICIT NONE
   !---RETURN THE DOT PRODUCT OF P(*,I) AND P(*,J)
   INTEGER I,J
   DOUBLE PRECISION DOT,P(4,7)
   DOT=P(4,I)*P(4,J)-P(3,I)*P(3,J)-P(2,I)*P(2,J)-P(1,I)*P(1,J)
 END FUNCTION DOT

end module mod_disent_interface
