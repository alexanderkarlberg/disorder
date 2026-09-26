!----------------------------------------------------------------------
! Unit test of DISENT's three-parton tree-level matrix element MATTHR
! (src/libdisent.f), generalised to gamma/Z and W exchange with the
! couplings from parton_couplings (src/mod_matrix_element.f90).
!
! Run with disorder's own command line flags (process selection). For
! three-parton configurations built here in the Breit frame, with
! DISENT's conventions, it checks
!   - for photon exchange: that MATTHR is bit-for-bit identical to the
!     original DISENT expression (written out again below);
!   - for all processes: that M(i) = c2(i) QQ + c3(i) QQ3 and
!     M(0) = sum_i cg(i) GQ, with QQ, GQ the original photon structures;
!   - that in the initial- and final-state collinear limits the
!     parity-violating structure satisfies QQ3/QQ -> Y-/Y+ =
!     (1-(1-y)^2)/(1+(1-y)^2), so that the three-parton matrix element
!     factorises onto the Born cross section Y+ c2 + Y- c3 (which
!     test_matrix_element checks against the quark-parton model).
!----------------------------------------------------------------------
program test_matthr
  use types, only: dp
  use mod_parameters
  use mod_matrix_element, only: parton_couplings
  use test_utils
  implicit none
  ! DISENT's /COLFAC/ (renamed locally, cutoff clashes with mod_parameters)
  integer :: SCHEME_D, NF
  double precision :: CF, CA, TR, PI_D, PISQ, HF, CUTOFF_D, EQ(-6:6), SCALE_D
  common /COLFAC/ CF, CA, TR, PI_D, PISQ, HF, CUTOFF_D, EQ, SCALE_D, SCHEME_D, NF
  real(dp) :: P(4,7), M(-6:6), c2(-6:6), c3(-6:6), cg(6)
  real(dp) :: Q, y, xi, cth, phi, QQ, QQ3, GQ, Mexp, r, yratio
  real(dp), parameter :: Qs(3) = [15.0_dp, 60.0_dp, 110.0_dp]
  real(dp), parameter :: ys(3) = [0.15_dp, 0.5_dp, 0.85_dp]
  logical :: photon, ok_bits, ok_decomp
  integer :: i, iq, iy, k, ntest
  character(len=80) :: tag

  call set_parameters()
  photon = NC .and. noZ .and. .not. CC

  ! the constants DISENTFULL sets up
  CF = 4.0_dp/3.0_dp; CA = 3; TR = 0.5_dp; NF = nflav
  PI_D = atan(1d0)*4; PISQ = PI_D**2; HF = 0.5_dp
  EQ(0) = 0
  EQ(1) = -1D0/3
  EQ(2) = EQ(1) + 1
  do i = 1, 6
     if (i > 2) EQ(i) = EQ(i-2)
     EQ(-i) = -EQ(i)
  enddo

  ok_bits = .true.; ok_decomp = .true.; ntest = 0
  r = 0.1234_dp
  do iq = 1, size(Qs)
     do iy = 1, size(ys)
        do k = 1, 20
           ! pseudo-random but reproducible configurations
           r = mod(r * 97.0_dp + 0.3137_dp, 1.0_dp); xi  = 0.05_dp + 0.94_dp * r
           r = mod(r * 97.0_dp + 0.3137_dp, 1.0_dp); cth = -0.99_dp + 1.98_dp * r
           r = mod(r * 97.0_dp + 0.3137_dp, 1.0_dp); phi = 6.2831853_dp * r
           call three_parton(Qs(iq), ys(iy), xi, cth, phi, P)
           call MATTHR(P, M)
           call structures(P, QQ, QQ3, GQ)
           call disent_couplings(-DOT(P,5,5), EQ, c2, c3, cg)
           ntest = ntest + 1
           if (photon) then
              ! the original DISENT expressions
              do i = -6, 6
                 if (i == 0) cycle
                 if (abs(i) <= nflav) ok_bits = ok_bits .and. M(i) == EQ(i)**2*QQ
              enddo
              Mexp = 0
              do i = 1, NF
                 Mexp = Mexp + EQ(i)**2*GQ
              enddo
              ok_bits = ok_bits .and. M(0) == Mexp
           endif
           do i = -6, 6
              if (i == 0) cycle
              ok_decomp = ok_decomp .and. abs(M(i) - (c2(i)*QQ + c3(i)*QQ3)) &
                   & <= 1e-13_dp * (abs(c2(i)*QQ) + abs(c3(i)*QQ3))
           enddo
           Mexp = sum(cg(1:NF)) * GQ
           ok_decomp = ok_decomp .and. abs(M(0) - Mexp) <= 1e-13_dp * abs(Mexp)
        enddo
     enddo
  enddo
  if (photon) call check_true('photon exchange: MATTHR identical to original DISENT', ok_bits)
  call check_true('MATTHR = c2 QQ + c3 QQ3 (quarks), sum cg GQ (gluon)', ok_decomp)
  call check_true('tested 180 configurations', ntest == 180)

  ! collinear limits: QQ3/QQ -> Y-/Y+
  do iy = 1, size(ys)
     y = ys(iy)
     yratio = (1 - (1-y)**2) / (1 + (1-y)**2)
     Q = 40.0_dp
     ! initial-state collinear: parton 3 along the incoming parton
     call three_parton(Q, y, 0.3_dp, 1.0_dp - 1e-10_dp, 0.7_dp, P)
     call structures(P, QQ, QQ3, GQ)
     write(tag,'(a,f5.2)') 'initial-state collinear limit, y =', y
     call check_close(trim(tag)//': QQ3/QQ = Y-/Y+', QQ3/QQ, yratio, 1e-4_dp)
     ! final-state collinear: xi -> 1 (partons 2 and 3 collinear)
     call three_parton(Q, y, 1.0_dp - 1e-9_dp, 0.3_dp, 0.7_dp, P)
     call structures(P, QQ, QQ3, GQ)
     write(tag,'(a,f5.2)') 'final-state collinear limit, y =', y
     call check_close(trim(tag)//': QQ3/QQ = Y-/Y+', QQ3/QQ, yratio, 1e-4_dp)
  enddo

  call finish_tests()

contains

  ! A three-parton configuration in DISENT's layout and Breit frame
  ! (P(1:4,i) = (px,py,pz,E); 1 incoming parton, 2,3 outgoing partons,
  ! 5 = q, 6/7 incoming/outgoing lepton), with the leptons as in GENTWO,
  ! incoming parton momentum 1/xi times that of the Born parton, and
  ! parton 3 at polar angle acos(cth) (w.r.t. the incoming parton) and
  ! azimuth phi in the rest frame of p2 + p3.
  subroutine three_parton(Q, y, xi, cth, phi, P)
    real(dp), intent(in)  :: Q, y, xi, cth, phi
    real(dp), intent(out) :: P(4,7)
    real(dp) :: E, E1, W, beta, gam, sth, k3(4)
    P = 0
    E = Q/2
    E1 = E/xi
    P(3,1) = E1; P(4,1) = E1
    P(3,5) = -2*E
    P(1,6) = E/y*2*sqrt(1-y); P(3,6) = -E; P(4,6) = E/y*(2-y)
    P(1,7) = E/y*2*sqrt(1-y); P(3,7) =  E; P(4,7) = E/y*(2-y)
    W = sqrt(4*E*E1 - 4*E**2)
    sth = sqrt(max(0.0_dp, 1 - cth**2))
    k3 = [W/2*sth*cos(phi), W/2*sth*sin(phi), W/2*cth, W/2]
    beta = (E1 - 2*E) / E1
    gam = 1 / sqrt(1 - beta**2)
    P(:,3) = [k3(1), k3(2), gam*(k3(3) + beta*k3(4)), gam*(k3(4) + beta*k3(3))]
    P(:,2) = P(:,1) + P(:,5) - P(:,3)
  end subroutine three_parton

  ! DISENT's photon structures QQ and GQ, and the parity-violating QQ3
  subroutine structures(P, QQ, QQ3, GQ)
    real(dp), intent(in)  :: P(4,7)
    real(dp), intent(out) :: QQ, QQ3, GQ
    QQ=8*(4*PI_D/137)**2* &
         & (DOT(P,1,6)**2+DOT(P,1,7)**2+DOT(P,2,7)**2+DOT(P,2,6)**2) &
         & *16*PISQ*CF/(-4*DOT(P,2,3)*DOT(P,1,3)*DOT(P,5,5))
    QQ3=8*(4*PI_D/137)**2* &
         & (DOT(P,1,6)**2-DOT(P,1,7)**2+DOT(P,2,7)**2-DOT(P,2,6)**2) &
         & *16*PISQ*CF/(-4*DOT(P,2,3)*DOT(P,1,3)*DOT(P,5,5))
    GQ=8*(4*PI_D/137)**2* &
         & (DOT(P,3,6)**2+DOT(P,3,7)**2+DOT(P,2,7)**2+DOT(P,2,6)**2) &
         & *16*PISQ*TR/(-4*DOT(P,2,1)*DOT(P,3,1)*DOT(P,5,5))
  end subroutine structures

  real(dp) function DOT(P, I, J)
    real(dp), intent(in) :: P(4,7)
    integer, intent(in) :: I, J
    DOT=P(4,I)*P(4,J)-P(3,I)*P(3,J)-P(2,I)*P(2,J)-P(1,I)*P(1,J)
  end function DOT

end program test_matthr
