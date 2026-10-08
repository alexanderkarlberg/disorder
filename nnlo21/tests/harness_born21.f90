!-----------------------------------------------------------------------
! The 2+1 Born of nnlo21 with photon + Z (8 Oct 2026) against DISENT's
! MATTHR, at random Born points, with disorder's flags (-includeZ,
! -positron, -neutrino; ew31 set from them):
!  1. hard21's helicity classes: hk(0,1) : hk(0,2) = S : O, with S =
!     (p1.p6)^2 + (p2.p7)^2, O = (p1.p7)^2 + (p2.p6)^2 (DISENT labels);
!  2. quark Born: MATTHR's M(f) = QQ/(S + O) (S w1 + O w2), w = ew31_w(f)
!     in units of the photon coupling (sliced21's weights);
!  3. gluon Born: M(0) = GQ sum_q (w1 + w2)/2 (the parity-odd part, odd
!     under quark <-> antiquark, dropped as in MATTHR);
!  4. one loop per class: H^(1)/H^(0) of class S (O) = the slicing's
!     hard_fact plus DISENT's non-factorising one-loop part of that class,
!     (QQ +- QQ3)/2 over its Born (VIRTHR: C2 QQ + C3 QQ3, QQ3 from
!     VIRT3PV, with C2 +- C3 = w1, w2). hard21 has a step of ~1e-5
!     relative within 1.3e-3 of the region boundaries s_ij = s45 (NNLOJET's
!     displacement at v -> 1; DISENT is smooth there): fixed seed, DIAG
!     prints such points.
!-----------------------------------------------------------------------
program harness_born21
  use types, only: dp
  use mod_parameters, only: set_parameters, CC, noZ, positron, neutrino
  use hard21
  use ew31
  use mod_slicing_scet, only: hard_fact, scet_set_colour
  implicit none
  integer :: SCHEME_D, NF
  double precision :: CF, CA, TR, PI_D, PISQ, HF, CUTOFF_D, EQ(-6:6), SCALE_D
  common /COLFAC/ CF, CA, TR, PI_D, PISQ, HF, CUTOFF_D, EQ, SCALE_D, SCHEME_D, NF
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
  real(dp), external :: DOT, LEIV, ERTV
  real(dp) :: qnf, q3, nx3, ny3, hf1, mk(2)
  real(dp) :: P(4,7), M(-6:6), h(0:2), hk(0:2,2), w(2), r(5), Q, Q2, S, O, QQ, GQ, ref, wg
  real(dp) :: worst(4)
  integer :: ip, f, i
  call set_parameters()
  if (CC) stop 'harness_born21: CC not implemented'
  ew31_mode = merge(0, 1, noZ)
  ew31_lepton = merge(1, 0, positron) + merge(2, 0, neutrino)
  CF = 4.0_dp/3; CA = 3; TR = 0.5_dp; NF = 5; PI_D = pi; PISQ = pi**2; HF = 0.5_dp
  CUTOFF_D = 1d-8; SCHEME_D = 0; SCALE_D = 1
  EQ = 0
  do i = 1, 5
     EQ(i) = merge(2.0_dp/3, -1.0_dp/3, mod(i, 2) == 0); EQ(-i) = -EQ(i)
  enddo
  call scet_set_colour(CF, CA, TR)
  call random_seed(put=[(2468 + 13*i, i = 1, 33)])
  worst = 0
  do ip = 1, 50
     call random_number(r)
     Q = 10 + 290*r(1); Q2 = Q*Q
     call breit_born(Q, 0.1_dp + 0.8_dp*r(2), 0.05_dp + 0.9_dp*r(3), acos(-1 + 2*r(4)), 2*pi*r(5), P)
     call MATTHR(P, M)
     S = DOT(P,1,6)**2 + DOT(P,2,7)**2; O = DOT(P,1,7)**2 + DOT(P,2,6)**2
     call hard21_eval(P, 1, -1.0_dp/3, h, hk)
     worst(1) = max(worst(1), abs(hk(0,1)/hk(0,2)*O/S - 1))
     QQ = 8*(4*pi/137)**2*(S + O)*16*pi**2*CF/(-4*DOT(P,2,3)*DOT(P,1,3)*DOT(P,5,5))
     GQ = 8*(4*pi/137)**2*(DOT(P,3,6)**2 + DOT(P,3,7)**2 + DOT(P,2,7)**2 + DOT(P,2,6)**2) &
          & *16*pi**2*TR/(-4*DOT(P,2,1)*DOT(P,3,1)*DOT(P,5,5))
     wg = 0
     do f = -5, 5
        if (f == 0) cycle
        call ew31_w(f, Q2, w)
        if (f > 0) wg = wg + (w(1) + w(2))/2
        ref = QQ/(S + O)*(S*w(1) + O*w(2))
        worst(2) = max(worst(2), abs(M(f)/ref - 1))
        if (ip == 1) write(*,'(a,i3,a,2es14.6)') ' f', f, '  MATTHR, sliced21 weights', M(f), ref
     enddo
     worst(3) = max(worst(3), abs(M(0)/(GQ*wg) - 1))
     qnf = -((4*pi/137)**2*4/Q2)*(2*LEIV(P, P(1,6), 2, -1, 3) - Q2/2*ERTV(P, 2, -1, 3))
     call VIRT3PV(P, nx3, ny3)
     q3 = -(4*pi/137)**2*4*CF*(CF*nx3 + CA*ny3)
     hf1 = hard_fact(.false., log(2*DOT(P,1,2)/Q2), log(2*DOT(P,1,3)/Q2), log(2*DOT(P,2,3)/Q2))
     mk = QQ/(S + O)*[S, O]
     worst(4) = max(worst(4), abs(hk(1,1) - hf1 - (qnf + q3)/2/mk(1))/max(1.0_dp, abs(hk(1,1))), &
          & abs(hk(1,2) - hf1 - (qnf - q3)/2/mk(2))/max(1.0_dp, abs(hk(1,2))))
     if (abs(hk(1,1) - hf1 - (qnf + q3)/2/mk(1)) > 1e-8_dp*max(1.0_dp, abs(hk(1,1)))) &
          & write(*,'(a,4es12.4,a,6es14.6)') ' DIAG s12 s13 s23 s45', 2*DOT(P,2,3), -2*DOT(P,1,2), -2*DOT(P,1,3), -Q2, &
          & ' h, S O sum; ref', hk(1,1), hk(1,2), h(1), hf1 + (qnf + q3)/2/mk(1), hf1 + (qnf - q3)/2/mk(2), hf1 + qnf/QQ
     if (ip <= 2) write(*,'(a,4es16.8)') ' H1/H0 classes S, O: hard21, DISENT', hk(1,1), hf1 + (qnf + q3)/2/mk(1), &
          & hk(1,2), hf1 + (qnf - q3)/2/mk(2)
  enddo
  write(*,'(a,4es10.2)') ' largest deviation: classes S:O, quark Born, gluon Born, one loop:', worst
  if (maxval(worst(1:3)) > 1d-10 .or. worst(4) > 1d-8) stop 1
contains
  ! as harness_hard21
  subroutine breit_born(Q, y, xp, th, ph, P)
    real(dp), intent(in) :: Q, y, xp, th, ph
    real(dp), intent(out) :: P(4,7)
    real(dp) :: Ea, W(4), sh, E, k(4), gam, bz
    integer :: n
    P = 0
    Ea = Q/(2*xp)
    P(:,1) = [0.0_dp, 0.0_dp, Ea, Ea]
    P(:,5) = [0.0_dp, 0.0_dp, -Q, 0.0_dp]
    W = P(:,1) + P(:,5)
    sh = W(4)**2 - W(3)**2
    E = sqrt(sh)/2
    k = [E*sin(th)*cos(ph), E*sin(th)*sin(ph), E*cos(th), E]
    bz = W(3)/W(4); gam = 1/sqrt(1 - bz**2)
    P(:,2) = k; P(:,3) = [-k(1), -k(2), -k(3), k(4)]
    do n = 2, 3
       E = P(4,n)
       P(4,n) = gam*(E + bz*P(3,n)); P(3,n) = gam*(P(3,n) + bz*E)
    enddo
    P(:,6) = Q/2*[2*sqrt(1 - y)/y, 0.0_dp, -1.0_dp, (2 - y)/y]
    P(:,7) = P(:,6) - P(:,5)
  end subroutine breit_born
end program harness_born21
