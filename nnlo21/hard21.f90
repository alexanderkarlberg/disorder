!-----------------------------------------------------------------------
! Hard function of DIS 2+1 (photon exchange) at one and two loops, in the
! SCET (MS-bar, minimal subtraction of the infrared poles) scheme, for the
! tau_2 slicing of NNLO DIS 2+1.
!
! Amplitudes: gamma* -> q qbar g helicity coefficients alpha, beta, gamma
! (Garland et al., hep-ph/0112081, 0206067) continued to DIS kinematics
! (Gehrmann, Remiddi; Gehrmann, Glover arXiv:0904.2665), as functions of the
! 16 kinematic regions in NNLOJET v1.0.2 (makecoef, makecoeftaylorU in
! src/process/Z/B1gNZ.f; GPL-3.0-or-later, linked from libnnlojet_core).
! The region selection below follows NNLOJET's helcoeff. NNLOJET's region
! coefficients equal the Catani-scheme finite remainders of 0904.2665
! (checked in all eight regions of its Fortran files, all 30 coefficients,
! to 1e-6 or better) up to a convention in the one-loop N and beta_0
! coefficients of alpha and beta:
!     a_GG = a_NJ + (11/24)(L13 + L23),   c_GG = c_NJ - (1/3)(L13 + L23),
! with L_ij = ln(-s_ij/mu^2 - i0), mu^2 = -s_123 = Q^2.
!
! Colour decomposition (0904.2665, eqs. (3.23), (3.26)), alpha_s/2pi:
!   Omega1 = N a + b/N + beta0 c,   beta0 = (11 N - 2 n_f)/6,
!   Omega2 = N^2 A + B + C/N^2 + N n_f D + n_f/N E + n_f^2 F + N_Fgam (4/N - N) G,
! N_Fgam = sum_q e_q / e_q (photon coupling to a closed quark loop).
! Catani -> SCET at mu^2 = Q^2 (derived in docs/nnlo21-plan.md; equal to
! MCFM 10.3's schemeconvC0, schemeconv2lM0 in the timelike region):
!   C1 = Omega1 + I1 Omega0,   C2 = Omega2 + I1 Omega1 + X2 Omega0,
! I1 = eps^0 coefficient of Catani's I^(1)(eps), X2 = I1^2/2 + R^(0).
!
! Spinor structures: MCFM's ampqqbgll (Eq. 22 of arXiv:1309.3245), all
! momenta outgoing: q(1) qbar(2) g(3) l-(4) l+(5); the other gluon
! helicity from the coefficients with q <-> qbar (as MCFM's Zampqqbgsq).
!
! hard21_eval(P, ichan, h): P(4,7) in DISENT's layout (1 incoming parton,
! 2, 3 outgoing, 6 incoming lepton, 7 outgoing lepton; px, py, pz, E);
! ichan = 1: incoming quark (or antiquark, same for photon exchange),
! parton 2 the outgoing quark, 3 the gluon; ichan = 2: incoming gluon,
! parton 2 the quark, 3 the antiquark. h(0) = sum over helicities of
! |M0|^2 (spinor level, no couplings), h(1), h(2) = H^(1)/H^(0) and
! H^(2)/H^(0) in powers of alpha_s(mu)/2pi at mu^2 = Q^2. scheme = 0: SCET,
! 1: Catani's finite remainders (for tests).
!-----------------------------------------------------------------------
module hard21
  implicit none
  private
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: mxpart = 14
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: zeta3 = 1.2020569031595942853997381615114499907650_dp
  real(dp), parameter :: xn = 3
  real(dp), save, public :: hard21_nf = 5
  ! N_Fgam (photon on a closed quark loop): sum_q e_q/e_q of the Born quark
  ! for the five light flavours (sum e_q = 1/3); 0 switches it off
  logical, save, public :: hard21_nfz = .true.
  integer, save, public :: hard21_scheme = 0
  public :: hard21_eval, hard21_coeffs
contains

  subroutine hard21_eval(P, ichan, eq, h)
    real(dp), intent(in) :: P(4,7), eq
    integer, intent(in) :: ichan
    real(dp), intent(out) :: h(0:2)
    complex(dp) :: za(mxpart,mxpart), zb(mxpart,mxpart)
    real(dp) :: pm(mxpart,4), s12, s13, s23, s45
    complex(dp) :: c(3,0:2,2), amp(0:2,8), al, be, ga, de
    complex(dp), external :: ampqqbgll
    integer :: k, j, iperm
    pm = 0
    if (ichan == 1) then
       pm(1,:) = P(:,2); pm(2,:) = -P(:,1); pm(3,:) = P(:,3)
    else
       pm(1,:) = P(:,2); pm(2,:) = P(:,3); pm(3,:) = -P(:,1)
    endif
    pm(4,:) = P(:,7); pm(5,:) = -P(:,6)
    call spinoru(5, pm, za, zb)
    s12 = mdot(pm(1,:) + pm(2,:)); s13 = mdot(pm(1,:) + pm(3,:))
    s23 = mdot(pm(2,:) + pm(3,:)); s45 = mdot(pm(4,:) + pm(5,:))
    call hard21_coeffs(s12, s13, s23, s45, eq, c)
    amp = 0
    do iperm = 1, 2
       do k = 0, 2
          al = c(1,k,iperm); be = c(2,k,iperm); ga = c(3,k,iperm)
          de = s12 * (al - be - ga) / (2 * (s12 + s13 + s23))
          if (iperm == 1) then
             amp(k,1) = ampqqbgll(1,2,3,5,4,al,be,ga,de,zb,za)
             amp(k,2) = ampqqbgll(1,2,3,5,4,al,be,ga,de,za,zb)
             amp(k,3) = ampqqbgll(1,2,3,4,5,al,be,ga,de,zb,za)
             amp(k,4) = ampqqbgll(1,2,3,4,5,al,be,ga,de,za,zb)
          else
             amp(k,5) = -ampqqbgll(2,1,3,4,5,al,be,ga,de,za,zb)
             amp(k,6) = -ampqqbgll(2,1,3,4,5,al,be,ga,de,zb,za)
             amp(k,7) = -ampqqbgll(2,1,3,5,4,al,be,ga,de,za,zb)
             amp(k,8) = -ampqqbgll(2,1,3,5,4,al,be,ga,de,zb,za)
          endif
       enddo
    enddo
    h = 0
    do j = 1, 8
       h(0) = h(0) + abs(amp(0,j))**2
       h(1) = h(1) + 2 * real(amp(0,j) * conjg(amp(1,j)), dp)
       h(2) = h(2) + 2 * real(amp(0,j) * conjg(amp(2,j)), dp) + abs(amp(1,j))**2
    enddo
    h(1:2) = h(1:2) / h(0)
  contains
    real(dp) function mdot(v)
      real(dp), intent(in) :: v(4)
      mdot = v(4)**2 - v(1)**2 - v(2)**2 - v(3)**2
    end function mdot
  end subroutine hard21_eval

  ! helicity coefficients c(alpha|beta|gamma, order, 1 = direct | 2 = q <->
  ! qbar) at mu^2 = -s45, in the scheme hard21_scheme; s12 = s(q,qbar),
  ! s13 = s(q,g), s23 = s(g,qbar), s45 = q^2 < 0; eq the Born quark charge
  subroutine hard21_coeffs(s12, s13, s23, s45, eq, c)
    real(dp), intent(in) :: s12, s13, s23, s45, eq
    complex(dp), intent(out) :: c(3,0:2,2)
    complex(dp) :: za(30), o1(3), o2(3), L12, L13, L23, S1, I10, X2
    real(dp) :: u(2), v(2), omv, deldis, adis, omumv, nf, b0, nfz
    integer :: kinregion(2), iregion, i
    nf = hard21_nf
    b0 = (11 * xn - 2 * nf) / 6
    nfz = 0
    if (hard21_nfz) nfz = (1.0_dp/3) / eq
    L12 = lgm(s12); L13 = lgm(s13); L23 = lgm(s23)
    ! kinematic region (NNLOJET helcoeff, s45 < 0)
    iregion = 0
    if (s12 > 0 .and. s13 > s45 .and. s23 > s45) iregion = 5   ! 1b
    if (s12 > 0 .and. s13 <= s45 .and. s23 <= s45) iregion = 14 ! 2d
    if (s12 > 0 .and. s13 <= s45 .and. s23 > s45) iregion = 11  ! 3c
    if (s12 > 0 .and. s13 > s45 .and. s23 <= s45) iregion = 12  ! 4c
    if (s13 > 0 .and. s12 > s45 .and. s23 > s45) iregion = 9    ! 1c
    if (s13 > 0 .and. s12 <= s45 .and. s23 > s45) iregion = 6   ! 2b
    if (s13 > 0 .and. s12 <= s45 .and. s23 <= s45) iregion = 15 ! 3d
    if (s13 > 0 .and. s12 > s45 .and. s23 <= s45) iregion = 8   ! 4b
    if (s23 > 0 .and. s12 > s45 .and. s13 > s45) iregion = 13   ! 1d
    if (s23 > 0 .and. s12 <= s45 .and. s13 > s45) iregion = 10  ! 2c
    if (s23 > 0 .and. s12 > s45 .and. s13 <= s45) iregion = 7   ! 3b
    if (s23 > 0 .and. s12 <= s45 .and. s13 <= s45) iregion = 16 ! 4d
    select case (iregion)
    case (5);  kinregion = [5, 5];   u = [(s12+s23)/s45, (s12+s13)/s45]; v = -s12/s45
    case (6);  kinregion = [6, 10];  u = -(s13+s23)/s12; v = s23/s12
    case (7);  kinregion = [7, 8];   u = (s13+s23)/s13;  v = s12/s13
    case (8);  kinregion = [8, 7];   u = (s13+s23)/s23;  v = s12/s23
    case (9);  kinregion = [9, 13];  u = (s12+s13)/s45;  v = -s13/s45
    case (10); kinregion = [10, 6];  u = -(s13+s23)/s12; v = s13/s12
    case (11); kinregion = [11, 12]; u = -(s12+s23)/s13; v = s23/s13
    case (12); kinregion = [12, 11]; u = -(s12+s13)/s23; v = s13/s23
    case (13); kinregion = [13, 9];  u = (s12+s23)/s45;  v = -s23/s45
    case (14); kinregion = [14, 14]; u = [(s12+s13)/s12, (s12+s23)/s12]; v = -s45/s12
    case (15); kinregion = [15, 16]; u = (s13+s23)/s13;  v = -s45/s13
    case (16); kinregion = [16, 15]; u = (s13+s23)/s23;  v = -s45/s23
    case default
       c = 0; return
    end select
    S1 = L13 + L23
    I10 = i1zero(L12, L13, L23, nf)
    X2 = x2zero(L12, L13, L23, nf)
    do i = 1, 2
       ! NNLOJET: displacement at v -> 1, Taylor series for u -> 0, 1-u-v -> 0
       omv = 1 - v(i)
       if (omv < 1e-3_dp) then
          v(i) = v(i) - 1e-3_dp; omv = 1 - v(i)
       endif
       deldis = 1e-3_dp * omv
       adis = min(1e-1_dp * omv, 1e-2_dp)
       omumv = 1 - u(i) - v(i)
       if (kinregion(i) >= 5 .and. u(i) < deldis .and. u(i) <= omumv) then
          call makecoeftaylorU(kinregion(i), u(i), v(i), adis, -2 * adis, za)
       elseif (kinregion(i) >= 5 .and. omumv < deldis) then
          call makecoeftaylorU(kinregion(i), u(i), v(i), adis, 2 * adis, za)
       else
          call makecoef(kinregion(i), u(i), v(i), za)
       endif
       ! NNLOJET -> 0904.2665 one-loop convention
       za(1) = za(1) + 11.0_dp/24 * S1; za(2) = za(2) + 11.0_dp/24 * S1
       za(7) = za(7) - S1/3;            za(8) = za(8) - S1/3
       o1 = xn * za(1:3) + za(4:6) / xn + b0 * za(7:9)
       o2 = xn**2 * za(10:12) + za(13:15) + za(16:18) / xn**2 + xn * nf * za(19:21) &
            + nf / xn * za(22:24) + nf**2 * za(25:27) + nfz * (4 / xn - xn) * za(28:30)
       c(:,0,i) = [(1.0_dp, 0.0_dp), (1.0_dp, 0.0_dp), (0.0_dp, 0.0_dp)]
       if (hard21_scheme == 1) then
          c(:,1,i) = o1; c(:,2,i) = o2
       else
          c(:,1,i) = o1 + I10 * c(:,0,i)
          c(:,2,i) = o2 + I10 * o1 + X2 * c(:,0,i)
       endif
    enddo
  contains
    ! L = ln(-s/mu^2 - i0), mu^2 = -s45
    complex(dp) function lgm(s)
      real(dp), intent(in) :: s
      lgm = cmplx(log(abs(s / s45)), merge(-pi, 0.0_dp, s > 0), dp)
    end function lgm
  end subroutine hard21_coeffs

  ! eps^0 coefficient of Catani's I^(1)(eps) (0904.2665 (3.19), N = 3, T_R = 1/2)
  complex(dp) function i1zero(L12, L13, L23, N_F) result(r)
    complex(dp), intent(in) :: L12, L13, L23
    real(dp), intent(in) :: N_F
    r = 0.083333333333333333333d0*L12**2 - 0.25d0*L12 - 0.75d0*L13**2 - &
      0.083333333333333333333d0*L13*N_F + 2.5d0*L13 - 0.75d0*L23**2 - &
      0.083333333333333333333d0*L23*N_F + 2.5d0*L23 + &
      2.330323261368320785d0
  end function i1zero

  ! X2 = I1^2/2 + R^(0) (see the header), N = 3, T_R = 1/2
  complex(dp) function x2zero(L12, L13, L23, N_F) result(r)
    complex(dp), intent(in) :: L12, L13, L23
    real(dp), intent(in) :: N_F
    include 'x2zero.inc'
  end function x2zero
end module hard21
