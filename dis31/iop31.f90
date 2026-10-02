!-----------------------------------------------------------------------
! Integrated Catani-Seymour dipoles for DIS 3+1 (one incoming parton;
! Catani, Seymour, hep-ph/9605323, section 8, MS-bar, K_F.S. = 0, n_f = 5),
! with the colour correlations of born31. Normalisation as born31/me31;
! all operators in units of alpha_s/2pi, with (4 pi)^eps/Gamma(1 - eps)
! factored out (as virt31_ren).
!
! iop31_i(P, fl, mu2, iv): <I(eps)> (CS 8.25 with 7.27, 7.28), Laurent
!   coefficients iv(-2:0), for the Born P(4,7), fl(4) (slot 1 the incoming
!   parton with its Born momentum), renormalisation scale mu2:
!     I = - sum_I 1/T_I^2 V_I(eps) sum_{J /= I} <T_I.T_J> (mu^2/(2 p_I.p_J))^eps,
!     V_I = T_I^2 (1/eps^2 - pi^2/3) + gamma_I/eps + gamma_I + K_I.
!   2 Re<M0|M1> (virt31_ren) + <I> is finite.
! iop31_kp(P, fl, muf2, b, g, lp): the Born-point factors of K + P for the
!   Born's incoming parton b (slot 1):
!     b  = |M_b|^2,
!     g  = sum_{i /= b} gamma_i/T_i^2 <T_i.T_b>,
!     lp = 1/T_b^2 sum_{i /= b} <T_i.T_b> ln(muf2/(2 p_b.p_i)),
!   so that, for an incoming parton a and momentum fraction x (CS 8.38, 8.39),
!     K^{a,b}(x) + P^{a,b}(x) = Kbar^{ab}(x) b
!        + delta^{ab} [(1/(1-x))_+ + delta(1-x)] g + P^{ab}(x) lp,
!   with Kbar^{ab} and P^{ab} from iop31_kernel.
! iop31_kernel(a, b, x, reg, plus, delta): Kbar^{ab}(x) and P^{ab}(x) split
!   into a regular part, the coefficient function of a plus distribution
!   [plus(x)]_+, and a delta(1-x) coefficient, for parton types a, b
!   (0 gluon, 1 quark); kind = 1: Kbar, kind = 2: P.
!-----------------------------------------------------------------------
module iop31
  use born31
  implicit none
  private
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: CF = 4.0_dp/3, CA = 3, TR = 0.5_dp, nf = 5
  real(dp), parameter :: gq = 1.5_dp*CF, gg = 11*CA/6 - 2*TR*nf/3
  real(dp), parameter :: Kq = (3.5_dp - pi**2/6)*CF, Kg = (67.0_dp/18 - pi**2/6)*CA - 10*TR*nf/9
  public :: iop31_i, iop31_kp, iop31_kernel
contains

  subroutine iop31_i(P, fl, mu2, iv)
    real(dp), intent(in) :: P(4,7), mu2
    integer, intent(in) :: fl(4)
    real(dp), intent(out) :: iv(-2:0)
    real(dp) :: msq, cc(4,4), CI, gam, KK, L
    integer :: i, j
    call born31_cc(P, fl, msq, cc)
    iv = 0
    do i = 1, 4
       if (fl(i) == 0) then
          CI = CA; gam = gg; KK = Kg
       else
          CI = CF; gam = gq; KK = Kq
       endif
       do j = 1, 4
          if (j == i) cycle
          L = log(mu2/(2*abs(mdot(P(:,i), P(:,j)))))
          ! - <T_i.T_j>/C_i [C_i (1/eps^2 - pi^2/3) + gam/eps + gam + K]
          !   (1 + eps L + eps^2 L^2/2)
          iv(-2) = iv(-2) - cc(i,j)
          iv(-1) = iv(-1) - cc(i,j)*(L + gam/CI)
          iv(0) = iv(0) - cc(i,j)*(L**2/2 - pi**2/3 + (gam*L + gam + KK)/CI)
       enddo
    enddo
  end subroutine iop31_i

  subroutine iop31_kp(P, fl, muf2, b, g, lp)
    real(dp), intent(in) :: P(4,7), muf2
    integer, intent(in) :: fl(4)
    real(dp), intent(out) :: b, g, lp
    real(dp) :: cc(4,4), Cb
    integer :: i
    call born31_cc(P, fl, b, cc)
    Cb = merge(CA, CF, fl(1) == 0)
    g = 0; lp = 0
    do i = 2, 4
       g = g + merge(gg/CA, gq/CF, fl(i) == 0)*cc(i,1)
       lp = lp + cc(i,1)/Cb*log(muf2/(2*abs(mdot(P(:,1), P(:,i)))))
    enddo
  end subroutine iop31_kp

  ! Kbar^{ab} (kind 1) and P^{ab} (kind 2): a = incoming parton type, b =
  ! the Born's incoming parton type (0 gluon, 1 quark/antiquark of the same
  ! flavour for a = b); the value at x is reg + plus (to be used as
  ! [plus]_+, i.e. integrated against h(x) - h(1)) + delta delta(1-x)
  subroutine iop31_kernel(kind, a, b, x, reg, plus, delta)
    integer, intent(in) :: kind, a, b
    real(dp), intent(in) :: x
    real(dp), intent(out) :: reg, plus, delta
    real(dp) :: l
    reg = 0; plus = 0; delta = 0
    l = log((1 - x)/x)
    if (kind == 1) then
       if (a == 1 .and. b == 0) then
          ! Kbar^{qg} = P^{qg} ln((1-x)/x) + CF x
          reg = CF*(1 + (1 - x)**2)/x*l + CF*x
       elseif (a == 0 .and. b == 1) then
          ! Kbar^{gq} = P^{gq} ln((1-x)/x) + TR 2x(1-x)
          reg = TR*(x**2 + (1 - x)**2)*l + TR*2*x*(1 - x)
       elseif (a == 1 .and. b == 1) then
          ! CF [2/(1-x) ln((1-x)/x)]_+ - CF (1+x) ln((1-x)/x) + CF (1-x)
          ! - delta(1-x) (5 - pi^2) CF
          plus = CF*2/(1 - x)*l
          reg = -CF*(1 + x)*l + CF*(1 - x)
          delta = -(5 - pi**2)*CF
       else
          ! 2CA [1/(1-x) ln((1-x)/x)]_+ + 2CA ((1-x)/x - 1 + x(1-x)) ln((1-x)/x)
          ! - delta(1-x) ((50/9 - pi^2) CA - 16/9 TR nf)
          plus = 2*CA/(1 - x)*l
          reg = 2*CA*((1 - x)/x - 1 + x*(1 - x))*l
          delta = -((50.0_dp/9 - pi**2)*CA - 16.0_dp/9*TR*nf)
       endif
    else
       if (a == 1 .and. b == 0) then
          reg = CF*(1 + (1 - x)**2)/x
       elseif (a == 0 .and. b == 1) then
          reg = TR*(x**2 + (1 - x)**2)
       elseif (a == 1 .and. b == 1) then
          ! CF ((1 + x^2)/(1 - x))_+
          plus = CF*(1 + x**2)/(1 - x)
       else
          ! 2CA [1/(1-x)]_+ + 2CA ((1-x)/x - 1 + x(1-x)) + delta(1-x) gamma_g
          plus = 2*CA/(1 - x)
          reg = 2*CA*((1 - x)/x - 1 + x*(1 - x))
          delta = gg
       endif
    endif
  end subroutine iop31_kernel

  pure real(dp) function mdot(a, b)
    real(dp), intent(in) :: a(4), b(4)
    mdot = a(4)*b(4) - a(1)*b(1) - a(2)*b(2) - a(3)*b(3)
  end function mdot
end module iop31
