!----------------------------------------------------------------------
! Unit tests of the Born phase-space generator gen_phsp_born and of
! p2bmomenta (src/mod_phase_space.f90): kinematics of the generated
! momenta, consistency between the two frames, and the Jacobian,
! whose integral over the unit square must be the phase-space volume.
program test_phase_space
  use types, only: dp
  use mod_parameters
  use mod_phase_space
  use test_utils
  implicit none

  El = 27.5_dp
  Eh = 820.0_dp
  s  = 4.0_dp * El * Eh

  ! x, y and Q2 ranges, with xmin low enough that the smallest x values
  ! have no phase space (Q2min > x*ymax*s), so that the vanishing
  ! Jacobian is exercised too.
  call set_ranges(1e-4_dp, 0.5_dp, 0.01_dp, 0.95_dp, 10.0_dp, 1e4_dp)
  call check_momenta('x,y,Q2 ranges')
  call check_close('x,y,Q2 ranges: integral of Jacobian = volume in (x,Q2)', &
       & integrate_jacobian(1000, 1000), volume_x_Q2(), 1e-4_dp)

  ! Fixed x: result is differential in x, integrated over Q2
  call set_ranges(0.01_dp, 0.01_dp, 0.01_dp, 0.95_dp, 10.0_dp, 1e4_dp)
  call check_momenta('fixed x')
  call check_close('fixed x: integral of Jacobian = Q2 range', &
       & integrate_jacobian(1, 1000), &
       & min(Q2max, xmin*ymax*s) - max(Q2min, xmin*ymin*s), 1e-6_dp)

  ! Fixed Q2: result is differential in Q2, integrated over x. The
  ! Jacobian then only depends on xborn(1) but has steps where y leaves
  ! its range, so the midpoint rule needs many points.
  call set_ranges(1e-3_dp, 0.5_dp, 0.01_dp, 0.95_dp, 100.0_dp, 100.0_dp)
  call check_momenta('fixed Q2')
  call check_close('fixed Q2: integral of Jacobian = x range', &
       & integrate_jacobian(200000, 1), &
       & min(xmax, Q2min/(ymin*s)) - max(xmin, Q2min/(ymax*s)), 1e-4_dp)

  call finish_tests()

contains

  subroutine set_ranges(x1, x2, y1, y2, q21, q22)
    real(dp), intent(in) :: x1, x2, y1, y2, q21, q22
    xmin = x1; xmax = x2
    ymin = y1; ymax = y2
    Q2min = q21; Q2max = q22
  end subroutine set_ranges

  ! Check the generated momenta at a grid of phase-space points.
  subroutine check_momenta(tag)
    character(len=*), intent(in) :: tag
    real(dp) :: xborn(2), x, y, Qsq, Qvec(0:3), jac
    real(dp) :: plab(0:3,4), pbreit(0:3,4), pb2(0:3,4), pl2(0:3,4)
    real(dp) :: P(0:3), k(0:3), Qval
    integer :: i, j, npoints
    logical :: ok_range, ok_cons, ok_mass, ok_inv, ok_breit, ok_boost, ok_p2b

    P = [Eh, 0.0_dp, 0.0_dp, Eh]   ! proton
    k = [El, 0.0_dp, 0.0_dp, -El]  ! incoming lepton
    ok_range = .true.; ok_cons = .true.; ok_mass = .true.; ok_inv = .true.
    ok_breit = .true.; ok_boost = .true.; ok_p2b = .true.
    npoints = 0
    do i = 1, 15
       do j = 1, 15
          xborn = [(i - 0.5_dp)/15, (j - 0.5_dp)/15]
          call gen_phsp_born(xborn, x, y, Qsq, Qvec, jac, plab, pbreit)
          if (jac == 0.0_dp) cycle
          npoints = npoints + 1
          Qval = sqrt(Qsq)

          ok_range = ok_range .and. x >= xmin*(1-1e-12_dp) .and. x <= xmax*(1+1e-12_dp) &
               & .and. y >= ymin*(1-1e-12_dp) .and. y <= ymax*(1+1e-12_dp) &
               & .and. Qsq >= Q2min*(1-1e-12_dp) .and. Qsq <= Q2max*(1+1e-12_dp)
          ! lab frame: beams, momentum conservation, massless particles,
          ! and x, y, Q2 recovered from the momenta
          ok_cons = ok_cons .and. close_vec(plab(:,1), k, El) &
               & .and. close_vec(plab(:,2), x*P, Eh) &
               & .and. close_vec(plab(:,1) + plab(:,2), plab(:,3) + plab(:,4), Eh) &
               & .and. close_vec(Qvec, plab(:,1) - plab(:,3), Eh)
          ok_mass = ok_mass .and. abs(dot(plab(:,3),plab(:,3))) < 1e-9_dp*Eh**2 &
               & .and. abs(dot(plab(:,4),plab(:,4))) < 1e-9_dp*Eh**2
          ok_inv = ok_inv .and. abs(-dot(Qvec,Qvec)/Qsq - 1) < 1e-9_dp &
               & .and. abs(Qsq/(2*dot(P,Qvec))/x - 1) < 1e-9_dp &
               & .and. abs(dot(P,Qvec)/dot(P,k)/y - 1) < 1e-9_dp &
               & .and. abs(Qsq/(x*y*s) - 1) < 1e-12_dp
          ! Breit frame: q = (0,0,0,-Q), incoming parton (Q/2)(1,0,0,1)
          ok_breit = ok_breit &
               & .and. close_vec(pbreit(:,1) - pbreit(:,3), [0.0_dp, 0.0_dp, 0.0_dp, -Qval], Qval) &
               & .and. close_vec(pbreit(:,2), 0.5_dp*Qval*[1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], Qval) &
               & .and. close_vec(pbreit(:,1) + pbreit(:,2), pbreit(:,3) + pbreit(:,4), Qval) &
               & .and. abs(dot(pbreit(:,1),pbreit(:,1))) < 1e-9_dp*Qsq/y**2 &
               & .and. abs(dot(pbreit(:,3),pbreit(:,3))) < 1e-9_dp*Qsq/y**2
          ! the hand-built Breit momenta are the boosted lab momenta
          call mlab2breit(4, Qvec, plab, pb2, .true.)
          ok_boost = ok_boost .and. close_vec(reshape(pb2,[16]), reshape(pbreit,[16]), Qval/y)
          ! p2bmomenta builds the same momenta from (x,y,Q2)
          call p2bmomenta(x, y, Qsq, pb2, pl2)
          ok_p2b = ok_p2b .and. all(pb2 == pbreit) .and. all(pl2 == plab)
       enddo
    enddo
    call check_true(tag//': found phase-space points', npoints > 20)
    call check_true(tag//': x, y, Q2 inside the requested ranges', ok_range)
    call check_true(tag//': lab frame beams and momentum conservation', ok_cons)
    call check_true(tag//': lab frame outgoing particles massless', ok_mass)
    call check_true(tag//': Q2, x, y recovered from lab momenta', ok_inv)
    call check_true(tag//': Breit frame momenta', ok_breit)
    call check_true(tag//': mlab2breit(plab) = pbreit', ok_boost)
    call check_true(tag//': p2bmomenta agrees with gen_phsp_born', ok_p2b)
  end subroutine check_momenta

  ! Midpoint rule for the integral of the Jacobian over [0,1]^2
  real(dp) function integrate_jacobian(n1, n2)
    integer, intent(in) :: n1, n2
    real(dp) :: xborn(2), x, y, Qsq, Qvec(0:3), jac
    real(dp) :: plab(0:3,4), pbreit(0:3,4)
    integer :: i, j
    integrate_jacobian = 0.0_dp
    do i = 1, n1
       do j = 1, n2
          xborn = [(i - 0.5_dp)/n1, (j - 0.5_dp)/n2]
          call gen_phsp_born(xborn, x, y, Qsq, Qvec, jac, plab, pbreit)
          integrate_jacobian = integrate_jacobian + jac
       enddo
    enddo
    integrate_jacobian = integrate_jacobian / (real(n1,dp) * n2)
  end function integrate_jacobian

  ! Area of the allowed region in the (x,Q2) plane,
  ! xmin<x<xmax, Q2min<Q2<Q2max, ymin<Q2/(x s)<ymax. The integrand
  ! in x is piecewise linear, so the trapezoidal rule on a fine grid
  ! is accurate far beyond the tolerance used.
  real(dp) function volume_x_Q2()
    integer, parameter :: n = 200000
    real(dp) :: xa, xb, fa, fb
    integer :: i
    volume_x_Q2 = 0.0_dp
    xa = xmin
    fa = width(xa)
    do i = 1, n
       xb = xmin + (xmax - xmin) * i / n
       fb = width(xb)
       volume_x_Q2 = volume_x_Q2 + 0.5_dp * (fa + fb) * (xb - xa)
       xa = xb; fa = fb
    enddo
  end function volume_x_Q2

  real(dp) function width(x)
    real(dp), intent(in) :: x
    width = max(0.0_dp, min(Q2max, x*ymax*s) - max(Q2min, x*ymin*s))
  end function width

  pure real(dp) function dot(a, b)
    real(dp), intent(in) :: a(0:3), b(0:3)
    dot = a(0)*b(0) - a(1)*b(1) - a(2)*b(2) - a(3)*b(3)
  end function dot

  pure logical function close_vec(a, b, scale)
    real(dp), intent(in) :: a(:), b(:), scale
    close_vec = all(abs(a - b) <= 1e-10_dp * scale)
  end function close_vec

end program test_phase_space
