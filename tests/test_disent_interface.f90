!----------------------------------------------------------------------
! Unit test of the interface to DISENT (src/mod_disent_interface.f90).
!
! DISENTFULL is run for a few thousand events up to O(alphas^2) with a
! checking callback instead of the real `user` routine. For every
! Born, real and double-real configuration DISENT hands over, the
! callback applies the same helpers `user` uses (disent_kinematics,
! disent_to_momenta, p2bmomenta, mbreit2lab) and checks that
!   - DISENT's momenta conserve momentum, are massless, lie inside
!     the cuts returned by dis_cuts,
!   - x, y, Q2 are reconstructed consistently,
!   - after the transformation to the lab frame the beams are the
!     physical ones (lepton (El,0,0,-El), parton eta*(Eh,0,0,Eh)),
!     momentum is conserved, and the outgoing lepton is the same as
!     in the projected Born event, as it must be in DIS.
! It also checks dis_cuts, disent_muf and the mapping of DISENT's
! three muF points onto our 7-point scale variation.
module disent_checks
  use types, only: dp
  use mod_parameters
  use mod_phase_space
  use mod_disent_interface
  implicit none
  integer, save :: ncalls(2:4) = 0
  logical, save :: ok_cons = .true., ok_mass = .true., ok_cuts = .true.
  logical, save :: ok_kin = .true., ok_lab = .true., ok_lepton = .true.
  real(dp), save :: cut_xmin, cut_xmax, cut_Q2min, cut_Q2max, cut_ymin, cut_ymax

contains

  subroutine checking_user(n, na, nt, p, sdis, weight, scale2)
    integer, intent(in)  :: n, na, nt
    real(dp), intent(in) :: sdis, p(4,7), weight(-6:6), scale2
    real(dp) :: eta, x, y, Q2, pbreit(0:3,6), plab(0:3,6), mom(4)
    real(dp) :: p2bbreit(0:3,4), p2blab(0:3,4), Qlab(0:3), escale
    integer :: i

    if (n == 0) return
    ncalls(n) = ncalls(n) + 1
    escale = p(4,6) + p(4,1)

    ! DISENT momenta: q = k - k' = p5, conservation, massless partons
    mom = p(:,1) + p(:,6) - p(:,7)
    do i = 2, n
       mom = mom - p(:,i)
    enddo
    ok_cons = ok_cons .and. all(abs(mom) < 1e-9_dp * escale) &
         & .and. all(abs(p(:,6) - p(:,7) - p(:,5)) < 1e-9_dp * escale)
    do i = 1, n
       ok_mass = ok_mass .and. abs(DOT(p,i,i)) < 1e-8_dp * escale**2
    enddo
    ok_mass = ok_mass .and. abs(DOT(p,6,6)) < 1e-8_dp * escale**2 &
         & .and. abs(DOT(p,7,7)) < 1e-8_dp * escale**2

    ! Born variables, as `user` computes them
    eta = 2 * DOT(p,1,6) / sdis
    call disent_kinematics(p, eta, x, y, Q2)
    ok_kin = ok_kin .and. abs(Q2 / (x * y * sdis) - 1) < 1e-9_dp &
         & .and. eta >= x * (1 - 1e-12_dp) .and. eta <= 1
    ok_cuts = ok_cuts .and. inside(x, cut_xmin, cut_xmax) &
         & .and. inside(Q2, cut_Q2min, cut_Q2max) .and. inside(y, cut_ymin, cut_ymax)

    ! To the lab frame, as `user` does it
    call p2bmomenta(x, y, Q2, p2bbreit, p2blab)
    Qlab(:) = p2blab(:,1) - p2blab(:,3)
    call disent_to_momenta(n, p, pbreit(:,1:n+2))
    call mbreit2lab(n+2, Qlab, pbreit(:,1:n+2), plab(:,1:n+2), .true.)
    ok_lab = ok_lab .and. close(plab(:,1), [El, 0.0_dp, 0.0_dp, -El], Eh) &
         & .and. close(plab(:,2), eta * [Eh, 0.0_dp, 0.0_dp, Eh], Eh) &
         & .and. close(plab(:,1) + plab(:,2), sum(plab(:,3:n+2), dim=2), Eh)
    ok_lepton = ok_lepton .and. close(plab(:,3), p2blab(:,3), Eh)
  end subroutine checking_user

  logical function inside(v, vmin, vmax)
    real(dp), intent(in) :: v, vmin, vmax
    inside = v >= vmin * (1 - 1e-9_dp) .and. v <= vmax * (1 + 1e-9_dp)
  end function inside

  logical function close(a, b, scale)
    real(dp), intent(in) :: a(0:3), b(0:3), scale
    close = all(abs(a - b) <= 1e-9_dp * scale)
  end function close

end module disent_checks

program test_disent_interface
  use types, only: dp
  use mod_parameters
  use mod_disent_interface
  use disent_checks
  use test_utils
  implicit none
  real(dp) :: xl, xu, Q2l, Q2u, yl, yu, muF2, p(4,7)
  integer :: isc
  real(dp), parameter :: disent_muf_factors(3) = [1.0_dp, 2.0_dp, 0.5_dp]

  ! Module state that set_parameters would otherwise provide
  El = 27.5_dp
  Eh = 820.0_dp
  s  = 4 * El * Eh
  Qmin = 1.0_dp
  scale_choice = 2
  xmur = 1.0_dp
  xmuf = 1.0_dp

  !---- dis_cuts
  xmin = 0.01_dp; xmax = 0.01_dp; Q2min = 100.0_dp; Q2max = 100.0_dp
  ymin = 0.3_dp; ymax = 0.3_dp
  call dis_cuts(s, xl, xu, Q2l, Q2u, yl, yu)
  call check_true('dis_cuts: y left free when x and Q2 are fixed', &
       & yl == 0 .and. yu == 1 .and. xl == xmin .and. Q2l == Q2min)
  xmin = 1e-3_dp; xmax = 0.5_dp; Q2min = 10.0_dp; Q2max = 1e4_dp
  ymin = 0.05_dp; ymax = 0.9_dp
  call dis_cuts(s, xl, xu, Q2l, Q2u, yl, yu)
  call check_true('dis_cuts: ranges passed through', &
       & xl == xmin .and. xu == xmax .and. Q2l == Q2min .and. Q2u == Q2max &
       & .and. yl == ymin .and. yu == ymax)

  !---- disent_muf returns (muF/Q)^2 for a Born configuration
  p = 0
  p(:,1) = [0.0_dp, 0.0_dp, 5.0_dp, 5.0_dp]
  p(:,5) = [0.0_dp, 0.0_dp, -10.0_dp, 0.0_dp]
  p(:,6) = [5.0_dp/0.4_dp*2*sqrt(0.6_dp), 0.0_dp, -5.0_dp, 5.0_dp/0.4_dp*1.6_dp]
  xmuf = 2.0_dp
  call disent_muf(p, s, muF2)
  call check_close('disent_muf: (muF/Q)^2 = xmuf^2 for muF = xmuf Q', muF2, 4.0_dp, 1e-12_dp)
  xmuf = 1.0_dp

  !---- mapping of DISENT's three muF points onto the 7-point variation
  do isc = 1, maxscales
     call check_close('DISENT muF point for scale variation point', &
          & disent_muf_factors(disent_imuf(isc)), scales_muf(isc), 1e-14_dp)
     call check_close('DISENT muR point for scale variation point', &
          & scales_mur(disent_imur(isc)), scales_mur(isc), 1e-14_dp)
  enddo

  !---- DISENT events, with x, y and Q2 ranges
  call dis_cuts(s, cut_xmin, cut_xmax, cut_Q2min, cut_Q2max, cut_ymin, cut_ymax)
  call DISENTFULL(2000, s, 5, checking_user, dis_cuts, 12345, 67890, &
       & 2.0_dp, 4.0_dp, 1e-8_dp, 2, disent_muf, 4.0_dp/3.0_dp, 3.0_dp, &
       & 0.5_dp, .false.)
  call report('x, y, Q2 ranges')

  !---- and with fixed x and Q2, as in the validation runs
  xmin = 0.01_dp; xmax = 0.01_dp; Q2min = 100.0_dp; Q2max = 100.0_dp
  ymin = Q2min / (s * xmin); ymax = ymin
  call dis_cuts(s, cut_xmin, cut_xmax, cut_Q2min, cut_Q2max, cut_ymin, cut_ymax)
  cut_ymin = ymin; cut_ymax = ymax
  ncalls = 0
  call DISENTFULL(2000, s, 5, checking_user, dis_cuts, 12345, 67890, &
       & 2.0_dp, 4.0_dp, 1e-8_dp, 2, disent_muf, 4.0_dp/3.0_dp, 3.0_dp, &
       & 0.5_dp, .false.)
  call report('fixed x and Q2')

  call finish_tests()

contains

  subroutine report(tag)
    character(len=*), intent(in) :: tag
    ! DISENTFULL (as modified for disorder) hands over real and
    ! double-real configurations; the Born ones come from the
    ! structure functions.
    call check_true(tag//': DISENT produced real and double-real events', &
         & ncalls(3) > 1000 .and. ncalls(4) > 1000)
    call check_true(tag//': DISENT momenta conserve momentum', ok_cons)
    call check_true(tag//': DISENT momenta massless', ok_mass)
    call check_true(tag//': x, Q2, y inside the cuts', ok_cuts)
    call check_true(tag//': Q2 = x y s and x <= eta <= 1', ok_kin)
    call check_true(tag//': lab frame beams and momentum conservation', ok_lab)
    call check_true(tag//': lab frame lepton same as in projected Born', ok_lepton)
  end subroutine report

end program test_disent_interface
