!----------------------------------------------------------------------
! Component checks of mod_slicing_scet:
!  1. cli2 against mpmath values;
!  2. the soft non-hemisphere integrals I0, I1 against mpmath (30 digits,
!     the same 1D representation, split at the kinks; a plain scipy 2D
!     quadrature was off by 3e-5 and fails where the log singularity at
!     y = 1, phi = 0 lies inside the region, as for (2, 0.3));
!  3. the complete one-loop 1-jettiness cumulant for a 1+1 Born (beam,
!     jet, two-direction soft and hard functions) against the validated
!     Python implementation of KLS (173) (tau_1^a; slicing/tau1b_nlo.py
!     --cumulant-only) at x = 0.01, Q = 20, y = 0.5, NNPDF30_nlo_as_0118.
!----------------------------------------------------------------------
program test_scet
  use types, only: dp
  use mod_slicing_scet
  implicit none
  real(dp), external :: alphasPDF
  complex(dp) :: z(6), ref(6)
  real(dp) :: I0, I1, x, Q, y, a, Yp, tau, f1(-6:6), c0(-6:6), c1(-6:6), c2(-6:6)
  real(dp) :: nhat(3,2), cas(2), tt(2,2), h2, s2, jq, res, eq2(-6:6), pyref(2), taus(2)
  integer :: i, it
  logical :: ok

  ok = .true.
  z = [ (0.5_dp,0.0_dp), (-1.0_dp,0.0_dp), (2.0_dp,0.0_dp), (0.3_dp,0.8_dp), &
       & (0.8104534588022096_dp,1.2622064772118446_dp), (0.9_dp,0.1_dp) ]
  ref = [ (0.5822405264650125_dp,0.0_dp), (-0.8224670334241132_dp,0.0_dp), &
       & (2.4674011002723395_dp,-2.177586090303602_dp), (0.12315425344338118_dp,0.8609654090425412_dp), &
       & (0.28850052402031634_dp,1.521257699097651_dp), (1.264186732338754_dp,0.24373567998101406_dp) ]
  do i = 1, 6
     ! Li2 at real z > 1: the branch with -i pi ln z (mpmath's) vs +i pi: compare real parts there
     if (aimag(z(i)) == 0 .and. real(z(i)) > 1) then
        call chk('cli2 real part (z > 1)', real(cli2(z(i))), real(ref(i)), 1e-13_dp)
     else
        call chk('cli2 re', real(cli2(z(i))), real(ref(i)), 1e-13_dp)
        call chk('cli2 im', aimag(cli2(z(i))), aimag(ref(i)), 1e-13_dp)
     endif
  enddo

  call soft_I0I1(0.5_dp, 0.7_dp, I0, I1)
  call chk('I0(0.5,0.7)', I0, 0.423937845282691_dp, 1e-10_dp)
  call chk('I1(0.5,0.7)', I1, -0.105451538601642_dp, 1e-10_dp)
  call soft_I0I1(1.0_dp, 1.0_dp, I0, I1)
  call chk('I0(1,1)', I0, 0.323065947219451_dp, 1e-10_dp)
  call chk('I1(1,1)', I1, -0.36554090374405_dp, 1e-10_dp)
  call soft_I0I1(2.0_dp, 0.3_dp, I0, I1)
  call chk('I0(2,0.3)', I0, 0.622381149336604_dp, 1e-10_dp)
  call chk('I1(2,0.3) (log point inside)', I1, -1.00628977909669_dp, 1e-10_dp)
  call soft_I0I1(0.3_dp, 2.5_dp, I0, I1)
  call chk('I0(0.3,2.5) empty region', I0, 0.0_dp, 1e-14_dp)
  ! extreme ratios (mpmath, 25 digits): narrow supports at large alpha
  call soft_I0I1(30.0_dp, 0.2_dp, I0, I1)
  call chk('I0 ln a + I1 (30, 0.2)', I0 * log(30.0_dp) + I1, -0.0336153051830945_dp, 1e-8_dp)
  call soft_I0I1(200.0_dp, 50.0_dp, I0, I1)
  call chk('I0 ln a + I1 (200, 50)', I0 * log(200.0_dp) + I1, -0.00500626392807682_dp, 1e-8_dp)
  call soft_I0I1(5000.0_dp, 4000.0_dp, I0, I1)
  call chk('I0 ln a + I1 (5000, 4000)', I0 * log(5000.0_dp) + I1, -0.000200010000888989_dp, 1e-8_dp)
  call soft_I0I1(0.01_dp, 0.02_dp, I0, I1)
  call chk('I0 ln a + I1 (0.01, 0.02)', I0 * log(0.01_dp) + I1, -7.66198710923967_dp, 1e-8_dp)
  call soft_I0I1(50.0_dp, 60.0_dp, I0, I1)
  call chk('I0 ln a + I1 (50, 60)', I0 * log(50.0_dp) + I1, -0.000502003785089934_dp, 1e-8_dp)

  call InitPDFsetByName('NNPDF30_nlo_as_0118')
  call InitPDF(0)
  x = 0.01_dp; Q = 20.0_dp; y = 0.5_dp
  a = alphasPDF(Q) / (2 * pi)
  Yp = 1 + (1 - y)**2
  eq2 = 0
  do i = 1, 5
     eq2(i) = merge(4.0_dp/9, 1.0_dp/9, mod(i,2) == 0); eq2(-i) = eq2(i)
  enddo
  nhat = 0; nhat(3,1) = 1; nhat(3,2) = -1
  cas = CF; tt = 0; tt(1,2) = -CF; tt(2,1) = -CF
  h2 = CF * (2 - pi**2) + 2 * (CF * pi**2 / 3 - 1.5_dp * CF - (3.5_dp - zeta2) * CF) + pi**2 / 12 * 2 * CF
  call chk('two-parton hard function = CF(-8 + pi^2/6)', h2, CF * (-8 + zeta2), 1e-14_dp)
  taus = [1e-2_dp, 1e-3_dp]
  pyref = [-1.6989899716e+02_dp, -3.7025011428e+02_dp]
  call EvolvePDF(x, Q, f1)
  call beam_coeffs(x, Q, c0, c1, c2)
  do it = 1, 2
     tau = taus(it)
     s2 = soft_cum(2, nhat, cas, tt, tau)
     jq = jet_cum(.false., tau)
     res = 0
     do i = -5, 5
        if (i == 0) cycle
        res = res + eq2(i) * ((h2 + jq + s2) * f1(i) + c0(i) + c1(i) * log(tau) + c2(i) * log(tau)**2)
     enddo
     res = a * Yp * res / x
     call chk('1+1 tau_1^a cumulant vs Python (KLS 173)', res, pyref(it), 2e-6_dp)
  enddo
  if (.not. ok) stop 1
  print *, 'all checks passed'
contains
  subroutine chk(tag, v, r, tol)
    character(len=*), intent(in) :: tag
    real(dp), intent(in) :: v, r, tol
    logical :: good
    good = abs(v - r) <= tol * max(1.0_dp, abs(r))
    write(*,'(a50,2es24.15,l3)') tag, v, r, good
    ok = ok .and. good
  end subroutine chk
end program test_scet
