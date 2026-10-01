! Leading-power tau_1 cumulant for DIS 1+1 (photon exchange) at NLO and NNLO,
! assembled with MCFM 10.3's SCET pieces (DY 0-jettiness structure):
!   tau_1 = min over partitions [2 x P.p(beam region) + m^2(jet region)] / Q^2
!   (recoil-free, thrust-like in the Breit frame), tau_hat = Q tau_1 =
!   t_B/Q + t_J/Q + k_s, Q_B = Q_J = Q, mu_R = mu_F = Q, so L = ln(tau_1,cut).
! beam a = quark beam function (xbeam1bis/xbeam2bis, e_q^2 weighted),
! "beam b" = quark jet function (jetq), soft = qqbar 0-jettiness soft function
! (equal to the DIS one at O(alpha_s^2), Kang-Labun-Lee 1501.04110),
! hard = hardqq at Qsq = -Q^2 (spacelike: real logs, |C_V(Q^2)|^2).
! Output: coefficients of (alpha_s/2pi)^k per Born (divided by sum e_q^2 f_q(x)).
subroutine fdist(ih, x, xmu, fx, ibeam)
  implicit none
  integer :: ih, ibeam
  double precision :: x, xmu, fx(-5:5), f(-6:6)
  integer :: pdfmask
  common/lpmask/pdfmask
  call evolvePDF(x, xmu, f)
  fx = f(-5:5) / x
  if (pdfmask == 1) fx(0) = 0
  if (pdfmask == 2) then
     fx(-5:-1) = 0; fx(1:5) = 0
  endif
end subroutine fdist

program dis_tau1_lp
  use SCET_Jet, only: jetq
  implicit none
  integer, parameter :: dp = kind(1d0)
  real(dp) :: x, Q2, Q, scale, musq, facscale, gsq, as, ason2pi, ason4pi
  integer :: nflav, pdfmask, kpart, i, j, it, nargs
  logical :: coeffonly
  common/mcfmscale/scale, musq
  common/facscale/facscale
  common/nflav/nflav
  common/qcdcouple/gsq, as, ason2pi, ason4pi
  common/coeffonly/coeffonly
  common/kpart/kpart
  common/lpmask/pdfmask
  real(dp) :: fx(-5:5), e2(-5:5), b0, b1(-1:1), b2(-1:3), J1(-1:1), J2(-1:3)
  real(dp) :: s1(-1:1), s2(-1:3), hard(2), taus(12), c1, c2, assemble
  real(dp) :: res(8), err
  character(len=100) :: pdfname, arg
  external assemble
  data taus /1e-1_dp, 3e-2_dp, 1e-2_dp, 3e-3_dp, 1e-3_dp, 3e-4_dp, 1e-4_dp, 3e-5_dp, 1e-5_dp, 1e-6_dp, 1e-7_dp, 1e-8_dp/
  pdfname = 'NNPDF30_nlo_as_0118'; x = 0.01_dp; Q2 = 400; pdfmask = 0
  nargs = command_argument_count()
  i = 1
  do while (i <= nargs)
     call get_command_argument(i, arg)
     select case (trim(arg))
     case ('-x'); call get_command_argument(i+1, arg); read(arg,*) x; i = i + 1
     case ('-Q2'); call get_command_argument(i+1, arg); read(arg,*) Q2; i = i + 1
     case ('-pdf'); call get_command_argument(i+1, pdfname); i = i + 1
     case ('-pdfmask'); call get_command_argument(i+1, arg); read(arg,*) pdfmask; i = i + 1
     end select
     i = i + 1
  enddo
  call initPDFsetByName(trim(pdfname)); call initPDF(0)
  Q = sqrt(Q2); scale = Q; musq = Q2; facscale = Q; nflav = 5
  ason2pi = 1; ason4pi = 0.5_dp; coeffonly = .true.; kpart = 0
  e2 = [1,4,1,4,1, 0, 1,4,1,4,1] / 9.0_dp
  call fdist(1, x, Q, fx, 1)
  b0 = sum(e2 * fx)
  call adapt(0.0_dp, 1.0_dp, res, err)
  b1 = res(1:3); b2 = res(4:8)
  call jetq(2, Q, J1, J2, Q)
  ! jetq is normalised to alpha_s/4pi (MCFM's assemblejet uses J1/2, J2/4)
  J1 = J1 / 2; J2 = J2 / 4
  call softqqbis(2, s1, s2)
  call hardqq(-Q2, Q2, hard)
  write(*,'(a,a,a,f8.5,a,f9.2,a,i2)') ' pdf ', trim(pdfname), '  x ', x, '  Q2 ', Q2, '  pdfmask ', pdfmask
  write(*,'(a,es14.6,a,es9.2)') ' Born sum e^2 f(x) ', b0, '   beam z-integration error estimate ', err
  write(*,'(a,3es14.6)') ' beam1 (-1:1)/b0 ', b1 / b0
  write(*,'(a,5es14.6)') ' beam2 (-1:3)/b0 ', b2 / b0
  write(*,'(a,3es14.6,a,5es14.6)') ' jet1 ', J1, '  jet2 ', J2
  write(*,'(a,3es14.6,a,5es14.6)') ' soft1 ', s1, '  soft2 ', s2
  write(*,'(a,2es14.6)') ' hard (spacelike) ', hard
  write(*,'(a)') '   tau_cut        LP O(as/2pi) per Born     LP O((as/2pi)^2) per Born'
  do it = 1, size(taus)
     c1 = assemble(1, taus(it) * Q, b0, 1.0_dp, b1, J1, b2, J2, s1, s2, hard) / b0
     c2 = assemble(2, taus(it) * Q, b0, 1.0_dp, b1, J1, b2, J2, s1, s2, hard) / b0
     write(*,'(es10.2,2es24.12)') taus(it), c1, c2
  enddo
contains
  subroutine integrand(z, f)
    real(dp), intent(in) :: z
    real(dp), intent(out) :: f(8)
    real(dp) :: bt1(-5:5,-1:1), bt2(-5:5,-1:3)
    integer :: k
    call xbeam1bis(1, z, x, Q, bt1, 1)
    call xbeam2bis(1, z, x, Q, bt2, 1)
    do k = -1, 1
       f(k+2) = sum(e2 * bt1(:,k))
    enddo
    do k = -1, 3
       f(k+5) = sum(e2 * bt2(:,k))
    enddo
  end subroutine integrand
  subroutine gk15(a, b, rk, rg)
    real(dp), intent(in) :: a, b
    real(dp), intent(out) :: rk(8), rg(8)
    real(dp), parameter :: xk(8) = [0.991455371120812639_dp, 0.949107912342758525_dp, &
         0.864864423359769073_dp, 0.741531185599394440_dp, 0.586087235467691130_dp, &
         0.405845151377397167_dp, 0.207784955007898468_dp, 0.0_dp]
    real(dp), parameter :: wk(8) = [0.022935322010529225_dp, 0.063092092629978553_dp, &
         0.104790010322250184_dp, 0.140653259715525919_dp, 0.169004726639267903_dp, &
         0.190350578064785410_dp, 0.204432940075298892_dp, 0.209482141084727828_dp]
    real(dp), parameter :: wg(4) = [0.129484966168869693_dp, 0.279705391489276668_dp, &
         0.381830050505118945_dp, 0.417959183673469388_dp]
    real(dp) :: c, h, fa(8), fb(8)
    integer :: jj
    c = 0.5_dp * (a + b); h = 0.5_dp * (b - a)
    call integrand(c, fa)
    rk = wk(8) * fa; rg = wg(4) * fa
    do jj = 1, 7
       call integrand(c - h * xk(jj), fa); call integrand(c + h * xk(jj), fb)
       rk = rk + wk(jj) * (fa + fb)
       if (mod(jj, 2) == 0) rg = rg + wg(jj/2) * (fa + fb)
    enddo
    rk = rk * h; rg = rg * h
  end subroutine gk15
  recursive subroutine adapt_rec(a, b, r, e, depth)
    real(dp), intent(in) :: a, b
    real(dp), intent(out) :: r(8), e
    integer, intent(in) :: depth
    real(dp) :: rk(8), rg(8), r1(8), r2(8), e1, e2_, m
    call gk15(a, b, rk, rg)
    e = maxval(abs(rk - rg)) / b0
    if (depth > 30 .or. e < 1e-12_dp) then
       r = rk; return
    endif
    m = 0.5_dp * (a + b)
    call adapt_rec(a, m, r1, e1, depth + 1); call adapt_rec(m, b, r2, e2_, depth + 1)
    r = r1 + r2; e = e1 + e2_
  end subroutine adapt_rec
  subroutine adapt(a, b, r, e)
    real(dp), intent(in) :: a, b
    real(dp), intent(out) :: r(8), e
    call adapt_rec(a, b, r, e, 0)
  end subroutine adapt
end program dis_tau1_lp
