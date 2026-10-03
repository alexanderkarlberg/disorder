!-----------------------------------------------------------------------
! Leading-power tau_2 cumulant of DIS 2+1 (photon exchange) at O(alpha_s)
! and O(alpha_s^2) relative to the Born, in the jets'-rest-frame
! geometric measure (slicing/mod_tau2_run.f90, measure 2), assembled with
! MCFM 10.3's 1-jettiness pieces (built by build.sh from an MCFM tree):
!   beam a  = beam function (xbeam1bis/xbeam2bis), z-integrated in a table in
!             xi at mu = mu_F = Q with Q_B = mu; Q_B = 2 E_a by the exact
!             cumulant shift s'(k) = sum_{m>=k} C(m,k) s(m) l^(m-k),
!             s'(-1) = s(-1) + sum_k s(k) l^(k+1)/(k+1), l = ln(Q_B/mu)
!             (the same as MCFM's jet functions do internally);
!   "beam b" = jet function of jet 2 (MCFM's jetq/jetg are in alpha_s/4pi:
!             J1/2, J2/4 in the beam slot; the jet slot converts itself);
!   jet     = jet function of jet 1;
!   soft    = soft_ab_* + soft_nab_* (MCFM's qgq: 1 = jet q, 2 = jet g,
!             3 = beam q; qag: 1 = jet q, 2 = jet qbar, 3 = beam g), with
!             all six I_ij,m from soft_G of slicing/mod_slicing_scet
!             (MCFM's computeIijm sets I13x2 = I23x1, not general);
!   hard    = hard21 (two loops, SCET, mu = Q).
! All energies and directions in the rest frame of the two jets (the
! Born's partonic CM frame), y_ij = n_i.n_j/2, Q_i = 2 E_i.
!
! lp21_born(P, xi, ntc, tcs, f0, c1, c2): P(4,7) the Born in DISENT's
! layout (1 incoming parton, 2 the outgoing quark, 3 the gluon for a quark
! Born; 2 quark, 3 antiquark for a gluon Born), xi its momentum fraction;
! tcs = tau_cut = T_cut/Q. Classes k = 1 up-type quarks (q + qbar, e^2
! weighted), 2 down-type, 3 gluon (weighted by sum e_q^2 over the produced
! flavour). f0(k): the Born luminosity sum e^2 f(xi) (MCFM's f = xf/x);
! c1, c2(it,k): the O(alpha_s/2pi) and O((alpha_s/2pi)^2) coefficients of the
! cumulant in the same units (divide by f0 for per-Born factors).
!-----------------------------------------------------------------------
subroutine fdist(ih, x, xmu, fx, ibeam)
  implicit none
  integer :: ih, ibeam
  double precision :: x, xmu, fx(-5:5), f(-6:6)
  call evolvePDF(x, xmu, f)
  fx = f(-5:5) / x
end subroutine fdist

module lp21
  use SCET_Jet, only: jetq, jetg
  use mod_slicing_scet, only: soft_G
  use hard21, only: hard21_eval
  implicit none
  private
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: eu = 2.0_dp/3, ed = -1.0_dp/3
  ! beam table: classes 1..3, coefficients b0, b1(-1:1), b2(-1:3) = 9 numbers
  integer, save :: nt = 0
  real(dp), save :: tQ = 0, t0, th
  real(dp), allocatable, save :: tab(:,:,:)       ! (0:nt, 9, 3)
  real(dp), save :: zx                              ! xi of the z integrand
  integer, save :: zlim = 30
  public :: lp21_init, lp21_born, lp21_beam_direct, lp21_shift
contains

  ! MCFM's common blocks at mu = mu_F = Q; table of the beam coefficients
  ! for ximin <= xi < 1, nodes equally spaced in ln(xi/(1-xi)) by h
  subroutine lp21_init(Q, ximin, h)
    real(dp), intent(in) :: Q, ximin, h
    real(dp) :: scale, musq, facscale, gsq, as, ason2pi, ason4pi, t1, t, b(9,3)
    integer :: nflav, kpart, i
    logical :: coeffonly
    common/mcfmscale/scale, musq
    common/facscale/facscale
    common/nflav/nflav
    common/qcdcouple/gsq, as, ason2pi, ason4pi
    common/coeffonly/coeffonly
    common/kpart/kpart
    scale = Q; musq = Q*Q; facscale = Q; nflav = 5
    ason2pi = 1; ason4pi = 0.5_dp; coeffonly = .true.; kpart = 0
    tQ = Q
    if (h <= 0) then
       nt = 0; return
    endif
    t0 = log(ximin/(1 - ximin)); t1 = log((1 - 1e-6_dp)/1e-6_dp); th = h
    nt = int((t1 - t0)/h) + 1
    if (allocated(tab)) deallocate(tab)
    allocate(tab(0:nt, 9, 3))
    do i = 0, nt
       t = t0 + h*i
       call lp21_beam_direct(1/(1 + exp(-t)), b)
       tab(i,:,:) = b
    enddo
  end subroutine lp21_init

  ! beam coefficients at xi (mu = Q, Q_B = mu) by adaptive z integration
  subroutine lp21_beam_direct(xi, b)
    real(dp), intent(in) :: xi
    real(dp), intent(out) :: b(9,3)
    real(dp) :: fx(-5:5), r(24), e
    call fdist(1, xi, tQ, fx, 1)
    b = 0
    b(1,1) = eu**2*(fx(2) + fx(-2) + fx(4) + fx(-4))
    b(1,2) = ed**2*(fx(1) + fx(-1) + fx(3) + fx(-3) + fx(5) + fx(-5))
    b(1,3) = fx(0)
    zx = xi
    call adapt_rec(0.0_dp, 1.0_dp, r, e, 0)
    b(2:9,1) = r(1:8); b(2:9,2) = r(9:16); b(2:9,3) = r(17:24)
  end subroutine lp21_beam_direct

  subroutine integrand(z, f)
    real(dp), intent(in) :: z
    real(dp), intent(out) :: f(24)
    real(dp) :: bt1(-5:5,-1:1), bt2(-5:5,-1:3), w(-5:5,3)
    integer :: k, c
    call xbeam1bis(1, z, zx, tQ, bt1, 1)
    call xbeam2bis(1, z, zx, tQ, bt2, 1)
    w = 0
    w([2, -2, 4, -4], 1) = eu**2
    w([1, -1, 3, -3, 5, -5], 2) = ed**2
    w(0, 3) = 1
    do c = 1, 3
       do k = -1, 1
          f(8*(c-1) + k + 2) = sum(w(:,c)*bt1(:,k))
       enddo
       do k = -1, 3
          f(8*(c-1) + k + 5) = sum(w(:,c)*bt2(:,k))
       enddo
    enddo
  end subroutine integrand

  subroutine gk15(a, b, rk, rg)
    real(dp), intent(in) :: a, b
    real(dp), intent(out) :: rk(24), rg(24)
    real(dp), parameter :: xk(8) = [0.991455371120812639_dp, 0.949107912342758525_dp, &
         0.864864423359769073_dp, 0.741531185599394440_dp, 0.586087235467691130_dp, &
         0.405845151377397167_dp, 0.207784955007898468_dp, 0.0_dp]
    real(dp), parameter :: wk(8) = [0.022935322010529225_dp, 0.063092092629978553_dp, &
         0.104790010322250184_dp, 0.140653259715525919_dp, 0.169004726639267903_dp, &
         0.190350578064785410_dp, 0.204432940075298892_dp, 0.209482141084727828_dp]
    real(dp), parameter :: wg(4) = [0.129484966168869693_dp, 0.279705391489276668_dp, &
         0.381830050505118945_dp, 0.417959183673469388_dp]
    real(dp) :: c, h, fa(24), fb(24)
    integer :: jj
    c = 0.5_dp*(a + b); h = 0.5_dp*(b - a)
    call integrand(c, fa)
    rk = wk(8)*fa; rg = wg(4)*fa
    do jj = 1, 7
       call integrand(c - h*xk(jj), fa); call integrand(c + h*xk(jj), fb)
       rk = rk + wk(jj)*(fa + fb)
       if (mod(jj, 2) == 0) rg = rg + wg(jj/2)*(fa + fb)
    enddo
    rk = rk*h; rg = rg*h
  end subroutine gk15

  recursive subroutine adapt_rec(a, b, r, e, depth)
    real(dp), intent(in) :: a, b
    real(dp), intent(out) :: r(24), e
    integer, intent(in) :: depth
    real(dp) :: rk(24), rg(24), r1(24), r2(24), e1, e2, m, sc
    call gk15(a, b, rk, rg)
    sc = max(maxval(abs(rk)), 1e-300_dp)
    e = maxval(abs(rk - rg))/sc
    if (depth > zlim .or. e < 1e-10_dp) then
       r = rk; return
    endif
    m = 0.5_dp*(a + b)
    call adapt_rec(a, m, r1, e1, depth + 1); call adapt_rec(m, b, r2, e2, depth + 1)
    r = r1 + r2; e = e1 + e2
  end subroutine adapt_rec

  ! beam coefficients at xi from the table (Catmull-Rom), or direct
  subroutine beam_at(xi, b)
    real(dp), intent(in) :: xi
    real(dp), intent(out) :: b(9,3)
    real(dp) :: t, x, fx, w(4)
    integer :: ix, k
    t = log(xi/(1 - xi))
    if (nt > 0 .and. t >= t0 + th .and. t <= t0 + th*(nt - 2)) then
       x = (t - t0)/th; ix = int(x); fx = x - ix
       w = [-0.5_dp*fx**3 + fx**2 - 0.5_dp*fx, 1.5_dp*fx**3 - 2.5_dp*fx**2 + 1, &
            -1.5_dp*fx**3 + 2*fx**2 + 0.5_dp*fx, 0.5_dp*fx**3 - 0.5_dp*fx**2]
       b = 0
       do k = 1, 4
          b = b + w(k)*tab(ix + k - 2,:,:)
       enddo
    else
       call lp21_beam_direct(xi, b)
    endif
  end subroutine beam_at

  ! shift of distribution coefficients s(-1:n) from logs of tau Q_i/mu^2 to
  ! logs of tau/mu: l = ln(Q_i/mu)
  subroutine lp21_shift(n, s, l)
    integer, intent(in) :: n
    real(dp), intent(inout) :: s(-1:n)
    real(dp), intent(in) :: l
    real(dp) :: t(-1:n)
    integer :: k, m
    t = 0
    do k = 0, n
       do m = k, n
          t(k) = t(k) + s(m)*binom(m, k)*l**(m - k)
       enddo
    enddo
    t(-1) = s(-1)
    do k = 0, n
       t(-1) = t(-1) + s(k)*l**(k + 1)/(k + 1)
    enddo
    s = t
  contains
    real(dp) function binom(m, k)
      integer, intent(in) :: m, k
      integer :: i
      binom = 1
      do i = 1, k
         binom = binom*(m - k + i)/i
      enddo
    end function binom
  end subroutine lp21_shift

  subroutine lp21_born(P, xi, ntc, tcs, f0, c1, c2)
    real(dp), intent(in) :: P(4,7), xi, tcs(:)
    integer, intent(in) :: ntc
    real(dp), intent(out) :: f0(3), c1(ntc,3), c2(ntc,3)
    real(dp) :: pb(4,3), y, ch, sh, Ea, E2, E3, nh(3,3), yy(3,3), Q
    real(dp) :: b(9,3), ba0, ba1(-1:1), ba2(-1:3), lB
    real(dp) :: Jq1(-1:1), Jq2(-1:3), Jg1(-1:1), Jg2(-1:3), Jb1(-1:1), Jb2(-1:3), Jc1(-1:1), Jc2(-1:3)
    real(dp) :: s1(-1:1), s2(-1:3), s2n(-1:3), Iq(6), Ig(6), hard(2), h(0:2), hz(0:2), tc
    real(dp), external :: assemblejet
    integer :: k, it, c
    Q = tQ
    ! the jets' rest frame: boost along z
    y = -atanh((P(3,2) + P(3,3))/(P(4,2) + P(4,3)))
    ch = cosh(y); sh = sinh(y)
    do k = 1, 3
       pb(:,k) = P(:,k)
       pb(4,k) = ch*P(4,k) + sh*P(3,k); pb(3,k) = sh*P(4,k) + ch*P(3,k)
    enddo
    Ea = pb(4,1); E2 = pb(4,2); E3 = pb(4,3)
    do k = 1, 3
       nh(:,k) = pb(1:3,k)/sqrt(sum(pb(1:3,k)**2))
    enddo
    do k = 1, 3
       do c = 1, 3
          yy(k,c) = 0.5_dp*(1 - dot_product(nh(:,k), nh(:,c)))
       enddo
    enddo
    ! beam
    call beam_at(xi, b)
    lB = log(2*Ea/Q)
    ! soft: MCFM labels (1, 2, 3) = (parton 2, parton 3, beam)
    call iijm6(yy(2,3), yy(3,1), yy(1,2), Iq)
    Ig = Iq
    ! jet functions (alpha_s/4pi, logs shifted to tau/mu by MCFM)
    call jetq(2, 2*E2, Jq1, Jq2, Q)
    call jetg(2, 2*E3, Jg1, Jg2, Q)
    do c = 1, 3
       f0(c) = b(1,c)
       ba0 = b(1,c); ba1 = b(2:4,c); ba2 = b(5:9,c)
       call lp21_shift(1, ba1, lB); call lp21_shift(3, ba2, lB)
       if (c <= 2) then
          ! quark Born: jet q (2) in the beam-b slot, jet g (3) in the jet slot
          call soft_ab_qgq(2, yy(2,3), yy(3,1), yy(1,2), Iq, 1, 2, 3, 4, 5, 6, s1, s2)
          call soft_nab_qgq(2, yy(2,3), yy(3,1), yy(1,2), Iq, 1, 2, 3, 4, 5, 6, s2n)
          Jb1 = Jq1/2; Jb2 = Jq2/4; Jc1 = Jg1; Jc2 = Jg2
          call hard21_eval(P, 1, merge(eu, ed, c == 1), h)
       else
          ! gluon Born: jets q (2), qbar (3); n_f,gamma term summed over the
          ! produced flavour: sum e^2 (sum e/e) = (sum e)^2 -> eq = sum e^2/sum e
          call soft_ab_qag(2, yy(2,3), yy(3,1), yy(1,2), Ig, 1, 2, 3, 4, 5, 6, s1, s2)
          call soft_nab_qag(2, yy(2,3), yy(3,1), yy(1,2), Ig, 1, 2, 3, 4, 5, 6, s2n)
          call jetq(2, 2*E3, Jc1, Jc2, Q)
          Jb1 = Jq1/2; Jb2 = Jq2/4
          call hard21_eval(P, 2, (11.0_dp/9)/(1.0_dp/3), h)
       endif
       s2 = s2 + s2n
       hard = h(1:2)
       do it = 1, ntc
          tc = tcs(it)*Q
          c1(it,c) = assemblejet(1, tc, ba0, 1.0_dp, ba1, Jb1, ba2, Jb2, s1, s2, Jc1, Jc2, hard)
          c2(it,c) = assemblejet(2, tc, ba0, 1.0_dp, ba1, Jb1, ba2, Jb2, s1, s2, Jc1, Jc2, hard)
       enddo
    enddo
  end subroutine lp21_born

  ! I_ij,m = I0 ln(alpha) + I1 for MCFM's six (i,j,m), alpha = y_jm/y_ij,
  ! beta = y_im/y_ij (soft_G of mod_slicing_scet, checked against mpmath)
  subroutine iijm6(y12, y23, y31, I)
    real(dp), intent(in) :: y12, y23, y31
    real(dp), intent(out) :: I(6)
    integer, parameter :: ii(6) = [1, 2, 2, 3, 3, 1], jj(6) = [2, 1, 3, 2, 1, 3], mm(6) = [3, 3, 1, 1, 2, 2]
    real(dp) :: yv(3,3)
    integer :: n
    yv = 0
    yv(1,2) = y12; yv(2,1) = y12; yv(2,3) = y23; yv(3,2) = y23; yv(1,3) = y31; yv(3,1) = y31
    do n = 1, 6
       I(n) = soft_G(yv(jj(n),mm(n))/yv(ii(n),jj(n)), yv(ii(n),mm(n))/yv(ii(n),jj(n)))
    enddo
  end subroutine iijm6
end module lp21
