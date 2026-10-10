!-----------------------------------------------------------------------
! Leading-power tau_2 cumulant of DIS 2+1 (photon; photon + Z) at O(alpha_s)
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
! lp21_born(P, xi, ntc, tcs, wk, nfk, f0, c1, c2): P(4,7) the Born in DISENT's
! layout (1 incoming parton, 2 the outgoing quark, 3 the gluon for a quark
! Born; 2 quark, 3 antiquark for a gluon Born), xi its momentum fraction;
! tcs = tau_cut = T_cut/Q. Beam classes c = 1 up-type quarks, 2 up-type
! antiquarks, 3 down-type quarks, 4 down-type antiquarks, 5 gluon (the
! beam table unweighted); wk(k,c) the couplings of class c for the two
! helicity classes k of hard21 (photon: e_q^2 for both; gluon: summed over
! the produced flavour), combined with hard21's Born fractions of k (with
! photon + Z the hard functions of the two classes differ); nfk(k,c) the
! two-loop N_F,V coupling ratio (the boson on a closed quark loop, hard21's
! G term: photon sum_q e_q / e_q, gluon (sum e)^2/sum e^2; with Z the
! vector couplings of the loop, sliced21), H^(2) linear in it. f0(c): the
! Born luminosity f(xi) (MCFM's f = xf/x) times the combined coupling;
! c1, c2(it,c): the O(alpha_s/2pi) and O((alpha_s/2pi)^2) coefficients of
! the cumulant in the same units (divide by f0 for per-Born factors).
!-----------------------------------------------------------------------
subroutine fdist(ih, x, xmu, fx, ibeam)
  implicit none
  integer :: ih, ibeam
  double precision :: x, xmu, fx(-5:5), f(-6:6)
  integer :: lpmask
  common/lpmask/lpmask
  call evolvePDF(x, xmu, f)
  fx = f(-5:5) / x
  if (lpmask == 1) fx(0) = 0
  if (lpmask == 2) then
     fx(-5:-1) = 0; fx(1:5) = 0
  endif
end subroutine fdist

module lp21
  use SCET_Jet, only: jetq, jetg
  use mod_slicing_scet, only: soft_G
  use hard21, only: hard21_eval, hard21_nfz
  implicit none
  private
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: eu = 2.0_dp/3, ed = -1.0_dp/3
  ! beam classes (lp21_born) and their flavours (q; qbar)
  integer, parameter, public :: nbc = 5
  integer, parameter :: cup(2) = [2, 4], cdn(3) = [1, 3, 5]
  ! beam table: classes 1..nbc, coefficients b0, b1(-1:1), b2(-1:3) = 9 numbers
  integer, save :: nt = 0
  real(dp), save :: tQ = 0, t0, th
  real(dp), allocatable, save :: tab(:,:,:)       ! (0:nt, 9, nbc)
  real(dp), save :: zx                              ! xi of the z integrand
  integer, save :: zlim = 30
  ! grid in Q (sliced21 mode 2, Q = mu varies per event): nodes ln Q =
  ! gq0 + gqh*i, i = 0..nqg-1, each a table as tab; cubic in ln Q
  integer, save :: nqg = 0
  real(dp), save :: gq0, gqh
  real(dp), allocatable, save :: tabg(:,:,:,:)    ! (0:nt, 9, nbc, 0:nqg-1)
  public :: lp21_init, lp21_born, lp21_beam_direct, lp21_shift, beam_at
  public :: lp21_grid_build, lp21_grid_load, lp21_setq
  ! W exchange (9 Oct): b without coupling, so not in the down-type classes;
  ! set before lp21_init/grid_build/grid_load. Grid tables built with it
  ! carry an extra flag word (the loader refuses a mismatch)
  logical, public :: lp21_nodn = .false.
contains

  ! MCFM's common blocks at mu = mu_F = Q; table of the beam coefficients
  ! for ximin <= xi < 1, nodes equally spaced in ln(xi/(1-xi)) by h
  subroutine lp21_init(Q, ximin, h)
    real(dp), intent(in) :: Q, ximin, h
    real(dp) :: t1, t, b(9,nbc)
    integer :: i
    call lp21_setq(Q)
    if (h <= 0) then
       nt = 0; return
    endif
    t0 = log(ximin/(1 - ximin)); t1 = log((1 - 1e-6_dp)/1e-6_dp); th = h
    nt = int((t1 - t0)/h) + 1
    if (allocated(tab)) deallocate(tab)
    allocate(tab(0:nt, 9, nbc))
    do i = 0, nt
       t = t0 + h*i
       call lp21_beam_direct(1/(1 + exp(-t)), b)
       tab(i,:,:) = b
    enddo
  end subroutine lp21_init

  ! MCFM's common blocks and the beam scale at mu = mu_F = Q (per event in
  ! sliced21 mode 2)
  subroutine lp21_setq(Q)
    real(dp), intent(in) :: Q
    real(dp) :: scale, musq, facscale, gsq, as, ason2pi, ason4pi
    integer :: nflav, kpart
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
  end subroutine lp21_setq

  ! Q grid: nodes ln Q = ln Qlo - dl + dl*i, i = 0..nq-1 (one node beyond each
  ! end for the cubic interpolation), nq = ceiling(ln(Qhi/Qlo)/dl) + 3.
  ! build: tables of nodes i0..i1 into files <prefix>_<i>.tab
  subroutine lp21_grid_build(Qlo, Qhi, dl, ximin, h, i0, i1, prefix)
    real(dp), intent(in) :: Qlo, Qhi, dl, ximin, h
    integer, intent(in) :: i0, i1
    character(*), intent(in) :: prefix
    integer :: i, u, nq
    character(16) :: num
    nq = ceiling(log(Qhi/Qlo)/dl) + 3
    do i = max(i0, 0), min(i1, nq - 1)
       call lp21_init(exp(log(Qlo) - dl + dl*i), ximin, h)
       write(num, '(i0)') i
       open(newunit=u, file=prefix//'_'//trim(num)//'.tab', form='unformatted', access='stream', status='replace')
       write(u) nt, t0, th, tQ
       write(u) tab
       if (lp21_nodn) write(u) 1
       close(u)
    enddo
  end subroutine lp21_grid_build

  ! load the grid <prefix>_0.tab, _1.tab, ... (spacing in ln Q and the number
  ! of nodes from the files; the Q range must be covered, Qlo/Qhi checked)
  subroutine lp21_grid_load(Qlo, Qhi, prefix)
    real(dp), intent(in) :: Qlo, Qhi
    character(*), intent(in) :: prefix
    integer :: i, u, n
    integer(8) :: sz
    real(dp) :: a, b, q, q0
    character(16) :: num
    logical :: ex
    nqg = 0
    do
       write(num, '(i0)') nqg
       inquire(file=prefix//'_'//trim(num)//'.tab', exist=ex)
       if (.not. ex) exit
       nqg = nqg + 1
    enddo
    if (nqg < 4) stop 'lp21_grid_load: fewer than 4 Q nodes'
    do i = 0, nqg - 1
       write(num, '(i0)') i
       open(newunit=u, file=prefix//'_'//trim(num)//'.tab', form='unformatted', access='stream', status='old')
       read(u) n, a, b, q
       inquire(unit=u, size=sz)
       if (sz /= 28 + 8_8*(n + 1)*9*nbc + merge(4, 0, lp21_nodn)) &
            & stop 'lp21_grid_load: not a table of nbc beam classes for this boson (old format, or b in the down classes: rebuild)'
       if (i == 0) then
          nt = n; t0 = a; th = b; q0 = q
          if (allocated(tabg)) deallocate(tabg)
          allocate(tabg(0:nt, 9, nbc, 0:nqg-1))
       elseif (n /= nt .or. a /= t0 .or. b /= th) then
          stop 'lp21_grid_load: tables differ in xi nodes'
       endif
       if (i == 1) then
          gq0 = log(q0); gqh = log(q) - gq0
       elseif (i > 1) then
          if (abs(log(q) - (gq0 + gqh*i)) > 1e-10_dp) stop 'lp21_grid_load: Q nodes not equally spaced'
       endif
       read(u) tabg(:,:,:,i)
       close(u)
    enddo
    if (log(Qlo) < gq0 + gqh .or. log(Qhi) > gq0 + gqh*(nqg - 2)) stop 'lp21_grid_load: Q range not covered'
  end subroutine lp21_grid_load

  ! beam coefficients at xi (mu = Q, Q_B = mu) by adaptive z integration
  subroutine lp21_beam_direct(xi, b)
    real(dp), intent(in) :: xi
    real(dp), intent(out) :: b(9,nbc)
    real(dp) :: fx(-5:5), r(8*nbc), e, w(-5:5,nbc)
    integer :: c
    call fdist(1, xi, tQ, fx, 1)
    call class_masks(w)
    zx = xi
    call adapt_rec(0.0_dp, 1.0_dp, r, e, 0)
    do c = 1, nbc
       b(1,c) = sum(w(:,c)*fx)
       b(2:9,c) = r(8*(c-1) + 1:8*c)
    enddo
  end subroutine lp21_beam_direct

  subroutine integrand(z, f)
    real(dp), intent(in) :: z
    real(dp), intent(out) :: f(8*nbc)
    real(dp) :: bt1(-5:5,-1:1), bt2(-5:5,-1:3), w(-5:5,nbc)
    integer :: k, c
    call xbeam1bis(1, z, zx, tQ, bt1, 1)
    call xbeam2bis(1, z, zx, tQ, bt2, 1)
    call class_masks(w)
    do c = 1, nbc
       do k = -1, 1
          f(8*(c-1) + k + 2) = sum(w(:,c)*bt1(:,k))
       enddo
       do k = -1, 3
          f(8*(c-1) + k + 5) = sum(w(:,c)*bt2(:,k))
       enddo
    enddo
  end subroutine integrand

  ! the Born flavours of each beam class
  subroutine class_masks(w)
    real(dp), intent(out) :: w(-5:5,nbc)
    w = 0
    w(cup, 1) = 1; w(-cup, 2) = 1
    w(cdn, 3) = 1; w(-cdn, 4) = 1
    if (lp21_nodn) then
       w(5, 3) = 0; w(-5, 4) = 0
    endif
    w(0, 5) = 1
  end subroutine class_masks

  subroutine gk15(a, b, rk, rg)
    real(dp), intent(in) :: a, b
    real(dp), intent(out) :: rk(8*nbc), rg(8*nbc)
    real(dp), parameter :: xk(8) = [0.991455371120812639_dp, 0.949107912342758525_dp, &
         0.864864423359769073_dp, 0.741531185599394440_dp, 0.586087235467691130_dp, &
         0.405845151377397167_dp, 0.207784955007898468_dp, 0.0_dp]
    real(dp), parameter :: wk(8) = [0.022935322010529225_dp, 0.063092092629978553_dp, &
         0.104790010322250184_dp, 0.140653259715525919_dp, 0.169004726639267903_dp, &
         0.190350578064785410_dp, 0.204432940075298892_dp, 0.209482141084727828_dp]
    real(dp), parameter :: wg(4) = [0.129484966168869693_dp, 0.279705391489276668_dp, &
         0.381830050505118945_dp, 0.417959183673469388_dp]
    real(dp) :: c, h, fa(8*nbc), fb(8*nbc)
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
    real(dp), intent(out) :: r(8*nbc), e
    integer, intent(in) :: depth
    real(dp) :: rk(8*nbc), rg(8*nbc), r1(8*nbc), r2(8*nbc), e1, e2, m, sc
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
    real(dp), intent(out) :: b(9,nbc)
    real(dp) :: t, x, fx, w(4)
    integer :: ix, k
    t = log(xi/(1 - xi))
    if (nqg > 0) then
       call grid_at(t, b)
       return
    endif
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

  ! Catmull-Rom in t = ln(xi/(1-xi)) and ln Q on the grid; direct outside it
  subroutine grid_at(t, b)
    real(dp), intent(in) :: t
    real(dp), intent(out) :: b(9,nbc)
    real(dp) :: x, fx, w(4), xq, fq, wq(4)
    integer :: ix, iq, k, l
    xq = (log(tQ) - gq0)/gqh
    if (t < t0 + th .or. t > t0 + th*(nt - 2) .or. xq < 1 .or. xq > nqg - 2) then
       call lp21_beam_direct(1/(1 + exp(-t)), b)
       return
    endif
    x = (t - t0)/th; ix = int(x); fx = x - ix
    iq = min(int(xq), nqg - 3); fq = xq - iq
    w = cr(fx); wq = cr(fq)
    b = 0
    do l = 1, 4
       do k = 1, 4
          b = b + wq(l)*w(k)*tabg(ix + k - 2,:,:,iq + l - 2)
       enddo
    enddo
  contains
    function cr(f) result(c)
      real(dp), intent(in) :: f
      real(dp) :: c(4)
      c = [-0.5_dp*f**3 + f**2 - 0.5_dp*f, 1.5_dp*f**3 - 2.5_dp*f**2 + 1, &
           -1.5_dp*f**3 + 2*f**2 + 0.5_dp*f, 0.5_dp*f**3 - 0.5_dp*f**2]
    end function cr
  end subroutine grid_at

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

  subroutine lp21_born(P, xi, ntc, tcs, wk, nfk, f0, c1, c2)
    real(dp), intent(in) :: P(4,7), xi, tcs(:), wk(2,nbc), nfk(2,nbc)
    integer, intent(in) :: ntc
    real(dp), intent(out) :: f0(nbc), c1(ntc,nbc), c2(ntc,nbc)
    real(dp) :: pb(4,3), y, ch, sh, Ea, E2, E3, nh(3,3), yy(3,3), Q
    real(dp) :: b(9,nbc), ba0, ba1(-1:1), ba2(-1:3), lB
    real(dp) :: Jq1(-1:1), Jq2(-1:3), Jg1(-1:1), Jg2(-1:3), Jb1(-1:1), Jb2(-1:3), Jc1(-1:1), Jc2(-1:3)
    real(dp) :: s1(-1:1), s2(-1:3), s2n(-1:3), Iq(6), Ig(6), h(0:2), hk(0:2,2,2), hn(0:2,2), gk(2,2), hc(2), tc, rk(2)
    real(dp), external :: assemblejet
    integer :: k, it, c, ht
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
    ! hard functions per helicity class, quark (1) and gluon (2) Born, with
    ! the N_F,V term off and with coupling ratio 1 (hard21: (1/3)/eq);
    ! H^(2) of class k, beam class c: hk(2,k) + nfk(k,c) hg(k)
    hard21_nfz = .false.
    call hard21_eval(P, 1, 1.0_dp/3, h, hk(:,:,1))
    call hard21_eval(P, 2, 1.0_dp/3, h, hk(:,:,2))
    hard21_nfz = .true.
    do ht = 1, 2
       call hard21_eval(P, ht, 1.0_dp/3, h, hn)
       gk(:,ht) = hn(2,:) - hk(2,:,ht)
    enddo
    do c = 1, nbc
       ht = merge(1, 2, c <= 4)
       rk = hk(0,:,ht)/sum(hk(0,:,ht))*wk(:,c)
       f0(c) = b(1,c)*sum(rk)
       ba0 = b(1,c); ba1 = b(2:4,c); ba2 = b(5:9,c)
       call lp21_shift(1, ba1, lB); call lp21_shift(3, ba2, lB)
       if (c <= 4) then
          ! quark Born: jet q (2) in the beam-b slot, jet g (3) in the jet slot
          call soft_ab_qgq(2, yy(2,3), yy(3,1), yy(1,2), Iq, 1, 2, 3, 4, 5, 6, s1, s2)
          call soft_nab_qgq(2, yy(2,3), yy(3,1), yy(1,2), Iq, 1, 2, 3, 4, 5, 6, s2n)
          Jb1 = Jq1/2; Jb2 = Jq2/4; Jc1 = Jg1; Jc2 = Jg2
       else
          ! gluon Born: jets q (2), qbar (3)
          call soft_ab_qag(2, yy(2,3), yy(3,1), yy(1,2), Ig, 1, 2, 3, 4, 5, 6, s1, s2)
          call soft_nab_qag(2, yy(2,3), yy(3,1), yy(1,2), Ig, 1, 2, 3, 4, 5, 6, s2n)
          call jetq(2, 2*E3, Jc1, Jc2, Q)
          Jb1 = Jq1/2; Jb2 = Jq2/4
       endif
       s2 = s2 + s2n
       c1(:,c) = 0; c2(:,c) = 0
       do it = 1, ntc
          tc = tcs(it)*Q
          do k = 1, 2
             hc = [hk(1,k,ht), hk(2,k,ht) + nfk(k,c)*gk(k,ht)]
             c1(it,c) = c1(it,c) + rk(k)*assemblejet(1, tc, ba0, 1.0_dp, ba1, Jb1, ba2, Jb2, s1, s2, Jc1, Jc2, hc)
             c2(it,c) = c2(it,c) + rk(k)*assemblejet(2, tc, ba0, 1.0_dp, ba1, Jb1, ba2, Jb2, s1, s2, Jc1, Jc2, hc)
          enddo
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
