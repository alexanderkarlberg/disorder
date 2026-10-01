!----------------------------------------------------------------------
! One-loop SCET ingredients for N-jettiness slicing in DIS (photon
! exchange, massless quarks, mu_R = mu_F = mu = Q), in units of
! alpha_s/(2 pi). Conventions as in Gaunt, Stahlhofen, Tackmann, Walsh,
! arXiv:1505.04794 (GSTW): geometric measure, Q_i = 2 E_i (Breit frame),
!   T_N = sum_k min_i n_i.p_k ,   n_i = (1, nhat_i),
! and with xi = T_cut the slicing cumulant at NLO is C_{-1}^(1)(xi):
!   hard + sum_jets J_{-1}(lambda_J) + B_{-1}(lambda_B) + S_{-1}(lambda_S),
!   lambda_J = 2 E_J T_cut/mu^2, lambda_B = 2 E_a T_cut/mu^2,
!   lambda_S = T_cut/mu.
! Jet: GSTW (A.8)-(A.12); beam: GSTW (A.19) with I^(1)_qq, I^(1)_qg of
! arXiv:1401.5478 and I^(1)_gg, I^(1)_gq of arXiv:1405.1044 (A.4);
! soft: GSTW (A.31)-(A.32) with the non-hemisphere integrals I0, I1 of
! Jouttenus, Stewart, Tackmann, Waalewijn, arXiv:1102.4344, Eqs. (56)-(60).
! The hard function is built from DISENT's factorising virtual + I
! operator (VIRTHR's QQ, GG), the finite part of the Catani-Seymour I
! operator and the conversion to MS-bar (see slicing/README.md).
!----------------------------------------------------------------------
module mod_slicing_scet
  use types, only: dp
  implicit none
  private
  public :: cli2, soft_I0I1, soft_cum, soft_geom, soft_from_geom, jet_cum, beam_coeffs, hard_fact
  public :: scet_set_nf, pi, zeta2

  real(dp), parameter :: pi = 3.141592653589793238462643383279502884_dp
  real(dp), parameter :: zeta2 = pi**2 / 6
  real(dp), parameter :: CF = 4.0_dp/3.0_dp, CA = 3.0_dp, TF = 0.5_dp
  integer, save :: nf = 5
  real(dp), public, save :: soft_tol = 1e-11_dp   ! GK tolerance of soft_I0I1
  ! diagnostics: 0 = all PDFs, 1 = quarks only (gluon PDF set to zero), 2 = gluon only
  integer, public, save :: pdf_mask = 0
  public :: mask_pdf

  ! Gauss-Legendre nodes/weights on [0,1] for the beam convolutions
  integer, parameter :: ngl = 64
  real(dp), save :: glx(ngl), glw(ngl)
  logical, save :: gl_init = .false.

  interface
     subroutine EvolvePDF(x, Q, res)
       import dp
       real(dp), intent(in) :: x, Q
       real(dp), intent(out) :: res(-6:6)
     end subroutine EvolvePDF
  end interface

contains

  subroutine mask_pdf(f)
    real(dp), intent(inout) :: f(-6:6)
    if (pdf_mask == 1) then
       f(0) = 0
    elseif (pdf_mask == 2) then
       f(1:6) = 0; f(-6:-1) = 0
    endif
  end subroutine mask_pdf

  subroutine scet_set_nf(n)
    integer, intent(in) :: n
    nf = n
  end subroutine scet_set_nf

  real(dp) function beta0()
    beta0 = 11.0_dp/3 * CA - 4.0_dp/3 * TF * nf
  end function beta0

  !--------------------------------------------------------------------
  ! Complex dilogarithm Li2(z) (Bernoulli series after mapping to
  ! |z| <= 1, Re z <= 1/2).
  recursive complex(dp) function cli2(z) result(res)
    complex(dp), intent(in) :: z
    complex(dp) :: u, u2, term
    real(dp), parameter :: b(10) = [ -0.25_dp, 1.0_dp/36, -1.0_dp/3600, &
         & 1.0_dp/211680, -1.0_dp/10886400, 1.0_dp/526901760, &
         & -4.064761645144225524e-11_dp, 8.921691020456452555e-13_dp, &
         & -1.993929586072107568e-14_dp, 4.518980029619918192e-16_dp ]
    integer :: k
    if (abs(z) < 1e-300_dp) then
       res = (0.0_dp, 0.0_dp); return
    endif
    if (abs(z - (1.0_dp, 0.0_dp)) < 1e-15_dp) then
       res = cmplx(zeta2, 0.0_dp, dp); return
    endif
    if (abs(z) > 1.0_dp) then
       res = -cli2(1 / z) - zeta2 - 0.5_dp * log(-z)**2
       return
    endif
    if (real(z) > 0.5_dp) then
       res = -cli2(1 - z) + zeta2 - log(z) * log(1 - z)
       return
    endif
    ! Li2(z) = sum_n B_n u^(n+1)/(n+1)!, u = -log(1-z), in the form
    ! u - u^2/4 + sum_k B_2k u^(2k+1)/(2k+1)!
    u = -log(1 - z)
    u2 = u * u
    res = u + b(1) * u2
    term = u
    do k = 2, 10
       term = term * u2
       res = res + b(k) * term
    enddo
  end function cli2

  !--------------------------------------------------------------------
  ! Non-hemisphere integrals of JSTW (60), as one-dimensional integrals
  ! over phi (inner y integral analytic):
  !   I0 = (2/pi) int_0^pi dphi ln(y_hi/y_lo),
  !   I1 = (2/pi) int_0^pi dphi (-2) Re[Li2(y_hi e^{i phi}) - Li2(y_lo e^{i phi})],
  ! over y > sqrt(beta/alpha) and 1 + y^2 - 2 y cos(phi) < 1/alpha.
  subroutine soft_I0I1(alpha, beta, I0, I1)
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(out) :: I0, I1
    real(dp) :: a, b, e0, e1
    a = 0; b = pi
    call adapt(a, b, I0, I1)
    I0 = 2 / pi * I0
    I1 = 2 / pi * I1
  contains
    subroutine inner(phi, f0, f1)
      real(dp), intent(in) :: phi
      real(dp), intent(out) :: f0, f1
      real(dp) :: disc, ym, yp, y0, lo, hi
      complex(dp) :: e
      f0 = 0; f1 = 0
      disc = 1 / alpha - sin(phi)**2
      if (disc <= 0) return
      ym = cos(phi) - sqrt(disc); yp = cos(phi) + sqrt(disc)
      y0 = sqrt(beta / alpha)
      lo = max(ym, y0); hi = yp
      if (hi <= lo) return
      f0 = log(hi / lo)
      e = cmplx(cos(phi), sin(phi), dp)
      f1 = -2 * real(cli2(hi * e) - cli2(lo * e))
    end subroutine inner
    ! adaptive 15-point Gauss-Kronrod on [a,b] for both integrands
    recursive subroutine adapt_rec(a, b, r0, r1, depth)
      real(dp), intent(in) :: a, b
      real(dp), intent(out) :: r0, r1
      integer, intent(in) :: depth
      real(dp) :: g0, g1, k0, k1, m, l0, l1, u0, u1
      call gk15(a, b, g0, g1, k0, k1)
      if (depth > 40 .or. (abs(k0 - g0) < soft_tol * max(1.0_dp, abs(k0)) .and. &
           & abs(k1 - g1) < soft_tol * max(1.0_dp, abs(k1)))) then
         r0 = k0; r1 = k1; return
      endif
      m = 0.5_dp * (a + b)
      call adapt_rec(a, m, l0, l1, depth + 1)
      call adapt_rec(m, b, u0, u1, depth + 1)
      r0 = l0 + u0; r1 = l1 + u1
    end subroutine adapt_rec
    subroutine adapt(a, b, r0, r1)
      real(dp), intent(in) :: a, b
      real(dp), intent(out) :: r0, r1
      call adapt_rec(a, b, r0, r1, 0)
    end subroutine adapt
    subroutine gk15(a, b, g0, g1, k0, k1)
      real(dp), intent(in) :: a, b
      real(dp), intent(out) :: g0, g1, k0, k1
      real(dp), parameter :: xk(8) = [0.991455371120812639_dp, 0.949107912342758525_dp, &
           & 0.864864423359769073_dp, 0.741531185599394440_dp, 0.586087235467691130_dp, &
           & 0.405845151377397167_dp, 0.207784955007898468_dp, 0.0_dp]
      real(dp), parameter :: wk(8) = [0.022935322010529225_dp, 0.063092092629978553_dp, &
           & 0.104790010322250184_dp, 0.140653259715525919_dp, 0.169004726639267903_dp, &
           & 0.190350578064785410_dp, 0.204432940075298892_dp, 0.209482141084727828_dp]
      real(dp), parameter :: wg(4) = [0.129484966168869693_dp, 0.279705391489276668_dp, &
           & 0.381830050505118945_dp, 0.417959183673469388_dp]
      real(dp) :: c, h, f0a, f1a, f0b, f1b
      integer :: j
      c = 0.5_dp * (a + b); h = 0.5_dp * (b - a)
      call inner(c, f0a, f1a)
      k0 = wk(8) * f0a; k1 = wk(8) * f1a
      g0 = wg(4) * f0a; g1 = wg(4) * f1a
      do j = 1, 7
         call inner(c - h * xk(j), f0a, f1a)
         call inner(c + h * xk(j), f0b, f1b)
         k0 = k0 + wk(j) * (f0a + f0b); k1 = k1 + wk(j) * (f1a + f1b)
         if (mod(j, 2) == 0) then
            g0 = g0 + wg(j/2) * (f0a + f0b); g1 = g1 + wg(j/2) * (f1a + f1b)
         endif
      enddo
      k0 = h * k0; k1 = h * k1; g0 = h * g0; g1 = h * g1
    end subroutine gk15
  end subroutine soft_I0I1

  !--------------------------------------------------------------------
  ! One-loop soft cumulant (alpha_s/2pi units) for the N-jettiness soft
  ! function with directions nhat(3,n), Casimirs cas(n), colour
  ! correlators tt(n,n) (T_i.T_j), at lambda = T_cut/mu. GSTW (A.24),
  ! (A.31), (A.32); s_ij = n_i.n_j/2 (geometric measure, Q_i = 2 E_i).
  real(dp) function soft_cum(n, nhat, cas, tt, lam) result(res)
    integer, intent(in) :: n
    real(dp), intent(in) :: nhat(3,n), cas(n), tt(n,n), lam
    real(dp) :: sij(n,n), s1, s0, sm1, I0, I1, L
    integer :: i, j, m
    do i = 1, n
       do j = 1, n
          sij(i,j) = 0.5_dp * (1 - dot_product(nhat(:,i), nhat(:,j)))
       enddo
    enddo
    s1 = -8 * sum(cas)
    s0 = 0; sm1 = 0
    do i = 1, n
       do j = 1, n
          if (i == j) cycle
          s0 = s0 - 4 * tt(i,j) * log(sij(i,j))
          sm1 = sm1 + tt(i,j) * (log(sij(i,j))**2 - zeta2)
          do m = 1, n
             if (m == i .or. m == j) cycle
             call soft_I0I1(sij(j,m) / sij(i,j), sij(i,m) / sij(i,j), I0, I1)
             sm1 = sm1 + tt(i,j) * 4 * (I0 * log(sij(j,m) / sij(i,j)) + I1)
          enddo
       enddo
    enddo
    L = log(lam)
    res = 0.5_dp * (sm1 + s0 * L + s1 * L * L / 2)
  end function soft_cum

  !--------------------------------------------------------------------
  ! Geometry-only part of the soft function for n directions (shared by
  ! all colour channels): for ordered pairs (i,j), g(i,j) = ln^2 s_ij -
  ! zeta2 + 4 sum_m I_ij,m and ls(i,j) = ln s_ij.
  subroutine soft_geom(n, nhat, g, ls)
    integer, intent(in) :: n
    real(dp), intent(in) :: nhat(3,n)
    real(dp), intent(out) :: g(n,n), ls(n,n)
    real(dp) :: sij(n,n), I0, I1
    integer :: i, j, m
    do i = 1, n
       do j = 1, n
          sij(i,j) = 0.5_dp * (1 - dot_product(nhat(:,i), nhat(:,j)))
       enddo
    enddo
    g = 0; ls = 0
    do i = 1, n
       do j = 1, n
          if (i == j) cycle
          ls(i,j) = log(sij(i,j))
          g(i,j) = ls(i,j)**2 - zeta2
          do m = 1, n
             if (m == i .or. m == j) cycle
             call soft_I0I1(sij(j,m) / sij(i,j), sij(i,m) / sij(i,j), I0, I1)
             g(i,j) = g(i,j) + 4 * (I0 * log(sij(j,m) / sij(i,j)) + I1)
          enddo
       enddo
    enddo
  end subroutine soft_geom

  ! soft cumulant (alpha_s/2pi) from the geometry part, for Casimirs cas
  ! and colour correlators tt, at lambda = T_cut/mu
  real(dp) function soft_from_geom(n, g, ls, cas, tt, lam) result(res)
    integer, intent(in) :: n
    real(dp), intent(in) :: g(n,n), ls(n,n), cas(n), tt(n,n), lam
    real(dp) :: L
    L = log(lam)
    res = 0.5_dp * (sum(tt * g) - 4 * sum(tt * ls) * L - 8 * sum(cas) * L * L / 2)
  end function soft_from_geom

  !--------------------------------------------------------------------
  ! One-loop jet cumulant J_{-1}(lambda) in alpha_s/2pi units.
  real(dp) function jet_cum(gluon, lam) result(res)
    logical, intent(in) :: gluon
    real(dp), intent(in) :: lam
    real(dp) :: L
    L = log(lam)
    if (gluon) then
       res = 0.5_dp * (CA * (4.0_dp/3 - pi**2) + 5.0_dp/3 * beta0() - beta0() * L + 2 * CA * L * L)
    else
       res = 0.5_dp * (CF * (7 - pi**2) - 3 * CF * L + 2 * CF * L * L)
    endif
  end function jet_cum

  !--------------------------------------------------------------------
  ! Beam-function coefficients at x = eta, scale Q (mu = mu_F = Q), in
  ! alpha_s/2pi units and x f(x) form: for each incoming flavour I
  ! (0 = gluon) B(I) = c0(I) + c1(I) ln(lambda_B) + c2(I) ln^2(lambda_B),
  ! to be used instead of x f_I(eta) at O(alpha_s).
  subroutine beam_coeffs(eta, Q, c0, c1, c2)
    real(dp), intent(in) :: eta, Q
    real(dp), intent(out) :: c0(-6:6), c1(-6:6), c2(-6:6)
    real(dp) :: f1(-6:6), fz(-6:6), z, w, dz, sq, l1m, lz, pqg, pgq, pggr, r1q, r1g
    real(dp) :: hq1(-6:6), hg1, aq0(-6:6), aq1(-6:6), ag0, ag1, sumqz
    integer :: k, i
    call init_gl()
    call EvolvePDF(eta, Q, f1)
    call mask_pdf(f1)
    c0 = 0; c1 = 0; c2 = 0
    aq0 = 0; aq1 = 0; ag0 = 0; ag1 = 0
    ! plus-distribution endpoint pieces: h(z) = r(z) xf(eta/z), r(1) = 2 for
    ! both (1+z^2) and 2(1-z+z^2)^2/z
    hq1 = 2 * f1
    hg1 = 2 * f1(0)
    do k = 1, ngl
       w = glx(k)
       z = 1 - (1 - eta) * w * w
       dz = 2 * (1 - eta) * w * glw(k)
       if (z <= eta) cycle
       call EvolvePDF(eta / z, Q, fz)
       call mask_pdf(fz)
       sumqz = 0
       do i = 1, nf
          sumqz = sumqz + fz(i) + fz(-i)
       enddo
       l1m = log(1 - z); lz = log(z)
       pqg = z * z + (1 - z)**2
       pgq = (1 + (1 - z)**2) / z
       r1q = 1 + z * z
       r1g = 2 * (1 - z + z * z)**2 / z
       pggr = r1g / (1 - z) * lz        ! P_gg(z) ln z (regular)
       ! quark: CF { (1+z^2)[ln(1-z)/(1-z)]_+ + (1-z) - (1+z^2) ln z/(1-z) } + TF { Pqg ln((1-z)/z) + 2z(1-z) }
       !        ln(lambda): CF (1+z^2)[1/(1-z)]_+ + TF Pqg
       do i = -nf, nf
          if (i == 0) cycle
          aq0(i) = aq0(i) + dz * (CF * (l1m / (1 - z) * (r1q * fz(i) - hq1(i)) &
               & + (1 - z - r1q / (1 - z) * lz) * fz(i)) &
               & + TF * (pqg * (l1m - lz) + 2 * z * (1 - z)) * fz(0))
          aq1(i) = aq1(i) + dz * (CF * (r1q * fz(i) - hq1(i)) / (1 - z) + TF * pqg * fz(0))
       enddo
       ! gluon: CA { 2(1-z+z^2)^2/z [ln(1-z)/(1-z)]_+ - Pgg ln z } + CF { Pgq ln((1-z)/z) + z } sum_q
       !        ln(lambda): CA 2(1-z+z^2)^2/z [1/(1-z)]_+ + CF Pgq sum_q
       ag0 = ag0 + dz * (CA * (l1m / (1 - z) * (r1g * fz(0) - hg1) - pggr * fz(0)) &
            & + CF * (pgq * (l1m - lz) + z) * sumqz)
       ag1 = ag1 + dz * (CA * (r1g * fz(0) - hg1) / (1 - z) + CF * pgq * sumqz)
    enddo
    ! endpoint pieces -h(1) int_0^eta [..] and delta terms
    ! [g]_+ endpoint: - h(1) int_0^eta g, with int_0^eta ln(1-z)/(1-z) = -ln^2(1-eta)/2
    ! and int_0^eta 1/(1-z) = -ln(1-eta); times the colour factor of the plus term
    do i = -nf, nf
       if (i == 0) cycle
       c0(i) = aq0(i) + CF * hq1(i) * 0.5_dp * log(1 - eta)**2 - CF * zeta2 * f1(i)
       c1(i) = aq1(i) + CF * hq1(i) * log(1 - eta)
       c2(i) = CF * f1(i)
    enddo
    c0(0) = ag0 + CA * hg1 * 0.5_dp * log(1 - eta)**2 - CA * zeta2 * f1(0)
    c1(0) = ag1 + CA * hg1 * log(1 - eta)
    c2(0) = CA * f1(0)
  end subroutine beam_coeffs

  subroutine init_gl()
    integer :: i, j
    real(dp) :: z, z1, p1, p2, p3, pp
    if (gl_init) return
    do i = 1, (ngl + 1) / 2
       z = cos(pi * (i - 0.25_dp) / (ngl + 0.5_dp))
       do
          p1 = 1; p2 = 0
          do j = 1, ngl
             p3 = p2; p2 = p1
             p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j
          enddo
          pp = ngl * (z * p1 - p2) / (z * z - 1)
          z1 = z; z = z1 - p1 / pp
          if (abs(z - z1) < 1e-15_dp) exit
       enddo
       glx(i) = 0.5_dp * (1 - z); glx(ngl + 1 - i) = 0.5_dp * (1 + z)
       glw(i) = 1 / ((1 - z * z) * pp * pp); glw(ngl + 1 - i) = glw(i)
    enddo
    gl_init = .true.
  end subroutine init_gl

  !--------------------------------------------------------------------
  ! Factorising part of the one-loop hard function (alpha_s/2pi units,
  ! relative to the Born, mu = Q) for the three-parton DIS Born, from
  ! DISENT's VIRTHR constants: H = QQ (or GG) + I0_CS + (pi^2/12) sum C_i,
  ! where I0_CS is the finite part of the Catani-Seymour I operator,
  !   sum_i [C_i pi^2/3 - gamma_i - K_i]
  ! + sum_{i/=k} T_i.T_k [ l_ik^2/2 - (gamma_i/C_i) l_ik ],  l_ik = ln(2 p_i.p_k/Q^2).
  ! Checked on the two-parton case: VIRTWO's QQ = CF(2 - pi^2) gives
  ! CF(-8 + pi^2/6), the DIS quark form factor at mu = Q.
  ! gluon = .false.: partons (1,2,3) = (q_in, q_out, g); .true.: (g_in, q, qbar).
  real(dp) function hard_fact(gluon, l12, l13, l23) result(res)
    logical, intent(in) :: gluon
    real(dp), intent(in) :: l12, l13, l23
    real(dp) :: cas(3), gam(3), kk(3), t12, t13, t23, gq, gg, kq, kg, qqv
    gq = 1.5_dp * CF
    gg = 11.0_dp/6 * CA - 2.0_dp/3 * TF * nf
    kq = (3.5_dp - zeta2) * CF
    kg = (67.0_dp/18 - zeta2) * CA - 10.0_dp/9 * TF * nf
    if (.not. gluon) then
       cas = [CF, CF, CA]; gam = [gq, gq, gg]; kk = [kq, kq, kg]
       qqv = CF*2 + CA*50.0_dp/9 - TF*nf*16.0_dp/9 - CF*pi**2 &
            & - 3*(CF - CA/2)*l12 - (5*CA - TF*nf)/3*(l13 + l23)
    else
       cas = [CA, CF, CF]; gam = [gg, gq, gq]; kk = [kg, kq, kq]
       qqv = CF*2 + CA*50.0_dp/9 - TF*nf*16.0_dp/9 - CA*pi**2 &
            & - 3*(CF - CA/2)*l23 - (5*CA - TF*nf)/3*(l12 + l13)
    endif
    t12 = 0.5_dp * (cas(3) - cas(1) - cas(2))
    t13 = 0.5_dp * (cas(2) - cas(1) - cas(3))
    t23 = 0.5_dp * (cas(1) - cas(2) - cas(3))
    res = qqv + sum(cas * pi**2 / 3 - gam - kk) &
         & + t12 * (l12**2 - (gam(1)/cas(1) + gam(2)/cas(2)) * l12) &
         & + t13 * (l13**2 - (gam(1)/cas(1) + gam(3)/cas(3)) * l13) &
         & + t23 * (l23**2 - (gam(2)/cas(2) + gam(3)/cas(3)) * l23) &
         & + pi**2 / 12 * sum(cas)
  end function hard_fact

end module mod_slicing_scet
