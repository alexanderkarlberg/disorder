!-----------------------------------------------------------------------
! Slicing-adapted multichannel phase space for DIS n+1 (n = 3, 4 final
! partons) at fixed (x_B, Q^2, y), for nlo31 mode 1 (the above-cut part of
! tau_2-sliced NNLO 2+1), where the flat sequential-decay generator leaves
! the near-2+1 region (relative measure ~1e-5) to VEGAS and biases it low.
!
! Channels (Catani-Seymour phase-space factorisation, hep-ph/9605323 sec. 5,
! d = 4, measure dPhi_n = prod d^3p/((2pi)^3 2E) (2pi)^4 delta^4):
!   n = 3: flat (sequential decays, as nlo31's breit), or a 2+1 Born (log eta,
!          isotropic two-body) times one emission FF(i,j;k), FI(i,j;a) or
!          IF(a,i;k), all labels: 3 + 3 + 6 = 12 channels;
!   n = 4: flat, or the full n = 3 density times a second emission:
!          12 + 6 + 12 = 30 channels.
! Emission variables: FF: y, z_i, phi; FI: x, z_i, phi; IF: x, u, phi; with
! y, 1-x and both ends of z, u sampled logarithmically down to 1e-10
! (psmc_set_edge).
!   FF: dPhi_n = dPhi_{n-1} (2 pij.pk)/(16 pi^2) (1 - y) dy dz dphi/(2pi)
!   FI: deta dPhi_n(eta P) = (deta~/x) dx dPhi_{n-1}(eta~ P) (2 pij.pa)/(16 pi^2) dz dphi/(2pi)
!   IF: deta dPhi_n(eta P) = (deta~/x) dx dPhi_{n-1}(eta~ P) (2 pk.pa)/(16 pi^2) du dphi/(2pi)
! (pij, pk the mapped Born momenta, pa the incoming parton of the n+1 event).
! The returned weight is 1/g(Phi), g = sum_c alpha_c g_c(Phi) the density
! with respect to deta dPhi_n(eta P + q), every g_c from the exact inverse
! maps; nlo31 then uses jac*dphi = (y/x_B) wps.
!
! Momentum layout as nlo31's breit: P(:,1) incoming parton (Breit frame,
! along +z), P(:,2:n+1) partons, P(:,n+2) q, P(:,n+3), P(:,n+4) leptons
! (px, py, pz, E).
!-----------------------------------------------------------------------
module psmc
  implicit none
  private
  integer, parameter :: dp = kind(1.0d0), qp = selected_real_kind(30)
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
  ! lower edge of the log maps (y, 1-x) and of both ends of the logistic map
  ! (z, u): default 1e-10; psmc_set_edge changes it (5 Oct: technical-cut tests)
  real(dp), save :: vmin = 1e-10_dp, tz = 23.025850929940457_dp   ! ln(1e10)
  ! fraction of points from the flat generator
  real(dp), save, public :: psmc_aflat = 0.1_dp
  public :: psmc_set_edge
  real(dp), save :: xB = 0, Q = 0, yl = 0, Lx = 0
  ! channel lists: (type 1 FF / 2 FI / 3 IF, i, j or 0, k or 0)
  integer, save :: nch(3:4) = 0, ch(4,30,3:4)
  public :: psmc_init, psmc_gen, psmc_density
  public :: psmc_born, psmc_nchan
contains

  subroutine psmc_set_edge(e)
    real(dp), intent(in) :: e
    vmin = e; tz = log(1/e)
  end subroutine psmc_set_edge

  subroutine psmc_init(xB_in, Q2_in, y_in)
    real(dp), intent(in) :: xB_in, Q2_in, y_in
    integer :: n
    xB = xB_in; Q = sqrt(Q2_in); yl = y_in; Lx = log(1/xB)
    do n = 3, 4
       call channels(n)
    enddo
  end subroutine psmc_init

  subroutine channels(n)
    integer, intent(in) :: n
    integer :: i, j, k, m
    m = 0
    do i = 2, n + 1
       do j = i + 1, n + 1
          do k = 2, n + 1
             if (k == i .or. k == j) cycle
             m = m + 1; ch(:,m,n) = [1, i, j, k]
          enddo
          m = m + 1; ch(:,m,n) = [2, i, j, 0]
       enddo
    enddo
    do i = 2, n + 1
       do k = 2, n + 1
          if (k == i) cycle
          m = m + 1; ch(:,m,n) = [3, i, 0, k]
       enddo
    enddo
    nch(n) = m
  end subroutine channels

  ! leptons, q and the incoming parton for eta (Breit frame)
  subroutine frame(n, eta, P)
    integer, intent(in) :: n
    real(dp), intent(in) :: eta
    real(dp), intent(inout) :: P(4,n+4)
    real(dp) :: E
    E = eta*Q/(2*xB)
    P(:,1) = [0.0_dp, 0.0_dp, E, E]
    P(:,n+2) = [0.0_dp, 0.0_dp, -Q, 0.0_dp]
    P(:,n+3) = [Q/(2*yl)*2*sqrt(1 - yl), 0.0_dp, -Q/2, Q/(2*yl)*(2 - yl)]
    P(:,n+4) = [Q/(2*yl)*2*sqrt(1 - yl), 0.0_dp, Q/2, Q/(2*yl)*(2 - yl)]
  end subroutine frame

  real(dp) function eta_of(P) result(eta)
    real(dp), intent(in) :: P(4,*)
    eta = 2*xB*P(4,1)/Q
  end function eta_of

  !---------------------------------------------------------------- generation
  ! r(1): channel choice; then the variables. n = 3: 7 numbers, n = 4: 11.
  subroutine psmc_gen(n, r, P, eta, wps, ok)
    integer, intent(in) :: n
    real(dp), intent(in) :: r(:)
    real(dp), intent(out) :: P(4,n+4), eta, wps
    logical, intent(out) :: ok
    real(dp) :: g
    P = 0; ok = .false.; wps = 0; eta = 0
    if (r(1) < psmc_aflat) then
       call gen_flat(n, r(2:), P, ok)
    else
       call gen_cs(n, (r(1) - psmc_aflat)/(1 - psmc_aflat), r(2:), P, ok)
    endif
    if (.not. ok) return
    eta = eta_of(P)
    g = psmc_density(n, P)
    if (.not. (g > 0) .or. g /= g) then
       ok = .false.; return
    endif
    wps = 1/g
  end subroutine psmc_gen

  integer function psmc_nchan(n)
    integer, intent(in) :: n
    psmc_nchan = nch(n)
  end function psmc_nchan

  ! the 2+1 Born of the CS channels for given unit numbers r(1:3) (in
  ! psmc_gen(3, r) these are r(2:4)) and its density (same normalisation as
  ! psmc_density): correlated sampling, the below-cut 2+1 term at the Born of
  ! the above-cut 3+1 event (5 Oct)
  subroutine psmc_born(r, B, eta, wps, ok)
    real(dp), intent(in) :: r(3)
    real(dp), intent(out) :: B(4,6), eta, wps
    logical, intent(out) :: ok
    real(dp) :: g
    wps = 0; eta = 0
    call gen_born2(r, B, ok)
    if (.not. ok) return
    eta = eta_of(B)
    g = dens_born2(B)
    ok = g > 0 .and. g == g
    if (ok) wps = 1/g
  end subroutine psmc_born

  ! flat: eta (log) and sequential decays (nlo31's breit); n-1 + 2(n-1) + 1 numbers
  subroutine gen_flat(n, r, P, ok)
    integer, intent(in) :: n
    real(dp), intent(in) :: r(:)
    real(dp), intent(out) :: P(4,n+4)
    logical, intent(out) :: ok
    real(dp) :: eta, W2, m2(4), rest(4), pa(4), k(4,5), E, pz, beta, gam
    integer :: i, ir
    eta = xB*(1/xB)**r(1)
    call frame(n, eta, P)
    W2 = Q*Q*(eta/xB - 1)
    ok = W2 > 0
    if (.not. ok) return
    m2(1) = W2; ir = 1
    do i = 2, n - 1
       ir = ir + 1; m2(i) = m2(i-1)*r(ir)
    enddo
    rest = [0.0_dp, 0.0_dp, 0.0_dp, sqrt(W2)]
    do i = 1, n - 1
       if (i < n - 1) then
          call twobody(rest, m2(i+1), r(ir+1), r(ir+2), k(:,i), pa)
       else
          call twobody(rest, 0.0_dp, r(ir+1), r(ir+2), k(:,i), k(:,i+1))
       endif
       ir = ir + 2; rest = pa
    enddo
    E = P(4,1); pz = E - Q; beta = pz/E; gam = E/sqrt(W2)
    do i = 1, n
       P(:,i+1) = [k(1,i), k(2,i), gam*(k(3,i) + beta*k(4,i)), gam*(k(4,i) + beta*k(3,i))]
    enddo
  end subroutine gen_flat

  ! CS channels: n = 3 from a 2+1 Born, n = 4 from an n = 3 event (recursively)
  recursive subroutine gen_cs(n, rc, r, P, ok)
    integer, intent(in) :: n
    real(dp), intent(in) :: rc, r(:)
    real(dp), intent(out) :: P(4,n+4)
    logical, intent(out) :: ok
    real(dp) :: B(4,n+3), rsub
    integer :: ic, nr
    ok = .false.
    ic = min(nch(n), 1 + int(rc*nch(n)))
    ! the Born (n-1 final partons): its own numbers come first
    if (n == 3) then
       call gen_born2(r(1:3), B, ok); nr = 3
    else
       ! n - 1 = 3 from the full n = 3 density: flat or CS, chosen by r(1)
       rsub = r(1)
       if (rsub < psmc_aflat) then
          call gen_flat(3, r(2:7), B, ok)
       else
          call gen_cs(3, (rsub - psmc_aflat)/(1 - psmc_aflat), r(2:7), B, ok)
       endif
       nr = 7
    endif
    if (.not. ok) return
    call emit(n, ch(:,ic,n), B, r(nr+1:nr+3), P, ok)
  end subroutine gen_cs

  ! 2+1 Born: eta~ (log), isotropic two-body in the W rest frame
  subroutine gen_born2(r, B, ok)
    real(dp), intent(in) :: r(3)
    real(dp), intent(out) :: B(4,6)
    logical, intent(out) :: ok
    real(dp) :: eta, W2, rest(4), k1(4), k2(4), E, pz, beta, gam
    B = 0
    eta = xB*(1/xB)**r(1)
    call frame(2, eta, B)
    W2 = Q*Q*(eta/xB - 1)
    ok = W2 > 0
    if (.not. ok) return
    rest = [0.0_dp, 0.0_dp, 0.0_dp, sqrt(W2)]
    call twobody(rest, 0.0_dp, r(2), r(3), k1, k2)
    E = B(4,1); pz = E - Q; beta = pz/E; gam = E/sqrt(W2)
    B(:,2) = [k1(1), k1(2), gam*(k1(3) + beta*k1(4)), gam*(k1(4) + beta*k1(3))]
    B(:,3) = [k2(1), k2(2), gam*(k2(3) + beta*k2(4)), gam*(k2(4) + beta*k2(3))]
  end subroutine gen_born2

  ! one CS emission c = (type, i, j, k) from the Born B (n-1 final partons)
  ! into the n+1 event P; v(1:3) the emission's unit numbers
  subroutine emit(n, c, B, v, P, ok)
    integer, intent(in) :: n, c(4)
    real(dp), intent(in) :: B(4,n+3), v(3)
    real(dp), intent(out) :: P(4,n+4)
    logical, intent(out) :: ok
    integer :: slot(5), lab, m, i, j, k
    real(dp) :: a, eta
    ! in quad precision (9 Oct): a second emission off a nearly collinear
    ! pair (pt, pk) has transverse vectors with components ~ 1/theta, so in
    ! double precision the new momenta had m^2/E^2 up to 1e-7
    real(qp) :: pt(4), pk(4), pa(4), b1(4), z, x, u, yv, kt(4), e1(4), e2(4), phi, ktn
    ok = .false.
    i = c(2); j = c(3); k = c(4)
    ! Born slots 2..n hold the labels {2..n+1} minus j (FF, FI) or minus i (IF), ascending
    m = 1
    do lab = 2, n + 1
       if ((c(1) <= 2 .and. lab == j) .or. (c(1) == 3 .and. lab == i)) cycle
       m = m + 1; slot(m) = lab
    enddo
    P = 0
    phi = 2*pi*v(3)
    select case (c(1))
    case (1)   ! FF: same eta
       call frame(n, eta_of(B), P)
       call copy_others()
       pt = B(:, bslot(i)); pk = B(:, bslot(k))
       yv = lmap(v(1), vmin, 1.0_dp); z = zmap(v(2))
       call tbasis(pt, pk, e1, e2)
       ktn = sqrt(2*qdot(pt, pk)*yv*z*(1 - z))
       kt = ktn*(cos(phi)*e1 + sin(phi)*e2)
       P(:,i) = z*pt + (1 - z)*yv*pk + kt
       P(:,j) = (1 - z)*pt + z*yv*pk - kt
       P(:,k) = (1 - yv)*pk
    case (2)   ! FI: pa = pa~/x
       eta = eta_of(B)
       a = 1 - eta                          ! 1 - x <= 1 - eta~
       if (a <= vmin) return
       x = 1 - lmap(v(1), vmin, a); z = zmap(v(2))
       call frame(n, real(eta/x, dp), P)
       call copy_others()
       pt = B(:, bslot(i)); b1 = B(:,1); pa = b1/x
       call tbasis(pt, pa, e1, e2)
       ktn = sqrt(2*qdot(pt, b1)*(1 - x)/x*z*(1 - z))
       kt = ktn*(cos(phi)*e1 + sin(phi)*e2)
       P(:,i) = z*pt + (1 - z)*(1 - x)/x*b1 + kt
       P(:,j) = (1 - z)*pt + z*(1 - x)/x*b1 - kt
    case (3)   ! IF: pa = pa~/x
       eta = eta_of(B)
       a = 1 - eta
       if (a <= vmin) return
       x = 1 - lmap(v(1), vmin, a); u = zmap(v(2))
       call frame(n, real(eta/x, dp), P)
       call copy_others()
       pk = B(:, bslot(k)); b1 = B(:,1); pa = b1/x
       call tbasis(pk, pa, e1, e2)
       ktn = sqrt(2*qdot(pk, b1)*(1 - x)/x*u*(1 - u))
       kt = ktn*(cos(phi)*e1 + sin(phi)*e2)
       P(:,i) = u*pk + (1 - u)*(1 - x)/x*b1 + kt
       P(:,k) = (1 - u)*pk + u*(1 - x)/x*b1 - kt
    end select
    ok = .true.
  contains
    integer function bslot(lab)
      integer, intent(in) :: lab
      integer :: s
      bslot = 0
      do s = 2, n
         if (slot(s) == lab) bslot = s
      enddo
    end function bslot
    subroutine copy_others()
      integer :: s
      do s = 2, n
         if (slot(s) /= i .and. slot(s) /= k) P(:,slot(s)) = B(:,s)
      enddo
    end subroutine copy_others
  end subroutine emit

  !---------------------------------------------------------------- densities
  ! density g(Phi) with respect to deta dPhi_n(eta P + q), all channels
  recursive real(dp) function psmc_density(n, P) result(g)
    integer, intent(in) :: n
    real(dp), intent(in) :: P(4,n+4)
    real(dp) :: gc
    integer :: ic
    g = psmc_aflat*dens_flat(n, P)
    gc = 0
    do ic = 1, nch(n)
       gc = gc + dens_ch(n, ch(:,ic,n), P)
    enddo
    g = g + (1 - psmc_aflat)*gc/nch(n)
  end function psmc_density

  real(dp) function dens_flat(n, P) result(g)
    integer, intent(in) :: n
    real(dp), intent(in) :: P(4,n+4)
    real(dp) :: eta, W2, m2(4), dphi, s(4)
    integer :: i, k
    eta = eta_of(P)
    W2 = Q*Q*(eta/xB - 1)
    m2(1) = W2
    ! m2(i) = (sum of partons i+1 .. n)^2 (the remainders of the decay chain)
    do i = 2, n - 1
       s = 0
       do k = i + 1, n + 1
          s = s + P(:,k)
       enddo
       m2(i) = mdot(s, s)
    enddo
    dphi = 1
    do i = 2, n - 1
       dphi = dphi*m2(i-1)/(2*pi)
    enddo
    do i = 1, n - 1
       if (i < n - 1) then
          dphi = dphi*(1 - m2(i+1)/m2(i))/(8*pi)
       else
          dphi = dphi/(8*pi)
       endif
    enddo
    g = 1/(eta*Lx*dphi)
  end function dens_flat

  ! density of the 2+1 Born B (eta~ log, isotropic two-body: dPhi_2 = 1/(8 pi) per unit)
  real(dp) function dens_born2(B) result(g)
    real(dp), intent(in) :: B(4,6)
    g = 8*pi/(eta_of(B)*Lx)
  end function dens_born2

  ! density of channel c at the n+1 event P (inverse map, Born density times
  ! the emission density over the CS factor)
  recursive real(dp) function dens_ch(n, c, P) result(g)
    integer, intent(in) :: n, c(4)
    real(dp), intent(in) :: P(4,n+4)
    real(dp) :: B(4,n+3), pij, pik, pjk, pia, pja, pka, yv, z, x, u, gB, a, fac
    integer :: i, j, k, lab, m, slot(5)
    g = 0
    i = c(2); j = c(3); k = c(4)
    m = 1
    do lab = 2, n + 1
       if ((c(1) <= 2 .and. lab == j) .or. (c(1) == 3 .and. lab == i)) cycle
       m = m + 1; slot(m) = lab
    enddo
    B = 0
    select case (c(1))
    case (1)
       pij = mdot(P(:,i), P(:,j)); pik = mdot(P(:,i), P(:,k)); pjk = mdot(P(:,j), P(:,k))
       yv = pij/(pij + pik + pjk); z = pik/(pik + pjk)
       if (yv < vmin .or. yv >= 1) return
       call frame(n - 1, eta_of(P), B)
       do m = 2, n
          lab = slot(m)
          if (lab == i) then
             B(:,m) = P(:,i) + P(:,j) - yv/(1 - yv)*P(:,k)
          elseif (lab == k) then
             B(:,m) = P(:,k)/(1 - yv)
          else
             B(:,m) = P(:,lab)
          endif
       enddo
       fac = 2*mdot(B(:,ms(i)), B(:,ms(k)))/(16*pi**2)*(1 - yv)
       g = born_dens(n - 1, B)*hl(yv, vmin, 1.0_dp)*hz(z)/(2*pi)/fac*(2*pi)
    case (2)
       pij = mdot(P(:,i), P(:,j)); pia = mdot(P(:,i), P(:,1)); pja = mdot(P(:,j), P(:,1))
       x = (pia + pja - pij)/(pia + pja); z = pia/(pia + pja)
       call frame(n - 1, x*eta_of(P), B)
       a = 1 - eta_of(B)
       if (1 - x < vmin .or. 1 - x > a) return
       do m = 2, n
          lab = slot(m)
          if (lab == i) then
             B(:,m) = P(:,i) + P(:,j) - (1 - x)*P(:,1)
          else
             B(:,m) = P(:,lab)
          endif
       enddo
       fac = 2*mdot(B(:,ms(i)), P(:,1))/(16*pi**2)/x
       g = born_dens(n - 1, B)*hl(1 - x, vmin, a)*hz(z)/fac
    case (3)
       pia = mdot(P(:,i), P(:,1)); pka = mdot(P(:,k), P(:,1)); pik = mdot(P(:,i), P(:,k))
       x = (pka + pia - pik)/(pka + pia); u = pia/(pia + pka)
       call frame(n - 1, x*eta_of(P), B)
       a = 1 - eta_of(B)
       if (1 - x < vmin .or. 1 - x > a) return
       do m = 2, n
          lab = slot(m)
          if (lab == k) then
             B(:,m) = P(:,k) + P(:,i) - (1 - x)*P(:,1)
          else
             B(:,m) = P(:,lab)
          endif
       enddo
       fac = 2*mdot(B(:,ms(k)), P(:,1))/(16*pi**2)/x
       g = born_dens(n - 1, B)*hl(1 - x, vmin, a)*hz(u)/fac
    end select
    if (.not. (g > 0)) g = 0
  contains
    integer function ms(lab)
      integer, intent(in) :: lab
      integer :: s
      ms = 0
      do s = 2, n
         if (slot(s) == lab) ms = s
      enddo
    end function ms
  end function dens_ch

  recursive real(dp) function born_dens(nb, B) result(g)
    integer, intent(in) :: nb
    real(dp), intent(in) :: B(4,nb+4)
    if (nb == 2) then
       g = dens_born2(B)
    else
       g = psmc_density(nb, B)
    endif
  end function born_dens

  !---------------------------------------------------------------- helpers
  ! log map on [lo, hi] and its density
  real(dp) function lmap(r, lo, hi)
    real(dp), intent(in) :: r, lo, hi
    lmap = lo*(hi/lo)**r
  end function lmap
  real(dp) function hl(v, lo, hi)
    real(dp), intent(in) :: v, lo, hi
    hl = 1/(v*log(hi/lo))
  end function hl
  ! logistic map of z, both ends logarithmic down to e^-tz, and its density
  real(dp) function zmap(r)
    real(dp), intent(in) :: r
    zmap = 1/(1 + exp(-tz*(2*r - 1)))
  end function zmap
  real(dp) function hz(z)
    real(dp), intent(in) :: z
    hz = 0
    if (z > 0 .and. z < 1) then
       if (abs(log(z/(1 - z))) <= tz) hz = 1/(z*(1 - z)*2*tz)
    endif
  end function hz

  ! unit spacelike e1, e2 orthogonal to the light-like a, b (and to each other)
  ! (9 Oct: e1 from eps as e2, the largest of the three axes; the projection
  ! rr - (rr.b) a/ab - (rr.a) b/ab lost ~1/theta_ab^2 in the orthogonality,
  ! so that a second emission off a nearly collinear pair had m^2/E^2 up to
  ! 1e-6, negative p_i.p_j and negative-energy mapped dipole Borns)
  subroutine tbasis(a, b, e1, e2)
    real(qp), intent(in) :: a(4), b(4)
    real(qp), intent(out) :: e1(4), e2(4)
    real(qp) :: rr(4), et(4), n2, nt
    integer :: t
    n2 = -1
    do t = 1, 3
       rr = 0; rr(t) = 1
       et = eps(a, b, rr)
       nt = -qdot(et, et)
       if (nt > n2) then
          e1 = et; n2 = nt
       endif
    enddo
    e1 = e1/sqrt(n2)
    e2 = eps(a, b, e1)
    e2 = e2/sqrt(-qdot(e2, e2))
  end subroutine tbasis

  ! e^mu = eps^{mu nu rho sigma} a_nu b_rho c_sigma with lowered a, b, c
  ! (components (x, y, z, t), metric diag(-1,-1,-1,1)): orthogonal to a, b, c
  function eps(a, b, c) result(e)
    real(qp), intent(in) :: a(4), b(4), c(4)
    real(qp) :: e(4), al(4), bl(4), cl(4), s
    integer :: mu, nu, ro, si
    al = [-a(1:3), a(4)]; bl = [-b(1:3), b(4)]; cl = [-c(1:3), c(4)]
    e = 0
    do mu = 1, 4
       do nu = 1, 4
          do ro = 1, 4
             do si = 1, 4
                s = perm_sign([mu, nu, ro, si])
                if (s /= 0) e(mu) = e(mu) + s*al(nu)*bl(ro)*cl(si)
             enddo
          enddo
       enddo
    enddo
  end function eps

  real(dp) function perm_sign(p)
    integer, intent(in) :: p(4)
    integer :: i, j
    perm_sign = 1
    do i = 1, 4
       do j = i + 1, 4
          if (p(i) == p(j)) then
             perm_sign = 0; return
          endif
          if (p(i) > p(j)) perm_sign = -perm_sign
       enddo
    enddo
  end function perm_sign

  subroutine twobody(p, m2b, rc, rphi, pa, pb)
    real(dp), intent(in) :: p(4), m2b, rc, rphi
    real(dp), intent(out) :: pa(4), pb(4)
    real(dp) :: M, pmod, c, st, phi, d(3), qa(4), bp, f
    M = sqrt(max(p(4)**2 - sum(p(1:3)**2), 0.0_dp))
    pmod = (M**2 - m2b)/(2*M)
    c = 2*rc - 1; st = sqrt(max(1 - c*c, 0.0_dp)); phi = 2*pi*rphi
    d = [st*cos(phi), st*sin(phi), c]
    qa = [pmod*d, pmod]
    bp = dot_product(qa(1:3), p(1:3))
    pa(4) = (p(4)*qa(4) + bp)/M
    f = (bp/(p(4) + M) + qa(4))/M
    pa(1:3) = qa(1:3) + f*p(1:3)
    pb = p - pa
  end subroutine twobody

  pure real(dp) function mdot(a, b)
    real(dp), intent(in) :: a(4), b(4)
    mdot = a(4)*b(4) - a(1)*b(1) - a(2)*b(2) - a(3)*b(3)
  end function mdot
  pure real(qp) function qdot(a, b)
    real(qp), intent(in) :: a(4), b(4)
    qdot = a(4)*b(4) - a(1)*b(1) - a(2)*b(2) - a(3)*b(3)
  end function qdot
end module psmc
