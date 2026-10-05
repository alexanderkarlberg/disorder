!-----------------------------------------------------------------------
! NLO DIS 3+1 (photon exchange) with Catani-Seymour subtraction: a first
! integrator for the validation of the dis31 ingredients.
!
! sigma(>= 3 jets) at O(alpha_s^2) (LO) and O(alpha_s^3) (NLO correction),
! e p -> e + >= 3 jets, inclusive kt jets (R = 1, E-scheme) in the Breit
! frame with p_T > ptmin; cuts on Q^2 and y. Parts:
!   lo : |M_{3+1}|^2
!   vi : 2 Re <M0|M1> + <I>   (virt31_ren + iop31_i, finite parts)
!   kp : K + P (iop31_kp, iop31_kernel; x-convolution with the PDFs)
!   r  : |M_{4+1}|^2 - sum of dipoles (me41, dip41), each with its own jet
!        function
! Phase space in the Breit frame: Q^2 (log), y, eta (log, x_B < eta < 1),
! the hadronic system by sequential two-body decays; mu_R = mu_F = Q.
!
! Mode 1 (nnlo21): the above-cut part of tau_2-sliced NNLO DIS 2+1 at fixed
! (x, Q^2): dsigma/dx dQ^2 [pb/GeV^2] with the jet function replaced by
! theta(tau_2 > tau_cut) for ntc values of tau_cut (tau_2 = T_2/Q in the
! jets' rest frame, slicing/mod_tau2_run.f90 measure 2) in bins of tau_zQ
! (as the NLO slicing tests); lo then is the O(alpha_s^2) part of NLO 2+1
! above the cut, vi + kp + r the O(alpha_s^3) part of NNLO 2+1 above it.
!
! Usage: nlo31 part ncall itmx [seed [techcut [mode x Q2 [logmap|uniform [pdfmask]]]]]
!        (part = lo, vi, kp, r; pdfmask 0 all, 1 quarks only, 2 gluon only;
!        9th argument: uniform | logmap | psmc (slicing-adapted multichannel, dis31/psmc.f90))
!        optional 11th: vslice k (VEGAS adapts to tcs(k) < tau_2 < tcs(k-1))
!-----------------------------------------------------------------------
module nlo31_mod
  use psmc
  use born31
  use me41
  use dip41
  use virt31
  use iop31
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: gev2pb = 0.3893793721e9_dp
  real(dp), parameter :: CF = 4.0_dp/3, CA = 3
  ! set-up (as NNLOJET's epLJJ run for the comparison)
  real(dp), parameter :: Ee = 27.5_dp, Ep = 920.0_dp, s = 4*Ee*Ep
  real(dp), parameter :: q2min = 150, q2max = 15000, ymin = 0.1_dp, ymax = 0.9_dp
  real(dp), parameter :: ptmin = 5, rjet = 1
  integer, parameter :: njmin = 3
  real(dp) :: techcut = 1d-9
  character(8) :: part
  ! flavour lists with symmetry factors
  integer :: nb, nr
  integer :: flb(4,200), flr(5,400)
  real(dp) :: symb(200), symr(400)
  ! histogram: Q^2 bins
  integer, parameter :: nq = 6
  real(dp), parameter :: qedge(0:nq) = [150.0_dp, 200.0_dp, 300.0_dp, 500.0_dp, 1000.0_dp, 3000.0_dp, 15000.0_dp]
  real(dp) :: hist(nq), hist2(nq), hacc(nq)
  ! mode 1: fixed (x, Q^2), cells (tau_cut k, tau_zQ bin b) -> k + ntc*(b - 1)
  integer :: mode = 0
  real(dp) :: xfix = 0.01_dp, Q2fix = 400
  integer, parameter :: ntc = 10, nzb = 6, ncell = ntc*nzb
  real(dp), parameter :: tcs(ntc) = [2e-2_dp, 1e-2_dp, 5e-3_dp, 2e-3_dp, 1e-3_dp, 5e-4_dp, 2e-4_dp, 1e-4_dp, &
       & 3e-5_dp, 1e-5_dp]
  real(dp), parameter :: zlo(nzb) = [0.05_dp, 0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, 0.05_dp]
  real(dp), parameter :: zhi(nzb) = [0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp, 0.5_dp]
  integer :: nv = 1                ! length of the acceptance vector: 1 (mode 0), ncell (mode 1)
  integer :: iv = 1                ! the cell VEGAS integrates (mode 1: smallest tau_cut, all tau_zQ)
  real(dp) :: hc(ncell), hc2, hcacc(ncell)
  ! mode 1: map the mass ratios and decay cosines logistically onto
  ! [1e-12, 1 - 1e-12] (both ends logarithmic), so that configurations near
  ! the 2+1 limits (tau_2 ~ 1e-5 needs two small invariants at once in the
  ! 4+1 real) are sampled; VEGAS adapts on top. Off in mode 0.
  logical :: logmap = .false.
  ! mode 1: slicing-adapted multichannel phase space (dis31/psmc.f90)
  logical :: usepsmc = .false.
  ! incoming-parton mask (diagnostics): 0 all, 1 quarks only, 2 gluon only
  integer :: pdfmask = 0
  ! mode 1: VEGAS adapts to the slice tcs(vslice) < tau_2 < tcs(vslice-1) of
  ! all tau_zQ bins (0: to tau_2 > tcs(ntc), the default); the cells are
  ! filled as before
  integer :: vslice = 0
  real(dp), parameter :: tmap = 27.631021115928547_dp     ! ln(1e12)
contains

  subroutine setup_flavours()
    integer :: f, Q, Q2, sg
    nb = 0; nr = 0
    do f = -5, 5
       if (f == 0) then
          do Q = 1, 5
             nb = nb + 1; flb(:,nb) = [0, Q, -Q, 0]; symb(nb) = 1
             nr = nr + 1; flr(:,nr) = [0, Q, -Q, 0, 0]; symr(nr) = 0.5_dp
             do Q2 = Q, 5
                nr = nr + 1; flr(:,nr) = [0, Q, -Q, Q2, -Q2]
                symr(nr) = merge(0.25_dp, 1.0_dp, Q2 == Q)
             enddo
          enddo
       else
          sg = sign(1, f)
          nb = nb + 1; flb(:,nb) = [f, f, 0, 0]; symb(nb) = 0.5_dp
          nr = nr + 1; flr(:,nr) = [f, f, 0, 0, 0]; symr(nr) = 1.0_dp/6
          do Q = 1, 5
             nb = nb + 1; flb(:,nb) = [f, f, sg*Q, -sg*Q]
             symb(nb) = merge(0.5_dp, 1.0_dp, Q == abs(f))
             nr = nr + 1; flr(:,nr) = [f, f, sg*Q, -sg*Q, 0]
             symr(nr) = merge(0.5_dp, 1.0_dp, Q == abs(f))
          enddo
       endif
    enddo
  end subroutine setup_flavours

  ! the integrand (pb) at the point r; wgt = VEGAS weight (for histograms)
  real(dp) function integrand(r, wgt) result(res)
    real(dp), intent(in) :: r(:), wgt
    select case (trim(part))
    case ('lo', 'vi', 'kp')
       res = born_part(r, wgt)
    case ('r')
       res = real_part(r, wgt)
    case default
       stop 'unknown part'
    end select
  end function integrand

  ! lepton variables and the incoming parton: returns Q2, y, xB, eta and
  ! the jacobian of dQ2 dy deta (the flux and 1/(16 pi^2) not included)
  subroutine lepton(r, Q2, y, xB, eta, jac, ok)
    real(dp), intent(in) :: r(3)
    real(dp), intent(out) :: Q2, y, xB, eta, jac
    logical, intent(out) :: ok
    if (mode == 1) then
       ! fixed (x, Q^2): dsigma/dx dQ^2 = (y/x) dsigma/dQ^2 dy
       Q2 = Q2fix; xB = xfix; y = Q2/(xB*s)
       ok = y < 1
       jac = 0; eta = 0
       if (.not. ok) return
       eta = xB*(1/xB)**r(3)
       jac = y/xB*eta*log(1/xB)
       return
    endif
    Q2 = q2min*(q2max/q2min)**r(1)
    y = ymin + (ymax - ymin)*r(2)
    xB = Q2/(y*s)
    ok = xB < 1
    jac = 0; eta = 0
    if (.not. ok) return
    eta = xB*(1/xB)**r(3)
    jac = Q2*log(q2max/q2min)*(ymax - ymin)*eta*log(1/xB)
  end subroutine lepton

  ! Breit-frame momenta in the dis31 layout: incoming parton eta P along +z,
  ! q along -z, leptons; the n outgoing partons from sequential decays of
  ! the hadronic system (W rest frame), boosted along z. Returns dPhi_n.
  subroutine breit(n, rr, Q2, y, xB, eta, Pk, dphi)
    integer, intent(in) :: n
    real(dp), intent(in) :: rr(:), Q2, y, xB, eta
    real(dp), intent(out) :: Pk(4,n+4), dphi
    real(dp) :: Q, W2, W, E, pz, beta, gam, k(4,4), m2(3), mtot2, pa(4), rest(4), rm(size(rr))
    integer :: i, j, ir
    Q = sqrt(Q2); E = eta*Q/(2*xB)
    ! mapped variables (mass ratios: the first n-2; cosines: every second after them)
    rm = rr
    dphi = 1
    if (logmap) then
       do i = 1, n - 2
          call lmap(rm(i))
       enddo
       do i = n - 1, size(rr), 2
          call lmap(rm(i))
       enddo
    endif
    Pk = 0
    Pk(:,1) = [0.0_dp, 0.0_dp, E, E]
    Pk(:,n+2) = [0.0_dp, 0.0_dp, -Q, 0.0_dp]
    Pk(:,n+3) = [Q/(2*y)*2*sqrt(1 - y), 0.0_dp, -Q/2, Q/(2*y)*(2 - y)]
    Pk(:,n+4) = [Q/(2*y)*2*sqrt(1 - y), 0.0_dp, Q/2, Q/(2*y)*(2 - y)]
    W2 = Q2*(eta/xB - 1); W = sqrt(W2)
    ! masses of the successive remainders: m2(1) = W2 > m2(2) > ... (last
    ! remainder massless pair)
    m2(1) = W2
    ir = 0
    do i = 2, n - 1
       ir = ir + 1
       m2(i) = m2(i-1)*rm(ir)
       dphi = dphi*m2(i-1)/(2*pi)
    enddo
    ! decays: remainder of mass^2 m2(i) -> parton i + remainder m2(i+1)
    rest = [0.0_dp, 0.0_dp, 0.0_dp, W]
    do i = 1, n - 1
       mtot2 = m2(i)
       if (i < n - 1) then
          call twobody(rest, 0.0_dp, m2(i+1), rm(ir+1), rm(ir+2), k(:,i), pa)
          dphi = dphi*(1 - m2(i+1)/mtot2)/(8*pi)
       else
          call twobody(rest, 0.0_dp, 0.0_dp, rm(ir+1), rm(ir+2), k(:,i), k(:,i+1))
          dphi = dphi/(8*pi)
       endif
       ir = ir + 2
       rest = pa
    enddo
    ! boost along z: hadronic system (E_h, 0, 0, pz)
    pz = E - Q
    beta = pz/E; gam = E/W
    do j = 1, n
       Pk(:,j+1) = [k(1,j), k(2,j), gam*(k(3,j) + beta*k(4,j)), gam*(k(4,j) + beta*k(3,j))]
    enddo
  contains
    ! logistic map of a unit variable, both ends logarithmic down to 1e-12;
    ! multiplies dphi by the jacobian
    subroutine lmap(v)
      real(dp), intent(inout) :: v
      real(dp) :: t
      t = tmap*(2*v - 1)
      v = 1/(1 + exp(-t))
      dphi = dphi*v*(1 - v)*2*tmap
    end subroutine lmap
  end subroutine breit

  ! massless + mass^2 m2b decay of p (mass^2 p.p) isotropically in its rest
  ! frame, then boosted to p's frame
  subroutine twobody(p, m2a, m2b, rc, rphi, pa, pb)
    real(dp), intent(in) :: p(4), m2a, m2b, rc, rphi
    real(dp), intent(out) :: pa(4), pb(4)
    real(dp) :: M, ea, pmod, c, st, phi, d(3), qa(4)
    M = sqrt(max(p(4)**2 - sum(p(1:3)**2), 0.0_dp))
    pmod = (M**2 - m2b)/(2*M)
    ea = pmod
    c = 2*rc - 1; st = sqrt(max(1 - c*c, 0.0_dp)); phi = 2*pi*rphi
    d = [st*cos(phi), st*sin(phi), c]
    qa = [pmod*d, ea]
    pa = boostv(qa, p, M)
    pb = p - pa
    if (.false.) print *, m2a
  end subroutine twobody

  function boostv(q, Pt, M) result(r)
    real(dp), intent(in) :: q(4), Pt(4), M
    real(dp) :: r(4), bp, f
    bp = dot_product(q(1:3), Pt(1:3))
    r(4) = (Pt(4)*q(4) + bp)/M
    f = (bp/(Pt(4) + M) + q(4))/M
    r(1:3) = q(1:3) + f*Pt(1:3)
  end function boostv

  ! PDF number densities f(-5:5) at x, scale mu
  subroutine pdfs(x, mu, f)
    real(dp), intent(in) :: x, mu
    real(dp), intent(out) :: f(-5:5)
    real(dp) :: xf(-6:6)
    call evolvePDF(x, mu, xf)
    f = xf(-5:5)/x
    if (pdfmask == 1) f(0) = 0
    if (pdfmask == 2) then
       f(-5:-1) = 0; f(1:5) = 0
    endif
  end subroutine pdfs

  ! charge structures: the matrix elements are bilinear in the quark
  ! charges and invariant under charge conjugation, so each channel topology
  ! needs one (e_q^2 only) or three (e_q^2, e_Q^2, e_q e_Q) evaluations per
  ! point; the flavour sums are done with the PDFs and charges
  pure real(dp) function ech(f)
    integer, intent(in) :: f
    ech = merge(2.0_dp/3, -1.0_dp/3, mod(abs(f), 2) == 0)
  end function ech

  ! coefficients (x, y, z) of ea^2, eb^2, ea eb from the values m at
  ! (ea, eb) = (-1/3, 2/3), (2/3, -1/3), (-1/3, -1/3) (flavours d u, u d, d s)
  subroutine solve3(m, c)
    real(dp), intent(in) :: m(3)
    real(dp), intent(out) :: c(3)
    real(dp) :: A(3,3), ea(3), eb(3), det, B(3,3)
    integer :: i
    ea = [-1.0_dp/3, 2.0_dp/3, -1.0_dp/3]; eb = [2.0_dp/3, -1.0_dp/3, -1.0_dp/3]
    A(:,1) = ea**2; A(:,2) = eb**2; A(:,3) = ea*eb
    det = A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2)) - A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1)) &
         & + A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))
    do i = 1, 3
       B = A; B(:,i) = m
       c(i) = (B(1,1)*(B(2,2)*B(3,3) - B(2,3)*B(3,2)) - B(1,2)*(B(2,1)*B(3,3) - B(2,3)*B(3,1)) &
            & + B(1,3)*(B(2,1)*B(3,2) - B(2,2)*B(3,1)))/det
    enddo
  end subroutine solve3

  ! Born-level value(s) of one flavour assignment at the point: lo: |M|^2;
  ! vi: finite V + I; kp: (b, g, lp)
  subroutine born_eval(Pk, fl, Q2, val)
    real(dp), intent(in) :: Pk(4,7), Q2
    integer, intent(in) :: fl(4)
    real(dp), intent(out) :: val(3)
    real(dp) :: cc(4,4), v(-2:0), t, iv(-2:0)
    val = 0
    select case (trim(part))
    case ('lo')
       call born31_cc(Pk, fl, val(1), cc)
    case ('vi')
       call virt31_ren(Pk, fl, Q2, v, t)
       call iop31_i(Pk, fl, Q2, iv)
       val(1) = v(0) + iv(0)
    case ('kp')
       call iop31_kp(Pk, fl, Q2, val(1), val(2), val(3))
    end select
  end subroutine born_eval

  real(dp) function born_part(r, wgt) result(res)
    real(dp), intent(in) :: r(:), wgt
    real(dp), external :: alphasPDF
    real(dp) :: Q2, y, xB, eta, jac, Pk(4,7), dphi, as, fpdf(-5:5), w, x, fa(-5:5)
    real(dp) :: u1(3), u2(3), u4(3), m3(3,3), c3(3,3), t(3), eq1, eq2, acc(nv)
    logical :: ok
    integer :: f, Q, k
    res = 0
    if (usepsmc) then
       Q2 = Q2fix; xB = xfix; y = Q2/(xB*s)
       call psmc_gen(3, r(1:7), Pk, eta, dphi, ok)
       if (.not. ok) return
       jac = y/xB
    else
       call lepton(r(1:3), Q2, y, xB, eta, jac, ok)
       if (.not. ok) return
       call breit(3, r(4:8), Q2, y, xB, eta, Pk, dphi)
    endif
    call accept(Pk(:,1), Pk(:,2:4), 3, sqrt(Q2), acc)
    if (all(acc == 0)) return
    as = alphasPDF(sqrt(Q2))
    call pdfs(eta, sqrt(Q2), fpdf)
    w = jac*dphi/(2*eta*s)/(16*pi**2)*gev2pb*(as/(2*pi))**2
    if (trim(part) /= 'lo') w = w*as/(2*pi)
    ! unit-charge values: q g g (d), g -> q qbar g (d), identical (d)
    call born_eval(Pk, [1, 1, 0, 0], Q2, u1); u1 = u1*9
    call born_eval(Pk, [0, 1, -1, 0], Q2, u2); u2 = u2*9
    call born_eval(Pk, [1, 1, 1, -1], Q2, u4); u4 = u4*9
    ! different flavours: (e_q, e_Q) = (d,u), (u,d), (d,s)
    call born_eval(Pk, [1, 1, 2, -2], Q2, m3(:,1))
    call born_eval(Pk, [2, 2, 1, -1], Q2, m3(:,2))
    call born_eval(Pk, [1, 1, 3, -3], Q2, m3(:,3))
    do k = 1, 3
       call solve3(m3(k,:), c3(k,:))
    enddo
    if (trim(part) == 'kp') then
       x = eta + (1 - eta)*r(9)
       call pdfs(eta/x, sqrt(Q2), fa)
    endif
    do f = -5, 5
       if (f == 0) then
          ! g -> q qbar g
          do Q = 1, 5
             res = res + bterm(0, ech(Q)**2*u2)
          enddo
          cycle
       endif
       eq1 = ech(f)
       res = res + 0.5_dp*bterm(f, eq1**2*u1)
       res = res + 0.5_dp*bterm(f, eq1**2*u4)
       do Q = 1, 5
          if (Q == abs(f)) cycle
          eq2 = ech(Q)
          t = eq1**2*c3(:,1) + eq2**2*c3(:,2) + eq1*eq2*c3(:,3)
          res = res + bterm(f, t)
       enddo
    enddo
    res = res*w
    if (mode == 1) then
       if (res /= res) then
          res = 0; return
       endif
       hcacc = hcacc + res*wgt*acc
       res = res*vtarget(acc)
    elseif (res /= 0) then
       call fill(Q2, res*wgt)
    endif
  contains
    ! one Born flavour class with incoming fb and values val(3)
    real(dp) function bterm(fb, val)
      integer, intent(in) :: fb
      real(dp), intent(in) :: val(3)
      real(dp) :: cp
      integer :: ia
      select case (trim(part))
      case ('lo', 'vi')
         bterm = fpdf(fb)*val(1)
      case default
         call kp_conv(merge(1, 0, fb /= 0), merge(1, 0, fb /= 0), .true., val(1), val(2), val(3), &
              & x, eta, fa(fb), fpdf(fb), bterm)
         if (fb /= 0) then
            call kp_conv(0, 1, .false., val(1), val(2), val(3), x, eta, fa(0), fpdf(0), cp)
            bterm = bterm + cp
         else
            do ia = -5, 5
               if (ia == 0) cycle
               call kp_conv(1, 0, .false., val(1), val(2), val(3), x, eta, fa(ia), fpdf(ia), cp)
               bterm = bterm + cp
            enddo
         endif
      end select
    end function bterm
  end function born_part

  ! the x-convolution of K + P for incoming type ta (1 quark, 0 gluon) into
  ! Born type tb, at the point x (uniform in (xi, 1), jacobian 1 - xi):
  !   int_xi^1 dx [R(x) f(xi/x)/x + Pl(x) (f(xi/x)/x - f(xi))]
  !   - f(xi) int_0^xi Pl(x) dx + f(xi) Delta
  ! fx = f(xi/x), fxi = f(xi) (number densities)
  subroutine kp_conv(ta, tb, diag, b, g, lp, x, xi, fx, fxi, res)
    integer, intent(in) :: ta, tb
    logical, intent(in) :: diag
    real(dp), intent(in) :: b, g, lp, x, xi, fx, fxi
    real(dp), intent(out) :: res
    real(dp), external :: ddilog
    real(dp) :: kr, kpl, kd, pr, ppl, pd, R, Pl, D, ipl, A, Bf, Cf, l1
    call iop31_kernel(1, ta, tb, x, kr, kpl, kd)
    call iop31_kernel(2, ta, tb, x, pr, ppl, pd)
    R = kr*b + pr*lp
    Pl = kpl*b + ppl*lp
    D = kd*b + pd*lp
    ipl = 0
    if (diag) then
       Pl = Pl + g/(1 - x)
       D = D + g
       ! int_0^xi of the plus parts
       l1 = log(1 - xi)
       A = -0.5_dp*l1**2 + log(xi)*l1 + ddilog(xi)      ! ln((1-x)/x)/(1-x)
       Bf = -l1                                          ! 1/(1-x)
       Cf = -2*l1 - xi - xi**2/2                         ! (1+x^2)/(1-x)
       if (tb == 1) then
          ipl = CF*2*A*b + g*Bf + CF*Cf*lp
       else
          ipl = 2*CA*A*b + g*Bf + 2*CA*Bf*lp
       endif
    endif
    res = (1 - xi)*(R*fx/x + Pl*(fx/x - fxi)) - fxi*ipl + fxi*D
  end subroutine kp_conv

  ! real minus dipoles for one flavour assignment (with jet functions)
  real(dp) function real_eval(Pk, fl, pass4) result(sub)
    real(dp), intent(in) :: Pk(4,8)
    integer, intent(in) :: fl(5)
    logical, intent(in) :: pass4
    real(dp) :: subv(nv), F4(nv)
    F4 = merge(1.0_dp, 0.0_dp, pass4)
    call real_evalv(Pk, fl, F4, subv)
    sub = subv(1)
  end function real_eval

  ! the same with acceptance vectors: F4 for the real, accept() for each
  ! dipole's mapped Born
  subroutine real_evalv(Pk, fl, F4, sub)
    real(dp), intent(in) :: Pk(4,8), F4(nv)
    integer, intent(in) :: fl(5)
    real(dp), intent(out) :: sub(nv)
    real(dp) :: m4, P3(4,7,dip41_max), val(dip41_max), F3(nv)
    integer :: fl3(4,dip41_max), nd, id
    sub = 0
    if (any(F4 /= 0)) then
       call me41_tree(Pk, fl, m4)
       sub = m4*F4
    endif
    call dip41_list(Pk, fl, nd, P3, fl3, val)
    do id = 1, nd
       call accept(P3(:,1,id), P3(:,2:4,id), 3, sqrt(-mdot(P3(:,5,id), P3(:,5,id))), F3)
       sub = sub - val(id)*F3
    enddo
  end subroutine real_evalv

  real(dp) function real_part(r, wgt) result(res)
    real(dp), intent(in) :: r(:), wgt
    real(dp), external :: alphasPDF
    real(dp) :: Q2, y, xB, eta, jac, Pk(4,8), dphi, as, fpdf(-5:5), w, smin, W2
    real(dp) :: u1(nv), u2(nv), u4(nv), u6(nv), m3(nv,3), c3(nv,3), m5(nv,3), c5(nv,3), eq1, eq2, sg(nv), F4(nv)
    integer :: i, j, f, Q, Q2i, ic
    logical :: ok
    res = 0
    if (usepsmc) then
       Q2 = Q2fix; xB = xfix; y = Q2/(xB*s)
       call psmc_gen(4, r(1:11), Pk, eta, dphi, ok)
       if (.not. ok) return
       jac = y/xB
    else
       call lepton(r(1:3), Q2, y, xB, eta, jac, ok)
       if (.not. ok) return
       call breit(4, r(4:11), Q2, y, xB, eta, Pk, dphi)
    endif
    W2 = Q2*(eta/xB - 1)
    smin = huge(1.0_dp)
    do i = 1, 5
       do j = i + 1, 5
          smin = min(smin, 2*abs(mdot(Pk(:,i), Pk(:,j))))
       enddo
    enddo
    if (smin < techcut*W2) return
    call accept(Pk(:,1), Pk(:,2:5), 4, sqrt(Q2), F4)
    as = alphasPDF(sqrt(Q2))
    call pdfs(eta, sqrt(Q2), fpdf)
    w = jac*dphi/(2*eta*s)/(16*pi**2)*gev2pb*(as/(2*pi))**3
    ! unit-charge values
    call real_evalv(Pk, [1, 1, 0, 0, 0], F4, u1); u1 = 9*u1        ! q -> q g g g
    call real_evalv(Pk, [0, 1, -1, 0, 0], F4, u2); u2 = 9*u2       ! g -> q qbar g g
    call real_evalv(Pk, [1, 1, 1, -1, 0], F4, u4); u4 = 9*u4       ! q -> q q qbar g
    call real_evalv(Pk, [0, 1, -1, 1, -1], F4, u6); u6 = 9*u6      ! g -> q qbar q qbar
    call real_evalv(Pk, [1, 1, 2, -2, 0], F4, m3(:,1)); call real_evalv(Pk, [2, 2, 1, -1, 0], F4, m3(:,2))
    call real_evalv(Pk, [1, 1, 3, -3, 0], F4, m3(:,3))
    call real_evalv(Pk, [0, 1, -1, 2, -2], F4, m5(:,1)); call real_evalv(Pk, [0, 2, -2, 1, -1], F4, m5(:,2))
    call real_evalv(Pk, [0, 1, -1, 3, -3], F4, m5(:,3))
    do ic = 1, nv
       call solve3(m3(ic,:), c3(ic,:)); call solve3(m5(ic,:), c5(ic,:))
    enddo
    sg = 0
    do f = -5, 5
       if (f == 0) cycle
       eq1 = ech(f)
       sg = sg + fpdf(f)*eq1**2*(u1/6 + u4/2)
       do Q = 1, 5
          if (Q == abs(f)) cycle
          eq2 = ech(Q)
          sg = sg + fpdf(f)*(eq1**2*c3(:,1) + eq2**2*c3(:,2) + eq1*eq2*c3(:,3))
       enddo
    enddo
    do Q = 1, 5
       eq1 = ech(Q)
       sg = sg + fpdf(0)*eq1**2*(u2/2 + u6/4)
       do Q2i = Q + 1, 5
          eq2 = ech(Q2i)
          sg = sg + fpdf(0)*(eq1**2*c5(:,1) + eq2**2*c5(:,2) + eq1*eq2*c5(:,3))
       enddo
    enddo
    if (mode == 1) then
       ! points where a matrix element is not finite (extreme configurations
       ! of the logmap sampling) are dropped, as VEGAS does for the integral
       if (any(sg /= sg)) then
          res = 0; return
       endif
       hcacc = hcacc + sg*w*wgt
       res = vtarget(sg)*w
    else
       res = sg(1)*w
       if (res /= 0) call fill(Q2, res*wgt)
    endif
  end function real_part

  ! the quantity VEGAS integrates in mode 1: the cell iv, or a tau_2 slice
  real(dp) function vtarget(a) result(t)
    real(dp), intent(in) :: a(nv)
    integer :: k
    if (vslice <= 1) then
       t = a(iv)
    else
       k = vslice + ntc*(nzb - 1)
       t = a(k) - a(k - 1)
    endif
  end function vtarget

  ! acceptance vector F(nv) of the n outgoing partons p(4,n) (incoming parton
  ! pin): mode 0: F(1) = (>= njmin jets); mode 1: F(k + ntc*(b-1)) =
  ! theta(tau_2 > tcs(k)) theta(zlo(b) <= tau_zQ < zhi(b))
  subroutine accept(pin, p, n, Q, F)
    integer, intent(in) :: n
    real(dp), intent(in) :: pin(4), p(4,n), Q
    real(dp), intent(out) :: F(nv)
    real(dp) :: t2, tz
    integer :: k, b, i
    if (mode == 0) then
       F(1) = merge(1.0_dp, 0.0_dp, njets(p, n) >= njmin)
       return
    endif
    F = 0
    tz = 1
    do i = 1, n
       if (p(3,i) < 0) tz = tz + 2*p(3,i)/Q
    enddo
    t2 = tau2cm(pin, p, n)/Q
    do b = 1, nzb
       if (tz < zlo(b) .or. tz >= zhi(b)) cycle
       do k = 1, ntc
          if (t2 > tcs(k)) F(k + ntc*(b - 1)) = 1
       enddo
    enddo
  end subroutine accept

  ! T_2 in the jets' rest frame (slicing/mod_tau2_run.f90, measure 2). The
  ! frame u = (P_J1 + P_J2)/m of the partition that minimises the geometric
  ! T_2 in the input (Breit) frame; then the exact minimum over partitions
  ! of sum_beam P.p/P.u + sum_jets m^2/(u.P + |P|_u), all in that one frame.
  ! (A frame per partition is not IR safe: two collinear partons as the two
  ! jets give a degenerate, infinitely boosted frame and T -> 0.)
  real(dp) function tau2cm(pin, p, n) result(T)
    integer, intent(in) :: n
    real(dp), intent(in) :: pin(4), p(4,n)
    real(dp) :: u(4), ub(4), tb, best
    integer :: code
    ! step 1: the Breit-frame partition (u = (1,0,0,0): energies in the input frame)
    ub = [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
    best = huge(1.0_dp); u = 0
    do code = 0, 3**n - 1
       tb = tpart(code, ub, .true.)
       if (tb < best) then
          best = tb; u = ujets(code)
       endif
    enddo
    T = 0
    if (best == huge(1.0_dp) .or. mdot(u, u) <= 0) return
    u = u/sqrt(mdot(u, u))
    ! step 2: exact minimisation in that frame
    T = huge(1.0_dp)
    do code = 0, 3**n - 1
       T = min(T, tpart(code, u, .false.))
    enddo
    T = max(T, 0.0_dp)
  contains
    ! the partition of code (0 beam, 1, 2 jets; both jets non-empty, jet 1
    ! holds the first jet parton), or huge if not allowed
    subroutine decode(code, a, okp)
      integer, intent(in) :: code
      integer, intent(out) :: a(6)
      logical, intent(out) :: okp
      integer :: m, k
      m = code
      do k = 1, n
         a(k) = mod(m, 3); m = m/3
      enddo
      okp = any(a(1:n) == 1) .and. any(a(1:n) == 2)
      if (.not. okp) return
      do k = 1, n
         if (a(k) > 0) then
            okp = a(k) == 1; return
         endif
      enddo
    end subroutine decode
    function ujets(code) result(w)
      integer, intent(in) :: code
      real(dp) :: w(4)
      integer :: a(6), k
      logical :: okp
      call decode(code, a, okp)
      w = 0
      do k = 1, n
         if (a(k) > 0) w = w + p(:,k)
      enddo
    end function ujets
    real(dp) function tpart(code, uu, breitf) result(tp)
      integer, intent(in) :: code
      real(dp), intent(in) :: uu(4)
      logical, intent(in) :: breitf
      integer :: a(6), k
      logical :: okp
      real(dp) :: PJ(4,2), uP, m2
      call decode(code, a, okp)
      tp = huge(1.0_dp)
      if (.not. okp) return
      PJ = 0
      do k = 1, n
         if (a(k) > 0) PJ(:,a(k)) = PJ(:,a(k)) + p(:,k)
      enddo
      tp = 0
      do k = 1, n
         if (a(k) == 0) tp = tp + mdot(pin, p(:,k))/mdot(pin, uu)
      enddo
      do k = 1, 2
         uP = mdot(uu, PJ(:,k)); m2 = max(mdot(PJ(:,k), PJ(:,k)), 0.0_dp)
         tp = tp + m2/(uP + sqrt(max(uP*uP - m2, 0.0_dp)))
      enddo
      if (.false. .and. breitf) tp = tp
    end function tpart
  end function tau2cm

  ! inclusive kt (R = rjet, E-scheme) on n massless partons in the Breit
  ! frame (z = the proton direction); number of jets with p_T > ptmin
  integer function njets(p, n) result(nj)
    integer, intent(in) :: n
    real(dp), intent(in) :: p(4,n)
    real(dp) :: q(4,n), dmin, d, pt2(n), yr(n), ph(n), dy, dphi
    integer :: m, i, j, ii, jj, act(n)
    logical :: beam
    q = p; m = n
    act = [(i, i = 1, n)]
    nj = 0
    do while (m > 0)
       do i = 1, m
          call ptyphi(q(:,act(i)), pt2(i), yr(i), ph(i))
       enddo
       dmin = huge(1.0_dp); ii = 0; jj = 0; beam = .true.
       do i = 1, m
          if (pt2(i) < dmin) then
             dmin = pt2(i); ii = i; beam = .true.
          endif
          do j = i + 1, m
             dy = yr(i) - yr(j)
             dphi = abs(ph(i) - ph(j)); if (dphi > pi) dphi = 2*pi - dphi
             d = min(pt2(i), pt2(j))*(dy**2 + dphi**2)/rjet**2
             if (d < dmin) then
                dmin = d; ii = i; jj = j; beam = .false.
             endif
          enddo
       enddo
       if (beam) then
          if (pt2(ii) > ptmin**2) nj = nj + 1
          act(ii) = act(m); m = m - 1
       else
          q(:,act(ii)) = q(:,act(ii)) + q(:,act(jj))
          act(jj) = act(m); m = m - 1
       endif
    enddo
  end function njets

  subroutine ptyphi(p, pt2, y, phi)
    real(dp), intent(in) :: p(4)
    real(dp), intent(out) :: pt2, y, phi
    pt2 = p(1)**2 + p(2)**2
    y = 0.5_dp*log(max(p(4) + p(3), 1d-300)/max(p(4) - p(3), 1d-300))
    phi = atan2(p(2), p(1))
  end subroutine ptyphi

  subroutine fill(Q2, w)
    real(dp), intent(in) :: Q2, w
    integer :: i
    do i = 1, nq
       if (Q2 >= qedge(i-1) .and. Q2 < qedge(i)) hacc(i) = hacc(i) + w
    enddo
  end subroutine fill

  pure real(dp) function mdot(a, b)
    real(dp), intent(in) :: a(4), b(4)
    mdot = a(4)*b(4) - a(1)*b(1) - a(2)*b(2) - a(3)*b(3)
  end function mdot

  ! VEGAS (Lepage), importance sampling with a factorised grid; per
  ! iteration estimates combined with inverse-variance weights. Histograms:
  ! per-iteration sums of w f, combined the same way as the integral.
  subroutine vegas(ndim, ncall, itmx, avg, err, chi2, fun)
    integer, intent(in) :: ndim, ncall, itmx
    interface
       real(dp) function fun(r, wgt)
         import :: dp
         real(dp), intent(in) :: r(:), wgt
       end function fun
    end interface
    real(dp), intent(out) :: avg, err, chi2
    integer, parameter :: nbin = 50
    real(dp) :: xi(0:nbin,ndim), d(nbin,ndim), r(ndim), x(ndim), jac, f, f2, s1, s2, wsum, sumw, sumwi
    real(dp) :: hit(nbin), dt, rc, xin(0:nbin), est(200), var(200), xo, xn
    integer :: ia(ndim), it, ic, j, k, i
    do j = 1, ndim
       xi(:,j) = [(real(i, dp)/nbin, i = 0, nbin)]
    enddo
    hist = 0; hist2 = 0; hc = 0; hc2 = 0
    sumw = 0; sumwi = 0
    do it = 1, itmx
       d = 0; s1 = 0; s2 = 0; hacc = 0; hcacc = 0
       do ic = 1, ncall
          call random_number(r)
          jac = 1
          do j = 1, ndim
             xn = r(j)*nbin
             ia(j) = min(int(xn) + 1, nbin)
             xo = xi(ia(j),j) - xi(ia(j)-1,j)
             x(j) = xi(ia(j)-1,j) + (xn - (ia(j) - 1))*xo
             jac = jac*xo*nbin
          enddo
          f = fun(x, jac/ncall)*jac
          if (f /= f) f = 0
          s1 = s1 + f; s2 = s2 + f*f
          do j = 1, ndim
             d(ia(j),j) = d(ia(j),j) + f*f
          enddo
       enddo
       s1 = s1/ncall
       s2 = max((s2/ncall - s1**2)/(ncall - 1), 1d-300)
       est(it) = s1; var(it) = s2
       ! histograms: this iteration's estimate, weighted like the integral
       if (it >= min(2, itmx)) then
          hist = hist + hacc/s2; hist2 = hist2 + 1/s2
          hc = hc + hcacc/s2; hc2 = hc2 + 1/s2
       endif
       write(*,'(a,i3,a,es16.8,a,es12.4)') ' iteration', it, ':', s1, ' +-', sqrt(s2)
       flush(6)
       ! refine the grid
       do j = 1, ndim
          ! smooth
          hit = d(:,j)
          dt = sum(hit)
          if (dt <= 0) cycle
          do k = 1, nbin
             hit(k) = d(k,j)
          enddo
          hit(1) = (d(1,j) + d(2,j))/2
          hit(nbin) = (d(nbin-1,j) + d(nbin,j))/2
          do k = 2, nbin - 1
             hit(k) = (d(k-1,j) + d(k,j) + d(k+1,j))/3
          enddo
          dt = sum(hit)
          do k = 1, nbin
             rc = hit(k)/dt
             if (rc > 0 .and. abs(rc - 1) > 1d-12) then
                hit(k) = ((rc - 1)/log(rc))**1.5_dp
             elseif (rc > 0) then
                hit(k) = 1
             endif
          enddo
          call rebin(hit, nbin, xi(:,j), xin)
          xi(:,j) = xin
       enddo
    enddo
    ! combine iterations 2..itmx (the first trains the grid)
    wsum = 0; avg = 0
    do it = min(2, itmx), itmx
       avg = avg + est(it)/var(it); wsum = wsum + 1/var(it)
    enddo
    avg = avg/wsum; err = sqrt(1/wsum)
    chi2 = 0
    do it = min(2, itmx), itmx
       chi2 = chi2 + (est(it) - avg)**2/var(it)
    enddo
    chi2 = chi2/max(itmx - min(2, itmx), 1)
    hist = hist/hist2
    if (hc2 > 0) hc = hc/hc2
    if (.false.) print *, f2, sumw, sumwi
  end subroutine vegas

  subroutine rebin(w, nbin, xold, xnew)
    integer, intent(in) :: nbin
    real(dp), intent(in) :: w(nbin), xold(0:nbin)
    real(dp), intent(out) :: xnew(0:nbin)
    real(dp) :: target, acc, wtot
    integer :: k, i
    wtot = sum(w)
    xnew(0) = 0; xnew(nbin) = 1
    acc = 0; k = 0
    do i = 1, nbin - 1
       target = wtot*i/nbin
       do while (acc + w(k+1) < target)
          acc = acc + w(k+1); k = k + 1
       enddo
       if (w(k+1) > 0) then
          xnew(i) = xold(k) + (target - acc)/w(k+1)*(xold(k+1) - xold(k))
       else
          xnew(i) = xold(k+1)
       endif
    enddo
  end subroutine rebin
  ! consistency check of the charge-structure sums against the explicit sum
  ! over flavour assignments (lists flb, flr with symmetry factors)
  subroutine check_sums()
    real(dp) :: r(11), a, b, Q2, y, xB, eta, jac, Pk(4,8), Pb(4,7), dphi, fpdf(-5:5), val(3), W2
    real(dp) :: P3(4,7,dip41_max), dv(dip41_max), m4
    integer :: ipt, ib, ir, fl3(4,dip41_max), nd, id
    logical :: ok, pass4
    do ipt = 1, 20
       call random_number(r)
       ! Born
       part = 'lo'
       a = born_part(r(1:8), 0.0_dp)
       call lepton(r(1:3), Q2, y, xB, eta, jac, ok)
       b = 0
       if (ok) then
          call breit(3, r(4:8), Q2, y, xB, eta, Pb, dphi)
          if (njets(Pb(:,2:4), 3) >= njmin) then
             call pdfs(eta, sqrt(Q2), fpdf)
             do ib = 1, nb
                call born_eval(Pb, flb(:,ib), Q2, val)
                b = b + symb(ib)*fpdf(flb(1,ib))*val(1)
             enddo
             b = b*jac*dphi/(2*eta*s)/(16*pi**2)*gev2pb*(alphasPDF_(sqrt(Q2))/(2*pi))**2
          endif
       endif
       write(*,'(a,2es18.10)') ' born: structures, explicit', a, b
       ! real
       part = 'r'
       a = real_part(r, 0.0_dp)
       b = 0
       if (ok) then
          call breit(4, r(4:11), Q2, y, xB, eta, Pk, dphi)
          pass4 = njets(Pk(:,2:5), 4) >= njmin
          call pdfs(eta, sqrt(Q2), fpdf)
          do ir = 1, nr
             b = b + symr(ir)*fpdf(flr(1,ir))*real_eval(Pk, flr(:,ir), pass4)
          enddo
          b = b*jac*dphi/(2*eta*s)/(16*pi**2)*gev2pb*(alphasPDF_(sqrt(Q2))/(2*pi))**3
          W2 = Q2*(eta/xB - 1)
       endif
       write(*,'(a,2es18.10)') ' real: structures, explicit', a, b
    enddo
    if (.false.) print *, P3, dv, m4, fl3, nd, id, W2
  end subroutine check_sums

  real(dp) function alphasPDF_(q)
    real(dp), intent(in) :: q
    real(dp), external :: alphasPDF
    alphasPDF_ = alphasPDF(q)
  end function alphasPDF_
end module nlo31_mod

program nlo31
  use nlo31_mod
  implicit none
  character(32) :: arg
  integer :: ncall, itmx, ndim, seed, nseed, i
  integer, allocatable :: sd(:)
  real(dp) :: avg, err, chi2
  call get_command_argument(1, part)
  call get_command_argument(2, arg); read(arg, *) ncall
  call get_command_argument(3, arg); read(arg, *) itmx
  seed = 1
  if (command_argument_count() > 3) then
     call get_command_argument(4, arg); read(arg, *) seed
  endif
  if (command_argument_count() > 4) then
     call get_command_argument(5, arg); read(arg, *) techcut
  endif
  if (command_argument_count() > 7) then
     call get_command_argument(6, arg); read(arg, *) mode
     call get_command_argument(7, arg); read(arg, *) xfix
     call get_command_argument(8, arg); read(arg, *) Q2fix
  endif
  if (command_argument_count() > 8) then
     call get_command_argument(9, arg); logmap = trim(arg) == 'logmap'; usepsmc = trim(arg) == 'psmc'
  endif
  if (command_argument_count() > 9) then
     call get_command_argument(10, arg); read(arg, *) pdfmask
  endif
  if (command_argument_count() > 10) then
     call get_command_argument(11, arg); read(arg, *) vslice
  endif
  if (mode == 1) then
     nv = ncell; iv = ntc + ntc*(nzb - 1)
     if (usepsmc) call psmc_init(xfix, Q2fix, Q2fix/(xfix*s))
  endif
  call random_seed(size=nseed); allocate(sd(nseed))
  sd = [(1000003*seed + 7919*i, i = 1, nseed)]
  call random_seed(put=sd)
  call initPDFSetByName('NNPDF30_nlo_as_0118')
  call initPDF(0)
  call setup_flavours()
  virt31_finite_only = .true.
  if (trim(part) == 'chk') then
     call check_sums(); stop
  endif
  select case (trim(part))
  case ('lo', 'vi'); ndim = 8
  case ('kp'); ndim = 9
  case ('r'); ndim = 11
  case default; stop 'part: lo, vi, kp or r'
  end select
  write(*,'(a,a,a,i10,a,i4,a,i6,a,es9.2,a,i2,a,f9.6,a,f10.2)') ' nlo31 part ', trim(part), ' ncall', ncall, &
       & ' itmx', itmx, ' seed', seed, ' techcut', techcut, ' mode', mode, ' x', xfix, ' Q2', Q2fix
  call vegas(ndim, ncall, itmx, avg, err, chi2, integrand)
  write(*,'(a,a,a,es16.8,a,es12.4,a,f8.3)') ' RESULT ', trim(part), ' sigma(>=3 jets) [pb] = ', avg, ' +- ', err, &
       & '   chi2/it', chi2
  if (mode == 1) then
     write(*,'(a)') ' CELLS dsigma/dx dQ2 [pb/GeV2] above the cut: tau_zQ bin, then tau_cut columns'
     write(*,'(a,10es12.3)') ' tau_cut     ', tcs
     do i = 1, nzb
        write(*,'(a,2f5.2,10es16.8)') ' CELL ', zlo(i), zhi(i), hc(1 + ntc*(i - 1):ntc*i)
     enddo
  else
     do i = 1, nq
        write(*,'(a,2f9.1,es16.8)') ' Q2bin', qedge(i-1), qedge(i), hist(i)
     enddo
  endif
end program nlo31
