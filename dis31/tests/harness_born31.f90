! Tests of born31 (colour- and spin-correlated DIS 3+1 Borns):
!  1. msq = me31, all flavour assignments;
!  2. colour conservation: sum_{k /= i} cc(i,k) = -C_i msq;
!  3. polarisation sum: msqv(e1) + msqv(e2) = msq, ccv(e1) + ccv(e2) = cc
!     (e1, e2 orthonormal, transverse to the gluon and a reference vector);
!  4. soft gluon: me41 -> -8 pi^2 sum_{a /= b} p_a.p_b/(p_a.q p_b.q) cc(a,b)
!     (gluon from the CS FF map with 1 - z = lambda a, y = lambda b);
!  5. collinear limits at fixed azimuth (no averaging), which need the spin
!     correlations: final g -> g g, g -> q qbar (FF map),
!     initial q -> g + q, g -> g + g (IF map), against
!       FF: 16 pi^2/s_ij <P^{mu nu}>, IF: (1/x) 16 pi^2/(2 pa.pi) <P^{mu nu}>,
!     P_gg = 2 CA [-g (z/(1-z) + (1-z)/z) - 2 z(1-z) kk/k^2],
!     P_qqbar = TR [-g + 4 z(1-z) kk/k^2],
!     initial q -> g: CF [-g x - 4 (1-x)/x kk/k^2] (this sign of the kk term
!     is the one that reproduces the azimuthal dependence; averaged it gives
!     CF (1 + (1-x)^2)/x),
!     initial g -> g: 2 CA [-g (x/(1-x) + x(1-x)) - 2 (1-x)/x kk/k^2];
!     at y (u) = 1e-8 and 1e-10 (three azimuths each), the criterion on the
!     smaller one (smaller values lose digits to round-off in double
!     precision; wrong spin correlations give O(1) deviations).
! Couplings: photon, or with the environment EW31 = "mode lepton" those of
! ew31 (e.g. "1 0": photon + Z, e-).
program harness_born31
  use me31
  use me41
  use born31
  use ew31
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: CF = 4.0_dp/3, CA = 3, TR = 0.5_dp
  real(dp) :: P3(4,7), r(8), worst(5)
  integer :: ipt
  character(32) :: arg
  call get_environment_variable('EW31', arg)
  if (len_trim(arg) > 0) read(arg, *) ew31_mode, ew31_lepton
  call seed_rng()
  worst = 0
  do ipt = 1, 4
     ! a resolved 3+1 point (all 2 p_i.p_j above 2% of W^2), as a jet
     ! function would require; near-singular Borns converge slowly
     do
        call random_number(r)
        call dis_point(r, P3)
        if (resolved3(P3, 0.02_dp)) exit
     enddo
     write(*,'(a,i2)') ' point', ipt
     call check_flavours(P3)
     if (ew31_mode == 2) then
        ! W exchange (EW31 = "2 0", 9 Oct): line u -> d, pair (d, ubar)
        call soft('W q;qgg ', P3, wf([2,1,0,0]), 3, 4)
        call soft('W q;qgg ', P3, wf([2,1,0,0]), 2, 3)
        call soft('W qb;qbgg', P3, wf([-1,-2,0,0]), 4, 2)
        call soft('W g;qqbg', P3, wf([0,1,-2,0]), 2, 4)
        call soft('W g;qqbg', P3, wf([0,1,-2,0]), 3, 2)
        call soft('W q;qQQb', P3, wf([2,1,3,-3]), 3, 4)
        call soft('W q;qqqb', P3, wf([2,1,1,-1]), 2, 4)
        call soft('W q;duub', P3, wf([2,1,2,-2]), 2, 4)
        call soft('W d;ddub', P3, wf([1,1,1,-2]), 3, 4)
        call soft('W qb;QQb', P3, wf([-1,-2,3,-3]), 3, 2)
        call ffspin('W q;qgg g>gg', P3, wf([2,1,0,0]), 3, 2, wf([2,1,0,0,0]), 'ggg', 0.3_dp)
        call ffspin('W g;qqbg g>gg', P3, wf([0,1,-2,0]), 4, 3, wf([0,1,-2,0,0]), 'ggg', 0.6_dp)
        call ffspin('W q;qgg g>qqb', P3, wf([2,1,0,0]), 4, 3, wf([2,1,0,2,-2]), 'gqq', 0.35_dp)
        call ffspin('W g;qqbg g>qqb', P3, wf([0,1,-2,0]), 4, 2, wf([0,1,-2,1,-1]), 'gqq', 0.5_dp)
        call fispin('W q;g IS q>g', P3, wf([0,1,-2,0]), 2, wf([2,1,-2,0,2]), 'qg', 0.6_dp)
        call fispin('W g;g IS g>g', P3, wf([0,1,-2,0]), 4, wf([0,1,-2,0,0]), 'gg', 0.55_dp)
        cycle
     endif
     call soft('q;qgg ', P3, [2,2,0,0], 3, 4)
     call soft('q;qgg ', P3, [2,2,0,0], 2, 3)
     call soft('qb;qbgg', P3, [-1,-1,0,0], 4, 2)
     call soft('g;qqbg', P3, [0,2,-2,0], 2, 4)
     call soft('g;qqbg', P3, [0,1,-1,0], 3, 2)
     call soft('q;qQQb', P3, [1,1,2,-2], 3, 4)
     call soft('q;qqqb', P3, [2,2,2,-2], 2, 4)
     call soft('qb;qbQQb', P3, [-2,-2,1,-1], 3, 2)
     call soft('qb;qbqbq', P3, [-1,-1,-1,1], 4, 3)
     call ffspin('q;qgg g>gg', P3, [2,2,0,0], 3, 2, [2,2,0,0,0], 'ggg', 0.3_dp)
     call ffspin('g;qqbg g>gg', P3, [0,1,-1,0], 4, 3, [0,1,-1,0,0], 'ggg', 0.6_dp)
     call ffspin('q;qgg g>qqb', P3, [2,2,0,0], 4, 3, [2,2,0,2,-2], 'gqq', 0.35_dp)
     call ffspin('qb;qbgg g>qqb', P3, [-1,-1,0,0], 3, 2, [-1,-1,1,0,-1], 'gqq', 0.7_dp)
     call ffspin('g;qqbg g>qqb', P3, [0,2,-2,0], 4, 2, [0,2,-2,1,-1], 'gqq', 0.5_dp)
     call fispin('q;g... IS q>g', P3, [0,2,-2,0], 2, [2,2,-2,0,2], 'qg', 0.6_dp)
     call fispin('qb;g.. IS qb>g', P3, [0,1,-1,0], 3, [-1,1,-1,0,-1], 'qg', 0.4_dp)
     call fispin('g;g... IS g>g', P3, [0,2,-2,0], 4, [0,2,-2,0,0], 'gg', 0.55_dp)
  enddo
  write(*,'(a,5es10.2)') ' largest deviations (me31, colour cons., pol. sum, soft, spin):', worst
  if (worst(1) > 1d-12 .or. worst(2) > 1d-10 .or. worst(3) > 1d-10 .or. worst(4) > 1d-3 .or. worst(5) > 2d-3) stop 1

contains

  logical function resolved3(P3, frac)
    real(dp), intent(in) :: P3(4,7), frac
    real(dp) :: W2, smin
    integer :: i, j
    W2 = mdot(P3(:,1) + P3(:,5), P3(:,1) + P3(:,5))
    smin = huge(1.0_dp)
    do i = 1, 4
       do j = i + 1, 4
          smin = min(smin, 2*abs(mdot(P3(:,i), P3(:,j))))
       enddo
    enddo
    resolved3 = smin > frac*W2
  end function resolved3

  ! reproducible random points: fixed seed, optionally shifted by the first
  ! command-line argument
  subroutine seed_rng()
    integer :: nseed, off, m
    integer, allocatable :: seed(:)
    character(32) :: arg
    off = 0
    if (command_argument_count() > 0) then
       call get_command_argument(1, arg); read(arg, *) off
    endif
    call random_seed(size=nseed); allocate(seed(nseed))
    seed = [(4711 + 7919*off + 104729*m, m = 1, nseed)]
    call random_seed(put=seed)
  end subroutine seed_rng

  subroutine check_flavours(P3)
    real(dp), intent(in) :: P3(4,7)
    integer :: f, Q, fl(4), i, ig, list(4,1000), nl, j
    real(dp) :: m, msq, cc(4,4), Ci, e1(4), e2(4), mv1, mv2, cv1(4,4), cv2(4,4)
    nl = 0
    do f = -5, 5
       if (ew31_mode == 2) then
          ! W exchange: every flavour multiset (zero ones skipped below)
          call cc_list(f, list, nl)
          cycle
       endif
       if (f == 0) then
          do Q = 1, 5
             nl = nl + 1; list(:,nl) = [0, Q, -Q, 0]
          enddo
       else
          nl = nl + 1; list(:,nl) = [f, f, 0, 0]
          do Q = 1, 5
             nl = nl + 1; list(:,nl) = [f, f, sign(Q, f), -sign(Q, f)]
          enddo
       endif
    enddo
    do j = 1, nl
       fl = list(:,j)
       call me31_tree(P3, fl, m)
       call born31_cc(P3, fl, msq, cc)
       if (m == 0 .and. msq == 0) cycle
       worst(1) = max(worst(1), abs(msq/m - 1))
       do i = 1, 4
          Ci = merge(CA, CF, fl(i) == 0)
          worst(2) = max(worst(2), abs(sum(cc(i,:)) + Ci*msq)/msq)
       enddo
       do ig = 1, 4
          if (fl(ig) /= 0) cycle
          call transverse(P3(:,ig), P3(:,merge(2, 1, ig == 1)), e1, e2)
          call born31_sc(P3, fl, ig, e1, mv1, cv1)
          call born31_sc(P3, fl, ig, e2, mv2, cv2)
          worst(3) = max(worst(3), abs((mv1 + mv2)/msq - 1), maxval(abs(cv1 + cv2 - cc))/msq)
       enddo
    enddo
  end subroutine check_flavours

  ! W exchange: every flavour multiset of the incoming parton f
  subroutine cc_list(f, list, nl)
    integer, intent(in) :: f
    integer, intent(inout) :: list(:,:), nl
    integer :: a, b, c, sg
    if (f == 0) then
       do a = 1, 5
          do b = 1, 5
             nl = nl + 1; list(:,nl) = [0, a, -b, 0]
          enddo
       enddo
       return
    endif
    sg = sign(1, f)
    do a = 1, 5
       nl = nl + 1; list(:,nl) = [f, sg*a, 0, 0]
       do b = a, 5
          do c = 1, 5
             nl = nl + 1; list(:,nl) = [f, sg*a, sg*b, -sg*c]
          enddo
       enddo
    enddo
  end subroutine cc_list

  subroutine soft(name, P3, fl3, i, k)
    character(*), intent(in) :: name
    real(dp), intent(in) :: P3(4,7)
    integer, intent(in) :: fl3(4), i, k
    real(dp) :: P4(4,8), pt(4), pk(4), e1(4), e2(4), kp(4), lam, z, y, m4, msq, cc(4,4), eik, q(4), rat
    integer :: il, a, b
    call born31_cc(P3, fl3, msq, cc)
    pt = P3(:,i); pk = P3(:,k)
    call perp(pt, pk, e1, e2)
    write(*,'(3x,a,i2,i2,a)', advance='no') 'soft  ' // name // ' (', i, k, ')'
    do il = 1, 4
       lam = 10.0_dp**(-il - 1)
       z = 1 - 0.7_dp*lam; y = 1.3_dp*lam
       kp = sqrt(z*(1 - z)*y*2*mdot(pt, pk))*(0.6_dp*e1 + 0.8_dp*e2)
       call fill41(P3, P4)
       P4(:,i) = z*pt + (1 - z)*y*pk + kp
       P4(:,5) = (1 - z)*pt + z*y*pk - kp
       P4(:,k) = (1 - y)*pk
       call me41_tree(P4, [fl3, 0], m4)
       q = P4(:,5)
       eik = 0
       do a = 1, 4
          do b = 1, 4
             if (a /= b) eik = eik + mdot(P4(:,a), P4(:,b))/(mdot(P4(:,a), q)*mdot(P4(:,b), q))*cc(a,b)
          enddo
       enddo
       rat = m4/(-8*pi**2*eik)
       write(*,'(f13.8)', advance='no') rat
    enddo
    write(*,*)
    worst(4) = max(worst(4), abs(rat - 1))
  end subroutine soft

  ! final-state gluon i of the Born split into slot i (fraction z) and slot 5,
  ! spectator k, at three azimuths
  subroutine ffspin(name, P3, fl3, i, k, fl4, kind, z)
    character(*), intent(in) :: name, kind
    real(dp), intent(in) :: P3(4,7), z
    integer, intent(in) :: fl3(4), i, k, fl4(5)
    real(dp) :: P4(4,8), pt(4), pk(4), e1(4), e2(4), kp(4), y, m4, msq, cc(4,4), mv, cv(4,4), sij, pred, rat, phi
    integer :: iphi, iy
    call born31_cc(P3, fl3, msq, cc)
    pt = P3(:,i); pk = P3(:,k)
    call perp(pt, pk, e1, e2)
    write(*,'(3x,a24)', advance='no') 'FF ' // name
    do iy = 1, 2
    y = merge(1d-8, 1d-10, iy == 1)
    do iphi = 0, 2
       phi = 0.4_dp + 1.1_dp*iphi
       kp = sqrt(z*(1 - z)*y*2*mdot(pt, pk))*(cos(phi)*e1 + sin(phi)*e2)
       call fill41(P3, P4)
       P4(:,i) = z*pt + (1 - z)*y*pk + kp
       P4(:,5) = (1 - z)*pt + z*y*pk - kp
       P4(:,k) = (1 - y)*pk
       call me41_tree(P4, fl4, m4)
       call born31_sc(P3, fl3, i, kp, mv, cv)
       sij = 2*mdot(P4(:,i), P4(:,5))
       select case (kind)
       case ('ggg'); pred = 2*CA*((z/(1 - z) + (1 - z)/z)*msq - 2*z*(1 - z)*mv/mdot(kp, kp))
       case ('gqq'); pred = TR*(msq + 4*z*(1 - z)*mv/mdot(kp, kp))
       end select
       rat = m4/(16*pi**2/sij*pred)
       write(*,'(f13.8)', advance='no') rat
       if (iy == 2) worst(5) = max(worst(5), abs(rat - 1))
    enddo
    enddo
    write(*,*)
  end subroutine ffspin

  ! initial state: incoming a (flavour fl4(1)) -> Born gluon (slot 1, fraction
  ! x) + slot 5, spectator k, three azimuths
  subroutine fispin(name, P3, fl3, k, fl4, kind, x)
    character(*), intent(in) :: name, kind
    real(dp), intent(in) :: P3(4,7), x
    integer, intent(in) :: fl3(4), k, fl4(5)
    real(dp) :: P4(4,8), pa(4), pk(4), e1(4), e2(4), kp(4), u, m4, msq, cc(4,4), mv, cv(4,4), pred, rat, phi
    integer :: iphi, iu
    call born31_cc(P3, fl3, msq, cc)
    pa = P3(:,1); pk = P3(:,k)
    call perp(pa, pk, e1, e2)
    write(*,'(3x,a24)', advance='no') 'IF ' // name
    do iu = 1, 2
    u = merge(1d-8, 1d-10, iu == 1)
    do iphi = 0, 2
       phi = 0.4_dp + 1.1_dp*iphi
       kp = sqrt(u*(1 - u)*(1 - x)/x*2*mdot(pa, pk))*(cos(phi)*e1 + sin(phi)*e2)
       call fill41(P3, P4)
       P4(:,1) = pa/x
       P4(:,5) = (1 - u)*(1 - x)/x*pa + u*pk + kp
       P4(:,k) = u*(1 - x)/x*pa + (1 - u)*pk - kp
       call me41_tree(P4, fl4, m4)
       call born31_sc(P3, fl3, 1, kp, mv, cv)
       select case (kind)
       case ('qg'); pred = CF*(x*msq - 4*(1 - x)/x*mv/mdot(kp, kp))
       case ('gg'); pred = 2*CA*((x/(1 - x) + x*(1 - x))*msq - 2*(1 - x)/x*mv/mdot(kp, kp))
       end select
       rat = m4/(16*pi**2/(x*2*mdot(P4(:,1), P4(:,5)))*pred)
       write(*,'(f13.8)', advance='no') rat
       if (iu == 2) worst(5) = max(worst(5), abs(rat - 1))
    enddo
    enddo
    write(*,*)
  end subroutine fispin

  ! unit vectors transverse to the light-like g and the reference ref
  subroutine transverse(g, ref, e1, e2)
    real(dp), intent(in) :: g(4), ref(4)
    real(dp), intent(out) :: e1(4), e2(4)
    call perp(g, ref, e1, e2)
  end subroutine transverse

  ! 3+1 layout (1-4 partons, 5 q, 6, 7 leptons) -> 4+1 (1-5, 6 q, 7, 8)
  subroutine fill41(P3, P4)
    real(dp), intent(in) :: P3(4,7)
    real(dp), intent(out) :: P4(4,8)
    P4 = 0
    P4(:,1:4) = P3(:,1:4); P4(:,6:8) = P3(:,5:7)
  end subroutine fill41

  ! two unit spacelike vectors orthogonal to the light-like a and b and to
  ! each other
  subroutine perp(a, b, e1, e2)
    real(dp), intent(in) :: a(4), b(4)
    real(dp), intent(out) :: e1(4), e2(4)
    real(dp) :: r1(4), r2(4)
    r1 = [1.0_dp, 0.3_dp, -0.2_dp, 0.0_dp]; r2 = [-0.4_dp, 1.0_dp, 0.5_dp, 0.0_dp]
    e1 = r1 - mdot(r1, b)/mdot(a, b)*a - mdot(r1, a)/mdot(a, b)*b
    e1 = e1/sqrt(-mdot(e1, e1))
    e2 = r2 - mdot(r2, b)/mdot(a, b)*a - mdot(r2, a)/mdot(a, b)*b
    e2 = e2 + mdot(e2, e1)*e1
    e2 = e2/sqrt(-mdot(e2, e2))
  end subroutine perp

  pure real(dp) function mdot(a, b)
    real(dp), intent(in) :: a(4), b(4)
    mdot = a(4)*b(4) - a(1)*b(1) - a(2)*b(2) - a(3)*b(3)
  end function mdot

  subroutine dis_point(r, P)
    real(dp), intent(in) :: r(8)
    real(dp), intent(out) :: P(4,7)
    real(dp) :: Q, y, xi, E, E1, W, beta, m23, pa(4), pb(4), pc(4), X(4), d(3)
    Q = 10 + 90*r(1); y = 0.1_dp + 0.8_dp*r(2); xi = 0.05_dp + 0.8_dp*r(3)
    E = Q/2; E1 = E/xi; W = sqrt(4*E*E1 - 4*E**2); beta = (E1 - 2*E)/E1
    P = 0
    P(:,1) = [0.0_dp, 0.0_dp, E1, E1]
    P(:,5) = [0.0_dp, 0.0_dp, -2*E, 0.0_dp]
    P(:,6) = [E/y*2*sqrt(1-y), 0.0_dp, -E, E/y*(2-y)]
    P(:,7) = [E/y*2*sqrt(1-y), 0.0_dp,  E, E/y*(2-y)]
    m23 = W*(0.1_dp + 0.8_dp*r(4))
    d = unitv(2*r(5) - 1, 6.283185307_dp*r(6))
    pa = [(W*W - m23*m23)/(2*W)*d, (W*W - m23*m23)/(2*W)]
    X = [-pa(1:3), W - pa(4)]
    call twobody(X, unitv(2*r(7) - 1, 6.283185307_dp*r(8)), pb, pc)
    P(:,2) = boostz(pa, beta); P(:,3) = boostz(pb, beta); P(:,4) = boostz(pc, beta)
  end subroutine dis_point

  function unitv(cth, phi) result(d)
    real(dp), intent(in) :: cth, phi
    real(dp) :: d(3)
    d = [sqrt(1-cth**2)*cos(phi), sqrt(1-cth**2)*sin(phi), cth]
  end function unitv

  subroutine twobody(Ptot, dir, pa, pb)
    real(dp), intent(in)  :: Ptot(4), dir(3)
    real(dp), intent(out) :: pa(4), pb(4)
    real(dp) :: Mm, h
    Mm = sqrt(Ptot(4)**2 - sum(Ptot(1:3)**2)); h = Mm/2
    pa = boostv([h*dir, h], Ptot); pb = boostv([-h*dir, h], Ptot)
  end subroutine twobody

  function boostv(p, Pt) result(q)
    real(dp), intent(in) :: p(4), Pt(4)
    real(dp) :: q(4), Mm, bp, f
    Mm = sqrt(Pt(4)**2 - sum(Pt(1:3)**2))
    bp = dot_product(p(1:3), Pt(1:3))
    q(4) = (Pt(4)*p(4) + bp)/Mm
    f = (bp/(Pt(4) + Mm) + p(4))/Mm
    q(1:3) = p(1:3) + f*Pt(1:3)
  end function boostv

  function boostz(p, beta) result(q)
    real(dp), intent(in) :: p(4), beta
    real(dp) :: q(4), g
    g = 1/sqrt(1 - beta**2)
    q = p
    q(3) = g*(p(3) + beta*p(4)); q(4) = g*(p(4) + beta*p(3))
  end function boostz
  ! W+ (e+, nu): the W- channel lists with up and down exchanged (u <-> d,
  ! c <-> s), 9 Oct
  function wf(f) result(g)
    integer, intent(in) :: f(:)
    integer :: g(size(f)), i
    integer, parameter :: sw(0:5) = [0, 2, 1, 4, 3, 5]
    g = f
    if (ew31_lepton /= 1 .and. ew31_lepton /= 2) return
    do i = 1, size(f)
       g(i) = sign(sw(abs(f(i))), f(i))
    enddo
  end function wf
end program harness_born31
