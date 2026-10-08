! me41 (4+1 trees) against me31 (3+1 trees, validated against DISENT's
! MATFOR) in single collinear limits; photon exchange, or with the
! environment EW31 = "mode lepton" the couplings of ew31 (e.g. "1 0"):
!  - final state i || j (spectator k, CS FF map with y -> 0):
!      me41 -> 16 pi^2/s_ij P(z) me31,
!    P = CF (1+z^2)/(1-z) (q -> q(z) g), 2 CA [z/(1-z) + (1-z)/z + z(1-z)]
!    (g -> g(z) g), TR [1 - 2 z(1-z)] (g -> q(z) qbar); the azimuthal
!    correlations of the gluon splittings cancel in the average over
!    phi and phi + pi/2;
!  - initial state a -> a~(x) + i (spectator k final, CS IF map with u -> 0):
!      me41 -> (1/x) 16 pi^2/(2 pa.pi) P(x) me31(a~)
!    (spin- and colour-averaged matrix elements, so the averaged
!    Altarelli-Parisi kernels): q -> q: CF (1+x^2)/(1-x); g -> q:
!    TR [1 - 2x(1-x)].
! Ratios printed for decreasing y (u); they must tend to 1 (deviation
! about linear in sqrt(y) or y).
program harness_lim41
  use me31
  use me41
  use ew31
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: CF = 4.0_dp/3, CA = 3, TR = 0.5_dp
  real(dp) :: P3(4,7), r(8), worst
  integer :: ipt
  character(32) :: arg
  call get_environment_variable('EW31', arg)
  if (len_trim(arg) > 0) read(arg, *) ew31_mode, ew31_lepton
  call seed_rng()
  worst = 0
  do ipt = 1, 3
     ! a resolved 3+1 point (all 2 p_i.p_j above 2% of W^2); near-singular
     ! Borns converge slowly
     do
        call random_number(r)
        call dis_point(r, P3)
        if (resolved3(P3, 0.02_dp)) exit
     enddo
     write(*,'(a,i2)') ' point', ipt
     ! q g g -> q g g g: final q || g, g || g; initial q -> q g
     call ff('q;qgg   q||g  ', P3, [2,2,0,0], 2, 3, [2,2,0,0,0], 'qqg', 0.3_dp)
     call ff('q;qgg   g||g  ', P3, [2,2,0,0], 3, 2, [2,2,0,0,0], 'ggg', 0.4_dp)
     call ff('qb;qbgg q||g  ', P3, [-1,-1,0,0], 2, 4, [-1,-1,0,0,0], 'qqg', 0.6_dp)
     call fi('q;qgg   IS q>q', P3, [2,2,0,0], 2, [2,2,0,0,0], 'qq', 0.7_dp)
     ! g -> q qbar g g
     call ff('g;qqbg  q||g  ', P3, [0,2,-2,0], 2, 3, [0,2,-2,0,0], 'qqg', 0.35_dp)
     call ff('g;qqbg  g||g  ', P3, [0,2,-2,0], 4, 2, [0,2,-2,0,0], 'ggg', 0.55_dp)
     call fi('g;gqqb  IS g>q', P3, [2,2,0,0], 3, [0,2,0,0,-2], 'gq', 0.6_dp)
     ! q -> q Q Qbar g
     call ff('q;qQQb  q||g  ', P3, [1,1,2,-2], 2, 4, [1,1,2,-2,0], 'qqg', 0.4_dp)
     call ff('q;qQQb  Qb||g ', P3, [1,1,2,-2], 4, 2, [1,1,2,-2,0], 'qqg', 0.5_dp)
     call ff('q;qgg>qQQb g>QQb', P3, [1,1,0,0], 3, 2, [1,1,2,0,-2], 'gqq', 0.45_dp)
     call fi('q;qQQb  IS q>q', P3, [1,1,2,-2], 2, [1,1,2,-2,0], 'qq', 0.5_dp)
     call ff('qb;qbQQb q||g ', P3, [-2,-2,1,-1], 3, 2, [-2,-2,1,-1,0], 'qqg', 0.3_dp)
     ! identical quarks
     call ff('q;qqqb  q||g  ', P3, [2,2,2,-2], 3, 2, [2,2,2,-2,0], 'qqg', 0.4_dp)
     call ff('q;qqqb  qb||g ', P3, [2,2,2,-2], 4, 3, [2,2,2,-2,0], 'qqg', 0.6_dp)
     call ff('q;qgg>qqqb g>qqb', P3, [2,2,0,0], 3, 2, [2,2,2,0,-2], 'gqq', 0.5_dp)
     call fi('q;qqqb  IS q>q', P3, [2,2,2,-2], 3, [2,2,2,-2,0], 'qq', 0.4_dp)
     ! g -> q qbar Q Qbar (from q Q Qbar by IS g -> q, and from g q qbar g by
     ! g -> Q Qbar)
     call fi('g;qQQb  IS g>q', P3, [1,1,2,-2], 2, [0,1,2,-2,-1], 'gq', 0.5_dp)
     call ff('g;qqbg>qqbQQb', P3, [0,1,-1,0], 4, 2, [0,1,-1,2,-2], 'gqq', 0.4_dp)
     call fi('g;qqqb  IS g>q', P3, [2,2,2,-2], 4, [0,2,2,-2,-2], 'gq', 0.3_dp)
     call ff('g;qqbg>qqbqqb', P3, [0,1,-1,0], 4, 3, [0,1,-1,1,-1], 'gqq', 0.6_dp)
  enddo
  write(*,'(a,es10.2)') ' largest |ratio - 1| at the smallest y:', worst
  if (worst > 3d-3) stop 1   ! g -> Q Qbar limits averaged over two azimuths converge like sqrt(y)

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

  ! final-state splitting of parton i of P3 (flavour fl3(i)) into i (momentum
  ! fraction z, flavour fl4 slot i) and j (new slot 5, flavour fl4(5)), with
  ! spectator k
  subroutine ff(name, P3, fl3, i, k, fl4, kind, z)
    character(*), intent(in) :: name, kind
    real(dp), intent(in) :: P3(4,7), z
    integer, intent(in) :: fl3(4), i, k, fl4(5)
    real(dp) :: P4(4,8), pt(4), pk(4), e1(4), e2(4), kp(4), y, m3, m4, sij, Pk_, rat, avr
    integer :: iy, iphi, nzero
    call me31_tree(P3, fl3, m3)
    pt = P3(:,i); pk = P3(:,k)
    call perp(pt, pk, e1, e2)
    write(*,'(3x,a20)', advance='no') name
    do iy = 1, 4
       y = 10.0_dp**(-2*iy - 2)
       avr = 0; nzero = 0
       do iphi = 0, 1
          kp = sqrt(z*(1 - z)*y*2*mdot(pt, pk))*merge(e1, e2, iphi == 0)
          call fill41(P3, P4)
          P4(:,i) = z*pt + (1 - z)*y*pk + kp
          P4(:,5) = (1 - z)*pt + z*y*pk - kp
          P4(:,k) = (1 - y)*pk
          call me41_tree(P4, fl4, m4)
          sij = 2*mdot(P4(:,i), P4(:,5))
          select case (kind)
          case ('qqg'); Pk_ = CF*(1 + z**2)/(1 - z)
          case ('ggg'); Pk_ = 2*CA*(z/(1 - z) + (1 - z)/z + z*(1 - z))
          case ('gqq'); Pk_ = TR*(1 - 2*z*(1 - z))
          end select
          avr = avr + 0.5_dp*m4
          rat = 16*pi**2/sij*Pk_*m3
       enddo
       rat = avr/rat
       write(*,'(f13.8)', advance='no') rat
    enddo
    write(*,*)
    worst = max(worst, abs(rat - 1))
  end subroutine ff

  ! initial-state splitting: incoming a (flavour fl4(1)) -> a~ (fl3(1),
  ! momentum x pa) + i (new slot 5, flavour fl4(5)), spectator k final
  subroutine fi(name, P3, fl3, k, fl4, kind, x)
    character(*), intent(in) :: name, kind
    real(dp), intent(in) :: P3(4,7), x
    integer, intent(in) :: fl3(4), k, fl4(5)
    real(dp) :: P4(4,8), pa(4), pk(4), e1(4), e2(4), kp(4), u, m3, m4, Pk_, rat, avr
    integer :: iu, iphi
    call me31_tree(P3, fl3, m3)
    pa = P3(:,1); pk = P3(:,k)
    call perp(pa, pk, e1, e2)
    write(*,'(3x,a20)', advance='no') name
    do iu = 1, 4
       u = 10.0_dp**(-2*iu - 2)
       avr = 0
       do iphi = 0, 1
          kp = sqrt(u*(1 - u)*(1 - x)/x*2*mdot(pa, pk))*merge(e1, e2, iphi == 0)
          call fill41(P3, P4)
          P4(:,1) = pa/x
          P4(:,5) = (1 - u)*(1 - x)/x*pa + u*pk + kp
          P4(:,k) = u*(1 - x)/x*pa + (1 - u)*pk - kp
          call me41_tree(P4, fl4, m4)
          avr = avr + 0.5_dp*m4
       enddo
       select case (kind)
       case ('qq'); Pk_ = CF*(1 + x**2)/(1 - x)
       case ('gq'); Pk_ = TR*(1 - 2*x*(1 - x))
       end select
       rat = avr/(16*pi**2/(x*2*mdot(P4(:,1), P4(:,5)))*Pk_*m3)
       write(*,'(f13.8)', advance='no') rat
    enddo
    write(*,*)
    worst = max(worst, abs(rat - 1))
  end subroutine fi

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
end program harness_lim41
