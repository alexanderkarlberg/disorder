! dip41 (Catani-Seymour dipoles for DIS 4+1 -> 3+1) against me41: the
! ratio me41 / sum of dipoles must tend to 1 in every single-unresolved
! limit (as tests/test_subtraction.f90 for DISENT):
!  - collinear final-state pairs (CS FF map from a 3+1 point, y -> 0, at a
!    fixed azimuth, no averaging), all splittings q -> q g, g -> g g,
!    g -> q qbar;
!  - collinear initial-state splittings (IF map, u -> 0, fixed azimuth):
!    q -> q g, q -> g q, g -> q qbar, g -> g g;
!  - soft gluons (FF map with 1 - z, y ~ lambda -> 0).
! Each 4+1 point is built from a random 3+1 point; ratios printed for four
! decreasing values of the limit parameter. Dipoles whose mapped Born is
! itself unresolved (some 2 p_i.p_j below 1e-3 W^2) are left out: there the
! Born still contains the limit pair, and in a 3+1 calculation the jet
! function of the mapped Born removes them (as test_subtraction does for
! DISENT).
! Couplings: photon, or with the environment EW31 = "mode lepton" those of
! ew31 (e.g. "1 0": photon + Z, e-).
program harness_dip41
  use me41
  use dip41
  use ew31
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
  real(dp) :: P3(4,7), r(8), worst
  integer :: ipt
  character(32) :: ewarg
  call get_environment_variable('EW31', ewarg)
  if (len_trim(ewarg) > 0) read(ewarg, *) ew31_mode, ew31_lepton
  call seed_rng()
  worst = 0
  do ipt = 1, 3
     ! a resolved 3+1 point (all 2 p_i.p_j above 2% of W^2), as a jet
     ! function would require; near-singular Borns converge slowly
     do
        call random_number(r)
        call dis_point(r, P3)
        if (resolved3(P3, 0.02_dp)) exit
     enddo
     write(*,'(a,i2)') ' point', ipt
     ! final-state collinear: Born slot i split into i and 5, spectator k
     call ffc('q;qggg   q||g ', P3, 2, 3, [2,2,0,0,0], 0.3_dp, 0.0_dp)
     call ffc('q;qggg   g||g ', P3, 3, 4, [2,2,0,0,0], 0.4_dp, 0.0_dp)
     call ffc('qb;qbggg qb||g', P3, 2, 1, [-1,-1,0,0,0], 0.6_dp, 0.0_dp)
     call ffc('g;qqbgg  q||g ', P3, 2, 4, [0,2,-2,0,0], 0.35_dp, 0.0_dp)
     call ffc('g;qqbgg  g||g ', P3, 4, 1, [0,2,-2,0,0], 0.55_dp, 0.0_dp)
     call ffc('q;qgQQb  Q||Qb', P3, 4, 2, [1,1,0,2,-2], 0.45_dp, 1.0_dp)
     call ffc('q;qQQbg  q||g ', P3, 2, 4, [1,1,2,-2,0], 0.4_dp, 0.0_dp)
     call ffc('q;qQQbg  Qb||g', P3, 4, 1, [1,1,2,-2,0], 0.5_dp, 0.0_dp)
     call ffc('q;qqqbg  q||g ', P3, 3, 2, [2,2,2,-2,0], 0.4_dp, 0.0_dp)
     call ffc('g;qqbQQb Q||Qb', P3, 4, 2, [0,1,-1,2,-2], 0.4_dp, 1.0_dp)
     call ffc('g;qqbqqb q||qb', P3, 4, 3, [0,1,-1,1,-1], 0.6_dp, 1.0_dp)
     ! initial-state collinear: incoming a -> Born slot 1 + slot 5, spectator k
     call ifc('q;qggg   IS q>q', P3, 2, [2,2,0,0,0], 0.7_dp)
     call ifc('q;qQQbg  IS q>q', P3, 3, [1,1,2,-2,0], 0.5_dp)
     call ifc('q;qqbgq  IS q>g', P3, 2, [2,2,-2,0,2], 0.6_dp)
     call ifc('qb;..    IS qb>g', P3, 4, [-1,1,-1,0,-1], 0.45_dp)
     call ifc('g;qqbgg  IS g>q', P3, 3, [0,2,0,0,-2], 0.6_dp)
     call ifc('g;qQQbqb IS g>q', P3, 2, [0,1,2,-2,-1], 0.5_dp)
     call ifc('g;qqbgg  IS g>g', P3, 4, [0,2,-2,0,0], 0.55_dp)
     ! soft gluon (slot 5) between i and k
     call softg('q;qggg  ', P3, [2,2,0,0], 3, 4)
     call softg('q;qggg  ', P3, [2,2,0,0], 2, 3)
     call softg('qb;qbggg', P3, [-1,-1,0,0], 4, 2)
     call softg('g;qqbgg ', P3, [0,2,-2,0], 2, 4)
     call softg('q;qQQbg ', P3, [1,1,2,-2], 3, 4)
     call softg('q;qqqbg ', P3, [2,2,2,-2], 2, 4)
     call softg('qb;qbQQb', P3, [-2,-2,1,-1], 3, 2)
  enddo
  write(*,'(a,es10.2)') ' largest |me41/sum(dipoles) - 1| at the smallest parameter:', worst
  if (worst > 2d-3) stop 1

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

  real(dp) function dsum(P4, fl4)
    real(dp), intent(in) :: P4(4,8)
    integer, intent(in) :: fl4(5)
    real(dp) :: P3(4,7,dip41_max), val(dip41_max)
    integer :: fl3(4,dip41_max), nd
    integer :: id
    call dip41_list(P4, fl4, nd, P3, fl3, val)
    dsum = 0
    do id = 1, nd
       if (resolved(P3(:,:,id))) dsum = dsum + val(id)
    enddo
  end function dsum

  ! all invariants 2 p_i.p_j of the Born partons above 1e-3 W^2
  logical function resolved(P3)
    real(dp), intent(in) :: P3(4,7)
    real(dp) :: W2, smin
    integer :: i, j
    W2 = mdot(P3(:,1) + P3(:,5), P3(:,1) + P3(:,5))
    smin = huge(1.0_dp)
    do i = 1, 4
       do j = i + 1, 4
          smin = min(smin, 2*abs(mdot(P3(:,i), P3(:,j))))
       enddo
    enddo
    resolved = smin > 1d-3*W2
  end function resolved

  ! Born slot i (of the 3+1 point, flavours from fl4 with slot 5 merged
  ! into i) split into i (fraction z) and 5, spectator k (k = 1: the
  ! incoming parton, FI-like map with the final-state recoil instead: use
  ! k = a final slot otherwise); phi offset dphi
  subroutine ffc(name, P3, i, k, fl4, z, dphi)
    character(*), intent(in) :: name
    real(dp), intent(in) :: P3(4,7), z, dphi
    integer, intent(in) :: i, k, fl4(5)
    real(dp) :: P4(4,8), pt(4), pk(4), e1(4), e2(4), kp(4), y, m4, rat, phi
    integer :: iy, kk
    kk = k
    if (kk == 1) kk = merge(3, 2, i == 2)
    pt = P3(:,i); pk = P3(:,kk)
    call perp(pt, pk, e1, e2)
    phi = 0.7_dp + dphi
    write(*,'(3x,a22)', advance='no') 'FF ' // name
    do iy = 1, 4
       y = 10.0_dp**(-2*iy - 2)
       kp = sqrt(z*(1 - z)*y*2*mdot(pt, pk))*(cos(phi)*e1 + sin(phi)*e2)
       call fill41(P3, P4)
       P4(:,i) = z*pt + (1 - z)*y*pk + kp
       P4(:,5) = (1 - z)*pt + z*y*pk - kp
       P4(:,kk) = (1 - y)*pk
       call me41_tree(P4, fl4, m4)
       rat = m4/dsum(P4, fl4)
       write(*,'(f13.8)', advance='no') rat
    enddo
    write(*,*)
    worst = max(worst, abs(rat - 1))
  end subroutine ffc

  subroutine ifc(name, P3, k, fl4, x)
    character(*), intent(in) :: name
    real(dp), intent(in) :: P3(4,7), x
    integer, intent(in) :: k, fl4(5)
    real(dp) :: P4(4,8), pa(4), pk(4), e1(4), e2(4), kp(4), u, m4, rat
    integer :: iu
    pa = P3(:,1); pk = P3(:,k)
    call perp(pa, pk, e1, e2)
    write(*,'(3x,a22)', advance='no') 'IF ' // name
    do iu = 1, 4
       u = 10.0_dp**(-2*iu - 2)
       kp = sqrt(u*(1 - u)*(1 - x)/x*2*mdot(pa, pk))*(cos(0.9_dp)*e1 + sin(0.9_dp)*e2)
       call fill41(P3, P4)
       P4(:,1) = pa/x
       P4(:,5) = (1 - u)*(1 - x)/x*pa + u*pk + kp
       P4(:,k) = u*(1 - x)/x*pa + (1 - u)*pk - kp
       call me41_tree(P4, fl4, m4)
       rat = m4/dsum(P4, fl4)
       write(*,'(f13.8)', advance='no') rat
    enddo
    write(*,*)
    worst = max(worst, abs(rat - 1))
  end subroutine ifc

  subroutine softg(name, P3, fl3, i, k)
    character(*), intent(in) :: name
    real(dp), intent(in) :: P3(4,7)
    integer, intent(in) :: fl3(4), i, k
    real(dp) :: P4(4,8), pt(4), pk(4), e1(4), e2(4), kp(4), lam, z, y, m4, rat
    integer :: il
    pt = P3(:,i); pk = P3(:,k)
    call perp(pt, pk, e1, e2)
    write(*,'(3x,a,i2,i2,a)', advance='no') 'soft ' // name // ' (', i, k, ')'
    do il = 1, 4
       lam = 10.0_dp**(-il - 1)
       z = 1 - 0.7_dp*lam; y = 1.3_dp*lam
       kp = sqrt(z*(1 - z)*y*2*mdot(pt, pk))*(0.6_dp*e1 + 0.8_dp*e2)
       call fill41(P3, P4)
       P4(:,i) = z*pt + (1 - z)*y*pk + kp
       P4(:,5) = (1 - z)*pt + z*y*pk - kp
       P4(:,k) = (1 - y)*pk
       call me41_tree(P4, [fl3, 0], m4)
       rat = m4/dsum(P4, [fl3, 0])
       write(*,'(f13.8)', advance='no') rat
    enddo
    write(*,*)
    worst = max(worst, abs(rat - 1))
  end subroutine softg

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
end program harness_dip41
