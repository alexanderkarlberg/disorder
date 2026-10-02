! me31 (MCFM crossed, absolute normalisation) against DISENT's MATFOR,
! photon exchange, all incoming partons -5..5, summed over the labellings of
! the outgoing partons with 1/2 for identical final-state pairs: ratio 1.
program time_me31
  use types, only: dp
  use mod_parameters
  use me31
  implicit none
  integer :: SCHEME_D, NF
  double precision :: CF, CA, TR, PI_D, PISQ, HF, CUTOFF_D, EQ(-6:6), SCALE_D
  common /COLFAC/ CF, CA, TR, PI_D, PISQ, HF, CUTOFF_D, EQ, SCALE_D, SCHEME_D, NF
  real(dp) :: P(4,7), PP(4,7), M4(-6:6), r(8), d(-5:5), m(-5:5), x, worst
  integer :: perm(3,6), ip, ipt, i, f, Q
  integer, parameter :: ncall = 20000
  real(dp) :: t0, t1, tsum
  integer :: n
  perm = reshape([2,3,4, 2,4,3, 3,2,4, 3,4,2, 4,2,3, 4,3,2], [3,6])
  call set_parameters()
  CF = 4.0_dp/3.0_dp; CA = 3; TR = 0.5_dp; NF = 5
  PI_D = atan(1d0)*4; PISQ = PI_D**2; HF = 0.5_dp; CUTOFF_D = 1d-8
  SCHEME_D = 0; SCALE_D = 1
  EQ(0) = 0; EQ(1) = -1D0/3; EQ(2) = EQ(1) + 1
  do i = 1, 6
     if (i > 2) EQ(i) = EQ(i-2)
     EQ(-i) = -EQ(i)
  enddo
  call random_number(r)
  call dis_point(r, P)
  call cpu_time(t0)
  tsum = 0
  do n = 1, ncall
     P(1:3,2) = P(1:3,2)*(1 + 1d-12*n)   ! defeat any caching by the compiler
     call MATFOR(P, M4); tsum = tsum + M4(1)
  enddo
  call cpu_time(t1)
  write(*,'(a,f8.2,a)') ' MATFOR (all flavours):  ', 1d6*(t1 - t0)/ncall, ' us per point'
  call cpu_time(t0)
  do n = 1, ncall
     do f = -5, 5
        if (f == 0) then
           do Q = 1, 5
              call me31_tree(P, [0, Q, -Q, 0], x); tsum = tsum + x
           enddo
        else
           call me31_tree(P, [f, f, 0, 0], x); tsum = tsum + x
           do Q = 1, 5
              if (Q == abs(f)) then
                 call me31_tree(P, [f, f, f, -f], x)
              else
                 call me31_tree(P, [f, f, sign(Q, f), -sign(Q, f)], x)
              endif
              tsum = tsum + x
           enddo
        endif
     enddo
  enddo
  call cpu_time(t1)
  write(*,'(a,f8.2,a,es10.2)') ' me31 (all 71 flavour assignments, one labelling):', 1d6*(t1 - t0)/ncall, ' us per point', tsum
contains

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
end program time_me31
