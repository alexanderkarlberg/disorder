! Tests of iop31 (integrated CS dipoles for DIS 3+1):
!  - the poles of <I(eps)> cancel those of virt31_ren (all channels, two
!    scales);
!  - P: d lp/d ln muf^2 = -|M_b|^2 (colour conservation, CS eq. 8.42);
!  - V + I finite part printed per channel (for the record).
program harness_iop31
  use virt31
  use iop31
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp) :: P3(4,7), r(8), v(-2:0), t, iv(-2:0), mu2, worst(3), b, g, lp1, lp2
  integer :: ipt, ic, imu
  integer, parameter :: nfl = 9
  integer, parameter :: fls(4,nfl) = reshape([2,2,0,0, -1,-1,0,0, 0,1,-1,0, 0,2,-2,0, 1,1,2,-2, &
       & 2,2,2,-2, -1,-1,-2,2, -2,-2,-2,2, 1,2,1,-2], [4,nfl])
  worst = 0
  do ipt = 1, 4
     call random_number(r); call dis_point(r, P3)
     do imu = 1, 2
        mu2 = merge(50.0_dp, 800.0_dp, imu == 1)
        do ic = 1, nfl
           call virt31_ren(P3, fls(:,ic), mu2, v, t)
           call iop31_i(P3, fls(:,ic), mu2, iv)
           worst(1) = max(worst(1), abs(v(-2) + iv(-2))/t)
           worst(2) = max(worst(2), abs(v(-1) + iv(-1))/t)
           call iop31_kp(P3, fls(:,ic), mu2, b, g, lp1)
           call iop31_kp(P3, fls(:,ic), mu2*exp(0.1_dp), b, g, lp2)
           worst(3) = max(worst(3), abs((lp2 - lp1)/0.1_dp + b)/b)
           if (ipt == 1 .and. imu == 1) write(*,'(a,4i3,a,f12.6,a,2f12.6)') ' fl', fls(:,ic), &
                & '  (V + I)/tree finite:', (v(0) + iv(0))/t, '   g/b, lp/b:', g/b, lp1/b
        enddo
     enddo
  enddo
  write(*,'(a,3es10.2)') ' largest |V + I| poles (1/eps^2, 1/eps) per tree, dlp/dln muf^2 + b:', worst
  if (any(worst > 1d-10)) stop 1
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
end program harness_iop31
