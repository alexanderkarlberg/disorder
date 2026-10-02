! Harness: MCFM's crossed q qbar g g + gamma* tree (z2jetsq) against DISENT's
! four-parton matrix element MATFOR, photon exchange, incoming gluon
! (gamma* g -> q qbar g, no four-quark channel at this order). Both summed
! over the six assignments of (q, qbar, g) to the outgoing slots 2, 3, 4;
! the ratio must be the same at every point.
program harness31
  use types, only: dp
  use mod_parameters
  implicit none
  integer :: SCHEME_D, NF
  double precision :: CF, CA, TR, PI_D, PISQ, HF, CUTOFF_D, EQ(-6:6), SCALE_D
  common /COLFAC/ CF, CA, TR, PI_D, PISQ, HF, CUTOFF_D, EQ, SCALE_D, SCHEME_D, NF
  integer, parameter :: mxpart = 14
  real(dp) :: P(4,7), PP(4,7), M4(-6:6), dsum, msum, msq(2,2), pm(mxpart,4), r(8), rat(20)
  complex(dp) :: za(mxpart,mxpart), zb(mxpart,mxpart)
  integer :: perm(3,6), ip, ipt, i
  perm = reshape([2,3,4, 2,4,3, 3,2,4, 3,4,2, 4,2,3, 4,3,2], [3,6])

  call set_parameters()
  CF = 4.0_dp/3.0_dp; CA = 3; TR = 0.5_dp; NF = nflav
  PI_D = atan(1d0)*4; PISQ = PI_D**2; HF = 0.5_dp; CUTOFF_D = 1d-8
  SCHEME_D = 0; SCALE_D = 1
  EQ(0) = 0; EQ(1) = -1D0/3; EQ(2) = EQ(1) + 1
  do i = 1, 6
     if (i > 2) EQ(i) = EQ(i-2)
     EQ(-i) = -EQ(i)
  enddo

  do ipt = 1, 20
     call random_number(r)
     call dis_point(r, P)
     if (maxval(abs(P(:,1) + P(:,6) - P(:,7) - P(:,2) - P(:,3) - P(:,4))) > 1d-9*P(4,1)) stop 'momentum not conserved'
     dsum = 0; msum = 0
     do ip = 1, 6
        ! DISENT: relabel the outgoing partons
        PP = P
        PP(:,2) = P(:,perm(1,ip)); PP(:,3) = P(:,perm(2,ip)); PP(:,4) = P(:,perm(3,ip))
        call MATFOR(PP, M4)
        dsum = dsum + M4(0)
        ! MCFM, all momenta incoming: 1 = -gluon in, 2 = -electron in (the
        ! 'antilepton' slot), 3 = electron out, 4 = q, 5 = qbar, 6 = g out
        pm = 0
        pm(1,:) = -P(:,1); pm(2,:) = -P(:,6); pm(3,:) = P(:,7)
        pm(4,:) = P(:,perm(1,ip)); pm(5,:) = P(:,perm(2,ip)); pm(6,:) = P(:,perm(3,ip))
        call spinoru(6, pm, za, zb)
        call z2jetsq(4, 5, 3, 2, 1, 6, za, zb, msq)
        msum = msum + sum(msq)
     enddo
     rat(ipt) = dsum/msum
     write(*,'(a,i3,a,es14.6,a,es14.6,a,es20.12)') ' point', ipt, '  DISENT', dsum, '  MCFM', msum, '  ratio', rat(ipt)
  enddo
  write(*,'(a,es10.2)') ' max relative spread of the ratio:', (maxval(rat) - minval(rat))/abs(sum(rat)/20)

contains

  ! a DIS point in DISENT's layout (Breit-like frame from test_subtraction:
  ! incoming parton along +z, q along -z), three outgoing partons from a
  ! random massless three-body decay of W = p1 + q
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
    ! W rest frame: pa back to back with X (mass m23), X -> pb pc
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

  ! boost p (rest frame of Pt) to the frame where Pt has its given momentum
  function boostv(p, Pt) result(q)
    real(dp), intent(in) :: p(4), Pt(4)
    real(dp) :: q(4), Mm, bp, f
    Mm = sqrt(Pt(4)**2 - sum(Pt(1:3)**2))
    bp = dot_product(p(1:3), Pt(1:3))
    q(4) = (Pt(4)*p(4) + bp)/Mm
    f = (bp/(Pt(4) + Mm) + p(4))/Mm
    q(1:3) = p(1:3) + f*Pt(1:3)
  end function boostv

  ! boost along z with velocity beta (as test_subtraction: W rest frame -> Breit frame)
  function boostz(p, beta) result(q)
    real(dp), intent(in) :: p(4), beta
    real(dp) :: q(4), g
    g = 1/sqrt(1 - beta**2)
    q = p
    q(3) = g*(p(3) + beta*p(4)); q(4) = g*(p(4) + beta*p(3))
  end function boostz
end program harness31
