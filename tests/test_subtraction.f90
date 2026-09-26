!----------------------------------------------------------------------
! Unit test of the consistency of DISENT's four-parton matrix element
! (MATFOR) with its dipole subtraction terms (SUBFOR), for the process
! selected by disorder's command-line flags (photon, Z, gamma/Z, W
! exchange; charged leptons and neutrinos).
!
! Four-parton configurations are built that approach single-unresolved
! limits: two outgoing partons collinear, one outgoing parton soft, one
! outgoing parton collinear to the incoming one. In each limit the
! matrix element must be reproduced by the sum of the dipoles, i.e.
! the ratio M(i) / sum S(i) -> 1 for every incoming parton i whose Born
! is non-zero. As DISENT uses fixed labels for some structures and
! relies on its permutation-symmetric phase space, both sides are
! symmetrised over the labellings of the outgoing partons 2, 3, 4, and
! dipoles whose three-parton Born contains the unresolved parton (as
! spectator or kept parton; in a calculation these are removed by the
! observable) are left out. This tests the parity-violating (Z, W)
! parts of MATFOR and CONTHR3 against those of MATTHR, including the
! spin correlations in the final-state gluon splittings.
!----------------------------------------------------------------------
program test_subtraction
  use types, only: dp
  use mod_parameters
  use test_utils
  implicit none
  ! DISENT's /COLFAC/ (renamed locally, cutoff clashes with mod_parameters)
  integer :: SCHEME_D, NF
  double precision :: CF, CA, TR, PI_D, PISQ, HF, CUTOFF_D, EQ(-6:6), SCALE_D
  common /COLFAC/ CF, CA, TR, PI_D, PISQ, HF, CUTOFF_D, EQ, SCALE_D, SCHEME_D, NF
  integer, parameter :: nkind = 7
  character(len=8), parameter :: kinds(nkind) = &
       & [character(len=8) :: 'fs34', 'fs23', 'fs24', 'soft4', 'soft2', 'is4', 'is2']
  integer, parameter :: unres(nkind) = [0, 0, 0, 4, 2, 4, 2]
  real(dp) :: P(4,7), Msym(-6:6), Ssym(-6:6), r, lam
  character(len=60) :: tag
  integer :: i, k, ifl

  call set_parameters()
  CF = 4.0_dp/3.0_dp; CA = 3; TR = 0.5_dp; NF = nflav
  PI_D = atan(1d0)*4; PISQ = PI_D**2; HF = 0.5_dp; CUTOFF_D = 1d-8
  SCHEME_D = 0; SCALE_D = 1
  EQ(0) = 0
  EQ(1) = -1D0/3
  EQ(2) = EQ(1) + 1
  do i = 1, 6
     if (i > 2) EQ(i) = EQ(i-2)
     EQ(-i) = -EQ(i)
  enddo

  lam = 1d-7
  do k = 1, nkind
     call limit_point(trim(kinds(k)), lam, P)
     call symmetrised(P, unres(k), Msym, Ssym)
     do ifl = -2, 2
        if (abs(Ssym(ifl)) < 1d-6 * maxval(abs(Ssym))) cycle   ! no Born for this parton
        r = Msym(ifl) / Ssym(ifl)
        write(tag,'(a,a,a,i2,a)') 'limit ', trim(kinds(k)), ', parton ', ifl, ': M/sum S -> 1'
        call check_close(trim(tag), r, 1.0_dp, 3d-3)
     enddo
  enddo
  call finish_tests()

contains

  ! M and the dipole sum, symmetrised over the labellings of 2,3,4
  subroutine symmetrised(P, iunres, Mt, St)
    real(dp), intent(in)  :: P(4,7)
    integer, intent(in)   :: iunres
    real(dp), intent(out) :: Mt(-6:6), St(-6:6)
    integer, parameter :: perm3(3,6) = reshape([2,3,4, 2,4,3, 3,2,4, 3,4,2, 4,2,3, 4,3,2], [3,6])
    integer, parameter :: kperm(6) = [4,3,1,1,2,1], lperm(6) = [1,1,4,3,1,2]   ! SUBFOR's K, L
    real(dp) :: PP(4,7), QM(4,7), M4(-6:6), S4(-6:6), JAC
    integer :: ip, j, upos
    Mt = 0; St = 0
    do ip = 1, 6
       PP = P
       PP(:,2) = P(:,perm3(1,ip)); PP(:,3) = P(:,perm3(2,ip)); PP(:,4) = P(:,perm3(3,ip))
       call MATFOR(PP, M4)
       Mt = Mt + M4
       upos = 0
       do j = 1, 3
          if (perm3(j,ip) == iunres) upos = j + 1
       enddo
       do j = 1, 6
          if (upos /= 0 .and. (kperm(j) == upos .or. lperm(j) == upos)) cycle
          call SUBFOR(j, 1d6, PP, QM, S4, JAC, *100)
          St = St + S4
100       continue
       enddo
    enddo
  end subroutine symmetrised

  ! A four-parton point near a limit, in DISENT's layout and Breit frame,
  ! built in the rest frame of W = p1 + q (p1 along +z there).
  subroutine limit_point(kind, lam, P)
    character(len=*), intent(in) :: kind
    real(dp), intent(in)  :: lam
    real(dp), intent(out) :: P(4,7)
    real(dp), parameter :: Q = 40, y = 0.4_dp, xi = 0.3_dp
    real(dp) :: E, E1, W, beta, f(4,2:4), Wv(4), rest(4), d(3), e4, pc, th
    integer :: a, b, c
    E = Q/2; E1 = E/xi; W = sqrt(4*E*E1 - 4*E**2); beta = (E1 - 2*E)/E1
    P = 0
    P(:,1) = [0.0_dp, 0.0_dp, E1, E1]
    P(:,5) = [0.0_dp, 0.0_dp, -2*E, 0.0_dp]
    P(:,6) = [E/y*2*sqrt(1-y), 0.0_dp, -E, E/y*(2-y)]
    P(:,7) = [E/y*2*sqrt(1-y), 0.0_dp,  E, E/y*(2-y)]
    Wv = [0.0_dp, 0.0_dp, 0.0_dp, W]
    select case (kind(1:2))
    case ('fs')
       read(kind(3:3),*) a; read(kind(4:4),*) b; c = 9 - a - b
       d = unit(0.37_dp, 1.9_dp)
       pc = (W*W - lam*W*W)/(2*W)
       f(:,c) = [pc*d, pc]
       rest = Wv - f(:,c)
       call twobody(rest, unit(0.61_dp, 0.4_dp), f(:,a), f(:,b))
    case ('so')
       read(kind(5:5),*) a
       d = unit(-0.25_dp, 2.7_dp); e4 = lam*W
       f(:,a) = [e4*d, e4]
       rest = Wv - f(:,a)
       call others(a, b, c)
       call twobody(rest, unit(0.52_dp, 1.1_dp), f(:,b), f(:,c))
    case ('is')
       read(kind(3:3),*) a
       th = sqrt(lam); e4 = 0.3_dp*W
       f(:,a) = [e4*sin(th)*cos(0.8_dp), e4*sin(th)*sin(0.8_dp), e4*cos(th), e4]
       rest = Wv - f(:,a)
       call others(a, b, c)
       call twobody(rest, unit(0.52_dp, 1.1_dp), f(:,b), f(:,c))
    end select
    do a = 2, 4
       P(:,a) = boostz(f(:,a), beta)
    enddo
  end subroutine limit_point

  subroutine others(a, b, c)
    integer, intent(in) :: a
    integer, intent(out) :: b, c
    select case (a)
    case (2); b = 3; c = 4
    case (3); b = 2; c = 4
    case default; b = 2; c = 3
    end select
  end subroutine others

  function unit(cth, phi) result(d)
    real(dp), intent(in) :: cth, phi
    real(dp) :: d(3)
    d = [sqrt(1-cth**2)*cos(phi), sqrt(1-cth**2)*sin(phi), cth]
  end function unit

  ! massless two-body decay of Ptot along +-dir in its rest frame
  subroutine twobody(Ptot, dir, pa, pb)
    real(dp), intent(in)  :: Ptot(4), dir(3)
    real(dp), intent(out) :: pa(4), pb(4)
    real(dp) :: Mm, h
    Mm = sqrt(Ptot(4)**2 - sum(Ptot(1:3)**2)); h = Mm/2
    pa = boost([h*dir, h], Ptot); pb = boost([-h*dir, h], Ptot)
  end subroutine twobody

  function boost(p, Pt) result(q)
    real(dp), intent(in) :: p(4), Pt(4)
    real(dp) :: q(4), Mm, bv(3), b2, g, bp
    Mm = sqrt(Pt(4)**2 - sum(Pt(1:3)**2)); bv = Pt(1:3)/Pt(4); b2 = sum(bv**2)
    if (b2 < 1d-300) then
       q = p; return
    endif
    g = Pt(4)/Mm; bp = sum(bv*p(1:3))
    q(1:3) = p(1:3) + ((g-1)*bp/b2 + g*p(4))*bv
    q(4) = g*(p(4) + bp)
  end function boost

  function boostz(p, beta) result(q)
    real(dp), intent(in) :: p(4), beta
    real(dp) :: q(4), g
    g = 1/sqrt(1 - beta**2)
    q = [p(1), p(2), g*(p(3) + beta*p(4)), g*(p(4) + beta*p(3))]
  end function boostz

end program test_subtraction
