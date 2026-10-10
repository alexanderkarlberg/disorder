! Size of the terms dropped with photon + Z exchange (8 Oct 2026): the
! interference of the boson on different quark lines (me31_keepint), i.e.
! for the symmetric observable here the axial pair-flavour sum, integrated
! over the full 3+1 phase space at fixed (x, Q^2): its contribution to
! dsigma/dx dQ^2 at O(alpha_s^2) [pb/GeV^2], for e- and e+ (their half
! difference is the parity-violating, xF3-like part). The integrand is
! finite (no singular limit of the boson-on-the-pair amplitude), so plain
! Monte Carlo with psmc's phase space; the outgoing labels symmetrised.
!   axsum31 x Q2 npoints [seed]
program axsum31
  use me31
  use ew31
  use psmc
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: s = 4*27.5_dp*920, gev2pb = 0.3893793721e9_dp
  integer, parameter :: perm(3,6) = reshape([2,3,4, 2,4,3, 3,2,4, 3,4,2, 4,2,3, 4,3,2], [3,6])
  real(dp), external :: alphasPDF
  real(dp) :: x, Q2, y, r(7), P(4,7), PP(4,7), eta, wps, xf(-6:6), f(-5:5), w, as, d(0:1), sm(0:1), sq(0:1)
  integer :: n, i, il, ip, fi, Q, seed, nseed
  integer, allocatable :: sd(:)
  logical :: ok
  character(32) :: arg
  call get_command_argument(1, arg); read(arg, *) x
  call get_command_argument(2, arg); read(arg, *) Q2
  call get_command_argument(3, arg); read(arg, *) n
  seed = 1
  if (command_argument_count() > 3) then
     call get_command_argument(4, arg); read(arg, *) seed
  endif
  call random_seed(size=nseed); allocate(sd(nseed))
  sd = [(1000003*seed + 7919*i, i = 1, nseed)]
  call random_seed(put=sd)
  call initPDFSetByName('NNPDF30_nlo_as_0118')
  call initPDF(0)
  y = Q2/(x*s)
  call psmc_init(x, Q2, y)
  as = alphasPDF(sqrt(Q2))
  ew31_mode = 1
  sm = 0; sq = 0
  do i = 1, n
     call random_number(r)
     call psmc_gen(3, r, P, eta, wps, ok)
     if (.not. ok) cycle
     call evolvePDF(eta, sqrt(Q2), xf); f = xf(-5:5)/eta
     w = y/x*wps/(2*eta*s)/(16*pi**2)*gev2pb*(as/(2*pi))**2/6
     do il = 0, 1
        ew31_lepton = il
        d(il) = 0
        do ip = 1, 6
           PP = P
           PP(:,2) = P(:,perm(1,ip)); PP(:,3) = P(:,perm(2,ip)); PP(:,4) = P(:,perm(3,ip))
           do fi = -5, 5
              if (fi == 0) cycle
              do Q = 1, 5
                 if (Q == abs(fi)) then
                    d(il) = d(il) + 0.5_dp*f(fi)*dropped([fi, fi, fi, -fi])
                 else
                    d(il) = d(il) + f(fi)*dropped([fi, fi, sign(Q, fi), -sign(Q, fi)])
                 endif
              enddo
           enddo
        enddo
        sm(il) = sm(il) + w*d(il); sq(il) = sq(il) + (w*d(il))**2
     enddo
  enddo
  sm = sm/n; sq = sqrt(max(sq/n - sm**2, 0.0_dp)/(n - 1))
  write(*,'(a,f8.4,a,f9.1,a,f6.3)') ' x', x, '  Q2', Q2, '  y', y
  write(*,'(a,es12.4,a,es10.2)') ' dropped, e-: dsigma/dx dQ2 [pb/GeV2] =', sm(0), ' +-', sq(0)
  write(*,'(a,es12.4,a,es10.2)') ' dropped, e+: dsigma/dx dQ2 [pb/GeV2] =', sm(1), ' +-', sq(1)
contains
  real(dp) function dropped(fl)
    integer, intent(in) :: fl(4)
    real(dp) :: a, b
    me31_keepint = .true.;  call me31_tree(PP, fl, a)
    me31_keepint = .false.; call me31_tree(PP, fl, b)
    dropped = a - b
  end function dropped
end program axsum31
