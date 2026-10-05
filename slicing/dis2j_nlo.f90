!----------------------------------------------------------------------
! DISENT reference for ZEUS-like dijets at NLO (photon exchange, mu_R =
! mu_F = Q): the selection of nnlo21/validation/nnlojet_epLJJ_zeus2j.run
! (as dis31/nlo31 mode 2): 125 < Q^2 < 20000 GeV^2, 0.2 < y < 0.6; inclusive
! kt jets (R = 1, E-scheme) in the Breit frame, E_T = sqrt(p_T^2 + m^2) > 8
! GeV, then jets outside -1 < eta_lab < 2.5 dropped, >= 2 jets, m12 > 20 GeV
! (two leading jets in Breit p_T). Lab frame from the proton (direction of
! the incoming parton) and the incoming lepton: t = P/(2Ep) + k/(2Ee),
! z = P/(2Ep) - k/(2Ee), proton along +z.
!
! Output: sigma per bin in pb at O(alpha_s) (LO) and the O(alpha_s^2)
! coefficient (NLO correction), bins as nlo31/sliced21 mode 2 (total, Q^2,
! ptavg_12, m12). Errors from the event-by-event scatter.
!
! Usage: dis2j_nlo -pdf NAME -nev N -seed1 I -seed2 J [-cutoff 1e-8]
!----------------------------------------------------------------------
module dis2j_run
  use types, only: dp
  implicit none
  integer, parameter :: nob = 15
  real(dp), parameter :: Ee = 27.5_dp, Ep = 920.0_dp
  real(dp), parameter :: zetmin = 8, zetalo = -1, zetahi = 2.5_dp, zm12min = 20, rjet = 1
  real(dp), parameter :: zq2e(0:6) = [125.0_dp, 250.0_dp, 500.0_dp, 1000.0_dp, 2000.0_dp, 5000.0_dp, 20000.0_dp]
  real(dp), parameter :: zpte(0:4) = [8.0_dp, 15.0_dp, 22.0_dp, 30.0_dp, 60.0_dp]
  real(dp), parameter :: zmje(0:4) = [20.0_dp, 30.0_dp, 45.0_dp, 65.0_dp, 120.0_dp]
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
  real(dp), save :: ev(nob,2) = 0, s1(nob,2) = 0, s2(nob,2) = 0
  integer(8), save :: nevt = 0
  logical, save :: checked = .false.
  real(dp), external :: DOT, alphasPDF

contains

  subroutine z_cuts(s, xmn, xmx, q2mn, q2mx, ymn, ymx)
    real(dp), intent(in) :: s
    real(dp), intent(out) :: xmn, xmx, q2mn, q2mx, ymn, ymx
    xmn = 0; xmx = 0; q2mn = zq2e(0); q2mx = zq2e(6); ymn = 0.2_dp; ymx = 0.6_dp
  end subroutine z_cuts

  subroutine z_muf(p, s, scale)
    real(dp), intent(in) :: p(4,7), s
    real(dp), intent(out) :: scale
    scale = 1   ! (mu_F/Q)^2
  end subroutine z_muf

  subroutine z_user(n, na, ntyp, p, s, weight, scale2)
    integer, intent(in) :: n, na, ntyp
    real(dp), intent(in) :: s, p(4,7), weight(-6:6), scale2
    real(dp) :: eta, Q2, Q, xf(-6:6), w
    logical :: inb(nob)
    if (n == 0) then
       nevt = nevt + 1
       s1 = s1 + ev; s2 = s2 + ev**2; ev = 0
       return
    endif
    Q2 = -DOT(p,5,5); Q = sqrt(Q2)
    if (.not. checked) then
       ! Breit frame with the incoming parton along +z and q = (0, 0, -Q, 0)
       if (abs(p(3,5) + Q) > 1e-6_dp*Q .or. abs(p(4,5)) > 1e-6_dp*Q .or. p(3,1) <= 0 &
            & .or. abs(p(1,1)) + abs(p(2,1)) > 1e-9_dp*p(4,1)) stop 'dis2j_nlo: not the expected Breit frame'
       checked = .true.
    endif
    call bins(n, p, Q2, s, inb)
    if (.not. inb(1)) return
    eta = 2 * DOT(p,1,6) / s
    call EvolvePDF(eta, Q, xf)
    w = dot_product(weight, xf) * (alphasPDF(Q) / (2 * pi))**na
    where (inb) ev(:,na) = ev(:,na) + w
  end subroutine z_user

  subroutine bins(n, p, Q2, s, inb)
    integer, intent(in) :: n
    real(dp), intent(in) :: p(4,7), Q2, s
    logical, intent(out) :: inb(nob)
    real(dp) :: jets(4,4), Pp(4), tl(4), zl(4), El, pzl, pabs, etal, pt2j(4), pj(4), m12, ptavg, k(4)
    integer :: nj, i, i1, i2, m, b
    inb = .false.
    call kt_jets(p(:,2:n), n - 1, jets, nj)
    k = p(:,6)
    Pp = p(:,1) * (s / 2) / mdot(p(:,1), k)
    tl = Pp / (2 * Ep) + k / (2 * Ee); zl = Pp / (2 * Ep) - k / (2 * Ee)
    m = 0
    do i = 1, nj
       pj = jets(:,i)
       if (pj(1)**2 + pj(2)**2 + max(mdot(pj, pj), 0.0_dp) <= zetmin**2) cycle
       El = mdot(pj, tl); pzl = -mdot(pj, zl)
       pabs = sqrt(max(El**2 - max(mdot(pj, pj), 0.0_dp), 0.0_dp))
       if (pabs - abs(pzl) <= 0) cycle
       etal = 0.5_dp * log((pabs + pzl) / (pabs - pzl))
       if (etal <= zetalo .or. etal >= zetahi) cycle
       m = m + 1; jets(:,m) = pj; pt2j(m) = pj(1)**2 + pj(2)**2
    enddo
    if (m < 2) return
    i1 = maxloc(pt2j(1:m), 1); pt2j(i1) = -1; i2 = maxloc(pt2j(1:m), 1)
    pj = jets(:,i1) + jets(:,i2)
    m12 = sqrt(max(mdot(pj, pj), 0.0_dp))
    if (m12 <= zm12min) return
    ptavg = (sqrt(jets(1,i1)**2 + jets(2,i1)**2) + sqrt(jets(1,i2)**2 + jets(2,i2)**2)) / 2
    inb(1) = .true.
    do b = 1, 6
       inb(1 + b) = Q2 >= zq2e(b-1) .and. Q2 < zq2e(b)
    enddo
    do b = 1, 4
       inb(7 + b) = ptavg >= zpte(b-1) .and. ptavg < zpte(b)
       inb(11 + b) = m12 >= zmje(b-1) .and. m12 < zmje(b)
    enddo
  end subroutine bins

  ! inclusive kt (R = rjet, E-scheme) in the Breit frame: all jets
  subroutine kt_jets(p, n, jets, nj)
    integer, intent(in) :: n
    real(dp), intent(in) :: p(4,n)
    real(dp), intent(out) :: jets(4,4)
    integer, intent(out) :: nj
    real(dp) :: q(4,n), dmin, d, pt2(n), yr(n), ph(n), dy, dphi
    integer :: m, i, j, ii, jj, act(n)
    logical :: beam
    q = p; m = n
    act = [(i, i = 1, n)]
    nj = 0
    do while (m > 0)
       do i = 1, m
          pt2(i) = q(1,act(i))**2 + q(2,act(i))**2
          yr(i) = 0.5_dp * log(max(q(4,act(i)) + q(3,act(i)), 1d-300) / max(q(4,act(i)) - q(3,act(i)), 1d-300))
          ph(i) = atan2(q(2,act(i)), q(1,act(i)))
       enddo
       dmin = huge(1.0_dp); ii = 0; jj = 0; beam = .true.
       do i = 1, m
          if (pt2(i) < dmin) then
             dmin = pt2(i); ii = i; beam = .true.
          endif
          do j = i + 1, m
             dy = yr(i) - yr(j)
             dphi = abs(ph(i) - ph(j)); if (dphi > pi) dphi = 2 * pi - dphi
             d = min(pt2(i), pt2(j)) * (dy**2 + dphi**2) / rjet**2
             if (d < dmin) then
                dmin = d; ii = i; jj = j; beam = .false.
             endif
          enddo
       enddo
       if (beam) then
          nj = nj + 1; jets(:,nj) = q(:,act(ii))
          act(ii) = act(m); m = m - 1
       else
          q(:,act(ii)) = q(:,act(ii)) + q(:,act(jj))
          act(jj) = act(m); m = m - 1
       endif
    enddo
  end subroutine kt_jets

  pure real(dp) function mdot(a, b)
    real(dp), intent(in) :: a(4), b(4)
    mdot = a(4) * b(4) - a(1) * b(1) - a(2) * b(2) - a(3) * b(3)
  end function mdot

  subroutine z_report()
    real(dp) :: m(nob,2), e(nob,2)
    integer :: i
    ! DISENT's weights include 1/NEV and the conversion to pb (NRM in
    ! DISENTFULL): the cross section is the sum over events
    m = s1
    e = sqrt(max(s2 - s1**2 / nevt, 0.0_dp))
    write(*,'(a,i12)') ' events ', nevt
    write(*,'(a)') ' DCELL bin  lo  hi   LO [pb]  +-   NLO coefficient [pb]  +-'
    write(*,'(a,4es16.8)') ' DCELL total     0     0', m(1,1), e(1,1), m(1,2), e(1,2)
    do i = 1, 6
       write(*,'(a,2f8.0,4es16.8)') ' DCELL q2   ', zq2e(i-1), zq2e(i), m(1+i,1), e(1+i,1), m(1+i,2), e(1+i,2)
    enddo
    do i = 1, 4
       write(*,'(a,2f8.0,4es16.8)') ' DCELL ptavg', zpte(i-1), zpte(i), m(7+i,1), e(7+i,1), m(7+i,2), e(7+i,2)
    enddo
    do i = 1, 4
       write(*,'(a,2f8.0,4es16.8)') ' DCELL m12  ', zmje(i-1), zmje(i), m(11+i,1), e(11+i,1), m(11+i,2), e(11+i,2)
    enddo
  end subroutine z_report
end module dis2j_run

program dis2j_nlo
  use types, only: dp
  use mod_parameters, only: nflav, NC, CC, noZ, Zonly, intonly, neutrino, positron
  use dis2j_run
  use sub_defs_io
  implicit none
  character(len=100) :: pdf
  real(dp) :: s, cutoff
  integer :: nev, seed1, seed2
  external :: DISENTFULL
  pdf = string_val_opt('-pdf', 'NNPDF30_nlo_as_0118')
  nev = int_val_opt('-nev', 100000)
  seed1 = int_val_opt('-seed1', 12345)
  seed2 = int_val_opt('-seed2', 67890)
  cutoff = dble_val_opt('-cutoff', 1e-8_dp)
  s = 4 * Ee * Ep
  nflav = 5; NC = .true.; CC = .false.; noZ = .true.; Zonly = .false.
  intonly = .false.; neutrino = .false.; positron = .false.
  call InitPDFsetByName(trim(pdf))
  call InitPDF(0)
  write(*,'(a,a,a,i12,a,2i10)') ' dis2j_nlo pdf ', trim(pdf), '  nev', nev, '  seeds', seed1, seed2
  call DISENTFULL(nev, s, 5, z_user, z_cuts, seed1, seed2, 2.0_dp, 4.0_dp, &
       & cutoff, 2, z_muf, 4.0_dp/3, 3.0_dp, 0.5_dp, .false.)
  call z_report()
end program dis2j_nlo
