!----------------------------------------------------------------------
! NNLO DIS 1+1 (O(alpha_s^2) coefficient of 1+1 observables) three ways,
! in one DISENT run (photon exchange, fixed x and Q^2, mu_R = mu_F = Q):
!
!  1. DISENT + P2B: the O(alpha_s^2) 2+1 weights of DISENT (VIRTHR, COLFOR,
!     MATFOR, SUBFOR), each with its Born projection subtracted,
!     E2 = sum w [O(event) - O(Born 1+1)]; the inclusive part (C2_incl O(Born))
!     is added in the analysis from hoppet's NNLO structure functions.
!  2. tau_2 slicing + P2B: the same, with DISENT's O(alpha_s^2) 2+1 weights
!     replaced by tau_2 slicing (invariant measure): the four-parton reals
!     with tau_2 > cut, and below the cut the one-loop SCET cumulant on the
!     three-parton Born (lp21_weights). Absolute cuts tau_2 > c and relative
!     cuts tau_2 > rho tau_1 (tau_1 of the Born, resp. of the four-parton event).
!  3. pure tau_1 slicing: the O(alpha_s^2) weights of DISENT with
!     tau_1 > tau_cut (above the cut); the NNLO leading-power cumulant below
!     the cut times the Born is added in the analysis (dis_tau1_lp).
!     tau_1 = min over partitions [2 x P.p(beam region) + m^2(jet region)]/Q^2
!     (recoil-free, thrust-like at leading power in the Breit frame).
!  The O(alpha_s) parts (P2B: E1; tau_1 slicing: A1) are kept as checks.
!
! Needs DISENT with its two-parton Born given to USER (USER(2,0,0)): the
! slicing study builds a copy of the fixed libdisent.f with that call.
!
! Observables (lab frame, anti-k_t R = 1, E scheme, jets with p_t > 5 GeV
! and -1 < eta < 2.5; proton along +z): 1 total, 2 at least one jet,
! 3-9 leading-jet p_t, 10-15 leading-jet eta, 16 at least two jets.
!
! Usage: nnlo11 -x X -Q2 Q2 -s S -Ep EP -nev N -seed1 I -seed2 J [-cutoff 1e-8]
!----------------------------------------------------------------------
module nnlo11_run
  use types, only: dp
  use mod_slicing_scet
  use tau2_run, only: xfix, Q2fix, measure, sdis, tau2_of, lp21_weights, nb
  implicit none
  integer, parameter :: nb11 = 16, nt1 = 9, nt2 = 11
  real(dp), parameter :: tc1(nt1) = [1e-1_dp, 3e-2_dp, 1e-2_dp, 3e-3_dp, 1e-3_dp, 3e-4_dp, 1e-4_dp, 3e-5_dp, 1e-5_dp]
  ! tau_2 cuts: absolute values, then rho for tau_2 > rho tau_1
  real(dp), parameter :: tc2(nt2) = [1e-2_dp, 1e-3_dp, 1e-4_dp, 1e-5_dp, &
       & 1e-1_dp, 3e-2_dp, 1e-2_dp, 3e-3_dp, 1e-3_dp, 3e-4_dp, 1e-4_dp]
  logical, parameter :: rel2(nt2) = [.false., .false., .false., .false., .true., .true., .true., .true., .true., .true., .true.]
  real(dp), parameter :: ptedge(8) = [5.0_dp, 10.0_dp, 14.0_dp, 15.0_dp, 16.0_dp, 17.0_dp, 20.0_dp, 40.0_dp]
  real(dp), parameter :: etaedge(7) = [-1.0_dp, -0.6_dp, -0.4_dp, -0.2_dp, 0.2_dp, 1.0_dp, 2.5_dp]
  real(dp), parameter :: ptjmin = 5.0_dp, etamin = -1.0_dp, etamax = 2.5_dp, Rjet = 1.0_dp
  real(dp), save :: Ep = 920.0_dp, Ylab = 0, t1min = 1e-7_dp
  ! diagnostics: number and CPU time of the below-cut (LP) evaluations
  integer(8), save :: nlp = 0
  real(dp), save :: tlp = 0
  ! per-event values
  real(dp), save :: eB0(nb11) = 0, eE1(nb11) = 0, eE2d(nb11) = 0, eE2s(nt2,nb11) = 0
  real(dp), save :: eA1(nt1,nb11) = 0, eA2(nt1,nb11) = 0, fbB(nb11) = 0
  ! sums and sums of squares; D2s = E2s - E2d, Dp = A2 - E2d (per event)
  real(dp), save :: sB0(2,nb11) = 0, sE1(2,nb11) = 0, sE2d(2,nb11) = 0, sE2s(2,nt2,nb11) = 0
  real(dp), save :: sA1(2,nt1,nb11) = 0, sA2(2,nt1,nb11) = 0, sD2s(2,nt2,nb11) = 0, sDp(2,nt1,nb11) = 0
  integer(8), save :: nevt = 0, nborn11 = 0, nno = 0
  logical, save :: haveborn = .false.
  real(dp), save :: as2pi_save = 0, fbB_save(nb11) = 0
  real(dp), external :: DOT, alphasPDF

contains

  ! tau_1 = min over partitions [2 x P.p(S_B) + (sum S_J)^2] / Q^2, S_J not empty
  real(dp) function tau1_of(n, p) result(t)
    integer, intent(in) :: n
    real(dp), intent(in) :: p(4,7)
    real(dp) :: eta, Q2, pp(2:4), v(4), m2, b
    integer :: mask, j
    Q2 = -DOT(p,5,5); eta = 2 * DOT(p,1,6) / sdis
    do j = 2, n
       pp(j) = 2 * xfix * DOT(p,1,j) / eta     ! 2 x P.p_j
    enddo
    t = huge(1.0_dp)
    do mask = 1, 2**(n-1) - 1                 ! bit j-2 set: parton j in the jet
       v = 0; b = 0
       do j = 2, n
          if (btest(mask, j-2)) then
             v = v + p(:,j)
          else
             b = b + pp(j)
          endif
       enddo
       m2 = v(4)**2 - v(1)**2 - v(2)**2 - v(3)**2
       t = min(t, b + max(m2, 0.0_dp))
    enddo
    t = t / Q2
  end function tau1_of

  ! lab-frame jet observables of the final-state partons 2..n
  subroutine obs11(n, p, fb)
    integer, intent(in) :: n
    real(dp), intent(in) :: p(4,7)
    real(dp), intent(out) :: fb(nb11)
    real(dp) :: eta, pprot(4), k(4), tot(4), bv(3), b2, g, z(3), e1(3), e2(3), q(4)
    real(dp) :: jv(4,3), pt(3), y(3), ph(3), dmin, d, ptl, etal
    integer :: i, j, nj, ia, ja, nq, nj5, il
    logical :: alive(3), isjet(3)
    eta = 2 * DOT(p,1,6) / sdis
    pprot = p(:,1) / eta; k = p(:,6)
    tot = pprot + k
    bv = tot(1:3) / tot(4); b2 = sum(bv * bv); g = 1 / sqrt(1 - b2)
    call boostv(pprot, q); z = q(1:3) / sqrt(sum(q(1:3)**2))
    e1 = [1.0_dp, 0.0_dp, 0.0_dp] - z(1) * z
    if (sum(e1 * e1) < 1e-6_dp) e1 = [0.0_dp, 1.0_dp, 0.0_dp] - z(2) * z
    e1 = e1 / sqrt(sum(e1 * e1))
    e2 = [z(2) * e1(3) - z(3) * e1(2), z(3) * e1(1) - z(1) * e1(3), z(1) * e1(2) - z(2) * e1(1)]
    ! particles: final partons 2..n, without the zero momenta of DISENT's
    ! counter-configurations (three partons in a four-parton array) and without
    ! partons exactly along the beam (collinear configurations: never in a jet)
    nq = 0
    do i = 2, n
       call boostv(p(:,i), q)
       if (q(4) <= 0) cycle
       q = [dot_product(q(1:3), e1), dot_product(q(1:3), e2), dot_product(q(1:3), z), q(4)]
       if (q(1)**2 + q(2)**2 <= (1e-12_dp * q(4))**2) cycle
       nq = nq + 1
       jv(:,nq) = q
    enddo
    ! anti-k_t clustering (E scheme) of up to three particles: the smallest of
    ! d_iB = 1/pt_i^2 and d_ij = min(1/pt_i^2, 1/pt_j^2) dR_ij^2/R^2; a beam
    ! distance makes i a final jet, a pair distance merges i and j
    alive = .false.; alive(1:nq) = .true.; isjet = .false.
    do
       nj = count(alive)
       if (nj == 0) exit
       call kin()
       dmin = huge(1.0_dp); ia = 0; ja = 0
       do i = 1, nq
          if (.not. alive(i)) cycle
          d = 1 / max(pt(i), 1e-300_dp)**2
          if (d < dmin) then
             dmin = d; ia = i; ja = 0
          endif
          do j = i + 1, nq
             if (.not. alive(j)) cycle
             d = min(1 / max(pt(i), 1e-300_dp)**2, 1 / max(pt(j), 1e-300_dp)**2) * &
                  & ((y(i) - y(j))**2 + dphi(ph(i), ph(j))**2) / Rjet**2
             if (d < dmin) then
                dmin = d; ia = i; ja = j
             endif
          enddo
       enddo
       if (ja == 0) then
          alive(ia) = .false.; isjet(ia) = .true.
       else
          jv(:,ia) = jv(:,ia) + jv(:,ja); alive(ja) = .false.
       endif
    enddo
    alive = isjet
    call kin()
    nj5 = 0; ptl = -1; etal = 0; il = 0
    do i = 1, nq
       if (.not. isjet(i)) cycle
       if (pt(i) > ptjmin .and. y(i) > etamin .and. y(i) < etamax) then
          nj5 = nj5 + 1
          if (pt(i) > ptl) then
             ptl = pt(i); etal = y(i); il = i
          endif
       endif
    enddo
    fb = 0
    fb(1) = 1
    if (nj5 >= 1) then
       fb(2) = 1
       do i = 1, 7
          if (ptl >= ptedge(i) .and. ptl < ptedge(i+1)) fb(2+i) = 1
       enddo
       do i = 1, 6
          if (etal >= etaedge(i) .and. etal < etaedge(i+1)) fb(9+i) = 1
       enddo
    endif
    if (nj5 >= 2) fb(16) = 1
  contains
    subroutine boostv(a, b)
      real(dp), intent(in) :: a(4)
      real(dp), intent(out) :: b(4)
      real(dp) :: bp
      bp = sum(bv * a(1:3))
      b(1:3) = a(1:3) + ((g - 1) * bp / b2 - g * a(4)) * bv
      b(4) = g * (a(4) - bp)
    end subroutine boostv
    subroutine kin()
      integer :: m
      do m = 1, nq
         if (.not. alive(m)) cycle
         pt(m) = sqrt(jv(1,m)**2 + jv(2,m)**2)
         y(m) = 0.5_dp * log((jv(4,m) + jv(3,m)) / max(jv(4,m) - jv(3,m), 1e-300_dp)) + Ylab
         ph(m) = atan2(jv(2,m), jv(1,m))
      enddo
    end subroutine kin
    real(dp) function dphi(a1, a2)
      real(dp), intent(in) :: a1, a2
      dphi = abs(a1 - a2)
      if (dphi > acos(-1.0_dp)) dphi = 2 * acos(-1.0_dp) - dphi
    end function dphi
  end subroutine obs11

  subroutine nnlo11_user(n, na, ntyp, p, s, weight, scale2)
    integer, intent(in) :: n, na, ntyp
    real(dp), intent(in) :: s, p(4,7), weight(-6:6), scale2
    real(dp) :: eta, Q, xf(-6:6), as2pi, w, fb(nb11), t1, t1r, T2, tcs(nt2), wl(nt2), d(nb11), tq0, tq1
    integer :: it
    if (n == 0) then
       call endevent()
       return
    endif
    sdis = s
    eta = 2 * DOT(p,1,6) / s
    Q = sqrt(-DOT(p,5,5))
    call EvolvePDF(eta, Q, xf)
    call mask_pdf(xf)
    as2pi = alphasPDF(Q) / (2 * pi); as2pi_save = as2pi
    w = dot_product(weight, xf) * as2pi**na
    call obs11(n, p, fb)
    if (n == 2) then
       nborn11 = nborn11 + 1
       fbB = fb; fbB_save = fb; haveborn = .true.
       eB0 = eB0 + w * fb
       return
    endif
    if (.not. haveborn) then
       nno = nno + 1          ! (does not happen: DISENT gives the Born first)
       return
    endif
    d = fb - fbB
    if (na == 1 .and. n == 3) then
       eE1 = eE1 + w * d
       t1 = tau1_of(3, p)
       do it = 1, nt1
          if (t1 > tc1(it)) eA1(it,:) = eA1(it,:) + w * fb
       enddo
       ! the below-cut weights only enter times d = O(event) - O(Born): skip
       ! them (and the costly soft function) when the Born is in the same bins,
       ! and for Borns at the 1+1 edge (tau_1 < t1min), where d -> 0 and the
       ! soft function's angular integrals become very slow
       if (any(d /= 0) .and. t1 > t1min) then
          tcs = merge(tc2 * t1, tc2, rel2)
          call cpu_time(tq0)
          call lp21_weights(p, s, weight, xf, eta, Q, as2pi, nt2, tcs, wl)
          call cpu_time(tq1); tlp = tlp + (tq1 - tq0); nlp = nlp + 1
          do it = 1, nt2
             eE2s(it,:) = eE2s(it,:) + wl(it) * d
          enddo
       endif
    elseif (na == 2) then
       eE2d = eE2d + w * d
       t1 = tau1_of(n, p)
       do it = 1, nt1
          if (t1 > tc1(it)) eA2(it,:) = eA2(it,:) + w * fb
       enddo
       if (n == 4 .and. ntyp == 0) then
          T2 = tau2_of(4, p) / Q
          t1r = t1
          do it = 1, nt2
             if (T2 > merge(tc2(it) * t1r, tc2(it), rel2(it))) eE2s(it,:) = eE2s(it,:) + w * d
          enddo
       endif
    endif
  end subroutine nnlo11_user

  subroutine endevent()
    integer :: it
    nevt = nevt + 1
    call acc1(sB0, eB0); call acc1(sE1, eE1); call acc1(sE2d, eE2d)
    do it = 1, nt2
       call acc1(sE2s(:,it,:), eE2s(it,:)); call acc1(sD2s(:,it,:), eE2s(it,:) - eE2d)
    enddo
    do it = 1, nt1
       call acc1(sA1(:,it,:), eA1(it,:)); call acc1(sA2(:,it,:), eA2(it,:))
       call acc1(sDp(:,it,:), eA2(it,:) - eE2d)
    enddo
    eB0 = 0; eE1 = 0; eE2d = 0; eE2s = 0; eA1 = 0; eA2 = 0; haveborn = .false.
  contains
    subroutine acc1(sm, v)
      real(dp), intent(inout) :: sm(2,nb11)
      real(dp), intent(in) :: v(nb11)
      sm(1,:) = sm(1,:) + v; sm(2,:) = sm(2,:) + v * v
    end subroutine acc1
  end subroutine endevent

  subroutine report11()
    integer :: u, ib
    open(newunit=u, file='nnlo11.dat', status='replace')
    write(u,*) nevt, nb11, nt1, nt2
    write(u,*) tc1
    write(u,*) tc2
    write(u,*) merge(1, 0, rel2)
    write(u,*) as2pi_save, fbB_save
    write(u,*) sB0, sE1, sE2d
    write(u,*) sE2s, sD2s
    write(u,*) sA1, sA2, sDp
    close(u)
    write(*,'(a,i12,a,i12,a,i6,a,es14.6)') ' events ', nevt, '  1+1 Borns ', nborn11, '  without Born ', nno, &
         & '  as/2pi ', as2pi_save
    write(*,'(a,16f4.0)') ' Born 1+1 bins ', fbB_save
    write(*,'(a,i12,a,f10.2,a)') ' below-cut evaluations ', nlp, '  CPU ', tlp, ' s'
    write(*,'(a)') ' bin   Born           E1 (P2B)       E2 DISENT+P2B   E2 slicing+P2B (tau2 1e-4 abs)   A2 (tau1 > 1e-3)'
    do ib = 1, nb11
       write(*,'(i4,5es16.6)') ib, sB0(1,ib), sE1(1,ib), sE2d(1,ib), sE2s(1,3,ib), sA2(1,5,ib)
    enddo
  end subroutine report11

end module nnlo11_run

program nnlo11
  use types, only: dp
  use mod_parameters, only: nflav, NC, CC, noZ, Zonly, intonly, neutrino, positron
  use tau2_run, only: xfix, Q2fix, measure, slice_cuts, slice_muf
  use nnlo11_run
  use mod_slicing_scet, only: soft_tol, pdf_mask, scet_set_colour, CF, CA, TF
  use sub_defs_io
  implicit none
  character(len=100) :: pdf
  real(dp) :: s, npow1, npow2, cutoff
  integer :: nev, seed1, seed2
  external :: DISENTFULL

  pdf = string_val_opt('-pdf', 'NNPDF30_nlo_as_0118')
  xfix = dble_val_opt('-x', 0.01_dp)
  Q2fix = dble_val_opt('-Q2', 400.0_dp)
  s = dble_val_opt('-s', 4 * 27.5_dp * 920.0_dp)
  Ep = dble_val_opt('-Ep', 920.0_dp)
  t1min = dble_val_opt('-t1min', 1e-7_dp)
  nev = int_val_opt('-nev', 100000)
  seed1 = int_val_opt('-seed1', 12345)
  seed2 = int_val_opt('-seed2', 67890)
  npow1 = dble_val_opt('-npow1', 2.0_dp)
  npow2 = dble_val_opt('-npow2', 4.0_dp)
  soft_tol = dble_val_opt('-softtol', 1e-9_dp)
  pdf_mask = int_val_opt('-pdfmask', 0)
  cutoff = dble_val_opt('-cutoff', 1e-8_dp)
  measure = 1                                    ! invariant tau_2 measure for the slicing
  Ylab = 0.5_dp * log(Ep / (s / (4 * Ep)))       ! lab rapidity of the (P + k) rest frame
  call scet_set_colour(4.0_dp/3.0_dp, 3.0_dp, 0.5_dp)

  nflav = 5; NC = .true.; CC = .false.; noZ = .true.; Zonly = .false.
  intonly = .false.; neutrino = .false.; positron = .false.
  call InitPDFsetByName(trim(pdf))
  call InitPDF(0)
  write(*,'(a,a,a,f8.5,a,f9.2,a,f10.1,a,f8.2,a,i10,a,i2,a,es9.2)') ' pdf ', trim(pdf), '  x ', xfix, '  Q2 ', Q2fix, &
       & '  s ', s, '  Ep ', Ep, '  nev ', nev, '  pdfmask ', pdf_mask, '  cutoff ', cutoff

  call DISENTFULL(nev, s, 5, nnlo11_user, slice_cuts, seed1, seed2, npow1, npow2, &
       & cutoff, 2, slice_muf, CF, CA, TF, .false.)
  call report11()
end program nnlo11
