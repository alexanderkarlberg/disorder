!-----------------------------------------------------------------------
! Tests of hard21 (two-loop DIS 2+1 hard function) at random Born points
! in the Breit frame:
!  1. tree: h(0) against DISENT's MATTHR (ratio constant over points and
!     channels up to the charge and colour factors);
!  2. one loop: H^(1)/H^(0) (SCET, mu = Q) against the hard function of the
!     NLO slicing (slicing/mod_tau2_run.f90): hard_fact (from DISENT's
!     VIRTHR and the CS I operator) plus DISENT's non-factorising one-loop
!     part (LEIV, ERTV) over the Born;
!  3. two loop: continuity across the region boundaries, and the
!     renormalisation-group structure is checked in harness_rg21 (todo).
!-----------------------------------------------------------------------
program harness_hard21
  use hard21
  use mod_slicing_scet, only: hard_fact, scet_set_colour
  use mod_parameters, only: nflav, NC, CC, noZ, Zonly, intonly, neutrino, positron
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
  integer :: SCHEME, NF
  real(dp) :: CF, CA, TR, PIc, PISQ, HF, CUTOFF, EQ(-6:6), SCALE
  common /COLFAC/ CF, CA, TR, PIc, PISQ, HF, CUTOFF, EQ, SCALE, SCHEME, NF
  real(dp) :: P(4,7), M(-6:6), h(0:2), hc(0:2), Q, Q2, y, xp, th, ph, r1, r2, rat, ref, qqnf
  real(dp) :: l12, l13, l23, eq1, worst1(2), ratmin(2), ratmax(2)
  real(dp), external :: DOT, LEIV, ERTV
  integer :: ip, ichan, n
  real(dp) :: rnd
  CF = 4.0_dp/3; CA = 3; TR = 0.5_dp; PIc = pi; PISQ = pi**2; HF = 0.5_dp; NF = 5; SCALE = 1
  EQ = 0
  do n = 1, 5
     EQ(n) = merge(2.0_dp/3, -1.0_dp/3, mod(n, 2) == 0); EQ(-n) = -EQ(n)
  enddo
  call scet_set_colour(CF, CA, TR)
  nflav = 5; NC = .true.; CC = .false.; noZ = .true.; Zonly = .false.
  intonly = .false.; neutrino = .false.; positron = .false.
  worst1 = 0; ratmin = huge(1.0_dp); ratmax = -huge(1.0_dp)
  call random_seed(put=[(12345 + 7*n, n = 1, 33)])
  do ip = 1, 40
     call random_number(r1); Q = 10 + 90*r1; Q2 = Q*Q
     call random_number(r1); y = 0.1_dp + 0.8_dp*r1
     call random_number(r1); xp = 0.05_dp + 0.9_dp*r1      ! x/xi of the Born
     call random_number(r1); th = acos(-1 + 2*r1)
     call random_number(r1); ph = 2*pi*r1
     call breit_born(Q, y, xp, th, ph, P)
     call MATTHR(P, M)
     l12 = log(2*DOT(P,1,2)/Q2); l13 = log(2*DOT(P,1,3)/Q2); l23 = log(2*DOT(P,2,3)/Q2)
     do ichan = 1, 2
        if (ichan == 1) then
           eq1 = -1.0_dp/3
           call hard21_eval(P, 1, eq1, h)
           rat = M(1)/eq1**2 / h(0) * Q2**2
           qqnf = -((4*pi/137)**2*4/Q2)*(2*LEIV(P, P(1,6), 2, -1, 3) - Q2/2*ERTV(P, 2, -1, 3))
           ref = hard_fact(.false., l12, l13, l23) + qqnf/(M(1)/eq1**2)
        else
           eq1 = 2.0_dp/3
           call hard21_eval(P, 2, eq1, h)
           rat = M(0)/(2*(1.0_dp/9 + 4.0_dp/9) + 1.0_dp/9) / h(0) * Q2**2
           qqnf = 0.5_dp/(4.0_dp/3)*((4*pi/137)**2*4/Q2)*(2*LEIV(P, P(1,6), 2, 3, -1) - Q2/2*ERTV(P, 2, 3, -1))
           ref = hard_fact(.true., l12, l13, l23) + qqnf/(M(0)/(2*(1.0_dp/9 + 4.0_dp/9) + 1.0_dp/9))
        endif
        ratmin(ichan) = min(ratmin(ichan), rat); ratmax(ichan) = max(ratmax(ichan), rat)
        worst1(ichan) = max(worst1(ichan), abs(h(1) - ref)/max(1.0_dp, abs(ref)))
        if (ip <= 3) print '(a,i2,a,i2,a,es14.6,a,2es16.8,a,es14.6)', 'pt', ip, ' chan', ichan, &
             '  tree ratio', rat, '  H1/H0 hard21, slicing', h(1), ref, '  H2/H0', h(2)
     enddo
  enddo
  print '(a,2(2es14.6,2x))', 'tree ratio min/max (quark, gluon):', ratmin(1), ratmax(1), ratmin(2), ratmax(2)
  print '(a,2es10.2)', 'one loop: max |hard21 - slicing| (quark, gluon):', worst1
  ! two loop: scan in the CM angle at fixed Q, y, xp across the region
  ! boundaries (s_ij = s45 and s_ij = 0 lines); the largest jump of H^(2)/H^(0)
  ! between neighbouring points relative to the local slope
  block
    integer, parameter :: ns = 4000
    real(dp) :: hs(0:ns,2), t, jmp, d1, d2
    integer :: i, ix, ic, imax
    do ix = 1, 3
       xp = merge(0.1_dp, merge(0.4_dp, 0.8_dp, ix == 2), ix == 1)
       do i = 0, ns
          t = pi*(0.002_dp + 0.996_dp*i/ns)
          call breit_born(30.0_dp, 0.5_dp, xp, t, 0.7_dp, P)
          call hard21_eval(P, 1, -1.0_dp/3, h); hs(i,1) = h(2)
          call hard21_eval(P, 2, 2.0_dp/3, h); hs(i,2) = h(2)
       enddo
       do ic = 1, 2
          jmp = 0; imax = 0
          do i = 2, ns - 1
             d1 = abs(hs(i,ic) - hs(i-1,ic)); d2 = 0.5_dp*(abs(hs(i-1,ic) - hs(i-2,ic)) + abs(hs(i+1,ic) - hs(i,ic)))
             if (d1 / max(d2, 1e-12_dp*max(1.0_dp, abs(hs(i,ic)))) > jmp) then
                jmp = d1 / max(d2, 1e-12_dp*max(1.0_dp, abs(hs(i,ic)))); imax = i
             endif
          enddo
          print '(a,f4.1,a,i2,a,es10.2,a,f8.4,a,2es14.6)', 'scan xp', xp, ' chan', ic, '  largest step/neighbour steps', jmp, &
               '  at theta/pi', 0.002_dp + 0.996_dp*imax/ns, '  H2/H0 there', hs(imax-1,ic), hs(imax,ic)
       enddo
    enddo
  end block
contains
  ! Born of DIS 2+1 in the Breit frame (DISENT layout, px py pz E):
  ! q = (0,0,-Q,0), incoming parton along +z with momentum fraction
  ! xi = x/xp of the photon's light-cone momentum, the two outgoing partons
  ! back to back in their CM frame at angles th, ph
  subroutine breit_born(Q, y, xp, th, ph, P)
    real(dp), intent(in) :: Q, y, xp, th, ph
    real(dp), intent(out) :: P(4,7)
    real(dp) :: Ea, W(4), sh, E, k(4), gam, bz
    integer :: n
    P = 0
    Ea = Q/(2*xp)
    P(:,1) = [0.0_dp, 0.0_dp, Ea, Ea]
    P(:,5) = [0.0_dp, 0.0_dp, -Q, 0.0_dp]
    W = P(:,1) + P(:,5)
    sh = W(4)**2 - W(3)**2
    E = sqrt(sh)/2
    k = [E*sin(th)*cos(ph), E*sin(th)*sin(ph), E*cos(th), E]
    ! boost along z from the CM frame to the Breit frame
    bz = W(3)/W(4); gam = 1/sqrt(1 - bz**2)
    P(:,2) = k; P(:,3) = [-k(1), -k(2), -k(3), k(4)]
    do n = 2, 3
       E = P(4,n)
       P(4,n) = gam*(E + bz*P(3,n)); P(3,n) = gam*(P(3,n) + bz*E)
    enddo
    P(:,6) = Q/2*[2*sqrt(1 - y)/y, 0.0_dp, -1.0_dp, (2 - y)/y]
    P(:,7) = P(:,6) - P(:,5)
  end subroutine breit_born
end program harness_hard21
