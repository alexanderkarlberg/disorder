!----------------------------------------------------------------------
! NLO DIS 2+1 (O(alpha_s^2) coefficient of 2+1-jet observables) with
! 2-jettiness slicing, against DISENT's dipole-subtracted NLO, in one
! DISENT run (photon exchange, fixed x and Q^2, mu_R = mu_F = Q).
!
! For every DISENT event:
!  - reference: all O(alpha_s^2) contributions (VIRTHR, COLFOR, MATFOR,
!    SUBFOR counter-events) -> the full NLO coefficient;
!  - above the cut: the real four-parton events (MATFOR) with
!    tau_2 = T_2/Q > tau_cut. Counter-events, collinear terms and the
!    three-parton virtual have Born kinematics (tau_2 = 0) and drop out;
!  - below the cut: the three-parton Born (MATTHR) reweighted by the
!    one-loop SCET cumulant at T_cut (hard, two jet, beam and soft
!    functions; mod_slicing_scet), O(alpha_s^2).
! T_2 (geometric measure, Breit frame, proton along +z):
!   T_2 = min( min_j (E_j - p_zj), min_{j<k} (E_j + E_k - |p_j + p_k|) ).
! Observable: tau_zQ = 1 - (2/Q) sum_{p_z<0} |p_z| (= tau_1^b) in bins
! between 0.05 and 0.5, which vanish at 1+1 Born kinematics and avoid the
! IR-unsafe edge tau_zQ = 1 (empty current hemisphere).
!
! Usage: tau2_nlo -pdf NAME -x X -Q2 Q2 -s S -nev N -seed1 I -seed2 J [-npow1 2 -npow2 4]
!----------------------------------------------------------------------
module tau2_run
  use types, only: dp
  use mod_slicing_scet
  implicit none
  integer, parameter :: nt = 10, nb = 6
  real(dp), parameter :: taus(nt) = [2e-2_dp, 1e-2_dp, 5e-3_dp, 2e-3_dp, 1e-3_dp, 5e-4_dp, 2e-4_dp, 1e-4_dp, &
       & 3e-5_dp, 1e-5_dp]
  ! tau_zQ bins below 0.5 only: a Born with an empty current hemisphere has
  ! tau_zQ = 1 exactly, and a soft gluon into the current hemisphere moves it
  ! below 1, so any bin edge at tau_zQ = 1 is not IR safe (first runs).
  real(dp), parameter :: blo(nb) = [0.05_dp, 0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, 0.05_dp]
  real(dp), parameter :: bhi(nb) = [0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp, 0.5_dp]
  real(dp), save :: xfix, Q2fix
  ! per-event values and global sums (s1 = sum, s2 = sum of squares)
  real(dp), save :: ev_ref(nb) = 0, ev_ab(nt,nb) = 0, ev_be(nt,nb) = 0
  real(dp), save :: s1_ref(nb) = 0, s2_ref(nb) = 0, s1_ab(nt,nb) = 0, s2_ab(nt,nb) = 0
  real(dp), save :: s1_be(nt,nb) = 0, s2_be(nt,nb) = 0, s1_d(nt,nb) = 0, s2_d(nt,nb) = 0
  real(dp), save :: tau2max_ct = 0
  real(dp), save :: ev_born(nb) = 0, s1_born(nb) = 0, s2_born(nb) = 0   ! O(alpha_s) Born rate
  integer(8), save :: nevt = 0, nborn = 0, nabort = 0
  ! DISENT aborts the rest of an event when VIRTHR's collinear X has
  ! 1 - X < CUTOFF (alternate return), after the three-parton Born has
  ! been given; such events have no O(alpha_s^2) part, so their Born and
  ! below-cut weights are dropped too (keepall = .true.: old behaviour)
  logical, save :: ev_na2 = .false., keepall = .false.
  integer, save :: ndebug = 0
  ! diagnostics: rapidity of a boost along z applied before computing T_2
  ! and the SCET ingredients (geometric measure in another frame); the
  ! observable stays in the Breit frame
  real(dp), save :: boostY = 0
  ! diagnostics: all pairs (beam, jet 1, jet 2) at s_ij > smin in the Breit frame (no cut for 0)
  real(dp), save :: smin = 0
  ! measure: 0 = geometric (Q_i = 2 E_i in the Breit frame), 1 = invariant
  ! (Q_i = Q for all regions, frame independent):
  !   tau_2 = min( min_j 2 eta_j P.p_j, min_{j<k} 2 p_j.p_k ) / Q^2,
  !   eta_j = x (1 + s_kl/Q^2) (Born momentum fraction of the other two)
  ! with the cumulant at lambda_B = lambda_J = tau_cut and the soft function
  ! of the invariants s_ij = 2 p_i.p_j/Q^2 of the Born
  integer, save :: measure = 0
  real(dp), save :: sdis = 0
  ! diagnostics: Born weight (total bin) with min_ij s_ij below 1e-2, 1e-3, 1e-4,
  ! for the Breit-frame directions (geo) and the invariants 2p_i.p_j/Q^2 (inv)
  real(dp), save :: wsmall_geo(3) = 0, wsmall_inv(3) = 0, wborn_tot = 0
  real(dp), external :: DOT, LEIV, ERTV, alphasPDF

contains

  subroutine slice_cuts(s, xmn, xmx, q2mn, q2mx, ymn, ymx)
    real(dp), intent(in) :: s
    real(dp), intent(out) :: xmn, xmx, q2mn, q2mx, ymn, ymx
    xmn = xfix; xmx = xfix; q2mn = Q2fix; q2mx = Q2fix; ymn = 0; ymx = 1
  end subroutine slice_cuts

  subroutine slice_muf(p, s, scale)
    real(dp), intent(in) :: p(4,7), s
    real(dp), intent(out) :: scale
    scale = 1   ! (mu_F/Q)^2
  end subroutine slice_muf

  ! tau_zQ bins for the outgoing partons 2..n
  subroutine observable(n, p, Q, fb)
    integer, intent(in) :: n
    real(dp), intent(in) :: p(4,7), Q
    real(dp), intent(out) :: fb(nb)
    real(dp) :: t
    integer :: i
    t = 1
    do i = 2, n
       if (p(3,i) < 0) t = t + 2 * p(3,i) / Q
    enddo
    fb = merge(1.0_dp, 0.0_dp, t >= blo .and. t < bhi)
    if (smin > 0) then
       if (jet_min_shat(n, p) < smin) fb = 0
    endif
  end subroutine observable

  ! diagnostics (-smin): smallest s_ij = (1 - cos theta_ij)/2 among the beam
  ! and the two jets of the minimising 2-jettiness partition (Breit frame): for three
  ! partons the partons 2, 3; for four, the two left when one parton goes
  ! to the beam, or the merged pair and the third parton
  real(dp) function jet_min_shat(n, p) result(sm)
    integer, intent(in) :: n
    real(dp), intent(in) :: p(4,7)
    real(dp) :: best, v, ja(4), jb(4), jj(4,2)
    integer :: j, k, m
    if (n == 3) then
       ja = p(:,2); jb = p(:,3)
    else
       best = huge(1.0_dp)
       do j = 2, 4
          v = p(4,j) - p(3,j)
          if (v < best) then
             best = v; m = 0
             do k = 2, 4
                if (k /= j) then
                   m = m + 1; jj(:,m) = p(:,k)
                endif
             enddo
          endif
       enddo
       do j = 2, 4
          do k = j + 1, 4
             v = p(4,j) + p(4,k) - sqrt(sum((p(1:3,j) + p(1:3,k))**2))
             if (v < best) then
                best = v
                jj(:,1) = p(:,j) + p(:,k); jj(:,2) = p(:,9 - j - k)
             endif
          enddo
       enddo
       ja = jj(:,1); jb = jj(:,2)
    endif
    sm = min(shat1(ja), shat1(jb), shatjj(ja, jb))
  contains
    ! the two jets with each other (Breit-frame angle)
    real(dp) function shatjj(qa, qb)
      real(dp), intent(in) :: qa(4), qb(4)
      real(dp) :: a, b
      a = sqrt(sum(qa(1:3)**2)); b = sqrt(sum(qb(1:3)**2))
      shatjj = merge(0.5_dp * (1 - dot_product(qa(1:3), qb(1:3)) / (a * b)), 1.0_dp, a > 0 .and. b > 0)
    end function shatjj
    real(dp) function shat1(q)
      real(dp), intent(in) :: q(4)
      real(dp) :: a
      a = sqrt(sum(q(1:3)**2))
      shat1 = merge(0.5_dp * (1 - q(3) / a), 1.0_dp, a > 0)
    end function shat1
  end function jet_min_shat

  subroutine zboost(pin, pout)
    real(dp), intent(in) :: pin(4,7)
    real(dp), intent(out) :: pout(4,7)
    real(dp) :: ch, sh
    integer :: i
    ch = cosh(boostY); sh = sinh(boostY)
    pout = pin
    do i = 1, 7
       pout(4,i) = ch * pin(4,i) + sh * pin(3,i)
       pout(3,i) = sh * pin(4,i) + ch * pin(3,i)
    enddo
  end subroutine zboost

  real(dp) function tau2_of(n, pin) result(T)
    integer, intent(in) :: n
    real(dp), intent(in) :: pin(4,7)
    integer :: j, k
    real(dp) :: v(3), p(4,7), Q2, eta, ppj
    if (measure == 1) then
       ! invariant measure, returned as T = tau_2 Q (Q = sqrt(Q2))
       Q2 = -DOT(pin,5,5); eta = 2 * DOT(pin,1,6) / sdis
       T = huge(1.0_dp)
       do j = 2, n
          ppj = DOT(pin,1,j) / eta          ! P.p_j
          T = min(T, 2 * xfix * (1 + spair(j) / Q2) * ppj)
          do k = j + 1, n
             T = min(T, 2 * DOT(pin,j,k))
          enddo
       enddo
       T = max(T, 0.0_dp) / sqrt(Q2)
       return
    endif
    call zboost(pin, p)
    T = huge(1.0_dp)
    do j = 2, n
       T = min(T, p(4,j) - p(3,j))
       do k = j + 1, n
          v = p(1:3,j) + p(1:3,k)
          T = min(T, p(4,j) + p(4,k) - sqrt(sum(v * v)))
       enddo
    enddo
    T = max(T, 0.0_dp)
  contains
    ! 2 p_k.p_l summed over the pairs of outgoing partons other than j (the
    ! Born left when j goes to the beam; one pair for four partons)
    real(dp) function spair(j)
      integer, intent(in) :: j
      integer :: k, l
      spair = 0
      do k = 2, n
         do l = k + 1, n
            if (k /= j .and. l /= j) spair = spair + 2 * DOT(pin,k,l)
         enddo
      enddo
    end function spair
  end function tau2_of

  subroutine slice_user(n, na, ntyp, p, s, weight, scale2)
    integer, intent(in) :: n, na, ntyp
    real(dp), intent(in) :: s, p(4,7), weight(-6:6), scale2
    real(dp) :: eta, Q2, Q, xf(-6:6), as2pi, w, fb(nb), T, d
    integer :: it, ib
    if (n == 0) then
       nevt = nevt + 1
       if (.not. ev_na2) then
          if (any(ev_born /= 0)) nabort = nabort + 1
          if (.not. keepall) then
             ev_be = 0; ev_born = 0
          endif
       endif
       ev_na2 = .false.
       s1_ref = s1_ref + ev_ref; s2_ref = s2_ref + ev_ref**2
       s1_ab = s1_ab + ev_ab; s2_ab = s2_ab + ev_ab**2
       s1_be = s1_be + ev_be; s2_be = s2_be + ev_be**2
       do ib = 1, nb
          do it = 1, nt
             d = ev_be(it,ib) + ev_ab(it,ib) - ev_ref(ib)
             s1_d(it,ib) = s1_d(it,ib) + d; s2_d(it,ib) = s2_d(it,ib) + d * d
          enddo
       enddo
       s1_born = s1_born + ev_born; s2_born = s2_born + ev_born**2
       ev_ref = 0; ev_ab = 0; ev_be = 0; ev_born = 0
       return
    endif
    if (na == 2) ev_na2 = .true.
    sdis = s
    eta = 2 * DOT(p,1,6) / s
    Q2 = -DOT(p,5,5); Q = sqrt(Q2)
    call observable(n, p, Q, fb)
    if (all(fb == 0)) then
       if (na == 2 .and. n == 4 .and. ntyp /= 0) tau2max_ct = max(tau2max_ct, tau2_of(n, p) / Q)
       return
    endif
    call EvolvePDF(eta, Q, xf)
    call mask_pdf(xf)
    as2pi = alphasPDF(Q) / (2 * pi)
    w = dot_product(weight, xf) * as2pi**na
    if (na == 2) then
       ev_ref = ev_ref + w * fb
       if (n == 4) then
          T = tau2_of(n, p) / Q
          if (ntyp == 0) then
             do it = 1, nt
                if (T > taus(it)) ev_ab(it,:) = ev_ab(it,:) + w * fb
             enddo
          else
             tau2max_ct = max(tau2max_ct, T)
          endif
       endif
    elseif (na == 1 .and. n == 3) then
       ev_born = ev_born + w * fb
       call below_cut(p, s, weight, xf, eta, Q, as2pi, fb)
    endif
  end subroutine slice_user

  ! O(alpha_s^2) below-cut weight of a three-parton Born event
  subroutine below_cut(p, s, weight, xf, eta, Q, as2pi, fb)
    real(dp), intent(in) :: p(4,7), s, weight(-6:6), xf(-6:6), eta, Q, as2pi, fb(nb)
    real(dp) :: M(-6:6), r, Q2, l12, l13, l23, hq, hg, qqnf, ggnf, qqb, gqb, eq2sum
    real(dp) :: nhat(3,3), g(3,3), ls(3,3), casq(3), casg(3), ttq(3,3), ttg(3,3)
    real(dp) :: c0(-6:6), c1(-6:6), c2(-6:6), lb, wq, wg, tot, Ea, E2, E3, tc
    real(dp) :: EQ(-6:6), pb(4,7), sij(3,3)
    integer :: i, j, it, imax
    nborn = nborn + 1
    Q2 = Q * Q
    call MATTHR(p, M)
    ! common factor DISENT weight / matrix element
    imax = maxloc(abs(M), dim=1) - 7
    if (M(imax) == 0) return
    r = weight(imax) / M(imax)
    l12 = log(2 * DOT(p,1,2) / Q2); l13 = log(2 * DOT(p,1,3) / Q2); l23 = log(2 * DOT(p,2,3) / Q2)
    hq = hard_fact(.false., l12, l13, l23)
    hg = hard_fact(.true., l12, l13, l23)
    ! non-factorising one-loop part, as in DISENT's VIRTHR (photon exchange)
    qqnf = -((4 * pi / 137)**2 * 4 / Q2) * (2 * LEIV(p, p(1,6), 2, -1, 3) - Q2 / 2 * ERTV(p, 2, -1, 3))
    ggnf = TF / CF * ((4 * pi / 137)**2 * 4 / Q2) * (2 * LEIV(p, p(1,6), 2, 3, -1) - Q2 / 2 * ERTV(p, 2, 3, -1))
    EQ = 0
    do i = 1, 5
       EQ(i) = merge(2.0_dp/3, -1.0_dp/3, mod(i,2) == 0); EQ(-i) = -EQ(i)
    enddo
    qqb = M(1) / EQ(1)**2
    eq2sum = sum(EQ(1:5)**2)
    gqb = M(0) / eq2sum
    hq = hq + qqnf / qqb
    if (gqb /= 0) hg = hg + ggnf / gqb   ! (no gluon Born for TR = 0, diagnostics)
    ! soft function: directions (incoming parton along +z, partons 2 and 3),
    ! in the frame of the measure (Breit frame unless -boostY)
    ! edge diagnostics (Born weight in the total bin)
    block
      real(dp) :: sg, si, wb, u(3,3)
      integer :: a, b2
      wb = r * dot_product(M, xf) * as2pi * fb(nb)
      do a = 1, 3
         u(:,a) = p(1:3,a) / sqrt(sum(p(1:3,a)**2))
      enddo
      sg = huge(1.0_dp); si = huge(1.0_dp)
      do a = 1, 3
         do b2 = a + 1, 3
            sg = min(sg, 0.5_dp * (1 - dot_product(u(:,a), u(:,b2))))
            si = min(si, 2 * DOT(p,a,b2) / Q2)
         enddo
      enddo
      wborn_tot = wborn_tot + wb
      wsmall_geo = wsmall_geo + merge(wb, 0.0_dp, sg < [1e-2_dp, 1e-3_dp, 1e-4_dp])
      wsmall_inv = wsmall_inv + merge(wb, 0.0_dp, si < [1e-2_dp, 1e-3_dp, 1e-4_dp])
    end block
    call zboost(p, pb)
    if (measure == 1) then
       ! invariant measure: directions n_i = 2 p_i/Q, s_ij = 2 p_i.p_j/Q^2,
       ! and all logs at tau_cut (the energies below enter as 2 E/Q = 1)
       do i = 1, 3
          do j = 1, 3
             sij(i,j) = merge(0.0_dp, 2 * DOT(p,i,j) / Q2, i == j)
          enddo
       enddo
       call soft_from_s(3, sij, g, ls)
       pb(4,1:3) = Q / 2
    else
       nhat(:,1) = [0.0_dp, 0.0_dp, 1.0_dp]
       nhat(:,2) = pb(1:3,2) / sqrt(sum(pb(1:3,2)**2))
       nhat(:,3) = pb(1:3,3) / sqrt(sum(pb(1:3,3)**2))
       call soft_geom(3, nhat, g, ls)
    endif
    casq = [CF, CF, CA]; casg = [CA, CF, CF]
    call ttmat(casq, ttq); call ttmat(casg, ttg)
    ! beam function
    call beam_coeffs(eta, Q, c0, c1, c2)
    Ea = pb(4,1); E2 = pb(4,2); E3 = pb(4,3)
    do it = 1, nt
       tc = taus(it)
       lb = log(2 * Ea * tc / Q)
       ! quark channel: jets q (2) and g (3); gluon channel: q (2), qbar (3)
       wq = hq + jet_cum(.false., 2 * E2 * tc / Q) + jet_cum(.true., 2 * E3 * tc / Q) &
            & + soft_from_geom(3, g, ls, casq, ttq, tc)
       wg = hg + jet_cum(.false., 2 * E2 * tc / Q) + jet_cum(.false., 2 * E3 * tc / Q) &
            & + soft_from_geom(3, g, ls, casg, ttg, tc)
       tot = 0
       do i = -5, 5
          if (i == 0) then
             tot = tot + M(0) * (wg * xf(0) + c0(0) + c1(0) * lb + c2(0) * lb * lb)
          else
             tot = tot + M(i) * (wq * xf(i) + c0(i) + c1(i) * lb + c2(i) * lb * lb)
          endif
       enddo
       ev_be(it,:) = ev_be(it,:) + r * tot * as2pi**2 * fb
       if (ndebug > 0 .and. it == 5) then
          ndebug = ndebug - 1
          write(*,'(a,3f8.4,a,2f9.4,a,3es11.3)') ' DEBUG x_p,E2/Q,E3/Q ', Q*Q/(2*DOT(p,1,5))/1, E2/Q, E3/Q, &
               & '  tau_zQ-bin/eta ', sum(fb), eta, '  l12 l13 l23 ', l12, l13, l23
          write(*,'(a,es12.4,a,4f10.4)') '   tau_cut ', tc, '  quark: H(fact+nf) Jq Jg S ', hq, &
               & jet_cum(.false., 2 * E2 * tc / Q), jet_cum(.true., 2 * E3 * tc / Q), soft_from_geom(3, g, ls, casq, ttq, tc)
          write(*,'(a,4f10.4)') '                        gluon: H Jq Jq S ', hg, &
               & jet_cum(.false., 2 * E2 * tc / Q), jet_cum(.false., 2 * E3 * tc / Q), soft_from_geom(3, g, ls, casg, ttg, tc)
          write(*,'(a,2f10.4,a,2f10.4,a,f10.4)') '   nonfact q,g ', qqnf / qqb, ggnf / gqb, &
               & '  beam/xf (u, g) ', (c0(2) + c1(2) * lb + c2(2) * lb * lb) / xf(2), &
               & (c0(0) + c1(0) * lb + c2(0) * lb * lb) / xf(0), '  soft geom I-part q ', &
               & 0.5_dp * sum(ttq * (g - (ls**2 - zeta2) * merge(1.0_dp, 0.0_dp, ls /= 0)))
       endif
    enddo
  end subroutine below_cut

  subroutine ttmat(cas, tt)
    real(dp), intent(in) :: cas(3)
    real(dp), intent(out) :: tt(3,3)
    tt = 0
    tt(1,2) = 0.5_dp * (cas(3) - cas(1) - cas(2)); tt(2,1) = tt(1,2)
    tt(1,3) = 0.5_dp * (cas(2) - cas(1) - cas(3)); tt(3,1) = tt(1,3)
    tt(2,3) = 0.5_dp * (cas(1) - cas(2) - cas(3)); tt(3,2) = tt(2,3)
  end subroutine ttmat

  subroutine report()
    integer :: ib, it, u
    real(dp) :: n, m, e
    n = real(nevt, dp)
    ! raw sums for combining runs (tools: slicing/tau2_combine.py)
    open(newunit=u, file='tau2_nlo.dat', status='replace')
    write(u,*) nevt, nt, nb
    write(u,*) taus
    write(u,*) blo, bhi
    write(u,*) s1_ref, s2_ref
    write(u,*) s1_ab, s2_ab
    write(u,*) s1_be, s2_be
    write(u,*) s1_d, s2_d
    write(u,*) s1_born, s2_born
    close(u)
    write(*,'(a,i12,a,i12,a,es10.2)') ' events ', nevt, '  Born events reweighted ', nborn, &
         & '  max tau_2 of counter-events/collinear terms ', tau2max_ct
    write(*,'(a,i12,a,l2)') ' events in the bins without O(alpha_s^2) part (DISENT abort) ', nabort, &
         & '  kept ', keepall
    write(*,'(a,3es11.3,a,3es11.3)') ' Born fraction with min s_ij < 1e-2,1e-3,1e-4: Breit-frame angles', &
         & wsmall_geo / wborn_tot, '  invariants', wsmall_inv / wborn_tot
    do ib = 1, nb
       write(*,'(/a,f5.2,a,f5.2,a,es14.6,a,es10.3,a,es14.6)') ' tau_zQ in [', blo(ib), ',', bhi(ib), &
            & ')   NLO coefficient (DISENT) ', s1_ref(ib), ' +- ', err(s1_ref(ib), s2_ref(ib)), &
            & '   O(alpha_s) Born ', s1_born(ib)
       write(*,'(a)') '   tau_cut       below            above            sum           (sum-ref)/ref     +-'
       do it = 1, nt
          m = s1_d(it,ib); e = err(s1_d(it,ib), s2_d(it,ib))
          write(*,'(es10.1,3es17.6,2es14.3)') taus(it), s1_be(it,ib), s1_ab(it,ib), &
               & s1_be(it,ib) + s1_ab(it,ib), m / s1_ref(ib), e / abs(s1_ref(ib))
       enddo
    enddo
  contains
    real(dp) function err(a, b)
      real(dp), intent(in) :: a, b
      err = sqrt(max(b - a * a / n, 0.0_dp))
    end function err
  end subroutine report

end module tau2_run

program tau2_nlo
  use types, only: dp
  use mod_parameters, only: nflav, NC, CC, noZ, Zonly, intonly, neutrino, positron
  use tau2_run
  use mod_slicing_scet, only: soft_tol, pdf_mask, scet_set_colour, CF, CA, TF
  use sub_defs_io
  implicit none
  character(len=100) :: pdf
  real(dp) :: s, npow1, npow2, cf_in, ca_in, tf_in, cutoff
  integer :: nev, seed1, seed2
  external :: DISENTFULL

  pdf = string_val_opt('-pdf', 'NNPDF30_nlo_as_0118')
  xfix = dble_val_opt('-x', 0.01_dp)
  Q2fix = dble_val_opt('-Q2', 400.0_dp)
  s = dble_val_opt('-s', 4 * 27.5_dp * 920.0_dp)
  nev = int_val_opt('-nev', 100000)
  seed1 = int_val_opt('-seed1', 12345)
  seed2 = int_val_opt('-seed2', 67890)
  npow1 = dble_val_opt('-npow1', 2.0_dp)
  npow2 = dble_val_opt('-npow2', 4.0_dp)
  soft_tol = dble_val_opt('-softtol', 1e-9_dp)
  pdf_mask = int_val_opt('-pdfmask', 0)   ! 1: quarks only, 2: gluon only (diagnostics)
  ndebug = int_val_opt('-debug', 0)       ! print the cumulant pieces for this many Born events
  boostY = dble_val_opt('-boostY', 0.0_dp) ! measure in a frame boosted along z (diagnostics)
  smin = dble_val_opt('-smin', 0.0_dp)     ! jets away from the beam: s_1J > smin (diagnostics)
  measure = merge(1, 0, string_val_opt('-measure', 'geo') == 'inv')   ! geo (default) or inv
  keepall = log_val_opt('-keepall')        ! keep Born/below-cut of aborted events (old behaviour)
  cutoff = dble_val_opt('-cutoff', 1e-8_dp)  ! DISENT's technical cutoff
  ! colour factors (diagnostics: e.g. -CA 0 or -TR 0 isolate colour structures; DISENT gets the same)
  cf_in = dble_val_opt('-CF', 4.0_dp/3.0_dp)
  ca_in = dble_val_opt('-CA', 3.0_dp)
  tf_in = dble_val_opt('-TR', 0.5_dp)
  call scet_set_colour(cf_in, ca_in, tf_in)

  nflav = 5; NC = .true.; CC = .false.; noZ = .true.; Zonly = .false.
  intonly = .false.; neutrino = .false.; positron = .false.
  call InitPDFsetByName(trim(pdf))
  call InitPDF(0)
  write(*,'(a,a,a,f8.5,a,f9.2,a,f10.1,a,i10,a,i2,a,3f9.5,a,i2)') ' pdf ', trim(pdf), '  x ', xfix, '  Q2 ', Q2fix, '  s ', s, &
       & '  nev ', nev, '  pdfmask ', pdf_mask, '  CF CA TR ', CF, CA, TF, '  measure ', measure

  call DISENTFULL(nev, s, 5, slice_user, slice_cuts, seed1, seed2, npow1, npow2, &
       & cutoff, 2, slice_muf, CF, CA, TF, .false.)
  call report()
end program tau2_nlo
