!-----------------------------------------------------------------------
! DIS 4+1 tree matrix elements (photon exchange; photon + Z with ew31,
! conventions as me31) from MCFM 10.3's Z+3
! parton amplitudes (xzqqggg: q qbar g g g; msq_gqqQQg, the photon version
! of msq_ZqqQQg: q qbar Q Qbar g) crossed to DIS (dis31/mcfm, crossing
! rules in its README).
!
! me41_tree(P, fl, msq): P(4,8) with 1 the incoming parton, 2-5 the
! outgoing partons, 6 q, 7 the incoming lepton, 8 the outgoing lepton
! (components px, py, pz, E); fl(1:5) the flavours of partons 1-5 (0
! gluon, 1-5 d u s c b, negative antiquarks). msq: |M|^2 averaged over the
! spins and colours of the incoming lepton and parton, summed over the
! final state, with alpha = 1/137 and divided by (alpha_s/2pi)^3 (as me31
! and DISENT's MATFOR, one power of alpha_s/2pi more); no symmetry factors
! for identical final-state partons. Zero for flavour assignments that do
! not occur. W exchange (9 Oct): as me31 (four_quark_cc, gluon_cc;
! msq_gqqQQg with the couplings of the exchanged pairing); checked in all
! single-unresolved limits against me31 (tests/harness_lim41, EW31 = "2 l").
!-----------------------------------------------------------------------
module me41
  use ew31
  implicit none
  private
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: mxpart = 14
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: xn = 3, V = xn**2 - 1
  ! (4 pi alpha)^2 g_s^6 / (alpha_s/2pi)^3 = (4 pi/137)^2 (8 pi^2)^3
  real(dp), parameter :: cnorm = (4*pi/137.0_dp)**2*(8*pi**2)**3
  ! MCFM's average over two incoming gluons in xzqqggg (spinave/V^2)
  real(dp), parameter :: avegg = 0.25_dp/V**2
  public :: me41_tree
contains

  subroutine me41_tree(P, fl, msq)
    real(dp), intent(in) :: P(4,8)
    integer, intent(in) :: fl(5)
    real(dp), intent(out) :: msq
    complex(dp) :: za(mxpart,mxpart), zb(mxpart,mxpart)
    common /zprods/ za, zb
    real(dp) :: gsq, as, ason2pi, ason4pi, Gf, gw, xw, gwsq, esq, vevsq
    common /qcdcouple/ gsq, as, ason2pi, ason4pi
    common /ewcouple/ Gf, gw, xw, gwsq, esq, vevsq
    integer :: colourchoice
    common /ColC/ colourchoice
    real(dp) :: pm(mxpart,4), m2(2,2), avg, MN, MI, c(2,2), cb(2,2), Q2
    integer :: ng, i, k, n, ig(3), iq, iqb, sg, nsame, ip(4), kqb(2), kg, ipos(2), ineg(2)
    msq = 0
    Q2 = sum(P(1:3,6)**2) - P(4,6)**2
    ! couplings of the MCFM routines: g_s = e = 1, all colour structures
    gsq = 1; esq = 1; colourchoice = 0
    ng = count(fl == 0)
    if (fl(1) == 0) then
       avg = 1.0_dp/(2*2*V)          ! lepton spin, gluon polarisation and colour
    else
       avg = 1.0_dp/(2*2*xn)         ! lepton spin, quark spin and colour
    endif
    pm = 0
    pm(3,:) = P(:,8); pm(4,:) = -P(:,7)
    if (ng == 3) then
       if (fl(1) /= 0) then
          ! q g g g (incoming quark or antiquark; the antiquark line with
          ! the quark-helicity label exchanged, ew31)
          k = 0
          do i = 2, 5
             if (fl(i) == ew31_out(fl(1))) k = i
          enddo
          if (k == 0 .or. ew31_out(fl(1)) == 0) return
          n = 0
          do i = 2, 5
             if (i /= k) then
                n = n + 1; ig(n) = i
             endif
          enddo
          pm(1,:) = -P(:,1); pm(2,:) = P(:,k)
          pm(5,:) = P(:,ig(1)); pm(6,:) = P(:,ig(2)); pm(7,:) = P(:,ig(3))
          call spinoru(7, pm, za, zb)
          call xzqqggg(2, 5, 6, 7, 1, 3, 4, m2)
          call ew31_cpl(fl(1), Q2, c)
       else
          ! g -> q qbar g g
          iq = 0; iqb = 0; n = 0
          do i = 2, 5
             if (fl(i) > 0) iq = i
             if (fl(i) < 0) iqb = i
             if (fl(i) == 0) then
                n = n + 1; ig(n) = i
             endif
          enddo
          if (iq == 0 .or. iqb == 0 .or. fl(iq) /= ew31_out(-fl(iqb))) return
          pm(1,:) = -P(:,1); pm(2,:) = P(:,iq); pm(5,:) = P(:,iqb)
          pm(6,:) = P(:,ig(1)); pm(7,:) = P(:,ig(2))
          call spinoru(7, pm, za, zb)
          call xzqqggg(2, 1, 6, 7, 5, 3, 4, m2)
          call ew31_cpl(fl(iq), Q2, c)
       endif
       msq = avg*cnorm*sum(c**2*m2)/avegg
       return
    endif
    if (ng /= 1) return
    if (fl(1) /= 0) then
       ! q -> q Q Qbar g (Q = q: identical quarks); an incoming antiquark by
       ! conjugating all flavours (couplings of the conjugated lines, ew31).
       ! The leptons enter as (4, 3): lepton-helicity label exchanged
       sg = merge(1, -1, fl(1) > 0)
       nsame = 0; n = 0; kg = 0
       do i = 2, 5
          if (fl(i) == 0) then
             kg = i
          elseif (sg*fl(i) > 0) then
             nsame = nsame + 1; ip(nsame) = i
          else
             n = n + 1; kqb(n) = i
          endif
       enddo
       if (nsame /= 2 .or. n /= 1) return
       pm(1,:) = -P(:,1); pm(6,:) = P(:,kqb(1)); pm(7,:) = P(:,kg)
       call ew31_cpl(fl(1), Q2, c); c = c(:,[2, 1])
       if (ew31_mode == 2) then
          call four_quark_cc()
       elseif (fl(ip(1)) == fl(1) .and. fl(ip(2)) == fl(1)) then
          if (fl(kqb(1)) /= -fl(1)) return
          pm(2,:) = P(:,ip(1)); pm(5,:) = P(:,ip(2))
          call spinoru(7, pm, za, zb)
          call msq_gqqQQg(2, 1, 5, 6, 7, 4, 3, c, c, c, c, ew31_mode /= 0, MN, MI)
          msq = MI
       else
          if (fl(ip(1)) == fl(1)) then
             k = ip(1); iq = ip(2)
          elseif (fl(ip(2)) == fl(1)) then
             k = ip(2); iq = ip(1)
          else
             return
          endif
          if (fl(kqb(1)) /= -fl(iq)) return
          pm(2,:) = P(:,k); pm(5,:) = P(:,iq)
          call spinoru(7, pm, za, zb)
          call ew31_cpl(fl(iq), Q2, cb); cb = cb(:,[2, 1])
          call msq_gqqQQg(2, 1, 5, 6, 7, 4, 3, c, cb, cb, c, ew31_mode /= 0, MN, MI)
          msq = MN
       endif
    else
       ! g -> q qbar Q Qbar (Q = q: identical quarks)
       nsame = 0; n = 0
       do i = 2, 5
          if (fl(i) > 0) then
             nsame = nsame + 1; if (nsame <= 2) ipos(nsame) = i
          elseif (fl(i) < 0) then
             n = n + 1; if (n <= 2) ineg(n) = i
          endif
       enddo
       if (nsame /= 2 .or. n /= 2) return
       if (ew31_mode == 2) then
          call gluon_cc()
          msq = avg*cnorm*32*msq
          return
       endif
       ! pair each quark with the antiquark of its flavour
       if (fl(ineg(1)) /= -fl(ipos(1))) ineg = ineg([2, 1])
       if (fl(ineg(1)) /= -fl(ipos(1)) .or. fl(ineg(2)) /= -fl(ipos(2))) return
       pm(1,:) = -P(:,1)
       pm(2,:) = P(:,ipos(1)); pm(5,:) = P(:,ineg(1))
       pm(6,:) = P(:,ipos(2)); pm(7,:) = P(:,ineg(2))
       call spinoru(7, pm, za, zb)
       call ew31_cpl(fl(ipos(1)), Q2, c); c = c(:,[2, 1])
       call ew31_cpl(fl(ipos(2)), Q2, cb); cb = cb(:,[2, 1])
       call msq_gqqQQg(2, 5, 6, 7, 1, 4, 3, c, cb, cb, c, ew31_mode /= 0, MN, MI)
       msq = merge(MI, MN, fl(ipos(1)) == fl(ipos(2)))
    endif
    ! MCFM's four-quark normalisation (qqb_z2jet_g): 4 g^6 e^4 8 MN
    msq = avg*cnorm*32*msq
  contains
    ! W exchange (9 Oct), as me31's four_quark_cc: assignment ia puts ip(ia)
    ! on the incoming line; the W on the line (cl) if that quark is
    ! ew31_out of the incoming one, on the pair (cp) if the incoming quark
    ! continues; both assignments: MI with the exchanged pairing's
    ! couplings (cl, cp of assignment 2 on the lines 5-1, 2-6)
    subroutine four_quark_cc()
      real(dp) :: cl(2,2,2), cp(2,2,2), M1, M2
      integer :: ia, kq
      logical :: on(2), bb
      cl = 0; cp = 0; on = .false.; bb = .true.
      do ia = 1, 2
         k = ip(ia); kq = ip(3 - ia)
         if (ew31_out(fl(1)) /= 0 .and. fl(k) == ew31_out(fl(1)) .and. fl(kqb(1)) == -fl(kq)) then
            cl(:,:,ia) = c; on(ia) = .true.; bb = .false.
         elseif (fl(k) == fl(1) .and. ew31_out(-fl(kqb(1))) == fl(kq) .and. fl(kq) /= 0) then
            call ew31_cpl(fl(kq), Q2, cp(:,:,ia)); cp(:,:,ia) = cp(:,[2, 1],ia); on(ia) = .true.
         endif
      enddo
      if (.not. any(on)) return
      if (.not. on(1)) then
         ip(1:2) = ip([2, 1]); cl = cl(:,:,[2, 1]); cp = cp(:,:,[2, 1]); on = on([2, 1])
      endif
      pm(2,:) = P(:,ip(1)); pm(5,:) = P(:,ip(2))
      call spinoru(7, pm, za, zb)
      call msq_gqqQQg(2, 1, 5, 6, 7, 4, 3, cl(:,:,1), cp(:,:,1), cl(:,:,2), cp(:,:,2), .false., MN, MI)
      if (.not. on(2)) then
         msq = MN
      elseif (bb .and. ew31_ccmatfor) then
         ! the two W-on-the-pair assignments without their interference
         M1 = MN
         call msq_gqqQQg(5, 1, 2, 6, 7, 4, 3, cl(:,:,2), cp(:,:,2), cl(:,:,1), cp(:,:,1), .false., M2, MI)
         msq = M1 + M2
      else
         msq = MI
      endif
    end subroutine four_quark_cc
    ! W exchange: g -> q qbar Q Qbar with the W on one pair, the other
    ! neutral; pairings (q1 qb1, q2 qb2) and (q2 qb1, q1 qb2)
    subroutine gluon_cc()
      real(dp) :: ce(2,2,4)
      integer :: kk
      logical :: on(2)
      ce = 0; on = .false.
      ! line couplings in msq_gqqQQg's order: q1-qb1, q2-qb2, q2-qb1, q1-qb2
      call wpair(ipos(1), ineg(1), ipos(2), ineg(2), ce(:,:,1))
      call wpair(ipos(2), ineg(2), ipos(1), ineg(1), ce(:,:,2))
      call wpair(ipos(2), ineg(1), ipos(1), ineg(2), ce(:,:,3))
      call wpair(ipos(1), ineg(2), ipos(2), ineg(1), ce(:,:,4))
      on(1) = any(ce(:,:,1:2) /= 0); on(2) = any(ce(:,:,3:4) /= 0)
      if (.not. any(on)) return
      if (.not. on(1)) then
         ineg = ineg([2, 1]); ce = ce(:,:,[4, 3, 2, 1]); on = on([2, 1])
      endif
      pm(1,:) = -P(:,1)
      pm(2,:) = P(:,ipos(1)); pm(5,:) = P(:,ineg(1))
      pm(6,:) = P(:,ipos(2)); pm(7,:) = P(:,ineg(2))
      call spinoru(7, pm, za, zb)
      call msq_gqqQQg(2, 5, 6, 7, 1, 4, 3, ce(:,:,1), ce(:,:,2), ce(:,:,3), ce(:,:,4), .false., MN, MI)
      kk = count(on)
      msq = merge(MI, MN, kk == 2)
    end subroutine gluon_cc
    ! the coupling of the W on the line (iq, iqb) if (jq, jqb) is a neutral
    ! pair, else zero
    subroutine wpair(iq, iqb, jq, jqb, cw)
      integer, intent(in) :: iq, iqb, jq, jqb
      real(dp), intent(out) :: cw(2,2)
      cw = 0
      if (fl(jqb) /= -fl(jq)) return
      if (ew31_out(-fl(iqb)) /= fl(iq)) return
      call ew31_cpl(fl(iq), Q2, cw); cw = cw(:,[2, 1])
    end subroutine wpair
  end subroutine me41_tree
end module me41
