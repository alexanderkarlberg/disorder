!-----------------------------------------------------------------------
! DIS 3+1 tree matrix elements (photon exchange; photon + Z with ew31) from
! MCFM 10.3's Z+2 jet amplitudes crossed to DIS (dis31/mcfm, crossing rules
! in its README).
!
! me31(P, fl, msq): P(4,7) in DISENT's layout (1 incoming parton, 2-4
! outgoing partons, 5 q, 6 incoming lepton, 7 outgoing lepton; components
! px, py, pz, E); fl(1:4) the flavours of partons 1-4 (0 gluon, 1-5 d u s c
! b, negative antiquarks). msq: |M|^2 averaged over the spins and colours of
! the incoming lepton and parton, summed over the final state, with
! alpha = 1/137 and divided by (alpha_s/2pi)^2 (DISENT's MATFOR
! normalisation); no symmetry factors for identical final-state partons.
! Zero for flavour assignments that do not occur. Couplings per helicity of
! the quark and lepton lines from ew31 (8 Oct 2026); with photon + Z the
! interference of the boson on different quark lines for different
! flavours is dropped, as in disorder's MATFOR and the structure functions
! (odd for vector couplings; for axial ones proportional to the sum of the
! axial couplings of the pair flavours), also its pair flavour = incoming
! flavour member in the identical-quark |D|^2, |E|^2; the direct-exchange
! interference of identical quarks is kept.
! W exchange (ew31_mode = 2, 9 Oct): the boson line changes flavour
! (ew31_out); four quarks: four_quark_cc (W on the incoming line or on the
! pair, the two assignments interfering; = disorder's CC MATFOR with
! ew31_ccmatfor, which drops the W-on-the-pair identical-quark
! interference, tests/harness_me31 -noNC -CC).
!-----------------------------------------------------------------------
module me31
  use ew31
  implicit none
  private
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: mxpart = 14
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: xn = 3, V = xn**2 - 1
  ! (4 pi alpha)^2 g_s^4 / (alpha_s/2pi)^2 = (4 pi/137)^2 64 pi^4
  real(dp), parameter :: cnorm = (4*pi/137.0_dp)**2*64*pi**4
  integer, parameter :: swp(2) = [2, 1]
  ! diagnostic: keep the interference of the boson on different quark lines
  ! with photon + Z (dis31/tests/axsum31: size of the dropped terms)
  logical, public :: me31_keepint = .false.
  public :: me31_tree
contains

  subroutine me31_tree(P, fl, msq)
    real(dp), intent(in) :: P(4,7)
    integer, intent(in) :: fl(4)
    real(dp), intent(out) :: msq
    complex(dp) :: za(mxpart,mxpart), zb(mxpart,mxpart)
    common /zprods/ za, zb
    real(dp) :: pm(mxpart,4), m2(2,2), avg, c(2,2), Q2
    integer :: nq, ng, i, iq, iqb, ig(2), iother(3), k
    msq = 0
    Q2 = sum(P(1:3,5)**2) - P(4,5)**2
    ng = count(fl == 0)
    if (fl(1) == 0) then
       avg = 1.0_dp/(2*2*V)          ! lepton spin, gluon polarisation and colour
    else
       avg = 1.0_dp/(2*2*xn)         ! lepton spin, quark spin and colour
    endif
    pm = 0
    pm(3,:) = P(:,7); pm(4,:) = -P(:,6)
    if (ng == 2) then
       ! q g g (incoming quark or antiquark) or g -> q qbar g
       if (fl(1) /= 0) then
          ! the outgoing (anti)quark of the incoming line: same flavour
          k = 0
          do i = 2, 4
             if (fl(i) == ew31_out(fl(1))) k = i
          enddo
          if (k == 0 .or. ew31_out(fl(1)) == 0) return
          call others(k, iother)
          pm(1,:) = -P(:,1); pm(2,:) = P(:,k)
          pm(5,:) = P(:,iother(1)); pm(6,:) = P(:,iother(2))
          call spinoru(6, pm, za, zb)
          ! an antiquark line: the quark-helicity label exchanged (ew31)
          call z2jetsq(2, 1, 3, 4, 5, 6, za, zb, m2)
          call ew31_cpl(fl(1), Q2, c)
          msq = avg*cnorm*V*xn/4*sum(c**2*m2)
       else
          ! g -> q qbar g: q, qbar, g among 2-4
          iq = 0; iqb = 0; ig = 0
          do i = 2, 4
             if (fl(i) > 0) iq = i
             if (fl(i) < 0) iqb = i
             if (fl(i) == 0) ig(1) = i
          enddo
          if (iq == 0 .or. iqb == 0 .or. fl(iq) /= ew31_out(-fl(iqb))) return
          pm(1,:) = -P(:,1); pm(2,:) = P(:,iq); pm(5,:) = P(:,iqb); pm(6,:) = P(:,ig(1))
          call spinoru(6, pm, za, zb)
          call z2jetsq(2, 5, 3, 4, 1, 6, za, zb, m2)
          call ew31_cpl(fl(iq), Q2, c)
          msq = avg*cnorm*V*xn/4*sum(c**2*m2)
       endif
       return
    endif
    if (ng /= 0 .or. fl(1) == 0) return
    ! four (anti)quarks: incoming q (or qbar) of flavour fl(1), outgoing
    ! line partner, and a pair Q Qbar
    nq = fl(1)
    if (ew31_mode == 2) then
       call four_quark_cc(P, fl, Q2, msq)
    else
       call four_quark(P, fl, Q2, msq)
    endif
    msq = avg*cnorm*msq
  end subroutine me31_tree

  ! four-quark |M|^2 without the average and normalisation: colour factor
  ! 4 V (MCFM's faclo); A(j1,j2,j3): boson on the line (2,1), B: on (5,6),
  ! j1, j2 their helicities, j3 the lepton's. The line (5,6) is read in
  ! MCFM's orientation (outgoing quark in its outgoing slot), so its
  ! couplings are ew31's with the quark-helicity label exchanged back
  ! (cr, cQ2; fixed against disorder's NC MATFOR, harness_me31)
  subroutine four_quark(P, fl, Q2, s)
    real(dp), intent(in) :: P(4,7), Q2
    integer, intent(in) :: fl(4)
    real(dp), intent(out) :: s
    complex(dp) :: za(mxpart,mxpart), zb(mxpart,mxpart)
    common /zprods/ za, zb
    real(dp) :: pm(mxpart,4), cq(2,2), cr(2,2), cQ2(2,2), sg
    complex(dp) :: A(2,2,2), B(2,2,2), Ae(2,2,2), Be(2,2,2), D(2,2,2), E(2,2,2)
    integer :: i, k, kq, kqb, j1, j2, j3, f0, nsame, ip(3)
    s = 0
    f0 = fl(1)
    ! work with an incoming quark: for an incoming antiquark conjugate all
    ! flavours (the helicity labels of the antiquark lines exchanged, ew31)
    sg = merge(1, -1, f0 > 0)
    ! the outgoing partons: two with the sign of the incoming parton, one
    ! with the opposite sign
    nsame = 0; kqb = 0
    do i = 2, 4
       if (sg*fl(i) > 0) then
          nsame = nsame + 1; ip(nsame) = i
       elseif (sg*fl(i) < 0) then
          kqb = i
       endif
    enddo
    if (nsame /= 2 .or. kqb == 0) return
    pm = 0
    pm(1,:) = -P(:,1); pm(3,:) = P(:,7); pm(4,:) = -P(:,6); pm(6,:) = P(:,kqb)
    call ew31_cpl(f0, Q2, cq)
    cr = cq([2, 1],:)
    if (fl(ip(1)) == f0 .and. fl(ip(2)) == f0) then
       ! identical quarks: direct D (lines 2-1, 5-6) and exchange E (5-1, 2-6)
       ! amplitudes, interfering for opposite helicity labels j1, j2 (fixed
       ! against Feynman diagrams, dis31/tests/fd31.py)
       if (fl(kqb) /= -f0) return
       pm(2,:) = P(:,ip(1)); pm(5,:) = P(:,ip(2))
       call spinoru(6, pm, za, zb)
       call ampqqb_qqb(2, 1, 5, 6, A, B)
       call ampqqb_qqb(5, 1, 2, 6, Ae, Be)
       do j1 = 1, 2; do j2 = 1, 2; do j3 = 1, 2
          D(j1,j2,j3) = cq(j1,j3)*A(j1,j2,j3) - cr(j2,j3)*B(j1,j2,j3)
          E(j1,j2,j3) = cq(j1,j3)*Ae(j1,j2,j3) - cr(j2,j3)*Be(j1,j2,j3)
       enddo; enddo; enddo
       if (ew31_mode == 0 .or. me31_keepint) then
          s = sum(abs(D)**2) + sum(abs(E)**2)
       else
          ! photon + Z: |D|^2, |E|^2 without the interference of the boson on
          ! the two lines (the pair flavour = q member of the dropped sum)
          s = 0
          do j1 = 1, 2; do j2 = 1, 2; do j3 = 1, 2
             s = s + abs(cq(j1,j3)*A(j1,j2,j3))**2 + abs(cr(j2,j3)*B(j1,j2,j3))**2 &
                  & + abs(cq(j1,j3)*Ae(j1,j2,j3))**2 + abs(cr(j2,j3)*Be(j1,j2,j3))**2
          enddo; enddo; enddo
       endif
       do j1 = 1, 2; do j3 = 1, 2
          s = s + 2/xn*real(D(j1,swp(j1),j3)*conjg(E(j1,swp(j1),j3)), dp)
       enddo; enddo
       s = 4*V*s
       return
    endif
    ! different flavours: one outgoing parton continues the incoming line
    if (fl(ip(1)) == f0) then
       k = ip(1); kq = ip(2)
    elseif (fl(ip(2)) == f0) then
       k = ip(2); kq = ip(1)
    else
       return
    endif
    if (fl(kqb) /= -fl(kq)) return
    pm(2,:) = P(:,k); pm(5,:) = P(:,kq)
    call spinoru(6, pm, za, zb)
    call ampqqb_qqb(2, 1, 5, 6, A, B)
    call ew31_cpl(fl(kq), Q2, cQ2)
    cQ2 = cQ2([2, 1],:)
    ! the incoming line is read as (2,1), against MCFM's orientation (1,2)
    ! (qqb_z2jet): the photon on the pair line (B) has the opposite sign
    ! (charge-odd e_q e_Q term; checked against Feynman diagrams, fd31.py)
    do j1 = 1, 2; do j2 = 1, 2; do j3 = 1, 2
       if (ew31_mode == 0 .or. me31_keepint) then
          s = s + abs(cq(j1,j3)*A(j1,j2,j3) - cQ2(j2,j3)*B(j1,j2,j3))**2
       else
          s = s + abs(cq(j1,j3)*A(j1,j2,j3))**2 + abs(cQ2(j2,j3)*B(j1,j2,j3))**2
       endif
    enddo; enddo; enddo
    s = 4*V*s
  end subroutine four_quark

  ! W exchange (9 Oct): the two same-sign outgoing quarks ip(1), ip(2) and
  ! the antiquark kqb; assignment ia puts ip(ia) on the incoming line. The W
  ! sits on the incoming line if that quark is ew31_out of the incoming one
  ! (A), or on the pair line if the incoming quark continues and the pair
  ! couples to the W (B). The two assignments interfere as the identical-
  ! quark D, E of four_quark: W on the incoming line in both (EXX of
  ! disorder's MATFOR, identical quarks), on the line in one and on the pair
  ! in the other (EXY), or on the pair in both (identical quarks; not in
  ! MATFOR: dropped with ew31_ccmatfor).
  subroutine four_quark_cc(P, fl, Q2, s)
    real(dp), intent(in) :: P(4,7), Q2
    integer, intent(in) :: fl(4)
    real(dp), intent(out) :: s
    complex(dp) :: za(mxpart,mxpart), zb(mxpart,mxpart)
    common /zprods/ za, zb
    real(dp) :: pm(mxpart,4), cq(2,2), cl(2,2,2), cp(2,2,2), sg
    complex(dp) :: A(2,2,2), B(2,2,2), Ae(2,2,2), Be(2,2,2), M(2,2,2,2)
    integer :: i, f0, nsame, ip(3), kqb, ia, k, kq, j1, j2, j3
    logical :: on(2), bb
    s = 0
    f0 = fl(1)
    sg = merge(1, -1, f0 > 0)
    nsame = 0; kqb = 0
    do i = 2, 4
       if (sg*fl(i) > 0) then
          nsame = nsame + 1; ip(nsame) = i
       elseif (sg*fl(i) < 0) then
          kqb = i
       endif
    enddo
    if (nsame /= 2 .or. kqb == 0) return
    call ew31_cpl(f0, Q2, cq)
    cl = 0; cp = 0; on = .false.; bb = .true.
    do ia = 1, 2
       k = ip(ia); kq = ip(3 - ia)
       if (ew31_out(f0) /= 0 .and. fl(k) == ew31_out(f0) .and. fl(kqb) == -fl(kq)) then
          cl(:,:,ia) = cq; on(ia) = .true.; bb = .false.
       elseif (fl(k) == f0 .and. ew31_out(-fl(kqb)) == fl(kq) .and. fl(kq) /= 0) then
          call ew31_cpl(fl(kq), Q2, cp(:,:,ia))
          cp(:,:,ia) = cp([2, 1],:,ia); on(ia) = .true.
       endif
    enddo
    if (.not. any(on)) return
    pm = 0
    pm(1,:) = -P(:,1); pm(3,:) = P(:,7); pm(4,:) = -P(:,6); pm(6,:) = P(:,kqb)
    pm(2,:) = P(:,ip(1)); pm(5,:) = P(:,ip(2))
    call spinoru(6, pm, za, zb)
    call ampqqb_qqb(2, 1, 5, 6, A, B)
    call ampqqb_qqb(5, 1, 2, 6, Ae, Be)
    do j1 = 1, 2; do j2 = 1, 2; do j3 = 1, 2
       M(j1,j2,j3,1) = cl(j1,j3,1)*A(j1,j2,j3) - cp(j2,j3,1)*B(j1,j2,j3)
       M(j1,j2,j3,2) = cl(j1,j3,2)*Ae(j1,j2,j3) - cp(j2,j3,2)*Be(j1,j2,j3)
    enddo; enddo; enddo
    s = sum(abs(M)**2)
    if (all(on) .and. .not. (bb .and. ew31_ccmatfor)) then
       do j1 = 1, 2; do j3 = 1, 2
          s = s + 2/xn*real(M(j1,swp(j1),j3,1)*conjg(M(j1,swp(j1),j3,2)), dp)
       enddo; enddo
    endif
    s = 4*V*s
  end subroutine four_quark_cc

  subroutine others(k, o)
    integer, intent(in) :: k
    integer, intent(out) :: o(3)
    select case (k)
    case (2); o = [3, 4, 0]
    case (3); o = [2, 4, 0]
    case default; o = [2, 3, 0]
    end select
  end subroutine others
end module me31
