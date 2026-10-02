!-----------------------------------------------------------------------
! DIS 3+1 tree matrix elements (photon exchange) from MCFM 10.3's Z+2 jet
! amplitudes crossed to DIS (dis31/mcfm, crossing rules in its README).
!
! me31(P, fl, msq): P(4,7) in DISENT's layout (1 incoming parton, 2-4
! outgoing partons, 5 q, 6 incoming lepton, 7 outgoing lepton; components
! px, py, pz, E); fl(1:4) the flavours of partons 1-4 (0 gluon, 1-5 d u s c
! b, negative antiquarks). msq: |M|^2 averaged over the spins and colours of
! the incoming lepton and parton, summed over the final state, with
! alpha = 1/137 and divided by (alpha_s/2pi)^2 (DISENT's MATFOR
! normalisation); no symmetry factors for identical final-state partons.
! Zero for flavour assignments that do not occur.
!-----------------------------------------------------------------------
module me31
  implicit none
  private
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: mxpart = 14
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: xn = 3, V = xn**2 - 1
  ! (4 pi alpha)^2 g_s^4 / (alpha_s/2pi)^2 = (4 pi/137)^2 64 pi^4
  real(dp), parameter :: cnorm = (4*pi/137.0_dp)**2*64*pi**4
  real(dp), parameter :: eq(5) = [-1.0_dp/3, 2.0_dp/3, -1.0_dp/3, 2.0_dp/3, -1.0_dp/3]
  integer, parameter :: swp(2) = [2, 1]
  public :: me31_tree
contains

  subroutine me31_tree(P, fl, msq)
    real(dp), intent(in) :: P(4,7)
    integer, intent(in) :: fl(4)
    real(dp), intent(out) :: msq
    complex(dp) :: za(mxpart,mxpart), zb(mxpart,mxpart)
    common /zprods/ za, zb
    real(dp) :: pm(mxpart,4), m2(2,2), avg, c1, c2
    complex(dp) :: A(2,2,2), B(2,2,2), Ad(2,2,2), Bd(2,2,2), Ae(2,2,2), Be(2,2,2), a1, b1
    integer :: nq, ng, i, iq, iqb, ig(2), j1, j2, j3, iother(3), k
    msq = 0
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
             if (fl(i) == fl(1)) k = i
          enddo
          if (k == 0) return
          call others(k, iother)
          pm(1,:) = -P(:,1); pm(2,:) = P(:,k)
          pm(5,:) = P(:,iother(1)); pm(6,:) = P(:,iother(2))
          call spinoru(6, pm, za, zb)
          ! for photon exchange the antiquark line has the same |M|^2
          ! (charge conjugation of the line)
          call z2jetsq(2, 1, 3, 4, 5, 6, za, zb, m2)
          c1 = eq(abs(fl(1)))
          msq = avg*cnorm*c1**2*V*xn/4*sum(m2)
       else
          ! g -> q qbar g: q, qbar, g among 2-4
          iq = 0; iqb = 0; ig = 0
          do i = 2, 4
             if (fl(i) > 0) iq = i
             if (fl(i) < 0) iqb = i
             if (fl(i) == 0) ig(1) = i
          enddo
          if (iq == 0 .or. iqb == 0 .or. fl(iq) /= -fl(iqb)) return
          pm(1,:) = -P(:,1); pm(2,:) = P(:,iq); pm(5,:) = P(:,iqb); pm(6,:) = P(:,ig(1))
          call spinoru(6, pm, za, zb)
          call z2jetsq(2, 5, 3, 4, 1, 6, za, zb, m2)
          c1 = eq(fl(iq))
          msq = avg*cnorm*c1**2*V*xn/4*sum(m2)
       endif
       return
    endif
    if (ng /= 0 .or. fl(1) == 0) return
    ! four (anti)quarks: incoming q (or qbar) of flavour fl(1), outgoing
    ! line partner, and a pair Q Qbar
    nq = fl(1)
    call four_quark(P, fl, msq)
    msq = avg*cnorm*msq
  end subroutine me31_tree

  ! four-quark |M|^2 without the average and normalisation: colour factor
  ! 4 V (MCFM's faclo), photon couplings
  subroutine four_quark(P, fl, s)
    real(dp), intent(in) :: P(4,7)
    integer, intent(in) :: fl(4)
    real(dp), intent(out) :: s
    complex(dp) :: za(mxpart,mxpart), zb(mxpart,mxpart)
    common /zprods/ za, zb
    real(dp) :: pm(mxpart,4), cq, cQ2, sg
    complex(dp) :: A(2,2,2), B(2,2,2), Ad(2,2,2), Bd(2,2,2), Ae(2,2,2), Be(2,2,2), aa, bb
    integer :: i, k, kq, kqb, j1, j2, j3, f0, nsame, ip(3)
    s = 0
    f0 = fl(1)
    ! work with an incoming quark: for an incoming antiquark conjugate all
    ! flavours (photon exchange: same |M|^2)
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
    cq = eq(abs(f0))
    if (fl(ip(1)) == f0 .and. fl(ip(2)) == f0) then
       ! identical quarks: the pair has the flavour of the incoming quark
       if (fl(kqb) /= -f0) return
       pm(2,:) = P(:,ip(1)); pm(5,:) = P(:,ip(2))
       call spinoru(6, pm, za, zb)
       call ampqqb_qqb(2, 1, 5, 6, A, B)
       call ampqqb_qqb(5, 1, 2, 6, Ae, Be)
       do j1 = 1, 2; do j2 = 1, 2; do j3 = 1, 2
          s = s + cq**2*(abs(A(j1,j2,j3))**2 + abs(B(j1,j2,j3))**2 + abs(Ae(j1,j2,j3))**2 + abs(Be(j1,j2,j3))**2)
       enddo; enddo; enddo
       ! interference from MCFM's own (phase-consistent) construction
       call ampqqb_qqb(1, 2, 6, 5, Ad, Bd)
       call ampqqb_qqb(1, 5, 2, 6, Ae, Be)
       do j1 = 1, 2; do j3 = 1, 2
          s = s + cq**2*2/xn*real((Ad(j1,swp(j1),j3) - Bd(j1,swp(j1),j3))*conjg(Ae(j1,j1,j3) + Be(j1,j1,j3)), dp)
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
    cQ2 = eq(abs(fl(kq)))
    do j1 = 1, 2; do j2 = 1, 2; do j3 = 1, 2
       s = s + abs(cq*A(j1,j2,j3) + cQ2*B(j1,j2,j3))**2
    enddo; enddo; enddo
    s = 4*V*s
  end subroutine four_quark

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
