!-----------------------------------------------------------------------
! DIS 4+1 tree matrix elements (photon exchange) from MCFM 10.3's Z+3
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
! not occur.
!-----------------------------------------------------------------------
module me41
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
  real(dp), parameter :: eq(5) = [-1.0_dp/3, 2.0_dp/3, -1.0_dp/3, 2.0_dp/3, -1.0_dp/3]
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
    real(dp) :: pm(mxpart,4), m2(2,2), avg, MN, MI, cq
    integer :: ng, i, k, n, ig(3), iq, iqb, sg, nsame, ip(4), kqb(2), kg, ipos(2), ineg(2)
    msq = 0
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
          ! q g g g (incoming quark or antiquark; photon exchange: the
          ! antiquark line has the same |M|^2)
          k = 0
          do i = 2, 5
             if (fl(i) == fl(1)) k = i
          enddo
          if (k == 0) return
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
          cq = eq(abs(fl(1)))
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
          if (iq == 0 .or. iqb == 0 .or. fl(iq) /= -fl(iqb)) return
          pm(1,:) = -P(:,1); pm(2,:) = P(:,iq); pm(5,:) = P(:,iqb)
          pm(6,:) = P(:,ig(1)); pm(7,:) = P(:,ig(2))
          call spinoru(7, pm, za, zb)
          call xzqqggg(2, 1, 6, 7, 5, 3, 4, m2)
          cq = eq(fl(iq))
       endif
       msq = avg*cnorm*cq**2*sum(m2)/avegg
       return
    endif
    if (ng /= 1) return
    if (fl(1) /= 0) then
       ! q -> q Q Qbar g (Q = q: identical quarks); an incoming antiquark by
       ! conjugating all flavours (photon exchange: same |M|^2)
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
       cq = eq(abs(fl(1)))
       if (fl(ip(1)) == fl(1) .and. fl(ip(2)) == fl(1)) then
          if (fl(kqb(1)) /= -fl(1)) return
          pm(2,:) = P(:,ip(1)); pm(5,:) = P(:,ip(2))
          call spinoru(7, pm, za, zb)
          call msq_gqqQQg(2, 1, 5, 6, 7, 4, 3, cq, cq, MN, MI)
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
          call msq_gqqQQg(2, 1, 5, 6, 7, 4, 3, cq, eq(abs(fl(iq))), MN, MI)
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
       ! pair each quark with the antiquark of its flavour
       if (fl(ineg(1)) /= -fl(ipos(1))) ineg = ineg([2, 1])
       if (fl(ineg(1)) /= -fl(ipos(1)) .or. fl(ineg(2)) /= -fl(ipos(2))) return
       pm(1,:) = -P(:,1)
       pm(2,:) = P(:,ipos(1)); pm(5,:) = P(:,ineg(1))
       pm(6,:) = P(:,ipos(2)); pm(7,:) = P(:,ineg(2))
       call spinoru(7, pm, za, zb)
       if (fl(ipos(1)) == fl(ipos(2))) then
          cq = eq(fl(ipos(1)))
          call msq_gqqQQg(2, 5, 6, 7, 1, 4, 3, cq, cq, MN, MI)
          msq = MI
       else
          call msq_gqqQQg(2, 5, 6, 7, 1, 4, 3, eq(fl(ipos(1))), eq(fl(ipos(2))), MN, MI)
          msq = MN
       endif
    endif
    ! MCFM's four-quark normalisation (qqb_z2jet_g): 4 g^6 e^4 8 MN
    msq = avg*cnorm*32*msq
  end subroutine me41_tree
end module me41
