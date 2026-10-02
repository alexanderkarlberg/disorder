!-----------------------------------------------------------------------
! DIS 3+1 one-loop matrix elements (photon exchange), from the one-loop
! amplitudes of Bern, Dixon, Kosower (hep-ph/9708239: q qbar g g + V;
! hep-ph/9610370: q qbar Q Qbar + V) as implemented in MCFM 10.3
! (dis31/mcfm/loop), crossed to DIS as me31.
!
! virt31_ren(P, fl, mu2, v, t): the one-loop interference 2 Re <M0|M1>,
!   UV-renormalised in MS-bar at the scale mu2 (= mu_R^2), 't Hooft-Veltman
!   scheme, with (4 pi)^eps/Gamma(1 - eps) factored out (the Catani-Seymour
!   normalisation of the I operator), as Laurent coefficients v(-2:0) (1/eps^2,
!   1/eps, finite), in units of me31 times alpha_s/2pi; t = the tree (= me31).
!   n_f = 5 massless flavours, no top loops; the photon couples to no closed
!   quark loop (Furry).
!   Its poles are those of -<I> (checked, tests/harness_virt31), and the
!   finite part agrees with NNLOJET v1.0.2 for every channel up to the
!   known constant (pi^2/12) sum_i C_i t of NNLOJET's normalisation
!   (dis31_nnlojet/harness_v31*, outside the repository).
! virt31_qqgg, virt31_4q: MCFM's raw results (BDK conventions: FDH, overall
!   c_Gamma = (4 pi)^eps/Gamma(1 - eps) + O(eps^3)); q qbar g g unrenormalised
!   with DRED coupling, four quarks UV-renormalised in a6routine (DRED
!   coupling converted to MS-bar). virt31_ren converts:
!     FDH -> HV: - sum_i gamma~_i t, gamma~_q = CF/2, gamma~_g = CA/6 (the
!       standard regularisation-scheme constants; confirmed numerically
!       against NNLOJET for both channel types, 2 Oct 2026);
!     q qbar g g: UV counterterm - 2 beta0/eps t and the DRED -> MS-bar
!       coupling for alpha_s^2, + N/3 t (as MCFM's qqb_z2jet_v, subuv).
!   Net: q qbar g g: v(-1) - 2 beta0 t, v(0) - CF t; four quarks: v(0) - 2 CF t.
!-----------------------------------------------------------------------
module virt31
  implicit none
  private
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: mxpart = 14
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: xn = 3, nadj = xn**2 - 1
  real(dp), parameter :: cnorm = (4*pi/137.0_dp)**2*64*pi**4
  real(dp), parameter :: avg4 = 1.0_dp/(2*2*xn)
  real(dp), parameter :: eq(5) = [-1.0_dp/3, 2.0_dp/3, -1.0_dp/3, 2.0_dp/3, -1.0_dp/3]
  public :: virt31_ren, virt31_qqgg, virt31_4q
contains

  subroutine virt31_ren(P, fl, mu2, v, t)
    real(dp), intent(in) :: P(4,7), mu2
    integer, intent(in) :: fl(4)
    real(dp), intent(out) :: v(-2:0), t
    real(dp), parameter :: CF = (xn**2 - 1)/(2*xn), b0 = (11*xn - 2*5)/6
    if (count(fl == 0) == 2) then
       call virt31_qqgg(P, fl, mu2, v, t)
       v(-1) = v(-1) - 2*b0*t
       v(0) = v(0) - CF*t
    elseif (count(fl == 0) == 0) then
       call virt31_4q(P, fl, mu2, v, t)
       v(0) = v(0) - 2*CF*t
    else
       v = 0; t = 0
    endif
  end subroutine virt31_ren

  subroutine virt31_qqgg(P, fl, mu2, v, t)
    real(dp), intent(in) :: P(4,7), mu2
    integer, intent(in) :: fl(4)
    real(dp), intent(out) :: v(-2:0), t
    complex(dp) :: za(mxpart,mxpart), zb(mxpart,mxpart)
    common /zprods/ za, zb
    real(dp) :: epinv, epinv2, scale, musq
    common /epinv/ epinv
    common /epinv2/ epinv2
    common /mcfmscale/ scale, musq
    integer :: toploops
    logical :: toplight, topvector, topaxial, onlyaxial
    common /toploops/ toploops, toplight, topvector, topaxial, onlyaxial
    real(dp) :: pm(mxpart,4), avg, ch, w(3,2), tr
    integer :: k, i, iq, iqb, igl(2), ie
    real(dp), parameter :: evals(3) = [0.0_dp, 1.0_dp, -1.0_dp]
    v = 0; t = 0
    if (count(fl == 0) /= 2) return
    ! MCFM labels (0 -> q(1) g(2) g(3) qbar(4) l(5) a(6), all outgoing)
    pm = 0
    pm(5,:) = P(:,7); pm(6,:) = -P(:,6)
    if (fl(1) /= 0) then
       k = 0
       do i = 2, 4
          if (fl(i) == fl(1)) k = i
       enddo
       if (k == 0) return
       igl = pack([2, 3, 4], [2, 3, 4] /= k)
       pm(1,:) = P(:,k); pm(4,:) = -P(:,1); pm(2,:) = P(:,igl(1)); pm(3,:) = P(:,igl(2))
       ch = eq(abs(fl(1)))**2
       avg = 1.0_dp/(2*2*xn)
    else
       iq = 0; iqb = 0
       do i = 2, 4
          if (fl(i) > 0) iq = i
          if (fl(i) < 0) iqb = i
          if (fl(i) == 0) igl(2) = i
       enddo
       if (iq == 0 .or. iqb == 0 .or. fl(iq) /= -fl(iqb)) return
       pm(1,:) = P(:,iq); pm(4,:) = P(:,iqb); pm(2,:) = -P(:,1); pm(3,:) = P(:,igl(2))
       ch = eq(fl(iq))**2
       avg = 1.0_dp/(2*2*nadj)
    endif
    call spinoru(6, pm, za, zb)
    toploops = 1; toplight = .false.; topvector = .false.; topaxial = .false.; onlyaxial = .false.
    musq = mu2; scale = sqrt(mu2)
    ! MCFM's poles: epinv = epinv2 = 1/eps (the double pole is
    ! epinv*epinv2), so the result is a quadratic in e = 1/eps: evaluate at
    ! e = 0, 1, -1
    do ie = 1, 3
       epinv = evals(ie); epinv2 = epinv
       call qqgg_photon(w(ie,:), tr)
    enddo
    ! w(ie,1) = interference (xzqqgg_v's colour combination, couplings
    ! stripped), w(ie,2) unused; tr = the tree, me31 = 96 avg cnorm e_q^2 tr;
    ! the virtual relative to the tree: N (xzqqgg_v: fac = xn ason2pi
    ! times xzqqgg's)
    t = 96*avg*cnorm*ch*tr
    v(0) = 96*xn*avg*cnorm*ch*w(1,1)
    v(-1) = 96*xn*avg*cnorm*ch*(w(2,1) - w(3,1))/2
    v(-2) = 96*xn*avg*cnorm*ch*((w(2,1) + w(3,1))/2 - w(1,1))
  end subroutine virt31_qqgg

  ! xzqqgg_v (MCFM src/Zbb) for colourchoice = 0, photon exchange, no
  ! boson-on-loop (Furry) or top-loop pieces; w(1) = sum_hel of the
  ! interference with the colour weights of xzqqgg_v, w(2) = the tree
  ! sum_hel [|m1|^2 + |m2|^2 - |m1 + m2|^2/N^2] times N^2 (V N/4 per
  ! colour... normalised below), both without couplings
  subroutine qqgg_photon(w, tr)
    real(dp), intent(out) :: w(2), tr
    complex(dp) :: za(mxpart,mxpart), zb(mxpart,mxpart)
    common /zprods/ za, zb
    include 'heldefs.f'
    complex(dp) :: m(2), ml1(2), ml2(2), ml3, ml4(2)
    complex(dp), external :: a6treeg1, a61g1lc, a61g1slc, a61g1nf, a63g1
    integer :: j, lh, h2, h3, hq, h(2:3)
    integer, parameter :: i1(2) = [1, 4], i2(2) = [2, 3], i3(2) = [3, 2], &
         & i4(2) = [4, 1], i5(2) = [6, 5], i6(2) = [5, 6]
    integer, parameter :: st1(2,2) = reshape([hqpgmgmqbm, hqpgmgpqbm, hqpgpgmqbm, hqpgpgpqbm], [2,2])
    integer, parameter :: st2(2,2) = reshape([hqpqbmgpgp, hqpqbmgpgm, hqpqbmgmgp, hqpqbmgmgm], [2,2])
    integer, parameter :: st3(2,2) = reshape([hqpqbmgmgm, hqpqbmgmgp, hqpqbmgpgm, hqpqbmgpgp], [2,2])
    real(dp), parameter :: xnsq = xn**2
    w = 0; tr = 0
    do hq = 1, 2
       do lh = 1, 2
          do h2 = 1, 2
             do h3 = 1, 2
                h(2) = h2; h(3) = h3
                do j = 1, 2
                   if (hq == 1) then
                      m(j) = a6treeg1(st1(3-h(i2(j)),3-h(i3(j))), i1(1), i2(j), i3(j), i4(1), i6(lh), i5(lh), zb, za)
                      ml1(j) = a61g1lc(st1(3-h(i2(j)),3-h(i3(j))), i1(1), i2(j), i3(j), i4(1), i6(lh), i5(lh), zb, za)
                      ml2(j) = a61g1slc(st2(3-h(i2(j)),3-h(i3(j))), i1(1), i2(j), i3(j), i4(1), i6(lh), i5(lh), zb, za)
                      ml4(j) = a61g1nf(st1(3-h(i2(j)),3-h(i3(j))), i1(1), i2(j), i3(j), i4(1), i6(lh), i5(lh), zb, za)
                   else
                      m(j) = a6treeg1(st1(h(i2(j)),h(i3(j))), i1(1), i2(j), i3(j), i4(1), i5(lh), i6(lh), za, zb)
                      ml1(j) = a61g1lc(st1(h(i2(j)),h(i3(j))), i1(1), i2(j), i3(j), i4(1), i5(lh), i6(lh), za, zb)
                      ml2(j) = a61g1slc(st2(h(i2(j)),h(i3(j))), i1(1), i2(j), i3(j), i4(1), i5(lh), i6(lh), za, zb)
                      ml4(j) = a61g1nf(st1(h(i2(j)),h(i3(j))), i1(1), i2(j), i3(j), i4(1), i5(lh), i6(lh), za, zb)
                   endif
                enddo
                if (hq == 1) then
                   ml3 = a63g1(st3(3-h2,3-h3), 1, 4, 2, 3, i6(lh), i5(lh), zb, za)
                else
                   ml3 = a63g1(st3(h2,h3), 1, 4, 2, 3, i5(lh), i6(lh), za, zb)
                endif
                w(1) = w(1) + real(conjg(m(1))*(ml1(1) - (ml1(1) + ml2(1) + ml1(2) - ml3)/xnsq &
                     & + (ml2(1) + ml2(2))/xnsq**2 + ml4(1)/xn - (ml4(1) + ml4(2))/xn**3), dp) &
                     & + real(conjg(m(2))*(ml1(2) - (ml1(2) + ml2(2) + ml1(1) - ml3)/xnsq &
                     & + (ml2(1) + ml2(2))/xnsq**2 + ml4(2)/xn - (ml4(1) + ml4(2))/xn**3), dp)
                tr = tr + abs(m(1))**2 + abs(m(2))**2 - abs(m(1) + m(2))**2/xnsq
             enddo
          enddo
       enddo
    enddo
  end subroutine qqgg_photon

  ! four (anti)quarks: incoming q (or qbar by charge conjugation), outgoing
  ! q of the same line, pair Q Qbar (Q = q: identical quarks). MCFM's q q ->
  ! q q channel of qqb_z2jet_v crossed: MCFM slot 1 = -incoming quark, 5 =
  ! outgoing quark of that line, 2 = the pair's antiquark, 6 = its quark;
  ! identical quarks: 5 <-> 6 exchanged (the two outgoing quarks). The
  ! a63z (vector loop) terms vanish for the photon.
  ! v(-2:0): Laurent coefficients of MCFM's (DRED, UV-renormalised in
  ! a6routine) interference, t = the tree (= me31), units as virt31_qqgg.
  subroutine virt31_4q(P, fl, mu2, v, t)
    real(dp), intent(in) :: P(4,7), mu2
    integer, intent(in) :: fl(4)
    real(dp), intent(out) :: v(-2:0), t
    complex(dp) :: za(mxpart,mxpart), zb(mxpart,mxpart)
    common /zprods/ za, zb
    real(dp) :: epinv, epinv2, scale, musq
    common /epinv/ epinv
    common /epinv2/ epinv2
    common /mcfmscale/ scale, musq
    integer :: toploops
    logical :: toplight, topvector, topaxial, onlyaxial
    common /toploops/ toploops, toplight, topvector, topaxial, onlyaxial
    real(dp) :: pm(mxpart,4), w(3), tr, cq, cQ2
    integer :: sg
    integer :: i, k, kq, kqb, nsame, ip(3), f0, ie
    logical :: ident
    real(dp), parameter :: evals(3) = [0.0_dp, 1.0_dp, -1.0_dp]
    v = 0; t = 0
    if (count(fl == 0) /= 0 .or. fl(1) == 0) return
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
    cq = eq(abs(f0))
    ident = fl(ip(1)) == f0 .and. fl(ip(2)) == f0
    if (ident) then
       if (fl(kqb) /= -f0) return
       k = ip(1); kq = ip(2); cQ2 = cq
    else
       if (fl(ip(1)) == f0) then
          k = ip(1); kq = ip(2)
       elseif (fl(ip(2)) == f0) then
          k = ip(2); kq = ip(1)
       else
          return
       endif
       if (fl(kqb) /= -fl(kq)) return
       cQ2 = eq(abs(fl(kq)))
    endif
    pm = 0
    pm(1,:) = -P(:,1); pm(5,:) = P(:,k); pm(2,:) = P(:,kqb); pm(6,:) = P(:,kq)
    pm(3,:) = P(:,7); pm(4,:) = -P(:,6)
    call spinoru(6, pm, za, zb)
    toploops = 1; toplight = .false.; topvector = .false.; topaxial = .false.; onlyaxial = .false.
    musq = mu2; scale = sqrt(mu2)
    do ie = 1, 3
       epinv = evals(ie); epinv2 = epinv
       call fourq_photon(cq, cQ2, ident, w(ie), tr)
    enddo
    t = avg4*cnorm*tr
    v(0) = avg4*cnorm*w(1)
    v(-1) = avg4*cnorm*(w(2) - w(3))/2
    v(-2) = avg4*cnorm*((w(2) + w(3))/2 - w(1))
  end subroutine virt31_4q

  ! qqb_z2jet_v's q q branch for the photon (couplings Q(j) q1 -> charges;
  ! the lepton charge is an overall sign): w = 2 Re(tree^* loop) summed,
  ! tr = the tree, both times MCFM's faclo/(couplings) = 4 V and xn/2 for
  ! the loop (fac = faclo xn/2 ason2pi)
  subroutine fourq_photon(cq, cQ2, ident, w, tr)
    real(dp), intent(in) :: cq, cQ2
    logical, intent(in) :: ident
    real(dp), intent(out) :: w, tr
    complex(dp) :: za(mxpart,mxpart), zb(mxpart,mxpart)
    common /zprods/ za, zb
    complex(dp), external :: atreez, a61z, a62z
    complex(dp) :: tamp, lamp, tamps, lamps, lampx, lampsx
    integer :: polq, polb, polz
    w = 0; tr = 0
    do polq = 1, 2
       do polz = 1, 2
          do polb = 1, 2
             tamp = atreez(polq,polb,polz,5,2,6,1,4,3,za,zb)*cq &
                  & - atreez(3-polb,3-polq,polz,2,5,1,6,4,3,za,zb)*cQ2
             lamp = a61z(polq,polb,polz,5,2,6,1,4,3,za,zb)*cq &
                  & - a61z(3-polb,3-polq,polz,2,5,1,6,4,3,za,zb)*cQ2
             tr = tr + abs(tamp)**2
             w = w + xn/2*2*real(tamp*conjg(lamp), dp)
             if (ident) then
                tamps = -(atreez(polq,polb,polz,6,2,5,1,4,3,za,zb)*cq &
                     & - atreez(3-polb,3-polq,polz,2,6,1,5,4,3,za,zb)*cQ2)
                lamps = -(a61z(polq,polb,polz,6,2,5,1,4,3,za,zb)*cq &
                     & - a61z(3-polb,3-polq,polz,2,6,1,5,4,3,za,zb)*cQ2)
                lampx = -(a62z(polq,polb,polz,6,2,5,1,4,3,za,zb)/xn*cq &
                     & - a62z(3-polb,3-polq,polz,2,6,1,5,4,3,za,zb)/xn*cQ2)
                lampsx = a62z(polq,polb,polz,5,2,6,1,4,3,za,zb)/xn*cq &
                     & - a62z(3-polb,3-polq,polz,2,5,1,6,4,3,za,zb)/xn*cQ2
                tr = tr + abs(tamps)**2
                w = w + xn/2*2*real(tamps*conjg(lamps), dp)
                if (polq == polb) then
                   tr = tr - 2/xn*real(tamp*conjg(tamps), dp)
                   w = w + xn/2*2*(real(tamp*conjg(lampx), dp) + real(tamps*conjg(lampsx), dp))
                endif
             endif
          enddo
       enddo
    enddo
    tr = 4*nadj*tr; w = 4*nadj*w
  end subroutine fourq_photon
end module virt31
