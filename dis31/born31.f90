!-----------------------------------------------------------------------
! Colour- and spin-correlated DIS 3+1 Born matrix elements (photon
! exchange) for Catani-Seymour dipoles, from MCFM 10.3's colour-ordered
! amplitudes crossed to DIS (as me31; dis31/mcfm/README.md).
!
! DIS layout and normalisation as me31 (P(4,7): 1 incoming parton, 2-4
! outgoing partons, 5 q, 6, 7 leptons; fl(1:4); averaged over the incoming
! lepton and parton, alpha = 1/137, divided by (alpha_s/2pi)^2).
!
! born31_cc(P, fl, msq, cc): msq = |M|^2 (= me31), cc(i,k) = <M|T_i.T_k|M>
!   for i /= k (cc(i,i) = 0), colour charges with all partons outgoing
!   (Catani-Seymour): sum_{k /= i} cc(i,k) = -C_i msq.
! born31_sc(P, fl, ig, n, msqv, ccv): the gluon in slot ig contracted with
!   the real vector n (n.p_ig = 0) instead of summed over polarisations:
!   msqv = M_mu M_nu^* n^mu n^nu summed over colours, ccv(i,k) the same
!   with T_i.T_k. For n = e1, e2 (transverse, orthonormal) msqv(e1) +
!   msqv(e2) = msq.
!
! Colour bases and matrices <c_m|T_i.T_k|c_n> from explicit SU(3) (dis31/
! tests/colour.py): q qbar g g: c1 = (T^A T^B)_{q qbar}, c2 = (T^B T^A);
! four quarks: D = T^a_{q1 qb2} T^a_{q3 qb4}, E = T^a_{q3 qb2} T^a_{q1 qb4}
! (all partons outgoing: an incoming quark is an outgoing antiquark).
!-----------------------------------------------------------------------
module born31
  implicit none
  private
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: mxpart = 14
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: xn = 3, V = xn**2 - 1
  real(dp), parameter :: cnorm = (4*pi/137.0_dp)**2*64*pi**4
  real(dp), parameter :: eq(5) = [-1.0_dp/3, 2.0_dp/3, -1.0_dp/3, 2.0_dp/3, -1.0_dp/3]
  integer, parameter :: swp(2) = [2, 1]
  ! q qbar g g: roles 1 = q, 2 = qbar, 3 = A, 4 = B; cqg(:,:,r1,r2) for r1 < r2
  real(dp) :: cqg(2,2,4,4), gqg(2,2)
  ! four quarks: roles 1 = q1, 2 = qb2, 3 = q3, 4 = qb4
  real(dp) :: c4q(2,2,4,4), g4q(2,2)
  logical :: init = .false.
  ! debug/validation: exchange the assignment of the colour-ordered
  ! amplitudes to c1, c2 (q qbar g g)
  logical, public :: born31_swapc = .false.
  public :: born31_cc, born31_sc
contains

  subroutine setup()
    integer :: nwz
    common /nwz/ nwz
    nwz = 0
    cqg = 0; c4q = 0
    gqg = reshape([16.0_dp/3, -2.0_dp/3, -2.0_dp/3, 16.0_dp/3], [2,2])
    call put(cqg, 1, 2, -1.0_dp/9, -10.0_dp/9, -1.0_dp/9)
    call put(cqg, 1, 3, -8.0_dp, 1.0_dp, 1.0_dp)
    call put(cqg, 1, 4, 1.0_dp, 1.0_dp, -8.0_dp)
    call put(cqg, 2, 3, 1.0_dp, 1.0_dp, -8.0_dp)
    call put(cqg, 2, 4, -8.0_dp, 1.0_dp, 1.0_dp)
    call put(cqg, 3, 4, -9.0_dp, 0.0_dp, -9.0_dp)
    g4q = reshape([2.0_dp, -2.0_dp/3, -2.0_dp/3, 2.0_dp], [2,2])
    call put(c4q, 1, 2, 1.0_dp/3, -1.0_dp/9, -7.0_dp/3)
    call put(c4q, 1, 3, -2.0_dp/3, 10.0_dp/9, -2.0_dp/3)
    call put(c4q, 1, 4, -7.0_dp/3, -1.0_dp/9, 1.0_dp/3)
    call put(c4q, 2, 3, -7.0_dp/3, -1.0_dp/9, 1.0_dp/3)
    call put(c4q, 2, 4, -2.0_dp/3, 10.0_dp/9, -2.0_dp/3)
    call put(c4q, 3, 4, 1.0_dp/3, -1.0_dp/9, -7.0_dp/3)
    init = .true.
  end subroutine setup

  subroutine put(c, i, k, a11, a12, a22)
    real(dp), intent(inout) :: c(2,2,4,4)
    integer, intent(in) :: i, k
    real(dp), intent(in) :: a11, a12, a22
    c(:,:,i,k) = reshape([a11, a12, a12, a22], [2,2])
    c(:,:,k,i) = c(:,:,i,k)
  end subroutine put

  subroutine born31_cc(P, fl, msq, cc)
    real(dp), intent(in) :: P(4,7)
    integer, intent(in) :: fl(4)
    real(dp), intent(out) :: msq, cc(4,4)
    call born31_any(P, fl, 0, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], msq, cc)
  end subroutine born31_cc

  subroutine born31_sc(P, fl, ig, n, msqv, ccv)
    real(dp), intent(in) :: P(4,7), n(4)
    integer, intent(in) :: fl(4), ig
    real(dp), intent(out) :: msqv, ccv(4,4)
    if (fl(ig) /= 0) stop 'born31_sc: slot ig is not a gluon'
    call born31_any(P, fl, ig, n, msqv, ccv)
  end subroutine born31_sc

  ! ig = 0: polarisation sums; ig > 0: gluon ig contracted with n
  subroutine born31_any(P, fl, ig, n, msq, cc)
    real(dp), intent(in) :: P(4,7), n(4)
    integer, intent(in) :: fl(4), ig
    real(dp), intent(out) :: msq, cc(4,4)
    complex(dp) :: za(mxpart,mxpart), zb(mxpart,mxpart)
    common /zprods/ za, zb
    real(dp) :: pm(mxpart,4), avg, ch, q(4,4)
    integer :: ng, i, k, iq, iqb, igl(2), role(4), sg, nsame, ip(3), kqb, kq, f0, ms(4)
    if (.not. init) call setup()
    msq = 0; cc = 0
    ng = count(fl == 0)
    avg = merge(1.0_dp/(2*2*V), 1.0_dp/(2*2*xn), fl(1) == 0)
    pm = 0
    pm(3,:) = P(:,7); pm(4,:) = -P(:,6)
    if (ng == 2) then
       if (fl(1) /= 0) then
          ! q g g: MCFM i1 = outgoing (anti)quark k, i2 = incoming one, gluons
          ! A, B in MCFM slots 5, 6
          k = 0
          do i = 2, 4
             if (fl(i) == fl(1)) k = i
          enddo
          if (k == 0) return
          igl = pack([2, 3, 4], [2, 3, 4] /= k)
          pm(1,:) = -P(:,1); pm(2,:) = P(:,k); pm(5,:) = P(:,igl(1)); pm(6,:) = P(:,igl(2))
          ! colour roles of the slots (q, qbar, A, B) for the colour-ordered
          ! amplitudes A1 <-> c1 = (T^A T^B)_{q qbar}: q = MCFM slot 2 (= k),
          ! qbar = slot 1. An incoming antiquark has the same amplitudes
          ! (charge conjugation), with q <-> qbar and the colour order
          ! reversed, which leaves these roles unchanged.
          role = 0; role(k) = 1; role(1) = 2
          role(igl(1)) = 3; role(igl(2)) = 4
          ms = [2, 1, 5, 6]
          ch = eq(abs(fl(1)))**2
       else
          ! g -> q qbar g: MCFM i1 = q, i2 = qbar, gluons A = incoming (slot 1),
          ! B = outgoing
          iq = 0; iqb = 0
          do i = 2, 4
             if (fl(i) > 0) iq = i
             if (fl(i) < 0) iqb = i
             if (fl(i) == 0) igl(2) = i
          enddo
          if (iq == 0 .or. iqb == 0 .or. fl(iq) /= -fl(iqb)) return
          igl(1) = 1
          pm(1,:) = -P(:,1); pm(2,:) = P(:,iq); pm(5,:) = P(:,iqb); pm(6,:) = P(:,igl(2))
          role = 0; role(iq) = 1; role(iqb) = 2; role(1) = 3; role(igl(2)) = 4
          ms = [2, 5, 1, 6]
          ch = eq(fl(iq))**2
       endif
       call spinoru(6, pm, za, zb)
       call qqgg_forms(pm, ms, ig, igl, n, q)
       msq = avg*cnorm*ch*q(1,1)
       do i = 1, 4
          do k = 1, 4
             if (i /= k) cc(i,k) = avg*cnorm*ch*q(role(i), role(k))
          enddo
       enddo
       return
    endif
    if (ng /= 0 .or. fl(1) == 0 .or. ig /= 0) return
    ! four (anti)quarks
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
    if (fl(ip(1)) == f0 .and. fl(ip(2)) == f0) then
       if (fl(kqb) /= -f0) return
       k = ip(1); kq = ip(2)
       call fourq_forms(P, k, kq, kqb, eq(abs(f0)), eq(abs(f0)), .true., q)
    else
       if (fl(ip(1)) == f0) then
          k = ip(1); kq = ip(2)
       elseif (fl(ip(2)) == f0) then
          k = ip(2); kq = ip(1)
       else
          return
       endif
       if (fl(kqb) /= -fl(kq)) return
       call fourq_forms(P, k, kq, kqb, eq(abs(f0)), eq(abs(fl(kq))), .false., q)
    endif
    ! roles: quark case q1 = k, qb2 = 1, q3 = kq, qb4 = kqb; antiquark case
    ! (all partons outgoing) q1 = 1, qb2 = k, q3 = kqb, qb4 = kq
    if (sg > 0) then
       role(k) = 1; role(1) = 2; role(kq) = 3; role(kqb) = 4
    else
       role(1) = 1; role(k) = 2; role(kqb) = 3; role(kq) = 4
    endif
    msq = avg*cnorm*q(1,1)
    do i = 1, 4
       do k = 1, 4
          if (i /= k) cc(i,k) = avg*cnorm*q(role(i), role(k))
       enddo
    enddo
  end subroutine born31_any

  ! q qbar g g: quadratic colour forms sum_hel A^+ C^{rs} A for all role
  ! pairs r /= s, and the metric (q(1,1) = sum A^+ G A); ig > 0: gluon ig
  ! (DIS slot) contracted with n. pm, za, zb set; ms = the MCFM slots of
  ! (q, qbar, A, B); igl = the DIS slots of (A, B).
  subroutine qqgg_forms(pm, ms, ig, igl, n, q)
    real(dp), intent(in) :: pm(mxpart,4), n(4)
    integer, intent(in) :: ms(4), ig, igl(2)
    real(dp), intent(out) :: q(4,4)
    complex(dp) :: za(mxpart,mxpart), zb(mxpart,mxpart)
    common /zprods/ za, zb
    complex(dp) :: q1(-1:1,-1:1), q2(-1:1,-1:1), zab(mxpart,mxpart), zba(mxpart,mxpart)
    complex(dp) :: ab(2,2,2), ba(2,2,2), A(2,16)
    integer, parameter :: pol(2) = [-1, 1]
    integer :: nh, j, k, l, r, s, pg, pq, pl, ia, ib
    real(dp) :: nDp, fac
    nh = 0
    if (ig == 0) then
       ! subqcd: LL (leptons 3, 4) and LR (4, 3); the right-handed quark line
       ! by parity (factor 2)
       do l = 1, 2
          if (l == 1) then
             call subqcd(ms(1), ms(2), 3, 4, ms(3), ms(4), za, zb, q1)
             call subqcd(ms(1), ms(2), 3, 4, ms(4), ms(3), za, zb, q2)
          else
             call subqcd(ms(1), ms(2), 4, 3, ms(3), ms(4), za, zb, q1)
             call subqcd(ms(1), ms(2), 4, 3, ms(4), ms(3), za, zb, q2)
          endif
          do j = 1, 2
             do k = 1, 2
                nh = nh + 1
                A(1,nh) = q1(pol(j),pol(k)); A(2,nh) = q2(pol(k),pol(j))
             enddo
          enddo
       enddo
       fac = 2
    else
       ! subqcdn: its second gluon (ib) contracted with n, the first (ia)
       ! summed; qcdab corresponds to the colour ordering (ib, ia) of
       ! subqcd's convention, qcdba to (ia, ib) (fixed by the polarisation
       ! sum ccv(e1) + ccv(e2) = cc, harness_born31)
       if (ig == igl(2)) then
          ia = ms(3); ib = ms(4)
       else
          ia = ms(4); ib = ms(3)
       endif
       nDp = n(4)*pm(ia,4) - n(1)*pm(ia,1) - n(2)*pm(ia,2) - n(3)*pm(ia,3)
       call spinork(6, pm, zab, zba, n)
       call subqcdn(ms(1), ms(2), 3, 4, ia, ib, nDp, za, zb, zab, zba, ab, ba)
       do pg = 1, 2
          do pq = 1, 2
             do pl = 1, 2
                nh = nh + 1
                if (ig == igl(2)) then
                   A(1,nh) = ba(pg,pq,pl); A(2,nh) = ab(pg,pq,pl)
                else
                   A(1,nh) = ab(pg,pq,pl); A(2,nh) = ba(pg,pq,pl)
                endif
             enddo
          enddo
       enddo
       ! MCFM's contraction with n is normalised such that n = e1, e2
       ! (orthonormal, transverse) sum to half the polarisation sum
       fac = 2
    endif
    if (born31_swapc) then
       do j = 1, nh
          A(:,j) = A([2, 1],j)
       enddo
    endif
    ! off-diagonal: the role pairs; q(1,1): the metric (|M|^2)
    q = 0
    do r = 1, 4
       do s = r + 1, 4
          q(r,s) = fac*form(cqg(:,:,r,s)); q(s,r) = q(r,s)
       enddo
    enddo
    q(1,1) = fac*form(gqg)
  contains
    real(dp) function form(c)
      real(dp), intent(in) :: c(2,2)
      integer :: h
      form = 0
      do h = 1, nh
         form = form + c(1,1)*abs(A(1,h))**2 + c(2,2)*abs(A(2,h))**2 &
              & + 2*c(1,2)*real(A(1,h)*conjg(A(2,h)), dp)
      enddo
    end function form
  end subroutine qqgg_forms

  ! four quarks: the forms for the direct (D) and, for identical quarks,
  ! exchange (E) colour structures, sum_hel; q(1,1) = metric
  subroutine fourq_forms(P, k, kq, kqb, cq, cQ2, ident, q)
    real(dp), intent(in) :: P(4,7), cq, cQ2
    integer, intent(in) :: k, kq, kqb
    logical, intent(in) :: ident
    real(dp), intent(out) :: q(4,4)
    complex(dp) :: za(mxpart,mxpart), zb(mxpart,mxpart)
    common /zprods/ za, zb
    real(dp) :: pm(mxpart,4), sDD, sEE, sDE
    complex(dp) :: A(2,2,2), B(2,2,2), Ae(2,2,2), Be(2,2,2)
    integer :: j1, j3, r, s
    pm = 0
    pm(1,:) = -P(:,1); pm(2,:) = P(:,k); pm(3,:) = P(:,7); pm(4,:) = -P(:,6)
    pm(5,:) = P(:,kq); pm(6,:) = P(:,kqb)
    call spinoru(6, pm, za, zb)
    call ampqqb_qqb(2, 1, 5, 6, A, B)
    if (.not. ident) then
       ! e_q A - e_Q B (me31), one colour structure D
       sDD = sum(abs(cq*A - cQ2*B)**2)
       sEE = 0; sDE = 0
    else
       call ampqqb_qqb(5, 1, 2, 6, Ae, Be)
       A = cq*(A - B); Ae = cq*(Ae - Be)
       sDD = sum(abs(A)**2); sEE = sum(abs(Ae)**2)
       ! E' = -E in the colour basis (me31's interference sign)
       sDE = 0
       do j1 = 1, 2; do j3 = 1, 2
          sDE = sDE - real(A(j1,swp(j1),j3)*conjg(Ae(j1,swp(j1),j3)), dp)
       enddo; enddo
    endif
    ! |M|^2 = 16 [C11 sDD + C22 sEE + 2 C12 sDE] (4V = 16 G11)
    q = 0
    do r = 1, 4
       do s = r + 1, 4
          q(r,s) = 16*(c4q(1,1,r,s)*sDD + c4q(2,2,r,s)*sEE + 2*c4q(1,2,r,s)*sDE)
          q(s,r) = q(r,s)
       enddo
    enddo
    q(1,1) = 16*(g4q(1,1)*sDD + g4q(2,2)*sEE + 2*g4q(1,2)*sDE)
  end subroutine fourq_forms
end module born31
