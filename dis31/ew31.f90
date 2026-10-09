!-----------------------------------------------------------------------
! Electroweak couplings of the DIS 3+1 / 4+1 matrix elements (8 Oct
! 2026): photon (default) or photon + Z exchange, for helicity-resolved
! MCFM squares m2(hq, hl) (hq: helicity of the quark line, hl: of the
! lepton line; 1 = L, 2 = R, as MCFM's z2jetsq).
!
! ew31_cpl(f, Q2, c): c(hq, hl) = Q_f q_l + z_f(hq) z_l(hl) Q^2/(Q^2 + MZ^2)
!   for a quark line of flavour f (1-5 d u s c b; negative: antiquark) at
!   the spacelike q^2 = -Q^2. The DIS crossing of me31/born31 puts the
!   outgoing quark of a line into MCFM's incoming-quark slot, so the
!   helicity label hq is that of an antiquark: z_f exchanged for quarks
!   (fixed against disorder's NC MATFOR, harness_me31). As MCFM's
!   qqb_z2jet with prop = s34/(s34 - MZ^2), no width; Z couplings as MCFM's zcouple,
!   MZ, MW and sin^2(theta_W) = 1 - MW^2/MZ^2 as disorder's defaults.
!   For the photon (ew31_mode = 0) c = Q_f q_l for all helicities.
!
! ew31_mode: 0 photon, 1 photon + Z, 2 W (CC, 9 Oct). ew31_lepton: 0 e-, 1
! e+ (hl exchanged), 2 neutrino, 3 antineutrino (one helicity: couplings
! times sqrt(2) against the lepton spin average 1/2 of the matrix elements).
!
! W exchange (ew31_mode = 2): c(L,L) = Q^2/(Q^2 + MW^2)/(2 sin^2 theta_W)
! (left-handed quark and lepton lines, labels as for Z), W- for e- and
! nubar, W+ for e+ and nu; unit CKM within the complete generations (u,d),
! (c,s), no coupling for b (as disorder's CC and HOPPET). The line changes
! flavour: ew31_out(f) is the outgoing flavour of the boson line for an
! incoming parton f (f itself for photon/Z; W-: u -> d, c -> s, dbar ->
! ubar, sbar -> cbar; W+ the reverse), 0 if the W does not couple to f.
!
! Basis couplings for the flavour sums (nlo31, ew31_mode = 1): ew31_basis
! = 1 (2) gives c = 1 on the diagonal (off-diagonal) helicity pairs for the
! line(s) of flavour |f| = ew31_bf (exchanged for quarks as above, so that
! dipoles with the line crossed keep the relative exchange) and c = 0 for
! all other lines. With option 1 every channel is linear in the squared
! couplings of each line, and at tree level m2(1,1) = m2(2,2), m2(1,2) =
! m2(2,1) (parity), so |M|^2 = [E1 (c11^2 + c22^2) + E2 (c12^2 + c21^2)]/2
! summed over the lines, E1, E2 the values with the basis couplings on the
! line of a quark f > 0 and c = ew31_cpl(-f): ew31_w(f, Q2, w) gives
! these weights w = (c11^2 + c22^2, c12^2 + c21^2)/2 (photon: e_f^2 both).
!
! The boson on a closed quark loop (virt31, vector coupling): ew31_cv(Q2,
! cv), cv(l) = sum_f (c_f(1,l) + c_f(2,l))/2 over the five flavours (photon:
! q_l sum e_f); 1 in basis mode (the term is then c(h,l) X(h,l)).
! ew31_wl(f, Q2, w): its weights for a line of flavour f, as ew31_w with
! c^2 -> c cv, w = (c11 cv1 + c22 cv2, c21 cv1 + c12 cv2)/2 (basis off).
!-----------------------------------------------------------------------
module ew31
  implicit none
  private
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: mz = 91.1876_dp, mw = 80.398_dp
  real(dp), parameter :: xw = 1 - (mw/mz)**2, sin2w = 2*sqrt(xw*(1 - xw))
  real(dp), parameter :: eq(5) = [-1.0_dp/3, 2.0_dp/3, -1.0_dp/3, 2.0_dp/3, -1.0_dp/3]
  real(dp), parameter :: tau(5) = [-1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp]
  integer, public :: ew31_mode = 0, ew31_lepton = 0, ew31_basis = 0, ew31_bf = 0
  ! W exchange: drop the interference of the two W-on-the-pair assignments
  ! of identical quarks (me31, me41), as disorder's MATFOR (harness checks)
  logical, public :: ew31_ccmatfor = .false.
  public :: ew31_cpl, ew31_w, ew31_cv, ew31_wl, ew31_out
contains

  subroutine ew31_cpl(f, Q2, c)
    integer, intent(in) :: f
    real(dp), intent(in) :: Q2
    real(dp), intent(out) :: c(2,2)
    real(dp) :: zq(2), zl(2), ql, prop
    integer :: a, hq, hl
    a = abs(f)
    if (ew31_basis > 0) then
       c = 0
       if (a /= ew31_bf) return
       if (ew31_basis == 1) then
          c(1,1) = 1; c(2,2) = 1
       else
          c(1,2) = 1; c(2,1) = 1
       endif
       if (f > 0) c = c([2, 1],:)
       return
    endif
    if (ew31_mode == 2) then
       ! W: left-handed lines only, labels as for Z
       c = 0
       if (a > 4) return
       zq = [1.0_dp, 0.0_dp]; zl = [1.0_dp, 0.0_dp]
       if (ew31_lepton >= 2) zl = zl*sqrt(2.0_dp)
       if (mod(ew31_lepton, 2) == 1) zl = zl([2, 1])
       if (f > 0) zq = zq([2, 1])
       prop = Q2/(Q2 + mw**2)/(2*xw)
       do hq = 1, 2
          do hl = 1, 2
             c(hq,hl) = zq(hq)*zl(hl)*prop
          enddo
       enddo
       return
    endif
    if (ew31_lepton <= 1) then
       ql = -1; zl = [(-1 + 2*xw)/sin2w, 2*xw/sin2w]
    else
       ql = 0; zl = [sqrt(2.0_dp)/sin2w, 0.0_dp]
    endif
    if (mod(ew31_lepton, 2) == 1) zl = zl([2, 1])
    zq = [(tau(a) - 2*eq(a)*xw)/sin2w, -2*eq(a)*xw/sin2w]
    if (f > 0) zq = zq([2, 1])
    prop = merge(Q2/(Q2 + mz**2), 0.0_dp, ew31_mode == 1)
    do hq = 1, 2
       do hl = 1, 2
          c(hq,hl) = eq(a)*ql + zq(hq)*zl(hl)*prop
       enddo
    enddo
  end subroutine ew31_cpl

  subroutine ew31_w(f, Q2, w)
    integer, intent(in) :: f
    real(dp), intent(in) :: Q2
    real(dp), intent(out) :: w(2)
    real(dp) :: c(2,2)
    call ew31_cpl(-f, Q2, c)
    w = [c(1,1)**2 + c(2,2)**2, c(1,2)**2 + c(2,1)**2]/2
  end subroutine ew31_w

  subroutine ew31_cv(Q2, cv)
    real(dp), intent(in) :: Q2
    real(dp), intent(out) :: cv(2)
    real(dp) :: c(2,2)
    integer :: f
    cv = 1
    if (ew31_basis > 0) return
    cv = 0
    ! W: no closed loop (a single W vertex changes the loop's flavour)
    if (ew31_mode == 2) return
    do f = 1, 5
       call ew31_cpl(f, Q2, c)
       cv = cv + (c(1,:) + c(2,:))/2
    enddo
  end subroutine ew31_cv

  subroutine ew31_wl(f, Q2, w)
    integer, intent(in) :: f
    real(dp), intent(in) :: Q2
    real(dp), intent(out) :: w(2)
    real(dp) :: c(2,2), cv(2)
    call ew31_cv(Q2, cv)
    call ew31_cpl(-f, Q2, c)
    w = [c(1,1)*cv(1) + c(2,2)*cv(2), c(2,1)*cv(1) + c(1,2)*cv(2)]/2
  end subroutine ew31_wl

  ! the outgoing flavour of the boson line for an incoming parton f (0: no
  ! coupling); f for photon/Z
  integer function ew31_out(f)
    integer, intent(in) :: f
    logical :: wm
    ew31_out = f
    if (ew31_mode /= 2) return
    ew31_out = 0
    if (f == 0 .or. abs(f) > 4) return
    wm = ew31_lepton == 0 .or. ew31_lepton == 3
    if (wm) then
       ! W- absorbed: up-type quark -> down-type, down-type antiquark -> up-type
       if (f > 0 .and. mod(f, 2) == 0) ew31_out = f - 1
       if (f < 0 .and. mod(-f, 2) == 1) ew31_out = f - 1
    else
       if (f > 0 .and. mod(f, 2) == 1) ew31_out = f + 1
       if (f < 0 .and. mod(-f, 2) == 0) ew31_out = f + 1
    endif
  end function ew31_out
end module ew31
