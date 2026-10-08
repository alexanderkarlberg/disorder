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
! ew31_mode: 0 photon, 1 photon + Z. ew31_lepton: 0 e-, 1 e+ (hl
! exchanged), 2 neutrino, 3 antineutrino (one helicity: couplings times
! sqrt(2) against the lepton spin average 1/2 of the matrix elements).
!-----------------------------------------------------------------------
module ew31
  implicit none
  private
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: mz = 91.1876_dp, mw = 80.398_dp
  real(dp), parameter :: xw = 1 - (mw/mz)**2, sin2w = 2*sqrt(xw*(1 - xw))
  real(dp), parameter :: eq(5) = [-1.0_dp/3, 2.0_dp/3, -1.0_dp/3, 2.0_dp/3, -1.0_dp/3]
  real(dp), parameter :: tau(5) = [-1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp]
  integer, public :: ew31_mode = 0, ew31_lepton = 0
  public :: ew31_cpl
contains

  subroutine ew31_cpl(f, Q2, c)
    integer, intent(in) :: f
    real(dp), intent(in) :: Q2
    real(dp), intent(out) :: c(2,2)
    real(dp) :: zq(2), zl(2), ql, prop
    integer :: a, hq, hl
    a = abs(f)
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
end module ew31
