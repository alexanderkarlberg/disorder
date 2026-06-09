! Small odule to carry information about the EW state to pass to disent
module mod_ew_state
  use types, only: dp
  implicit none
  private

  public :: set_ew_state, set_ew_constants
  public :: ew_nc, ew_cc, ew_positron, ew_neutrino
  public :: ew_noz, ew_zonly, ew_intonly, ew_nc_mode
  public :: ew_alpha_em, ew_sin_thw_sq, ew_sin_2thw_sq
  public :: ew_ve, ew_ae, ew_ve2, ew_ae2, ew_ve2_ae2, ew_two_ve_ae
  public :: ew_mz, ew_mw, ew_gf

  logical, save :: ew_nc = .true.
  logical, save :: ew_cc = .false.
  logical, save :: ew_positron = .false.
  logical, save :: ew_neutrino = .false.
  logical, save :: ew_noz = .true.
  logical, save :: ew_zonly = .false.
  logical, save :: ew_intonly = .false.
  integer, save :: ew_nc_mode = 0

  real(dp), save :: ew_alpha_em = 0._dp
  real(dp), save :: ew_sin_thw_sq = 0._dp
  real(dp), save :: ew_sin_2thw_sq = 0._dp
  real(dp), save :: ew_ve = 0._dp
  real(dp), save :: ew_ae = 0._dp
  real(dp), save :: ew_ve2 = 0._dp
  real(dp), save :: ew_ae2 = 0._dp
  real(dp), save :: ew_ve2_ae2 = 0._dp
  real(dp), save :: ew_two_ve_ae = 0._dp
  real(dp), save :: ew_mz = 0._dp
  real(dp), save :: ew_mw = 0._dp
  real(dp), save :: ew_gf = 0._dp

contains

  subroutine set_ew_state(nc, cc, positron, neutrino, noZ, Zonly, intonly)
    logical, intent(in) :: nc, cc, positron, neutrino, noZ, Zonly, intonly
    ew_nc = nc
    ew_cc = cc
    ew_positron = positron
    ew_neutrino = neutrino
    ew_noz = noZ
    ew_zonly = Zonly
    ew_intonly = intonly

    if (.not. nc) then
      ew_nc_mode = -1
    elseif (noZ) then
      ew_nc_mode = 1
    elseif (Zonly) then
      ew_nc_mode = 3
    elseif (intonly) then
      ew_nc_mode = 2
    else
      ew_nc_mode = 0
    endif
  end subroutine set_ew_state

  subroutine set_ew_constants(alpha_em, sin_thw_sq, sin_2thw_sq, ve, ae, ve2, ae2, ve2_ae2, two_ve_ae, mz, mw, gf)
    real(dp), intent(in) :: alpha_em, sin_thw_sq, sin_2thw_sq
    real(dp), intent(in) :: ve, ae, ve2, ae2, ve2_ae2, two_ve_ae
    real(dp), intent(in) :: mz, mw, gf
    ew_alpha_em = alpha_em
    ew_sin_thw_sq = sin_thw_sq
    ew_sin_2thw_sq = sin_2thw_sq
    ew_ve = ve
    ew_ae = ae
    ew_ve2 = ve2
    ew_ae2 = ae2
    ew_ve2_ae2 = ve2_ae2
    ew_two_ve_ae = two_ve_ae
    ew_mz = mz
    ew_mw = mw
    ew_gf = gf
  end subroutine set_ew_constants

end module mod_ew_state
