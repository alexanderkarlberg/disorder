! Small odule to carry information about the EW state to pass to disent
module mod_ew_state
  implicit none
  private

  public :: set_ew_state
  public :: ew_nc, ew_cc, ew_positron, ew_neutrino
  public :: ew_noz, ew_zonly, ew_intonly, ew_nc_mode

  logical, save :: ew_nc      = .true.
  logical, save :: ew_cc      = .false.
  logical, save :: ew_positron = .false.
  logical, save :: ew_neutrino  = .false.
  logical, save :: ew_noz     = .true.
  logical, save :: ew_zonly   = .false.
  logical, save :: ew_intonly = .false.
  integer, save :: ew_nc_mode = 0
  ! 0 = full NC
  ! 1 = photon only
  ! 2 = gamma/Z interference only
  ! 3 = Z only

contains

  subroutine set_ew_state(nc, cc, positron, neutrino, noZ, Zonly, intonly)
    logical, intent(in) :: nc, cc, positron, neutrino, noZ, Zonly, intonly

    ew_nc       = nc
    ew_cc       = cc
    ew_positron = positron
    ew_neutrino = neutrino
    ew_noz      = noZ
    ew_zonly    = Zonly
    ew_intonly  = intonly

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

end module mod_ew_state
