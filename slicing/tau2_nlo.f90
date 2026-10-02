!----------------------------------------------------------------------
! NLO DIS 2+1 (O(alpha_s^2) coefficient of 2+1-jet observables) with
! 2-jettiness slicing, against DISENT's dipole-subtracted NLO, in one
! DISENT run (photon exchange, fixed x and Q^2, mu_R = mu_F = Q).
!
! For every DISENT event:
!  - reference: all O(alpha_s^2) contributions (VIRTHR, COLFOR, MATFOR,
!    SUBFOR counter-events) -> the full NLO coefficient;
!  - above the cut: the real four-parton events (MATFOR) with
!    tau_2 = T_2/Q > tau_cut. Counter-events, collinear terms and the
!    three-parton virtual have Born kinematics (tau_2 = 0) and drop out;
!  - below the cut: the three-parton Born (MATTHR) reweighted by the
!    one-loop SCET cumulant at T_cut (hard, two jet, beam and soft
!    functions; mod_slicing_scet), O(alpha_s^2).
! T_2 (geometric measure, Breit frame, proton along +z):
!   T_2 = min( min_j (E_j - p_zj), min_{j<k} (E_j + E_k - |p_j + p_k|) ).
! Observable: tau_zQ = 1 - (2/Q) sum_{p_z<0} |p_z| (= tau_1^b) in bins
! between 0.05 and 0.5, which vanish at 1+1 Born kinematics and avoid the
! IR-unsafe edge tau_zQ = 1 (empty current hemisphere).
!
! Usage: tau2_nlo -pdf NAME -x X -Q2 Q2 -s S -nev N -seed1 I -seed2 J [-npow1 2 -npow2 4]
!----------------------------------------------------------------------
program tau2_nlo
  use types, only: dp
  use mod_parameters, only: nflav, NC, CC, noZ, Zonly, intonly, neutrino, positron
  use tau2_run
  use mod_slicing_scet, only: soft_tol, pdf_mask, scet_set_colour, CF, CA, TF, soft_table_init, soft_ncalls, &
       & beam_table_init, beam_ncalls
  use sub_defs_io
  implicit none
  character(len=100) :: pdf
  character(len=300) :: softtab
  real(dp) :: hw
  logical :: beamtab
  real(dp) :: s, npow1, npow2, cf_in, ca_in, tf_in, cutoff
  integer :: nev, seed1, seed2
  external :: DISENTFULL

  pdf = string_val_opt('-pdf', 'NNPDF30_nlo_as_0118')
  xfix = dble_val_opt('-x', 0.01_dp)
  Q2fix = dble_val_opt('-Q2', 400.0_dp)
  s = dble_val_opt('-s', 4 * 27.5_dp * 920.0_dp)
  nev = int_val_opt('-nev', 100000)
  seed1 = int_val_opt('-seed1', 12345)
  seed2 = int_val_opt('-seed2', 67890)
  npow1 = dble_val_opt('-npow1', 2.0_dp)
  npow2 = dble_val_opt('-npow2', 4.0_dp)
  soft_tol = dble_val_opt('-softtol', 1e-9_dp)
  ! tabulated soft function: table file (built there if missing), node spacing in w (in v: twice that)
  softtab = string_val_opt('-softtable', '')
  hw = dble_val_opt('-softtable-h', 0.025_dp)
  beamtab = log_val_opt('-beamtable')      ! tabulate the beam coefficients at the fixed Q
  pdf_mask = int_val_opt('-pdfmask', 0)   ! 1: quarks only, 2: gluon only (diagnostics)
  ndebug = int_val_opt('-debug', 0)       ! print the cumulant pieces for this many Born events
  boostY = dble_val_opt('-boostY', 0.0_dp) ! measure in a frame boosted along z (diagnostics)
  smin = dble_val_opt('-smin', 0.0_dp)     ! jets away from the beam: s_1J > smin (diagnostics)
  measure = merge(1, 0, string_val_opt('-measure', 'geo') == 'inv')   ! geo (default) or inv
  keepall = log_val_opt('-keepall')        ! keep Born/below-cut of aborted events (old behaviour)
  cutoff = dble_val_opt('-cutoff', 1e-8_dp)  ! DISENT's technical cutoff
  ! colour factors (diagnostics: e.g. -CA 0 or -TR 0 isolate colour structures; DISENT gets the same)
  cf_in = dble_val_opt('-CF', 4.0_dp/3.0_dp)
  ca_in = dble_val_opt('-CA', 3.0_dp)
  tf_in = dble_val_opt('-TR', 0.5_dp)
  call scet_set_colour(cf_in, ca_in, tf_in)

  nflav = 5; NC = .true.; CC = .false.; noZ = .true.; Zonly = .false.
  intonly = .false.; neutrino = .false.; positron = .false.
  if (softtab /= '') call soft_table_init(trim(softtab), hw, 2 * hw)
  call InitPDFsetByName(trim(pdf))
  call InitPDF(0)
  if (beamtab) call beam_table_init(sqrt(Q2fix), 0.999_dp * xfix, 0.001_dp)
  write(*,'(a,a,a,f8.5,a,f9.2,a,f10.1,a,i10,a,i2,a,3f9.5,a,i2)') ' pdf ', trim(pdf), '  x ', xfix, '  Q2 ', Q2fix, '  s ', s, &
       & '  nev ', nev, '  pdfmask ', pdf_mask, '  CF CA TR ', CF, CA, TF, '  measure ', measure

  call DISENTFULL(nev, s, 5, slice_user, slice_cuts, seed1, seed2, npow1, npow2, &
       & cutoff, 2, slice_muf, CF, CA, TF, .false.)
  write(*,'(a,2i14)') ' soft function G: table / direct evaluations ', soft_ncalls
  write(*,'(a,2i14)') ' beam coefficients: table / direct evaluations ', beam_ncalls
  call report()
end program tau2_nlo
