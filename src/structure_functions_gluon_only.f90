!======================================================================
!
! Summary of our understanding of the coefficient functions
! ---------------------------------------------------------
!
! - F1 = (F2 - FL)/2x
!
! - hep-ph/0504042 (MMV): gives F2 and FL results (electromagnetic)
!
!   (1/x) F_a = C_{a,ns} \otimes q_{ns}
!                   +  <e^2> (C_{a,q} \otimes q_s  + C_{a,g} \otimes g)
!
!   with z = 2,L, where:
!
!   * <e^2> is the average electric charge
!   * q_{ns} = ???? (we're supposed to deduce it from Eq.(4.2))
!   * C_{a,q} = C_{a,ns} + C_{a,ps};
!   * C_{2,ps}^{(0)} = C_{2,ps}^{(1)} = 0 (and presumably for FL too)
!
! - http://www.sciencedirect.com/science/article/pii/055032139290087R#
!   (Zijlstra & van Neerven) has the second-order coefficient
!   functions. 
!
! - from 1109.3717 (Zaro & Co):
!
!   * q_{ns,i}^+ = (q_i + qbar_i) - q_s
!   * q_{ns,i}^- = (q_i - qbar_i) - q_{ns}^v
!   * q_{ns}^v   = \sum_{i=1}^{n_f} (q_i - qbar_i)
!   * q_s        = \sum_{i=1}^{n_f} (q_i + qbar_i)
!
!   That, together with [C_{a,q} = C_{a,ns} + C_{a,ps}] means that the combination
! 
!       q_{ns,j}^+ * C_{i,ns}^+ + q_s * C_{i,q}
!
!   reduces to 
!
!       (q_j+qbar_j) * C_{i,ns}^+ + q_s * C_{i,ps}
!
module structure_functions_gluon_only
  use pdf_representation
  use streamlined_interface
  use splitting_functions
  use coefficient_functions_holder
  use qcd
  use structure_functions
  implicit none

  private
  public :: InitStrFct_gluon_only, StartStrFct_gluon_only
  public :: F_LO_gluon_only, F_NLO_gluon_only
!  public :: muR, muF, quark_masses_sf
  
  ! holds the coefficient functions
  type(coef_holder), save :: ch
  
  ! indices for the different structure functions
  !integer, parameter :: F1Wp= 1, F2Wp= 2, F3Wp= 3
  !integer, parameter :: F1Wm=-1, F2Wm=-2, F3Wm=-3
  !integer, parameter :: F1Z = 4, F2Z = 5, F3Z = 6
  !integer, parameter :: F1EM = -4, F2EM = -5
  !integer, parameter :: F1gZ = 0, F2gZ = -6, F3gZ = 7
  
  ! constants and fixed parameters
  !real(dp), parameter :: quark_masses_sf(3:6) = &
  !     & (/ 0.95_dp, 1.275_dp, 4.18_dp, 173.21_dp/)
  logical,  save      :: use_mass_thresholds
  real(dp), parameter :: viW = 1/sqrt(two), aiW = viW
  real(dp), save      :: vi2_ai2_avg_W, two_vi_ai_avg_W
  real(dp), save      :: vi2_ai2_Z_down, vi2_ai2_Z_up
  real(dp), save      :: two_vi_ai_Z_down, two_vi_ai_Z_up
  real(dp), save      :: two_vi_Z_down, two_vi_Z_up
  real(dp), save      :: two_ai_Z_down, two_ai_Z_up
  real(dp), parameter :: e_up = 2.0_dp/3.0_dp, e_down = -1.0_dp/3.0_dp
  real(dp), parameter :: e2_up = 4.0_dp/9.0_dp, e2_down = 1.0_dp/9.0_dp
  ! these log terms are only used for scale_choice = 0,1, to speed up the code
  real(dp), save      :: log_muF2_over_Q2, log_muR2_over_Q2
  real(dp), save      :: Qmin
  !real(dp), public, save :: toy_alphas_Q0 = 0.1185_dp
  !type(running_coupling), public, save :: toy_coupling

  ! scale choices
  real(dp), save :: cst_mu, xmuR, xmuF
  integer, save  :: scale_choice

  logical, save  :: use_sep_orders
  logical, save  :: exact_coefs
  integer :: nf_lcl, order_setup


contains

  !----------------------------------------------------------------------
  ! Setup of constants and parameters needed for structure functions
  subroutine StartStrFct_gluon_only(rts, order_max, nflav, xR, xF, sc_choice, cmu, param_coefs, &
       & Qmin_PDF, wmass, zmass)
    real(dp), intent(in) :: rts
    integer, intent(in)  :: order_max
    integer, optional    :: nflav, sc_choice
    real(dp), optional   :: xR, xF, cmu, Qmin_PDF, wmass, zmass
    logical, optional    :: param_coefs
    !----------------------------------------------------------------------
    real(dp) :: sin_thw_sq, ymax, dy, minQval, maxQval, dlnlnQ, mw, mz
    integer  :: nloop, order

    ! take sensible default value for mw and mz
    mw = 80.398_dp
    mz = 91.187_dp
    ! otherwise use user input
    if(present(wmass)) mw = wmass
    if(present(zmass)) mz = zmass
    ! compute sin(\theta_w) from W/Z mass
    sin_thw_sq = 1.0_dp - (mw/mz)**2
    
    ! evaluate parameters needed for the structure functions
    ! cf. Eqs. (3.10+11) of 1109.3717
    vi2_ai2_Z_up     = one/four + (half - (four/three) * sin_thw_sq)**2
    vi2_ai2_Z_down   = one/four + (half - (two/three)  * sin_thw_sq)**2
    two_vi_ai_Z_up   = half - (four/three) * sin_thw_sq
    two_vi_ai_Z_down = half - (two/three)  * sin_thw_sq

    ! cf. Eq.3.20 + 3.17 & 3.18
    ! 1/nf \sum_j=1^nf (vj^2 + aj^2)
    vi2_ai2_avg_W = viW**2 + aiW**2
    ! 1/nf \sum_j=1^nf 2*vj*aj
    two_vi_ai_avg_W = two * viW * aiW

    ! these parameters are needed explicitly for the gamma-Z structure function
    two_vi_Z_up = one - (8.0_dp/three) * sin_thw_sq
    two_vi_Z_down   = -one + (four/three) * sin_thw_sq
    two_ai_Z_up = one
    two_ai_Z_down   = -one
    
    ! default settings
    xmuR = one
    xmuF = one
    scale_choice = 1
    exact_coefs  = .false.
    cst_mu       = zero
    Qmin         = one
    
    ! change to user input if specified
    if(present(xR)) xmuR=xR
    if(present(xF)) xmuF=xF
    if(present(sc_choice)) scale_choice=sc_choice
    if(present(cmu)) cst_mu=cmu
    if(present(param_coefs)) exact_coefs = .not.param_coefs
    if(present(Qmin_PDF)) Qmin=Qmin_PDF

    ! Streamlined initialization
    ! including  parameters for x-grid
    order = -6 
    ymax  = 16.0_dp
    dy    = 0.05_dp  ! dble_val_opt("-dy",0.1_dp)
    dlnlnQ = dy/4.0_dp
    nloop = 3
    minQval = min(xmuF*Qmin, Qmin)
    maxQval = max(xmuF*rts, rts)

    ! initialise the grid and dglap holder
    call hoppetStartExtended(ymax,dy,minQval,maxQval,dlnlnQ,nloop,&
         &         order,factscheme_MSbar)

    if (present(nflav)) then
       ! if nflav is present, then use a fixed flavour number
       call qcd_SetNf(nflav)
       call InitCoefHolder(grid, ch, order_max, exact_coefs)
       use_mass_thresholds = .false.
       nf_lcl = nf_int
       write(6,*) "nf_lcl = ", nf_lcl
    else
       ! otherwise, use variable flavour number scheme
       call InitCoefHolder(grid, ch, order_max, exact_coefs, nflo=3, nfhi=5)
       use_mass_thresholds = .true.
       ! and start with a sensible local nf (which will be 5 here) 
       nf_lcl = nf_int
    endif
    
    ! if mass thresholds are used in the structure functions, only scale_choice 0 and 1 are allowed
    if ((use_mass_thresholds).and.(scale_choice.gt.1)) then
       call wae_error('StartStrFct', 'illegal value for scale_choice with mass thresholds turned on', intval = scale_choice)
    end if

  end subroutine StartStrFct_gluon_only

  
  !----------------------------------------------------------------------
  ! Initialize the structure functions up to specified order
  ! this requires the PDF to have been set up beforehand, and filled in tables(0)
  subroutine InitStrFct_gluon_only(order, separate_orders)
    integer, intent(in) :: order
    logical, optional :: separate_orders

    ! To turn off b quarks completely (only for testing and comparison)
    ! uncomment following two lines:
    ! tables(0)%tab(:,-5,:) = zero
    ! tables(0)%tab(:,+5,:) = zero
! gluon only
     tables(0)%tab(:,:-1,:) = zero
     tables(0)%tab(:,+1:,:) = zero

    ! default is to sum each order before convolutions (it's faster)
    if (present(separate_orders)) then
       use_sep_orders = separate_orders
    else
       use_sep_orders = .false.
    endif

    if (use_sep_orders) then
       ! First we treat the case where we want to separate out each order in
       ! the final structure functions.
       ! This is slower, but is needed e.g. for VBFH
          ! if scale_choice = 0,1 use fast implementation
          ! tables is saved as an array in Q, and only tables(0),
          ! tables(1), tables(2), tables(4), tables(8) are non zero
          if (order.ge.1) call set_LO_structure_functions()
          if (order.ge.2) call set_NLO_structure_functions()
    else
       ! Now set up the default case where we sum up everything right away.
       if (scale_choice.gt.1) & ! only allow for scale_choice = 0,1
            & call wae_error('InitStrFct', 'illegal value for scale_choice', intval = scale_choice)
       call set_structure_functions_upto(order)
    endif

    ! To rescale PDF by the N3LO F2 structure function (evaluated at 8 GeV)
    ! as a check of the size of missing N3LO PDFs, uncomment following line:
    ! call rescale_pdf_nnlo(8.0_dp,order)
    ! call rescale_pdf_n3lo(8.0_dp,order)

    setup_done(1:) = .true.
    order_setup = order
    
  end subroutine InitStrFct_gluon_only

  !----------------------------------------------------------------------
  ! set up structure functions for scale_choice = 0, 1, adding up each order
  subroutine set_structure_functions_upto (order)
    integer, intent(in) :: order
    integer :: iQ
    real(dp) :: Q, as2pi
    real(dp) :: f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f2(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: fL(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f3(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f2_fl11(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: fL_fl11(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO_f (0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PNLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO2_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PNNLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLONLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PNLOLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO3_f(0:grid%ny,ncompmin:ncompmax)
    
    do iQ = 0, tables(0)%nQ
       
       Q = tables(0)%Q_vals(iQ)
       call EvalPdfTable_Q(tables(0),muF(Q),f)

       call set_scale_logs(Q)
       
       as2pi = alphasLocal(muR(Q)) / (twopi)
       
       if (use_mass_thresholds) then
          call use_vfns(f, Q)
       endif

       ! start with LO
       if (order.ge.1) then
          tables(14)%tab(:,:,iQ) = structure_function_general(ch%C2LO*f, ch%CLLO*f, ch%C3LO*f)
       endif
       
       ! now add NLO terms
       if (order.ge.2) then
          if ((scale_choice.eq.1).and.(xmuR.eq.one).and.(xmuF.eq.one)) then
             ! For central scale with scale_choice=1, we don't need to do all
             ! the convolutions of the splitting functions
             f2 = ch%C2NLO * f
             fL = ch%CLNLO * f
             f3 = ch%C3NLO * f
          else
             ! do the convolution with the coefficient functions and also the
             ! corresponding splitting-function contributions when scales
             ! are not equal to Q
             PLO_f = dh%P_LO * f
             f2 = CxNLO_with_logs(ch%C2LO, ch%C2NLO, f, PLO_f)
             fL = CxNLO_with_logs(ch%CLLO, ch%CLNLO, f, PLO_f)
             f3 = CxNLO_with_logs(ch%C3LO, ch%C3NLO, f, PLO_f)
          endif

          tables(14)%tab(:,:,iQ) = tables(14)%tab(:,:,iQ) + &
               & as2pi * structure_function_general(f2, fL, f3)
       endif

    end do
  end subroutine set_structure_functions_upto

  
  !----------------------------------------------------------------------
  ! set up LO structure functions for scale_choice = 0, 1
  subroutine set_LO_structure_functions()
    integer :: iQ
    real(dp) :: f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: Q

    ! all coefficient functions at LO are delta functions (F2, FL and F3),
    ! so simply pass table(0) for each of the pieces
    do iQ = 0, tables(0)%nQ
       Q = tables(0)%Q_vals(iQ)
       ! explicitly evaluate the PDF at scale muF(Q)
       call EvalPdfTable_Q(tables(0),muF(Q),f)
       if (use_mass_thresholds) then
          call use_vfns(f, Q)
       endif
       tables(14)%tab(:,:,iQ) = structure_function_general(ch%C2LO*f, ch%CLLO*f, ch%C3LO*f)
    end do
    
  end subroutine set_LO_structure_functions
  
  !----------------------------------------------------------------------
  ! set up NLO structure functions for scale_choice = 0, 1
  subroutine set_NLO_structure_functions()
    integer :: iQ
    real(dp) :: f (0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f2(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: fL(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: f3(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: PLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: Q
    
    do iQ = 0, tables(0)%nQ
       
       Q = tables(0)%Q_vals(iQ)
       call EvalPdfTable_Q(tables(0),muF(Q),f)
       call set_scale_logs(Q)
       
       if (use_mass_thresholds) then
          call use_vfns(f, Q)
       endif

       if ((scale_choice.eq.1).and.(xmuR.eq.one).and.(xmuF.eq.one)) then
       ! For central scale with scale_choice=1, we don't need to do all
       ! the convolutions of the splitting functions
          f2 = ch%C2NLO * f
          fL = ch%CLNLO * f
          f3 = ch%C3NLO * f
       else
       ! do the convolution with the coefficient functions and also the
       ! corresponding splitting-function contributions when scales
       ! are not equal to Q
          PLO_f = dh%P_LO * f
          f2 = CxNLO_with_logs(ch%C2LO, ch%C2NLO, f, PLO_f)
          fL = CxNLO_with_logs(ch%CLLO, ch%CLNLO, f, PLO_f)
          f3 = CxNLO_with_logs(ch%C3LO, ch%C3NLO, f, PLO_f)
       endif
          
       tables(15)%tab(:,:,iQ) = structure_function_general(f2, fL, f3)

    end do

  end subroutine set_NLO_structure_functions

  !----------------------------------------------------------------------
  ! returns the convolution of coefficient and splitting functions
  ! for NLO, with scale dependence; leading factor of as2pi left out.
  !
  ! This routine assumes that set_scale_logs(Q) has been called
  ! beforehand.
  function CxNLO_with_logs(CxLO, CxNLO, f, PLO_f) result(res)
    real(dp),        intent(in) :: CxLO
    type(split_mat), intent(in) :: CxNLO
    real(dp),        intent(in) :: f    (0:grid%ny,ncompmin:ncompmax)
    real(dp),        intent(in) :: PLO_f(0:grid%ny,ncompmin:ncompmax)
    real(dp)                    :: res  (0:grid%ny,ncompmin:ncompmax)
    !----------------------------------------------------------------

    res = CxNLO * f - (CxLO * log_muF2_over_Q2) * PLO_f
  end function CxNLO_with_logs
  
 !----------------------------------------------------------------------
 ! Structure function valid up to NNLO
  function structure_function_general(C2_f, CL_f, C3_f) result(res)
    real(dp), intent(in) :: C2_f(0:,ncompmin:), CL_f(0:,ncompmin:), C3_f(0:,ncompmin:)
    real(dp) :: C2_f_fl11(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: CL_f_fl11(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: res(0:ubound(C2_f,dim=1), -6:7)
    C2_f_fl11 = zero
    CL_f_fl11 = zero
    res = structure_function_general_full(C2_f, CL_f, C3_f, C2_f_fl11, CL_f_fl11)
  end function structure_function_general


 !----------------------------------------------------------------------
 ! Structure function including N3LO fl_11 piece for neutral current
  function structure_function_general_full(C2_f, CL_f, C3_f, C2_f_fl11, CL_f_fl11) result(res)
    real(dp), intent(in) :: C2_f(0:,ncompmin:), CL_f(0:,ncompmin:), C3_f(0:,ncompmin:)
    real(dp), intent(in) :: C2_f_fl11(0:,ncompmin:), CL_f_fl11(0:,ncompmin:)
    real(dp) :: C2_f_NC(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: CL_f_NC(0:grid%ny,ncompmin:ncompmax)
    real(dp)             :: res(0:ubound(C2_f,dim=1), -6:7)
    !----------------------
    ! not just up and down, but the sum 
    real(dp) :: u(0:grid%ny), d(0:grid%ny), ubar(0:grid%ny), dbar(0:grid%ny) 
    real(dp) :: two_xvals(0:grid%ny)
    integer  :: nf_save
    
    two_xvals = two*xValues(grid) ! move this outside at some point
    C2_f_NC = C2_f + C2_f_fl11
    CL_f_NC = CL_f + CL_f_fl11

    res = zero
    !--- deal with Z case -----------------------------------------
    res(:,F2Z) = (dlike(C2_f_NC) + dbarlike(C2_f_NC))*vi2_ai2_Z_down + &
         &       (ulike(C2_f_NC) + ubarlike(C2_f_NC))*vi2_ai2_Z_up

    ! temporarily put FL into F1;
    res(:,F1Z) = (dlike(CL_f_NC) + dbarlike(CL_f_NC))*vi2_ai2_Z_down + &
         &       (ulike(CL_f_NC) + ubarlike(CL_f_NC))*vi2_ai2_Z_up
    ! then convert to F1
    res(:,F1Z) = (res(:,F2Z) - res(:,F1Z)) / two_xvals

    res(:,F3Z) = (dlike(C3_f) - dbarlike(C3_f))*two_vi_ai_Z_down + &
         &       (ulike(C3_f) - ubarlike(C3_f))*two_vi_ai_Z_up
    res(:,F3Z) = res(:,F3Z)/xValues(grid)

    !--- deal with EM case -----------------------------------------
    res(:,F2EM) = (dlike(C2_f_NC) + dbarlike(C2_f_NC))*e2_down + &
         &        (ulike(C2_f_NC) + ubarlike(C2_f_NC))*e2_up
    
    ! temporarily put FL into F1;
    res(:,F1EM) = (dlike(CL_f_NC) + dbarlike(CL_f_NC))*e2_down + &
         &        (ulike(CL_f_NC) + ubarlike(CL_f_NC))*e2_up
    ! then convert to F1
    res(:,F1EM) = (res(:,F2EM) - res(:,F1EM)) / two_xvals
    
    !--- deal with gamma-Z case -----------------------------------------
    ! equations taken from section 19 of PDG (19.18)
    res(:,F2gZ) = (dlike(C2_f_NC) + dbarlike(C2_f_NC))*e_down*two_vi_Z_down + &
         &        (ulike(C2_f_NC) + ubarlike(C2_f_NC))*e_up * two_vi_Z_up
    
    ! temporarily put FL into F1;
    res(:,F1gZ) = (dlike(CL_f_NC) + dbarlike(CL_f_NC))*e_down*two_vi_Z_down + &
         &        (ulike(CL_f_NC) + ubarlike(CL_f_NC))*e_up * two_vi_Z_up
    ! then convert to F1
    res(:,F1gZ) = (res(:,F2gZ) - res(:,F1gZ)) / two_xvals
    
    res(:,F3gZ) = (dlike(C3_f) - dbarlike(C3_f))*e_down*two_ai_Z_down + &
         &        (ulike(C3_f) - ubarlike(C3_f))*e_up * two_ai_Z_up
    res(:,F3gZ) = res(:,F3gZ)/xValues(grid)

    ! for the W cases, it only makes sense to sum over an even number
    ! of light flavours; so save the actual number of flavours, switch
    ! the module-local nf_lcl variable to the nearest (lower) event number
    ! for our W calculations; switch back later
    nf_save = nf_lcl
    nf_lcl = (nf_lcl/2) * 2
    
    !--- deal with Wp case -----------------------------------------
    res(:,F2Wp) = (ulike(C2_f) + dbarlike(C2_f))*vi2_ai2_avg_W
    ! temporarily put FL into F1;
    res(:,F1Wp) = (ulike(CL_f) + dbarlike(CL_f))*vi2_ai2_avg_W
    ! then convert to F1
    res(:,F1Wp) = (res(:,F2Wp) - res(:,F1Wp)) / two_xvals
    res(:,F3Wp) = (ulike(C3_f) - dbarlike(C3_f))*two_vi_ai_avg_W
    res(:,F3Wp) = res(:,F3Wp)/xValues(grid)

    !--- deal with Wm case -----------------------------------------
    res(:,F2Wm) = (dlike(C2_f) + ubarlike(C2_f))*vi2_ai2_avg_W
    ! temporarily put FL into F1;
    res(:,F1Wm) = (dlike(CL_f) + ubarlike(CL_f))*vi2_ai2_avg_W
    ! then convert to F1
    res(:,F1Wm) = (res(:,F2Wm) - res(:,F1Wm)) / two_xvals
    res(:,F3Wm) = (dlike(C3_f) - ubarlike(C3_f))*two_vi_ai_avg_W
    res(:,F3Wm) = res(:,F3Wm)/xValues(grid)

    ! reset nf_lcl to the full (possibly odd) saved value
    nf_lcl = nf_save
    
    ! overall factor of two that we still haven't fully looked into as
    ! of 2015-02-24 [GPS TMP]
    res = res * two
    !! GPS+AK TMP: we have just included a factor of 2 but this plainly
    !! should not be present for the electromagnetic part; so here
    !! we eliminate it again...
    res(:,F2EM) = half * res(:,F2EM)
    res(:,F1EM) = half * res(:,F1EM)

  end function structure_function_general_full
  

  !----------------------------------------------------------------------
  function ulike(f) result(res)
    real(dp), intent(in) :: f(0:,ncompmin:)
    real(dp)             :: res(0:ubound(f,dim=1))
    res = sum(f(:, 2: nf_lcl: 2), dim=2)
  end function ulike
  !----------------------------------------------------------------------
  function dlike(f) result(res)
    real(dp), intent(in) :: f(0:,ncompmin:)
    real(dp)             :: res(0:ubound(f,dim=1))
    res = sum(f(:, 1: nf_lcl: 2), dim=2)
  end function dlike
  !----------------------------------------------------------------------
  function ubarlike(f) result(res)
    real(dp), intent(in) :: f(0:,ncompmin:)
    real(dp)             :: res(0:ubound(f,dim=1))
    res = sum(f(:,-2:-nf_lcl:-2), dim=2)
  end function ubarlike
  !----------------------------------------------------------------------
  function dbarlike(f) result(res)
    real(dp), intent(in) :: f(0:,ncompmin:)
    real(dp)             :: res(0:ubound(f,dim=1))
    res = sum(f(:,-1:-nf_lcl:-2), dim=2)
  end function dbarlike
  
  
  !----------------------------------------------------------------------
  real(dp) function alphasLocal(muR)
    real(dp), intent(in) :: muR
    real(dp) :: muR_lcl
    muR_lcl = max(muR,Qmin)
    ! we use alphas from the LHAPDF PDF
    ! alphasLocal = alphasPDF(muR_lcl)
    ! we use alphas from HOPPET
    alphasLocal = Value(coupling, muR_lcl)
  end function alphasLocal
  

  !----------------------------------------------------------------------
  ! F_LO
  ! calculate the leading order structure function at x, muF
  !
  function F_LO_gluon_only (y, Q, muR, muF) result(res)
    real(dp), intent(in)  :: y, Q, muR, muF
    real(dp) :: res(-6:7)
    real(dp) :: Q_or_muF
    
    ! if scale_choice = 0,1, then evaluate tables at Q
    ! (since it is saved as an array in Q)
    Q_or_muF = Q
    ! if scale_choice >= 2, then evaluate tables at muF
    ! (since it is saved as an array in muF)
    if (scale_choice.ge.2) Q_or_muF = muF

    call EvalPdfTable_yQ(tables(14), y, Q_or_muF, res)
    
  end function F_LO_gluon_only

  !----------------------------------------------------------------------
  ! F_NLO
  ! calculate the next-to-leading order structure function at x, muF
  !
  ! LRQ2 == ln muR^2/Q^2
  ! LFQ2 == ln muF^2/Q^2
  !
  function F_NLO_gluon_only (y, Q, muR, muF) result(res)
    real(dp), intent(in)  :: y, Q, muR, muF
    real(dp) :: res(-6:7), as2pi, LFQ2
    real(dp) :: C1f(-6:7), C0P0f(-6:7)
    real(dp) :: Q_or_muF

    as2pi = alphasLocal(muR) / (twopi)

    ! if scale_choice = 0,1, then evaluate tables at Q
    ! (since it is saved as an array in Q)
    Q_or_muF = Q
    ! if scale_choice >= 2, then evaluate tables at muF
    ! (since it is saved as an array in muF)
    if (scale_choice.ge.2) Q_or_muF = muF
    
    ! C_NLO x f (x) in C1f(:) 
    call EvalPdfTable_yQ(tables(15), y, Q_or_muF, C1f)
    res = C1f
    
    ! if scale_choice = 0,1 then this term is already taken care of
    if (scale_choice.ge.2) then
       LFQ2 = two*log(muF/Q)
       ! C_LO x P_LO x f (x) in C0P0f(:)
       call EvalPdfTable_yQ(tables(3), y, Q_or_muF, C0P0f)
       res = res - C0P0f * LFQ2
    endif
    
    res = res * as2pi
    
  end function F_NLO_gluon_only

  !----------------------------------------------------------------------
  ! set_scale_logs is only used for scale_choice = 0,1
  subroutine set_scale_logs(Q)
    real(dp), intent(in) :: Q

    log_muF2_over_Q2 = two * log(muF(Q) / Q)
    log_muR2_over_Q2 = two * log(muR(Q) / Q)
  end subroutine set_scale_logs
  
  subroutine use_vfns (f, Q)
    real(dp), intent(in) :: Q
    real(dp), intent(inout) :: f(0:grid%ny,ncompmin:ncompmax)
    integer :: inf, current_nf

    ! set the max value for nf
    current_nf = ch%nfhi
    
    ! set pdf to zero if below threshold
    do inf = 3, 6
       if (Q.lt.quark_masses_sf(inf)) then
          f(:,-inf) = zero
          f(:, inf) = zero
          if ((inf.lt.current_nf).and.(inf.ge.ch%nflo)) current_nf = inf
       endif
    enddo

    ! set nf in coeff fct to current number of flavours above threshold
    call SetNfCoefHolder(ch, current_nf)
    
  end subroutine use_vfns

end module structure_functions_gluon_only
