module mod_dsigma
  use hoppet_v1, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF
  use mod_parameters
  use mod_analysis
  use types
  
  implicit none
  
contains
  
  !------------------------------------------------------------
  !     dsigma function
  double precision function dsigma(xrand, vegas_weight)
    use mod_matrix_element
    use mod_phase_space
    implicit none
    !     xrand contains a vector of random numbers in [0,1]
    real(dp) :: x, y, Qsq, jacborn, Qvec(0:3)
    real(dp) :: dsigma_all_scales(maxscales)
    real(dp) :: xrand(4), vegas_weight
    integer vegas_ncall
    common/vegas_ncall/vegas_ncall
    
    dsigma = 0d0
    dsigma_all_scales = 0d0
    NC_reduced_dsigma = 0d0
    CC_reduced_dsigma = 0d0
    

    !     generate phase space 
    call gen_phsp_born(xrand(1:2),x,y,Qsq,Qvec,jacborn,pbornlab,pbornbreit)
    
    !     skip phase space points with vanishing jacobian or with Q < Qmin
    if (.not.(jacborn.ne.0d0).and.(Qsq.gt.(Qmin**2))) return
    if(separate_orders) then
       dsigma_all_scales = eval_matrix_element(order_min,order_max, x, y, Qsq)
    else
       dsigma_all_scales = eval_matrix_element(1,1, x, y, Qsq)
    endif
    dsigma_all_scales = dsigma_all_scales * gev2pb * jacborn 
    dsigma = dsigma_all_scales(1)
    dsigma_all_scales = dsigma_all_scales * vegas_weight / itmx2
    sigma_all_scales = sigma_all_scales + dsigma_all_scales

    ! Do reduced cross sections
    NC_reduced_dsigma = NC_reduced_dsigma * jacborn * vegas_weight / itmx2
    CC_reduced_dsigma = CC_reduced_dsigma * jacborn * vegas_weight / itmx2

    NC_reduced_sigma = NC_reduced_sigma + NC_reduced_dsigma
    CC_reduced_sigma = CC_reduced_sigma + CC_reduced_dsigma

    !     remove the rare outliers where we get dsigma = NaN (never happens)
    if (dsigma.ne.dsigma) then
       dsigma = 0d0
       dsigma_all_scales = 0d0
    endif

    ! Do the analysis
    if(fillplots) then
       call analysis(4,dsigma_all_scales*vegas_ncall,x,y,Qsq)
       call pwhgaccumup
    endif
  end function dsigma
  
  
end module mod_dsigma
