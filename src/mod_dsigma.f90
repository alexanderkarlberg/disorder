module mod_dsigma
  use hoppet_v1, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF
  use mod_parameters
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
    !     generate phase space 
    call gen_phsp_born(xrand(1:2),x,y,Qsq,Qvec,jacborn,pbornlab,pbornbreit)
    
    !     skip phase space points with vanishing jacobian or with Q < Qmin
    if (.not.(jacborn.ne.0d0).and.(Qsq.gt.(Qmin**2))) return
    dsigma_all_scales = eval_matrix_element(order_min,order_max, x, y, Qsq)
    dsigma_all_scales = dsigma_all_scales * gev2pb * jacborn 
    dsigma = dsigma_all_scales(1)
    dsigma_all_scales = dsigma_all_scales * vegas_weight 
    sigma_all_scales = sigma_all_scales + dsigma_all_scales
    
    !      if(scaleuncert) then
    !         dsigma_all_scales = eval_matrix_element(order_min,order_max, x, y, Qsq, nscales)
    !         dsigma_all_scales = dsigma_all_scales * gev2pb * jacborn 
    !         dsigma = dsigma_all_scales(1)
    !         dsigma_all_scales = dsigma_all_scales * vegas_weight 
    !         sigma_all_scales = sigma_all_scales + dsigma_all_scales
    !      else
    !         dsigma_all_scales = eval_matrix_element(order_min,order_max, x, y, Qsq, 1)
    !         dsigma = dsigma_all_scales(1)
    !         dsigma = dsigma * gev2pb * jacborn
    !         dsigma_all_scales = dsigma * vegas_weight
    !      endif
    
    !     remove the rare outliers where we get dsigma = NaN (never happens)
    if (dsigma.ne.dsigma) then
       dsigma = 0d0
       dsigma_all_scales = 0d0
    endif
    
    ! Do the analysis
    if(fillplots) then
       call user_analysis(4,dsigma_all_scales*vegas_ncall,x,y,Qsq)
       call pwhgaccumup
    endif
  end function dsigma
  
  
end module mod_dsigma
