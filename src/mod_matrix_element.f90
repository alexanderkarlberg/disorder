! This small routine implements the DIS cross section differential in
! x and Q (and y) as given in eq. 4.5 in the pink book (Ellis,
! Stirling, Webber). Equations also taken from PDG Review Chapter 18.
module mod_matrix_element
  use hoppet_v1!, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF
  use streamlined_interface
  use mod_parameters
  implicit none

  private
  public :: eval_matrix_element

contains
  !----------------------------------------------------------------------
  function eval_matrix_element(order_start,order_stop, x, yDIS, Qsq) result(res)
    integer , intent(in) :: order_start,order_stop
    real(dp), intent(in) :: x, yDIS, Qsq
    real(dp)             :: res(maxscales)
    !----------------------------------------------------------------------
    real(dp) :: y, Qval
    real(dp) :: muRval, muFval
    real(dp) :: Fx(-6:7,4)
    real(dp) :: F1NC, F2NC, F3NC, F1CC, F2CC, F3CC, sigma(4)
    real(dp) :: overall_norm, propZ, propW, propgZ
    integer  :: i, iscale

    y = -log(x)

    Qval = sqrt(Qsq)

    res = zero

    ! Eq. 4.5 in pink book and 4.19-4.21
    ! Note that there are typos in 4.19-4.21
    ! https://www.hep.phy.cam.ac.uk/theory/webber/QCDupdates.html
    overall_norm = 4.0_dp * pi * alpha_em**2 / Qsq**2 / x

    do iscale = 1,Nscales
       F1NC  = zero
       F2NC  = zero
       F3NC  = zero
       F1CC  = zero
       F2CC  = zero
       F3CC  = zero
       Fx  = zero
       muRval = Qval * scales_mur(iscale)
       muFval = Qval * scales_muf(iscale)
       
       ! Compute the structure functions
       Fx(:,1) = F_LO(y, Qval, muRval, muFval)
       if (order_stop.ge.2) Fx(:,2) = F_NLO(y, Qval, muRval, muFval)
       if (order_stop.ge.3) Fx(:,3) = F_NNLO(y, Qval, muRval, muFval)
       if (order_stop.ge.4) Fx(:,4) = F_N3LO(y, Qval, muRval, muFval)

       if(NC) then
          propgZ = Qsq / (Qsq + MZ**2) / sin_2thw_sq! Z propagator
          !    propgZ = (GF*MZ**2/(two*sqrt(two)*pi*alpha_em)) * Qsq / (Qsq + MZ**2)
          propZ  = propgZ**2
          do i = order_start,order_stop
             if(noZ) then
                F1NC = F1NC + Fx(F1EM,i)
                F2NC = F2NC + Fx(F2EM,i)
             elseif(Zonly) then ! No interference or γ
                F1NC = F1NC +   Ve2_Ae2 * propZ * Fx(F1Z,i) 
                F2NC = F2NC +   Ve2_Ae2 * propZ * Fx(F2Z,i) 
                F3NC = F3NC + two_Ve_Ae * propZ * Fx(F3Z,i) 
             elseif(intonly) then ! No γ/Z
                F1NC = F1NC  - Ve * propgZ * Fx(F1gZ,i) 
                F2NC = F2NC  - Ve * propgZ * Fx(F2gZ,i) 
                F3NC = F3NC  - Ae * propgZ * Fx(F3gZ,i) 
             else
                F1NC = F1NC + Fx(F1EM,i) - (Ve * propgZ * Fx(F1gZ,i) -   Ve2_Ae2 * propZ * Fx(F1Z,i)) 
                F2NC = F2NC + Fx(F2EM,i) - (Ve * propgZ * Fx(F2gZ,i) -   Ve2_Ae2 * propZ * Fx(F2Z,i)) 
                F3NC = F3NC              - (Ae * propgZ * Fx(F3gZ,i) - two_Ve_Ae * propZ * Fx(F3Z,i)) 
             endif
          enddo
       endif
       
       if(CC) then
          !    propW = half * (GF * MW**2/(four * pi * alpha_em) * Qsq / (Qsq + MW**2))**2 ! W propagator 
          !          propW = half * (1/(sqrt(two) * four * sin_thw_sq) * Qsq / (Qsq + MW**2))**2 ! W propagator
          ! AK There seems to be some factors of 2 not fully accounted for. 
          propW = two * (1/(sqrt(two) * four * sin_thw_sq) * Qsq / (Qsq + MW**2))**2 ! W propagator 
          do i = order_start,order_stop
             if(positron) then ! Always W+
                F1CC = F1CC + propW * Fx(F1Wp,i) 
                F2CC = F2CC + propW * Fx(F2Wp,i) 
                F3CC = F3CC - propW * Fx(F3Wp,i)
             else
                F1CC = F1CC + propW * Fx(F1Wm,i) 
                F2CC = F2CC + propW * Fx(F2Wm,i) 
                F3CC = F3CC + propW * Fx(F3Wm,i)
             endif
          enddo
          
       endif

       sigma = compute_sigmas(x,yDIS,Qsq,F1NC,F2NC,F3NC,F1CC,F2CC,F3CC) * overall_norm

       res(iscale) = sigma(1) + sigma(2) ! NC + CC contributions

       NC_reduced_dsigma(iscale) = sigma(3)
       CC_reduced_dsigma(iscale) = sigma(4)

    enddo
  end function eval_matrix_element
    
  !----------------------------------------------------------------------
  ! mu_R as a function of Q
  real(dp) function muRlcl(Q)
    real(dp), intent(in) :: Q
    muRlcl = zero
    if (scale_choice.le.1) then
       ! if scale_choice = 0,1 then muR1(Q1,Q2) = muR(Q1)
       muRlcl = muR(Q)
    elseif (scale_choice.eq.2) then
       ! else if scale_choice=2, use Q
       muRlcl = xmur * Q
    endif
  end function muRlcl
  
  !----------------------------------------------------------------------
  ! mu_R as a function of Q
  real(dp) function muFlcl(Q)
    real(dp), intent(in) :: Q
    muFlcl = zero
    if (scale_choice.le.1) then
       ! if scale_choice = 0,1 then muF1(Q1,Q2) = muF(Q1)
       muFlcl = muF(Q)
    elseif (scale_choice.eq.2) then
       ! else if scale_choice=2, use Q
       muFlcl = xmur * Q
    endif
  end function muFlcl

  ! These are basically taken from 1206.7007 eq.1, 6, 7, 11
  function compute_sigmas(x,y,Qsq,F1NC,F2NC,F3NC,F1CC,F2CC,F3CC) result(res)
    implicit none
    real(dp), intent(in) :: x, y, Qsq, F1NC, F2NC, F3NC, F1CC, F2CC,&
         & F3CC
    real(dp) :: res(4) ! 1: NC inclusive, 2: CC inclusive, 3: NC
    ! reduced, 4: CC reduced
    real(dp) :: yp, ym, FLNC, FLCC

    FLNC = F2NC - two * x * F1NC
    FLCC = F2CC - two * x * F1CC
    
    yp = one + (one - y)**2
    ym = one - (one - y)**2

    res = zero
    
    if(NC) then
       res(1) = yp * F2NC + x * ym * F3NC - y**2 * FLNC
       res(3) = x * Qsq**2 / (two * pi * alpha_em**2 * yp) * res(1)
    endif
    if(CC) then
       res(2) = yp * F2CC + x * ym * F3CC - y**2 * FLCC
       res(4) = 4.0_dp * pi * x  / GF**2 * ((mw**2 + Qsq) / mw**2)**2 * res(2)
    endif

    ! Before we used
    ! res = y**2 * x * F1 + (one - y) * F2 + y * (one - half * y) * x * F3)
    !
    ! Which differs by a factor 2
    
    res = half * res 
  end function compute_sigmas
  
end module mod_matrix_element
