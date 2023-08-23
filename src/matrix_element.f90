! This small routine implements the DIS cross section differential in
! x and Q (and y) as given in eq. 4.5 in the pink book (Ellis,
! Stirling, Webber). Equations also taken from PDG Review Chapter 18.
module matrix_element
  use hoppet_v1, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF
  use streamlined_interface
  use parameters
  implicit none

  private
  public :: eval_matrix_element, eval_matrix_element_scale_variation

contains

  !----------------------------------------------------------------------
  function eval_matrix_element(order_start,order_stop, x, yDIS, Qsq) result(res)
    integer , intent(in) :: order_start,order_stop
    real(dp), intent(in) :: x, yDIS, Qsq
    real(dp)             :: res
    !----------------------------------------------------------------------
    real(dp) :: y, Qval
    real(dp) :: muRval, muFval
    real(dp) :: Fx(-6:7,4)
    real(dp) :: F1, F2, F3
    real(dp) :: overall_norm, propgZ, propZ, propW
    integer  :: i

    y = -log(x)

    Qval = sqrt(Qsq)

    muRval = muRlcl(Qval)
    muFval = muFlcl(Qval)

    F1 = zero
    F2 = zero
    F3 = zero
    Fx = zero

    ! Eq. 4.5 in pink book and 4.19-4.21
    ! Note that there are typos in 4.19-4.21
    ! https://www.hep.phy.cam.ac.uk/theory/webber/QCDupdates.html
    overall_norm = 4.0_dp * pi * alpha_em**2 / Qsq**2

    propgZ = Qsq / (Qsq + MZ**2) / sin_2thw_sq! Z propagator
    !    propgZ = (GF*MZ**2/(two*sqrt(two)*pi*alpha_em)) * Qsq / (Qsq + MZ**2)
    propZ  = propgZ**2
    !    propW = half * (GF * MW**2/(four * pi * alpha_em) * Qsq / (Qsq + MW**2))**2 ! W propagator 
    propW = half * (1/(sqrt(two) * four * sin_thw_sq) * Qsq / (Qsq + MW**2))**2 ! W propagator 

    ! Compute the structure funtion by adding all the pieces from
    ! tables
    Fx(:,1) = F_LO(y, Qval, muRval, muFval)
    if (order_stop.ge.2) Fx(:,2) = F_NLO(y, Qval, muRval, muFval)
    if (order_stop.ge.3) Fx(:,3) = F_NNLO(y, Qval, muRval, muFval)
    if (order_stop.ge.4) Fx(:,4) = F_N3LO(y, Qval, muRval, muFval)

    if(NC) then
       do i = order_start,order_stop
          if(noZ) then
             F1 = F1 + Fx(F1EM,i)
             F2 = F2 + Fx(F2EM,i)
          elseif(Zonly) then ! No interference or γ
             F1 = F1 +   Ve2_Ae2 * propZ * Fx(F1Z,i) 
             F2 = F2 +   Ve2_Ae2 * propZ * Fx(F2Z,i) 
             F3 = F3 + two_Ve_Ae * propZ * Fx(F3Z,i) 
          elseif(intonly) then ! No γ/Z
             F1 = F1  - Ve * propgZ * Fx(F1gZ,i) 
             F2 = F2  - Ve * propgZ * Fx(F2gZ,i) 
             F3 = F3  - Ae * propgZ * Fx(F3gZ,i) 
          else
             F1 = F1 + Fx(F1EM,i) - (Ve * propgZ * Fx(F1gZ,i) -   Ve2_Ae2 * propZ * Fx(F1Z,i)) 
             F2 = F2 + Fx(F2EM,i) - (Ve * propgZ * Fx(F2gZ,i) -   Ve2_Ae2 * propZ * Fx(F2Z,i)) 
             F3 = F3              - (Ae * propgZ * Fx(F3gZ,i) - two_Ve_Ae * propZ * Fx(F3Z,i)) 
          endif
       enddo
    endif
    
    if(CC) then
       do i = order_start,order_stop
          if(positron) then ! Always W+
             F1 = F1 + propW * Fx(F1Wp,i) 
             F2 = F2 + propW * Fx(F2Wp,i) 
             F3 = F3 - propW * Fx(F3Wp,i)
          else
             F1 = F1 + propW * Fx(F1Wm,i) 
             F2 = F2 + propW * Fx(F2Wm,i) 
             F3 = F3 + propW * Fx(F3Wm,i)
          endif
       enddo
    endif

    res = (                   x * yDIS**2 * F1 &
         +                   (one - yDIS) * F2 &
         + x * yDIS * (one - half * yDIS) * F3)

    res = res * overall_norm / x
    
  end function eval_matrix_element

  !----------------------------------------------------------------------
  function eval_matrix_element_scale_variation(order_start,order_stop, x, yDIS, Qsq, iscales) result(res)
    integer , intent(in) :: order_start,order_stop,iscales
    real(dp), intent(in) :: x, yDIS, Qsq
    real(dp)             :: res(maxscales)
    !----------------------------------------------------------------------
    real(dp) :: y, Qval
    real(dp) :: muRval, muFval
    real(dp) :: Fx(-6:7,4)
    real(dp) :: F1, F2, F3
    real(dp) :: overall_norm, propZ, propW, propgZ
    integer  :: i, iscale

    if(iscales.gt.maxscales.or.iscales.lt.1) stop 'Need number of scales between 1 and 7'
    
    y = -log(x)

    Qval = sqrt(Qsq)

    res = zero

    ! Eq. 4.5 in pink book and 4.19-4.21
    ! Note that there are typos in 4.19-4.21
    ! https://www.hep.phy.cam.ac.uk/theory/webber/QCDupdates.html
    overall_norm = 4.0_dp * pi * alpha_em**2 / Qsq**2

    propgZ = Qsq / (Qsq + MZ**2) / sin_2thw_sq! Z propagator
    !    propgZ = (GF*MZ**2/(two*sqrt(two)*pi*alpha_em)) * Qsq / (Qsq + MZ**2)
    propZ  = propgZ**2
    !    propW = half * (GF * MW**2/(four * pi * alpha_em) * Qsq / (Qsq + MW**2))**2 ! W propagator 
    propW = half * (1/(sqrt(two) * four * sin_thw_sq) * Qsq / (Qsq + MW**2))**2 ! W propagator 

    do iscale = 1,iscales
       F1  = zero
       F2  = zero
       F3  = zero
       Fx  = zero
       muRval = Qval * scales_mur(iscale)
       muFval = Qval * scales_muf(iscale)
       
       ! Compute the structure functions
       Fx(:,1) = F_LO(y, Qval, muRval, muFval)
       if (order_stop.ge.2) Fx(:,2) = F_NLO(y, Qval, muRval, muFval)
       if (order_stop.ge.3) Fx(:,3) = F_NNLO(y, Qval, muRval, muFval)
       if (order_stop.ge.4) Fx(:,4) = F_N3LO(y, Qval, muRval, muFval)

       if(NC) then
          do i = order_start,order_stop
             if(noZ) then
                F1 = F1 + Fx(F1EM,i)
                F2 = F2 + Fx(F2EM,i)
             elseif(Zonly) then ! No interference or γ
                F1 = F1 +   Ve2_Ae2 * propZ * Fx(F1Z,i) 
                F2 = F2 +   Ve2_Ae2 * propZ * Fx(F2Z,i) 
                F3 = F3 + two_Ve_Ae * propZ * Fx(F3Z,i) 
             elseif(intonly) then ! No γ/Z
                F1 = F1  - Ve * propgZ * Fx(F1gZ,i) 
                F2 = F2  - Ve * propgZ * Fx(F2gZ,i) 
                F3 = F3  - Ae * propgZ * Fx(F3gZ,i) 
             else
                F1 = F1 + Fx(F1EM,i) - (Ve * propgZ * Fx(F1gZ,i) -   Ve2_Ae2 * propZ * Fx(F1Z,i)) 
                F2 = F2 + Fx(F2EM,i) - (Ve * propgZ * Fx(F2gZ,i) -   Ve2_Ae2 * propZ * Fx(F2Z,i)) 
                F3 = F3              - (Ae * propgZ * Fx(F3gZ,i) - two_Ve_Ae * propZ * Fx(F3Z,i)) 
             endif
          enddo
       endif
       
       if(CC) then
          do i = order_start,order_stop
             if(positron) then ! Always W+
                F1 = F1 + propW * Fx(F1Wp,i) 
                F2 = F2 + propW * Fx(F2Wp,i) 
                F3 = F3 - propW * Fx(F3Wp,i)
             else
                F1 = F1 + propW * Fx(F1Wm,i) 
                F2 = F2 + propW * Fx(F2Wm,i) 
                F3 = F3 + propW * Fx(F3Wm,i)
             endif
          enddo
          
       endif
       
       res(iscale) = (                   x * yDIS**2 * F1 &
            +                   (one - yDIS) * F2 &
            + x * yDIS * (one - half * yDIS) * F3)
       
       res(iscale) = res(iscale) * overall_norm / x
    enddo
  end function eval_matrix_element_scale_variation
    
  !----------------------------------------------------------------------
  real(dp) function alphasLocal(muR)
    real(dp), intent(in) :: muR
    real(dp) :: muR_lcl, alphasPDF
    muR_lcl = max(muR,Qmin)
    ! we use alphas from the LHAPDF PDF
     alphasLocal = alphasPDF(muR_lcl)
    ! we use alphas from HOPPET
    !alphasLocal = Value(coupling, muR_lcl)
   end function alphasLocal

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

end module matrix_element
