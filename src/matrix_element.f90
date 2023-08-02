! This small routine implements the DIS cross section differential in
! x and Q (and y) as given in eq. 4.5 in the pink book (Ellis,
! Stirling, Webber).
module matrix_element
  use hoppet_v1, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF
  use streamlined_interface
  use parameters
  use structure_functions_gluon_only
  implicit none

  private
  public :: eval_matrix_element, real_emission_matrix_element

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
    real(dp) :: overall_norm, prop
    real(dp) :: P1q1
    integer  :: i, j
    integer :: iorder

    y = -log(x)

    Qval = sqrt(Qsq)

    muRval = muR(Qval)
    muFval = muF(Qval)


    F1 = zero
    F2 = zero
    F3 = zero
    Fx = zero

    ! Eq. 4.5 in pink book and 4.19-4.21
    ! Note that there are typos in 4.19-4.21
    ! https://www.hep.phy.cam.ac.uk/theory/webber/QCDupdates.html
    overall_norm = 4.0_dp * pi * alpha_em**2 / Qsq**2
    ! Not sure how to treat the Breit-Wigner in the interference here..
    prop = Qsq / (Qsq + MZ**2) / sin_2thw_sq! Z propagator wo Breit-Wigner

    ! Compute the LO structure funtion by adding all the pieces
    ! from tables
    if(gluon_only) then
      Fx(:,1) = F_LO_gluon_only(y, Qval, muRval, muFval)
    else
      Fx(:,1) = F_LO(y, Qval, muRval, muFval)
    endif
    if (order_stop.ge.2) then
       ! Compute the NLO structure funtion by adding all the pieces
       ! from tables
       if(gluon_only) then
        Fx(:,2) = F_NLO_gluon_only(y, Qval, muRval, muFval)
       else
        Fx(:,2) = F_NLO(y, Qval, muRval, muFval)
       endif
    endif

    if (order_stop.ge.3) then
       ! Compute the NNLO structure funtion by adding all the pieces
       ! from tables
       Fx(:,3) = F_NNLO(y, Qval, muRval, muFval)
    endif

    if (order_stop.ge.4) then
       ! Compute the N3LO structure funtion by adding all the pieces
       ! from tables
       Fx(:,4) = F_N3LO(y, Qval, muRval, muFval)
    endif

    if(noZ) then
       do i = order_start,order_stop
          F1 = F1 + Fx(F1EM,i)
          F2 = F2 + Fx(F2EM,i)
       end do
              
       res = ( one + (one - yDIS)**2 ) * F1 + (1 - yDIS) / x * (F2 - two * x * F1)      

       res = res * overall_norm
    else
       ! Careful. In HOPPET a factor 2 is included for the structure
       ! functions as defined in VBF. For the EM piece this is already
       ! divided out, but we divide it out by hand for the Z and γZ
       ! pieces here.
       do i = order_start,order_stop
          if(Zonly) then ! No interference or γ
             F1 = F1 +   Ve2_Ae2 * prop**2 * Fx(F1Z,i) * half
             F2 = F2 +   Ve2_Ae2 * prop**2 * Fx(F2Z,i) * half
             F3 = F3 + two_Ve_Ae * prop**2 * Fx(F3Z,i) * half
          else if(intonly) then ! No γ/Z
             F1 = F1  - Ve * prop * Fx(F1gZ,i) * half
             F2 = F2  - Ve * prop * Fx(F2gZ,i) * half
             F3 = F3  - Ae * prop * Fx(F3gZ,i) * half
          else
             F1 = F1 + Fx(F1EM,i) - (Ve * prop * Fx(F1gZ,i) -   Ve2_Ae2 * prop**2 * Fx(F1Z,i)) * half
             F2 = F2 + Fx(F2EM,i) - (Ve * prop * Fx(F2gZ,i) -   Ve2_Ae2 * prop**2 * Fx(F2Z,i)) * half
             F3 = F3              - (Ae * prop * Fx(F3gZ,i) - two_Ve_Ae * prop**2 * Fx(F3Z,i)) * half
          endif
       end do

       res = (                  x * yDIS**2 * F1 &
           +                   (one - yDIS) * F2 &
           + x * yDIS * (one - half * yDIS) * F3)
       
       res = res * overall_norm / x
    endif
  end function eval_matrix_element

  ! Eq. 3.1, 3.2 and 3.3 of hep-ph/9704297
  function real_emission_matrix_element(x, y, Qsq, csi, z) result(res)
  implicit none
  real(dp), intent(in) :: x, y, Qsq, csi, z
  real(dp) :: res
  real(dp) :: overall_norm, q, g, Qval, as, pdf(-6:6), FT, FL
  real(dp) :: muRval, muFval

  Qval = sqrt(Qsq)
  ! Returns x / ξ * PDF(x/ξ)
  muRval = muR(Qval)
  muFval = muF(Qval)
  call EvolvePDF(x/csi,muFval,pdf)   
  pdf = pdf / x 
  !as = alphasPDF(Qval)
  as = alphasLocal(muRval)
  overall_norm = 2.0_dp * pi * alpha_em**2 / Qsq**2


  q = Qdn**2 * (pdf(1) + pdf(3) + pdf(5) + pdf(-1) + pdf(-3) + pdf(-5)) &
      & + Qup**2 * (pdf(2) + pdf(4) + pdf(-2) + pdf(-4))  

  g = pdf(0)
  Tf = Tr * (3.0_dp * Qdn**2 + 2.0_dp * Qup**2)
  
  if(gluon_only) q = 0.0_dp

  FT = as / 2.0_dp / pi * (CF * CTq(csi,z) * q + Tf * CTg(csi,z) * g)
  FL = as / 2.0_dp / pi * (CF * CLq(csi,z) * q + Tf * CLg(csi,z) * g)

  res = overall_norm * ((1.0_dp + (1.0_dp - y)**2) * FT &
      &                   + 2.0_dp * (1.0_dp - y) * FL)
  end function

!  Eq. 3.4 of hep-ph/9704297
  function CTq(csi,z) result(res)
  implicit none
  real(dp), intent(in) :: csi, z
  real(dp) :: res

  res = (csi**2 + z**2)/ (1.0_dp - csi) / (1.0_dp - z) + 2.0_dp * csi * z + 2.0_dp
  end

  function CLq(csi,z) result(res)
  implicit none
  real(dp), intent(in) :: csi, z
  real(dp) :: res

  res = 4.0_dp * csi * z
  end

  function CTg(csi,z) result(res)
  implicit none
  real(dp), intent(in) :: csi, z
  real(dp) :: res

  res = (csi**2 + (1.0_dp - csi)**2) * (z**2 + (1.0_dp - z)**2) / (z * (1.0_dp - z))
  end

  function CLg(csi,z) result(res)
  implicit none
  real(dp), intent(in) :: csi, z
  real(dp) :: res

  res = 8.0_dp * csi * (1.0_dp - csi)
  end

    
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

end module matrix_element
