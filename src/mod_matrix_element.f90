! This small routine implements the DIS cross section differential in
! x and Q (and y) as given in eq. 4.5 in the pink book (Ellis,
! Stirling, Webber). Equations also taken from PDG Review Chapter 18.
module mod_matrix_element
  use hoppet!, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF
  use streamlined_interface
  use mod_parameters
  implicit none

  private
  public :: eval_matrix_element, eval_matrix_element_new, muR_muF
  public :: parton_couplings

contains
  !----------------------------------------------------------------------
  function eval_matrix_element(order_start,order_stop, x, yDIS, Qsq) result(res)
    integer , intent(in) :: order_start,order_stop
    real(dp), intent(in) :: x, yDIS, Qsq
    real(dp)             :: res(maxscales)
    !----------------------------------------------------------------------
    real(dp) :: Qval
    real(dp) :: muRval, muFval
    real(dp) :: Fx(-6:7,4)
    real(dp) :: F1NC, F2NC, F3NC, F1CC, F2CC, F3CC, sigma(4)
    real(dp) :: overall_norm, propZ, propW, propgZ
    integer  :: i, iscale


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
       call muR_muF(x,yDIS,Qval,muRval,muFval)
       muRval = muRval * scales_mur(iscale)
       muFval = muFval * scales_muf(iscale)
!       muRval = Qval * scales_mur(iscale)
!       muFval = Qval * scales_muf(iscale)
       
       ! Compute the structure functions
       if(separate_orders) then ! Fill individual structure functions
          Fx(:,1) = F_LO(x, Qval, muRval, muFval)
          if (order_stop.ge.2) Fx(:,2) = F_NLO(x, Qval, muRval, muFval)
          if (order_stop.ge.3) Fx(:,3) = F_NNLO(x, Qval, muRval, muFval)
          if (order_stop.ge.4) Fx(:,4) = F_N3LO(x, Qval, muRval, muFval)
       else
          Fx(:,1) = StrFct(x, Qval, muRval, muFval)
          if(order_stop.ne.1) stop 'eval_matrix_element called with wrong order_stop'
       endif

       if(NC) then
          propgZ = Qsq / (Qsq + MZ**2) / sin_2thw_sq! Z propagator
          !    propgZ = (GF*MZ**2/(two*sqrt(two)*pi*alpha_em)) * Qsq / (Qsq + MZ**2)
          propZ  = propgZ**2
          do i = order_start,order_stop
             if(noZ) then
                F1NC = F1NC + Fx(iF1EM,i)
                F2NC = F2NC + Fx(iF2EM,i)
             elseif(Zonly) then ! No interference or γ
                F1NC = F1NC +   Ve2_Ae2 * propZ * Fx(iF1Z,i) 
                F2NC = F2NC +   Ve2_Ae2 * propZ * Fx(iF2Z,i) 
                F3NC = F3NC + two_Ve_Ae * propZ * Fx(iF3Z,i) 
             elseif(intonly) then ! No γ/Z
                F1NC = F1NC  - Ve * propgZ * Fx(iF1gZ,i) 
                F2NC = F2NC  - Ve * propgZ * Fx(iF2gZ,i) 
                F3NC = F3NC  - Ae * propgZ * Fx(iF3gZ,i) 
             else
                F1NC = F1NC + Fx(iF1EM,i) - (Ve * propgZ * Fx(iF1gZ,i) -   Ve2_Ae2 * propZ * Fx(iF1Z,i)) 
                F2NC = F2NC + Fx(iF2EM,i) - (Ve * propgZ * Fx(iF2gZ,i) -   Ve2_Ae2 * propZ * Fx(iF2Z,i)) 
                F3NC = F3NC               - (Ae * propgZ * Fx(iF3gZ,i) - two_Ve_Ae * propZ * Fx(iF3Z,i)) 
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
                F1CC = F1CC + propW * Fx(iF1Wp,i) 
                F2CC = F2CC + propW * Fx(iF2Wp,i) 
                F3CC = F3CC - propW * Fx(iF3Wp,i)
             else
                F1CC = F1CC + propW * Fx(iF1Wm,i) 
                F2CC = F2CC + propW * Fx(iF2Wm,i) 
                F3CC = F3CC + propW * Fx(iF3Wm,i)
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
  ! Rewritten using notation closer to PDF section 18 on Structure functions
  ! https://pdg.lbl.gov/2019/reviews/rpp2019-rev-structure-functions.pdf
  function eval_matrix_element_new(order_start,order_stop, x, yDIS, Qsq) result(res)
    integer , intent(in) :: order_start,order_stop
    real(dp), intent(in) :: x, yDIS, Qsq
    real(dp)             :: res(maxscales)
    !----------------------------------------------------------------------
    real(dp) :: Qval
    real(dp) :: muRval, muFval
    real(dp) :: Fx(-6:7,4)
    real(dp) :: F1NC, F2NC, F3NC, F1CC, F2CC, F3CC, sigma(4)
    real(dp) :: overall_norm, propZ, propW, propgZ
    ! eq. 18.9 and below. Only one of the helicities (+1,-1) will
    ! contribute to the electron giving a factor 4, which is divided
    ! by 2 since we average for unpolarised beams.
    real(dp), parameter :: etanu = four, etael = two
    integer  :: i, iscale


    Qval = sqrt(Qsq)

    res = zero

    ! Eq. 18.6 (dy -> Q2/y dQ2)
    overall_norm = twopi * alpha_em**2 / Qsq**2 / x

    do iscale = 1,Nscales
       F1NC  = zero
       F2NC  = zero
       F3NC  = zero
       F1CC  = zero
       F2CC  = zero
       F3CC  = zero
       Fx  = zero
       call muR_muF(x,yDIS,Qval,muRval,muFval)
       muRval = muRval * scales_mur(iscale)
       muFval = muFval * scales_muf(iscale)
       
       ! Compute the structure functions
       if(separate_orders) then ! Fill individual structure functions
          Fx(:,1) = F_LO(x, Qval, muRval, muFval)
          if (order_stop.ge.2) Fx(:,2) = F_NLO(x, Qval, muRval, muFval)
          if (order_stop.ge.3) Fx(:,3) = F_NNLO(x, Qval, muRval, muFval)
          if (order_stop.ge.4) Fx(:,4) = F_N3LO(x, Qval, muRval, muFval)
       else
          Fx(:,1) = StrFct(x, Qval, muRval, muFval)
          if(order_stop.ne.1) stop 'eval_matrix_element called with wrong order_stop'
       endif

       if(NC) then
          propgZ = Qsq / (Qsq + MZ**2) / sin_2thw_sq! Z propagator from eq. 18.4 using
          ! GF =  pi * alpha_em / sqrt(2) / mw**2/sin_thw_sq
          ! sin_2thw_sq = 4 * (1 - sin_thw_sq) * sin_thw_sq
          ! propgZ = (GF*MZ**2/(2*sqrt(2)*pi*alpha_em)) * Qsq / (Qsq + MZ**2)
          propZ  = propgZ**2
          if(neutrino) propZ = two * propZ ! Polarisation sum is 1, not 1/2.
          do i = order_start,order_stop
             ! Eqs. 18.11-18.12
             if(noZ) then
                F1NC = F1NC + Fx(iF1EM,i)
                F2NC = F2NC + Fx(iF2EM,i)
             elseif(Zonly) then ! No interference or γ
                F1NC = F1NC +   Ve2_Ae2 * propZ * Fx(iF1Z,i) 
                F2NC = F2NC +   Ve2_Ae2 * propZ * Fx(iF2Z,i) 
                F3NC = F3NC + two_Ve_Ae * propZ * Fx(iF3Z,i) 
             elseif(intonly) then ! No γ/Z
                F1NC = F1NC  - Ve * propgZ * Fx(iF1gZ,i) 
                F2NC = F2NC  - Ve * propgZ * Fx(iF2gZ,i) 
                F3NC = F3NC  - Ae * propgZ * Fx(iF3gZ,i) 
             else
                F1NC = F1NC + Fx(iF1EM,i) - (Ve * propgZ * Fx(iF1gZ,i) -   Ve2_Ae2 * propZ * Fx(iF1Z,i)) 
                F2NC = F2NC + Fx(iF2EM,i) - (Ve * propgZ * Fx(iF2gZ,i) -   Ve2_Ae2 * propZ * Fx(iF2Z,i)) 
                F3NC = F3NC               - (Ae * propgZ * Fx(iF3gZ,i) - two_Ve_Ae * propZ * Fx(iF3Z,i)) 
             endif
          enddo
       endif
       
       if(CC) then
          ! eq. 18.4 using
          ! propW = half * (GF * MW**2/(four * pi * alpha_em) * Qsq / (Qsq + MW**2))**2 ! W propagator 
          !       = half * (1/(sqrt(two) * four * sin_thw_sq) * Qsq / (Qsq + MW**2))**2 ! W propagator
          ! AK There seems to be some factors of 2 not fully accounted for. 
          propW = half * (1/(sqrt(two) * four * sin_thw_sq) * Qsq / (Qsq + MW**2))**2 ! W propagator
          if(neutrino) then
             propw = etanu * propW
          else
             propW = etael * propW
          endif
          do i = order_start,order_stop
             ! Eq. 18.9 seems to suggest that the W structure
             ! functions is defined with a factor two that we do not
             ! have in hoppet.
             if(neutrino) then
                if(positron) then ! Always W-
                   F1CC = F1CC + propW * two * Fx(iF1Wm,i) 
                   F2CC = F2CC + propW * two * Fx(iF2Wm,i) 
                   F3CC = F3CC - propW * two * Fx(iF3Wm,i)
                else
                   F1CC = F1CC + propW * two * Fx(iF1Wp,i) 
                   F2CC = F2CC + propW * two * Fx(iF2Wp,i) 
                   F3CC = F3CC + propW * two * Fx(iF3Wp,i)
                endif
             else
                if(positron) then ! Always W+
                   F1CC = F1CC + propW * two * Fx(iF1Wp,i) 
                   F2CC = F2CC + propW * two * Fx(iF2Wp,i) 
                   F3CC = F3CC - propW * two * Fx(iF3Wp,i)
                else
                   F1CC = F1CC + propW * two * Fx(iF1Wm,i) 
                   F2CC = F2CC + propW * two * Fx(iF2Wm,i) 
                   F3CC = F3CC + propW * two * Fx(iF3Wm,i)
                endif
             endif
          enddo
          
       endif

       sigma = compute_sigmas_new(x,yDIS,Qsq,F1NC,F2NC,F3NC,F1CC,F2CC,F3CC) * overall_norm

       res(iscale) = sigma(1) + sigma(2) ! NC + CC contributions

       NC_reduced_dsigma(iscale) = sigma(3)
       CC_reduced_dsigma(iscale) = sigma(4)

    enddo
  end function eval_matrix_element_new

  !----------------------------------------------------------------------
  ! Electroweak couplings of each incoming parton at LO, in units of
  ! the photon-exchange coupling e_q^2 (and with the same propagator
  ! factors, lepton species and gamma/Z selection as
  ! eval_matrix_element_new above, which must be kept in step):
  !   F2/x  = sum_i c2(i) f_i(x),    F3 = sum_i c3(i) f_i(x),
  ! i = -6..6 as in the PDF arrays (1 = d, 2 = u, negative = antiquarks,
  ! c2(0) = c3(0) = 0). Equivalently, for helicity-conserving
  ! lepton-parton scattering, c2 + c3 (c2 - c3) is the coupling of the
  ! same-helicity (opposite-helicity) lepton-parton configuration.
  ! Used by MATTHR in DISENT (src/libdisent.f) through disent_couplings
  ! (src/mod_disent_interface.f90). eq(-6:6) are the quark charges; the
  ! photon term is eq(i)**2 exactly, so that pure photon exchange
  ! reproduces DISENT's original arithmetic bit for bit.
  subroutine parton_couplings(Qsq, eq, c2, c3)
    real(dp), intent(in)  :: Qsq, eq(-6:6)
    real(dp), intent(out) :: c2(-6:6), c3(-6:6)
    real(dp) :: propgZ, propZ, propW, cw, sw2, c2q, c3q
    real(dp) :: two_vq, two_aq, vq2_aq2, two_vq_aq
    integer  :: i, iu, id, ngen

    c2 = zero
    c3 = zero

    if(NC) then
       propgZ = Qsq / (Qsq + MZ**2) / sin_2thw_sq
       propZ  = propgZ**2
       if(neutrino) propZ = two * propZ ! Polarisation sum is 1, not 1/2.
       sw2 = sin_thw_sq
       do i = 1, nflav
          ! Z couplings in HOPPET's normalisation (structure_functions.f90)
          if (mod(i,2) == 0) then ! up type
             two_vq    = one - (8.0_dp/three) * sw2
             two_aq    = one
             vq2_aq2   = one/four + (half - (four/three) * sw2)**2
             two_vq_aq = half - (four/three) * sw2
          else                    ! down type
             two_vq    = -one + (four/three) * sw2
             two_aq    = -one
             vq2_aq2   = one/four + (half - (two/three)  * sw2)**2
             two_vq_aq = half - (two/three)  * sw2
          endif
          if(noZ) then
             c2q = eq(i)**2
             c3q = zero
          elseif(Zonly) then
             c2q =   Ve2_Ae2 * propZ * vq2_aq2
             c3q = two_Ve_Ae * propZ * two_vq_aq
          elseif(intonly) then
             c2q = - Ve * propgZ * eq(i) * two_vq
             c3q = - Ae * propgZ * eq(i) * two_aq
          else
             c2q = eq(i)**2 - (Ve * propgZ * eq(i) * two_vq - Ve2_Ae2 * propZ * vq2_aq2)
             c3q =          - (Ae * propgZ * eq(i) * two_aq - two_Ve_Ae * propZ * two_vq_aq)
          endif
          c2( i) = c2( i) + c2q
          c2(-i) = c2(-i) + c2q
          c3( i) = c3( i) + c3q
          c3(-i) = c3(-i) - c3q
       enddo
    endif

    if(CC) then
       ! As in eval_matrix_element_new, including the factor two
       ! multiplying HOPPET's W structure functions (whose couplings are
       ! vi^2 + ai^2 = 2 vi ai = 1 per flavour).
       propW = half * (1/(sqrt(two) * four * sin_thw_sq) * Qsq / (Qsq + MW**2))**2
       if(neutrino) then
          propW = four * propW
       else
          propW = two * propW
       endif
       cw = two * propW
       ! HOPPET only uses complete generations (u,d), (c,s) for W exchange
       ngen = nflav / 2
       do i = 1, ngen
          id = 2*i - 1
          iu = 2*i
          if(neutrino .neqv. positron) then ! e+ or nu: W+ absorbed, d -> u, ubar -> dbar
             c2( id) = c2( id) + cw
             c2(-iu) = c2(-iu) + cw
             if(neutrino) then ! nu (left-handed): same helicity as d
                c3( id) = c3( id) + cw
                c3(-iu) = c3(-iu) - cw
             else              ! e+ (right-handed)
                c3( id) = c3( id) - cw
                c3(-iu) = c3(-iu) + cw
             endif
          else                              ! e- or nubar: W-, u -> d, dbar -> ubar
             c2( iu) = c2( iu) + cw
             c2(-id) = c2(-id) + cw
             if(.not.neutrino) then ! e- (left-handed): same helicity as u
                c3( iu) = c3( iu) + cw
                c3(-id) = c3(-id) - cw
             else                   ! nubar (right-handed)
                c3( iu) = c3( iu) - cw
                c3(-id) = c3(-id) + cw
             endif
          endif
       enddo
    endif
  end subroutine parton_couplings
    
  subroutine muR_muF(x,y,Q,muR,muF)
    implicit none
    real(dp), intent(in)  :: x,y,Q
    real(dp), intent(out) :: muR, muF
    real(dp) :: mu
    if (scale_choice.le.1) then
       ! if scale_choice = 0,1 then muF1(Q) = muF(Q)
       muR = sf_muR(Q)
       muF = sf_muF(Q)
    elseif (scale_choice.eq.2) then
       ! else if scale_choice=2, use Q - this is needed for scale variations
       mu  = Q
       muR = xmur * mu
       muF = xmuf * mu
    elseif (scale_choice.eq.3) then
       ! else if scale_choice=2, use pt lepton 
       mu  = Q * sqrt(1 - y) 
       muR = xmur * mu
       muF = xmuf * mu
    elseif (scale_choice.eq.4) then
       ! else if scale_choice=4, use mu^2 = Q^2 * (1-x)/x advocated by Stefano Forte
       mu  = Q * sqrt((1 - x) / x)
       muR = xmur * mu
       muF = xmuf * mu
    else
       stop 'Wrong scale choice in muR_muF'
    endif
    muR = max(Qmin,muR)
    muF = max(Qmin,muF)
  end subroutine muR_muF

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
  
  ! These are basically taken from 1206.7007 eq.1, 6, 7, 11 and eq. 18.16 in the PDG
  function compute_sigmas_new(x,y,Qsq,F1NC,F2NC,F3NC,F1CC,F2CC,F3CC) result(res)
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
       ! Reduced sigma
       res(3) = x * Qsq**2 / (two * pi * alpha_em**2 * yp) * res(1)
    endif
    if(CC) then
       res(2) = yp * F2CC + x * ym * F3CC - y**2 * FLCC
       ! Reduced sigma
       res(4) = 4.0_dp * pi * x  / GF**2 * ((mw**2 + Qsq) / mw**2)**2 * res(2)
    endif

    res = res 
  end function compute_sigmas_new
  
end module mod_matrix_element
