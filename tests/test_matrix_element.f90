!----------------------------------------------------------------------
! Unit test of eval_matrix_element_new (src/mod_matrix_element.f90).
!
! Run with disorder's own command line flags, e.g.
!   test_matrix_element -lo -toyQ0 2.0 -includeZ -positron [...]
! It evaluates the LO cross section d^2sigma/dx/dQ^2 at a set of (x,Q)
! points and compares it with the quark-parton-model expressions,
! written out independently here from the PDFs:
!   - charged-lepton NC (gamma, Z, gamma/Z) and the reduced cross
!     section, as in H1, arXiv:1206.7007, eqs. 1-7,
!   - CC and neutrino NC, from the V-A four-fermion interaction,
!   - the reduced CC cross section, as in 1206.7007, eq. 11, i.e.
!     4 pi x/GF^2 ((MW^2+Q^2)/MW^2)^2 d2sigma/dx/dQ2 (note that e.g. the
!     HERA combination, arXiv:1506.06042, uses half of this).
! Beyond LO (for photon exchange only) it instead checks how the
! structure functions from HOPPET are combined, which tests the F_L
! term that vanishes at LO.
! It also checks the central-scale choices of muR_muF and the scale
! labels used in the output file names.
program test_matrix_element
  use types, only: dp
  use hoppet, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF
  use mod_parameters
  use mod_matrix_element
  use mod_pdf_setup
  use test_utils
  implicit none
  real(dp), parameter :: xs(5) = [0.002_dp, 0.02_dp, 0.1_dp, 0.3_dp, 0.6_dp]
  real(dp), parameter :: Qs(5) = [3.0_dp, 10.0_dp, 30.0_dp, 80.0_dp, 150.0_dp]
  real(dp) :: x, y, Q, res(maxscales), expected, muR, muF
  integer :: ix, iQ, npoints
  character(len=80) :: tag

  call set_parameters()
  if (order_max > 1 .and. (CC .or. .not. noZ .or. neutrino .or. separate_orders)) &
       & stop 'test_matrix_element: beyond LO only photon exchange is tested'
  call setup_structure_functions()

  npoints = 0
  do ix = 1, size(xs)
     do iQ = 1, size(Qs)
        x = xs(ix)
        Q = Qs(iQ)
        y = Q**2 / (x * s)
        if (y > 0.95_dp .or. y < 0.01_dp) cycle
        npoints = npoints + 1
        write(tag,'(a,f6.3,a,f6.1)') 'x =', x, ', Q =', Q

        res = eval_matrix_element_new(1, 1, x, y, Q**2)
        call muR_muF(x, y, Q, muR, muF)
        if (order_max > 1) then
           call check_close(trim(tag)//': d2sigma/dx/dQ2 from F2, FL', res(1), &
                & photon_cross_section(x, y, Q, muR, muF), 1e-12_dp)
           cycle
        endif
        expected = lo_cross_section(x, y, Q, muF)
        call check_close(trim(tag)//': LO d2sigma/dx/dQ2', res(1), expected, 1e-5_dp)

        ! reduced cross sections as defined by H1, 1206.7007 (charged leptons)
        if (.not. neutrino) then
           if (NC) call check_close(trim(tag)//': NC reduced cross section', &
                & NC_reduced_dsigma(1), &
                & x * Q**4 / (2 * pi * alpha_em**2 * (1 + (1-y)**2)) &
                & * lo_cross_section(x, y, Q, muF, only='NC'), 1e-5_dp)
           if (CC) call check_close(trim(tag)//': CC reduced cross section', &
                & CC_reduced_dsigma(1), &
                & 4 * pi * x / GF**2 * ((mw**2 + Q**2) / mw**2)**2 &
                & * lo_cross_section(x, y, Q, muF, only='CC'), 1e-5_dp)
        endif
     enddo
  enddo
  call check_true('tested at least 8 (x,Q) points', npoints >= 8)

  call check_central_scales()
  call check_scale_labels()

  call finish_tests()

contains

  ! Photon-exchange cross section from HOPPET's structure functions,
  ! d2sigma/dx/dQ2 = 2 pi alpha^2/(x Q^4) [Y+ F2 - y^2 FL]
  real(dp) function photon_cross_section(x, y, Q, muR, muF) result(sig)
    real(dp), intent(in) :: x, y, Q, muR, muF
    real(dp) :: Fx(-6:7), F2, FL
    Fx = StrFct(x, Q, muR, muF)
    F2 = Fx(iF2EM)
    FL = F2 - 2 * x * Fx(iF1EM)
    sig = 2 * pi * alpha_em**2 / (x * Q**4) * ((1 + (1-y)**2) * F2 - y**2 * FL)
  end function photon_cross_section

  ! LO cross section d^2sigma/dx/dQ^2 in GeV^-4, for the process
  ! selected on the command line (NC and/or CC, lepton species).
  real(dp) function lo_cross_section(x, y, Q, muF, only) result(sig)
    real(dp), intent(in) :: x, y, Q, muF
    character(len=*), intent(in), optional :: only
    real(dp) :: xf(-6:6), f(-6:6), q2, yp, ym, sw2, kappa, propZ, propW
    real(dp) :: ve, ae, F2t, xF3t, F2g, F2gZ, F2Z, xF3gZ, xF3Z, eq, vq, aq
    real(dp) :: epsL, epsR, qsum, qbsum, sgn
    real(dp) :: u, ub, d, db, st, sb, c, cb
    integer :: i
    logical :: doNC, doCC

    doNC = NC; doCC = CC
    if (present(only)) then
       doNC = only == 'NC'
       doCC = only == 'CC'
    endif

    call hoppetEval(x, muF, xf)
    f = xf / x          ! number densities
    q2 = Q**2
    yp = 1 + (1-y)**2
    ym = 1 - (1-y)**2
    sw2 = 1 - (mw/mz)**2
    sig = 0.0_dp

    if (doNC .and. .not. neutrino) then
       ! kappa_Z of 1206.7007 eq. 5
       kappa = q2 / (q2 + mz**2) / (4 * sw2 * (1 - sw2))
       ve = -0.5_dp + 2 * sw2
       ae = -0.5_dp
       F2g = 0; F2gZ = 0; F2Z = 0; xF3gZ = 0; xF3Z = 0
       do i = 1, 5
          call quark_couplings(i, sw2, eq, vq, aq)
          F2g   = F2g   + x * eq**2           * (f(i) + f(-i))
          F2gZ  = F2gZ  + x * 2 * eq * vq     * (f(i) + f(-i))
          F2Z   = F2Z   + x * (vq**2 + aq**2) * (f(i) + f(-i))
          xF3gZ = xF3gZ + x * 2 * eq * aq     * (f(i) - f(-i))
          xF3Z  = xF3Z  + x * 2 * vq * aq     * (f(i) - f(-i))
       enddo
       if (noZ) then
          F2t = F2g
          xF3t = 0
       elseif (Zonly) then
          F2t  = (ve**2 + ae**2) * kappa**2 * F2Z
          xF3t = 2 * ve * ae * kappa**2 * xF3Z
       elseif (intonly) then
          F2t  = - ve * kappa * F2gZ
          xF3t = - ae * kappa * xF3gZ
       else
          F2t  = F2g - ve * kappa * F2gZ + (ve**2 + ae**2) * kappa**2 * F2Z
          xF3t = - ae * kappa * xF3gZ + 2 * ve * ae * kappa**2 * xF3Z
       endif
       ! 1206.7007 eqs. 1-2: sigma_r = F2 -/+ Y-/Y+ xF3 for e+/e-
       sgn = merge(-1.0_dp, 1.0_dp, positron)
       sig = sig + 2 * pi * alpha_em**2 / (x * q2**2) * (yp * F2t + sgn * ym * xF3t)
    endif

    if (doNC .and. neutrino) then
       ! nu q -> nu q through Z exchange:
       ! d2sigma/dx/dQ2 = GF^2/pi (MZ^2/(MZ^2+Q^2))^2
       !   sum_q [(epsL^2 + epsR^2 (1-y)^2) q + (epsR^2 + epsL^2 (1-y)^2) qbar]
       ! with epsL = T3 - e sw2, epsR = -e sw2, and L <-> R for nubar.
       propZ = (mz**2 / (mz**2 + q2))**2
       do i = 1, 5
          call quark_couplings(i, sw2, eq, vq, aq)
          epsL = aq - eq * sw2    ! aq = T3
          epsR = - eq * sw2
          if (positron) call swap(epsL, epsR)
          qsum  = (epsL**2 + epsR**2 * (1-y)**2) * f(i)
          qbsum = (epsR**2 + epsL**2 * (1-y)**2) * f(-i)
          sig = sig + GF**2 / pi * propZ * (qsum + qbsum)
       enddo
    endif

    if (doCC) then
       ! V-A: same-helicity lepton-quark scattering is flat in y, the
       ! opposite-helicity one goes as (1-y)^2. No top quark, so b
       ! (bbar) cannot be converted by a W. A charged lepton is averaged
       ! over two helicities, a neutrino has only one.
       d = f(1); u = f(2); st = f(3); c = f(4)
       db = f(-1); ub = f(-2); sb = f(-3); cb = f(-4)
       propW = (mw**2 / (mw**2 + q2))**2
       if (.not. neutrino) then
          if (.not. positron) then ! e- p -> nu X
             sig = sig + GF**2 / (2*pi) * propW * (u + c + (1-y)**2 * (db + sb))
          else                     ! e+ p -> nubar X
             sig = sig + GF**2 / (2*pi) * propW * (ub + cb + (1-y)**2 * (d + st))
          endif
       else
          if (.not. positron) then ! nu p -> e- X
             sig = sig + GF**2 / pi * propW * (d + st + (1-y)**2 * (ub + cb))
          else                     ! nubar p -> e+ X
             sig = sig + GF**2 / pi * propW * (db + sb + (1-y)**2 * (u + c))
          endif
       endif
    endif
  end function lo_cross_section

  ! Electric charge and the Z vector and axial couplings, in the
  ! normalisation v = T3 - 2 e sw2, a = T3, of quark flavour i
  ! (1=d, 2=u, 3=s, 4=c, 5=b).
  subroutine quark_couplings(i, sw2, eq, vq, aq)
    integer, intent(in)   :: i
    real(dp), intent(in)  :: sw2
    real(dp), intent(out) :: eq, vq, aq
    if (mod(i,2) == 0) then
       eq = 2.0_dp/3.0_dp;  aq = 0.5_dp
    else
       eq = -1.0_dp/3.0_dp; aq = -0.5_dp
    endif
    vq = aq - 2 * eq * sw2
  end subroutine quark_couplings

  subroutine swap(a, b)
    real(dp), intent(inout) :: a, b
    real(dp) :: t
    t = a; a = b; b = t
  end subroutine swap

  ! muR_muF for the central scale choices that do not go through
  ! HOPPET: 2: Q, 3: Q*sqrt(1-y), 4: mu^2 = Q^2 (1-x)/x.
  subroutine check_central_scales()
    real(dp) :: muR, muF, xmur_save, xmuf_save
    integer :: sc_save

    sc_save = scale_choice
    xmur_save = xmur; xmuf_save = xmuf
    xmur = 2.0_dp; xmuf = 0.5_dp

    scale_choice = 2
    call muR_muF(0.1_dp, 0.3_dp, 20.0_dp, muR, muF)
    call check_close('scale choice 2: muR = xmur Q', muR, 40.0_dp, 1e-14_dp)
    call check_close('scale choice 2: muF = xmuf Q', muF, 10.0_dp, 1e-14_dp)
    scale_choice = 3
    call muR_muF(0.1_dp, 0.36_dp, 20.0_dp, muR, muF)
    call check_close('scale choice 3: muR = xmur Q sqrt(1-y)', muR, 2 * 20 * 0.8_dp, 1e-14_dp)
    call check_close('scale choice 3: muF = xmuf Q sqrt(1-y)', muF, 0.5_dp * 20 * 0.8_dp, 1e-14_dp)
    scale_choice = 4
    call muR_muF(0.2_dp, 0.3_dp, 20.0_dp, muR, muF)
    call check_close('scale choice 4: muR^2 = xmur^2 Q^2 (1-x)/x', muR**2, 4 * 400 * 4.0_dp, 1e-14_dp)
    call check_close('scale choice 4: muF^2 = xmuf^2 Q^2 (1-x)/x', muF**2, 0.25_dp * 400 * 4.0_dp, 1e-14_dp)
    ! scales are never taken below Qmin
    scale_choice = 2
    call muR_muF(0.1_dp, 0.3_dp, 0.1_dp, muR, muF)
    call check_close('scales are cut off at Qmin (muR)', muR, Qmin, 1e-14_dp)
    call check_close('scales are cut off at Qmin (muF)', muF, Qmin, 1e-14_dp)

    scale_choice = sc_save
    xmur = xmur_save; xmuf = xmuf_save
  end subroutine check_central_scales

  ! The 7-point variation: central point first, all points distinct,
  ! 1/2 <= muR/muF <= 2, and output labels matching the scale factors.
  subroutine check_scale_labels()
    character(len=17) :: label
    integer :: i, j
    logical :: distinct
    call check_true('scale variation: central scale first', &
         & scales_mur(1) == 1.0_dp .and. scales_muf(1) == 1.0_dp)
    call check_true('scale variation: 1/2 <= muR/muF <= 2', &
         & all(scales_mur/scales_muf >= 0.5_dp .and. scales_mur/scales_muf <= 2.0_dp))
    distinct = .true.
    do i = 1, maxscales
       do j = i+1, maxscales
          distinct = distinct .and. .not. (scales_mur(i) == scales_mur(j) &
               & .and. scales_muf(i) == scales_muf(j))
       enddo
    enddo
    call check_true('scale variation: all points distinct', distinct)
    do i = 1, maxscales
       write(label,'(a,f3.1,a,f3.1)') '_μR_', scales_mur(i), '_μF_', scales_muf(i)
       call check_true('output label of scale point '//trim(scalestr(i)), &
            & trim(label) == trim(scalestr(i)))
    enddo
  end subroutine check_scale_labels

end program test_matrix_element
