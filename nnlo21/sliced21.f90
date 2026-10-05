!-----------------------------------------------------------------------
! Below-cut part of tau_2-sliced NNLO DIS 2+1 at fixed (x, Q^2): the 2+1
! Born (photon exchange) in bins of tau_zQ, times the leading-power
! cumulant of nnlo21/lp21 for each tau_cut. With the above-cut part
! (dis31/nlo31 in mode 1) this gives dsigma/dx dQ^2 [pb/GeV^2]:
!   b0: O(alpha_s)   Born
!   b1: O(alpha_s^2) below the cut   (+ nlo31 lo above it = NLO 2+1 correction)
!   b2: O(alpha_s^3) below the cut   (+ nlo31 vi + kp + r above it = NNLO 2+1 correction)
! Born matrix elements: DISENT's MATTHR (photon exchange; alpha_s/2pi
! factored out, the normalisation of born31/me41), parton 2 the quark and 3
! the gluon (quark Born) or the antiquark (gluon Born).
!
! Usage: sliced21 part ncall itmx seed x Q2 [softtable [h [pdfmask]]]   (part = b0, b1, b2;
! tabchk: beam table at node spacing h (default 0.02) against direct evaluation)
!
! Mode 2 (ZEUS-like dijets, as dis31/nlo31 mode 2; sigma per bin in pb):
!   sliced21 part ncall itmx seed zeus <table prefix> [softtable]
! integrated over Q^2 and y with the selection on the 2+1 Born jets; the beam
! tables on a grid in Q (mu = Q per event), built once by
!   sliced21 mktab i0 i1 <prefix> [dlnQ [h]]   (nodes i0..i1, in parallel;
!                                    defaults 0.1 in ln Q, 0.02 in ln(xi/(1-xi)))
!   sliced21 tabchkg <prefix>          (grid against direct evaluation)
!
! Mode 3 (P2B, as dis31/nlo31 mode 3): sliced21 part ncall itmx seed p2b <table prefix> [softtable]
! 1+1 lab-frame jet bins with O(2+1 Born) - O(projected 1+1 Born) per event
!-----------------------------------------------------------------------
module sliced21_mod
  use nlo31_mod
  use lp21
  implicit none
  character(8) :: bpart
  integer :: nemit = 1      ! correlated sampling: emissions per Born (C1_NEMIT)
  integer :: vbal = 0       ! c1: balanced VEGAS target (C1_VBAL)
  integer, parameter :: ivb = 8   ! tau_cut column of the balanced target (1e-4)
  ! mode 2: beam-table grid in Q (lp21_grid_build/load)
  real(dp), parameter :: gqlo = 11.180339887498949_dp, gqhi = 141.42135623730951_dp, gdl = 0.1_dp
  real(dp), parameter :: gximin = 2e-3_dp, gh = 0.02_dp

contains

  real(dp) function b21_part(r, wgt) result(res)
    real(dp), intent(in) :: r(:), wgt
    real(dp) :: Q2, y, xB, eta, jac, Pk(4,6), dphi
    logical :: ok
    res = 0
    call lepton(r(1:3), Q2, y, xB, eta, jac, ok)
    if (.not. ok) return
    call breit(2, r(4:5), Q2, y, xB, eta, Pk, dphi)
    res = b21_eval(Pk, Q2, xB, eta, jac, dphi, wgt)
  end function b21_part

  ! correlated sampling (5 Oct): the below-cut 2+1 term at the 2+1 Born that
  ! psmc builds the above-cut 3+1 event from (same unit numbers: lepton2 from
  ! r(1:2), the Born from r(4:6) = psmc_gen(3, r(3:9))'s r(2:4)); with
  ! born_part (nlo31, part lo) on the same r this is the NLO 2+1 coefficient
  ! b1 + lo sampled point by point together
  real(dp) function b21c_part(r, wgt) result(res)
    real(dp), intent(in) :: r(:), wgt
    real(dp) :: Q2, y, xB, eta, jac, Pk(4,6), dphi
    logical :: ok
    res = 0
    call lepton2(r(1:2), Q2, y, xB, jac, ok)
    if (.not. ok) return
    call psmc_born(r(4:6), Pk, eta, dphi, ok)
    if (.not. ok) return
    res = b21_eval(Pk, Q2, xB, eta, jac, dphi, wgt)
  end function b21c_part

  ! the correlated NLO 2+1 integrand: lo (nlo31) + b1 at the same point; with
  ! nemit > 1 (environment C1_NEMIT) the lo term is averaged over nemit
  ! emissions from the same Born (fresh unit numbers for the channel and the
  ! emission, r(3) and r(7:9), for all but the first): each is unbiased, b1
  ! (the expensive part) is evaluated once
  ! VEGAS target for c1 with vbal (environment C1_VBAL = 1): the size of the
  ! point's whole contribution vector, sqrt(sum_b f_b^2) over the jet bins at
  ! tau_cut = tcs(ivb) (f = b1 + lo per bin, from the hcacc increments), so the
  ! grid adapts to the variance of the cancelling sum in all bins at once
  real(dp) function c1_part(r, wgt) result(res)
    real(dp), intent(in) :: r(:), wgt
    real(dp) :: h0(ncell), f
    integer :: b
    if (vbal == 0) then
       res = c1_core(r, wgt)
       return
    endif
    h0 = hcacc
    f = c1_core(r, wgt)
    res = 0
    do b = 2, nob
       res = res + ((hcacc(ivb + ntc*(b - 1)) - h0(ivb + ntc*(b - 1)))/wgt)**2
    enddo
    res = sqrt(res)
  end function c1_part

  real(dp) function c1_core(r, wgt) result(res)
    real(dp), intent(in) :: r(:), wgt
    real(dp) :: rr(size(r)), u(4), a
    integer :: m, nc
    if (nemit < 0) then
       ! stratified over psmc's channels (nemit = -1): one emission per channel
       ! from the same Born, weighted by the channel probabilities (psmc_gen uses
       ! r(1) = rr(3) only to select the channel); fresh emission numbers except
       ! for the flat channel, which takes r(4:9)
       nc = psmc_nchan(3)
       res = 0
       do m = 0, nc
          rr = r
          if (m == 0) then
             a = psmc_aflat; rr(3) = 0.5_dp*psmc_aflat
          else
             a = (1 - psmc_aflat)/nc; rr(3) = psmc_aflat + (1 - psmc_aflat)*(m - 0.5_dp)/nc
             call random_number(u(1:3)); rr(7:9) = u(1:3)
          endif
          res = res + a*born_part(rr, wgt*a)
       enddo
       res = res + b21c_part(r, wgt)
       return
    endif
    res = born_part(r, wgt/nemit)/nemit
    do m = 2, nemit
       rr = r
       call random_number(u)
       rr(3) = u(1); rr(7:9) = u(2:4)
       res = res + born_part(rr, wgt/nemit)/nemit
    enddo
    res = res + b21c_part(r, wgt)
  end function c1_core

  real(dp) function b21_eval(Pk, Q2, xB, eta, jac, dphi, wgt) result(res)
    real(dp), intent(in) :: Pk(4,6), Q2, xB, eta, jac, dphi, wgt
    real(dp), external :: alphasPDF
    real(dp) :: P(4,7), as, fpdf(-5:5), w, Q, tz
    real(dp) :: QQ, GQ, born(3), f0(3), c1(ntc,3), c2(ntc,3), val(ntc), mu(3)
    logical :: inb(nobmax)
    real(dp) :: dp2b(nobmax)
    integer :: b, k, c
    res = 0
    Q = sqrt(Q2)
    if (mode == 2) then
       klep = Pk(:,5); Q2cur = Q2
       call zeus_bins(Pk(:,1), Pk(:,2:3), 2, inb)
       if (.not. inb(1)) return
       if (trim(bpart) /= 'b0') call lp21_setq(Q)
    elseif (mode == 3) then
       klep = Pk(:,5); kout = Pk(:,6); Q2cur = Q2; xcur = xB
       call p2b_diff(Pk(:,1), Pk(:,2:3), 2, dp2b)
       if (all(dp2b == 0)) return
       if (trim(bpart) /= 'b0') call lp21_setq(Q)
    else
       tz = 1
       do k = 2, 3
          if (Pk(3,k) < 0) tz = tz + 2*Pk(3,k)/Q
       enddo
       if (tz < zlo(nzb) .or. tz >= zhi(nzb)) return
    endif
    P = 0
    P(:,1:3) = Pk(:,1:3); P(:,5) = Pk(:,4); P(:,6) = Pk(:,5); P(:,7) = Pk(:,6)
    ! DISENT's MATTHR per unit charge^2 (QQ: quark Born, GQ: gluon Born, one flavour)
    QQ = 8*(4*pi/137)**2*(dd(1,6)**2 + dd(1,7)**2 + dd(2,7)**2 + dd(2,6)**2) &
         & *16*pi**2*CF/(-4*dd(2,3)*dd(1,3)*dd(5,5))
    GQ = 8*(4*pi/137)**2*(dd(3,6)**2 + dd(3,7)**2 + dd(2,7)**2 + dd(2,6)**2) &
         & *16*pi**2*0.5_dp/(-4*dd(2,1)*dd(3,1)*dd(5,5))
    as = alphasPDF(Q)
    call pdfs(eta, Q, fpdf)
    w = jac*dphi/(2*eta*s)/(16*pi**2)*gev2pb*(as/(2*pi))
    born(1) = QQ*4.0_dp/9*(fpdf(2) + fpdf(-2) + fpdf(4) + fpdf(-4))
    born(2) = QQ*1.0_dp/9*(fpdf(1) + fpdf(-1) + fpdf(3) + fpdf(-3) + fpdf(5) + fpdf(-5))
    born(3) = GQ*11.0_dp/9*fpdf(0)
    select case (trim(bpart))
    case ('b0')
       val = sum(born)
    case ('b1', 'b2')
       call lp21_born(P, eta, ntc, tcs, f0, c1, c2)
       ! unit-charge matrix element times the cumulant coefficients (which
       ! carry the e^2-weighted PDFs themselves; with a PDF mask a class can
       ! have f0 = 0 but c1, c2 /= 0 through the off-diagonal beam functions)
       mu(1:2) = QQ; mu(3) = GQ*11.0_dp/9
       val = 0
       do c = 1, 3
          if (trim(bpart) == 'b1') then
             val = val + mu(c)*c1(:,c)*(as/(2*pi))
          else
             val = val + mu(c)*c2(:,c)*(as/(2*pi))**2
          endif
       enddo
    end select
    val = val*w
    if (val(ntc) /= val(ntc)) then
       res = 0; return
    endif
    do b = 1, nob
       if (mode == 3) then
          if (dp2b(b) /= 0) hcacc(1 + ntc*(b - 1):ntc*b) = hcacc(1 + ntc*(b - 1):ntc*b) + val*wgt*dp2b(b)
          cycle
       elseif (mode == 2) then
          if (.not. inb(b)) cycle
       elseif (tz < zlo(b) .or. tz >= zhi(b)) then
          cycle
       endif
       hcacc(1 + ntc*(b - 1):ntc*b) = hcacc(1 + ntc*(b - 1):ntc*b) + val*wgt
    enddo
    res = val(ntc)
    ! mode 3: the VEGAS target is the >= 1 jet row (the total row vanishes in P2B)
    if (mode == 3) res = val(ntc)*dp2b(2)
  contains
    real(dp) function dd(i, j)
      integer, intent(in) :: i, j
      dd = mdot(P(:,i), P(:,j))
    end function dd
  end function b21_eval
end module sliced21_mod

program sliced21
  use sliced21_mod
  use mod_slicing_scet, only: soft_table_init, scet_set_colour
  implicit none
  character(256) :: arg, softtab, gprefix
  integer :: ncall, itmx, seed, nseed, i
  integer, allocatable :: sd(:)
  integer :: lpmask
  common/lpmask/lpmask
  real(dp) :: avg, err, chi2, hb, bt(9,3), bd(9,3), xi, rr, em(3)
  integer :: k
  call get_command_argument(1, bpart)
  if (trim(bpart) == 'mktab' .or. trim(bpart) == 'tabchkg') then
     call initPDFSetByName('NNPDF30_nlo_as_0118')
     call initPDF(0)
     call scet_set_colour(4.0_dp/3, 3.0_dp, 0.5_dp)
     if (trim(bpart) == 'mktab') then
        call get_command_argument(2, arg); read(arg, *) i
        call get_command_argument(3, arg); read(arg, *) k
        call get_command_argument(4, softtab)
        rr = gdl; hb = gh
        if (command_argument_count() > 4) then
           call get_command_argument(5, arg); read(arg, *) rr
        endif
        if (command_argument_count() > 5) then
           call get_command_argument(6, arg); read(arg, *) hb
        endif
        call lp21_grid_build(gqlo, gqhi, rr, gximin, hb, i, k, trim(softtab))
     else
        call get_command_argument(2, softtab)
        call lp21_grid_load(gqlo, gqhi, trim(softtab))
        em = 0
        do k = 1, 300
           call random_number(rr); xi = gximin*(0.95_dp/gximin)**rr
           call random_number(rr); call lp21_setq(gqlo*(gqhi/gqlo)**rr)
           call beam_at(xi, bt); call lp21_beam_direct(xi, bd)
           do i = 1, 3
              em(i) = max(em(i), maxval(abs(bt(:,i) - bd(:,i)))/maxval(abs(bd(:,i))))
           enddo
        enddo
        write(*,'(a,3es10.2)') ' beam grid: max deviation / max coefficient (up, down, gluon)', em
     endif
     stop
  endif
  call get_command_argument(2, arg); read(arg, *) ncall
  call get_command_argument(3, arg); read(arg, *) itmx
  call get_command_argument(4, arg); read(arg, *) seed
  call get_command_argument(5, arg)
  if (trim(arg) == 'zeus' .or. trim(arg) == 'p2b') then
     mode = merge(2, 3, trim(arg) == 'zeus')
     call get_command_argument(6, gprefix)
  else
     read(arg, *) xfix
     call get_command_argument(6, arg); read(arg, *) Q2fix
  endif
  softtab = ''; hb = 0.02_dp
  if (command_argument_count() > 6) call get_command_argument(7, softtab)
  if (softtab == '-') softtab = ''
  if (command_argument_count() > 7 .and. mode /= 2) then
     call get_command_argument(8, arg); read(arg, *) hb
  endif
  if (command_argument_count() > 8) then
     ! incoming-parton mask on the PDFs everywhere (Born luminosity of
     ! nlo31_mod and the beam functions of lp21): 1 quarks only, 2 gluon only
     call get_command_argument(9, arg); read(arg, *) pdfmask
     lpmask = pdfmask
  endif
  if (mode >= 2) then
     nob = merge(15, 16, mode == 2); nv = ntc*nob; iv = merge(ntc, 2*ntc, mode == 2)
     q2lo = gqlo**2; q2hi = gqhi**2; ylo = 0.2_dp; yhi = 0.6_dp
  else
     mode = 1; nob = nzb; nv = ntc*nob; iv = ntc + ntc*(nzb - 1)
  endif
  call random_seed(size=nseed); allocate(sd(nseed))
  sd = [(1000003*seed + 7919*i, i = 1, nseed)]
  call random_seed(put=sd)
  call initPDFSetByName('NNPDF30_nlo_as_0118')
  call initPDF(0)
  call scet_set_colour(4.0_dp/3, 3.0_dp, 0.5_dp)
  if (softtab /= '') call soft_table_init(trim(softtab), 0.025_dp, 0.05_dp)
  if (mode >= 2) then
     if (trim(bpart) /= 'b0') call lp21_grid_load(gqlo, gqhi, trim(gprefix))
  elseif (trim(bpart) /= 'b0') then
     call lp21_init(sqrt(Q2fix), xfix, hb)
  endif
  if (trim(bpart) == 'tabchk') then
     ! largest deviation per class, relative to the largest coefficient of the class
     em = 0
     do k = 1, 200
        call random_number(rr); xi = xfix*(0.95_dp/xfix)**rr
        call beam_at(xi, bt); call lp21_beam_direct(xi, bd)
        do i = 1, 3
           em(i) = max(em(i), maxval(abs(bt(:,i) - bd(:,i)))/maxval(abs(bd(:,i))))
        enddo
        if (maxval(abs(bt(:,3) - bd(:,3)))/maxval(abs(bd(:,3))) > 3e-4_dp) &
             write(*,'(a,es12.4,a,9es11.3)') ' xi', xi, '  gluon table-direct / max:', (bt(:,3) - bd(:,3))/maxval(abs(bd(:,3)))
     enddo
     write(*,'(a,f7.4,a,3es10.2)') ' beam table h', hb, ': max deviation (up, down, gluon)', em
     stop
  endif
  write(*,'(a,a,a,i10,a,i4,a,i6,a,f9.6,a,f10.2)') ' sliced21 part ', trim(bpart), ' ncall', ncall, ' itmx', itmx, &
       & ' seed', seed, ' x', xfix, ' Q2', Q2fix
  if (trim(bpart) == 'c1') then
     ! correlated b1 + lo (mode >= 2 with psmc): nlo31's lo on the same points
     if (mode < 2) stop 'c1: mode 2 or 3 (zeus | p2b) only'
     part = 'lo'; bpart = 'b1'; usepsmc = .true.
     call get_environment_variable('C1_NEMIT', arg)
     if (len_trim(arg) > 0) read(arg, *) nemit
     call get_environment_variable('C1_VBAL', arg)
     if (len_trim(arg) > 0) read(arg, *) vbal
     write(*,'(a,i4,a,i2)') ' correlated sampling c1: emissions per Born', nemit, '  balanced VEGAS target', vbal
     call setup_flavours()
     call vegas(9, ncall, itmx, avg, err, chi2, c1_part)
     bpart = 'c1'
  else
     call vegas(5, ncall, itmx, avg, err, chi2, b21_part)
  endif
  write(*,'(a,a,a,es16.8,a,es12.4,a,f8.3)') ' RESULT ', trim(bpart), ' dsigma/dx dQ2 [pb/GeV2] (smallest tau_cut, all bins) = ', &
       & avg, ' +- ', err, '   chi2/it', chi2
  if (mode == 3) then
     write(*,'(a)') ' LCELLS sigma per bin [pb] below the cut, O(event) - O(Born): lab11 bin, then tau_cut columns'
     write(*,'(a,10es16.8)') ' tau_cut     ', tcs
     do i = 1, 16
        write(*,'(a,i3,10es16.8)') ' LCELL ', i, hc(1 + ntc*(i - 1):ntc*i)
     enddo
  elseif (mode == 2) then
     write(*,'(a)') ' ZCELLS sigma per bin [pb] below the cut: observable bin, then tau_cut columns'
     write(*,'(a,10es12.3)') ' tau_cut     ', tcs
     write(*,'(a,10es16.8)') ' ZCELL total     0     0', hc(1:ntc)
     do i = 1, 6
        write(*,'(a,2f8.0,10es16.8)') ' ZCELL q2   ', zq2e(i-1), zq2e(i), hc(1 + ntc*i:ntc*(i + 1))
     enddo
     do i = 1, 4
        write(*,'(a,2f8.0,10es16.8)') ' ZCELL ptavg', zpte(i-1), zpte(i), hc(1 + ntc*(6 + i):ntc*(7 + i))
     enddo
     do i = 1, 4
        write(*,'(a,2f8.0,10es16.8)') ' ZCELL m12  ', zmje(i-1), zmje(i), hc(1 + ntc*(10 + i):ntc*(11 + i))
     enddo
  else
     write(*,'(a)') ' CELLS dsigma/dx dQ2 [pb/GeV2] below the cut: tau_zQ bin, then tau_cut columns'
     write(*,'(a,10es12.3)') ' tau_cut     ', tcs
     do i = 1, nzb
        write(*,'(a,2f5.2,10es16.8)') ' CELL ', zlo(i), zhi(i), hc(1 + ntc*(i - 1):ntc*i)
     enddo
  endif
end program sliced21
