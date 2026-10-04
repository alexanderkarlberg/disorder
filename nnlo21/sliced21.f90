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
!-----------------------------------------------------------------------
module sliced21_mod
  use nlo31_mod
  use lp21
  implicit none
  character(8) :: bpart

contains

  real(dp) function b21_part(r, wgt) result(res)
    real(dp), intent(in) :: r(:), wgt
    real(dp), external :: alphasPDF
    real(dp) :: Q2, y, xB, eta, jac, Pk(4,6), P(4,7), dphi, as, fpdf(-5:5), w, Q, tz
    real(dp) :: QQ, GQ, born(3), f0(3), c1(ntc,3), c2(ntc,3), val(ntc)
    logical :: ok
    integer :: b, k, c
    res = 0
    call lepton(r(1:3), Q2, y, xB, eta, jac, ok)
    if (.not. ok) return
    call breit(2, r(4:5), Q2, y, xB, eta, Pk, dphi)
    Q = sqrt(Q2)
    tz = 1
    do k = 2, 3
       if (Pk(3,k) < 0) tz = tz + 2*Pk(3,k)/Q
    enddo
    if (tz < zlo(nzb) .or. tz >= zhi(nzb)) return
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
       val = 0
       do c = 1, 3
          if (f0(c) == 0) cycle
          if (trim(bpart) == 'b1') then
             val = val + born(c)*c1(:,c)/f0(c)*(as/(2*pi))
          else
             val = val + born(c)*c2(:,c)/f0(c)*(as/(2*pi))**2
          endif
       enddo
    end select
    val = val*w
    do b = 1, nzb
       if (tz >= zlo(b) .and. tz < zhi(b)) hcacc(1 + ntc*(b - 1):ntc*b) = hcacc(1 + ntc*(b - 1):ntc*b) + val*wgt
    enddo
    res = val(ntc)
  contains
    real(dp) function dd(i, j)
      integer, intent(in) :: i, j
      dd = mdot(P(:,i), P(:,j))
    end function dd
  end function b21_part
end module sliced21_mod

program sliced21
  use sliced21_mod
  use mod_slicing_scet, only: soft_table_init, scet_set_colour
  implicit none
  character(256) :: arg, softtab
  integer :: ncall, itmx, seed, nseed, i
  integer, allocatable :: sd(:)
  integer :: lpmask
  common/lpmask/lpmask
  real(dp) :: avg, err, chi2, hb, bt(9,3), bd(9,3), xi, rr, em(3)
  integer :: k
  call get_command_argument(1, bpart)
  call get_command_argument(2, arg); read(arg, *) ncall
  call get_command_argument(3, arg); read(arg, *) itmx
  call get_command_argument(4, arg); read(arg, *) seed
  call get_command_argument(5, arg); read(arg, *) xfix
  call get_command_argument(6, arg); read(arg, *) Q2fix
  softtab = ''; hb = 0.02_dp
  if (command_argument_count() > 6) call get_command_argument(7, softtab)
  if (softtab == '-') softtab = ''
  if (command_argument_count() > 7) then
     call get_command_argument(8, arg); read(arg, *) hb
  endif
  if (command_argument_count() > 8) then
     ! incoming-parton mask on the PDFs everywhere (Born luminosity of
     ! nlo31_mod and the beam functions of lp21): 1 quarks only, 2 gluon only
     call get_command_argument(9, arg); read(arg, *) pdfmask
     lpmask = pdfmask
  endif
  mode = 1; nv = ncell; iv = ntc + ntc*(nzb - 1)
  call random_seed(size=nseed); allocate(sd(nseed))
  sd = [(1000003*seed + 7919*i, i = 1, nseed)]
  call random_seed(put=sd)
  call initPDFSetByName('NNPDF30_nlo_as_0118')
  call initPDF(0)
  call scet_set_colour(4.0_dp/3, 3.0_dp, 0.5_dp)
  if (softtab /= '') call soft_table_init(trim(softtab), 0.025_dp, 0.05_dp)
  if (trim(bpart) /= 'b0') call lp21_init(sqrt(Q2fix), xfix, hb)
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
  call vegas(5, ncall, itmx, avg, err, chi2, b21_part)
  write(*,'(a,a,a,es16.8,a,es12.4,a,f8.3)') ' RESULT ', trim(bpart), ' dsigma/dx dQ2 [pb/GeV2] (smallest tau_cut, all bins) = ', &
       & avg, ' +- ', err, '   chi2/it', chi2
  write(*,'(a)') ' CELLS dsigma/dx dQ2 [pb/GeV2] below the cut: tau_zQ bin, then tau_cut columns'
  write(*,'(a,10es12.3)') ' tau_cut     ', tcs
  do i = 1, nzb
     write(*,'(a,2f5.2,10es16.8)') ' CELL ', zlo(i), zhi(i), hc(1 + ntc*(i - 1):ntc*i)
  enddo
end program sliced21
