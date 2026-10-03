!-----------------------------------------------------------------------
! test_lp21: the O(alpha_s) part of lp21 (MCFM pieces) against the NLO
! tau_2 cumulant of the slicing code (mod_slicing_scet: beam_coeffs,
! jet_cum, soft_geom/soft_from_geom; the hard part is hard21's h(1) in both,
! which equals the DISENT-based hard function, nnlo21/tests). Also MCFM's
! one-loop soft function alone against soft_from_geom, and the O(alpha_s^2)
! coefficients printed for inspection.
!-----------------------------------------------------------------------
program test_lp21
  use lp21
  use hard21, only: hard21_eval
  use mod_slicing_scet, only: beam_coeffs, jet_cum, soft_geom, soft_from_geom, scet_set_colour, &
       & scet_set_nf, CF, CA
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
  integer, parameter :: ntc = 4
  real(dp) :: tcs(ntc) = [1e-2_dp, 1e-3_dp, 1e-4_dp, 1e-5_dp]
  real(dp) :: P(4,7), Q, x, xp, th, ph, r1, xi, f0(3), c1(ntc,3), c2(ntc,3), h(0:2)
  real(dp) :: c0b(-6:6), c1b(-6:6), c2b(-6:6), xf(-6:6), pb(4,3), y, ch, sh, nh(3,3), g(3,3), ls(3,3)
  real(dp) :: tq(3,3), tg(3,3), ref, lb, beam, worst(3), Ea, E2, E3, wsum, wbeam
  real(dp) :: e2c(-5:5)
  integer :: ip, it, c, k, n, i
  call scet_set_colour(4.0_dp/3, 3.0_dp, 0.5_dp)
  call InitPDFsetByName('NNPDF30_nlo_as_0118'); call InitPDF(0)
  Q = 20; x = 0.01_dp
  call lp21_init(Q, 0.5_dp*x, 0.0_dp)        ! no table: direct beam integrals
  e2c = [1, 4, 1, 4, 1, 0, 1, 4, 1, 4, 1]/9.0_dp
  call ttm([CF, CF, CA], tq); call ttm([CA, CF, CF], tg)
  worst = 0
  call random_seed(put=[(4321 + 11*n, n = 1, 33)])
  do ip = 1, 6
     call random_number(r1); xp = 0.1_dp + 0.85_dp*r1
     call random_number(r1); th = acos(-0.9_dp + 1.8_dp*r1)
     call random_number(r1); ph = 2*pi*r1
     call breit_born(Q, 0.5_dp, xp, th, ph, P)
     xi = x/xp
     call lp21_born(P, xi, ntc, tcs, f0, c1, c2)
     ! reference pieces in the jets' frame
     y = -atanh((P(3,2) + P(3,3))/(P(4,2) + P(4,3))); ch = cosh(y); sh = sinh(y)
     do k = 1, 3
        pb(:,k) = P(:,k); pb(4,k) = ch*P(4,k) + sh*P(3,k); pb(3,k) = sh*P(4,k) + ch*P(3,k)
        nh(:,k) = pb(1:3,k)/sqrt(sum(pb(1:3,k)**2))
     enddo
     Ea = pb(4,1); E2 = pb(4,2); E3 = pb(4,3)
     call soft_geom(3, nh, g, ls)
     call beam_coeffs(xi, Q, c0b, c1b, c2b)
     call evolvePDF(xi, Q, xf)
     do c = 1, 3
        if (c <= 2) then
           call hard21_eval(P, 1, merge(2.0_dp/3, -1.0_dp/3, c == 1), h)
        else
           call hard21_eval(P, 2, 11.0_dp/3, h)
        endif
        do it = 1, ntc
           lb = log(2*Ea*tcs(it)/Q)
           wsum = 0; wbeam = 0
           do i = -5, 5
              if ((c == 1 .and. (abs(i) == 2 .or. abs(i) == 4)) .or. &
                  (c == 2 .and. (abs(i) == 1 .or. abs(i) == 3 .or. abs(i) == 5)) .or. (c == 3 .and. i == 0)) then
                 wsum = wsum + merge(1.0_dp, e2c(i), c == 3)*xf(i)
                 wbeam = wbeam + merge(1.0_dp, e2c(i), c == 3)*(c0b(i) + c1b(i)*lb + c2b(i)*lb*lb)
              endif
           enddo
           if (c <= 2) then
              ref = h(1) + jet_cum(.false., 2*E2*tcs(it)/Q) + jet_cum(.true., 2*E3*tcs(it)/Q) &
                   + soft_from_geom(3, g, ls, [CF, CF, CA], tq, tcs(it)) + wbeam/wsum
           else
              ref = h(1) + jet_cum(.false., 2*E2*tcs(it)/Q) + jet_cum(.false., 2*E3*tcs(it)/Q) &
                   + soft_from_geom(3, g, ls, [CA, CF, CF], tg, tcs(it)) + wbeam/wsum
           endif
           worst(c) = max(worst(c), abs(c1(it,c)/f0(c) - ref)/max(1.0_dp, abs(ref)))
           if (ip <= 2) print '(a,i2,a,i2,a,es9.1,a,2es16.8,a,es14.6)', 'pt', ip, ' class', c, ' tau_cut', tcs(it), &
                '  O(as): lp21, slicing', c1(it,c)/f0(c), ref, '  O(as^2) lp21', c2(it,c)/f0(c)
        enddo
     enddo
  enddo
  print '(a,3es10.2)', 'O(alpha_s): max |lp21 - slicing| (up, down, gluon):', worst
contains
  subroutine ttm(cas, tt)
    real(dp), intent(in) :: cas(3)
    real(dp), intent(out) :: tt(3,3)
    tt = 0
    tt(1,2) = 0.5_dp*(cas(3) - cas(1) - cas(2)); tt(2,1) = tt(1,2)
    tt(1,3) = 0.5_dp*(cas(2) - cas(1) - cas(3)); tt(3,1) = tt(1,3)
    tt(2,3) = 0.5_dp*(cas(1) - cas(2) - cas(3)); tt(3,2) = tt(2,3)
  end subroutine ttm
  subroutine breit_born(Q, y, xp, th, ph, P)
    real(dp), intent(in) :: Q, y, xp, th, ph
    real(dp), intent(out) :: P(4,7)
    real(dp) :: Ea, W(4), sh2, E, k(4), gam, bz
    integer :: n
    P = 0
    Ea = Q/(2*xp)
    P(:,1) = [0.0_dp, 0.0_dp, Ea, Ea]
    P(:,5) = [0.0_dp, 0.0_dp, -Q, 0.0_dp]
    W = P(:,1) + P(:,5)
    sh2 = W(4)**2 - W(3)**2
    E = sqrt(sh2)/2
    k = [E*sin(th)*cos(ph), E*sin(th)*sin(ph), E*cos(th), E]
    bz = W(3)/W(4); gam = 1/sqrt(1 - bz**2)
    P(:,2) = k; P(:,3) = [-k(1), -k(2), -k(3), k(4)]
    do n = 2, 3
       E = P(4,n)
       P(4,n) = gam*(E + bz*P(3,n)); P(3,n) = gam*(P(3,n) + bz*E)
    enddo
    P(:,6) = Q/2*[2*sqrt(1 - y)/y, 0.0_dp, -1.0_dp, (2 - y)/y]
    P(:,7) = P(:,6) - P(:,5)
  end subroutine breit_born
end program test_lp21
