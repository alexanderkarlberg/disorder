! W exchange (9 Oct 2026): the flavour sums of nlo31 from representative
! channels (generation 1: line q1 -> q1o, W pair of generation 2 (q2o, -q2);
! unit CKM in (u,d), (c,s), no coupling for b) against brute-force sums over
! every final-state flavour multiset (1/prod n_i! for identical partons),
! both summed over the permutations of the outgoing momenta, for the 3+1
! (me31) and 4+1 (me41) trees and each incoming parton -5..5. W- and W+ from
! EW31 = "2 lepton" (default "2 0").
program harness_ccsum
  use me31
  use me41
  use ew31
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
  real(dp) :: P3(4,7), P4(4,8), r(12), b3(-5:5), f3(-5:5), b4(-5:5), f4(-5:5), worst
  integer :: ipt, f, q1, q1o, q2, q2o
  character(32) :: arg
  ew31_mode = 2; ew31_lepton = 0
  call get_environment_variable('EW31', arg)
  if (len_trim(arg) > 0) read(arg, *) ew31_mode, ew31_lepton
  q1 = merge(2, 1, ew31_lepton == 0 .or. ew31_lepton == 3)
  q1o = ew31_out(q1); q2 = q1 + 2; q2o = q1o + 2
  worst = 0
  do ipt = 1, 3
     call random_number(r)
     call point(r, 3, P3, P4)
     call random_number(r)
     call point(r, 4, P3, P4)
     do f = -5, 5
        b3(f) = brute3(f); f3(f) = form3(f)
        b4(f) = brute4(f); f4(f) = form4(f)
     enddo
     write(*,'(a,i2,a,11f9.5)') ' point', ipt, '  3+1 formula/brute', merge(f3/b3, 1.0_dp, b3 /= 0)
     write(*,'(a,i2,a,11f9.5)') ' point', ipt, '  4+1 formula/brute', merge(f4/b4, 1.0_dp, b4 /= 0)
     worst = max(worst, maxval(abs(f3 - b3)/max(abs(b3), tiny(1.0_dp))), maxval(abs(f4 - b4)/max(abs(b4), tiny(1.0_dp))))
  enddo
  write(*,'(a,es10.2)') ' largest |formula/brute - 1|:', worst
  if (worst > 1d-12) stop 1

contains

  ! 3+1: me31 summed over the 6 permutations of the outgoing momenta
  real(dp) function m3(fl)
    integer, intent(in) :: fl(4)
    integer, parameter :: perm(3,6) = reshape([2,3,4, 2,4,3, 3,2,4, 3,4,2, 4,2,3, 4,3,2], [3,6])
    real(dp) :: PP(4,7), x
    integer :: ip
    m3 = 0
    do ip = 1, 6
       PP = P3; PP(:,2:4) = P3(:,perm(:,ip))
       call me31_tree(PP, fl, x); m3 = m3 + x
    enddo
  end function m3

  real(dp) function m4(fl)
    integer, intent(in) :: fl(5)
    real(dp) :: PP(4,8), x
    integer :: a, b, c, d
    m4 = 0
    do a = 2, 5; do b = 2, 5; do c = 2, 5; do d = 2, 5
       if (a == b .or. a == c .or. a == d .or. b == c .or. b == d .or. c == d) cycle
       PP = P4; PP(:,2) = P4(:,a); PP(:,3) = P4(:,b); PP(:,4) = P4(:,c); PP(:,5) = P4(:,d)
       call me41_tree(PP, fl, x); m4 = m4 + x
    enddo; enddo; enddo; enddo
  end function m4

  ! the representative-channel sums (nlo31's weights)
  real(dp) function form3(f)
    integer, intent(in) :: f
    form3 = 0
    if (f == 0) then
       form3 = 2*m3([0, q1o, -q1, 0])
    elseif (f == q1 .or. f == q2) then
       form3 = 0.5_dp*m3([q1, q1o, 0, 0]) + 3*m3([q1, q1o, 5, -5]) + 0.5_dp*m3([q1, q1o, q1o, -q1o]) &
            & + m3([q1, q1o, q1, -q1]) + m3([q1, q1, q2o, -q2])
    elseif (f == q1o .or. f == q2o) then
       form3 = 0.5_dp*m3([q1o, q1o, q1o, -q1]) + m3([q1o, q1o, q2o, -q2])
    elseif (f == 5) then
       form3 = 2*m3([q1o, q1o, q2o, -q2])
    elseif (f == -q1o .or. f == -q2o) then
       form3 = 0.5_dp*m3([-q1o, -q1, 0, 0]) + 3*m3([-q1o, -q1, -5, 5]) + 0.5_dp*m3([-q1o, -q1, -q1, q1]) &
            & + m3([-q1o, -q1, -q1o, q1o]) + m3([-q1o, -q1o, -q2, q2o])
    elseif (f == -q1 .or. f == -q2) then
       form3 = 0.5_dp*m3([-q1, -q1, -q1, q1o]) + m3([-q1, -q1, -q2, q2o])
    elseif (f == -5) then
       form3 = 2*m3([-q1, -q1, -q2, q2o])
    endif
  end function form3

  real(dp) function form4(f)
    integer, intent(in) :: f
    form4 = 0
    if (f == 0) then
       form4 = 2*(0.5_dp*m4([0, q1o, -q1, 0, 0]) + 3*m4([0, q1o, -q1, 5, -5]) &
            & + 0.5_dp*m4([0, q1o, -q1, q1o, -q1o]) + 0.5_dp*m4([0, q1o, -q1, q1, -q1]))
    elseif (f == q1 .or. f == q2) then
       form4 = m4([q1, q1o, 0, 0, 0])/6 + 3*m4([q1, q1o, 5, -5, 0]) + 0.5_dp*m4([q1, q1o, q1o, -q1o, 0]) &
            & + m4([q1, q1o, q1, -q1, 0]) + m4([q1, q1, q2o, -q2, 0])
    elseif (f == q1o .or. f == q2o) then
       form4 = 0.5_dp*m4([q1o, q1o, q1o, -q1, 0]) + m4([q1o, q1o, q2o, -q2, 0])
    elseif (f == 5) then
       form4 = 2*m4([q1o, q1o, q2o, -q2, 0])
    elseif (f == -q1o .or. f == -q2o) then
       form4 = m4([-q1o, -q1, 0, 0, 0])/6 + 3*m4([-q1o, -q1, -5, 5, 0]) + 0.5_dp*m4([-q1o, -q1, -q1, q1, 0]) &
            & + m4([-q1o, -q1, -q1o, q1o, 0]) + m4([-q1o, -q1o, -q2, q2o, 0])
    elseif (f == -q1 .or. f == -q2) then
       form4 = 0.5_dp*m4([-q1, -q1, -q1, q1o, 0]) + m4([-q1, -q1, -q2, q2o, 0])
    elseif (f == -5) then
       form4 = 2*m4([-q1, -q1, -q2, q2o, 0])
    endif
  end function form4

  ! brute force: every multiset of outgoing flavours (non-decreasing codes),
  ! weight 1/prod n_i!
  real(dp) function brute3(f)
    integer, intent(in) :: f
    integer :: a, b, c
    brute3 = 0
    do a = -5, 5; do b = a, 5; do c = b, 5
       brute3 = brute3 + m3([f, a, b, c])/sym([a, b, c])
    enddo; enddo; enddo
  end function brute3

  real(dp) function brute4(f)
    integer, intent(in) :: f
    integer :: a, b, c, d
    brute4 = 0
    do a = -5, 5; do b = a, 5; do c = b, 5; do d = c, 5
       ! charge (quark number) conservation per generation is left to me41
       brute4 = brute4 + m4([f, a, b, c, d])/sym([a, b, c, d])
    enddo; enddo; enddo; enddo
  end function brute4

  real(dp) function sym(c)
    integer, intent(in) :: c(:)
    integer :: i, n
    sym = 1
    n = 1
    do i = 2, size(c)
       if (c(i) == c(i-1)) then
          n = n + 1; sym = sym*n
       else
          n = 1
       endif
    enddo
  end function sym

  ! a DIS 3+1 or 4+1 point (Breit-like frame from flat sequential decays)
  subroutine point(r, n, P3, P4)
    real(dp), intent(in) :: r(12)
    integer, intent(in) :: n
    real(dp), intent(inout) :: P3(4,7), P4(4,8)
    real(dp) :: E1, y, Q, W, ks(4,4)
    integer :: i
    Q = 30; y = 0.3_dp + 0.4_dp*r(1)
    W = Q*(1.5_dp + 3*r(2))
    E1 = (W**2 + Q**2)/(2*Q)
    call rambo(n, W, r(3:), ks)
    if (n == 3) then
       P3 = 0
       P3(:,1) = [0.0_dp, 0.0_dp, E1, E1]
       P3(:,5) = [0.0_dp, 0.0_dp, -Q, 0.0_dp]
       P3(:,6) = [Q/(2*y)*2*sqrt(1 - y), 0.0_dp, -Q/2, Q/(2*y)*(2 - y)]
       P3(:,7) = [Q/(2*y)*2*sqrt(1 - y), 0.0_dp, Q/2, Q/(2*y)*(2 - y)]
       do i = 1, 3
          P3(:,1 + i) = boost(ks(:,i), P3(:,1) + P3(:,5))
       enddo
    else
       P4 = 0
       P4(:,1) = [0.0_dp, 0.0_dp, E1, E1]
       P4(:,6) = [0.0_dp, 0.0_dp, -Q, 0.0_dp]
       P4(:,7) = [Q/(2*y)*2*sqrt(1 - y), 0.0_dp, -Q/2, Q/(2*y)*(2 - y)]
       P4(:,8) = [Q/(2*y)*2*sqrt(1 - y), 0.0_dp, Q/2, Q/(2*y)*(2 - y)]
       do i = 1, 4
          P4(:,1 + i) = boost(ks(:,i), P4(:,1) + P4(:,6))
       enddo
    endif
  end subroutine point

  ! massless n-body momenta in the rest frame of mass W (RAMBO-like: random
  ! directions and energies, rescaled; the weight is irrelevant here)
  subroutine rambo(n, W, r, k)
    integer, intent(in) :: n
    real(dp), intent(in) :: W, r(:)
    real(dp), intent(out) :: k(4,4)
    real(dp) :: q(4,4), Qs(4), M, b(3), g, a, x, c, s, ph, e, bq
    integer :: i
    k = 0
    do i = 1, n
       c = 2*r(3*i - 2) - 1; s = sqrt(1 - c*c); ph = 2*pi*r(3*i - 1)
       e = -log(max(r(3*i)*0.999_dp + 1d-3*r(3*i - 2), 1d-12))
       q(:,i) = e*[s*cos(ph), s*sin(ph), c, 1.0_dp]
    enddo
    Qs = sum(q(:,1:n), dim=2)
    M = sqrt(Qs(4)**2 - sum(Qs(1:3)**2))
    b = -Qs(1:3)/M; g = Qs(4)/M; a = 1/(1 + g); x = W/M
    do i = 1, n
       bq = dot_product(b, q(1:3,i))
       k(4,i) = x*(g*q(4,i) + bq)
       k(1:3,i) = x*(q(1:3,i) + b*q(4,i) + a*bq*b)
    enddo
  end subroutine rambo

  ! boost k from the rest frame of P (timelike) to the frame of P
  function boost(k, P) result(kb)
    real(dp), intent(in) :: k(4), P(4)
    real(dp) :: kb(4), M, bp, f
    M = sqrt(P(4)**2 - sum(P(1:3)**2))
    bp = dot_product(k(1:3), P(1:3))
    kb(4) = (P(4)*k(4) + bp)/M
    f = (bp/(P(4) + M) + k(4))/M
    kb(1:3) = k(1:3) + f*P(1:3)
  end function boost
end program harness_ccsum
