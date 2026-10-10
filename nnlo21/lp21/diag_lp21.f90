! diagnostics of the O(alpha_s^2) cumulant: polynomial fit in L = ln(tau_cut)
! of c1 and c2 at one Born point (class 1, up quarks)
program diag_lp21
  use lp21
  implicit none
  integer, parameter :: dp = kind(1.0d0), n = 9
  real(dp) :: P(4,7), tcs(n), f0(nbc), c1(n,nbc), c2(n,nbc), L(n), A(n,5), b(n), cf1(5), cf2(5), wk(2,nbc), nfk(2,nbc)
  integer :: i
  call InitPDFsetByName('NNPDF30_nlo_as_0118'); call InitPDF(0)
  call lp21_init(20.0_dp, 0.005_dp, 0.0_dp)
  P = 0
  P(:,1) = [0.0_dp, 0.0_dp, 25.0_dp, 25.0_dp]; P(:,5) = [0.0_dp, 0.0_dp, -20.0_dp, 0.0_dp]
  P(:,2) = [8.0_dp, 3.0_dp, 0.0_dp, 0.0_dp]; P(4,2) = sqrt(sum(P(1:3,2)**2))
  P(:,3) = P(:,1) + P(:,5) - P(:,2)
  P(:,6) = 10*[2*sqrt(0.5_dp)/0.5_dp, 0.0_dp, -1.0_dp, 1.5_dp/0.5_dp]; P(:,7) = P(:,6) - P(:,5)
  print *, 'p3^2 =', P(4,3)**2 - sum(P(1:3,3)**2)
  do i = 1, n
     tcs(i) = 10.0_dp**(-1 - 0.5_dp*i)
  enddo
  wk = 1; nfk = 0.5_dp
  call lp21_born(P, 0.4_dp, n, tcs, wk, nfk, f0, c1, c2)
  L = log(tcs)
  do i = 1, 5
     A(:,i) = L**(i - 1)
  enddo
  call lsq(A(:,1:3), c1(:,1)/f0(1), cf1(1:3), 3)
  call lsq(A, c2(:,1)/f0(1), cf2, 5)
  print '(a,3es14.6)', ' c1: L^0, L^1, L^2  ', cf1(1:3)
  print '(a,5es14.6)', ' c2: L^0 .. L^4     ', cf2
  print '(a,es14.6)', ' (c1_L2)^2/2        ', cf1(3)**2/2
contains
  subroutine lsq(M, y, c, k)
    integer, intent(in) :: k
    real(dp), intent(in) :: M(:,:), y(:)
    real(dp), intent(out) :: c(k)
    real(dp) :: N(k,k), r(k), piv
    integer :: i, j
    N = matmul(transpose(M(:,1:k)), M(:,1:k)); r = matmul(transpose(M(:,1:k)), y)
    do i = 1, k
       piv = N(i,i)
       N(i,:) = N(i,:)/piv; r(i) = r(i)/piv
       do j = 1, k
          if (j /= i) then
             r(j) = r(j) - N(j,i)*r(i); N(j,:) = N(j,:) - N(j,i)*N(i,:)
          endif
       enddo
    enddo
    c = r
  end subroutine lsq
end program diag_lp21
