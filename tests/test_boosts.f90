!----------------------------------------------------------------------
! Unit tests of the lab <-> Breit frame transformations mlab2breit and
! mbreit2lab (src/mod_phase_space.f90).
program test_boosts
  use types, only: dp
  use mod_phase_space
  use test_utils
  implicit none
  integer, parameter :: npart = 5
  real(dp) :: q(0:3), p(0:3,npart), pb(0:3,npart), pl(0:3,npart)
  real(dp) :: hadron(0:3,1), hb(0:3,1), qb(0:3,1), qq(0:3,1)
  real(dp) :: Qval, xB, yB, phi, r(3), tol, trans(0:3,4)
  real(dp), parameter :: El = 27.5_dp, Eh = 820.0_dp
  integer :: itry, i, j
  character(len=80) :: tag
  logical :: isPlus

  call random_seed(put=[(12345 + i, i = 1, 64)])

  do itry = 1, 20
     isPlus = mod(itry,2) == 1
     write(tag,'(a,i0,a,l1)') 'config ', itry, ' isPlus=', isPlus

     ! q = k - k' for HERA beams (lepton along -z, 27.5 GeV on 820 GeV),
     ! random x, y and lepton azimuth; mirrored in z for isPlus=.false.
     call random_number(r)
     xB = 1e-4_dp * (0.9_dp/1e-4_dp)**r(1)
     yB = 0.01_dp + 0.98_dp * r(2)
     phi = 6.283185307179586_dp * r(3)
     Qval = sqrt(xB * yB * 4.0_dp * El * Eh)
     q = [yB * (El - xB * Eh), -Qval * sqrt(1 - yB) * cos(phi), &
          & -Qval * sqrt(1 - yB) * sin(phi), -yB * (El + xB * Eh)]
     hadron(:,1) = [1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
     if (.not. isPlus) then
        q(3) = -q(3)
        hadron(3,1) = -1.0_dp
     endif
     ! Rounding errors grow with the square of the largest factor by
     ! which the transformation enlarges a momentum component.
     call mlab2breit(4, q, reshape([1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1]*1.0_dp, [4,4]), &
          & trans, isPlus)
     tol = 1e-15_dp * (1.0_dp + maxval(abs(trans)))**2

     do j = 1, npart
        call random_number(r)
        p(1:3,j) = 100.0_dp * (r - 0.5_dp)
        p(0,j) = sqrt(sum(p(1:3,j)**2)) ! massless
     enddo

     ! Breit -> lab undoes lab -> Breit
     call mlab2breit(npart, q, p, pb, isPlus)
     call mbreit2lab(npart, q, pb, pl, isPlus)
     do j = 1, npart
        call check_vec_close(trim(tag)//': lab->Breit->lab round trip', pl(:,j), p(:,j), tol)
     enddo
     call mbreit2lab(npart, q, p, pl, isPlus)
     call mlab2breit(npart, q, pl, pb, isPlus)
     do j = 1, npart
        call check_vec_close(trim(tag)//': Breit->lab->Breit round trip', pb(:,j), p(:,j), tol)
     enddo

     ! Lorentz transformation: all scalar products preserved
     call mlab2breit(npart, q, p, pb, isPlus)
     do i = 1, npart
        do j = i, npart
           ! (absolute tolerance for the vanishing p_i^2 of massless momenta)
           call check_close(trim(tag)//': invariant p_i.p_j', minkowski(pb(:,i),pb(:,j)), &
                & minkowski(p(:,i),p(:,j)), tol, &
                & 1e-14_dp*max(p(0,i)*p(0,j), pb(0,i)*pb(0,j)))
        enddo
     enddo

     ! Defining properties of the Breit frame: q = (0,0,0,-Q) and the
     ! incoming parton x*P = (Q/2)(1,0,0,1), with x = Q^2/(2 P.q).
     call check_close(trim(tag)//': Q^2 = x y s', -minkowski(q,q), Qval**2, 1e-12_dp)
     qq(:,1) = q
     call mlab2breit(1, q, qq, qb, isPlus)
     call check_vec_close(trim(tag)//': q in Breit frame', qb(:,1), &
          & [0.0_dp, 0.0_dp, 0.0_dp, merge(-Qval, Qval, isPlus)], tol)
     xB = Qval**2 / (2.0_dp * minkowski(hadron(:,1), q))
     call mlab2breit(1, q, xB*hadron, hb, isPlus)
     call check_vec_close(trim(tag)//': incoming parton in Breit frame', hb(:,1), &
          & 0.5_dp*Qval*[1.0_dp, 0.0_dp, 0.0_dp, merge(1.0_dp, -1.0_dp, isPlus)], tol)
  enddo

  ! isPlus = .false. is the mirror image (z -> -z) of isPlus = .true.
  q = [3.0_dp, -4.0_dp, 1.5_dp, -12.0_dp]
  do j = 1, npart
     call random_number(r)
     p(1:3,j) = 40.0_dp * (r - 0.5_dp)
     p(0,j) = sqrt(sum(p(1:3,j)**2))
  enddo
  call mlab2breit(npart, q, p, pb, .true.)
  call mlab2breit(npart, mirror(q), mirror_all(p), pl, .false.)
  do j = 1, npart
     call check_vec_close('isPlus=.false. is the z-mirror of isPlus=.true.', &
          & pl(:,j), mirror(pb(:,j)), 1e-12_dp)
  enddo

  call finish_tests()

contains

  pure function minkowski(a, b)
    real(dp), intent(in) :: a(0:3), b(0:3)
    real(dp) :: minkowski
    minkowski = a(0)*b(0) - a(1)*b(1) - a(2)*b(2) - a(3)*b(3)
  end function minkowski

  pure function mirror(a)
    real(dp), intent(in) :: a(0:3)
    real(dp) :: mirror(0:3)
    mirror = [a(0), a(1), a(2), -a(3)]
  end function mirror

  pure function mirror_all(a)
    real(dp), intent(in) :: a(0:,:)
    real(dp) :: mirror_all(0:3,size(a,2))
    integer :: k
    do k = 1, size(a,2)
       mirror_all(:,k) = mirror(a(:,k))
    enddo
  end function mirror_all

end program test_boosts
