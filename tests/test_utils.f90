!----------------------------------------------------------------------
! Minimal assertion helpers for the unit tests. Each check prints one
! line; finish_tests() stops with a non-zero exit code if any failed,
! which is what ctest looks at.
module test_utils
  use types, only: dp
  implicit none

  private
  public :: check_close, check_true, check_vec_close, finish_tests

  integer, save :: nchecks = 0, nfailed = 0

contains

  subroutine check_close(name, got, expected, rtol, atol)
    character(len=*), intent(in)   :: name
    real(dp), intent(in)           :: got, expected
    real(dp), intent(in), optional :: rtol, atol
    real(dp) :: rt, at, tol

    rt = 1e-10_dp
    at = 0.0_dp
    if (present(rtol)) rt = rtol
    if (present(atol)) at = atol
    tol = rt * max(abs(got), abs(expected)) + at
    ! Written so that a NaN fails
    call record(name, abs(got - expected) <= tol)
    if (.not. (abs(got - expected) <= tol)) then
       write(*,'(a,es24.16,a,es24.16,a,es10.3)') '      got ', got, &
            & ' expected ', expected, ' rel diff ', &
            & abs(got - expected) / max(abs(expected), tiny(1.0_dp))
    endif
  end subroutine check_close

  ! Elementwise comparison, with the tolerance set by the largest
  ! component (appropriate for four-momenta).
  subroutine check_vec_close(name, got, expected, rtol, atol)
    character(len=*), intent(in)   :: name
    real(dp), intent(in)           :: got(:), expected(:)
    real(dp), intent(in), optional :: rtol, atol
    real(dp) :: rt, at, tol

    rt = 1e-10_dp
    at = 0.0_dp
    if (present(rtol)) rt = rtol
    if (present(atol)) at = atol
    tol = rt * max(maxval(abs(got)), maxval(abs(expected))) + at
    call record(name, all(abs(got - expected) <= tol))
    if (.not. all(abs(got - expected) <= tol)) then
       write(*,'(a,*(es15.7))') '      got      ', got
       write(*,'(a,*(es15.7))') '      expected ', expected
    endif
  end subroutine check_vec_close

  subroutine check_true(name, cond)
    character(len=*), intent(in) :: name
    logical, intent(in)          :: cond
    call record(name, cond)
  end subroutine check_true

  subroutine record(name, ok)
    character(len=*), intent(in) :: name
    logical, intent(in)          :: ok
    nchecks = nchecks + 1
    if (ok) then
       write(*,'(a,a)') '  PASS  ', name
    else
       nfailed = nfailed + 1
       write(*,'(a,a)') '  FAIL  ', name
    endif
  end subroutine record

  subroutine finish_tests()
    write(*,'(/,i0,a,i0,a)') nchecks - nfailed, ' of ', nchecks, ' checks passed'
    if (nfailed > 0) error stop 1
  end subroutine finish_tests

end module test_utils
