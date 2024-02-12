module test_suite_faddeeva
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use dawson, only: dawsn
  use ecerad_utils, only: linspace
  use iso_c_binding, only: wp => c_double
  implicit none

  interface 
    real(c_double) function Faddeeva_Dawson_re(x) result(r) bind(C,name="Faddeeva_Dawson_re")
    use iso_c_binding, only: c_double
    implicit none
    real(c_double), value :: x

    end function
  end interface

contains
  subroutine collect_faddeeva_test(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
      new_unittest("Test on some values", test_on_some_values), &
      new_unittest("Test on some negative values", test_on_some_negative_values), &
      new_unittest("Test on linspaced values", test_on_linspaced_values) &
    ]
  end subroutine

  subroutine test_on_some_values(error)
    type(error_type), allocatable,intent(out) :: error 
    real(wp) :: x(10) = [real(wp):: 0.02_wp, 0.01_wp, 0.1_wp, 0.2_wp, 10, 23, 50, 70, 300, 6e7]
    call check(error,test_equal(x))
  end subroutine
  
  subroutine test_on_some_negative_values(error)
    type(error_type), allocatable,intent(out) :: error 
    real(wp) :: x(10) = [real(wp):: 0.02_wp, 0.01_wp, 0.1_wp, 0.2_wp, 10, 23, 50, 70, 300, 6e7]
    call check(error,test_equal(-x))
  end subroutine

  subroutine test_on_linspaced_values(error)
    type(error_type), allocatable,intent(out) :: error 
    real(wp) :: x(400)
    
    call linspace(x, 0.0_wp, 45.0_wp)
    
    call check(error,test_equal(x))
    if (allocated(error)) return

    call check(error,test_equal(-x))

  end subroutine
 
  logical function test_equal(x)
    real(wp), intent(in) :: x(:)
    real(wp) :: r(size(x)), ro(size(x))
    integer :: i
    do i=1, size(x)
        ro(i) = Faddeeva_Dawson_re(x(i))
    enddo
    r = dawsn(x)
    if (any(ro /= r)) then
      test_equal = .FALSE.
    else
      test_equal = .TRUE.
    endif
  
  end function
end module
