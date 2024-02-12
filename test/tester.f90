program tester
    use, intrinsic :: iso_fortran_env, only : error_unit
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type
    use test_suite_faddeeva, only: collect_faddeeva_test
    use test_suite_phi, only: collect_phi_test
    !use test_suite2, only : collect_suite2
    implicit none
    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'
  
    stat = 0
  
    testsuites = [ &
       new_testsuite("Faddeeva. Test the Fortran version is like the original C version", collect_faddeeva_test), &
       new_testsuite("Calculation of the phi", collect_phi_test) &
      ]
  
    do is = 1, size(testsuites)
      write(error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
    end do
  
    if (stat > 0) then
      write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      !error stop
    end if
  
  end program tester