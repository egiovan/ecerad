module test_suite_phi
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use ecerad_base, only: phi_integral, phi_taylor
  use ecerad_parameters, only: wp
  use ecerad_utils, only: linspace
  implicit none
  real(wp), parameter :: MU_START = 0.001
  real(wp), parameter :: MAX_ERROR_ON_MU = 1e-4

  public collect_phi_test

contains
  subroutine collect_phi_test(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    character(len=90) :: diff_tailor_message
    write(diff_tailor_message,"(A,F8.5)") "Test Taylor almost equal original within:", MAX_ERROR_ON_MU
    testsuite = [ &
      new_unittest("Test on some precalculated values", test_on_some_precalculated_values_one), &
      new_unittest("Test phi integral close to 1", test_phi_integral), & 
      new_unittest(trim(diff_tailor_message), test_taylor_almost_equal_original) &
    ]
  end subroutine

  subroutine test_taylor_almost_equal_original(error)
    ! Test that where both formula "Taylor expansion" and original are
    ! valid they are both close within a prearranged value
    type(error_type), allocatable,intent(out) :: error
    real(wp), parameter :: te = 4e3 ! eV
    real(wp), parameter :: theta = 0.01 ! it should be ok for both taylor and integral  
    integer, parameter ::  N = 800
    real(wp) :: mu(N), phit(N), phio(N), r, dr(N)
    integer :: i
    call linspace(mu, MU_START, 1.0_wp)
    phio = phi_integral(mu, te, theta)
    phit = phi_taylor(mu, te, theta)
    r = maxval(abs(phio - phit))/maxval(phio)
    call check(error, r < MAX_ERROR_ON_MU)
    if (allocated(error)) return
    
    dr = 0
    where (phio>1e-3_wp) dr = abs(phit/phio - 1)
    ! Except for mu close to 1 the difference is far smaller.
    ! mu equal 1 is where things can break.
    call check(error, maxval(dr) < 1e-2_wp)
    if (allocated(error)) return

    ! Just some code to check if something is going wrong
    ! it shoudn't print anything with this settings
    do i=1,N
      if (dr(i)> 1e-2_wp) then
        print *, i, mu(i), dr(i), phio(i), phit(i)
      endif
    enddo

  end subroutine
!----------------------------------------
  subroutine test_phi_integral(error)
    ! Test that the integral of the phi is one. Well, close to one. 
    type(error_type), allocatable,intent(out) :: error
    real(wp), parameter :: te = 4e3 ! eV
    real(wp), parameter :: theta = 0.01 ! it should be ok for both taylor and integral
    integer, parameter ::  N = 400
    real(wp) :: mu(N), phit(N), r

    call linspace(mu, MU_START, 1.0_wp)
    phit = phi_integral(mu, te, theta)
    r = trapz(mu, phit)

    call check(error, abs(r-1) < MAX_ERROR_ON_MU)

  end subroutine
!-------------------------------
  function trapz(x, y) result(r)
    real(wp), intent(in) :: x(:), y(:)
    real(wp) :: r
    integer :: n
    n = size(x)
    r = sum((x(2:n) - x(1:n-1))*(y(2:n) + y(1:n-1))/2)
  end function

  subroutine test_on_some_precalculated_values_one(error)
    ! The comparison has been done on some precalculated data from a Python code
    ! That does mean that the Fortran version is giving the same results of
    ! the Python code. The Taylor version of the Python code had the same error 
    ! of a factor 2 that was previosly corrected.
    ! This is why "re" is divide by two.
    type(error_type), allocatable,intent(out) :: error
    real(wp) :: mu(40)
    real(wp) :: r(size(mu)), dr(size(mu)), re(size(mu)), re2(size(mu))
    re = [ &
    5.79064369d-13, 1.75016620d-12, 5.27804891d-12, 1.58799671d-11, &
    4.76587357d-11, 1.42652809d-10, 4.25778575d-10, 1.26696690d-09, &
    3.75774556d-09, 1.11061000d-08, 3.26999004d-08, 9.58838284d-08, &
    2.79900800d-07, 8.13106837d-07, 2.34950194d-06, 6.74930381d-06, &
    1.92633307d-05, 5.45862613d-05, 1.53444403d-04, 4.27467448d-04, &
    1.17875204d-03, 3.21279320d-03, 8.64001479d-03, 2.28748349d-02, &
    5.94559287d-02, 1.51163723d-01, 3.74126369d-01, 8.95440849d-01, &
    2.05316884d+00, 4.44738107d+00, 8.90091170d+00, 1.58390263d+01, &
    2.32329361d+01, 2.33311693d+01, 7.67748271d+00, 0.00000000d+00, &
    0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00]/2.0_wp
    re2 = [ &
    3.62840187d-13, 1.14162172d-12, 3.58672516d-12, 1.12444357d-11, &
    3.51511568d-11, 1.09498010d-10, 3.39658378d-10, 1.04846532d-09, &
    3.21845290d-09, 9.81804132d-09, 2.97430522d-08, 8.94169471d-08, &
    2.66567614d-07, 7.87431645d-07, 2.30293374d-06, 6.66243905d-06, &
    1.90483197d-05, 5.37650561d-05, 1.49644788d-04, 4.10181113d-04, &
    1.10559518d-03, 2.92532104d-03, 7.58270900d-03, 1.92082974d-02, &
    4.74100724d-02, 1.13592263d-01, 2.62932793d-01, 5.84268833d-01, &
    1.23564131d+00, 2.45640398d+00, 4.50475364d+00, 7.38986895d+00, &
    1.02507543d+01, 1.06391454d+01, 5.71692616d+00, 4.02252456d-02, &
    0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00]
    call linspace(mu, 0.2_wp, 1.10_wp)

    r = phi_taylor(mu, 5.0e3_wp, 0.0_wp)
    dr = 0.0
    where (re>0) dr = abs(r/re - 1)
    call check(error, maxval(dr) < 1e-8)
    if (allocated(error)) return

    r = phi_integral(mu, 5.0e3_wp, 0.1_wp)
    dr = 0.0
    where (re2>0) dr = abs(r/re2 - 1)

    call check(error, maxval(dr) < 1e-8)
    if (allocated(error)) return

  end subroutine
end module