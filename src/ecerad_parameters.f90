module ecerad_parameters
  use iso_c_binding, only: wp => c_double
  implicit none

  real(wp), parameter :: RELATIVE_TOLERANCE = 1.0e-5_wp
  real(wp), parameter :: ABSOLUTE_TOLERANCE = 1.0e-7_wp

  real(wp), parameter :: THETA_TAYLOR_MAX = 0.01_wp 

end module