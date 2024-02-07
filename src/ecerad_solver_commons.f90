module ecerad_solver_commons
  use ecerad_base, only: wp
  implicit none

  real(wp), parameter :: RELATIVE_TOLERANCE = 1.0e-5_wp
  real(wp), parameter :: ABSOLUTE_TOLERANCE = 1.0e-7_wp

  abstract interface
    subroutine solver_func(r, y, f)
      import wp
      !! Right-hand side of van der Pol's equation
      implicit none
      real(wp),intent(in)               :: r
      real(wp),dimension(:),intent(in)  :: y
      real(wp),dimension(:),intent(out) :: f
    end subroutine
  end interface

end module