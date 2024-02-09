module ecerad_solver_commons
  use ecerad_parameters, only: wp
  implicit none
  PRIVATE
  
  public solver_func


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