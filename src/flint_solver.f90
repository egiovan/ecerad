module flint_solver
  use FLINT, only: DiffEqSys, ERK_class, ERK_DOP54
  use ecerad_solver_commons, only: solver_func
  use ecerad_parameters, only: wp, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE
  implicit none
  
  PRIVATE
  public solve_flint

  type , extends(DiffEqSys) :: flint_sys
    procedure(solver_func), pointer, nopass :: step
  contains
    procedure :: F => flint_func     ! Differential equations
  end type flint_sys

contains

  function flint_func(me, X, Y, Params)
    implicit none
    class(flint_sys), intent(inout) :: me !< Differential Equation object
    real(WP), intent(in) :: X
    real(WP), intent(in), dimension(me%n) :: Y
    real(WP), intent(in), dimension(:), optional :: Params
    real(WP), dimension(size(Y)) :: flint_func
    call me%step(X,Y,flint_func)
  end function

  
subroutine solve_flint(func_flint, rr, te, i_start, rs, te_rad_profile, te_rad)
  procedure(solver_func) :: func_flint

  real(wp),intent(in) :: rr(:), te(:), rs
  integer,intent(in) :: i_start 
  real(wp) :: te_rad, te_rad_profile(:), y0(1), yf(1), h

  real(wp) :: r_start, r_end
  type(ERK_class) :: erk
  type(flint_sys),target :: prop

  prop%n = 1
  prop%step => func_flint

  call erk%Init(prop, 10000, Method=ERK_DOP54, &
     ATol=[ABSOLUTE_TOLERANCE], RTol=[RELATIVE_TOLERANCE], &
   InterpOn=.TRUE., eventsOn=.False., maxstepsize=0.02_wp)
  r_start = rr(i_start)
  r_end = rr(size(rr))
  y0 = 0.0
  h = 0.0
  call erk%Integrate(r_start, y0, r_end, yf, StepSz=h, IntStepsOn=.FALSE.)
  te_rad = yf(1)
  te_rad_profile = 0.0
  call erk%Interpolate(rr(i_start:), te_rad_profile(i_start:))
end subroutine
end module
