module dop853_solver
  use dop853_module, only: dop853_class
  use ecerad_utils
  use ecerad_solver_commons
  implicit none

  private
  public solve_dop853

contains

subroutine solve_dop853(func, rr, te, i_start, rs, te_rad_profile, te_rad)
  procedure(solver_func) :: func

  real(wp),intent(in) :: rr(:), te(:), rs
  integer,intent(in) :: i_start 
  real(wp) :: te_rad, te_rad_profile(:)

  type(dop853_class) :: prop
  logical :: status_ok
  integer :: idid, ix
  real(wp) :: y_sol(1), rtol(1), atol(1), r, f_sol(1)
  real(wp) :: r_start, r_end
  call prop%initialize(n=1, fcn=fvpol, status_ok=status_ok, solout=solout_func, icomp=[1], hmax=0.02_wp)
  if (.not. status_ok) error stop "Error initializing solver"
  
  y_sol = 0.0
  rtol = RELATIVE_TOLERANCE
  atol = ABSOLUTE_TOLERANCE
  r_start = rr(i_start)
  r_end = rr(size(rr))
  r = r_start
  ix = i_start
  te_rad_profile = 0.0
  call prop%integrate(r,y_sol,r_end,rtol,atol,iout=3,idid=idid)

  if (idid < 0) then
    print *,'--------------------------------------------------------------------'
    print *,'IDID:', idid, i_start,  y_sol
    print *,'IDID:', r_start, r, rs, rr(i_start + 1)
    print *,'---- ', np_interp(r_start, rr,te),np_interp(r, rr,te), np_interp(rs, rr,te)
    y_sol = 0.0
    call fvpol(prop, r, y_sol, f_sol)
    print *,'     ', f_sol
    print *,'---------------------------------------------------------------------'
  endif

  te_rad = y_sol(1)
contains
  subroutine fvpol(this, r, y, f)
    implicit none

    class(dop853_class),intent(inout) :: this
    real(wp),intent(in)               :: r
    real(wp),dimension(:),intent(in)  :: y
    real(wp),dimension(:),intent(out) :: f
    call func(r,y,f)
  end subroutine

  subroutine solout_func(this,nr,xold,x,y,irtrn,xout)
    !! `solout` furnishes the solution `y` at the `nr`-th
    !! grid-point `x` (thereby the initial value is
    !! the first grid-point).

    class(dop853_class),intent(inout) :: this
    integer,intent(in)                :: nr    !! grid point (0,1,...)
    real(wp),intent(in)               :: xold  !! the preceding grid point
    real(wp),intent(in)               :: x     !! current grid point
    real(wp),dimension(:),intent(in)  :: y     !! state vector \( y(x) \) [size n]
    integer,intent(inout)             :: irtrn !! serves to interrupt the integration. if
                                              !! `irtrn` is set `<0`, [[dop853]] will return to
                                              !! the calling program. if the numerical solution
                                              !! is altered in `solout`, set `irtrn = 2`.
    real(wp),intent(out)              :: xout  !! `xout` can be used for efficient intermediate output
                                              !! if one puts `iout=3`. when `nr=1` define the first
                                              !! output point `xout` in `solout`. the subroutine
                                              !! `solout` will be called only when `xout` is in the
                                              !! interval `[xold,x]`; during this call
                                              !! a new value for `xout` can be defined, etc.
    !print *,'SOLOUT ---------------------------------------------'
    do 
      if (ix>size(rr)) exit
      if (rr(ix)>= xold .and. rr(ix) <= x) then
        te_rad_profile(ix) = this%contd8(1, rr(ix))
        !print *,'SOLOUT', ix, pr%r(ix), this%contd8(1, pr%r(ix))
        ix = ix + 1
      else
        xout = rr(ix)
        exit
      endif
    enddo

  end subroutine solout_func
end subroutine
end module