module ecerad
  use ecerad_base
  use dop853_module, only: dop853_class
  use FLINT, only: DiffEqSys
  implicit none

  real(wp), parameter :: kb = e ! Temperature in eV
  logical, parameter :: USE_FLINT = .TRUE.
  real(wp), parameter :: RELATIVE_TOLERANCE = 1.0e-5_wp
  real(wp), parameter :: ABSOLUTE_TOLERANCE = 1.0e-7_wp

  type Profile
    real(wp), allocatable :: r(:)
    real(wp), allocatable :: ne(:)
    real(wp), allocatable :: te(:)
    real(wp), allocatable :: b(:)
  end type

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

elemental function refractor_index(X, Y, theta) result(r)
  real(wp), intent(in) :: X, Y, theta
  real(wp) :: xy2, y2m, x2m, xym, s2, c2, b, delta, a
  real(wp) :: r

  xy2 = X**2 * Y**2
  y2m = 1 - Y**2
  x2m = 1 - X**2
  xym = 1 - X**2 - Y**2
    
  s2 = cos(theta)**2 ! sin -> cos per theta -> pi/2 - theta
  c2 = sin(theta)**2 ! cos -> sin per theta -> pi/2 - theta
  
  b = -xy2*s2 + 2*xym*y2m
  delta = xy2**2*s2**2 + 4*y2m**2*xy2**2*Y**2*c2
  a = xym*s2 + y2m*x2m*c2
  r = (b - sqrt(delta))/2/a
  !return , (b + np.sqrt(delta))/2/a
end function


subroutine radiation_temperature(r, ne, te, b, rs, theta, te_rad, te_rad_profile)
  real(wp),intent(in) :: r(:), ne(:), te(:), b(:)
  real(wp),intent(in) :: rs(:), theta
  real(wp), intent(out) :: te_rad(:), te_rad_profile(:,:)
  type(Profile) :: pr
  integer :: i

  if (.not. check_equal([size(r), size(ne), size(te), size(b)])) error stop
  if (.not. check_equal([size(rs), size(te_rad)])) error stop
  if (size(r)/= size(te_rad_profile,1) .or. size(rs)/= size(te_rad_profile,2))  error stop

  pr%r = r
  pr%ne=ne
  pr%te=te
  pr%b=b

  do i = 1,size(rs)
    call radiation_evolution(pr, theta, rs(i), te_rad(i), te_rad_profile(:,i))
  enddo
end subroutine

subroutine radiation_evolution(pr, theta, rs, te_rad, te_rad_profile)
  type(Profile), intent(in) :: pr
  real(wp), intent(in) :: theta, rs
  real(wp), intent(out) :: te_rad, te_rad_profile(:)
  real(wp) :: eta_2x, Bs, omega

  real(wp),dimension(size(pr%B)) :: omega_x, omega_p, X, Y, n2
  integer  :: i_start
  real(wp) :: r_start, r_end
  real(wp) :: te_rad_flint, te_rad_profile_flint(size(te_rad_profile))
  real(wp) :: t0, t1, t2 
  Bs = np_interp(rs, pr%r, pr%B)
  omega = 2*e*Bs/me
  eta_2x = 0.5_wp + &
              (1.0_wp/8_wp*cos(theta)**4 + sin(theta)**2)/ &
              ((1 + sin(theta)**2)*sqrt(sin(theta)**2 + 1.0_wp/16_wp*cos(theta)**4)) 

  omega_x = e*pr%B/me
  omega_p = sqrt(pr%ne*e**2/epsilon_0/me)
  
  X = omega_x/omega
  Y = omega_p/omega 
  n2 = refractor_index(X, Y, theta)
  
  i_start = last_greater_than_zero(n2)
  r_start = pr%r(i_start)
  r_end = pr%r(size(pr%r))

  if (USE_FLINT) then 
    call solve_flint(func, pr%r, pr%te, i_start, rs, te_rad_profile, te_rad)
  else
    call solve_dop853(func, pr%r, pr%te, i_start, rs, te_rad_profile, te_rad)
  endif
  !
  !print *, te_rad_profile - te_rad_profile_flint
  !print *,'DOP853,FLINT te_rad', te_rad, te_rad_flint, te_rad_flint/te_rad - 1, t1-t0, t2-t1
contains
!------------------------------------------------------------------------------------------
  subroutine func(r, y, f)
    !! Right-hand side of van der Pol's equation
    implicit none
    real(wp),intent(in)               :: r
    real(wp),dimension(:),intent(in)  :: y
    real(wp),dimension(:),intent(out) :: f

    real(wp) :: trad, te, ne, B, zeta, omega_2x, mu
    real(wp) :: total_emission, sqrt_mu, termine, j2x, dy
    trad = y(1)
    te = np_interp(r, pr%r, pr%te)
    ne = np_interp(r, pr%r, pr%ne)
    B = np_interp(r, pr%r, pr%B)
          
    zeta = mec2/te/2
          
    omega_2x = 2*e*B/me
    mu = (omega/omega_2x)**2
          
    total_emission = (e*omega_2x/4/pi/zeta)**2*ne/epsilon_0/c * cos(theta)**2*(1+sin(theta)**2)
    sqrt_mu = sqrt(mu) 
    termine = total_emission*8*pi**3*c**2/(kb*omega**2) * 2*sqrt_mu/omega_2x
          
    j2x = eta_2x * termine * phi(mu, te, theta)
    dy = j2x*(1 - trad/te)
    f(1) = dy
  end subroutine
end subroutine

subroutine solve_flint(func_flint, rr, te, i_start, rs, te_rad_profile, te_rad)
  use FLINT
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
    !! Right-hand side of van der Pol's equation
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


!===============================================================================================

!---------------------------------------------------------------------------------------
integer function last_greater_than_zero(v) result(r)
  real(wp), intent(in) :: v(:)
  integer :: i, n
  logical :: flag
  n = size(v)
  r = 1
  do i=n, 1, -1
    if (v(i) < 0) then
      r = i + 1
      exit
    endif
  enddo

end function

logical function check_equal(sh) result(r)
  integer :: sh(:)
  r = all(sh(1)==sh(2:))
end function

function np_interp(x, xf, vf) result(r)
  real(wp),intent(in) :: x, xf(:), vf(:)
  real(wp) :: r
  integer :: n, i,is, ie, im
  n = size(xf)
  
  if (x <= xf(1)) then
    r = vf(1)
    RETURN
  endif
  if (x >= xf(n)) then
    r = vf(n)
    return
  endif
! I'm doing it in a stupid way
  if (n < 16) then ! stupid way
    do i=1, size(xf)-1
      if (x < xf(i+1)) then ! is between i, i+1
        r = linear_interp(x, i)
        exit
      endif
    enddo
  else
    is = 1
    ie = n
    do 
      if (is + 1 == ie) exit
      im = (is+ie)/2
      if (x >= xf(im)) then
        is = im
      else
        ie = im
      endif
    enddo
    r = linear_interp(x, is)
  endif
contains
  real(wp) function linear_interp(x, i)
    real(wp), intent(in) :: x
    integer,intent(in) :: i
    linear_interp = (vf(i+1) - vf(i))/(xf(i+1) - xf(i))*(x - xf(i)) + vf(i)
  end function  
end function
end module ecerad
