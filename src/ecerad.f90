module ecerad
  use ecerad_base
  use dop853_solver, only: solve_dop853
  use flint_solver, only: solve_flint
  use ecerad_utils
  implicit none

  real(wp), parameter :: kb = e ! Temperature in eV
  logical, parameter :: USE_FLINT = .TRUE.

  type Profile
    real(wp), allocatable :: r(:)
    real(wp), allocatable :: ne(:)
    real(wp), allocatable :: te(:)
    real(wp), allocatable :: b(:)
  end type

contains

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
              (1.0_wp/8.0_wp*cos(theta)**4 + sin(theta)**2)/ &
              ((1 + sin(theta)**2)*sqrt(sin(theta)**2 + 1.0_wp/16.0_wp*cos(theta)**4)) 

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

contains
!------------------------------------------------------------------------------------------
  subroutine func(r, y, f)
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
end module ecerad
