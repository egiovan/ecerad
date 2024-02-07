module ecerad_utils
  use iso_c_binding, only: wp => c_double
  implicit none



contains

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

subroutine linspace(v, vmin, vmax)
  real(wp), intent(out) :: v(:)
  real(wp), intent(in) :: vmin, vmax
  real(wp) :: dv
  integer :: i, n

  n = size(v)
  dv = (vmax - vmin)/(n-1)
  v = [(vmin + dv*(i-1), i=1,n)]
end subroutine
end module