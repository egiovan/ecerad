program check
    use dawson
    use iso_c_binding, only: wp => c_double
implicit none
interface 
    real(c_double) function Faddeeva_Dawson_re(x) result(r) bind(C,name="Faddeeva_Dawson_re")
    use iso_c_binding, only: c_double
    implicit none
    real(c_double), value :: x

    end function
end interface

real(wp) :: x(10) = [real(wp):: 0.02_wp, 0.01_wp, 0.1_wp, 0.2_wp, 10, 23, 50, 70, 300, 6e7]
real(wp) :: xb(400)
call test(x)
x = -x
call test(x)
call linspace(xb, 0.0_wp, 45.0_wp)
call test(xb)
xb = -xb
call test(xb)
contains
    subroutine linspace(v, vmin, vmax)
        real(wp), intent(out) :: v(:)
        real(wp), intent(in) :: vmin, vmax
        real(wp) :: dv
        integer :: i, n

        n = size(v)
        dv = (vmax - vmin)/(n-1)
        v = [(vmin + dv*(i-1), i=1,n)]
    end subroutine
    subroutine test(x)
        real(wp), intent(in) :: x(:)
        real(wp) :: r(size(x)), ro(size(x))
        integer :: i
        do i=1, size(x)
            ro(i) = Faddeeva_Dawson_re(x(i))
        enddo
        r = dawsn(x)
        print *, "Testing Faddeeva"
        do i=1, size(x)
            if (r(i) /= ro(i)) then
                print *,'Error', x(i)
            else
                write(*,"(A)",advance='no') "."
            endif
        enddo
        print *
    end subroutine
end program check
