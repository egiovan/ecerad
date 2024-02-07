program check_integration
    use ecerad
    use ecerad_utils, only: linspace
    implicit none
    integer, parameter :: NR = 100
    integer, parameter :: NS = 4
    real(wp) :: r(NR), te(NR), ne(NR), b(NR)
    real(wp) :: rs(NS), te_rad(NS), te_rad_profile(NR,NS)
    real(wp) :: te_rad_solution(NS)
    real :: t1, t2

    call basic_profiles(8.0e3_wp, 0.9e20_wp, 2.9_wp, r, te, ne, b)
    rs = [3.5_wp, 3.75_wp, 3.80_wp, 3.85_wp]
    te_rad_solution = [4.64560774_wp, 1.60564565_wp, 0.76468464_wp, 0.98271209_wp]*1e3
    call cpu_time(t1)
    call radiation_temperature(r, ne, te, b, rs, 0.0_wp, te_rad, te_rad_profile)
    call cpu_time(t2)
    print *, "Time 1", NR, NS, t2-t1, (t2-t1)/NS
    print *,te(NR/2) ==  8009.9839679599672_wp
    !print *, rs, te_rad, np_interp(rs(1), r, te), te(500)
    print *, te_rad / te_rad_solution - 1
    print *, "massimo errore", maxval(abs(te_rad / te_rad_solution - 1))
    block
      real(wp) :: te_rad(NR), te_rad_profile(NR,NR)
      integer :: i
      te_rad = 0
      call cpu_time(t1)
      call radiation_temperature(r, ne, te, b, r(2:), 0.0_wp, te_rad(2:), te_rad_profile(:,2:))
      call cpu_time(t2)
      print *, "Time 2",NR, NR-1, t2-t1, (t2-t1)/(NR-1)
      do i=1, NR
        write(700, *) r(i), te(i), te_rad(i)
      enddo
    end block
    block
      real(wp) :: te_rad(3), te_rad_profile(NR,3), rs(3)
      integer :: i
      te_rad = 0
      print *,'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
      rs = [2.25_wp, 3.5_wp, 3.75_wp]
      call radiation_temperature(r, ne, te, b, rs, 0.0_wp, te_rad, te_rad_profile)
      print *,'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
      print *,'TERAD', te_rad
      do i=1, NR
        write(800, *) r(i), te(i), te_rad_profile(i,:)
      enddo
    end block
contains
  subroutine basic_profiles(te0, ne0, b0, r, te, ne, b)
    real(wp),intent(in) :: te0, ne0, b0
    real(wp), intent(out) :: r(:), te(:), ne(:), b(:)
    real(wp),parameter :: R0 = 3.0_wp
    real(wp) :: te_edge, ne_edge

    call linspace(r, 2.0_wp, 4.0_wp)

    te_edge = 10.0_wp
    ne_edge= 1e17_wp
    te = te0 * (1-(r-R0)**2)**2 + te_edge
    ne = ne0 * sqrt(1-(r-R0)**2) + ne_edge
    b = b0*R0/r
    
    call linear_edge(r, te, 3.8_wp, 0.02_wp, te_edge)
    call linear_edge(r, ne, 3.8_wp, 0.02_wp, ne_edge)
  end subroutine
  subroutine linear_edge(r, ve, r_edge, dr, v_edge)
    real(wp),intent(in) :: r(:), r_edge, dr, v_edge
    real(wp),intent(inout) :: ve(:) 
    real(wp) :: r_min, r_max, ve_up

    r_min = r_edge - dr
    r_max = r_edge + dr
    ve_up = np_interp(r_min, r, ve)

    where (r> r_edge) ve = v_edge
    where (r >= r_min .and. r <= r_max)
      ve = (ve_up - v_edge)/(r_min - r_max)*(r - r_max) + v_edge
    endwhere
  end subroutine

!----------------------------------------------------------------

end program
