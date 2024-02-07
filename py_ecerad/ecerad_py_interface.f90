module ecerad_py_interface
    use iso_c_binding
    use ecerad, only: rad_te => radiation_temperature
    implicit none
contains
subroutine radiation_temperature(nr, ns, r, ne, te, b, rs, theta, te_rad, te_rad_profile) bind(C)
    use iso_c_binding
    integer(c_int), value :: nr, ns
    real(c_double),intent(in) :: r(nr), ne(nr), te(nr), b(nr)
    real(c_double),intent(in) :: rs(ns)
    real(c_double),value :: theta
    real(c_double), intent(out) :: te_rad(ns), te_rad_profile(nr,ns)
  
    call rad_te(r, ne, te, b, rs, theta, te_rad, te_rad_profile)
    
  end subroutine
end module
