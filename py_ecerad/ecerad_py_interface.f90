module ecerad_py_interface
    use iso_c_binding
    use ecerad, only: rad_te => radiation_temperature, j2x_for_extern_call
    use ecerad_base, only: phi_int => phi_integral, phi_tayl => phi_taylor
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

  subroutine j2x_python(nr, rp, ne, te, b, theta, rs, n, r, j2x, taylor) bind(C)
    integer,value :: nr
    real(c_double),intent(in) :: rp(nr)
    real(c_double),intent(in) :: ne(nr)
    real(c_double),intent(in) :: te(nr)
    real(c_double),intent(in) :: b(nr)

    real(c_double),value :: theta
    real(c_double),value :: rs
    integer,value :: n
    real(c_double),intent(in) :: r(n)
    real(c_double),intent(out) :: j2x(n)
    integer, value :: taylor

    call j2x_for_extern_call(rp, ne, te, b, theta, rs, r, j2x, taylor==1)
  end subroutine

  subroutine phi_taylor(n, mu, te, theta, r) bind(C)
    integer, value :: n
    real(c_double), intent(in) :: mu(n)
    real(c_double), value :: te
    real(c_double), value :: theta
    real(c_double), intent(out) :: r(n)
    r = phi_tayl(mu, te, theta)
  end subroutine

end module
