module ecerad_base
    use iso_c_binding, only: wp => c_double
    use dawson
    implicit none

    real(wp), parameter :: mev_to_ev = 1e6_wp
    real(wp), parameter :: c = 299792458.0_wp
    real(wp), parameter :: e = 1.602176634e-19_wp
    real(wp), parameter :: me = 9.1093837015e-31_wp
    real(wp), parameter :: mec2 = 0.51099895000_wp*mev_to_ev ! in eV
    real(wp), parameter :: pi = 3.141592653589793_wp
    real(wp), parameter :: epsilon_0 = 8.8541878128e-12_wp
    
contains
!-------------------------------------------------------------------
elemental function zeta(te) result(r)
    real(wp),intent(in) :: te
    real(wp) :: r
    r = mec2/te/2
end function
!-------------------------------------------------------------------
elemental subroutine alpha(mu, theta, ap, am)
    real(wp), intent(in) :: mu, theta
    real(wp), intent(out) :: ap, am

    real(wp) :: eta

    eta = 1 + mu*sin(theta)**2
    ap = (1 + sqrt(1 - mu*cos(theta)**2)*sin(theta))**2/eta**2
    am = (1 - sqrt(1 - mu*cos(theta)**2)*sin(theta))**2/eta**2
end subroutine

elemental subroutine approx(z, mu, ap5, ap7, ap9)
    real(wp),intent(in) :: z, mu
    real(wp), intent(out) :: ap5, ap7, ap9

    ap5 = -(16*z**2*(mu-1)**3*mu**2*exp(z*mu-z))/(15*sqrt(1-mu))
    ap7 = (8*z**2*(mu-1)**2*mu**2*(12*z**2*mu**4-24*z**2*mu**3+102*z*mu**3+ &
           12*z**2*mu**2-120*z*mu**2+168*mu**2+18*z*mu-28*mu-35)*exp(z*mu-z))/(315*sqrt(1-mu))
    ap9 = -(2*z**2*(mu-1)*mu**2*(16*z**4*mu**8-64*z**4*mu**7+368*z**3*mu**7+96*z**4*mu**6-1184*z**3*mu**6+ &
            2652*z**2*mu**6-64*z**4*mu**5+1344*z**3*mu**5-6336*z**2*mu**5+6840*z*mu**5+16*z**4*mu**4 &
            -608*z**3*mu**4+4608*z**2*mu**4-10272*z*mu**4+5040*mu**4+80*z**3*mu**3-816*z**2*mu**3 &
            +2292*z*mu**3-2688*mu**3-108*z**2*mu**2+1392*z*mu**2-2128*mu**2-252*z*mu+560*mu+161)*exp(z*mu-z))/(945*sqrt(1-mu))
end subroutine

elemental subroutine sin5_approx(s5, s7, s9)
    real(wp), intent(out) :: s5, s7, s9
    s5 = 1.0_wp
    s7 = -5.0_wp/6.0_wp
    s9 = 23.0_wp/72.0_wp
    ! x^5-(5*x^7)/6+(23*x^9)/72-(227*x^11)/3024
end subroutine

elemental function tayrat(a, b, c, aa, bb, cc, x) result(r)
    real(wp), intent(in) :: a, b, c, aa, bb, cc, x
    real(wp) :: r
    r = a/aa -((a*bb-aa*b)*x**2)/aa**2-((a*aa*cc-aa**2*c-a*bb**2+aa*b*bb)*x**4)/aa**3

end function

elemental function main_approx(z, mu, theta) result(r)
    real(wp), intent(in) :: z, mu, theta
    real(wp) :: r
    real(wp) :: s5, s7, s9
    real(wp) :: ap5, ap7, ap9
    real(wp) :: f1
    
    if (mu >= 1) then
        r = 0.0
        return
    endif
    f1 = z**1.5/(sqrt(pi) *  mu**2)
    call sin5_approx(s5, s7, s9)
    call approx(z, mu, ap5, ap7, ap9)
    r = tayrat(ap5, ap7, ap9, s5, s7, s9, theta)*f1
end function

elemental function interno(theta, alpha, z, mu, eta) result(r)
    real(wp),intent(in) :: theta, alpha, z, mu, eta
    real(wp) :: r
    real(wp) :: f2, fp_1, fp_2, fp_3

    f2 = exp(-z*(1-mu*alpha))
    fp_1 = - cos(theta)**2/sqrt(alpha)
    fp_2 = eta**2/z/mu
    fp_3 = -1/sqrt(z*mu)*(3*eta - 2*z*mu*cos(theta)**2)*dawsn(sqrt(z*mu*alpha))

    r = f2*(fp_1 + fp_2 + fp_3)
end function

elemental function phi(mu, te, theta) result(r)
    real(wp), intent(in) :: mu, te
    real(wp), intent(in), optional :: theta
    real(wp) :: r
    real(wp) :: f1, eta, alpha1, alpha2, int1, int2, z, theta_loc
    z = zeta(te)
    theta_loc = 0.0
    if (present(theta)) theta_loc = theta

    f1 = z**1.5/(sqrt(pi) *  mu**2 * sin(theta_loc)**5)/2
    eta = 1 + mu*sin(theta_loc)**2
    !alpha1 = (1 + np.sqrt(1 - mu*np.cos(theta)**2)*np.sin(theta))**2/eta**2
    !alpha2 = (1 - np.sqrt(1 - mu*np.cos(theta)**2)*np.sin(theta))**2/eta**2
    if (abs(theta_loc) < 0.01) then
        r = main_approx(z, mu, theta_loc)
        return
    end if
    call alpha(mu, theta_loc, alpha1, alpha2)
    
    int1 = interno(theta_loc, alpha1, z, mu, eta)
    int2 = interno(theta_loc, alpha2, z, mu, eta)
    r = f1*(int1 - int2)

end function


end module