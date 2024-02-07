# ecerad
Calculate ECE emission from a tokamak plasma according to the approximations found in:

"Rathgeber,Plasma Phys. Control. Fusion 55 (2013) 025004 (15pp)"

In particular it assumes the radiation at double the extraordinary wave in an almost perpendicular settings. The plasma refraction index is assumed to be one.

Differently to the article the *theta* angular variable is referred to the perpendicular to the magnetic field: *theta_rathgeber* = *pi*/2  - *theta*.

It can use two kind of solver:

- https://github.com/egiovan/FLINT (original https://github.com/princemahajan/FLINT)
- https://github.com/jacobwilliams/dop853

Ther is just a routine with the following interface:
```fortran
subroutine radiation_temperature(r, ne, te, b, rs, theta, te_rad, te_rad_profile)
  real(wp),intent(in) :: r(:), ne(:), te(:), b(:)
  real(wp),intent(in) :: rs(:), theta
  real(wp), intent(out) :: te_rad(:), te_rad_profile(:,:)
end subroutine
```
Description of the plasma profiles:

r:     Major radius (m)
ne:    electron density (m^-3)
te:    electron temperatue (eV)
b:     magnetic field (T)

Description of the receiving apparatus:

theta:  angle to the perpendicular of the magnetic field 
rs:     major radius where 2*omega_x will be resonant
te_rad: radiation temperature (eV)
te_rad_profile: evolution of the radiation temperature along the major radius. 

For thick plasma generally te is close to te_rad, this is typically not true in case of a cutoff, or at the edge of the plsama where the density is too low and the so called shine-thorugh happens. 

Example:
```fortran
use ecerad, only: radiation_temperature
use iso_c_binding, only: wp => c_double

integer, parameter :: NR = ..., NS = ...
real(wp) :: r(NR), ne(NR), te(NR), b(NR)
real(wp) :: rs(NS), te_rad(NS), te_rad_profile(NR,NS)
real(wp) :: theta 
...

call radiation_temperature(r, ne, te, b, rs, theta, te_rad, te_rad_profile)

...
```





