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
  real(wp),intent(out) :: te_rad(:), te_rad_profile(:,:)
end subroutine
```

<table>
<caption style="text-align:left"> Description of the plasma profiles </caption>
<tr><td> r: </td><td> Major radius (m) </td></tr>
<tr><td> ne: </td><td> electron density (m^-3)</td></tr>
<tr><td> te: </td><td> electron temperatue (eV)</td></tr>
<tr><td> b: </td><td>  magnetic field (T)</td></tr>
</table>


<table>
<caption style="text-align:left"> Description of the receiving apparatus </caption>
<tr><td> theta: </td><td> angle to the perpendicular of the magnetic field </td></tr>

<tr><td> rs: </td><td>   major radius where 2*omega_x will be resonant </td></tr>
</table>

<table>
<caption style="text-align:left"> Output </caption>
<tr><td> te_rad: </td><td> radiation temperature (eV) </td></tr>
<tr><td> te_rad_profile: </td><td> evolution of the radiation temperature along the major radius.  (eV) </td></tr>
</table>

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

Compilation:

A Makefile is used to both launch `fpm` and to compile a shared library to be used by Python codes.

The shared library is in `./py_ecerad` and the module is `ecerad.py`


