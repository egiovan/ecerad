# ecerad
Calculate the ECE emission of second harmonic of the extraordinary wave from a tokamak plasma according to the approximations found in[^1].

In particular it assumes the radiation at double the extraordinary wave in an almost perpendicular settings. The plasma refraction index is assumed to be one.

Differently to the article the $\theta$ angular variable is referred to the perpendicular to the magnetic field: $\theta_{Rathgeber} = \pi/2  - \theta$.
Moreover the software assume a propagation with $\theta$ close to zero (perpendicular to the magnetic field) as it generally happens in actual ECE diagnostic. 

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

## Notes on the equations

The equations are those present in Rathgeber[^1], with all $\sin$'s converted to $\cos$'s and viceversa, so that the angle $\theta$ will be equal to zero cloe to the perpedicular to the magnetic field (it is adifferent convention compared to what found in the literature).

Following[^1] we are solving the following equation:

$$ \frac d {ds} T_{rad}(s,\omega)= \frac{8\pi^3c^2}{k_B\omega^2} j_{2X}(s,\omega) \left( 1- \frac{T_{rad}(s,\omega)}{T_e(s)} \right)$$

That is the same equation (7)[^1] with $N^2=1$ having transformed the intensity to temperature usign (9)[^1]. Moreover the $s$ cohordinate is assumend to be coincident with the $r$, the cohordinate used to specify the profiles.

$$j_{2X} = \eta_{2X} \frac{e^2\omega_{2X}^2}{16\pi^2\zeta^2} \frac{n_e}{\epsilon_0 c}  \cos^2 \theta (1+\sin^2 \theta)
\frac {2\omega}{\omega_{2X}^2}\bar\Phi(\mu)$$

with:

$$\omega_{2X} = 2 \frac{eB}{m_e}\text{,    } 
\zeta= \frac{m_ec^2}{2k_BT_e}\text{,    }
\mu=\frac{\omega^2}{\omega^2_{2X}}$$

and:

$$\eta_{2X}=\frac 1 2 + \frac{\frac 1 8 \cos^4 \theta + \sin^2 \theta}{(1+\sin^2 \theta)\sqrt{\sin^2 \theta + \frac 1 {16} \cos^4 \theta}}$$

$$\bar\Phi(\mu) = \frac{\zeta^{3/2}}{2\sqrt\pi \mu^2 \sin^5\theta}\left(I(\alpha_1) -I(\alpha_2)\right)$$

$$I(\alpha) = \exp({-\zeta[1-\mu\alpha]})\left( -\frac {\cos^2\theta}{\sqrt\alpha} + \frac {\eta^2}{\zeta\mu} - \frac 1 {\sqrt{\zeta\mu}}\left[ 3\eta -2\zeta\cos^2\theta F(\sqrt{\zeta\mu\alpha})\right]
\right)$$

Where we have:

$$
\alpha_{1,2} = \frac {(1 \pm \sqrt{1 -\mu\cos^2\theta} \sin\theta)^2}{(1+\mu\sin^2\theta)^2}
$$

and

$$
\eta = 1+\mu\sin^2\theta
$$

Moreover $F(x)=\exp(-x^2) \int_0^x \exp y^2 dy$ is the Dawson integral[^2], which has been translated to Fortran from a C source.
To compare this definition of $\bar\Phi(\mu)$ with that in [^1] one have to remind that $\Phi(\omega)=2\omega/\omega_{2X}^2\bar\Phi(\mu)$, just a change of coordinates.  

Close to $\theta=0$ a taylor approximation of $\bar\Phi(\mu)$ is used. That exact Taylor approximation has been obtained with Maxima [^3], and it is not reported here due to its complexity but can be read on the sources.


[^1]: Rathgeber,Plasma Phys. Control. Fusion 55 (2013) 025004 (15pp)
[^2]: http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package
[^3]:  Maxima, a Computer Algebra System. https://maxima.sourceforge.io/