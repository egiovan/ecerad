from ctypes import CDLL, c_int, c_size_t, c_double, POINTER
import os
import numpy as np
from numpy.ctypeslib import ndpointer

_path = os.path.dirname(__file__)
_library = CDLL(os.path.join(_path, "ecerad_py_interface.so"))

#subroutine radiation_temperature(nr, ns, r, ne, te, b, rs, theta, te_rad, te_rad_profile) bind(C)
#    use iso_c_binding
#    integer(c_int), value :: nr, ns
#    real(c_double),intent(in) :: r(nr), ne(nr), te(nr), b(nr)
#    real(c_double),intent(in) :: rs(ns), theta
#    real(c_double), intent(out) :: te_rad(ns), te_rad_profile(nr,ns)

_library.radiation_temperature.restype = None
_library.radiation_temperature.argtypes = [
        c_int,                          # NR
        c_int,                          # NS
        ndpointer(dtype=float),    # r
        ndpointer(dtype=np.float64),    # ne
        ndpointer(dtype=np.float64),    # te
        ndpointer(dtype=np.float64),    # b
        ndpointer(dtype=np.float64),    # rs
        c_double,                       # theta
        ndpointer(dtype=np.float64),    # te_rad
        ndpointer(dtype=np.float64),    # te_rad_profile
    ]

def radiation_temperature(r, ne, te, b, rs, theta):
    """ ECE Radiation temperature.
    r    [nr]   major radius (m)
    ne   [nr]   electron density (m^-3)
    te   [nr]   electron temperature (eV)
    b    [nr]   magnetic field (T)
    rs   [ns]   radii where a resonant radiation temperature is calculated
    theta       angle of the viewing line respect to the perpendicular to the magnetic field

    The array r, ne, te, and b should have the same lenght and specify the plasma profiles.
    the r array should be strictly increasing.

    rs is used to specify the frequency of the radiation temperature. In a thick plasma the 
    radiation temperature at the frequency corresponding to rs will be equal to the actual
    temerature at rs.
    return:
    te_rad [ns]      radiation temperature
    te_rad_profile [ns, nr]  radiation temperature evolution.
    typically te_rad == te_rad_profile[:,-1]

    """
    r = np.ascontiguousarray(r, dtype=np.float64)
    ne = np.ascontiguousarray(ne, dtype=np.float64)
    te = np.ascontiguousarray(te, dtype=np.float64)
    b = np.ascontiguousarray(b, dtype=np.float64)
    rs = np.ascontiguousarray(rs, dtype=np.float64)
    nr = r.size
    assert nr == ne.size and nr == te.size and nr == b.size
    ns = rs.size
    te_rad  = np.zeros(ns, dtype=float)
    te_rad_profile = np.zeros((ns,nr), dtype=float)
    _library.radiation_temperature(nr, ns, r, ne, te, b, rs, theta, te_rad, te_rad_profile)
    return te_rad, te_rad_profile