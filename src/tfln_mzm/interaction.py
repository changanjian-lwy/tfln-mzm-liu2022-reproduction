"""Stable exact integral of Liu 2022 Eq.(1), uniform effective line."""
import numpy as np
from .microwave import propagation, C0
from .termination import reflection_coefficients


def _exprel(z):
    """(exp(z)-1)/z, continued at zero; callers have Re(z)<=0."""
    z=np.asarray(z,complex)
    small=np.abs(z)<1e-5
    safe=np.where(small,1,z)
    regular=np.expm1(safe)/safe
    return np.where(small,1+z/2+z*z/6+z**3/24+z**4/120,regular)


def average_voltage(f_hz, *, length_m, z0, zg, zl, n_m, n_g, alpha_np_m, vg=1.0):
    """Complex Vavg in volts. n_g is optical group index, not phase index.

    Derived by x=L-z in printed Eq.(1), beta_e=omega*n_m/c-j*alpha.
    No fitting, clipping, or inferred electrode parameters.
    """
    if not np.isfinite(length_m) or length_m<=0 or not np.isfinite(n_g) or n_g<=0 or not np.isfinite(vg):
        raise ValueError('Require L>0, n_g>0 and finite source amplitude')
    g=propagation(f_hz,n_m,alpha_np_m)
    bo=2*np.pi*np.asarray(f_hz)*n_g/C0
    r1,r2=reflection_coefficients(z0,zg,zl)
    u=(g-1j*bo)*length_m
    w=(g+1j*bo)*length_m
    return vg*(1+r1)/2 * (_exprel(-u)+r2*np.exp(-u)*_exprel(-w))/(1+r1*r2*np.exp(-2*g*length_m))
