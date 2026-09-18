"""Paper reflection conventions and S11 at the electrode input plane."""
import numpy as np


def reflection_coefficients(z0, zg, zl):
    z0,zg,zl=np.broadcast_arrays(np.asarray(z0,float),np.asarray(zg,float),np.asarray(zl,float))
    if np.any(~np.isfinite(z0)) or np.any(z0<=0) or np.any(~np.isfinite(zg)) or np.any(zg<=0):
        raise ValueError('Require real finite positive Z0 and Zg')
    if np.any(np.isnan(zl)) or np.any(zl<0):
        raise ValueError('Require real passive load; +inf means open')
    rho1=(z0-zg)/(z0+zg)
    with np.errstate(invalid='ignore'):
        rho2=np.where(np.isposinf(zl),1,(zl-z0)/(zl+z0))
    return rho1,rho2


def s11(gamma, length_m, z0, zl, reference_ohm=50.0):
    if not np.isfinite(length_m) or length_m<=0:
        raise ValueError('Length must be finite and positive')
    g=np.asarray(gamma,complex)
    if np.any(~np.isfinite(g)) or np.any(g.real<0):
        raise ValueError('Require finite passive propagation')
    r,rl=reflection_coefficients(z0,reference_ohm,zl)
    q=rl*np.exp(-2*g*length_m)
    return (r+q)/(1+r*q)
