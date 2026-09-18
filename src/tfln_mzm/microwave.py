"""Uniform effective transmission line; SI units, real positive impedance."""
import numpy as np
C0 = 299792458.0


def propagation(f_hz, n_m, alpha_np_m):
    """gamma=alpha+j*omega*n_m/c; forward voltage exp(-gamma*z)."""
    f, n, a = np.broadcast_arrays(np.asarray(f_hz, float), np.asarray(n_m, float), np.asarray(alpha_np_m, float))
    if not all(np.all(np.isfinite(x)) for x in (f,n,a)) or np.any(f<0) or np.any(n<=0) or np.any(a<0):
        raise ValueError('Require finite f>=0, n_m>0, alpha>=0')
    return a + 2j*np.pi*f*n/C0


def db_per_mm_to_np_per_m(loss):
    loss=np.asarray(loss,float)
    if np.any(~np.isfinite(loss)) or np.any(loss<0):
        raise ValueError('Loss must be finite and nonnegative')
    return loss*1000*np.log(10)/20
