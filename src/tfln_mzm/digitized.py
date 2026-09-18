"""Bounded interpolation for independent digitized input providers."""
import numpy as np


def interpolate_supported(query, frequency, values, max_gap):
    """Return values and mask; no extrapolation or bridging of long gaps."""
    q=np.asarray(query,float);f=np.asarray(frequency,float);v=np.asarray(values,float)
    if q.ndim!=1 or f.ndim!=1 or v.shape!=f.shape or len(f)<2 or np.any(np.diff(f)<=0) or not np.isfinite(max_gap) or max_gap<=0:
        raise ValueError('Invalid ordered source or gap limit')
    if any(np.any(~np.isfinite(a)) for a in (q,f,v)):
        raise ValueError('Finite samples required')
    right=np.searchsorted(f,q,side='left');safe=np.clip(right,0,len(f)-1)
    exact=(f[safe]==q)
    lo=np.clip(right-1,0,len(f)-2);hi=lo+1
    supported=(q>=f[0])&(q<=f[-1])&(exact|((f[hi]-f[lo])<=max_gap))
    values=np.interp(q,f,v)
    return np.where(supported,values,np.nan),supported
