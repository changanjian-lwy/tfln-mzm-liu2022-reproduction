"""Magnitude interpretation of Eq.(4) and first -3 dB threshold extraction."""
import numpy as np


def magnitude_db(value):
    value=np.asarray(value,complex)
    if np.any(~np.isfinite(value)):
        raise ValueError('Require finite phasors')
    with np.errstate(divide='ignore'):
        return 20*np.log10(np.abs(value))


def eo_response_db(voltage, reference_voltage):
    reference_voltage=np.asarray(reference_voltage,complex)
    if np.any(~np.isfinite(reference_voltage)) or np.any(np.abs(reference_voltage)==0):
        raise ValueError('Reference voltage is zero or invalid; normalization undefined')
    return magnitude_db(voltage/reference_voltage)


def bandwidth_3db(frequency_hz, response_db):
    """First -3 dB contact/downcrossing, relative to supplied f0 normalization.

    Requires a strictly increasing grid starting at f0 and a zero-dB first
    sample. Censoring is conditional on sampling; callers must refine grids.
    """
    f=np.asarray(frequency_hz,float)
    y=np.asarray(response_db,float)
    if f.ndim!=1 or y.shape!=f.shape or len(f)<2 or np.any(~np.isfinite(f)) or np.any(f<0) or np.any(np.diff(f)<=0):
        raise ValueError('Require equal 1D arrays and finite increasing frequencies')
    if np.any(np.isnan(y)) or np.any(np.isposinf(y)) or not np.isclose(y[0],0,atol=1e-8,rtol=0):
        raise ValueError('Response must start at 0 dB at f0; NaN/+inf forbidden')
    events=[]
    for i in range(1,len(f)):
        if y[i-1]>-3 and y[i]<=-3:
            if y[i]==-3:
                value=float(f[i]);kind='threshold_contact'
            elif np.isneginf(y[i]):
                value=None;kind='bracket_requires_refinement'
            else:
                value=float(f[i-1]+(-3-y[i-1])/(y[i]-y[i-1])*(f[i]-f[i-1]));kind='downcrossing'
            events.append({'frequency_hz':value,'bracket_hz':[float(f[i-1]),float(f[i])],'kind':kind})
    first=events[0] if events else None
    return {'bandwidth_hz':first['frequency_hz'] if first else None,
            'censored':not bool(events),'lower_bound_hz':float(f[-1]) if not events else None,
            'status':first['kind'] if first else 'above_scan_limit','events':events}
