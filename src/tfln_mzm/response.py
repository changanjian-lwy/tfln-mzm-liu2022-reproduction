"""Magnitude interpretation of Eq.(4). Bandwidth extraction pending A02."""
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
