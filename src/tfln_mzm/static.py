"""Static, ideal MZM. SI units; V is the voltage across one electrode gap.

Single arm: one active arm, the other fixed.
Push-pull: both arms see equal and opposite fields of magnitude V/gap.
This is NOT a definition of differential driver voltage across two independent ports.
"""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class StaticParameters:
    wavelength_m: float = 1550e-9
    n_e: float = 2.2
    r33_m_per_v: float = 30e-12
    gap_m: float = 10e-6
    length_m: float = 0.02
    overlap: float = 1.0

    def __post_init__(self):
        values = (self.wavelength_m, self.n_e, self.r33_m_per_v,
                  self.gap_m, self.length_m, self.overlap)
        if not all(np.isfinite(v) and v > 0 for v in values):
            raise ValueError("All parameters must be finite and positive.")
        if self.overlap > 1:
            raise ValueError("This simplified overlap must be in (0, 1].")


def _arm_factor(drive):
    if drive not in ("single_arm", "push_pull"):
        raise ValueError("drive must be single_arm or push_pull")
    return 1 if drive == "single_arm" else 2


def delta_n(voltage_v, p=StaticParameters()):
    """Effective index perturbation for a +z electric field, first order."""
    return -0.5 * p.n_e**3 * p.r33_m_per_v * p.overlap * np.asarray(voltage_v) / p.gap_m


def phase_difference(voltage_v, p=StaticParameters(), drive="single_arm"):
    """Signed driven-arm minus reference-arm phase, radians; no static bias."""
    return _arm_factor(drive) * 2 * np.pi * p.length_m / p.wavelength_m * delta_n(voltage_v, p)


def half_wave_voltage(p=StaticParameters(), drive="single_arm"):
    """Positive gap voltage causing |differential phase| = pi."""
    return p.wavelength_m * p.gap_m / (
        _arm_factor(drive) * p.n_e**3 * p.r33_m_per_v * p.overlap * p.length_m)


def transmission(voltage_v, p=StaticParameters(), drive="single_arm", bias_rad=0.0):
    """Ideal bright-port optical power / input power (balanced, lossless)."""
    return 0.5 * (1 + np.cos(bias_rad + phase_difference(voltage_v, p, drive)))


def dc_gap_voltage(source_v, load_ohm, source_ohm=50.0):
    """Zero-series-resistance line at DC. source_v is Thevenin open-circuit Vg.

    +infinity denotes an open load. This is a DC limiting check, not an RF model.
    """
    load = np.asarray(load_ohm, dtype=float)
    if not np.isfinite(source_ohm) or source_ohm <= 0:
        raise ValueError("source_ohm must be finite and positive")
    if np.any(np.isnan(load)) or np.any(load < 0):
        raise ValueError("load_ohm must be nonnegative or +infinity")
    return np.asarray(source_v) / (1 + source_ohm / np.where(load == 0, np.inf, load)) * (load != 0)
