"""Preflight for scalar paper-baseline manifests, not a physics validator.

Future sampled-data providers need their own coverage/interpolation contracts.
A successful preflight never labels a run as a paper reproduction success.
"""
import math

UNITS = {"length":"m", "source_impedance":"ohm", "source_voltage":"V",
         "loads":"ohm", "reference_frequency":"Hz", "optical_group_index":"1",
         "characteristic_impedance":"ohm", "rf_attenuation":"Np/m", "microwave_index":"1"}
EVIDENCE = {"paper_explicit", "digitized_paper_calculation", "digitized_paper_measurement",
            "external_reference", "assumption", "sensitivity_only", "unknown"}


def preflight(config):
    """Reject malformed/conflicting metadata; return unknown inputs explicitly."""
    if config.get("track") not in {"A", "B", "C"}:
        raise ValueError("Unknown research track")
    if not config.get("experiment_id") or not config.get("scope"):
        raise ValueError("Experiment and claim scope are required")
    params = config["parameters"]
    if set(params) != set(UNITS):
        raise ValueError("Parameter roles missing or unrecognized")
    missing, assumptions = [], []
    for name, unit in UNITS.items():
        row = params[name]
        if row.get("unit") != unit:
            raise ValueError(f"{name}: expected {unit}")
        evidence = row.get("evidence")
        if evidence not in EVIDENCE or not row.get("source", "").strip():
            raise ValueError(f"{name}: evidence and source are required")
        value = row.get("value")
        if evidence == "unknown":
            if value is not None:
                raise ValueError(f"{name}: unknown input cannot carry an invented value")
            missing.append(name)
            continue
        if evidence == "sensitivity_only" and config["track"] == "A":
            raise ValueError(f"{name}: sensitivity-only values cannot be the reproduction baseline")
        if evidence in {"assumption", "external_reference", "sensitivity_only"}:
            assumptions.append(name)
        if name == "loads":
            if not isinstance(value, list) or not value:
                raise ValueError("loads must be a nonempty list")
            values = [v for v in value if v != "open"]
        else:
            values = [value]
        for v in values:
            if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
                raise ValueError(f"{name}: finite scalar required; use 'open' for open load")
            minimum_inclusive = name in {"loads", "rf_attenuation", "reference_frequency"}
            if v < 0 or (v == 0 and not minimum_inclusive):
                raise ValueError(f"{name}: outside baseline parameter domain")
    return {"experiment_id":config["experiment_id"], "scope":config["scope"],
            "status":"BLOCKED_BY_MISSING_DATA" if missing else "INPUTS_READY_NOT_VALIDATED",
            "missing":missing, "nonpaper_inputs":assumptions,
            "physics_validated":False}
