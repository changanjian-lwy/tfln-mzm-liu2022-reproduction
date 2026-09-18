"""Reproduce Day 1 figures and data from the repository root or any directory."""
from pathlib import Path
from dataclasses import asdict
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tfln_mzm.static import StaticParameters, delta_n, half_wave_voltage, transmission, dc_gap_voltage

ROOT = Path(__file__).resolve().parents[1]


def main():
    p = StaticParameters()
    vpi = half_wave_voltage(p)
    v = np.linspace(-2*vpi, 2*vpi, 1001)
    frame = pd.DataFrame({"voltage_gap_V":v, "delta_n_effective":delta_n(v,p),
        "T_single_arm":transmission(v,p), "T_push_pull":transmission(v,p,"push_pull"),
        "T_quadrature":transmission(v,p,bias_rad=np.pi/2), "provenance":"simulation"})
    frame.to_csv(ROOT/"data/simulated/day01_static.csv", index=False)
    loads = [20,40,50,80,np.inf]
    pd.DataFrame({"load_ohm":[str(r) for r in loads], "Vg_V":1.0,
        "gap_voltage_V":dc_gap_voltage(1,loads), "provenance":"derived_dc_limit"}).to_csv(ROOT/"data/simulated/day01_dc_limits.csv", index=False)
    summary = {"provenance":"simulation_from_lecture_parameters", "parameters_SI":asdict(p),
        "single_arm_Vpi_V":vpi, "push_pull_gap_Vpi_V":half_wave_voltage(p,"push_pull"),
        "delta_n_at_single_arm_Vpi":float(delta_n(vpi,p)), "single_arm_VpiL_V_cm":vpi*p.length_m*100,
        "quadrature_slope_per_V":float(np.pi/(2*vpi)),
        "boundary":"Ideal balanced lossless static model, no traveling-wave or device-field solution"}
    (ROOT/"data/simulated/day01_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    plt.rcParams.update({"font.size":11,"axes.spines.top":False,"axes.spines.right":False})
    fig, axes = plt.subplots(1,2,figsize=(11.4,4.3),layout="constrained")
    axes[0].plot(v,frame.T_single_arm,label="Single active arm",color="#195f9c",lw=2)
    axes[0].plot(v,frame.T_push_pull,label="Ideal push-pull (gap voltage)",color="#e68b2e",ls="--")
    axes[0].axvline(vpi,color="gray",ls=":")
    axes[0].set_title(f"Interference switching | single-arm Vpi = {vpi:.3f} V")
    axes[0].set_xlabel("Voltage across electrode gap (V)")
    axes[0].set_ylabel("Normalized optical power")
    axes[0].legend(loc="lower left",fontsize=9)
    small = np.linspace(-0.4*vpi,0.4*vpi,301)
    axes[1].plot(small/vpi,transmission(small,p,bias_rad=np.pi/2),lw=2,color="#195f9c",label="Exact interference")
    axes[1].plot(small/vpi,0.5+np.pi*small/(2*vpi),ls="--",color="#e68b2e",label="Small-signal approximation")
    axes[1].scatter([0],[0.5],color="#195f9c",s=30,zorder=3)
    axes[1].set_title("Quadrature bias: phase becomes power change")
    axes[1].set_xlabel("Applied voltage / single-arm Vpi")
    axes[1].set_ylabel("Normalized optical power")
    axes[1].legend(fontsize=9)
    for ax in axes:
        ax.grid(alpha=0.2)
        ax.set_ylim(-0.15,1.15)
    fig.suptitle("Day 1 | Ideal static MZM simulation - lecture parameters",fontsize=14)
    fig.savefig(ROOT/"figures/day01_static_mzm.png",dpi=180)
    plt.close(fig)
    print(json.dumps(summary,indent=2))
    return summary


if __name__ == "__main__":
    main()
