"""Generate a small paper-linked DC-limit table; not Figure 3 samples."""
from pathlib import Path
import json
import numpy as np
from tfln_mzm.interaction import average_voltage
from tfln_mzm.termination import reflection_coefficients

root=Path(__file__).resolve().parents[1]
rows=[]
for load in [20,40,50,80,np.inf]:
    voltage=average_voltage(0,length_m=.006,z0=50,zg=50,zl=load,n_m=2.25,n_g=2.25,alpha_np_m=0)
    expected=1 if np.isinf(load) else load/(50+load)
    rows.append({'load_ohm':'open' if np.isinf(load) else load,
        'rho2':float(reflection_coefficients(50,50,load)[1]),
        'Vavg_dc_V':float(voltage.real),'imaginary_residual_V':float(voltage.imag),
        'divider_reference_V':expected,'absolute_error_V':float(abs(voltage-expected))})
result={'scope':'A01_DC_limit_only','provenance':'simulation_with_declared_idealizations',
        'parameters':{'Vg_V':1,'Zg_ohm':50,'Z0_ohm':50,'length_m':.006,'alpha_Np_per_m':0,'n_m':2.25,'n_g':2.25,'frequency_Hz':0},
        'source_note':'Vg/Zg/L/loads from paper; constant Z0, zero loss and matching indices are ideal-limit assumptions; f=0 is a limit test, not a change to the Figure 3 f0 convention',
        'rows':rows}
(root/'experiments/A01_analytic_limits/dc_results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
