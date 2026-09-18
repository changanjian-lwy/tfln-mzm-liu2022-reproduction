"""A02 assumed-input mechanism study, not digitized Figure 3 reproduction."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from tfln_mzm.interaction import average_voltage
from tfln_mzm.microwave import propagation,db_per_mm_to_np_per_m
from tfln_mzm.termination import s11
from tfln_mzm.response import eo_response_db,magnitude_db,bandwidth_3db

OUT=Path(__file__).resolve().parents[1]/'experiments/track_A_reproduction/A02_frequency_response'
CONFIG={'length_m':.006,'z0':50,'zg':50,'n_m':2.25,'n_g':2.25,'vg':1}
F0=1e6

def alpha(f):return db_per_mm_to_np_per_m(.05*np.sqrt(np.asarray(f)/1e9))

def run():
    all_rows=[];summaries=[]
    fig,axes=plt.subplots(3,1,figsize=(9,10),layout='constrained',sharex=True)
    for load,color in zip([20,40,50,80,np.inf],['#d98b00','#d54b40','#262626','#3e70b5','#4c973e']):
        label='open' if np.isinf(load) else str(load)+' ohm'
        results=[]
        for count in [2001,4001]:
            f=np.linspace(F0,200e9,count)
            v=average_voltage(f,**CONFIG,zl=load,alpha_np_m=alpha(f))
            ref=average_voltage(F0,**CONFIG,zl=load,alpha_np_m=alpha(F0))
            response=eo_response_db(v,ref)
            results.append(bandwidth_3db(f,response))
        a,b=results
        same=a['status']==b['status'] and len(a['events'])==len(b['events'])
        shift=None if b['censored'] else abs(a['bandwidth_hz']-b['bandwidth_hz']) if a['bandwidth_hz'] is not None and b['bandwidth_hz'] is not None else None
        converged=same and (b['censored'] or (shift is not None and shift<=.1e9))
        reflection=s11(propagation(f,2.25,alpha(f)),.006,50,load)
        sd=magnitude_db(reflection)
        all_rows.append(pd.DataFrame({'f_Hz':f,'load':label,'Vavg_real_V':v.real,'Vavg_imag_V':v.imag,'Vavg_abs_V':abs(v),'EO_dB':response,'S11_dB':sd,'provenance':'simulation_assumed_provider'}))
        summaries.append({'load':label,'reference_voltage_abs_V':float(abs(ref)),'bandwidth':b,'refinement_shift_Hz':shift,'grid_check_pass':converged})
        axes[0].plot(f/1e9,abs(v),label=label,color=color)
        axes[1].plot(f/1e9,response,label=label,color=color)
        if load!=50:axes[2].plot(f/1e9,sd,label=label,color=color)
    for ax,ylabel in zip(axes,['Average drive magnitude (V)','EO response relative to 1 MHz (dB)','Input S11 (dB)']):
        ax.set_ylabel(ylabel);ax.grid(alpha=.2);ax.legend(ncol=3,fontsize=9)
    axes[1].axhline(-3,color='gray',ls='--',lw=1)
    axes[2].text(.02,.07,'50 ohm: S11 = 0, hence -infinity dB (not plotted)',transform=axes[2].transAxes,fontsize=9)
    axes[2].set_xlabel('Frequency (GHz)')
    fig.suptitle('A02 | Assumed-input simulation, NOT Figure 3 reproduction\nZ0 = 50 ohm; matched velocity; loss = 0.05 sqrt(f/GHz) dB/mm',fontsize=12)
    fig.savefig(OUT/'response.png',dpi=160);plt.close(fig)
    pd.concat(all_rows).to_csv(OUT/'curves.csv',index=False)
    result={'grade':'QUALITATIVE_ONLY','provider_evidence':'assumption','parameters':CONFIG,'f0_Hz':F0,'loss_law':'0.05 sqrt(f/GHz) dB/mm','frequency_range_Hz':[F0,200e9],'grid_points':[2001,4001],'results':summaries}
    (OUT/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    assert all(r['grid_check_pass'] for r in summaries),'Grid convergence failed; see results.json'
    for r in summaries:print(r['load'],r['bandwidth']['bandwidth_hz'],r['bandwidth']['censored'],r['refinement_shift_Hz'])

if __name__=='__main__':run()
