"""A11: C-02 reference-frequency sensitivity on the A04 L1 provider. See A11 BOUNDARY.

Only the Eq.(4) reference changes. Below the first Fig.2(c) sample, Z0 is held at
that sample and loss comes from A08's input-side fit (both labelled assumptions).
A05's implied offsets are compared with, never used to choose or fit anything.
"""
from pathlib import Path
import json,sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parent))
from run_a08 import loss_provider
from run_a09 import series,zl,vref
from tfln_mzm.digitized import interpolate_supported
from tfln_mzm.interaction import average_voltage
from tfln_mzm.microwave import db_per_mm_to_np_per_m
from tfln_mzm.response import eo_response_db

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'experiments/track_A_reproduction/A11_c02_reference_frequency'
LOADS=['20','40','50','80','open'];CANDIDATES_GHZ=[0.0,0.5,1.0,None];TOL_DB=.15  # None = first Fig.2(c) sample


def reference(d,label,f0,alpha_db,hold_shift=0.0):
    if f0==0.0:return complex(vref(label))
    zs=series(d,'fig2c');z,ok=interpolate_supported(np.array([f0]),zs.frequency_GHz,zs.value,1)
    z0=float(z[0]) if ok[0] else float(zs.value.iloc[0])+hold_shift
    v=average_voltage(np.array([f0*1e9]),length_m=.006,z0=z0,zg=50,zl=zl(label),n_m=2.25,n_g=2.25,alpha_np_m=db_per_mm_to_np_per_m(alpha_db(np.array([f0]))))
    return complex(v[0])


def fig3b_metrics(d,label,ref):
    """A04 fig3b steps with an arbitrary reference phasor."""
    t=series(d,'fig3b',label);f=t.frequency_GHz.to_numpy();y=t.value.to_numpy()
    (loss,ml),(z0,mz)=[interpolate_supported(f,s.frequency_GHz,s.value,1) for s in (series(d,'fig2a'),series(d,'fig2c'))]
    s=ml&mz;p=np.full(len(f),np.nan)
    v=average_voltage(f[s]*1e9,length_m=.006,z0=z0[s],zg=50,zl=zl(label),n_m=2.25,n_g=2.25,alpha_np_m=db_per_mm_to_np_per_m(loss[s]))
    p[s]=eo_response_db(v,ref)
    slope=np.gradient(y,f);nb=np.ones(len(f),bool);g=np.diff(f)>1;nb[:-1]&=~g;nb[1:]&=~g
    allowance=t.uncertainty_value.to_numpy()+abs(slope)*t.uncertainty_frequency_GHz.to_numpy()
    sc=s&nb&np.isfinite(p);e=p[sc]-y[sc]
    return {'MAE_dB':float(np.mean(abs(e))),'within_readout_fraction':float(np.mean(abs(e)<=allowance[sc])),'scored':int(sc.sum())}


def run():
    d=pd.read_csv(ROOT/'data/digitized/liu2022/curves.csv',dtype={'series':str})
    alpha_db,_,fit=loss_provider(d)
    first=float(series(d,'fig2c').frequency_GHz.iloc[0])
    a05={r['load']:r['equivalent_normalization_offset_dB'] for r in json.loads((ROOT/'experiments/track_A_reproduction/A05_convention_audit/results.json').read_text())['cross_panel_reference']}
    a04={r['load']:r['MAE'] for r in json.loads((ROOT/'experiments/track_A_reproduction/A04_digitized_comparison/results.json').read_text())['results'] if r['panel']=='fig3b'}
    table={};finite=True
    for c in CANDIDATES_GHZ:
        f0=first if c is None else c;key=f'{f0:.4f}_GHz' if f0 else 'DC'
        rows={}
        for l in LOADS:
            ref=reference(d,l,f0,alpha_db);finite&=bool(np.isfinite(ref) and abs(ref)>0)
            delta=float(20*np.log10(abs(ref)/vref(l)))
            rows[l]={'offset_dB':delta,'A05_implied_offset_dB':a05[l],'difference_dB':delta-a05[l],**fig3b_metrics(d,l,ref)}
        table[key]={'f0_GHz':f0,'loads':rows,'compatible_with_A05':all(abs(r['difference_dB'])<=TOL_DB for r in rows.values())}
    # Reported-only (added after the first run): hold-value sensitivity at f0'=1 GHz, +-2x median Z0 readout uncertainty.
    u=2*float(series(d,'fig2c').uncertainty_value.median())
    hold={f'{s:+.2f}_ohm':{l:float(20*np.log10(abs(reference(d,l,1.0,alpha_db,s))/vref(l))) for l in LOADS} for s in (-u,u)}
    criteria={'1_DC_reproduces_A04':all(abs(table['DC']['loads'][l]['MAE_dB']-a04[l])<1e-9 for l in LOADS),'2_finite_nonzero':finite,
              '3_compatible_candidates':[k for k,v in table.items() if v['compatible_with_A05']]}
    result={'grade':'SENSITIVITY_ONLY','scope':'Eq.(4) reference only; L1 provider of A04; low-frequency hold below the first Z0 sample',
            'first_Z0_sample_GHz':first,'first_Z0_sample_ohm':float(series(d,'fig2c').value.iloc[0]),'reported_only_hold_sensitivity_1GHz_offsets_dB':hold,'loss_fit':fit,'tolerance_dB':TOL_DB,'criteria':criteria,'candidates':table}
    (OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    fig,axs=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    x=np.arange(len(LOADS))
    for i,(k,v) in enumerate(table.items()):
        axs[0].plot(x,[v['loads'][l]['offset_dB'] for l in LOADS],'o-',label=f"model, f0'={k}")
        axs[1].plot(x,[v['loads'][l]['MAE_dB'] for l in LOADS],'o-',label=k)
    axs[0].errorbar(x,[a05[l] for l in LOADS],yerr=TOL_DB,fmt='ks',capsize=4,label='A05 implied offset ±0.15 dB')
    for ax,yl in zip(axs,['offset vs DC divider (dB)','Fig.3(b) MAE (dB)']):ax.set_xticks(x,LOADS);ax.set_xlabel('load');ax.set_ylabel(yl);ax.grid(alpha=.25);ax.legend(fontsize=8)
    fig.suptitle("A11 | C-02 reference frequency: fixed candidates only, no search; sensitivity",fontsize=12)
    fig.savefig(OUT/'reference_frequency.png',dpi=140);plt.close(fig)
    return result


if __name__=='__main__':
    r=run()
    for k,v in r['candidates'].items():print(k,v['compatible_with_A05'],{l:(round(x['offset_dB'],3),round(x['A05_implied_offset_dB'],3),round(x['MAE_dB'],3)) for l,x in v['loads'].items()})
    print(r['criteria'])
