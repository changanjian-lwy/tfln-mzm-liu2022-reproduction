"""A13: A04 criteria applied to the A06 Fig.3(a) main-panel readout (50-200 GHz). See A13 BOUNDARY."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from tfln_mzm.digitized import interpolate_supported
from tfln_mzm.interaction import average_voltage
from tfln_mzm.microwave import db_per_mm_to_np_per_m

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'experiments/track_A_reproduction/A13_fig3a_main_comparison'
SERIES=['20','40','open']


def score_voltage(curves,f,y,uf,uv,label):
    """A04 steps for |Vavg| targets; returns summary and prediction."""
    order=np.argsort(f);f,y,uf,uv=f[order],y[order],uf[order],uv[order]
    inputs=[curves[curves.panel==p].sort_values('frequency_GHz') for p in ('fig2a','fig2c')]
    (loss,ml),(z0,mz)=[interpolate_supported(f,d.frequency_GHz,d.value,1) for d in inputs]
    s=ml&mz;p=np.full(len(f),np.nan);load=np.inf if label=='open' else float(label)
    p[s]=abs(average_voltage(f[s]*1e9,length_m=.006,z0=z0[s],zg=50,zl=load,n_m=2.25,n_g=2.25,alpha_np_m=db_per_mm_to_np_per_m(loss[s])))
    slope=np.gradient(y,f);nb=np.ones(len(f),bool);g=np.diff(f)>1;nb[:-1]&=~g;nb[1:]&=~g
    allowance=uv+abs(slope)*uf;sc=s&nb&np.isfinite(p);e=p[sc]-y[sc]
    within=float(np.mean(abs(e)<=allowance[sc])) if sc.any() else None;nonfinite=int((s&~np.isfinite(p)).sum())
    return {'load':label,'target_count':int(len(f)),'scored_count':int(sc.sum()),'unsupported_input_count':int((~s).sum()),
            'target_gap_excluded_count':int((s&~nb).sum()),'nonfinite_prediction_count':nonfinite,
            'MAE_V':float(np.mean(abs(e))) if sc.any() else None,'RMSE_V':float(np.sqrt(np.mean(e**2))) if sc.any() else None,
            'max_abs_error_V':float(np.max(abs(e))) if sc.any() else None,'mean_signed_error_V':float(np.mean(e)) if sc.any() else None,
            'within_readout_fraction':within,'diagnostic_pass':bool(sc.sum()>=100 and within is not None and within>=.9 and nonfinite==0)},(f,y,p)


def run():
    curves=pd.read_csv(ROOT/'data/digitized/liu2022/curves.csv',dtype={'series':str})
    a06=pd.read_csv(ROOT/'data/digitized/liu2022/a06_curves.csv',dtype={'series':str})
    archived=[r for r in json.loads((ROOT/'experiments/track_A_reproduction/A04_digitized_comparison/results.json').read_text())['results'] if r['panel']=='fig3a_inset' and r['load']=='40'][0]['MAE']
    t=curves[(curves.panel=='fig3a_inset')&(curves.series=='40')]
    g1,_=score_voltage(curves,t.frequency_GHz.to_numpy(),t.value.to_numpy(),t.uncertainty_frequency_GHz.to_numpy(),t.uncertainty_value.to_numpy(),'40')
    rows=[];plots={}
    for label in SERIES:
        t=a06[(a06.panel=='fig3a_main')&(a06.series==label)]
        r,pl=score_voltage(curves,t.x_value.to_numpy(),t.value.to_numpy(),t.uncertainty_x.to_numpy(),t.uncertainty_value.to_numpy(),label)
        rows.append(r);plots[label]=pl
    allpass=all(r['diagnostic_pass'] for r in rows)
    criteria={'1_inset40_reproduces_A04':abs(g1['MAE_V']-archived)<1e-9,'2_all_finite':all(r['nonfinite_prediction_count']==0 for r in rows),
              '3_per_series_pass':{r['load']:r['diagnostic_pass'] for r in rows}}
    result={'grade':'DIGITIZED_COMPARISON_PASS' if allpass else 'PARTIAL_COMPARISON_FAIL',
            'scope':'Fig.3(a) main panel 50-200 GHz, 20/40 ohm and open only; L1 provider of A04',
            'a06_curves_sha256':__import__('hashlib').sha256((ROOT/'data/digitized/liu2022/a06_curves.csv').read_bytes()).hexdigest(),
            'criteria':criteria,'results':rows}
    (OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    col={'20':'#d98b00','40':'red','open':'green'}
    fig,axs=plt.subplots(3,1,figsize=(10,9),layout='constrained',sharex=True)
    for ax,label in zip(axs,SERIES):
        f,y,p=plots[label];ax.plot(f,y,'.',ms=2,color=col[label],alpha=.5,label=f'readout {label}');ax.plot(f,p,lw=.9,color='k',label='L1 model')
        r=[x for x in rows if x['load']==label][0]
        ax.set_title(f"{label}: MAE {r['MAE_V']:.4f} V, within {r['within_readout_fraction']:.1%}, pass {r['diagnostic_pass']}",fontsize=10)
        ax.set_ylabel('|Vavg| (V)');ax.grid(alpha=.25);ax.legend(fontsize=8)
    axs[-1].set_xlabel('Frequency (GHz)')
    fig.suptitle('A13 | Fig.3(a) main panel 50-200 GHz vs L1 model, A04 criteria; no fit',fontsize=12)
    fig.savefig(OUT/'fig3a_main_comparison.png',dpi=140);plt.close(fig)
    return result


if __name__=='__main__':
    r=run();print(json.dumps(r['criteria'],indent=1))
    for x in r['results']:print(x)
