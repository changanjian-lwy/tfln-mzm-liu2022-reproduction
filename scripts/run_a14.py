"""A14: input readout-uncertainty corners for the A04/A13 series. Diagnostic only. See A14 BOUNDARY.

Corners shift the whole Fig.2(c) readout by +-u_Z and the whole Fig.2(a) readout by
+-u_alpha (per-sample uncertainties; loss clipped at 0): systematic, fully correlated
calibration-type errors. Scoring points and target allowances are A04/A13's.
"""
from pathlib import Path
import itertools,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from tfln_mzm.digitized import interpolate_supported
from tfln_mzm.interaction import average_voltage
from tfln_mzm.microwave import propagation,db_per_mm_to_np_per_m
from tfln_mzm.response import magnitude_db,eo_response_db
from tfln_mzm.termination import s11

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'experiments/track_A_reproduction/A14_input_uncertainty_corners'
CORNERS=[(0,0)]+list(itertools.product((-1,1),(-1,1)))  # (Z0 shift sign, loss shift sign)


def predict(src,panel,label,f,sz,sa):
    loss_src,z_src=src
    loss,ml=interpolate_supported(f,loss_src.frequency_GHz,np.maximum(loss_src.value+sa*loss_src.uncertainty_value,0),1)
    z0,mz=interpolate_supported(f,z_src.frequency_GHz,z_src.value+sz*z_src.uncertainty_value,1)
    s=ml&mz;p=np.full(len(f),np.nan);load=np.inf if label=='open' else float(label);a=db_per_mm_to_np_per_m(loss[s]);fh=f[s]*1e9
    if panel=='fig3c':p[s]=magnitude_db(s11(propagation(fh,2.25,a),.006,z0[s],load))
    else:
        v=average_voltage(fh,length_m=.006,z0=z0[s],zg=50,zl=load,n_m=2.25,n_g=2.25,alpha_np_m=a)
        p[s]=eo_response_db(v,1 if label=='open' else load/(50+load)) if panel=='fig3b' else abs(v)
    return p,s


def assess(src,panel,label,f,y,uf,uv):
    o=np.argsort(f);f,y,uf,uv=f[o],y[o],uf[o],uv[o]
    preds=[];s=None
    for sz,sa in CORNERS:
        p,s=predict(src,panel,label,f,sz,sa);preds.append(p)
    nominal=preds[0];lo=np.min(preds,axis=0);hi=np.max(preds,axis=0)
    slope=np.gradient(y,f);nb=np.ones(len(f),bool);g=np.diff(f)>1;nb[:-1]&=~g;nb[1:]&=~g
    allowance=uv+abs(slope)*uf;sc=s&nb&np.isfinite(nominal)
    e=nominal[sc]-y[sc];orig=float(np.mean(abs(e)<=allowance[sc]))
    overlap=(y[sc]+allowance[sc]>=lo[sc])&(y[sc]-allowance[sc]<=hi[sc])
    finite=bool(all(np.all(np.isfinite(p[s])) for p in preds))
    return {'panel':panel,'load':label,'scored':int(sc.sum()),'MAE_nominal':float(np.mean(abs(e))),'within_original':orig,
            'within_combined':float(np.mean(overlap)),'median_envelope_width':float(np.median(hi[sc]-lo[sc])),'all_finite':finite},(f,y,nominal,lo,hi,sc)


def run():
    d=pd.read_csv(ROOT/'data/digitized/liu2022/curves.csv',dtype={'series':str})
    a6=pd.read_csv(ROOT/'data/digitized/liu2022/a06_curves.csv',dtype={'series':str})
    src=[d[d.panel==p].sort_values('frequency_GHz') for p in ('fig2a','fig2c')]
    a04={(r['panel'],r['load']):(r['MAE'],r['diagnostic_pass']) for r in json.loads((ROOT/'experiments/track_A_reproduction/A04_digitized_comparison/results.json').read_text())['results']}
    a13={('fig3a_main',r['load']):(r['MAE_V'],r['diagnostic_pass']) for r in json.loads((ROOT/'experiments/track_A_reproduction/A13_fig3a_main_comparison/results.json').read_text())['results']}
    archived={**a04,**a13};rows=[];plots={}
    for panel in ('fig3a_inset','fig3b','fig3c'):
        for label,t in d[d.panel==panel].groupby('series',sort=False):
            r,pl=assess(src,panel,label,t.frequency_GHz.to_numpy(),t.value.to_numpy(),t.uncertainty_frequency_GHz.to_numpy(),t.uncertainty_value.to_numpy())
            rows.append(r);plots[(panel,label)]=pl
    for label in ('20','40','open'):
        t=a6[(a6.panel=='fig3a_main')&(a6.series==label)]
        r,pl=assess(src,'fig3a_main',label,t.x_value.to_numpy(),t.value.to_numpy(),t.uncertainty_x.to_numpy(),t.uncertainty_value.to_numpy())
        rows.append(r);plots[('fig3a_main',label)]=pl
    for r in rows:
        mae,passed=archived[(r['panel'],r['load'])];r['archived_pass']=passed;r['mae_reproduced']=abs(r['MAE_nominal']-mae)<1e-9
        r['verdict']='passed_originally' if passed else ('failure_within_input_readout_precision' if r['within_combined']>=.9 else 'failure_beyond_input_readout_precision')
    criteria={'1_nominal_reproduces_A04_A13':all(r['mae_reproduced'] for r in rows),'2_all_finite':all(r['all_finite'] for r in rows),
              '3_failed_series_verdicts':{f"{r['panel']}:{r['load']}":r['verdict'] for r in rows if not r['archived_pass']},
              '4_passed_series_combined_fraction':{f"{r['panel']}:{r['load']}":r['within_combined'] for r in rows if r['archived_pass']}}
    result={'grade':'DIAGNOSTIC_ONLY','scope':'systematic input readout corners (Fig.2(a), Fig.2(c)); A04/A13 scoring','criteria':criteria,'series':rows}
    (OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    failed=[r for r in rows if not r['archived_pass']]
    n=len(failed);cols=3;rws=int(np.ceil(n/cols))
    fig,axs=plt.subplots(rws,cols,figsize=(15,3.2*rws),layout='constrained');axs=np.atleast_1d(axs).ravel()
    for ax,r in zip(axs,failed):
        f,y,nom,lo,hi,sc=plots[(r['panel'],r['load'])]
        ax.fill_between(f,lo,hi,color='0.8',label='corner envelope');ax.plot(f,nom,lw=.7,color='k',label='nominal');ax.plot(f[sc],y[sc],'.',ms=1.5,color='tab:red',label='readout')
        ax.set_title(f"{r['panel']} {r['load']}: {r['within_original']:.0%} -> {r['within_combined']:.0%}",fontsize=9);ax.grid(alpha=.25)
    for ax in axs[n:]:ax.axis('off')
    axs[0].legend(fontsize=7)
    fig.suptitle('A14 | failed series: nominal model, input-readout corner envelope and readout; diagnostic only',fontsize=12)
    fig.savefig(OUT/'corner_envelopes.png',dpi=110);plt.close(fig)
    return result


if __name__=='__main__':
    r=run();print(json.dumps({k:v for k,v in r['criteria'].items()},indent=1))
    for x in r['series']:print(x['panel'],x['load'],'orig',round(x['within_original'],3),'comb',round(x['within_combined'],3),'env',round(x['median_envelope_width'],4),x['verdict'])
