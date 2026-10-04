"""A15: A14 corners plus a frequency-axis shift of the Fig.2(c) readout. Diagnostic only.

27 combinations: Z0 vertical {-1,0,+1}u_Z x Z0 frequency {-1,0,+1}u_f x loss vertical {-1,0,+1}u_alpha.
Scoring points and allowances are A04/A13's. A14's 5 corners are recomputed as the G1 check.
"""
from pathlib import Path
import itertools,json,warnings
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
OUT=ROOT/'experiments/track_A_reproduction/A15_frequency_readout_corners'
ALL=list(itertools.product((-1,0,1),(-1,0,1),(-1,0,1)))           # (z vertical, z frequency, loss vertical)
A14=[(0,0,0)]+[(sz,0,sa) for sz in (-1,1) for sa in (-1,1)]


def predict(src,panel,label,f,sz,sf,sa):
    ls,zs=src
    loss,ml=interpolate_supported(f,ls.frequency_GHz,np.maximum(ls.value+sa*ls.uncertainty_value,0),1)
    z0,mz=interpolate_supported(f,zs.frequency_GHz+sf*zs.uncertainty_frequency_GHz,zs.value+sz*zs.uncertainty_value,1)
    s=ml&mz;p=np.full(len(f),np.nan);load=np.inf if label=='open' else float(label);a=db_per_mm_to_np_per_m(loss[s]);fh=f[s]*1e9
    if panel=='fig3c':p[s]=magnitude_db(s11(propagation(fh,2.25,a),.006,z0[s],load))
    else:
        v=average_voltage(fh,length_m=.006,z0=z0[s],zg=50,zl=load,n_m=2.25,n_g=2.25,alpha_np_m=a)
        p[s]=eo_response_db(v,1 if label=='open' else load/(50+load)) if panel=='fig3b' else abs(v)
    return p


def assess(src,panel,label,f,y,uf,uv):
    o=np.argsort(f);f,y,uf,uv=f[o],y[o],uf[o],uv[o]
    preds={c:predict(src,panel,label,f,*c) for c in ALL}
    nominal=preds[(0,0,0)];s=np.isfinite(nominal)
    slope=np.gradient(y,f);nb=np.ones(len(f),bool);g=np.diff(f)>1;nb[:-1]&=~g;nb[1:]&=~g
    allowance=uv+abs(slope)*uf;sc=s&nb
    def frac(corners):
        stack=np.array([preds[c] for c in corners])
        with warnings.catch_warnings():  # unscored frequencies may have no supported corner at all
            warnings.simplefilter('ignore',RuntimeWarning);lo=np.nanmin(stack,axis=0);hi=np.nanmax(stack,axis=0)
        return float(np.mean((y[sc]+allowance[sc]>=lo[sc])&(y[sc]-allowance[sc]<=hi[sc]))),lo,hi
    a14,_,_=frac(A14);full,lo,hi=frac(ALL)
    finite=all(np.all(np.isfinite(p[sc])|np.isnan(p[sc])) for p in preds.values())
    return {'panel':panel,'load':label,'scored':int(sc.sum()),'within_A14_corners':a14,'within_with_frequency_shift':full,'all_finite':bool(finite)},(f,y,nominal,lo,hi,sc)


def run():
    d=pd.read_csv(ROOT/'data/digitized/liu2022/curves.csv',dtype={'series':str})
    a6=pd.read_csv(ROOT/'data/digitized/liu2022/a06_curves.csv',dtype={'series':str})
    src=[d[d.panel==p].sort_values('frequency_GHz') for p in ('fig2a','fig2c')]
    a14={(r['panel'],r['load']):r for r in json.loads((ROOT/'experiments/track_A_reproduction/A14_input_uncertainty_corners/results.json').read_text())['series']}
    rows=[];plots={}
    for panel in ('fig3a_inset','fig3b','fig3c'):
        for label,t in d[d.panel==panel].groupby('series',sort=False):
            r,pl=assess(src,panel,label,t.frequency_GHz.to_numpy(),t.value.to_numpy(),t.uncertainty_frequency_GHz.to_numpy(),t.uncertainty_value.to_numpy());rows.append(r);plots[(panel,label)]=pl
    for label in ('20','40','open'):
        t=a6[(a6.panel=='fig3a_main')&(a6.series==label)]
        r,pl=assess(src,'fig3a_main',label,t.x_value.to_numpy(),t.value.to_numpy(),t.uncertainty_x.to_numpy(),t.uncertainty_value.to_numpy());rows.append(r);plots[('fig3a_main',label)]=pl
    for r in rows:
        ref=a14[(r['panel'],r['load'])];r['A14_verdict']=ref['verdict'];r['A14_reproduced']=abs(r['within_A14_corners']-ref['within_combined'])<1e-12
        if ref['verdict']=='passed_originally':r['verdict']='passed_originally'
        else:r['verdict']='failure_within_input_readout_precision' if r['within_with_frequency_shift']>=.9 else 'failure_beyond_input_readout_precision'
        r['changed_from_A14']=r['verdict']!=ref['verdict']
    criteria={'1_reproduces_A14':all(r['A14_reproduced'] for r in rows),'2_all_finite':all(r['all_finite'] for r in rows),
              '3_reassessed_A14_beyond':{f"{r['panel']}:{r['load']}":{'A14':a14[(r['panel'],r['load'])]['within_combined'],'A15':r['within_with_frequency_shift'],'verdict':r['verdict']}
                                         for r in rows if r['A14_verdict']=='failure_beyond_input_readout_precision'}}
    result={'grade':'DIAGNOSTIC_ONLY','scope':'A14 corners plus Fig.2(c) frequency-axis shift; 27 combinations','criteria':criteria,'series':rows}
    (OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    beyond=[r for r in rows if r['A14_verdict']=='failure_beyond_input_readout_precision']
    fig,axs=plt.subplots(int(np.ceil(len(beyond)/3)),3,figsize=(15,3.2*int(np.ceil(len(beyond)/3))),layout='constrained');axs=axs.ravel()
    for ax,r in zip(axs,beyond):
        f,y,nom,lo,hi,sc=plots[(r['panel'],r['load'])]
        ax.fill_between(f,lo,hi,color='0.8');ax.plot(f,nom,lw=.7,color='k');ax.plot(f[sc],y[sc],'.',ms=1.5,color='tab:red')
        ax.set_title(f"{r['panel']} {r['load']}: A14 {a14[(r['panel'],r['load'])]['within_combined']:.0%} -> A15 {r['within_with_frequency_shift']:.0%}",fontsize=9);ax.grid(alpha=.25)
    for ax in axs[len(beyond):]:ax.axis('off')
    fig.suptitle('A15 | series A14 judged beyond precision: envelope now includes Fig.2(c) frequency shift; diagnostic only',fontsize=12)
    fig.savefig(OUT/'frequency_corners.png',dpi=110);plt.close(fig)
    return result


if __name__=='__main__':
    r=run();print(json.dumps(r['criteria'],indent=1))
