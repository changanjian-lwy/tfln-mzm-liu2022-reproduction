"""A09: decidable checks of the paper's qualitative statements S1-S3. See A09 BOUNDARY.

Each check is applied to the A03 readout of the author's Fig.3(b,c) and to the L1
model with A04's inputs and conventions. No fitting.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from tfln_mzm.digitized import interpolate_supported
from tfln_mzm.interaction import average_voltage
from tfln_mzm.microwave import propagation,db_per_mm_to_np_per_m
from tfln_mzm.response import eo_response_db,magnitude_db
from tfln_mzm.termination import s11

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'experiments/track_A_reproduction/A09_qualitative_statements'
LOADS=['20','40','50','80','open'];RISE=.3;DROP=-.3;S2_ABS=1.0;S2_REL=.2;S11_MAX=-10.0;NEAR=2.0;MIN_PTS=5


def zl(label):return np.inf if label=='open' else float(label)
def vref(label):return 1.0 if label=='open' else float(label)/(50+float(label))


def series(d,panel,label=None):
    t=d[d.panel==panel]
    if label is not None:t=t[t.series==label]
    return t.sort_values('frequency_GHz')


def model_at(d,f,label):
    (loss,ml),(z0,mz)=[interpolate_supported(f,s.frequency_GHz,s.value,1) for s in (series(d,'fig2a'),series(d,'fig2c'))]
    ok=ml&mz;f=f[ok];a=db_per_mm_to_np_per_m(loss[ok])
    v=average_voltage(f*1e9,length_m=.006,z0=z0[ok],zg=50,zl=zl(label),n_m=2.25,n_g=2.25,alpha_np_m=a)
    sd=magnitude_db(s11(propagation(f*1e9,2.25,a),.006,z0[ok],zl(label))) if label!='open' else None
    return f,eo_response_db(v,vref(label)),sd,ok


def a04_fig3b_mae(d):
    out={}
    for label in LOADS:
        t=series(d,'fig3b',label);f=t.frequency_GHz.to_numpy();y=t.value.to_numpy()
        ff,m,_,ok=model_at(d,f,label);p=np.full(len(f),np.nan);p[ok]=m
        gaps=np.diff(f)>1;nb=np.ones(len(f),bool);nb[:-1]&=~gaps;nb[1:]&=~gaps
        sc=ok&nb&np.isfinite(p);out[label]=float(np.mean(abs(p[sc]-y[sc])))
    return out


def window(f,y,lo,hi,fn):
    w=(f>lo)&(f<=hi)&np.isfinite(y)
    return float(fn(y[w])) if w.sum()>=MIN_PTS else None


def checks(eo,sp,z0ref):
    """eo/sp: {label: (f, value)}; returns per-statement results."""
    cls={l:('above' if l=='open' or zl(l)>z0ref+NEAR else 'below' if zl(l)<z0ref-NEAR else 'near') for l in LOADS}
    s1a={l:window(*eo[l],1e-3,20,np.max) for l in LOADS if cls[l]=='below'}
    s1b={l:window(*eo[l],1e-3,10,np.min) for l in LOADS if cls[l]=='above'}
    above=sorted([l for l in s1b],key=zl)
    deeper=[s1b[l] for l in above]
    s1a_ok=all(v is not None and v>=RISE for v in s1a.values())
    s1b_ok=all(v is not None and v<=DROP for v in s1b.values()) and all(b<a for a,b in zip(deeper,deeper[1:]))
    roll={l:None if (a:=window(*eo[l],180,195,np.mean)) is None or (b:=window(*eo[l],95,110,np.mean)) is None else a-b for l in LOADS}
    at10={}
    for l in LOADS:
        f,y=eo[l];v,ok=interpolate_supported(np.array([10.0]),f,y,1);at10[l]=float(v[0]) if ok[0] else None
    vals=[v for v in roll.values() if v is not None]
    spread=max(vals)-min(vals) if len(vals)==len(LOADS) else None
    low=[v for v in at10.values() if v is not None];low_spread=max(low)-min(low) if len(low)==len(LOADS) else None
    # None = undecidable (BOUNDARY: too few supported points), never silently False.
    s2_ok=None if spread is None or low_spread is None else bool(spread<=S2_ABS and spread<=S2_REL*low_spread)
    s2_abs=None if spread is None else bool(spread<=S2_ABS)
    s3={l:window(*sp[l],0,200.1,np.max) for l in sp}
    s3_ok=all(s3[l] is not None and s3[l]<=S11_MAX for l in ('40','50','80'))
    return {'classification':cls,'S1a_max_dB_0_20GHz':s1a,'S1a':s1a_ok,'S1b_min_dB_0_10GHz':s1b,'S1b_order':above,'S1b':s1b_ok,
            'S2_rolloff_dB':roll,'S2_spread_dB':spread,'S2_absolute_part':s2_abs,'S2_M_at_10GHz':at10,'S2_low_spread_dB':low_spread,'S2':s2_ok,
            'S3_max_S11_dB':s3,'S3':s3_ok}


def run():
    d=pd.read_csv(ROOT/'data/digitized/liu2022/curves.csv',dtype={'series':str})
    z=series(d,'fig2c');z0ref=float(np.median(z[z.frequency_GHz<=20].value))
    mae=a04_fig3b_mae(d)
    archived={r['load']:r['MAE'] for r in json.loads((ROOT/'experiments/track_A_reproduction/A04_digitized_comparison/results.json').read_text())['results'] if r['panel']=='fig3b'}
    grid=np.arange(.05,200+1e-9,.05)
    eo_m={};sp_m={};nonfinite={};s11_zero={}
    for l in LOADS:
        f,m,sd,_=model_at(d,grid,l);eo_m[l]=(f,m);nonfinite[l]=int((~np.isfinite(m)).sum())
        if sd is not None:sp_m[l]=(f,sd);s11_zero[l]=int(np.isneginf(sd).sum());nonfinite[l]+=int(np.isnan(sd).sum()+np.isposinf(sd).sum())
    eo_r={l:(series(d,'fig3b',l).frequency_GHz.to_numpy(),series(d,'fig3b',l).value.to_numpy()) for l in LOADS}
    sp_r={l:(series(d,'fig3c',l).frequency_GHz.to_numpy(),series(d,'fig3c',l).value.to_numpy()) for l in ('20','40','50','80')}
    readout=checks(eo_r,sp_r,z0ref);model=checks(eo_m,sp_m,z0ref)
    statements=['S1a','S1b','S2','S3']
    def verdict(s):
        if not readout[s]:return 'not_scored_readout_inconsistent'
        return 'undecidable' if model[s] is None else 'reproduced' if model[s] else 'not_reproduced'
    scored={s:verdict(s) for s in statements}
    criteria={'1_A04_mae_reproduced':all(abs(mae[l]-archived[l])<1e-9 for l in LOADS),'2_model_finite':all(v==0 for v in nonfinite.values()),
              '3_readout_consistent':{s:readout[s] for s in statements},'4_model_verdict':scored}
    result={'grade':'QUALITATIVE_ONLY','scope':'three textual statements of p.855-856; L1 provider of A04','Z0ref_ohm_median_0_20GHz':z0ref,
            'A04_fig3b_mae_recomputed':mae,'model_s11_exact_zero_points':s11_zero,'criteria':criteria,'readout':readout,'model':model}
    (OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    fig,axs=plt.subplots(2,1,figsize=(10,9),layout='constrained',sharex=True)
    col={'20':'#d98b00','40':'red','50':'black','80':'blue','open':'green'}
    for l in LOADS:
        axs[0].plot(*eo_r[l],'.',ms=1.5,color=col[l],alpha=.35);axs[0].plot(*eo_m[l],lw=.9,color=col[l],label=l)
        if l in sp_m:axs[1].plot(*sp_r[l],'.',ms=1.5,color=col[l],alpha=.35);axs[1].plot(*sp_m[l],lw=.9,color=col[l],label=l)
    for a,b in ((95,110),(180,195)):axs[0].axvspan(a,b,color='0.9',zorder=0)
    axs[0].axvspan(0,20,color='#fff3d6',zorder=0);axs[1].axhline(S11_MAX,color='r',ls=':',lw=1)
    axs[0].set_ylabel('EO response (dB)');axs[1].set_ylabel('S11 (dB)');axs[1].set_ylim(-50,0);axs[1].set_xlabel('Frequency (GHz)')
    for a in axs:a.grid(alpha=.25);a.legend(ncol=5,fontsize=8)
    fig.suptitle('A09 | dots: author readout; lines: L1 model. Shaded: S1 (0-20 GHz) and S2 windows; dotted: S3 limit',fontsize=11)
    fig.savefig(OUT/'qualitative_checks.png',dpi=140);plt.close(fig)
    return result


if __name__=='__main__':
    r=run();print(json.dumps({k:r[k] for k in ('Z0ref_ohm_median_0_20GHz','criteria','readout','model')},indent=1))
