"""A10: is Fig.3(c)'s 40/80-ohm labelling swapped? Diagnostic only. See A10 BOUNDARY.

T1 reuses A05's phase-independent bounds (imported, not re-derived); T2 repeats A04's
fig3c comparison steps. Readout data and legends are never modified.
"""
from pathlib import Path
import json,sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parent))
from run_a05 import bounds
from tfln_mzm.digitized import interpolate_supported
from tfln_mzm.microwave import propagation,db_per_mm_to_np_per_m
from tfln_mzm.response import magnitude_db
from tfln_mzm.termination import s11

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'experiments/track_A_reproduction/A10_c07_fig3c_legend_swap'
ASSIGNMENTS={'original':{'40':40.0,'80':80.0},'swapped':{'40':80.0,'80':40.0}}
LOW_GHZ=3.0


def curve(d,panel,label=None):
    t=d[d.panel==panel]
    if label is not None:t=t[t.series==label]
    return t.sort_values('frequency_GHz')


def phase_bounds(d,label,load):
    """A05 section-2 loop, with the readout curve `label` tested against termination `load`."""
    zs=curve(d,'fig2c');ls=curve(d,'fig2a');t=curve(d,'fig3c',label);f=t.frequency_GHz.to_numpy()
    z,okz=interpolate_supported(f,zs.frequency_GHz,zs.value,1);loss,okl=interpolate_supported(f,ls.frequency_GHz,ls.value,1)
    rows=[]
    for i in np.flatnonzero(okz&okl):
        low,high=bounds(z[i],loss[i],load);err=t.uncertainty_value.iloc[i]
        targetlo=10**((t.value.iloc[i]-err)/20);targethi=10**((t.value.iloc[i]+err)/20)
        radius=t.uncertainty_frequency_GHz.iloc[i]+max(zs.uncertainty_frequency_GHz.max(),ls.uncertainty_frequency_GHz.max())
        zz=zs[abs(zs.frequency_GHz-f[i])<=radius];ll=ls[abs(ls.frequency_GHz-f[i])<=radius]
        if len(zz)==0 or len(ll)==0:continue
        zmin=float(min((zz.value-zz.uncertainty_value).min(),z[i]));zmax=float(max((zz.value+zz.uncertainty_value).max(),z[i]))
        lmin=max(0,float(min((ll.value-ll.uncertainty_value).min(),loss[i])));lmax=float(max((ll.value+ll.uncertainty_value).max(),loss[i]))
        bl,bh=bounds(np.linspace(zmin,zmax,101)[:,None],np.linspace(lmin,lmax,21)[None,:],load)
        rows.append({'frequency_GHz':float(f[i]),'target_dB':float(t.value.iloc[i]),'low_dB':float(20*np.log10(low)),'high_dB':float(20*np.log10(high)),
                     'outside_nominal':bool(targethi<low or targetlo>high),'outside_envelope':bool(targethi<float(bl.min()) or targetlo>float(bh.max()))})
    return pd.DataFrame(rows)


def a04_metrics(d,label,load):
    """A04 fig3c steps for readout `label` predicted with termination `load`."""
    t=curve(d,'fig3c',label);f=t.frequency_GHz.to_numpy();y=t.value.to_numpy()
    (loss,ml),(z0,mz)=[interpolate_supported(f,s.frequency_GHz,s.value,1) for s in (curve(d,'fig2a'),curve(d,'fig2c'))]
    s=ml&mz;p=np.full(len(f),np.nan);a=db_per_mm_to_np_per_m(loss[s])
    p[s]=magnitude_db(s11(propagation(f[s]*1e9,2.25,a),.006,z0[s],load))
    slope=np.gradient(y,f);nb=np.ones(len(f),bool);g=np.diff(f)>1;nb[:-1]&=~g;nb[1:]&=~g
    allowance=t.uncertainty_value.to_numpy()+abs(slope)*t.uncertainty_frequency_GHz.to_numpy()
    sc=s&nb&np.isfinite(p);e=p[sc]-y[sc]
    return {'scored':int(sc.sum()),'MAE_dB':float(np.mean(abs(e))),'within_readout_fraction':float(np.mean(abs(e)<=allowance[sc])),
            'nonfinite':int((s&~np.isfinite(p)).sum())},(f,p)


def run():
    d=pd.read_csv(ROOT/'data/digitized/liu2022/curves.csv',dtype={'series':str})
    a05={r['load']:r for r in json.loads((ROOT/'experiments/track_A_reproduction/A05_convention_audit/results.json').read_text())['phase_bounds']}
    a04={r['load']:r['MAE'] for r in json.loads((ROOT/'experiments/track_A_reproduction/A04_digitized_comparison/results.json').read_text())['results'] if r['panel']=='fig3c'}
    res={};curves={}
    for name,assign in ASSIGNMENTS.items():
        res[name]={}
        for label,load in assign.items():
            pb=phase_bounds(d,label,load);m,cv=a04_metrics(d,label,load);low=pb[pb.frequency_GHz<=LOW_GHZ]
            res[name][label]={'assumed_load_ohm':load,'points':int(len(pb)),'outside_nominal':int(pb.outside_nominal.sum()),
                'outside_envelope':int(pb.outside_envelope.sum()),'outside_nominal_fraction':float(pb.outside_nominal.mean()),
                'low_freq_points':int(len(low)),'low_freq_outside_nominal':int(low.outside_nominal.sum()),'A04_metrics':m,
                'all_finite':bool(np.all(np.isfinite(pb[['low_dB','high_dB']].to_numpy())) and m['nonfinite']==0)}
            curves[(name,label)]=(pb,cv)
    o,s=res['original'],res['swapped']
    g1=all(o[l]['points']==a05[l]['points'] and o[l]['outside_nominal']==a05[l]['outside_nominal'] and o[l]['outside_envelope']==a05[l]['outside_sampled_uncertainty_envelope'] for l in ('40','80')) \
        and all(abs(o[l]['A04_metrics']['MAE_dB']-a04[l])<1e-9 for l in ('40','80'))
    supported=all(s[l]['outside_nominal_fraction']<.1 and o[l]['outside_nominal_fraction']>.5 and s[l]['A04_metrics']['MAE_dB']<.5*o[l]['A04_metrics']['MAE_dB'] for l in ('40','80'))
    criteria={'1_reproduces_A05_A04_archives':bool(g1),'2_all_finite':all(r[l]['all_finite'] for r in res.values() for l in r),'3_swap_supported':bool(supported)}
    result={'grade':'DIAGNOSTIC_ONLY','scope':'Fig.3(c) 40/80-ohm label assignment; L1 inputs of A04/A05; readout unchanged',
            'dc_limit_dB':{'40':float(20*np.log10(10/90)),'80':float(20*np.log10(30/130))},'assignments':res,'criteria':criteria}
    (OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained',sharex=True,sharey=True)
    for i,name in enumerate(ASSIGNMENTS):
        for j,label in enumerate(('40','80')):
            ax=axs[j,i];pb,(f,p)=curves[(name,label)];load=ASSIGNMENTS[name][label]
            ax.fill_between(pb.frequency_GHz,pb.low_dB,pb.high_dB,color='0.85',label='any-phase bound (nominal)')
            ax.plot(f,p,lw=.8,color='k',label=f'L1 model, ZL={load:.0f} ohm')
            ax.plot(pb.frequency_GHz,pb.target_dB,'.',ms=2,color='red' if label=='40' else 'blue',label=f'readout labelled {label} ohm')
            out=pb[pb.outside_nominal];ax.plot(out.frequency_GHz,out.target_dB,'x',ms=3,color='m',label='outside bound')
            ax.set_title(f'{name}: readout "{label}" tested as {load:.0f} ohm | outside {res[name][label]["outside_nominal"]}/{res[name][label]["points"]}',fontsize=10)
            ax.set_ylim(-50,0);ax.grid(alpha=.25);ax.legend(fontsize=7,loc='lower right')
    for ax in axs[-1]:ax.set_xlabel('Frequency (GHz)')
    for ax in axs[:,0]:ax.set_ylabel('S11 (dB)')
    fig.suptitle('A10 | C-07: original vs swapped 40/80-ohm labels in Fig.3(c); diagnostic only, readout unchanged',fontsize=12)
    fig.savefig(OUT/'legend_swap.png',dpi=130);plt.close(fig)
    return result


if __name__=='__main__':
    print(json.dumps(run(),indent=1))
