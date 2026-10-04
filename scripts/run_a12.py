"""A12: physically consistent complex Z0 (G=0) for C-01, re-scored with A04's rules. See A12 BOUNDARY.

C-01b-Re : readout is Re Z0      -> Z0 = Zr (1 - j alpha/beta)
C-01b-Abs: readout is |Z0|       -> Z0 = Zr (1 - j alpha/beta)/sqrt(1+(alpha/beta)^2)
All other inputs, scoring points, allowances and criteria are A04's.
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
from tfln_mzm.microwave import propagation,db_per_mm_to_np_per_m,C0
from tfln_mzm.response import magnitude_db,eo_response_db
from tfln_mzm.termination import s11

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'experiments/track_A_reproduction/A12_c01_complex_z0'
NM=2.25


def build_z0(zr,alpha,f_hz,branch):
    if branch=='real':return zr
    ratio=alpha/(2*np.pi*f_hz*NM/C0)
    z=zr*(1-1j*ratio)
    return z if branch=='C-01b-Re' else z/np.sqrt(1+ratio**2)


def score_all(d,branch,zero_loss=False):
    inputs=[d[d.panel==p].sort_values('frequency_GHz') for p in ('fig2a','fig2c')]
    out=[];curves={}
    for panel in ('fig3a_inset','fig3b','fig3c'):
        for label,t in d[d.panel==panel].groupby('series',sort=False):
            t=t.sort_values('frequency_GHz');f=t.frequency_GHz.to_numpy();y=t.value.to_numpy()
            (loss,ml),(zr,mz)=[interpolate_supported(f,s.frequency_GHz,s.value,1) for s in inputs]
            s=ml&mz;load=np.inf if label=='open' else float(label)
            a=np.zeros(s.sum()) if zero_loss else db_per_mm_to_np_per_m(loss[s]);fh=f[s]*1e9
            z0=build_z0(zr[s],a,fh,branch)
            v=average_voltage(fh,length_m=.006,z0=z0,zg=50,zl=load,n_m=NM,n_g=2.25,alpha_np_m=a)
            ref=1 if label=='open' else load/(50+load)
            p=np.full(len(f),np.nan)
            p[s]=abs(v) if panel=='fig3a_inset' else eo_response_db(v,ref) if panel=='fig3b' else magnitude_db(s11(propagation(fh,NM,a),.006,z0,load))
            slope=np.gradient(y,f);nb=np.ones(len(f),bool);g=np.diff(f)>1;nb[:-1]&=~g;nb[1:]&=~g
            allowance=t.uncertainty_value.to_numpy()+abs(slope)*t.uncertainty_frequency_GHz.to_numpy()
            sc=s&nb&np.isfinite(p);e=p[sc]-y[sc];within=float(np.mean(abs(e)<=allowance[sc]))
            out.append({'panel':panel,'load':label,'scored':int(sc.sum()),'MAE':float(np.mean(abs(e))),'within_readout_fraction':within,
                        'diagnostic_pass':bool(sc.sum()>=100 and within>=.9 and int((s&~np.isfinite(p)).sum())==0),'nonfinite':int((s&~np.isfinite(p)).sum())})
            curves[(panel,label)]=(f,y,p)
    return out,curves


def run():
    d=pd.read_csv(ROOT/'data/digitized/liu2022/curves.csv',dtype={'series':str})
    a04={(r['panel'],r['load']):r['MAE'] for r in json.loads((ROOT/'experiments/track_A_reproduction/A04_digitized_comparison/results.json').read_text())['results']}
    res={};curves={}
    for b in ('real','C-01b-Re','C-01b-Abs'):res[b],curves[b]=score_all(d,b)
    zero={b:score_all(d,b,True)[1] for b in ('real','C-01b-Re','C-01b-Abs')}
    zero_diff=max(float(np.nanmax(abs(zero[b][k][2]-zero['real'][k][2]))) for b in ('C-01b-Re','C-01b-Abs') for k in zero['real'])
    base={(r['panel'],r['load']):r for r in res['real']}
    def delta(b):
        return [{'panel':r['panel'],'load':r['load'],'MAE':r['MAE'],'MAE_change':r['MAE']-base[(r['panel'],r['load'])]['MAE'],
                 'within_readout_fraction':r['within_readout_fraction'],'within_change':r['within_readout_fraction']-base[(r['panel'],r['load'])]['within_readout_fraction'],
                 'diagnostic_pass':r['diagnostic_pass']} for r in res[b]]
    branches={b:delta(b) for b in ('C-01b-Re','C-01b-Abs')}
    fig3c={b:[r for r in branches[b] if r['panel']=='fig3c'] for b in branches}
    criteria={'1_real_reproduces_A04':all(abs(r['MAE']-a04[(r['panel'],r['load'])])<1e-9 for r in res['real']),
              '1_zero_loss_branches_equal_real_max_abs':zero_diff,'1_zero_loss_pass':bool(zero_diff<1e-12),
              '2_all_finite':all(r['nonfinite']==0 for b in res.values() for r in b),
              '3_explains_S11':{b:sum(r['diagnostic_pass'] for r in fig3c[b])>=2 for b in fig3c},
              '3_material_effect':{b:any(abs(r['MAE_change'])>=1 for r in fig3c[b]) for b in fig3c},
              '3_max_abs_fig3c_MAE_change_dB':{b:max(abs(r['MAE_change']) for r in fig3c[b]) for b in fig3c}}
    result={'grade':'DIAGNOSTIC_ONLY','scope':'C-01b complex Z0 from G=0; A04 inputs, scoring and criteria','criteria':criteria,
            'baseline_real':res['real'],'branches':branches}
    (OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained',sharex=True,sharey=True)
    for ax,label in zip(axs.ravel(),['20','40','50','80']):
        f,y,_=curves['real'][('fig3c',label)]
        ax.plot(f,y,'.',ms=2,color='0.6',label='readout')
        for b,c in (('real','k'),('C-01b-Re','tab:red'),('C-01b-Abs','tab:blue')):ax.plot(f,curves[b][('fig3c',label)][2],lw=.8,color=c,label=b)
        ax.set_title(f'Fig.3(c) {label} ohm',fontsize=10);ax.set_ylim(-50,0);ax.grid(alpha=.25);ax.legend(fontsize=7)
    for ax in axs[-1]:ax.set_xlabel('Frequency (GHz)')
    for ax in axs[:,0]:ax.set_ylabel('S11 (dB)')
    fig.suptitle('A12 | C-01b: G=0 complex Z0 vs real Z0, A04 scoring; diagnostic only',fontsize=12)
    fig.savefig(OUT/'complex_z0.png',dpi=130);plt.close(fig)
    return result


if __name__=='__main__':
    r=run();print(json.dumps(r['criteria'],indent=1))
    for b,rows in r['branches'].items():
        print(b);[print(' ',x['panel'],x['load'],round(x['MAE'],4),'dMAE',round(x['MAE_change'],4),'within',round(x['within_readout_fraction'],3),'pass',x['diagnostic_pass']) for x in rows]
