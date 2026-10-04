"""A07: C-06 bandwidth-definition diagnostic on the A04 L1 provider. See A07 BOUNDARY.

Same inputs and conventions as A04; no fitting. Definitions D1-D3 are applied
identically to the L1 model (50 ohm) and to the A03 Fig.3(b) 50-ohm readout.
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
from tfln_mzm.response import eo_response_db,bandwidth_3db
from tfln_mzm.termination import s11

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'experiments/track_A_reproduction/A07_c06_bandwidth_definition'
DATA=ROOT/'data/digitized/liu2022/curves.csv'
F0_GHZ=1e-3;LOAD=50.0;VREF=LOAD/(50+LOAD);WINDOW_GHZ=10.0;MIN_COVER=.8;GAP_GHZ=1.0
FIG2B_STAR_GHZ=202.5996984545797  # A06 results.json, star_bw3dB_50ohm_GHz
EXPLAIN_GHZ=195.0


def series(data,panel,label=None):
    t=data[data.panel==panel]
    if label is not None:t=t[t.series==label]
    return t.sort_values('frequency_GHz')


def a04_mae_check(data):
    """Recompute the archived A04 fig3b 50-ohm MAE with A04's own steps."""
    t=series(data,'fig3b','50');f=t.frequency_GHz.to_numpy();y=t.value.to_numpy()
    (loss,ml),(z0,mz)=[interpolate_supported(f,d.frequency_GHz,d.value,1) for d in (series(data,'fig2a'),series(data,'fig2c'))]
    s=ml&mz
    v=average_voltage(f[s]*1e9,length_m=.006,z0=z0[s],zg=50,zl=LOAD,n_m=2.25,n_g=2.25,alpha_np_m=db_per_mm_to_np_per_m(loss[s]))
    p=np.full(len(f),np.nan);p[s]=eo_response_db(v,VREF)
    gaps=np.diff(f)>1;ok=np.ones(len(f),bool);ok[:-1]&=~gaps;ok[1:]&=~gaps
    scored=s&ok&np.isfinite(p)
    return float(np.mean(abs(p[scored]-y[scored])))


def model_response(data,step):
    f=np.arange(step,200+1e-9,step)
    (loss,ml),(z0,mz)=[interpolate_supported(f,d.frequency_GHz,d.value,GAP_GHZ) for d in (series(data,'fig2a'),series(data,'fig2c'))]
    s=ml&mz;f=f[s]
    v=average_voltage(f*1e9,length_m=.006,z0=z0[s],zg=50,zl=LOAD,n_m=2.25,n_g=2.25,alpha_np_m=db_per_mm_to_np_per_m(loss[s]))
    m=eo_response_db(v,VREF)
    return f,m,int((~np.isfinite(m)).sum()),int(len(s)-s.sum())


def with_reference(f,m):
    return np.concatenate([[F0_GHZ],f]),np.concatenate([[0.0],m])


def smoothed(f,m):
    out=np.full(len(f),np.nan)
    for i,fc in enumerate(f):
        w=(f>=fc-WINDOW_GHZ/2)&(f<=fc+WINDOW_GHZ/2)
        fw=f[w];d=np.diff(fw)
        cover=(d[d<=GAP_GHZ].sum())/WINDOW_GHZ if len(fw)>1 else 0
        if cover>=MIN_COVER:out[i]=m[w].mean()
    return out


def definitions(f,m):
    """f (GHz) excludes f0; returns D1-D3 in GHz. bandwidth_3db is unit-agnostic;
    its *_hz keys therefore hold GHz here."""
    def event_info(b):
        if b['censored']:return {'value_GHz':None,'censored':True,'lower_bound_GHz':b['lower_bound_hz']}
        e=b['events'][0];br=e['bracket_hz']
        return {'value_GHz':e['frequency_hz'],'censored':False,'bracket_GHz':br,'bracket_in_gap':(br[1]-br[0])>GAP_GHZ,'kind':e['kind']}
    ff,mm=with_reference(f,m)
    b=bandwidth_3db(ff,mm)
    d1=event_info(b)|{'event_count':len(b['events'])}
    if mm[-1]>-3:d2={'value_GHz':None,'censored':True,'lower_bound_GHz':float(ff[-1]),'final_dB':float(mm[-1])}
    else:
        e=b['events'][-1];d2={'value_GHz':e['frequency_hz'],'censored':False,'bracket_GHz':e['bracket_hz'],'final_dB':float(mm[-1])}
    sm=smoothed(f,m);ok=np.isfinite(sm)
    fs,ms=with_reference(f[ok],sm[ok])
    b3=bandwidth_3db(fs,ms)
    d3=event_info(b3)|{'event_count':len(b3['events']),'smoothed_min_dB':float(np.min(sm[ok])),'smoothed_min_at_GHz':float(f[ok][np.argmin(sm[ok])])}
    return {'D1_first':d1,'D2_stays_below':d2,'D3_smoothed_first':d3},sm


def explains(model,target):
    def high(d):return d['censored'] or (d['value_GHz'] is not None and d['value_GHz']>=EXPLAIN_GHZ)
    d1_low=not high(model['D1_first']) and not high(target['D1_first'])
    return {k:bool(d1_low and high(model[k]) and high(target[k])) for k in ('D2_stays_below','D3_smoothed_first')}


def run():
    data=pd.read_csv(DATA,dtype={'series':str})
    mae=a04_mae_check(data)
    archived=[r for r in json.loads((ROOT/'experiments/track_A_reproduction/A04_digitized_comparison/results.json').read_text())['results'] if r['panel']=='fig3b' and r['load']=='50'][0]['MAE']
    grids={}
    for step in (0.05,0.025):
        f,m,nonfinite,unsupported=model_response(data,step)
        defs,sm=definitions(f,m)
        grids[str(step)]={'definitions':defs,'nonfinite':nonfinite,'unsupported_grid_points':unsupported,'f':f,'m':m,'sm':sm}
    t=series(data,'fig3b','50');tf=t.frequency_GHz.to_numpy();tv=t.value.to_numpy();tu=t.uncertainty_value.to_numpy()
    target,tsm=definitions(tf,tv)
    sens={'plus_uncertainty':definitions(tf,tv+tu)[0]['D1_first'],'minus_uncertainty':definitions(tf,tv-tu)[0]['D1_first']}
    coarse,fine=grids['0.05']['definitions'],grids['0.025']['definitions']
    def shift(k):
        a,b=coarse[k],fine[k]
        if a['censored'] or b['censored']:return None if a['censored']==b['censored'] else float('inf')
        return abs(a['value_GHz']-b['value_GHz'])
    criteria={'1_A04_mae_reproduced':abs(mae-archived)<1e-9,
        '2_grid_events_equal':coarse['D1_first']['event_count']==fine['D1_first']['event_count'] and coarse['D3_smoothed_first']['event_count']==fine['D3_smoothed_first']['event_count'],
        '2_grid_shift_le_0.1GHz':{k:shift(k) for k in coarse},
        '3_nonfinite_zero':grids['0.05']['nonfinite']==0 and grids['0.025']['nonfinite']==0,
        '4_definition_explains_C06':explains(fine,target)}
    shifts=[v for k,v in criteria['2_grid_shift_le_0.1GHz'].items() if v is not None and not fine[k].get('bracket_in_gap',False)]
    criteria['2_grid_shift_pass']=all(v<=.1 for v in shifts)
    result={'grade':'DIAGNOSTIC_ONLY','scope':'L1 provider of A04, 50 ohm only; bandwidth definitions D1-D3','fig2b_star_GHz':FIG2B_STAR_GHZ,
        'explain_threshold_GHz':EXPLAIN_GHZ,'A04_mae_recomputed_dB':mae,'A04_mae_archived_dB':archived,
        'model_fine_grid':fine,'model_coarse_grid':coarse,'target_fig3b_50':target,'target_readout_sensitivity_D1':sens,'criteria':criteria}
    (OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    g=grids['0.025']
    fig,ax=plt.subplots(figsize=(10,5.5),layout='constrained')
    ax.plot(tf,tv,'.',ms=2,color='0.55',label='Fig.3(b) 50 ohm readout (A03)')
    ax.plot(g['f'],g['m'],lw=.9,color='k',label='L1 model, A04 inputs')
    ax.plot(g['f'],g['sm'],lw=1.6,color='tab:blue',label='L1 model, 10 GHz smoothed (D3)')
    ax.plot(tf,tsm,lw=1.2,ls='--',color='tab:orange',label='readout, 10 GHz smoothed (D3)')
    ax.axhline(-3,color='r',lw=.8,ls=':');ax.axvline(FIG2B_STAR_GHZ,color='tab:red',lw=1,ls='--',label='Fig.2(b) star, 202.6 GHz')
    d1=fine['D1_first']['value_GHz']
    if d1:ax.axvline(d1,color='k',lw=.8,ls='-.',label=f'model D1 = {d1:.1f} GHz')
    ax.set_xlim(0,205);ax.set_ylim(-5,1.5);ax.set_xlabel('Frequency (GHz)');ax.set_ylabel('EO response re 1 MHz (dB)');ax.grid(alpha=.25);ax.legend(fontsize=8,loc='lower left')
    ax.set_title('A07 | 50-ohm bandwidth under definitions D1/D3; diagnostic only, no fit')
    fig.savefig(OUT/'bandwidth_definitions.png',dpi=150);plt.close(fig)
    return result


if __name__=='__main__':
    r=run()
    print(json.dumps({k:r[k] for k in ('A04_mae_recomputed_dB','model_fine_grid','target_fig3b_50','target_readout_sensitivity_D1','criteria')},indent=1))
