"""B02-B04: single-factor scans (length, loss multiplier, mismatch) on the frozen A04 L1 provider.

Usage: python scripts/run_b02_b04.py --study B02|B03|B04. See each experiment's BOUNDARY.md.
Only the named factor changes; ZL is fixed at 40 and 50 ohm. D3 smoothing is A07's.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parent))
from run_a07 import smoothed
from run_b01 import inputs
from tfln_mzm.interaction import average_voltage
from tfln_mzm.microwave import propagation
from tfln_mzm.response import eo_response_db,magnitude_db,bandwidth_3db
from tfln_mzm.termination import s11

ROOT=Path(__file__).resolve().parents[1]
STUDIES={'B02':('length_mm',[3.0,6.0,9.0,12.0],'B02_length_scan'),
         'B03':('loss_multiplier',[.5,1.0,1.5,2.0],'B03_loss_scan'),
         'B04':('delta_n',[-.10,-.05,0.0,.05,.10],'B04_mismatch_scan')}
LOADS=[40.0,50.0];NG=2.25;F0=1e-3;BASE={'length_mm':6.0,'loss_multiplier':1.0,'delta_n':0.0}


def first_crossing(f,m):
    b=bandwidth_3db(np.concatenate([[F0],f]),np.concatenate([[0.0],m]))
    e=b['events'][0] if b['events'] else None
    return {'GHz':b['bandwidth_hz'],'censored':b['censored'],'bracket_GHz':e['bracket_hz'] if e else None}


def evaluate(f,a,z0,zl,factor,value):
    L=(value if factor=='length_mm' else 6.0)*1e-3
    alpha=a*(value if factor=='loss_multiplier' else 1.0)
    nm=NG+(value if factor=='delta_n' else 0.0)
    v=average_voltage(f*1e9,length_m=L,z0=z0,zg=50,zl=zl,n_m=nm,n_g=NG,alpha_np_m=alpha)
    m=eo_response_db(v,zl/(50+zl));low=f<=50
    sd=magnitude_db(s11(propagation(f*1e9,nm,alpha),L,z0,zl))
    sm=smoothed(f,m);ok=np.isfinite(sm)
    near=lambda g:float(m[np.argmin(abs(f-g))])
    out={'value':value,'ZL_ohm':zl,'P50_dB':float(max(0.0,m[low].max())-min(0.0,m[low].min())),'D1':first_crossing(f,m),
         'D3':first_crossing(f[ok],sm[ok]),'M_100GHz_dB':near(100),'M_200GHz_dB':near(200),'maxS11_dB':float(sd.max()),
         'nonfinite':int((~np.isfinite(m)).sum()+np.isnan(sd).sum())}
    if factor=='length_mm':
        vdc=zl/(50+zl)
        out['G_dB']={str(g):float(20*np.log10((value/6.0)*abs(v[np.argmin(abs(f-g))])/vdc)) for g in (50,100,200)}
    return out,m


def run(study):
    factor,values,folder=STUDIES[study];out_dir=ROOT/'experiments/track_B_extensions'/folder
    d=pd.read_csv(ROOT/'data/digitized/liu2022/curves.csv',dtype={'series':str})
    grids={}
    for step in (.025,.05):
        f,a,z0=inputs(d,np.arange(step,200+1e-9,step),False)
        grids[step]={(zl,v):evaluate(f,a,z0,zl,factor,v) for zl in LOADS for v in values}
    fine={k:r for k,(r,_) in grids[.025].items()};coarse={k:r for k,(r,_) in grids[.05].items()}
    b01=pd.read_csv(ROOT/'experiments/track_B_extensions/B01_termination_sweep/sweep.csv')
    g1=[]
    for zl in LOADS:
        r=fine[(zl,BASE[factor])];row=b01[b01.ZL_ohm.astype(str)==str(zl)].iloc[0]
        d1_same=(r['D1']['censored'] and bool(row.D1_censored)) or (r['D1']['GHz'] is not None and abs(r['D1']['GHz']-row.D1_GHz)<1e-9)
        g1.append(abs(r['P50_dB']-row.P50_dB)<1e-9 and d1_same)
    dp=max(abs(fine[k]['P50_dB']-coarse[k]['P50_dB']) for k in fine)
    def d1shift(k):
        a,b=fine[k]['D1'],coarse[k]['D1']
        if a['censored']!=b['censored']:return float('inf')
        if a['censored'] or a['bracket_GHz'] is None or a['bracket_GHz'][1]-a['bracket_GHz'][0]>1:return 0.0
        return abs(a['GHz']-b['GHz'])
    dd=max(d1shift(k) for k in fine)
    criteria={'1_baseline_matches_B01':all(g1),'2_P50_grid_shift_max_dB':float(dp),'2_D1_grid_shift_max_GHz':float(dd),
              '2_pass':bool(dp<=.02 and dd<=.1),'3_finite':all(r['nonfinite']==0 for r in fine.values())}
    rows=list(fine.values())
    result={'grade':'SENSITIVITY_ONLY','track':'B','study':study,'factor':factor,'values':values,'fixed_loads_ohm':LOADS,
            'parent':'A04 L1 provider','criteria':criteria,'results':rows}
    (out_dir/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    fig,axs=plt.subplots(1,2,figsize=(12,4.5),layout='constrained',sharey=True)
    for ax,zl in zip(axs,LOADS):
        for v in values:
            _,m=grids[.025][(zl,v)];f=inputs(d,np.arange(.025,200+1e-9,.025),False)[0]
            ax.plot(f,m,lw=.8,label=f'{factor}={v:g}')
        ax.axhline(-3,color='r',ls=':',lw=1);ax.set_title(f'ZL = {zl:.0f} ohm',fontsize=10);ax.set_xlabel('Frequency (GHz)');ax.grid(alpha=.25);ax.legend(fontsize=8)
    axs[0].set_ylabel('EO response re 1 MHz (dB)')
    fig.suptitle(f'{study} | {factor} scan on the A04 L1 provider; model sensitivity, not measurement',fontsize=12)
    fig.savefig(out_dir/f'{folder}.png',dpi=130);plt.close(fig)
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--study',choices=sorted(STUDIES),required=True)
    r=run(ap.parse_args().study);print(json.dumps(r['criteria']))
    for x in r['results']:
        print(x['ZL_ohm'],x['value'],'P50',round(x['P50_dB'],3),'D1',x['D1']['GHz'] and round(x['D1']['GHz'],1),x['D1']['censored'],'D3',x['D3']['GHz'] and round(x['D3']['GHz'],1),x['D3']['censored'],
              'M100',round(x['M_100GHz_dB'],2),'M200',round(x['M_200GHz_dB'],2),'S11max',round(x['maxS11_dB'],1),x.get('G_dB',''))
