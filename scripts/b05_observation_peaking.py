"""B05 post-hoc observation (written after the run; not a criterion): peaking of the J-optimal designs.

The B05 constraints bound the normalized response only from below (c1 >= -3 dB). This script reports,
for each grid optimum in front.csv and for the paper design (40 ohm, 6 mm), the maximum normalized
response up to B, the peak-to-peak ripple and the absolute response G(f) = 20log10[(L/6 mm)|Vavg|/0.5 V].
It reads B05 outputs and recomputes curves with the same provider; it changes no B05 output.
"""
from pathlib import Path
import json,sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parent))
from run_b01 import inputs
from run_b05 import line,INK,INK2,SERIES,PAPER
from tfln_mzm.interaction import average_voltage

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'experiments/track_B_extensions/B05_joint_zl_length_bo'


def curves(f,a,z0,zl,lmm):
    m,_=line(f,a,z0,np.array([[zl]]),lmm)
    v=average_voltage(f*1e9,length_m=lmm/1000,z0=z0,zg=50,zl=zl,n_m=2.25,n_g=2.25,alpha_np_m=a)
    return m[0],20*np.log10((lmm/6)*np.abs(v)/.5)


def run():
    d=pd.read_csv(ROOT/'data/digitized/liu2022/curves.csv',dtype={'series':str})
    f,a,z0=inputs(d,np.arange(.025,200+1e-9,.025),False)
    front=pd.read_csv(OUT/'front.csv',dtype={'B_GHz':str})
    rows=[]
    for _,r in front.iterrows():
        b=f[-1] if r.B_GHz=='B_full' else float(r.B_GHz);k=f<=b
        for name,zl,lmm in (('optimum',r.ZL_ohm,r.L_mm),('paper',*PAPER)):
            m,g=curves(f,a,z0,zl,lmm)
            rows.append({'B_GHz':r.B_GHz,'design':name,'ZL_ohm':float(zl),'L_mm':float(lmm),'max_M_dB':float(m[k].max()),
                         'f_at_max_GHz':float(f[k][np.argmax(m[k])]),'min_M_dB':float(m[k].min()),
                         'ripple_pp_dB':float(max(0.0,m[k].max())-min(0.0,m[k].min())),'min_G_dB':float(g[k].min())})
    (OUT/'observation_peaking.json').write_text(json.dumps({'note':'post-hoc observation, not a B05 criterion','rows':rows},indent=2)+'\n')
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
    fig,(a1,a2)=plt.subplots(1,2,figsize=(12,4.4),layout='constrained',sharex=True)
    picks=[('100',SERIES[0]),('150',SERIES[1]),('B_full',SERIES[2])]
    for key,col in picks:
        r=front[front.B_GHz==key].iloc[0];m,g=curves(f,a,z0,r.ZL_ohm,r.L_mm)
        lab=f'J* optimum for B={key}: {r.ZL_ohm:g} ohm, {r.L_mm:g} mm'
        a1.plot(f,np.where(np.diff(f,prepend=f[0])>.03,np.nan,m),color=col,lw=1.2,label=lab);a2.plot(f,np.where(np.diff(f,prepend=f[0])>.03,np.nan,g),color=col,lw=1.2,label=lab)
    m,g=curves(f,a,z0,*PAPER)
    a1.plot(f,np.where(np.diff(f,prepend=f[0])>.03,np.nan,m),color=INK,lw=1.2,ls='--',label='paper design: 40 ohm, 6 mm')
    a2.plot(f,np.where(np.diff(f,prepend=f[0])>.03,np.nan,g),color=INK,lw=1.2,ls='--',label='paper design: 40 ohm, 6 mm')
    a1.axhline(-3,color=INK2,lw=.8,ls=':');a1.axhline(0,color=INK2,lw=.6)
    a1.set_ylabel('normalized EO response M re DC divider (dB)');a2.set_ylabel('absolute response G = 20log10[(L/6 mm)|Vavg|/0.5 V] (dB)')
    for ax in (a1,a2):ax.set_xlabel('frequency (GHz); gaps = unsupported Z0 readout');ax.grid(alpha=.25);ax.legend(fontsize=7,frameon=False)
    fig.suptitle('B05 post-hoc observation | J-optimal designs equalize by low ZL and peak above 0 dB; model sensitivity, not measurement',fontsize=10,color=INK)
    fig.savefig(OUT/'observation_response.png',dpi=140);plt.close(fig)
    return rows


if __name__=='__main__':
    for r in run():print({k:(round(v,3) if isinstance(v,float) else v) for k,v in r.items()})
