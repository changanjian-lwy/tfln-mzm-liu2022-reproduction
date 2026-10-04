"""A08: does a matched 50-ohm provider reproduce the Fig.2(b) bandwidth curve? See A08 BOUNDARY.

Branch C-06b (assumptions): Z0=50 ohm, loss independent of BCB, nm(BCB) frequency
independent. Loss outside Fig.2(a) coverage uses an input-side fit a*sqrt(f)+b*f;
no output target is fitted or adjusted.
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
from tfln_mzm.microwave import db_per_mm_to_np_per_m
from tfln_mzm.response import eo_response_db,bandwidth_3db

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'experiments/track_A_reproduction/A08_c06_matched_provider'
A06=ROOT/'experiments/track_A_reproduction/A06_fig2b_fig3a_main_digitization/results.json'
L=.006;Z=50.0;F0_GHZ=1e-3;BCB_GAP=.05;MAX_POINT_GHZ=205.0


def loss_provider(data):
    t=data[data.panel=='fig2a'].sort_values('frequency_GHz')
    f=t.frequency_GHz.to_numpy();v=t.value.to_numpy();u=t.uncertainty_value.to_numpy()
    A=np.column_stack([np.sqrt(f),f]);w=1/u
    coef,*_=np.linalg.lstsq(A*w[:,None],v*w,rcond=None)
    fit_rms=float(np.sqrt(np.mean((A@coef-v)**2)))
    def alpha_db(q):
        q=np.asarray(q,float);fitted=coef[0]*np.sqrt(q)+coef[1]*q
        inside,ok=interpolate_supported(np.clip(q,f[0],f[-1]),f,v,1)
        return np.where((q>=f[0])&(q<=f[-1]),inside,fitted)
    def unc_db(q):
        return np.interp(np.asarray(q,float),f,u)  # constant end values outside coverage
    return alpha_db,unc_db,{'a_dB_per_mm_per_sqrtGHz':float(coef[0]),'b_dB_per_mm_per_GHz':float(coef[1]),'fit_rms_dB_per_mm':fit_rms,
                            'coverage_GHz':[float(f[0]),float(f[-1])]}


def response(f_ghz,nm,ng,alpha_db):
    f=np.asarray(f_ghz,float);nm=np.broadcast_to(np.asarray(nm,float),f.shape)
    v=average_voltage(f*1e9,length_m=L,z0=Z,zg=Z,zl=Z,n_m=nm,n_g=ng,alpha_np_m=db_per_mm_to_np_per_m(alpha_db(f)))
    ref=average_voltage(np.full(f.shape,F0_GHZ*1e9),length_m=L,z0=Z,zg=Z,zl=Z,n_m=nm,n_g=ng,alpha_np_m=db_per_mm_to_np_per_m(alpha_db(np.full(f.shape,F0_GHZ))))
    return eo_response_db(v,ref)


def interval(f,uf,nm,ng,alpha_db,unc_db):
    vals=[response(f+sf*uf,nm,ng,lambda q,s=sa:np.maximum(alpha_db(q)+s*unc_db(q),0)) for sa in (-1,0,1) for sf in (-1,0,1)]
    return np.min(vals,axis=0),np.max(vals,axis=0)


def run():
    data=pd.read_csv(ROOT/'data/digitized/liu2022/curves.csv',dtype={'series':str})
    a6=pd.read_csv(ROOT/'data/digitized/liu2022/a06_curves.csv',dtype={'series':str})
    stars=json.loads(A06.read_text())['derived']
    ng=stars['optical_index_dashed_line']['value']
    alpha_db,unc_db,fit=loss_provider(data)
    blue=a6[(a6.panel=='fig2b')&(a6.series=='nm_at_100GHz')].sort_values('x_value')
    red=a6[(a6.panel=='fig2b')&(a6.series=='bw3dB_50ohm_GHz')].sort_values('x_value')
    bcb=blue.x_value.to_numpy();nm=blue.value.to_numpy()
    bw,ok=interpolate_supported(bcb,red.x_value,red.value,BCB_GAP)
    ubw,_=interpolate_supported(bcb,red.x_value,red.uncertainty_value,BCB_GAP)
    bcb,nm,bw,ubw=bcb[ok],nm[ok],bw[ok],ubw[ok]
    m=response(bw,nm,ng,alpha_db);lo,hi=interval(bw,ubw,nm,ng,alpha_db,unc_db)
    contains=(lo<=-3)&(hi>=-3)
    s_nm=stars['star_nm_at_100GHz']['value'];s_bw=stars['star_bw3dB_50ohm_GHz']['value'];s_u=stars['star_bw3dB_50ohm_GHz']['uncertainty_value']
    sm=float(response(np.array([s_bw]),s_nm,ng,alpha_db)[0]);slo,shi=[float(x[0]) for x in interval(np.array([s_bw]),np.array([s_u]),s_nm,ng,alpha_db,unc_db)]
    # Reported only: matched-model first crossing per BCB, long extrapolation.
    grid=np.arange(.1,400+1e-9,.1);info=[]
    for b,n in zip(bcb[::10],nm[::10]):
        r=bandwidth_3db(np.concatenate([[F0_GHZ],grid]),np.concatenate([[0.0],response(grid,n,ng,alpha_db)]))
        info.append({'BCB_um':float(b),'nm':float(n),'matched_bw_GHz':r['bandwidth_hz'],'censored':r['censored']})
    star_bw=bandwidth_3db(np.concatenate([[F0_GHZ],grid]),np.concatenate([[0.0],response(grid,s_nm,ng,alpha_db)]))
    criteria={'1_star_values_from_A06':True,'2_points_le_205GHz':bool(np.all(bw<=MAX_POINT_GHZ)) and s_bw<=MAX_POINT_GHZ,
        '3_all_finite':bool(np.all(np.isfinite(m)) and np.isfinite(sm)),
        '4_star_interval_contains_-3dB':bool(slo<=-3<=shi),'4_fraction_points_contain_-3dB':float(np.mean(contains)),
        '4_matched_provider_explains':bool(slo<=-3<=shi and np.mean(contains)>=.8)}
    result={'grade':'DIAGNOSTIC_ONLY','branch':'C-06b matched provider: Z0=Zg=ZL=50 ohm; loss independent of BCB; nm frequency independent',
        'ng':ng,'loss_fit_outside_coverage':fit,'compared_points':int(len(bcb)),
        'star':{'BCB_um':stars['star_nm_at_100GHz']['BCB_um'],'nm':s_nm,'fig2b_bw_GHz':s_bw,'matched_response_dB_at_fig2b_bw':sm,'interval_dB':[slo,shi],
                'matched_first_crossing_GHz_long_extrapolation':star_bw['bandwidth_hz']},
        'points_response_dB':{'min':float(m.min()),'median':float(np.median(m)),'max':float(m.max()),'median_interval_width_dB':float(np.median(hi-lo))},
        'criteria':criteria,'reported_only_matched_bw_by_BCB':info}
    (OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    fig,axs=plt.subplots(2,1,figsize=(9,8),layout='constrained',sharex=True)
    axs[0].plot(red.x_value,red.value,'.',ms=2,color='tab:red',label='Fig.2(b) readout (A06)')
    axs[0].plot([r['BCB_um'] for r in info],[r['matched_bw_GHz'] for r in info],'-',color='k',label='matched provider, first crossing (long loss extrapolation; reported only)')
    axs[0].plot([stars['star_bw3dB_50ohm_GHz']['BCB_um']],[s_bw],'*',ms=12,color='tab:red')
    axs[0].set_ylabel('3-dB EO bandwidth (GHz)');axs[0].legend(fontsize=8);axs[0].grid(alpha=.25)
    axs[1].fill_between(bcb,lo,hi,color='0.8',label='interval from loss and readout uncertainty')
    axs[1].plot(bcb,m,color='k',lw=1,label='matched model at the Fig.2(b) bandwidth')
    axs[1].errorbar([1.5],[sm],yerr=[[sm-slo],[shi-sm]],fmt='*',ms=10,color='tab:red',label='design point (stars)')
    axs[1].axhline(-3,color='r',ls=':',lw=1);axs[1].set_ylabel('response at Fig.2(b) bandwidth (dB)');axs[1].set_xlabel('BCB thickness (um)')
    axs[1].legend(fontsize=8);axs[1].grid(alpha=.25)
    fig.suptitle('A08 | C-06b: matched 50-ohm provider vs Fig.2(b); diagnostic only',fontsize=12)
    fig.savefig(OUT/'matched_provider.png',dpi=150);plt.close(fig)
    return result


if __name__=='__main__':
    r=run();print(json.dumps({k:v for k,v in r.items() if k!='reported_only_matched_bw_by_BCB'},indent=1))
