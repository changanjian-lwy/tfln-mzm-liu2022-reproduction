"""B01: termination-resistance sweep on the frozen A04 L1 provider. See B01 BOUNDARY.

Only ZL changes. Primary results use supported frequencies only; a Z0-gap-bridged
variant is reported for comparison and never enters the criteria.
"""
from pathlib import Path
import json,sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parent))
from run_a09 import a04_fig3b_mae,series
from tfln_mzm.digitized import interpolate_supported
from tfln_mzm.interaction import average_voltage
from tfln_mzm.microwave import propagation,db_per_mm_to_np_per_m
from tfln_mzm.response import eo_response_db,magnitude_db,bandwidth_3db
from tfln_mzm.termination import s11

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'experiments/track_B_extensions/B01_termination_sweep'
LOADS=[float(z) for z in range(20,81)]+[np.inf];F0=1e-3;FLAT_GHZ=50.0


def inputs(d,f,bridged):
    ls=series(d,'fig2a');zs=series(d,'fig2c')
    loss,ml=interpolate_supported(f,ls.frequency_GHz,ls.value,1)
    if bridged:
        zf=zs.frequency_GHz.to_numpy();z0=np.interp(f,zf,zs.value);mz=(f>=zf[0])&(f<=zf[-1])
    else:
        z0,mz=interpolate_supported(f,zs.frequency_GHz,zs.value,1)
    ok=ml&mz
    return f[ok],db_per_mm_to_np_per_m(loss[ok]),z0[ok]


def metrics(f,a,z0,zl):
    ref=1.0 if np.isinf(zl) else zl/(50+zl)
    v=average_voltage(f*1e9,length_m=.006,z0=z0,zg=50,zl=zl,n_m=2.25,n_g=2.25,alpha_np_m=a)
    m=eo_response_db(v,ref)
    low=f<=FLAT_GHZ;p50=float(max(0.0,m[low].max())-min(0.0,m[low].min()))
    b=bandwidth_3db(np.concatenate([[F0],f]),np.concatenate([[0.0],m]))
    e=b['events'][0] if b['events'] else None
    out={'ZL_ohm':'open' if np.isinf(zl) else zl,'Vref_V':ref,'P50_dB':p50,'D1_GHz':b['bandwidth_hz'],'D1_censored':b['censored'],
         'D1_bracket_GHz':e['bracket_hz'] if e else None,'M_at_50GHz_nearest_dB':float(m[np.argmin(abs(f-50))]),'nonfinite_M':int((~np.isfinite(m)).sum())}
    if not np.isinf(zl):
        sd=magnitude_db(s11(propagation(f*1e9,2.25,a),.006,z0,zl))
        out.update(maxS11_0_50GHz_dB=float(sd[low].max()),S11_exact_zeros=int(np.isneginf(sd).sum()),nonfinite_S11=int(np.isnan(sd).sum()+np.isposinf(sd).sum()))
    return out,m


def sweep(d,step,bridged):
    f,a,z0=inputs(d,np.arange(step,200+1e-9,step),bridged)
    rows=[];curves={}
    for zl in LOADS:
        r,m=metrics(f,a,z0,zl);rows.append(r);curves[r['ZL_ohm']]=(f,m)
    return pd.DataFrame(rows),curves


def run():
    d=pd.read_csv(ROOT/'data/digitized/liu2022/curves.csv',dtype={'series':str})
    mae=a04_fig3b_mae(d)
    archived={r['load']:r['MAE'] for r in json.loads((ROOT/'experiments/track_A_reproduction/A04_digitized_comparison/results.json').read_text())['results'] if r['panel']=='fig3b'}
    fine,curves=sweep(d,.025,False);coarse,_=sweep(d,.05,False);bridged,_=sweep(d,.025,True)
    finite=fine[fine.ZL_ohm!='open']
    dp=abs(fine.P50_dB-coarse.P50_dB).max()
    both=fine.D1_GHz.notna()&coarse.D1_GHz.notna()
    narrow=np.array([b is not None and (b[1]-b[0])<=1 for b in fine.D1_bracket_GHz])
    dd=float(abs(fine.D1_GHz[both&narrow]-coarse.D1_GHz[both&narrow]).max()) if (both&narrow).any() else 0.0
    criteria={'1_A04_mae_reproduced':all(abs(mae[l]-archived[l])<1e-9 for l in archived),
              '2_P50_grid_shift_max_dB':float(dp),'2_D1_grid_shift_max_GHz':dd,'2_pass':bool(dp<=.02 and dd<=.1 and (fine.D1_censored==coarse.D1_censored).all()),
              '3_finite':bool((fine.nonfinite_M==0).all() and (finite.nonfinite_S11==0).all()),'3_S11_exact_zero_points':int(finite.S11_exact_zeros.sum())}
    def band(df):
        ok=df[(df.ZL_ohm!='open')&(df.P50_dB<1)]
        return None if ok.empty else {'ZL_min':float(ok.ZL_ohm.min()),'ZL_max':float(ok.ZL_ohm.max()),'contiguous':bool(len(ok)==ok.ZL_ohm.max()-ok.ZL_ohm.min()+1)}
    best=finite.loc[finite.P50_dB.astype(float).idxmin()]
    inband=finite[finite.P50_dB<1]
    summary={'P50_lt_1dB_ZL_range':band(fine),'P50_minimum':{'ZL_ohm':float(best.ZL_ohm),'P50_dB':float(best.P50_dB)},
             'in_range':None if inband.empty else {'D1_GHz_min':float(inband.D1_GHz.min()) if inband.D1_GHz.notna().any() else None,
                'D1_censored_any':bool(inband.D1_censored.any()),'maxS11_dB_range':[float(inband.maxS11_0_50GHz_dB.min()),float(inband.maxS11_0_50GHz_dB.max())],
                'Vref_V_range':[float(inband.Vref_V.min()),float(inband.Vref_V.max())]},
             'bridged_P50_lt_1dB_ZL_range':band(bridged),
             'bridged_minus_primary_P50_dB':{'median':float((bridged.P50_dB-fine.P50_dB).median()),'max':float((bridged.P50_dB-fine.P50_dB).max())}}
    # Reported-only observation (added after the first run): same P50 on the author's Fig.3(b) readout.
    readout_p50={}
    for l in ('20','40','50','80','open'):
        s=series(d,'fig3b',l);m=s.value.to_numpy()[s.frequency_GHz.to_numpy()<=FLAT_GHZ]
        readout_p50[l]=float(max(0.0,m.max())-min(0.0,m.min()))
    summary['observation_readout_P50_dB']=readout_p50
    fine.to_csv(OUT/'sweep.csv',index=False)
    result={'grade':'SENSITIVITY_ONLY','track':'B','parent':'A04 L1 provider','changed':'ZL 20-80 ohm step 1, plus open',
            'criteria':criteria,'summary':summary,'A04_fig3b_mae_recomputed':mae}
    (OUT/'results.json').write_text(json.dumps(result,indent=2,default=float)+'\n')
    z=finite.ZL_ohm.astype(float)
    fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    axs[0,0].plot(z,finite.P50_dB,label='supported points only');axs[0,0].plot(z,bridged[bridged.ZL_ohm!='open'].P50_dB,'--',label='Z0 gaps bridged (assumption)')
    axs[0,0].axhline(1,color='r',ls=':');axs[0,0].set_ylabel('0-50 GHz peak-to-peak (dB)');axs[0,0].legend(fontsize=8)
    axs[0,1].plot(z,finite.D1_GHz);axs[0,1].set_ylabel('first -3 dB crossing (GHz)')
    axs[1,0].plot(z,finite.maxS11_0_50GHz_dB);axs[1,0].axhline(-10,color='r',ls=':');axs[1,0].set_ylabel('max S11, 0-50 GHz (dB)')
    axs[1,1].plot(z,finite.Vref_V);axs[1,1].set_ylabel('low-frequency drive Vref (V)')
    for ax in axs.ravel():ax.set_xlabel('ZL (ohm)');ax.grid(alpha=.25)
    fig.suptitle('B01 | termination sweep on the A04 L1 provider; model sensitivity, not measurement',fontsize=12)
    fig.savefig(OUT/'termination_sweep.png',dpi=140);plt.close(fig)
    return result


if __name__=='__main__':
    r=run();print(json.dumps({k:r[k] for k in ('criteria','summary')},indent=1,default=float))
