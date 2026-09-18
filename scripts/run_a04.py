"""Frozen input comparison; never fits parameters to the output targets."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from tfln_mzm.digitized import interpolate_supported
from tfln_mzm.interaction import average_voltage
from tfln_mzm.microwave import propagation,db_per_mm_to_np_per_m
from tfln_mzm.termination import s11
from tfln_mzm.response import magnitude_db,eo_response_db

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'experiments/track_A_reproduction/A04_digitized_comparison'
SOURCE=ROOT/'data/digitized/liu2022/curves.csv'

def run():
    data=pd.read_csv(SOURCE,dtype={'series':str})
    inputs=[data[data.panel==p].sort_values('frequency_GHz') for p in ['fig2a','fig2c']]
    summaries=[];tables=[]
    fig,axs=plt.subplots(3,1,figsize=(11,11),layout='constrained')
    colors={'20':'#d98b00','40':'red','50':'black','80':'blue','open':'green'}
    for ax,panel in zip(axs,['fig3a_inset','fig3b','fig3c']):
        for label,target in data[data.panel==panel].groupby('series',sort=False):
            target=target.sort_values('frequency_GHz').copy()
            f=target.frequency_GHz.to_numpy();y=target.value.to_numpy()
            (loss,ml),(z0,mz)=[interpolate_supported(f,d.frequency_GHz,d.value,1) for d in inputs]
            supported=ml&mz
            predicted=np.full(len(f),np.nan)
            load=np.inf if label=='open' else float(label)
            alpha=db_per_mm_to_np_per_m(loss[supported]);freq=f[supported]*1e9
            v=average_voltage(freq,length_m=.006,z0=z0[supported],zg=50,zl=load,n_m=2.25,n_g=2.25,alpha_np_m=alpha)
            ref=1 if label=='open' else load/(50+load)
            if panel=='fig3a_inset':p=abs(v)
            elif panel=='fig3b':p=eo_response_db(v,ref)
            else:p=magnitude_db(s11(propagation(freq,2.25,alpha),.006,z0[supported],load))
            predicted[supported]=p
            slope=np.gradient(y,f)
            neighbor_ok=np.ones(len(f),bool)
            gaps=np.diff(f)>1
            neighbor_ok[:-1]&=~gaps;neighbor_ok[1:]&=~gaps
            allowance=target.uncertainty_value.to_numpy()+abs(slope)*target.uncertainty_frequency_GHz.to_numpy()
            scored=supported&neighbor_ok&np.isfinite(predicted)
            e=predicted[scored]-y[scored]
            within=float(np.mean(abs(e)<=allowance[scored])) if len(e) else None
            unresolved=int(np.sum(supported&~np.isfinite(predicted)))
            passed= len(e)>=100 and within>=.9 and unresolved==0
            summary={'panel':panel,'load':label,'target_count':len(f),'scored_count':int(scored.sum()),'unsupported_input_count':int((~supported).sum()),'nonfinite_prediction_count':unresolved,'target_gap_excluded_count':int((supported&~neighbor_ok).sum()),'MAE':float(np.mean(abs(e))) if len(e) else None,'RMSE':float(np.sqrt(np.mean(e**2))) if len(e) else None,'max_abs_error':float(max(abs(e))) if len(e) else None,'within_readout_fraction':within,'diagnostic_pass':bool(passed),'unit':target.unit.iloc[0]}
            summaries.append(summary)
            target['prediction']=predicted;target['supported']=supported;target['scored']=scored;target['allowance']=allowance
            tables.append(target)
            ax.scatter(f,y,s=4,color=colors[label],alpha=.35)
            ax.plot(f,predicted,color=colors[label],lw=1,label=label)
        ax.set_title(panel+' | dots: digitized author calculation; lines: model');ax.set_xlabel('Frequency (GHz)');ax.set_ylabel('V' if panel=='fig3a_inset' else 'dB');ax.legend(ncol=5,fontsize=8);ax.grid(alpha=.2)
    fig.suptitle('A04 | Digitized loss + impedance; assumed constant microwave index\nMissing source intervals omitted; no parameter fit',fontsize=13)
    fig.savefig(OUT/'comparison.png',dpi=150);plt.close(fig)
    pd.concat(tables).to_csv(OUT/'comparison.csv',index=False)
    result={'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'claim':'Partial independent comparison under explicit assumptions; not full reproduction','criterion':'At least 100 scored points; >=90% within readout allowance; no unresolved nonfinite predictions','results':summaries}
    (OUT/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for r in summaries:print(r['panel'],r['load'],'MAE',round(r['MAE'],4),'within',round(r['within_readout_fraction'],3),'pass',r['diagnostic_pass'])

if __name__=='__main__':run()
