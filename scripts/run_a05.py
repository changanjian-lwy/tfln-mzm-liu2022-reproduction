"""Independent cross-panel and phase-invariant diagnostics. No model fitting."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from tfln_mzm.digitized import interpolate_supported
from tfln_mzm.microwave import db_per_mm_to_np_per_m
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'experiments/track_A_reproduction/A05_convention_audit'

def bounds(z,loss,load):
 r=(z-50)/(z+50);rl=(load-z)/(load+z)
 a=abs(rl)*np.exp(-2*db_per_mm_to_np_per_m(loss)*.006)
 p=abs((r+a)/(1+r*a));m=abs((r-a)/(1-r*a))
 return np.minimum(p,m),np.maximum(p,m)

def run():
 d=pd.read_csv(ROOT/'data/digitized/liu2022/curves.csv',dtype={'series':str})
 def curve(panel,label=None):
  t=d[d.panel==panel]
  if label is not None:t=t[t.series==label]
  return t.sort_values('frequency_GHz')
 def interp(t,f,col='value'):return interpolate_supported(f,t.frequency_GHz,t[col],1)
 refs=[];rows=[];zs=curve('fig2c');ls=curve('fig2a')
 for load in ['20','40','50','80','open']:
  v=curve('fig3a_inset',load);m=curve('fig3b',load)
  f=v.frequency_GHz.to_numpy();mv,ok=interp(m,f);mu,_=interp(m,f,'uncertainty_value')
  ok&=(f>=2)&(f<=49)
  ref=v.value.to_numpy()[ok]/10**(mv[ok]/20)
  # Approximate vertical uncertainty; frequency mismatch not included here.
  u=ref*(v.uncertainty_value.to_numpy()[ok]/v.value.to_numpy()[ok]+np.log(10)/20*mu[ok])
  dc=1 if load=='open' else float(load)/(50+float(load))
  refs.append({'load':load,'points':len(ref),'divider_V':dc,'implied_reference_median_V':float(np.median(ref)),'p05_V':float(np.quantile(ref,.05)),'p95_V':float(np.quantile(ref,.95)),'median_vertical_uncertainty_V':float(np.median(u)),'equivalent_normalization_offset_dB':float(20*np.log10(np.median(ref)/dc))})
  if load=='open':continue
  t=curve('fig3c',load);f=t.frequency_GHz.to_numpy();z,okz=interp(zs,f);loss,okl=interp(ls,f)
  for i in np.flatnonzero(okz&okl):
   low,high=bounds(z[i],loss[i],float(load));err=t.uncertainty_value.iloc[i]
   targetlo=10**((t.value.iloc[i]-err)/20);targethi=10**((t.value.iloc[i]+err)/20)
   outside=targethi<low or targetlo>high
   # Local neighborhood envelope includes source ordinate error and x error.
   radius=t.uncertainty_frequency_GHz.iloc[i]+max(zs.uncertainty_frequency_GHz.max(),ls.uncertainty_frequency_GHz.max())
   zz=zs[abs(zs.frequency_GHz-f[i])<=radius];ll=ls[abs(ls.frequency_GHz-f[i])<=radius]
   if len(zz)==0 or len(ll)==0:continue
   zmin=float(min((zz.value-zz.uncertainty_value).min(),z[i]));zmax=float(max((zz.value+zz.uncertainty_value).max(),z[i]))
   lmin=max(0,float(min((ll.value-ll.uncertainty_value).min(),loss[i])));lmax=float(max((ll.value+ll.uncertainty_value).max(),loss[i]))
   # Dense bounded sensitivity envelope, not a mathematically certified interval.
   gridz=np.linspace(zmin,zmax,101)[:,None];gridl=np.linspace(lmin,lmax,21)[None,:]
   bl,bh=bounds(gridz,gridl,float(load));elo=float(bl.min());ehi=float(bh.max())
   rows.append({'load':load,'frequency_GHz':float(f[i]),'target_magnitude':float(10**(t.value.iloc[i]/20)),'phase_bound_low':float(low),'phase_bound_high':float(high),'outside_nominal_bounds':bool(outside),'envelope_low':elo,'envelope_high':ehi,'outside_sampled_uncertainty_envelope':bool(targethi<elo or targetlo>ehi)})
 frame=pd.DataFrame(rows);frame.to_csv(OUT/'phase_bounds.csv',index=False)
 summary=[]
 for load,t in frame.groupby('load'):
  summary.append({'load':load,'points':len(t),'outside_nominal':int(t.outside_nominal_bounds.sum()),'outside_sampled_uncertainty_envelope':int(t.outside_sampled_uncertainty_envelope.sum())})
 result={'scope':'diagnosis_only_no_fit','cross_panel_reference':refs,'phase_bounds':summary,'uncertainty_note':'reference estimates propagate vertical readout only; phase envelope samples source frequency/ordinate errors but is not a certified interval or complete input-model uncertainty'}
 (OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
 comparison=pd.read_csv(ROOT/'experiments/track_A_reproduction/A04_digitized_comparison/comparison.csv',dtype={'series':str})
 diagnostic=json.loads(json.dumps(result))
 for ref in diagnostic['cross_panel_reference']:
  t=comparison[(comparison.panel=='fig3b')&(comparison.series==ref['load'])&comparison.scored]
  residual=t.prediction-t.value
  ref['A04_EO_MAE_dB']=float(abs(residual).mean())
  ref['cross_panel_renormalized_MAE_dB']=float(abs(residual-ref['equivalent_normalization_offset_dB']).mean())
 diagnostic['renormalization_warning']='Output-calibrated diagnostic only. Never replace A04 reference or claim independent validation.'
 (OUT/'normalization_diagnostic.json').write_text(json.dumps(diagnostic,indent=2)+'\n')
 print(json.dumps(result,indent=2))
if __name__=='__main__':run()
