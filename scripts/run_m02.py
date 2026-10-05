"""M02: check the paper's three Fig.6 measurement statements against its own figure. See M02 BOUNDARY.

Usage: python scripts/run_m02.py
Inputs: M01 m01_curves.csv (digitized_paper_measurement) and the A03 Fig.2(a) 'this work'
loss (digitized_paper_calculation). No transmission-line model is run; the T2 sinc penalty
is the L0 analytic expression of BOUNDARIES B04 and is reported, not scored.
"""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
M01=ROOT/'data/digitized/liu2022/m01_curves.csv'
M01_RESULTS=ROOT/'experiments/track_M_measurement/M01_fig6_digitization/results.json'
A03=ROOT/'data/digitized/liu2022/curves.csv'
MANIFEST=ROOT/'experiments/EVIDENCE_MANIFEST.json'
OUT=ROOT/'experiments/track_M_measurement/M02_fig6_statement_checks'
# Locked before running (M02 BOUNDARY).
L_MM=6.0;NG=2.25;RED_TOL=0.004;C=299792458.0
LEVEL_DB=-6.4;CLAIM_GHZ=160.0
BANDS=[(0,10),(10,50),(50,170)];BAND_PASS=0.9
T3_FROM_GHZ=20.0;T3_HOLDS=0.9;T3_FAILS=0.1
MARKS=[50,100,150,160];MAX_GAP=1.0
MEASURED_FLAGS=('single','merged')


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def gate_inputs():
    """G1 input part: M01 hash as registered; M01 criteria 2, 3, 4 passed for the series used here."""
    m=json.loads(MANIFEST.read_text())
    entry=next((e for e in m['experiments'] if e['id']=='M01_fig6_digitization'),None)
    r=json.loads(M01_RESULTS.read_text())['criteria']
    return {'m01_registered':entry is not None,
            'm01_curves_hash_matches':entry is not None and entry['outputs'].get(str(M01.relative_to(ROOT)))==sha(M01),
            'm01_criterion_2':bool(r['2_dashed_line_-6.40pm0.13dB']),'m01_criterion_3':bool(r['3_optical_line_2.250pm0.004']),
            'm01_criterion_4':{k:bool(r['4_retained_ge_70pct'][k]) for k in ('S21','S11','nm_microwave')}}


def interp_supported(x,xs,ys,us):
    """A04 rule: linear between neighbours at most 1 GHz apart; no extrapolation. Uncertainty: larger neighbour."""
    i=int(np.searchsorted(xs,x))
    if i<len(xs) and xs[i]==x:return float(ys[i]),float(us[i])
    if i==0 or i==len(xs) or xs[i]-xs[i-1]>MAX_GAP:return None
    w=(x-xs[i-1])/(xs[i]-xs[i-1])
    return float(ys[i-1]+w*(ys[i]-ys[i-1])),float(max(us[i-1],us[i]))


def t1(s21):
    f,v,u,uf=(s21[k].to_numpy() for k in ('x_value','value','uncertainty_value','uncertainty_x'))
    below=np.flatnonzero(v+u<LEVEL_DB)
    out={'level_dB':LEVEL_DB,'claim_GHz':CLAIM_GHZ,'points':len(f)}
    if not len(below):return {**out,'verdict':'censored_no_crossing_within_readout'}
    b=int(below[0]);above=np.flatnonzero((v[:b]-u[:b])>LEVEL_DB)
    a=int(above[-1]) if len(above) else None
    lo=(f[a]-uf[a]) if a is not None else 0.0;hi=f[b]+uf[b]
    median_first=f[np.flatnonzero(v<LEVEL_DB)[0]]
    return {**out,'f_last_definitely_above_GHz':float(f[a]) if a is not None else None,'f_first_definitely_below_GHz':float(f[b]),
            'interval_GHz':[float(lo),float(hi)],'first_median_below_GHz':float(median_first),
            'verdict':'consistent' if lo<=CLAIM_GHZ<=hi else 'inconsistent'}


def sinc_penalty_db(f_ghz,dn):
    x=np.pi*f_ghz*1e9*dn*L_MM*1e-3/C
    s=np.where(x==0,1.0,np.sin(x)/np.where(x==0,1,x))
    return -20*np.log10(np.abs(s))


def t2(nm):
    f,v,u=(nm[k].to_numpy() for k in ('x_value','value','uncertainty_value'))
    occ=(nm['flag']=='occluded').to_numpy();match=np.abs(v-NG)<=u+RED_TOL
    bands=[]
    for lo,hi in BANDS:
        s=(f>lo)&(f<=hi);n=int(s.sum());k=s&~occ
        bands.append({'band_GHz':[lo,hi],'points':n,'occluded':int((s&occ).sum()),
            'match_fraction':float(match[s].mean()) if n else None,
            'match_fraction_without_occluded':float(match[k].mean()) if k.any() else None,
            'max_abs_deviation':float(np.max(np.abs(v[s]-NG))) if n else None,
            'max_sinc_penalty_dB':float(np.max(sinc_penalty_db(f[s],np.abs(v[s]-NG)))) if n else None,
            'verdict':('consistent' if match[s].mean()>=BAND_PASS else 'inconsistent') if n else 'no_points'})
    return {'ng':NG,'tolerance':'|nm-2.25| <= u_nm + 0.004','L_mm':L_MM,'bands':bands},match


def t3(s21,s11,fig2a):
    s11max=float((s11['value']+s11['uncertainty_value']).max());g=10**(s11max/20)
    M=float(-10*np.log10(1-g**2))
    xs,ys,us=(fig2a[k].to_numpy() for k in ('frequency_GHz','value','uncertainty_value'))
    rows=[];unsupported=0
    for f,v,u in s21[['x_value','value','uncertainty_value']].itertuples(index=False):
        sim=interp_supported(f,xs,ys,us)
        if sim is None:unsupported+=1;continue
        a_sim,u_sim=L_MM*sim[0],L_MM*sim[1];a_meas,u_meas=-v,u
        cls='higher' if a_meas-u_meas-M>a_sim+u_sim else 'lower' if a_meas+u_meas<a_sim-u_sim else 'undecided'
        rows.append({'frequency_GHz':f,'A_meas_dB':a_meas,'u_meas_dB':u_meas,'A_sim_dB':a_sim,'u_sim_dB':u_sim,
                     'mismatch_allowance_dB':M,'class':cls})
    df=pd.DataFrame(rows)
    hi=df[df.frequency_GHz>=T3_FROM_GHZ]
    frac_higher=float((hi['class']=='higher').mean());n_lower=int((df['class']=='lower').sum());frac_lower=n_lower/len(df)
    verdict='holds' if frac_higher>=T3_HOLDS and n_lower==0 else 'does_not_hold' if frac_lower>=T3_FAILS else 'undecided'
    marks=[]
    for m in MARKS:
        j=int(np.argmin(np.abs(df.frequency_GHz-m)))
        if abs(df.frequency_GHz.iloc[j]-m)<=MAX_GAP:
            r=df.iloc[j];marks.append({'target_GHz':m,'frequency_GHz':float(r.frequency_GHz),'A_meas_dB':float(r.A_meas_dB),
                'A_sim_dB':float(r.A_sim_dB),'ratio_upper_bound':float(r.A_meas_dB/r.A_sim_dB),'class':r['class']})
        else:marks.append({'target_GHz':m,'frequency_GHz':None})
    counts={c:int((df['class']==c).sum()) for c in ('higher','undecided','lower')}
    return {'L_mm':L_MM,'s11_max_plus_u_dB':s11max,'mismatch_allowance_dB':M,'comparable_points':len(df),
            'unsupported_points':unsupported,'counts_all':counts,'fraction_higher_from_20GHz':frac_higher,
            'points_from_20GHz':len(hi),'fraction_lower_all':frac_lower,'verdict':verdict,'marks':marks},df


def run():
    gate=gate_inputs()
    m=pd.read_csv(M01);a=pd.read_csv(A03)
    s21=m[(m.series=='S21')&m.flag.isin(MEASURED_FLAGS)].sort_values('x_value')
    s11=m[(m.series=='S11')&m.flag.isin(MEASURED_FLAGS)].sort_values('x_value')
    nm=m[m.series=='nm_microwave'].sort_values('x_value')
    fig2a=a[(a.panel=='fig2a')&(a.series=='this_work')].sort_values('frequency_GHz')
    used=[s21[['x_value','value','uncertainty_value']],s11[['value','uncertainty_value']],nm[['x_value','value','uncertainty_value']],
          fig2a[['frequency_GHz','value','uncertainty_value']]]
    nonfinite=int(sum((~np.isfinite(d.to_numpy(dtype=float))).sum() for d in used))
    res={'grade_cap':'DIAGNOSTIC_ONLY','inputs':{str(M01.relative_to(ROOT)):sha(M01),str(A03.relative_to(ROOT)):sha(A03)},
         'gate_inputs':gate,'nonfinite_values':nonfinite}
    ok_a=gate['m01_criterion_2'] and gate['m01_criterion_4']['S21'];ok_b=gate['m01_criterion_3'] and gate['m01_criterion_4']['nm_microwave']
    if not (gate['m01_registered'] and gate['m01_curves_hash_matches']):
        res['stopped']='M01 not registered or hash differs';(OUT/'results.json').write_text(json.dumps(res,indent=2)+'\n');return res
    res['T1']=t1(s21) if ok_a else 'not run (M01 criterion 2 or S21 coverage failed)'
    if ok_b:res['T2'],match=t2(nm)
    else:res['T2']='not run (M01 criterion 3 or nm coverage failed)'
    if ok_a and gate['m01_criterion_4']['S11']:
        res['T3'],df=t3(s21,s11,fig2a);df.to_csv(OUT/'t3_loss_comparison.csv',index=False)
    else:res['T3']='not run (M01 S21/S11 unusable)'
    (OUT/'results.json').write_text(json.dumps(res,indent=2)+'\n')
    fig,axs=plt.subplots(1,2,figsize=(13,4.8),layout='constrained')
    if isinstance(res['T3'],dict):
        colors={'higher':'tab:red','undecided':'tab:gray','lower':'tab:blue'}
        axs[0].fill_between(df.frequency_GHz,df.A_sim_dB-df.u_sim_dB,df.A_sim_dB+df.u_sim_dB,color='k',alpha=.25,lw=0,
                            label='6 mm x Fig.2(a) simulated loss (A03)')
        for c,col in colors.items():
            s=df[df['class']==c]
            axs[0].errorbar(s.frequency_GHz,s.A_meas_dB,yerr=s.u_meas_dB,fmt='.',ms=2,lw=.4,color=col,label=f'-S21 measured (M01): {c} ({len(s)})')
        axs[0].set_title(f'T3 | measured -S21 vs simulated line loss (mismatch allowance {res["T3"]["mismatch_allowance_dB"]:.2f} dB)')
        axs[0].set_xlabel('Frequency (GHz)');axs[0].set_ylabel('loss (dB)');axs[0].legend(fontsize=7);axs[0].grid(alpha=.2)
    if isinstance(res['T2'],dict):
        dev=nm.value.to_numpy()-NG;tol=nm.uncertainty_value.to_numpy()+RED_TOL;f=nm.x_value.to_numpy()
        axs[1].fill_between(f,-tol,tol,color='tab:red',alpha=.2,lw=0,label='match tolerance u_nm + 0.004')
        for flag,mk in (('single','.'),('merged','x'),('occluded','+')):
            s=(nm.flag==flag).to_numpy()
            axs[1].plot(f[s],dev[s],mk,ms=2.5,color='tab:blue',alpha=1 if flag=='single' else .5,label=f'nm - 2.25 ({flag})')
        for lo,hi in BANDS[1:]:axs[1].axvline(lo,color='k',lw=.6,ls=':')
        axs[1].set_ylim(-0.06,0.4);axs[1].set_title('T2 | extracted microwave index minus ng = 2.25')
        axs[1].set_xlabel('Frequency (GHz)');axs[1].legend(fontsize=7);axs[1].grid(alpha=.2)
    fig.suptitle('M02 | measured Fig.6 readouts vs the paper\'s statements (evidence class: measurement; not Figure 3)',fontsize=11)
    fig.savefig(OUT/'m02_checks.png',dpi=150);plt.close(fig)
    return res


if __name__=='__main__':print(json.dumps(run(),indent=1))
