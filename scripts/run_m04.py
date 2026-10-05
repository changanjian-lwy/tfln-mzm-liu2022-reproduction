"""M04: check the 38-ohm flatness and S11 statements about Fig.7 against the author's own measured
curves. See M04 BOUNDARY.

Usage: python scripts/run_m04.py
Input: M03 m03_curves.csv (digitized_paper_measurement, 'single' and 'merged' points). No model is run
and the author's dashed calculations are not used.
"""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
M03=ROOT/'data/digitized/liu2022/m03_curves.csv'
M03_RESULTS=ROOT/'experiments/track_M_measurement/M03_fig7_digitization/results.json'
M03_META=ROOT/'data/digitized/liu2022/m03_metadata.json'
MANIFEST=ROOT/'experiments/EVIDENCE_MANIFEST.json'
OUT=ROOT/'experiments/track_M_measurement/M04_fig7_statement_checks'
# Locked before running (M04 BOUNDARY).
RISE_MAX_GHZ=10.0;FLAT_MAX_GHZ=50.0;FLAT_DB=1.0
S11_LEVEL=-10.0;S11_APPROX=0.10;S11_SUB_BAND=(0,50)
FLAGS=('single','merged')


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def gate_inputs():
    m=json.loads(MANIFEST.read_text());entry=next((e for e in m['experiments'] if e['id']=='M03_fig7_digitization'),None)
    c=json.loads(M03_RESULTS.read_text())['criteria']
    return {'m03_registered':entry is not None,
            'm03_curves_hash_matches':entry is not None and entry['outputs'].get(str(M03.relative_to(ROOT)))==sha(M03),
            'm03_criterion_2':bool(c['2_reference_lines_pm0.09dB']),'m03_criterion_3':bool(c['3_tick_residual_le_1px']),
            'eo_38_usable':bool(c['4_retained_ge_70pct']['EO_response_38ohm'] and c['5_strict_median_le_1px']['EO_response_38ohm']),
            's11_38_usable':bool(c['4_retained_ge_70pct']['S11_38ohm'] and c['5_strict_median_le_1px']['S11_38ohm'])}


def flatness(s):
    s=s[s.x_value<=FLAT_MAX_GHZ]
    i,j=s.value.idxmax(),s.value.idxmin()
    vmax,umax,vmin,umin=s.value[i],s.uncertainty_value[i],s.value[j],s.uncertainty_value[j]
    hi=(vmax+umax)-(vmin-umin);lo=(vmax-umax)-(vmin+umin)
    return {'points':len(s),'last_point_GHz':float(s.x_value.max()),'max_dB':float(vmax),'max_at_GHz':float(s.x_value[i]),
            'min_dB':float(vmin),'min_at_GHz':float(s.x_value[j]),'peak_to_peak_dB':float(vmax-vmin),
            'P_lo_dB':float(lo),'P_hi_dB':float(hi),'verdict':'holds' if hi<FLAT_DB else 'does_not_hold' if lo>=FLAT_DB else 'undecided'}


def rise(s):
    s=s[s.x_value<=RISE_MAX_GHZ];above=s[(s.value-s.uncertainty_value)>0]
    return {'points':len(s),'points_definitely_above_0dB':len(above),'max_dB':float(s.value.max()),
            'max_at_GHz':float(s.x_value[s.value.idxmax()]),'verdict':'initial_rise' if len(above) else 'not_seen'}


def s11_below(s,clipped,band=None):
    if band:s=s[(s.x_value>band[0])&(s.x_value<=band[1])];clipped=[c for c in clipped if band[0]<c['x_value']<=band[1]]
    above=s[(s.value-s.uncertainty_value)>S11_LEVEL];n=len(s)+len(clipped);frac=len(above)/n
    i=s.value.idxmax()
    return {'points':len(s),'clipped_counted_below':len(clipped),'definitely_above':len(above),'fraction_above':frac,
            'max_dB':float(s.value[i]),'max_at_GHz':float(s.x_value[i]),
            'verdict':'strictly_below' if len(above)==0 else 'approximately_below' if frac<=S11_APPROX else 'not_below'}


def run():
    gate=gate_inputs();d=pd.read_csv(M03)
    clipped={k:[c for c in v if c['edge']=='bottom'] for k,v in json.loads(M03_META.read_text())['clipped'].items()}  # below -35 dB counts as below
    d=d[d.flag.isin(FLAGS)]
    nonfinite=int((~np.isfinite(d[['x_value','value','uncertainty_value']].to_numpy(dtype=float))).sum())
    res={'grade_cap':'DIAGNOSTIC_ONLY','inputs':{str(M03.relative_to(ROOT)):sha(M03)},'gate_inputs':gate,'nonfinite_values':nonfinite}
    if not (gate['m03_registered'] and gate['m03_curves_hash_matches']):
        res['stopped']='M03 not registered or hash differs';(OUT/'results.json').write_text(json.dumps(res,indent=2)+'\n');return res
    eo={l:d[d.series==f'EO_response_{l}'].sort_values('x_value') for l in ('38ohm','57ohm','83ohm','open')}
    s11={l:d[d.series==f'S11_{l}'].sort_values('x_value') for l in ('38ohm','57ohm','83ohm')}
    if gate['m03_criterion_2'] and gate['m03_criterion_3'] and gate['eo_38_usable']:
        res['U1a_initial_rise']=rise(eo['38ohm']);res['U1b_flatness_38ohm']=flatness(eo['38ohm'])
        res['U1b_reported_other_loads']={l:flatness(s) for l,s in eo.items() if l!='38ohm' and len(s)}
    else:res['U1']='not run (M03 Fig.7(a) calibration or 38-ohm EO series unusable)'
    if gate['m03_criterion_3'] and gate['s11_38_usable']:
        res['U2_s11_38ohm']=s11_below(s11['38ohm'],clipped['S11_38ohm'])
        res['U2_reported']={'38ohm_0_50GHz':s11_below(s11['38ohm'],clipped['S11_38ohm'],S11_SUB_BAND),
                            **{l:s11_below(s,clipped[f'S11_{l}']) for l,s in s11.items() if l!='38ohm'}}
    else:res['U2']='not run (M03 Fig.7(b) calibration or 38-ohm S11 series unusable)'
    (OUT/'results.json').write_text(json.dumps(res,indent=2)+'\n')
    fig,axs=plt.subplots(1,2,figsize=(13,4.8),layout='constrained')
    col={'38ohm':'purple','57ohm':'goldenrod','83ohm':'red','open':'blue'}
    if 'U1b_flatness_38ohm' in res:
        s=eo['38ohm'];f=res['U1b_flatness_38ohm']
        axs[0].errorbar(s.x_value,s.value,yerr=s.uncertainty_value,fmt='.',ms=2,lw=.4,color='purple',label='38 ohm measured (M03)')
        for v in (f['max_dB'],f['min_dB']):axs[0].axhline(v,color='k',lw=.6,ls=':')
        axs[0].set_title(f"U1b | 38 ohm EO: peak-to-peak {f['peak_to_peak_dB']:.2f} dB (bounds {f['P_lo_dB']:.2f}-{f['P_hi_dB']:.2f}) -> {f['verdict']}",fontsize=9)
        axs[0].set_xlabel('Frequency (GHz)');axs[0].set_ylabel('dB');axs[0].grid(alpha=.2);axs[0].legend(fontsize=7)
    if 'U2_s11_38ohm' in res:
        for l,s in s11.items():
            axs[1].errorbar(s.x_value,s.value,yerr=s.uncertainty_value,fmt='.',ms=1.5,lw=.3,color=col[l],alpha=1 if l=='38ohm' else .35,label=f'S11 {l}')
        axs[1].axhline(S11_LEVEL,color='k',lw=.8,ls='--')
        u=res['U2_s11_38ohm'];axs[1].set_title(f"U2 | 38 ohm S11: {u['fraction_above']:.1%} of points definitely above -10 dB -> {u['verdict']}",fontsize=9)
        axs[1].set_xlabel('Frequency (GHz)');axs[1].set_ylabel('dB');axs[1].grid(alpha=.2);axs[1].legend(fontsize=7)
    fig.suptitle("M04 | Fig.7 measured curves vs the paper's statements (measurement class; not Figure 3)",fontsize=11)
    fig.savefig(OUT/'m04_checks.png',dpi=150);plt.close(fig)
    return res


if __name__=='__main__':print(json.dumps(run(),indent=1))
