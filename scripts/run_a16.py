"""A16: statement S2 (similar high-frequency roll-off) with a supported 5-20 GHz reference. See A16 BOUNDARY.

Roll-off windows and thresholds are A09's; only the low-frequency reference changes.
"""
from pathlib import Path
import json,sys
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parent))
from run_a09 import LOADS,S2_ABS,S2_REL,MIN_PTS,series,model_at,window

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'experiments/track_A_reproduction/A16_s2_supported_reference'
LOW=(5.0,20.0)


def s2(eo):
    roll={l:window(*eo[l],180,195,np.mean)-window(*eo[l],95,110,np.mean) for l in LOADS}
    low={l:window(*eo[l],*LOW,np.mean) for l in LOADS}
    spread=max(roll.values())-min(roll.values());low_spread=max(low.values())-min(low.values())
    return {'rolloff_dB':roll,'rolloff_spread_dB':spread,'low_mean_5_20GHz_dB':low,'low_spread_dB':low_spread,
            'S2':bool(spread<=S2_ABS and spread<=S2_REL*low_spread)}


def run():
    d=pd.read_csv(ROOT/'data/digitized/liu2022/curves.csv',dtype={'series':str})
    a09=json.loads((ROOT/'experiments/track_A_reproduction/A09_qualitative_statements/results.json').read_text())
    grid=np.arange(.05,200+1e-9,.05)
    eo_m={l:model_at(d,grid,l)[:2] for l in LOADS}
    eo_r={l:(series(d,'fig3b',l).frequency_GHz.to_numpy(),series(d,'fig3b',l).value.to_numpy()) for l in LOADS}
    r,m=s2(eo_r),s2(eo_m)
    g1=all(abs(r['rolloff_dB'][l]-a09['readout']['S2_rolloff_dB'][l])<1e-12 and abs(m['rolloff_dB'][l]-a09['model']['S2_rolloff_dB'][l])<1e-12 for l in LOADS)
    verdict='not_scored_readout_inconsistent' if not r['S2'] else 'reproduced' if m['S2'] else 'not_reproduced'
    result={'grade':'QUALITATIVE_ONLY','scope':'statement S2 relative part; low-frequency reference = mean M over (5,20] GHz','min_points':MIN_PTS,
            'criteria':{'1_rolloff_reproduces_A09':g1,'2_readout_consistent':r['S2'],'3_model_verdict':verdict},'readout':r,'model':m}
    (OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    print(json.dumps(run(),indent=1))
