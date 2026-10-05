"""M02 post-hoc observation (written after the run; not a criterion).

1. Where does the author's '160 GHz' arrow in Fig.6(a) point? The shaft column is located from
   dark pixels (max(RGB)<140, the M01 occluder test) in rows 178-190, columns 880-910, and mapped
   with the M01 calibration; uncertainty +-2 px.
2. How do the T3 classes split by band, and above which frequency is every point 'higher'?
It reads the M01 raster and M02 outputs; it changes no M01 or M02 output.

Usage: python scripts/m02_observation.py --images ../tmp/pdfs/a03
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
import pandas as pd
from PIL import Image
sys.path.insert(0,str(Path(__file__).resolve().parent))
from digitize_m01 import PANELS,IMAGE,DARK_MAX,to_phys

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'experiments/track_M_measurement/M02_fig6_statement_checks'
ARROW_ROWS=(178,191);ARROW_COLS=(880,911)
BANDS=[(20,50),(50,100),(100,170)]


def arrow(images):
    rgb=np.array(Image.open(images/IMAGE).convert('RGB')).astype(int)
    dark=rgb.max(axis=2)<DARK_MAX
    ys,xs=np.nonzero(dark[ARROW_ROWS[0]:ARROW_ROWS[1],ARROW_COLS[0]:ARROW_COLS[1]])
    x=float(np.mean(xs)+ARROW_COLS[0]);P=PANELS['fig6a']
    dx=abs((P['x_phys'][1]-P['x_phys'][0])/(P['x_px'][1]-P['x_px'][0]))
    return {'shaft_column_px':x,'pixels':int(len(xs)),'frequency_GHz':to_phys(x,P['x_px'],P['x_phys']),'uncertainty_GHz':2*dx}


def t3_bands():
    df=pd.read_csv(OUT/'t3_loss_comparison.csv')
    bands=[]
    for lo,hi in BANDS:
        s=df[(df.frequency_GHz>lo)&(df.frequency_GHz<=hi)]
        bands.append({'band_GHz':[lo,hi],'points':len(s),**{c:int((s['class']==c).sum()) for c in ('higher','undecided','lower')},
                      'fraction_higher':float((s['class']=='higher').mean())})
    not_higher=df.frequency_GHz[df['class']!='higher']
    above=df.frequency_GHz[df.frequency_GHz>not_higher.max()]
    return {'bands':bands,'last_not_higher_GHz':float(not_higher.max()),
            'all_higher_from_GHz':float(above.min()) if len(above) else None}


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--images',type=Path,required=True)
    res={'note':'post-hoc observation, written after M02 ran; not a criterion','arrow_160GHz_label':arrow(ap.parse_args().images),'t3_by_band':t3_bands()}
    (OUT/'observation.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=1))
