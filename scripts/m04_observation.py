"""M04 post-hoc observation (written after the run; not a criterion).

Where are the two black lines of the author's '1 dB' annotation in Fig.7(a), and how do they compare
with the 38-ohm readout maximum and minimum used in U1b? Line rows are dark pixels (max(RGB)<100)
with more than 20 px in columns 270-390, rows 85-134, mapped with the M03 calibration; +-2 px per line.
It reads the M03 raster and the M04 results; it changes no M03 or M04 output.

Usage: python scripts/m04_observation.py --images ../tmp/pdfs/a03
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
from PIL import Image
sys.path.insert(0,str(Path(__file__).resolve().parent))
from digitize_m03 import PANELS,IMAGE
from digitize_m01 import to_phys,groups_of

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'experiments/track_M_measurement/M04_fig7_statement_checks'
ROWS=(85,135);COLS=(270,390);MIN_PX=20;DARK_MAX=100

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--images',type=Path,required=True);images=ap.parse_args().images
    a=np.array(Image.open(images/IMAGE).convert('RGB')).astype(int);dark=a.max(axis=2)<DARK_MAX
    rows=np.array([y for y in range(*ROWS) if dark[y,COLS[0]:COLS[1]].sum()>MIN_PX])
    P=PANELS['fig7a'];px_per_dB=abs(P['y_px'][1]-P['y_px'][0])/18
    lines=[float(np.mean(g)) for g in groups_of(rows)]
    vals=[to_phys(y,P['y_px'],P['y_phys']) for y in lines]
    u=json.loads((OUT/'results.json').read_text())['U1b_flatness_38ohm']
    res={'note':'post-hoc observation, written after M04 ran; not a criterion',
         'annotation_line_rows':lines,'annotation_line_dB':vals,'annotation_span_dB':abs(vals[0]-vals[-1]) if len(vals)==2 else None,
         'span_uncertainty_dB':4/px_per_dB,'readout_38ohm_max_dB':u['max_dB'],'readout_38ohm_min_dB':u['min_dB'],
         'readout_38ohm_peak_to_peak_dB':u['peak_to_peak_dB']}
    (OUT/'observation.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=1))
