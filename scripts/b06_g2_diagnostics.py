"""B06 G2 failure diagnosis (written after the failure; prints only, no registered output).

For each (B, R) that failed criterion 2 or 3: feasible-set connected components (8-connectivity) on the
fine grid, the size of the component holding the fine-grid optimum and how many coarse-subgrid points it
contains, and the J* shift. Reports counts and shifts only, not the front values.
"""
import json,sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import ndimage
sys.path.insert(0,str(Path(__file__).resolve().parent))
from run_b05 import BANDS,label
from run_b06 import ROOT,OUT,J,grid,feasible

r=json.loads((OUT/'results.json').read_text())['criteria']
bad2={(x['B_GHz'],x['R_dB']):x['dJ_freq_dB'] for x in r['2_freq_grid']['rows'] if x['dJ_freq_dB'] is None or x['dJ_freq_dB']>.05}
bad3={(x['B_GHz'],x['R_dB']):x['dJ_param_dB'] for x in r['3_param_grid']['rows'] if x['dJ_param_dB'] is None or x['dJ_param_dB']>.1}
d=pd.read_csv(ROOT/'data/digitized/liu2022/curves.csv',dtype={'series':str})
fine=grid(d,.025);sub=np.zeros(J.shape,bool);sub[::2,::2]=True
for key in sorted(set(bad2)|set(bad3),key=lambda k:(str(k[0]),k[1])):
    b=[label(x) for x in BANDS].index(key[0]);ok=feasible(fine['P'][b],fine['c2'][b],key[1])
    lab,n=ndimage.label(ok,structure=np.ones((3,3)))
    i,k=np.unravel_index(np.argmax(np.where(ok,J,-np.inf)),J.shape);comp=lab==lab[i,k]
    print({'B_GHz':key[0],'R_dB':key[1],'crit2_dJ_freq':bad2.get(key),'crit3_dJ_param':bad3.get(key),'components':int(n),
           'optimum_component_points':int(comp.sum()),'of_which_on_coarse_subgrid':int((comp&sub).sum())})
