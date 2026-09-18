"""Read native figure rasters, independently of model output. See A03 boundary.

Usage: python scripts/digitize_a03.py --images ../tmp/pdfs/a03
Rasters stay local; committed CSVs + metadata can be read without this dependency.
"""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
import pandas as pd
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/digitized/liu2022'
OUT=ROOT/'experiments/track_A_reproduction/A03_digitization'
PANELS=[
 ('fig2a','p2_2.jpeg',(120,687,526,106),(0,200,0,1.5),'dB/mm',[('this_work','red')]),
 ('fig2c','p2_2.jpeg',(495,1035,1150,732),(0,200,0,100),'ohm',[('Z0','black')]),
 ('fig3a_inset','p3_0.png',(650,1107,258,24),(0,50,.2,1),'V',[('20','orange'),('40','red'),('50','black'),('80','blue'),('open','green')]),
 ('fig3b','p3_0.png',(195,1111,858,522),(0,200,-12,6),'dB',[('20','orange'),('40','red'),('50','black'),('80','blue'),('open','green')]),
 ('fig3c','p3_0.png',(195,1111,1320,974),(0,200,-50,0),'dB',[('20','orange'),('40','red'),('50','black'),('80','blue')])]

def color_mask(rgb,color,strict=False):
 r,g,b=np.moveaxis(rgb.astype(float),-1,0)
 delta=100 if strict else 65
 if color=='red':return (r>160)&(r-g>delta)&(r-b>delta)&(g<120)
 if color=='orange':return (r>180)&(g>100)&(g<210)&(b<110)&(r-g>35)
 if color=='blue':return (b>150)&(b-r>delta)&(b-g>delta)
 if color=='green':return (g>70)&(g-r>20)&(g-b>35)&(r<195)
 return (np.max(rgb,axis=2)<(95 if strict else 140))

def extract(rgb,panel,pixels,physical,color,strict=False):
 x0,x1,y0,y1=pixels;f0,f1,v0,v1=physical
 mask=color_mask(rgb,color,strict)
 rows=[];missing=0;ambiguous=0
 for x in range(x0+4,x1-3):
  ys=np.flatnonzero(mask[y1+4:y0-3,x])+y1+4
  if panel=='fig2a':ys=ys[ys>280] # below legend (signal lies below this throughout)
  if panel=='fig2c':ys=ys[(ys<939)|(ys>947)];ys=ys[(ys>900)&(ys<1020)]
  if panel=='fig3c' and x>=740:ys=ys[ys<1180]
  if panel=='fig3a_inset' and x>=840:ys=ys[ys>(184 if color=='black' else 172)]
  if panel=='fig3a_inset' and color=='black':ys=ys[ys<220] # remove bottom-axis tick fragments; visible black signal is above this
  groups=np.split(ys,np.flatnonzero(np.diff(ys)>1)+1) if len(ys) else []
  groups=[g for g in groups if len(g)]
  if not groups:missing+=1;continue
  if len(groups)>1:ambiguous+=1;continue
  g=groups[0];y=float(np.median(g))
  rows.append({'x_pixel':x,'y_pixel':y,'frequency_GHz':f0+(x-x0)*(f1-f0)/(x1-x0),
   'value':v0+(y0-y)*(v1-v0)/(y0-y1),'uncertainty_frequency_GHz':2*(f1-f0)/(x1-x0),
   'uncertainty_value':(2+(g[-1]-g[0])/2)*(v1-v0)/(y0-y1)})
 return rows,missing,ambiguous

def run(images):
 DATA.mkdir(exist_ok=True,parents=True)
 fig,axs=plt.subplots(3,2,figsize=(12,11),layout='constrained');axs=axs.ravel()
 records=[];quality=[];sources={}
 for ax,(panel,filename,pixels,physical,unit,series) in zip(axs,PANELS):
  p=images/filename;rgb=np.array(Image.open(p).convert('RGB'));sources[filename]={'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'width':rgb.shape[1],'height':rgb.shape[0]}
  overlay=Image.fromarray(rgb).convert('RGB')
  from PIL import ImageDraw
  draw=ImageDraw.Draw(overlay)
  for label,color in series:
   rows,missing,ambig=extract(rgb,panel,pixels,physical,color)
   strict,_,_=extract(rgb,panel,pixels,physical,color,True)
   lookup={r['x_pixel']:r['value'] for r in strict}
   diffs=[abs(r['value']-lookup[r['x_pixel']]) for r in rows if r['x_pixel'] in lookup]
   for r in rows:
    r.update(panel=panel,series=label,unit=unit,provenance='digitized_paper_calculation',source=filename)
    records.append(r)
    if r['x_pixel']%4==0:draw.ellipse((r['x_pixel']-1,r['y_pixel']-1,r['x_pixel']+1,r['y_pixel']+1),fill=(255,0,255))
   ax.scatter([r['frequency_GHz'] for r in rows],[r['value'] for r in rows],s=2,label=label,color='goldenrod' if color=='orange' else color)
   xs=[r['frequency_GHz'] for r in rows]
   quality.append({'panel':panel,'series':label,'retained':len(rows),'missing_columns':missing,'ambiguous_columns':ambig,'max_retained_gap_GHz':float(max(np.diff(xs))) if len(xs)>1 else None,'strict_threshold_common_points':len(diffs),'strict_median_difference':float(np.median(diffs)) if diffs else None,'strict_max_difference':float(max(diffs)) if diffs else None})
  overlay.save(images/(panel+'_overlay.png'))
  ax.set_title(panel+' | digitized author calculation');ax.set_xlabel('Frequency (GHz)');ax.set_ylabel(unit);ax.legend(fontsize=8);ax.grid(alpha=.2)
 axs[-1].axis('off')
 fig.suptitle('A03 | Raster extraction only; missing/ambiguous columns omitted',fontsize=14)
 fig.savefig(OUT/'digitized_preview.png',dpi=150);plt.close(fig)
 pd.DataFrame(records).to_csv(DATA/'curves.csv',index=False)
 # PDF path may differ; image fingerprints are always present.
 pdf=images.parent/'reference_tfln_lpt2022.pdf'
 metadata={'pdf_sha256':hashlib.sha256(pdf.read_bytes()).hexdigest() if pdf.exists() else None,'doi':'10.1109/LPT.2022.3178214','sources':sources,'calibration':[{'panel':p,'image':im,'pixel_axes':px,'physical_axes':ph,'unit':u} for p,im,px,ph,u,s in PANELS],'method':'Per-column unique contiguous color cluster; masks in script; no model or interpolation','uncertainty':'2px axis localization plus half stroke height; not author error bars','quality':quality}
 (DATA/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
 print(json.dumps(quality,indent=2))

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--images',type=Path,required=True);run(ap.parse_args().images)
