"""M03: read the measured (solid) curves of Fig.7: EO response for 38/57/83 ohm and open, S11 for
38/57/83 ohm. See M03 BOUNDARY.

Usage: python scripts/digitize_m03.py --images ../tmp/pdfs/a03
Per-column cluster, merge-across-occluder and clipping rules of M01; orange/red/blue thresholds from
digitize_a03, purple defined from the legend colour. In Fig.7(a) the author's dashed calculations
share the colours of the measured curves, so only pixels in wide (>=30 px) connected components,
after bridging across occluders, are read; the dashes are never read. Evidence class
digitized_paper_measurement (SC-03). Writes new files only.
"""
from pathlib import Path
import argparse,hashlib,json,sys
import numpy as np
import pandas as pd
from PIL import Image,ImageDraw
from scipy import ndimage
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parent))
from digitize_a03 import color_mask
from digitize_m01 import to_phys,groups_of

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/digitized/liu2022'
OUT=ROOT/'experiments/track_M_measurement/M03_fig7_digitization'
IMAGE='p4_1.png'
# Locked before running (M03 BOUNDARY): native pixels -> physical axes; frames; scan windows.
X_PX=(206.5,621);X_PHYS=(20,100)
PANELS={
 'fig7a':dict(x_px=X_PX,x_phys=X_PHYS,y_px=(117.5,510.5),y_phys=(0,-18),frame=(104,672,74,510),scan_x=(108,372),
   quantity='EO_response',x_ticks_GHz=20,y_tick_step=3),
 'fig7b':dict(x_px=X_PX,x_phys=X_PHYS,y_px=(649.5,1085.5),y_phys=(0,-35),frame=(104,672,650,1085),scan_x=(108,669),
   quantity='S11',x_ticks_GHz=20,y_tick_step=5)}
LEGEND={'fig7a':{c:[(329,490,370,484)] for c in ('purple','orange','red','blue')},
        'fig7b':{'purple':[(426,490,681,690)],'orange':[(426,490,712,721)],'red':[(426,490,742,751)]}}
COLOURS=('purple','orange','red','blue')
SERIES=[('fig7a','38ohm','purple'),('fig7a','57ohm','orange'),('fig7a','83ohm','red'),('fig7a','open','blue'),
        ('fig7b','38ohm','purple'),('fig7b','57ohm','orange'),('fig7b','83ohm','red')]
BRIDGE_PX=3;SOLID_MIN_EXTENT=30
REF_LINES={'-1':(130,150),'-3':(175,192)};REF_COLS=(400,660);REF_MIN=80


def cmask(rgb,colour,strict=False):
    if colour!='purple':return color_mask(rgb,colour,strict)
    r,g,b=np.moveaxis(rgb.astype(int),-1,0)
    if strict:return (r>=105)&(r<=150)&(b>=125)&(b<=165)&(g<=70)
    return (r>=90)&(r<=200)&(b>=110)&(g<=110)&(b-g>=50)&(r-g>=40)


def grey(rgb):
    a=rgb.astype(int);return ((a.max(axis=2)-a.min(axis=2))<40)&(a.max(axis=2)<215)


def masked(rgb,panel,colour,strict):
    m=cmask(rgb,colour,strict).copy()
    for x0,x1,y0,y1 in LEGEND[panel].get(colour,[]):m[y0:y1+1,x0:x1+1]=False
    return m


def shift(m,dy,dx):
    """o[y,x] = m[y-dy, x-dx], zero outside."""
    o=np.zeros_like(m);H,W=m.shape
    o[max(dy,0):H+min(dy,0),max(dx,0):W+min(dx,0)]=m[max(-dy,0):H+min(-dy,0),max(-dx,0):W+min(-dx,0)]
    return o


def solid_only(m,occ,window):
    """Bridge across occluders, keep connected components at least SOLID_MIN_EXTENT px wide."""
    x0,x1,y0,y1=window;sub=m[y0:y1+1,x0:x1+1];o=occ[y0:y1+1,x0:x1+1];K=range(1,BRIDGE_PX+1)
    up=np.logical_or.reduce([shift(sub,k,0) for k in K]);dn=np.logical_or.reduce([shift(sub,-k,0) for k in K])
    lf=np.logical_or.reduce([shift(sub,0,k) for k in K]);rt=np.logical_or.reduce([shift(sub,0,-k) for k in K])
    lab,_=ndimage.label(sub|(o&((up&dn)|(lf&rt))),structure=np.ones((3,3)))
    keep=[i+1 for i,s in enumerate(ndimage.find_objects(lab)) if s[1].stop-s[1].start>=SOLID_MIN_EXTENT]
    solid=np.isin(lab,keep)&sub
    out=np.zeros_like(m);out[y0:y1+1,x0:x1+1]=solid
    cols=np.flatnonzero(solid.any(axis=0))
    return out,int((sub&~solid).sum()),(int(cols[-1])+x0 if len(cols) else None)


def extract(rgb,panel,colour,strict=False):
    P=PANELS[panel];xl,xr,yt,yb=P['frame'];r0,r1=yt+4,yb-3;c0,c1=P['scan_x']
    m=masked(rgb,panel,colour,strict)
    occ=grey(rgb)|np.logical_or.reduce([cmask(rgb,c,strict) for c in COLOURS if c!=colour])
    dash_px=None;end_col=c1
    if panel=='fig7a':
        m,dash_px,end_col=solid_only(m,occ,(c0,c1,r0,r1))
        if end_col is None:end_col=c0-1
    dy=abs((P['y_phys'][1]-P['y_phys'][0])/(P['y_px'][1]-P['y_px'][0]))
    dx=abs((P['x_phys'][1]-P['x_phys'][0])/(P['x_px'][1]-P['x_px'][0]))
    out=[];clipped=[];count=dict(single=0,merged=0,missing=0,ambiguous=0,clipped=0)
    for x in range(c0,end_col+1):
        gs=groups_of(np.flatnonzero(m[r0:r1+1,x])+r0)
        if len(gs)==1:g=gs[0];flag='single'
        elif len(gs)>1 and all(occ[a[-1]+1:b[0],x].all() for a,b in zip(gs[:-1],gs[1:])):
            g=np.arange(gs[0][0],gs[-1][-1]+1);flag='merged'
        elif len(gs)>1:count['ambiguous']+=1;continue
        else:count['missing']+=1;continue
        if g[0]<=r0 or g[-1]>=r1:
            count['clipped']+=1;clipped.append({'x_pixel':x,'x_value':to_phys(x,P['x_px'],P['x_phys']),'edge':'top' if g[0]<=r0 else 'bottom'});continue
        count[flag]+=1;y=float(np.median(g))
        out.append({'x_pixel':x,'y_pixel':y,'x_value':to_phys(x,P['x_px'],P['x_phys']),'value':to_phys(y,P['y_px'],P['y_phys']),
                    'uncertainty_x':2*dx,'uncertainty_value':(2+(g[-1]-g[0])/2)*dy,'flag':flag})
    return out,clipped,count,end_col-c0+1,dash_px,end_col,m


def runs(idx):
    idx=np.asarray(idx)
    return [(int(a[0]),int(a[-1])) for a in np.split(idx,np.flatnonzero(np.diff(idx)>1)+1)] if len(idx) else []


def reference_lines(rgb):
    P=PANELS['fig7a'];gr=grey(rgb);out={}
    for label,(y0,y1) in REF_LINES.items():
        rows=[y for y in range(y0,y1+1) if gr[y,REF_COLS[0]:REF_COLS[1]+1].sum()>REF_MIN]
        out[label]={'rows':rows,'value':to_phys(float(np.mean(rows)),P['y_px'],P['y_phys']) if rows else None,
                    'uncertainty':2*abs(18/(P['y_px'][1]-P['y_px'][0]))}
    return out


def ticks(rgb,ref_rows):
    """Detect tick centres (grey runs 2-10 px inside the frame) and their distance to the calibration grid."""
    gr=grey(rgb);report={}
    for panel,P in PANELS.items():
        xl,xr,yt,yb=P['frame']
        sx=(P['x_px'][1]-P['x_px'][0])/(P['x_phys'][1]-P['x_phys'][0])
        x_grid=[P['x_px'][0]+(f-P['x_phys'][0])*sx for f in range(0,111,P['x_ticks_GHz'])]
        sy=(P['y_px'][1]-P['y_px'][0])/(P['y_phys'][1]-P['y_phys'][0])
        y_grid=[P['y_px'][0]+(v-P['y_phys'][0])*sy for v in np.arange(P['y_phys'][0],P['y_phys'][1]-1e-9,-P['y_tick_step'])]
        found=[]
        for band in (gr[yb-10:yb-1,xl+3:xr-2],gr[yt+2:yt+11,xl+3:xr-2]):
            for a,b in runs(np.flatnonzero(band.sum(axis=0)>=3)+xl+3):
                c=(a+b)/2;found.append(('x',c,min(abs(c-g) for g in x_grid)))
        excluded=[]
        for band in (gr[yt+3:yb-2,xl+2:xl+11],gr[yt+3:yb-2,xr-10:xr-1]):
            for a,b in runs(np.flatnonzero(band.sum(axis=1)>=3)+yt+3):
                c=(a+b)/2
                if panel=='fig7a' and any(abs(c-np.mean(r))<=3 for r in ref_rows if r):excluded.append(c);continue
                found.append(('y',c,min(abs(c-g) for g in y_grid)))
        report[panel]={'ticks':[{'axis':ax,'centre_px':c,'residual_px':r} for ax,c,r in found],
                       'excluded_reference_line_rows':excluded,'max_residual_px':max(r for _,_,r in found) if found else None}
    return report


def run(images):
    p=images/IMAGE;raw=p.read_bytes();rgb=np.array(Image.open(p).convert('RGB'))
    pdf=images.parent/'reference_tfln_lpt2022.pdf'
    overlay=Image.fromarray(rgb).convert('RGB');draw=ImageDraw.Draw(overlay)
    records=[];quality=[];clipped_all={}
    fig,axs=plt.subplots(1,2,figsize=(13,4.8),layout='constrained')
    plot_colour={'purple':'purple','orange':'goldenrod','red':'red','blue':'blue'}
    for panel,load,colour in SERIES:
        rows_,clipped,count,eligible,dash_px,end_col,solid=extract(rgb,panel,colour)
        strict=extract(rgb,panel,colour,True)[0]
        lookup={r['x_pixel']:r['y_pixel'] for r in strict}
        diffs=[abs(r['y_pixel']-lookup[r['x_pixel']]) for r in rows_ if r['x_pixel'] in lookup]
        xs=[r['x_value'] for r in rows_];P=PANELS[panel];series=f"{P['quantity']}_{load}"
        q=dict(panel=panel,series=series,colour=colour,eligible_columns=eligible,retained=len(rows_),
               retained_fraction=len(rows_)/eligible if eligible>0 else 0.0,**count,dash_classified_px=dash_px,
               end_column=end_col,end_GHz=to_phys(end_col,P['x_px'],P['x_phys']),strict_common_points=len(diffs),
               strict_median_difference_px=float(np.median(diffs)) if diffs else None,
               max_retained_gap_GHz=float(np.max(np.diff(xs))) if len(xs)>1 else None)
        quality.append(q);clipped_all[series]=clipped
        if panel=='fig7a':
            raw_m=masked(rgb,panel,colour,False);c0,c1=P['scan_x'];r0,r1=P['frame'][2]+4,P['frame'][3]-3
            ys,xs_=np.nonzero(raw_m[r0:r1+1,c0:c1+1]&~solid[r0:r1+1,c0:c1+1])
            for y,x in zip(ys[::2],xs_[::2]):draw.point((x+c0,y+r0),fill=(0,255,255))
        for r in rows_:
            r.update(panel=panel,series=series,load=load,x_quantity='frequency_GHz',unit='dB',
                     provenance='digitized_paper_measurement',source=IMAGE)
            records.append(r)
            if r['x_pixel']%3==0:draw.ellipse((r['x_pixel']-1,r['y_pixel']-1,r['x_pixel']+1,r['y_pixel']+1),
                fill=(255,0,255) if r['flag']=='single' else (0,200,0))
        ax=axs[0] if panel=='fig7a' else axs[1]
        ax.plot(xs,[r['value'] for r in rows_],'.',ms=2,color=plot_colour[colour],label=f'{load} ({len(rows_)})')
    refs=reference_lines(rgb)
    for v in refs.values():
        for y in v['rows']:draw.line((REF_COLS[0],y,REF_COLS[1],y),fill=(0,255,0))
    tick=ticks(rgb,[v['rows'] for v in refs.values()])
    overlay.save(images/'fig7_overlay.png')
    axs[0].set_title('Fig.7(a) | digitized measured EO response (solid curves only)');axs[0].set_ylabel('dB')
    axs[1].set_title('Fig.7(b) | digitized measured S11');axs[1].set_ylabel('dB')
    for ax in axs:ax.set_xlabel('Frequency (GHz)');ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.suptitle('M03 | Raster extraction only; dashed author calculations not read',fontsize=12)
    fig.savefig(OUT/'digitized_preview.png',dpi=150);plt.close(fig)
    cols=['x_pixel','y_pixel','x_value','value','uncertainty_x','uncertainty_value','flag','panel','series','load',
          'x_quantity','unit','provenance','source']
    pd.DataFrame(records)[cols].to_csv(DATA/'m03_curves.csv',index=False)
    metadata={'doi':'10.1109/LPT.2022.3178214','evidence':'digitized_paper_measurement','scope_change':'SC-03 (D7, D13)',
        'sources':{IMAGE:hashlib.sha256(raw).hexdigest()},'pdf_sha256':hashlib.sha256(pdf.read_bytes()).hexdigest() if pdf.exists() else None,
        'extraction':'pdfimages -all -f 4 -l 4 reference_tfln_lpt2022.pdf (image 1 of page 4)',
        'calibration':PANELS,'legend_masks':LEGEND,'series':[dict(zip(('panel','load','colour'),s)) for s in SERIES],
        'rules':{'purple':'90<=R<=200, B>=110, G<=110, B-G>=50, R-G>=40 (strict 105<=R<=150, 125<=B<=165, G<=70)',
                 'occluder':'grey (max-min<40, max<215) or another curve colour','bridge_px':BRIDGE_PX,
                 'solid_min_extent_px':SOLID_MIN_EXTENT,'scan':'frame+4 .. frame-3; fig7a to x=372 (52 GHz)'},
        'method':'M01 per-column cluster with merge across occluders; fig7a solid-component filter; no interpolation',
        'uncertainty':'2px axis localization plus half stroke height; not author error bars',
        'quality':quality,'clipped':clipped_all,'reference_lines':refs,'ticks':tick}
    (DATA/'m03_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    criteria={'2_reference_lines_pm0.09dB':all(v['value'] is not None and abs(v['value']-float(k))<=0.09 for k,v in refs.items()),
        '3_tick_residual_le_1px':all(t['max_residual_px'] is not None and t['max_residual_px']<=1 for t in tick.values()),
        '4_retained_ge_70pct':{q['series']:q['retained_fraction']>=.7 for q in quality},
        '5_strict_median_le_1px':{q['series']:(q['strict_median_difference_px'] is None or q['strict_median_difference_px']<=1) for q in quality}}
    result={'grade':'PARTIAL_EXTRACTION_COMPLETE','scope':'raster readout of measured solid curves only; dashed calculations not read',
        'evidence':'digitized_paper_measurement','m03_curves_sha256':hashlib.sha256((DATA/'m03_curves.csv').read_bytes()).hexdigest(),
        'criteria':criteria,'reference_lines':refs,'tick_max_residual_px':{k:v['max_residual_px'] for k,v in tick.items()},'quality':quality}
    (OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--images',type=Path,required=True)
    print(json.dumps(run(ap.parse_args().images),indent=1))
