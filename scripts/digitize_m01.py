"""M01: read Fig.6(a) measured S21/S11 and Fig.6(b) extracted microwave index. See M01 BOUNDARY.

Usage: python scripts/digitize_m01.py --images ../tmp/pdfs/a03
Same per-column cluster method and uncertainty rule as A03/A06, colour thresholds imported
from digitize_a03. Two rules predeclared in the M01 boundary are added: (i) clusters of one
colour separated only by a declared occluder are merged ('merged'); (ii) a Fig.6(b) column
with no blue pixels but a red optical line is 'occluded' (blue indistinguishable from the
red line). Evidence class digitized_paper_measurement (SC-03). Writes new files only; the
A03/A06 data are never modified.
"""
from pathlib import Path
import argparse,hashlib,json,sys
import numpy as np
import pandas as pd
from PIL import Image,ImageDraw
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parent))
from digitize_a03 import color_mask

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/digitized/liu2022'
OUT=ROOT/'experiments/track_M_measurement/M01_fig6_digitization'
IMAGE='p4_0.jpeg'
# Locked before running (M01 BOUNDARY): native pixels -> physical axes; frame edges in pixels.
PANELS={
 'fig6a':dict(x_px=(209,837.5),x_phys=(0,150),y_px=(113.5,661),y_phys=(0,-36),frame=(209,921,113,661),
   legend={'blue':(541,601,250,255),'red':(717,775,582,586)}),
 'fig6b':dict(x_px=(1156,1790.5),x_phys=(0,150),y_px=(660.5,113),y_phys=(2.0,3.0),frame=(1156,1875,113,660),
   legend={'blue':(1471,1531,288,293),'red':(1472,1531,336,340)})}
SERIES=[  # panel, series, colour, unit, declared occluder
 ('fig6a','S21','blue','dB','dark'),
 ('fig6a','S11','red','dB',None),
 ('fig6b','nm_microwave','blue','1','red'),
 ('fig6b','optical_mode_line','red','1',None)]
LEGEND_PAD=3
DARK_MAX=140      # occluder 'dark': dashed line and arrow, max(RGB)<140
DASH_DARK_MAX=100 # criterion 2 dashed-line rows, same darkness test as A06's dashed line
DASH_ROWS=(195,225);DASH_COLS=(215,915);DASH_MIN=150


def to_phys(px,px_range,phys_range):
    (p0,p1),(v0,v1)=px_range,phys_range
    return v0+(px-p0)*(v1-v0)/(p1-p0)


def groups_of(ys):
    return [g for g in np.split(ys,np.flatnonzero(np.diff(ys)>1)+1) if len(g)] if len(ys) else []


def masks(rgb,panel,color,strict):
    m=color_mask(rgb,color,strict).copy()
    for c,(x0,x1,y0,y1) in PANELS[panel]['legend'].items():
        if c==color:m[y0-LEGEND_PAD:y1+LEGEND_PAD+1,x0-LEGEND_PAD:x1+LEGEND_PAD+1]=False
    return m


def occluder_mask(rgb,panel,name,strict):
    if name=='dark':return rgb.max(axis=2)<DARK_MAX
    if name=='red':return masks(rgb,panel,'red',strict)
    return None


def extract(rgb,panel,color,occluder,strict=False):
    P=PANELS[panel];xl,xr,yt,yb=P['frame'];r0,r1=yt+4,yb-3
    mask=masks(rgb,panel,color,strict);occ=occluder_mask(rgb,panel,occluder,strict) if occluder else None
    red=masks(rgb,panel,'red',strict) if (panel=='fig6b' and color=='blue') else None
    dy=abs((P['y_phys'][1]-P['y_phys'][0])/(P['y_px'][1]-P['y_px'][0]))
    dx=abs((P['x_phys'][1]-P['x_phys'][0])/(P['x_px'][1]-P['x_px'][0]))
    out=[];clipped=[];count=dict(single=0,merged=0,occluded=0,missing=0,ambiguous=0,clipped=0)
    for x in range(xl+4,xr-2):
        gs=groups_of(np.flatnonzero(mask[r0:r1+1,x])+r0)
        flag=None
        if len(gs)==1:g=gs[0];flag='single'
        elif len(gs)>1 and occ is not None and all(occ[a[-1]+1:b[0],x].all() for a,b in zip(gs[:-1],gs[1:])):
            g=np.arange(gs[0][0],gs[-1][-1]+1);flag='merged'
        elif len(gs)>1:count['ambiguous']+=1;continue
        elif red is not None:
            rg=groups_of(np.flatnonzero(red[r0:r1+1,x])+r0)
            if len(rg)==1:g=rg[0];flag='occluded'
            else:count['missing']+=1;continue
        else:count['missing']+=1;continue
        if g[0]<=r0 or g[-1]>=r1:
            count['clipped']+=1;clipped.append({'x_pixel':x,'x_value':to_phys(x,P['x_px'],P['x_phys']),'edge':'top' if g[0]<=r0 else 'bottom'});continue
        count[flag]+=1;y=float(np.median(g))
        out.append({'x_pixel':x,'y_pixel':y,'x_value':to_phys(x,P['x_px'],P['x_phys']),'value':to_phys(y,P['y_px'],P['y_phys']),
                    'uncertainty_x':2*dx,'uncertainty_value':(2+(g[-1]-g[0])/2)*dy,'flag':flag})
    return out,clipped,count,xr-2-(xl+4)


def dashed_line(rgb):
    P=PANELS['fig6a'];dark=rgb.max(axis=2)<DASH_DARK_MAX
    rows=[y for y in range(*DASH_ROWS) if dark[y,DASH_COLS[0]:DASH_COLS[1]].sum()>DASH_MIN]
    dy=abs(36/(P['y_px'][1]-P['y_px'][0]))
    return {'rows':rows,'value':to_phys(float(np.mean(rows)),P['y_px'],P['y_phys']) if rows else None,'uncertainty':2*dy}


def run(images):
    p=images/IMAGE;raw=p.read_bytes();rgb=np.array(Image.open(p).convert('RGB'))
    sources={IMAGE:hashlib.sha256(raw).hexdigest()}
    pdf=images.parent/'reference_tfln_lpt2022.pdf'
    overlay=Image.fromarray(rgb).convert('RGB');draw=ImageDraw.Draw(overlay)
    records=[];quality=[];clipped_all={}
    fig,axs=plt.subplots(1,2,figsize=(13,4.8),layout='constrained')
    for panel,series,color,unit,occ in SERIES:
        rows_,clipped,count,eligible=extract(rgb,panel,color,occ)
        strict,_,_,_=extract(rgb,panel,color,occ,True)
        lookup={r['x_pixel']:r['y_pixel'] for r in strict}
        diffs=[abs(r['y_pixel']-lookup[r['x_pixel']]) for r in rows_ if r['x_pixel'] in lookup]
        xs=[r['x_value'] for r in rows_]
        q=dict(panel=panel,series=series,colour=color,occluder=occ,eligible_columns=eligible,retained=len(rows_),
               retained_fraction=len(rows_)/eligible,**count,strict_common_points=len(diffs),
               strict_median_difference_px=float(np.median(diffs)) if diffs else None,
               max_retained_gap_GHz=float(np.max(np.diff(xs))) if len(xs)>1 else None)
        quality.append(q);clipped_all[f'{panel}:{series}']=clipped
        for r in rows_:
            r.update(panel=panel,series=series,x_quantity='frequency_GHz',unit=unit,provenance='digitized_paper_measurement',source=IMAGE)
            records.append(r)
            if r['x_pixel']%3==0:draw.ellipse((r['x_pixel']-1,r['y_pixel']-1,r['x_pixel']+1,r['y_pixel']+1),
                fill=(255,0,255) if r['flag']=='single' else (0,200,0))
        ax=axs[0] if panel=='fig6a' else axs[1]
        for flag,mk in (('single','.'),('merged','x'),('occluded','+')):
            sel=[r for r in rows_ if r['flag']==flag]
            if sel:ax.plot([r['x_value'] for r in sel],[r['value'] for r in sel],mk,ms=2.5,color=color,
                           alpha=1 if flag=='single' else .5,label=f'{series} ({flag}, {len(sel)})')
    dash=dashed_line(rgb)
    for y in dash['rows']:draw.line((215,y,915,y),fill=(0,255,255))
    overlay.save(images/'fig6_overlay.png')
    axs[0].axhline(-6.4,color='k',ls='--',lw=.8)
    axs[0].set_title('Fig.6(a) | digitized measurement (author probes, open terminal)');axs[0].set_ylabel('dB')
    axs[1].set_title('Fig.6(b) | digitized extracted microwave index');axs[1].set_ylabel('index')
    for ax in axs:ax.set_xlabel('Frequency (GHz)');ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.suptitle('M01 | Raster extraction only; missing/ambiguous/clipped columns omitted',fontsize=12)
    fig.savefig(OUT/'digitized_preview.png',dpi=150);plt.close(fig)
    cols=['x_pixel','y_pixel','x_value','value','uncertainty_x','uncertainty_value','flag','panel','series','x_quantity','unit','provenance','source']
    pd.DataFrame(records)[cols].to_csv(DATA/'m01_curves.csv',index=False)
    optical=[r['value'] for r in records if r['series']=='optical_mode_line' and r['flag']=='single']
    derived={'dashed_line_minus_6p4':dash,'optical_mode_line_median':float(np.median(optical)) if optical else None,
             'optical_mode_line_columns':len(optical),'optical_mode_line_uncertainty':2/547.5}
    metadata={'doi':'10.1109/LPT.2022.3178214','evidence':'digitized_paper_measurement','scope_change':'SC-03 (D7)',
        'sources':sources,'pdf_sha256':hashlib.sha256(pdf.read_bytes()).hexdigest() if pdf.exists() else None,
        'extraction':'pdfimages -all -f 4 -l 4 reference_tfln_lpt2022.pdf (native JPEG stream, image 0 of page 4)',
        'calibration':PANELS,'series':[dict(zip(('panel','series','colour','unit','occluder'),s)) for s in SERIES],
        'rules':{'legend_pad_px':LEGEND_PAD,'dark_occluder_max_rgb':DARK_MAX,'scan':'frame+4 .. frame-3',
                 'merged':'same-colour clusters separated only by declared occluder pixels',
                 'occluded':'fig6b blue absent, unique red optical-line cluster; value = red centre',
                 'clipped':'cluster touches the scan window edge; counted, not a value'},
        'method':'A03 per-column unique contiguous colour cluster; legends masked by pixel rectangles; no interpolation',
        'uncertainty':'2px axis localization plus half stroke height; not author error bars',
        'quality':quality,'clipped':clipped_all,'derived':derived}
    (DATA/'m01_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    # Predeclared criteria 2-5; visual QA (6) and evidence (1, 7) are recorded in RESULTS.md.
    measured=[q for q in quality if q['series']!='optical_mode_line']
    criteria={'2_dashed_line_-6.40pm0.13dB':dash['value'] is not None and abs(dash['value']+6.40)<=0.13,
        '3_optical_line_2.250pm0.004':derived['optical_mode_line_median'] is not None and abs(derived['optical_mode_line_median']-2.25)<=0.004,
        '4_retained_ge_70pct':{q['series']:q['retained_fraction']>=.7 for q in measured},
        '5_strict_median_le_1px':{q['series']:(q['strict_median_difference_px'] is None or q['strict_median_difference_px']<=1) for q in quality}}
    result={'grade':'PARTIAL_EXTRACTION_COMPLETE','scope':'raster readout of measured curves only; no model comparison',
        'evidence':'digitized_paper_measurement','m01_curves_sha256':hashlib.sha256((DATA/'m01_curves.csv').read_bytes()).hexdigest(),
        'criteria':criteria,'derived':derived,'quality':quality}
    (OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--images',type=Path,required=True)
    print(json.dumps(run(ap.parse_args().images),indent=1))
