"""A06: read Fig.2(b) and the 50-200 GHz part of Fig.3(a) main panel. See A06 BOUNDARY.

Usage: python scripts/digitize_a06.py --images ../tmp/pdfs/a03
Same per-column cluster method and uncertainty rule as A03; colour thresholds are
imported from digitize_a03 so both experiments use identical colour definitions.
Writes new files only; A03's curves.csv is never modified.
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
OUT=ROOT/'experiments/track_A_reproduction/A06_fig2b_fig3a_main_digitization'
# Locked before running (A06 BOUNDARY): native pixels -> physical axes.
FIG2B=dict(image='p2_2.jpeg',x_px=(855,1360),x_phys=(0,4),y_px=(526,106),rows=(110,522))
FIG2B_SERIES=[  # series, colour, physical y axis, x column range, excluded column ranges
 ('nm_at_100GHz','blue',(2.2,2.3),'1',(866,1356),[(1031,1059),(1165,1244)]),
 ('bw3dB_50ohm_GHz','red',(160,220),'GHz',(859,1349),[(1028,1057),(1186,1264)])]
STARS=[('nm_at_100GHz','blue',(2.2,2.3),(300,332)),('bw3dB_50ohm_GHz','red',(160,220),(208,248))]
FIG3A=dict(image='p3_0.png',x_px=(195,1111),x_phys=(0,200),y_px=(400,101),y_phys=(0,1.5),rows=(296,397),columns=(432,1108))
LOADS=[('20','orange'),('40','red'),('50','black'),('80','blue'),('open','green')]
OPTICAL_INDEX_ROWS=(300,332)


def to_phys(px,px_range,phys_range):
    (p0,p1),(v0,v1)=px_range,phys_range
    return v0+(px-p0)*(v1-v0)/(p1-p0)


def extract(rgb,color,columns,rows,excluded,x_px,x_phys,y_px,y_phys,strict=False):
    mask=color_mask(rgb,color,strict)
    out=[];missing=ambiguous=masked=0
    dx=abs((x_phys[1]-x_phys[0])/(x_px[1]-x_px[0]));dy=abs((y_phys[1]-y_phys[0])/(y_px[1]-y_px[0]))
    for x in range(columns[0],columns[1]+1):
        if any(a<=x<=b for a,b in excluded):masked+=1;continue
        ys=np.flatnonzero(mask[rows[0]:rows[1]+1,x])+rows[0]
        groups=[g for g in np.split(ys,np.flatnonzero(np.diff(ys)>1)+1) if len(g)] if len(ys) else []
        if not groups:missing+=1;continue
        if len(groups)>1:ambiguous+=1;continue
        g=groups[0];y=float(np.median(g))
        out.append({'x_pixel':x,'y_pixel':y,'x_value':to_phys(x,x_px,x_phys),'value':to_phys(y,y_px,y_phys),
                    'uncertainty_x':2*dx,'uncertainty_value':(2+(g[-1]-g[0])/2)*dy})
    return out,{'missing_columns':missing,'ambiguous_columns':ambiguous,'masked_columns':masked}


def star_center(rgb,color,rows):
    """Centroid of the filled star; the curve passes through the star's own centre."""
    m=color_mask(rgb,color)[rows[0]:rows[1]+1,1025:1064]
    ys,xs=np.nonzero(m)
    return float(np.mean(xs)+1025),float(np.mean(ys)+rows[0]),int(len(xs))


def run(images):
    records=[];quality=[];derived={};sources={}
    fig,axs=plt.subplots(2,1,figsize=(10,9),layout='constrained')
    # Fig.2(b)
    p=images/FIG2B['image'];rgb=np.array(Image.open(p).convert('RGB'))
    sources[FIG2B['image']]=hashlib.sha256(p.read_bytes()).hexdigest()
    overlay=Image.fromarray(rgb).convert('RGB');draw=ImageDraw.Draw(overlay)
    dark=rgb.max(axis=2)<100
    rows=[y for y in range(*OPTICAL_INDEX_ROWS) if dark[y,870:1340].sum()>100]
    derived['optical_index_dashed_line']={'rows':rows,'value':to_phys(float(np.mean(rows)),FIG2B['y_px'],(2.2,2.3)),
        'uncertainty':2*0.1/420}
    ax2=axs[0].twinx()
    for (series,color,yphys,unit,cols,excl),ax in zip(FIG2B_SERIES,[axs[0],ax2]):
        rows_,q=extract(rgb,color,cols,FIG2B['rows'],excl,FIG2B['x_px'],FIG2B['x_phys'],FIG2B['y_px'],yphys)
        strict,_=extract(rgb,color,cols,FIG2B['rows'],excl,FIG2B['x_px'],FIG2B['x_phys'],FIG2B['y_px'],yphys,True)
        lookup={r['x_pixel']:r['y_pixel'] for r in strict}
        diffs=[abs(r['y_pixel']-lookup[r['x_pixel']]) for r in rows_ if r['x_pixel'] in lookup]
        eligible=cols[1]-cols[0]+1-q['masked_columns']
        q.update(panel='fig2b',series=series,retained=len(rows_),eligible_columns=eligible,retained_fraction=len(rows_)/eligible,
                 strict_common_points=len(diffs),strict_median_difference_px=float(np.median(diffs)) if diffs else None,
                 excluded_column_ranges=excl,max_retained_gap_um=float(np.max(np.diff([r['x_value'] for r in rows_]))))
        quality.append(q)
        for r in rows_:
            r.update(panel='fig2b',series=series,x_quantity='BCB_thickness_um',unit=unit,provenance='digitized_paper_calculation',source=FIG2B['image'])
            records.append(r)
            if r['x_pixel']%4==0:draw.ellipse((r['x_pixel']-1,r['y_pixel']-1,r['x_pixel']+1,r['y_pixel']+1),fill=(255,0,255))
        ax.plot([r['x_value'] for r in rows_],[r['value'] for r in rows_],'.',ms=2,color=color,label=series)
    for series,color,yphys,rows_ in STARS:
        x,y,n=star_center(rgb,color,rows_)
        derived['star_'+series]={'x_pixel':x,'y_pixel':y,'pixels':n,'BCB_um':to_phys(x,FIG2B['x_px'],FIG2B['x_phys']),
            'value':to_phys(y,FIG2B['y_px'],yphys),'uncertainty_value':2*abs(yphys[1]-yphys[0])/420,'uncertainty_BCB_um':2*4/505}
        draw.rectangle((x-3,y-3,x+3,y+3),outline=(0,200,0))
    overlay.save(images/'fig2b_overlay.png')
    axs[0].set_title('Fig.2(b) | digitized author calculation (stars masked from curves)');axs[0].set_xlabel('BCB thickness (um)')
    axs[0].set_ylabel('nm @ 100 GHz');ax2.set_ylabel('3-dB EO bandwidth, 50 ohm (GHz)');axs[0].grid(alpha=.2)
    # Fig.3(a) main panel, 50-200 GHz
    p=images/FIG3A['image'];rgb=np.array(Image.open(p).convert('RGB'))
    sources[FIG3A['image']]=hashlib.sha256(p.read_bytes()).hexdigest()
    overlay=Image.fromarray(rgb).convert('RGB');draw=ImageDraw.Draw(overlay)
    for series,color in LOADS:
        rows_,q=extract(rgb,color,FIG3A['columns'],FIG3A['rows'],[],FIG3A['x_px'],FIG3A['x_phys'],FIG3A['y_px'],FIG3A['y_phys'])
        strict,_=extract(rgb,color,FIG3A['columns'],FIG3A['rows'],[],FIG3A['x_px'],FIG3A['x_phys'],FIG3A['y_px'],FIG3A['y_phys'],True)
        lookup={r['x_pixel']:r['y_pixel'] for r in strict}
        diffs=[abs(r['y_pixel']-lookup[r['x_pixel']]) for r in rows_ if r['x_pixel'] in lookup]
        xs=[r['x_value'] for r in rows_]
        q.update(panel='fig3a_main',series=series,retained=len(rows_),eligible_columns=FIG3A['columns'][1]-FIG3A['columns'][0]+1,
                 retained_fraction=len(rows_)/(FIG3A['columns'][1]-FIG3A['columns'][0]+1),strict_common_points=len(diffs),
                 strict_median_difference_px=float(np.median(diffs)) if diffs else None,
                 max_retained_gap_GHz=float(np.max(np.diff(xs))) if len(xs)>1 else None,usable_for_scoring=len(rows_)>=100)
        quality.append(q)
        for r in rows_:
            r.update(panel='fig3a_main',series=series,x_quantity='frequency_GHz',unit='V',provenance='digitized_paper_calculation',source=FIG3A['image'])
            records.append(r)
            if r['x_pixel']%4==0:draw.ellipse((r['x_pixel']-1,r['y_pixel']-1,r['x_pixel']+1,r['y_pixel']+1),fill=(255,0,255))
        axs[1].plot(xs,[r['value'] for r in rows_],'.',ms=2,color='goldenrod' if color=='orange' else color,label=series)
    overlay.save(images/'fig3a_main_overlay.png')
    axs[1].set_title('Fig.3(a) main panel, 50-200 GHz | digitized author calculation');axs[1].set_xlabel('Frequency (GHz)')
    axs[1].set_ylabel('Average voltage (V)');axs[1].legend(ncol=5,fontsize=8);axs[1].grid(alpha=.2)
    fig.suptitle('A06 | Raster extraction only; missing/ambiguous/masked columns omitted',fontsize=13)
    fig.savefig(OUT/'digitized_preview.png',dpi=150);plt.close(fig)
    cols=['x_pixel','y_pixel','x_value','value','uncertainty_x','uncertainty_value','panel','series','x_quantity','unit','provenance','source']
    pd.DataFrame(records)[cols].to_csv(DATA/'a06_curves.csv',index=False)
    metadata={'doi':'10.1109/LPT.2022.3178214','sources':sources,
        'calibration':{'fig2b':{k:v for k,v in FIG2B.items()},'fig2b_series':[{'series':s,'colour':c,'y_physical':y,'unit':u,'columns':cl,'excluded':ex} for s,c,y,u,cl,ex in FIG2B_SERIES],'fig3a_main':FIG3A},
        'method':'A03 per-column unique contiguous colour cluster; stars and annotations masked; no interpolation',
        'uncertainty':'2px axis localization plus half stroke height; not author error bars','quality':quality,'derived':derived}
    (DATA/'a06_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    # Predeclared A06 criteria 2-4 and reported-only item 6; visual QA (5) is recorded in RESULTS.md.
    fig2b=[q for q in quality if q['panel']=='fig2b']
    criteria={'2_optical_index_line_2.250pm0.001':abs(derived['optical_index_dashed_line']['value']-2.25)<=0.001,
        '3_fig2b_retained_ge_70pct':all(q['retained_fraction']>=.7 for q in fig2b),
        '4_strict_median_le_1px':all(q['strict_median_difference_px'] is None or q['strict_median_difference_px']<=1 for q in quality),
        '6_fig3a_main_usable_series':[q['series'] for q in quality if q['panel']=='fig3a_main' and q['usable_for_scoring']]}
    result={'grade':'PARTIAL_EXTRACTION_COMPLETE','scope':'raster readout only; no model comparison',
        'a06_curves_sha256':hashlib.sha256((DATA/'a06_curves.csv').read_bytes()).hexdigest(),'criteria':criteria,'derived':derived}
    (OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    return quality,derived


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--images',type=Path,required=True)
    q,d=run(ap.parse_args().images)
    print(json.dumps({'quality':q,'derived':d},indent=1))
