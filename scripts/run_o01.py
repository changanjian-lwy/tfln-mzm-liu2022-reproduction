"""O01: optical TE0 group index from the Liu 2022 cross-section (SC-05, Track O). See O01 BOUNDARY.

G1 anisotropic slab gate -> G2 mesh ladder (N and N without T-rails), domain and difference-step checks ->
one-factor configurations -> two corners, solved directly -> G2 mesh check of the corners -> verdict.
Every section solve uses the x >= 0 half domain with a PEC symmetry plane (BOUNDARY section 5).
BUDGET.json caps are enforced; on a cap the partial record is written and the script exits 2.
results.json holds 6-significant-digit values; timings go to run_log.json.
"""
from collections import OrderedDict
from dataclasses import replace
from importlib import metadata
from pathlib import Path
import json,resource,sys,time
import numpy as np
from tfln_mzm.optical_fem import (Section,build_mesh,build_slab_mesh,domain,geometry,group_index,slab_exact,
                                  solve_slab,solve_te0,LN_FILES,SILICA,index)

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'experiments/track_O_optical/O01_optical_mode_ng'
LAM0,DELTA,DELTAS_CHECK=1.55,0.005,(0.0025,0.01)       # um (boundary 0.2 and criterion 2)
M_START,M_MIN=2.0,0.25                              # ladder start per BOUNDARY section 5
G1={'neff':1e-4,'ng':5e-4,'te':0.99,'tm':0.01}          # criterion 1
CONV={'neff':1e-4,'ng':5e-4};DOMAIN_NG,STEP_NG,TE_SPREAD=2e-4,1e-4,0.05   # criteria 2, 2b
ROUND,WIDE=(2.245,2.255),(2.20,2.30)                   # criterion 4
N=Section()
# factor -> (nominal, endpoints); lam_nm is the wavelength, every other key is a Section field (boundary 1)
FACTORS=OrderedDict([('lam_nm',(1550.0,[1530.0,1570.0])),('theta_deg',(75.0,[60.0,90.0])),
                     ('width_ref',('top',['bottom'])),('clad_dn',(0.0,[0.02])),
                     ('substrate',('quartz',['fused'])),('rail_t',(0.5,[None,0.2,0.8]))])


class BudgetStop(Exception):pass


class Budget:
    def __init__(s,caps):s.caps=caps;s.t0=time.perf_counter();s.solves=0;s.log=[]
    def elapsed(s):return time.perf_counter()-s.t0
    def rss_gb(s):return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e9  # bytes on macOS
    def before_solve(s):
        if s.elapsed()>s.caps['wall_seconds']:raise BudgetStop(f"wall {s.elapsed():.0f} s > {s.caps['wall_seconds']}")
        if s.solves>=s.caps['max_solves']:raise BudgetStop(f"solves {s.solves} reached cap")
    def mesh(s,n,label,sec):
        s.log.append({'mesh':label,'elements':n,'seconds':round(sec,2)})
        if n>s.caps['max_mesh_elements']:raise BudgetStop(f'mesh {label} has {n} elements')
    def solved(s,seconds,label):
        s.solves+=1;s.log.append({'solve':label,'seconds':round(seconds,2),'rss_gb':round(s.rss_gb(),3)})
        if seconds>s.caps['per_solve_seconds']:raise BudgetStop(f"solve took {seconds:.0f} s")
        if s.rss_gb()>s.caps['max_rss_gb']:raise BudgetStop(f'peak memory {s.rss_gb():.1f} GB')


def sig(x,n=6):
    if x is None:return None
    if isinstance(x,complex):return [sig(x.real,n),sig(x.imag,n)]
    return float(f'{x:.{n}g}')


def config(**changes):
    """(Section, wavelength um) for N with the given factor values."""
    lam=changes.pop('lam_nm',1550.0)/1000
    return replace(N,**changes),lam


def label(sec,lam):
    diff={k:v for k,v in vars(sec).items() if v!=getattr(N,k)}
    if abs(lam-LAM0)>1e-12:diff['lam_nm']=round(lam*1000,3)
    return diff


class Runner:
    def __init__(s,budget):s.b=budget;s.meshes={};s.cache={}
    def mesh(s,sec,m,d):
        key=(sec,m,d)
        if key not in s.meshes:
            t=time.perf_counter();mesh=build_mesh(sec,m,d,half=True)
            s.b.mesh(int(mesh.nelements),{'m':m,'d':d,**label(sec,LAM0)},time.perf_counter()-t);s.meshes[key]=mesh
        return s.meshes[key]
    def te0(s,sec,lam,m,d):
        key=(sec,round(lam,9),m,d)
        if key not in s.cache:
            mesh=s.mesh(sec,m,d);s.b.before_solve();r=solve_te0(sec,mesh,lam)
            if not (sec==N and abs(lam-LAM0)<1e-12 and d==1.0):r.pop('mode',None)   # keep fields only for the figure
            s.cache[key]=r;s.b.solved(r['solve_s'],{'m':m,'d':d,'lam_um':lam,**label(sec,LAM0)})
        return s.cache[key]
    def ng(s,sec,lam,m,d=1.0,delta=DELTA):
        rs=[s.te0(sec,lam+k*delta,m,d) for k in (-1,0,1)]
        found=all(r['found'] for r in rs)
        out={'m':m,'d':d,'delta_nm':round(delta*1000,4),'found':found,'n_elements':rs[1]['n_elements'],
             'te_fractions':[r.get('te_fraction') for r in rs],
             'nonfinite':sum(not np.isfinite(r['neff']) for r in rs if r['found'])}
        if found:
            n=[r['neff'].real for r in rs]
            out.update(neff=rs[1]['neff'],ng=group_index(n[0],n[1],n[2],lam,delta),loss_dB_per_cm=rs[1]['loss_dB_per_cm'],
                       te_spread=max(out['te_fractions'])-min(out['te_fractions']))
        out['ok_2b']=bool(found and out['nonfinite']==0 and out['te_spread']<=TE_SPREAD)
        return out


def slab_gate(run,m):
    """Criterion 1 (FEM part): anisotropic CLN slab vs the exact dispersion relation at mesh scale m."""
    mesh=build_slab_mesh(m);run.b.mesh(int(mesh.nelements),{'slab':m},0.0)
    res={}
    for pol in ('TE','TM'):
        rs=[]
        for k in (-1,0,1):
            run.b.before_solve();r=solve_slab(mesh,LAM0+k*DELTA,pol);run.b.solved(r['solve_s'],{'slab':m,'pol':pol,'k':k});rs.append(r)
        if any(r['neff'] is None for r in rs):
            res[pol]={'found':False,'pass':False};continue
        fem=[r['neff'].real for r in rs];ex=[r['exact'] for r in rs]
        ng_f=group_index(*fem,LAM0,DELTA);ng_e=group_index(*ex,LAM0,DELTA)
        te=rs[1]['te_fraction']
        ok=abs(fem[1]-ex[1])<=G1['neff'] and abs(ng_f-ng_e)<=G1['ng'] and (te>=G1['te'] if pol=='TE' else te<=G1['tm'])
        res[pol]={'found':True,'neff_fem':sig(fem[1],9),'neff_exact':sig(ex[1],9),'dneff':sig(fem[1]-ex[1]),
                  'ng_fem':sig(ng_f,9),'ng_exact':sig(ng_e,9),'dng':sig(ng_f-ng_e),'te_fraction':sig(te),
                  'n_elements':rs[1]['n_elements'],'pass':bool(ok)}
    return {'m':m,**res,'pass':all(res[p]['pass'] for p in res)}


def pub(r):
    """JSON view of one ng record."""
    keys=('m','d','delta_nm','found','n_elements','nonfinite','ok_2b')
    out={k:r[k] for k in keys}
    out['te_fractions']=[sig(t) for t in r['te_fractions']]
    if r['found']:out.update(neff=sig(r['neff'],9),ng=sig(r['ng'],9),loss_dB_per_cm=sig(r['loss_dB_per_cm']),te_spread=sig(r['te_spread']))
    return out


def mesh_change(a,b):
    return {'dneff':sig(abs(a['neff'].real-b['neff'].real)),'dng':sig(abs(a['ng']-b['ng'])),
            'pass':bool(abs(a['neff'].real-b['neff'].real)<=CONV['neff'] and abs(a['ng']-b['ng'])<=CONV['ng'])}


def run():
    caps=json.loads((OUT/'BUDGET.json').read_text());b=Budget(caps);R=Runner(b)
    dist=metadata.distribution('femwell')
    res={'experiment':'O01','layer':'O-L2 (optical cross-section eigenmode; not part of the RF L0-L3 stack)',
         'scope':'TE0 group index of the Liu 2022 cross-section under the declared branches; compared only with ng ~2.25 (p.857)',
         'femwell':{'version':dist.version,'source':json.loads(dist.read_text('direct_url.json') or '{}')},
         'nominal':{k:v for k,v in vars(N).items()},'half_domain_pec_symmetry':True,'lam0_um':LAM0,'factors':{k:[v[0],v[1]] for k,v in FACTORS.items()},
         'budget':caps}
    status='COMPLETE'
    try:
        # G1: slab gate at m = 1, else 0.5
        gates=[slab_gate(R,M_START)]
        while not gates[-1]['pass'] and gates[-1]['m']/2>=0.5-1e-12:gates.append(slab_gate(R,gates[-1]['m']/2))
        res['G1_slab']=gates
        if not gates[-1]['pass']:
            status='G1_FAIL';raise StopIteration
        m=gates[-1]['m'];rounds=[]
        norail=config(rail_t=None)
        # G2 mesh ladder on N and N without rails
        while True:
            pair={name:(R.ng(*c,m),R.ng(*c,m/2)) for name,c in (('N',config()),('N_no_rails',norail))}
            ch={k:mesh_change(*v) if all(x['ok_2b'] for x in v) else {'pass':False,'reason':'2b'} for k,v in pair.items()}
            rounds.append({'m':m,'fine':m/2,'change':ch,'records':{k:[pub(x) for x in v] for k,v in pair.items()}})
            if all(c['pass'] for c in ch.values()):break
            if m/2<=M_MIN+1e-12:
                res['G2_mesh']={'rounds':rounds,'accepted':None};status='NUMERICAL_FAIL_MESH';raise StopIteration
            m/=2
        res['G2_mesh']={'rounds':rounds,'accepted':m}
        n0=R.ng(*config(),m)
        dom=R.ng(*config(),m,d=1.5)
        steps=[R.ng(*config(),m,delta=dl) for dl in DELTAS_CHECK]
        d_dom=abs(dom['ng']-n0['ng']) if dom['ok_2b'] else float('inf')
        ngs=[n0['ng']]+[s['ng'] for s in steps if s['ok_2b']]
        d_step=max(ngs)-min(ngs) if len(ngs)==3 else float('inf')
        res['G2_domain']={'record':pub(dom),'dng':sig(d_dom),'pass':bool(d_dom<=DOMAIN_NG)}
        res['G2_step']={'records':[pub(s) for s in steps],'max_dng':sig(d_step),'pass':bool(d_step<=STEP_NG)}
        if not (res['G2_domain']['pass'] and res['G2_step']['pass']):
            status='NUMERICAL_FAIL_DOMAIN_OR_STEP';raise StopIteration
        # one-factor configurations and the MgO sensitivity branch
        oat=[];excluded=[]
        for f,(nom,ends) in FACTORS.items():
            for v in ends:
                r=R.ng(*config(**{f:v}),m)
                row={'factor':f,'value':v,**pub(r)}
                if r['ok_2b']:row['dng_vs_N']=sig(r['ng']-n0['ng'])
                else:excluded.append({'factor':f,'value':v,'reason':'no TE0' if not r['found'] else 'nonfinite or TE-fraction spread'})
                oat.append(row)
        mgo=R.ng(*config(ln='MgO'),m)
        res['N']=pub(n0);res['one_factor']=oat;res['excluded_no_te0']=excluded
        res['MgO_sensitivity']={**pub(mgo),'dng_vs_N':sig(mgo['ng']-n0['ng']) if mgo['ok_2b'] else None,'in_interval':False}
        # corners: per factor choose among {nominal, endpoints} the value giving min / max ng
        def ng_of(f,v):
            if v==FACTORS[f][0]:return n0['ng']
            row=next(x for x in oat if x['factor']==f and x['value']==v)
            return row['ng'] if row['ok_2b'] else None
        corners={}
        for side,pick in (('lo',min),('hi',max)):
            choice={}
            for f,(nom,ends) in FACTORS.items():
                cand=[(ng_of(f,v),v) for v in [nom]+ends if ng_of(f,v) is not None]
                choice[f]=pick(cand,key=lambda t:t[0])[1]
            c=config(**{k:v for k,v in choice.items() if v!=FACTORS[k][0]})
            coarse,fine=R.ng(*c,m),R.ng(*c,m/2)
            additive=sum(ng_of(f,v)-n0['ng'] for f,v in choice.items())
            corners[side]={'choice':choice,'record':pub(coarse),'fine':pub(fine),
                           'mesh_check':mesh_change(coarse,fine) if coarse['ok_2b'] and fine['ok_2b'] else {'pass':False,'reason':'2b'},
                           'sum_of_one_factor_changes':sig(additive),
                           'solved_change':sig(coarse['ng']-n0['ng']) if coarse['ok_2b'] else None}
        res['corners']=corners
        if not all(c['mesh_check']['pass'] for c in corners.values()):
            status='NUMERICAL_FAIL_CORNER';raise StopIteration
        # interval, numerical uncertainty, verdict (criteria 3, 4)
        vals=[n0['ng']]+[x['ng'] for x in oat if x['ok_2b']]+[c['record']['ng'] for c in corners.values()]
        lo,hi=min(vals),max(vals)
        mesh_u=max([float(ch['dng']) for ch in rounds[-1]['change'].values()]+[float(c['mesh_check']['dng']) for c in corners.values()])
        u=mesh_u+d_dom+d_step
        clo,chi=lo-u,hi+u
        if chi>=ROUND[0] and clo<=ROUND[1]:verdict='compatible'
        elif chi<WIDE[0] or clo>WIDE[1]:verdict='incompatible'
        else:verdict='undetermined'
        res['interval']={'ng_lo':sig(lo,9),'ng_hi':sig(hi,9),'u':sig(u),'u_parts':{'mesh':sig(mesh_u),'domain':sig(d_dom),'step':sig(d_step)},
                         'compared':[sig(clo,9),sig(chi,9)],'n_values':len(vals),'target':2.25,'rounding_interval':list(ROUND),
                         'wide_interval':list(WIDE),'verdict':verdict,
                         'N_within_rounding_interval':bool(n0['ng']+u>=ROUND[0] and n0['ng']-u<=ROUND[1])}
        res['bulk_reference_1550nm']={'ne':sig(index(LN_FILES['CLN'][0],LAM0).real),'no':sig(index(LN_FILES['CLN'][1],LAM0).real),
                                      'silica':sig(index(SILICA,LAM0).real)}
        plot(R,res,m)
    except StopIteration:
        pass
    except BudgetStop as e:
        status=f'BUDGET_STOP: {e}'
    res['status']=status;res['solves']=b.solves
    res['nonfinite_total']=sum(not np.isfinite(r['neff']) for r in R.cache.values() if r['found'])
    (OUT/'results.json').write_text(json.dumps(res,indent=1,default=str)+'\n')
    (OUT/'run_log.json').write_text(json.dumps({'elapsed_s':round(b.elapsed(),1),'peak_rss_gb':round(b.rss_gb(),3),'events':b.log},indent=1,default=str)+'\n')
    return res


INK,INK2,GRID='#0b0b0b','#52514e','#e4e3df'
SERIES,OTHER='#2a78d6','#eb6834'


def plot(R,res,m):
    import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
    from matplotlib.tri import Triangulation
    plt.rcParams.update({'font.size':9,'text.color':INK,'axes.labelcolor':INK,'xtick.color':INK2,'ytick.color':INK2,
                         'axes.edgecolor':INK2,'axes.linewidth':0.6})
    fig,(a,bx)=plt.subplots(1,2,figsize=(12,4.6),gridspec_kw={'width_ratios':[1,1.15]})
    sec,lam=config()
    mode=R.cache[(sec,round(lam,9),m,1.0)]['mode']
    basis2,I=mode.calculate_intensity()
    mesh=basis2.mesh;face=I[basis2.element_dofs].mean(axis=0);face=face/face.max()
    for sgn in (1,-1):   # half domain, mirrored for display
        a.tripcolor(Triangulation(sgn*mesh.p[0],mesh.p[1],mesh.t.T),facecolors=face,cmap='Blues',vmin=0,vmax=1,shading='flat',rasterized=True)
    for name,shp in geometry(sec).items():   # physical interfaces only (slab_near is a mesh region)
        if name not in ('ridge','clad','rail_l','rail_r'):continue
        for poly in getattr(shp,'geoms',[shp]):
            x,y=poly.exterior.xy;a.plot(x,y,color=INK2,lw=0.6)
    for yy in (0,-sec.t_box):a.axhline(yy,color=INK2,lw=0.6)
    a.set_xlim(-3,3);a.set_ylim(-1.2,1.6);a.set_aspect('equal')
    a.set_xlabel('x (um, crystal Z)');a.set_ylabel('y (um, crystal X)')
    a.set_title(f"(a) N: TE0 Poynting intensity (half domain, mirrored)\nn_eff = {res['N']['neff'][0]:.4f}, TE fraction {res['N']['te_fractions'][1]:.3f}",loc='left',fontsize=9)
    a.text(2.3,0.95,'T-rail',color=INK2,fontsize=8,ha='center');a.text(-2.3,0.95,'T-rail',color=INK2,fontsize=8,ha='center')
    rows=[('N (nominal)',res['N']['ng'])]
    for x in res['one_factor']:
        if x['ok_2b']:rows.append((f"{x['factor']} = {x['value']}",x['ng']))
    for side in ('lo','hi'):rows.append((f'corner {side}',res['corners'][side]['record']['ng']))
    u=res['interval']['u'];y=np.arange(len(rows))[::-1]
    bx.axvspan(*WIDE,color='#f1f0ec',lw=0);bx.axvspan(*ROUND,color='#d9d8d2',lw=0)
    bx.errorbar([r[1] for r in rows],y,xerr=u,fmt='o',ms=6,color=SERIES,ecolor=SERIES,elinewidth=1,capsize=2)
    if res['MgO_sensitivity'].get('ng') is not None:
        bx.plot(res['MgO_sensitivity']['ng'],-1,'D',mfc='none',mec=OTHER,ms=6)
        rows_lbl=[r[0] for r in rows]+['MgO (not in interval)'];yy=list(y)+[-1]
    else:
        rows_lbl=[r[0] for r in rows];yy=list(y)
    bx.set_yticks(yy);bx.set_yticklabels(rows_lbl,fontsize=8)
    bx.grid(axis='x',color=GRID,lw=0.6);bx.set_axisbelow(True)
    for s in ('top','right'):bx.spines[s].set_visible(False)
    bx.set_xlabel('group index ng (TE0); dark band 2.245-2.255 (rounding of 2.25), light band 2.20-2.30')
    bx.set_title(f"(b) ng per configuration; bars = numerical u ({u:.1e})\nverdict: {res['interval']['verdict']}",loc='left',fontsize=9)
    fig.tight_layout();fig.savefig(OUT/'o01_ng.png',dpi=150);plt.close(fig)


if __name__=='__main__':
    r=run()
    print(json.dumps({'status':r['status'],'interval':r.get('interval'),'N':r.get('N',{}).get('ng')},indent=1))
    sys.exit(0 if r['status']=='COMPLETE' else 2)
