"""O02: electro-optic overlap and VpiL from the Liu 2022 cross-section (SC-05, Track O). See O02 BOUNDARY.

Two meshes per configuration (BOUNDARY section 5): the optical TE0 on O01's half-domain mesh, the potential on an
enlarged half-domain mesh with the same rules; E_x is carried over by point evaluation (eo_overlap.FieldX).
G1 two-layer capacitor gate and parent reproduction -> G2 mesh ladder (N), domain and linearity checks ->
13 one-factor configurations -> two corners by convention-A VpiL, solved directly -> G2 mesh check of the
corners -> verdicts for conventions A and B. Optical solves keep O01's quadrature; the Gamma integral uses a
6th-order rule (BOUNDARY section 5, ERRATA E-03). Corner choices
use full-precision values; rounding happens only in results.json. BUDGET.json caps are enforced.
"""
from collections import OrderedDict
from dataclasses import replace
from importlib import metadata
from pathlib import Path
import json,resource,sys,time
import numpy as np
from shapely.geometry import box
from tfln_mzm.optical_fem import Section,build_mesh,solve_te0,index,LN_FILES
from tfln_mzm.eo_overlap import (DC_LN,FieldX,dc_permittivity,deps_xx,drive_dofs,element_ex,gamma_a,solve_potential,vpil_a,vpil_b)

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'experiments/track_O_optical/O02_eo_overlap_vpil'
LAM0,GAP,V0,DUTY,INTORDER=1.55,3.0,10.0,0.9,None   # None = O01's solver quadrature
M_START,M_MIN=2.0,0.5
G1_TOL,PARENT_NEFF,PARENT_TOL=1e-6,1.8976164,1e-4              # criterion 1
CONV,DOMAIN,LIN,TE_DV=0.005,0.005,0.001,0.01                     # criteria 2, 2b
BAND,WIDE=(1.59,1.65),(1.32,1.92)                                # criteria 3, 4
N=Section()
FACTORS=OrderedDict([('lam_nm',(1550.0,[1530.0,1570.0])),('theta_deg',(75.0,[60.0,90.0])),('width_ref',('top',['bottom'])),
                     ('rail_t',(0.5,[0.2,0.8])),('rail_w',(2.0,[1.0,3.0])),('state',('clamped',['unclamped'])),
                     ('eps_sio2',(3.9,[3.8,4.5])),('substrate',('quartz',['fused']))])
EXTRA=('lam_nm','state','eps_sio2')


class BudgetStop(Exception):pass


class Budget:
    def __init__(s,caps):s.caps=caps;s.t0=time.perf_counter();s.solves=0;s.log=[]
    def elapsed(s):return time.perf_counter()-s.t0
    def rss_gb(s):return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e9
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


def extent(d):return (20*d,-2-18*d,0.6+20*d)


def config(**ch):
    """O02 configuration: Section plus wavelength (um), LN mechanical state and silica DC permittivity."""
    lam=ch.pop('lam_nm',1550.0)/1000;state=ch.pop('state','clamped');eps=ch.pop('eps_sio2',3.9)
    return (replace(N,**ch),lam,state,eps)


def label(c):
    sec,lam,state,eps=c
    out={k:v for k,v in vars(sec).items() if v!=getattr(N,k)}
    if abs(lam-LAM0)>1e-12:out['lam_nm']=round(lam*1000,3)
    if state!='clamped':out['state']=state
    if eps!=3.9:out['eps_sio2']=eps
    return out


class Runner:
    def __init__(s,b):s.b=b;s.meshes={};s.pots={};s.modes={};s.keep={}
    def mesh(s,sec,m,d=None):
        """d=None: O01's optical mesh; otherwise the enlarged electrostatic mesh of scale d."""
        k=(sec,m,d)
        if k not in s.meshes:
            t=time.perf_counter()
            mesh=build_mesh(sec,m,half=True) if d is None else build_mesh(sec,m,half=True,extent=extent(d))
            s.b.mesh(int(mesh.nelements),{'m':m,'d':d if d is not None else 'optical',**label((sec,LAM0,'clamped',3.9))},
                     time.perf_counter()-t);s.meshes[k]=mesh
        return s.meshes[k]
    def potential(s,c,m,d):
        """(FieldX at 1 V gap voltage, element-mean E_x on the optical mesh)."""
        sec,lam,state,eps=c;k=(sec,state,eps,m,d)
        if k not in s.pots:
            from skfem import Basis,ElementTriP2
            mesh=s.mesh(sec,m,d)
            basis,phi=solve_potential(mesh,dc_permittivity(sec,state,eps),drive_dofs(Basis(mesh,ElementTriP2()),1.0))
            f=FieldX(basis,phi);s.pots[k]=(f,element_ex(s.mesh(sec,m),f),basis,phi)
        return s.pots[k]
    def te0(s,c,m,d,V):
        sec,lam,state,eps=c
        k=(sec,round(lam,9),m) if V==0 else (sec,round(lam,9),state,eps,m,d,V)
        if k not in s.modes:
            mesh=s.mesh(sec,m);dep=None
            if V!=0:
                _,ex,_,_=s.potential(c,m,d)
                dep=deps_xx(mesh,ex,index(LN_FILES['CLN'][0],lam).real,DC_LN[state]['r33'],V)
            s.b.before_solve();r=solve_te0(sec,mesh,lam,deps_xx=dep,intorder=INTORDER)
            s.b.solved(r['solve_s'],{'m':m,'d':d,'V':V,**label(c)})
            mode=r.pop('mode',None)
            if V==0 and mode is not None:s.keep[k]=mode   # needed for convention A
            s.modes[k]=r
        return s.modes[k],s.keep.get((sec,round(lam,9),m))
    def evaluate(s,c,m,d=1.0,v0=V0):
        """Raw (full-precision) record for one configuration."""
        sec,lam,state,eps=c
        r0,mode=s.te0(c,m,d,0)
        rec={'cfg':c,'m':m,'d':d,'V0':v0,'found0':r0['found'],'n_elements':r0['n_elements']}
        if not r0['found']:rec['ok']=False;return rec
        field,_,_,_=s.potential(c,m,d)
        ne=index(LN_FILES['CLN'][0],lam).real;r33=DC_LN[state]['r33']
        g=gamma_a(mode,field,GAP,1.0)
        rp,_=s.te0(c,m,d,v0);rm,_=s.te0(c,m,d,-v0)
        found=rp['found'] and rm['found']
        rec.update(neff0=r0['neff'],te0=r0['te_fraction'],gamma=g,vpil_a=vpil_a(lam,GAP,ne,r33,g),ne=ne,r33=r33)
        if found:
            dn=(rp['neff'].real-rm['neff'].real)/(2*v0)
            rec.update(dn_dv=dn,vpil_b=vpil_b(lam,dn),te_pm=[rp['te_fraction'],rm['te_fraction']],
                       dte=max(abs(rp['te_fraction']-r0['te_fraction']),abs(rm['te_fraction']-r0['te_fraction'])))
        vals=[rec.get(k) for k in ('gamma','vpil_a','dn_dv','vpil_b')]
        rec['nonfinite']=sum(v is None or not np.isfinite(v) for v in vals)+sum(not np.isfinite(x['neff']) for x in (r0,rp,rm) if x['found'])
        rec['ok']=bool(found and rec['nonfinite']==0 and rec['dte']<=TE_DV)
        return rec


def pub(rec):
    out={'config':label(rec['cfg']),'m':rec['m'],'d':rec['d'],'V0':rec['V0'],'n_elements':rec['n_elements'],'ok':rec['ok']}
    for k in ('neff0','te0','gamma','vpil_a','dn_dv','vpil_b','dte','ne','r33'):
        if k in rec:out[k]=sig(rec[k],9) if k in ('neff0','gamma','vpil_a','dn_dv','vpil_b') else sig(rec[k])
    if 'te_pm' in rec:out['te_pm']=[sig(t) for t in rec['te_pm']]
    if 'vpil_a' in rec:out['vpil_a_device']=sig(rec['vpil_a']/DUTY,9)
    if 'vpil_b' in rec:out['vpil_b_device']=sig(rec['vpil_b']/DUTY,9)
    return out


def rel(a,b):return abs(a-b)/abs(b)


def two_layer_gate(axis,a=0.7,d=2.0,e1=(27.9,44.0),e2=(3.9,3.9)):
    """Criterion 1 static gate: plates normal to `axis`, layer 1 (LN clamped tensor) on [0, a], silica on [a, d]."""
    from femwell.mesh import mesh_from_OrderedDict
    from skfem import Basis,ElementTriP2
    from skfem.io.meshio import from_meshio
    sh=OrderedDict(l1=box(0,0,a,2.0),l2=box(a,0,d,2.0)) if axis==0 else OrderedDict(l1=box(0,0,2.0,a),l2=box(0,a,2.0,d))
    mesh=from_meshio(mesh_from_OrderedDict(sh,{},default_resolution_max=0.25))
    lo=mesh.facets_satisfying(lambda p:np.abs(p[axis])<1e-9,boundaries_only=True)
    hi=mesh.facets_satisfying(lambda p:np.abs(p[axis]-d)<1e-9,boundaries_only=True)
    b=Basis(mesh,ElementTriP2())
    basis,phi=solve_potential(mesh,{'l1':e1,'l2':e2},[(b.get_dofs(facets=lo).all(),0.0),(b.get_dofs(facets=hi).all(),1.0)])
    gr=basis.interpolate(phi).grad[axis]
    k1,k2=e1[axis],e2[axis];E1=1/(a+(d-a)*k1/k2);exact=[E1,E1*k1/k2]
    fem=[float(np.max(np.abs(gr[mesh.subdomains[n]]/x-1))) for n,x in zip(('l1','l2'),exact)]
    return {'axis':'x' if axis==0 else 'y','exact_V_per_um':[sig(x) for x in exact],'max_rel_err':[sig(e) for e in fem],'pass':bool(max(fem)<=G1_TOL)}


def run():
    caps=json.loads((OUT/'BUDGET.json').read_text());b=Budget(caps);R=Runner(b)
    dist=metadata.distribution('femwell')
    res={'experiment':'O02','layer':'O-L2 (optical mode + quasi-static field, 2D cross-section)',
         'scope':'Gamma and VpiL of the Liu 2022 rail cross-section under the declared branches; compared only with VpiL 1.62 V*cm (p.856)',
         'femwell':{'version':dist.version,'source':json.loads(dist.read_text('direct_url.json') or '{}')},
         'conventions':{'A':'Gamma=(g/V) int_LN Ex|Ex,opt|^2/int|Ex,opt|^2; VpiL=lam g/(2 ne^3 r33 Gamma)',
                        'B':'re-solve TE0 with deps_xx=-ne^4 r33 <Ex> at +-V0; VpiL=lam/(4|dn/dV|) (BOUNDARY section 5)'},
         'gap_um':GAP,'duty':DUTY,'V0':V0,'intorder':INTORDER,'factors':{k:[v[0],v[1]] for k,v in FACTORS.items()},
         'dc_ln':DC_LN,'budget':caps}
    status='COMPLETE'
    try:
        gates=[two_layer_gate(0),two_layer_gate(1)]
        n2=R.evaluate(config(),M_START)
        parent=abs(n2['neff0'].real-PARENT_NEFF) if n2['found0'] else float('inf')
        res['G1']={'two_layer':gates,'parent_neff':sig(n2['neff0'].real,9) if n2['found0'] else None,'parent_dneff':sig(parent),
                   'pass':bool(all(g['pass'] for g in gates) and parent<=PARENT_TOL)}
        if not res['G1']['pass']:status='G1_FAIL';raise StopIteration
        m=M_START;rounds=[]
        while True:
            a,f=R.evaluate(config(),m),R.evaluate(config(),m/2)
            ok=a['ok'] and f['ok']
            ch={'gamma':sig(rel(a['gamma'],f['gamma'])),'dn_dv':sig(rel(a['dn_dv'],f['dn_dv']))} if ok else {}
            passed=ok and rel(a['gamma'],f['gamma'])<=CONV and rel(a['dn_dv'],f['dn_dv'])<=CONV
            rounds.append({'m':m,'fine':m/2,'change':ch,'pass':bool(passed),'records':[pub(a),pub(f)]})
            if passed:break
            if m/2<=M_MIN+1e-12:
                res['G2_mesh']={'rounds':rounds,'accepted':None};status='NUMERICAL_FAIL_MESH';raise StopIteration
            m/=2
        res['G2_mesh']={'rounds':rounds,'accepted':m}
        n0=R.evaluate(config(),m);dom=R.evaluate(config(),m,d=2.0);lin=R.evaluate(config(),m,v0=2*V0)
        dg=rel(n0['gamma'],dom['gamma']) if dom['ok'] else float('inf')
        dd=rel(n0['dn_dv'],dom['dn_dv']) if dom['ok'] else float('inf')
        dl=rel(n0['dn_dv'],lin['dn_dv']) if lin['ok'] else float('inf')
        res['G2_domain']={'record':pub(dom),'gamma':sig(dg),'dn_dv':sig(dd),'pass':bool(dg<=DOMAIN and dd<=DOMAIN)}
        res['G2_linearity']={'record':pub(lin),'dn_dv':sig(dl),'pass':bool(dl<=LIN)}
        if not (res['G2_domain']['pass'] and res['G2_linearity']['pass'] and n0['ok']):
            status='NUMERICAL_FAIL_DOMAIN_OR_LINEARITY';raise StopIteration
        oat={}
        for fac,(nom,ends) in FACTORS.items():
            for v in ends:oat[(fac,v)]=R.evaluate(config(**{fac:v}),m)
        res['N']=pub(n0);res['one_factor']=[{'factor':k[0],'value':k[1],**pub(r)} for k,r in oat.items()]
        bad=[{'factor':k[0],'value':k[1]} for k,r in oat.items() if not r['ok']]
        res['excluded']=bad
        def va(fac,v):
            if v==FACTORS[fac][0]:return n0['vpil_a']
            r=oat[(fac,v)];return r['vpil_a'] if r['ok'] else None
        corners={};craw={}
        for side,pick in (('lo',min),('hi',max)):
            choice={}
            for fac,(nom,ends) in FACTORS.items():
                cand=[(va(fac,v),v) for v in [nom]+ends if va(fac,v) is not None]
                choice[fac]=pick(cand,key=lambda t:t[0])[1]
            c=config(**{k:v for k,v in choice.items() if v!=FACTORS[k][0]})
            a,f=R.evaluate(c,m),R.evaluate(c,m/2)
            mc={'gamma':sig(rel(a['gamma'],f['gamma'])),'dn_dv':sig(rel(a['dn_dv'],f['dn_dv'])),
                'pass':bool(a['ok'] and f['ok'] and rel(a['gamma'],f['gamma'])<=CONV and rel(a['dn_dv'],f['dn_dv'])<=CONV)} \
               if a['ok'] and f['ok'] else {'pass':False,'reason':'2b'}
            corners[side]={'choice':choice,'record':pub(a),'fine':pub(f),'mesh_check':mc};craw[side]=(a,f)
        res['corners']=corners
        if not all(c['mesh_check']['pass'] for c in corners.values()):status='NUMERICAL_FAIL_CORNER';raise StopIteration
        recs=[n0]+[r for r in oat.values() if r['ok']]+[craw['lo'][0],craw['hi'][0]]
        last=rounds[-1]['records']
        mesh_u={'A':max([rel(R.evaluate(config(),m)['gamma'],R.evaluate(config(),m/2)['gamma'])]+[rel(a['gamma'],f['gamma']) for a,f in craw.values()]),
                'B':max([rel(R.evaluate(config(),m)['dn_dv'],R.evaluate(config(),m/2)['dn_dv'])]+[rel(a['dn_dv'],f['dn_dv']) for a,f in craw.values()])}
        verdicts={}
        for conv,key,u in (('A','vpil_a',mesh_u['A']+dg),('B','vpil_b',mesh_u['B']+dd+dl)):
            vals=[r[key] for r in recs];lo,hi=min(vals),max(vals)
            dlo,dhi=lo/DUTY*(1-u),hi/DUTY*(1+u)
            if dhi>=BAND[0] and dlo<=BAND[1]:v='compatible'
            elif dhi<WIDE[0] or dlo>WIDE[1]:v='incompatible'
            else:v='undetermined'
            dneed=[lo/BAND[1],hi/BAND[0]]
            verdicts[conv]={'vpil_rail':[sig(lo,9),sig(hi,9)],'u_rel':sig(u),'vpil_device_compared':[sig(dlo,9),sig(dhi,9)],
                            'measured_band':list(BAND),'wide_band':list(WIDE),'verdict':v,'duty_needed':[sig(x) for x in dneed],
                            'duty_needed_contains_0.9':bool(dneed[0]<=DUTY<=dneed[1]),'duty_needed_exceeds_1':bool(dneed[0]>1),
                            'n_values':len(vals)}
        res['verdicts']=verdicts;res['u_parts']={'mesh_A':sig(mesh_u['A']),'mesh_B':sig(mesh_u['B']),'domain_gamma':sig(dg),
                                                  'domain_dn_dv':sig(dd),'linearity':sig(dl)}
        plot(R,res,m,n0)
    except StopIteration:
        pass
    except BudgetStop as e:
        status=f'BUDGET_STOP: {e}'
    res['status']=status;res['solves']=b.solves
    (OUT/'results.json').write_text(json.dumps(res,indent=1,default=str)+'\n')
    (OUT/'run_log.json').write_text(json.dumps({'elapsed_s':round(b.elapsed(),1),'peak_rss_gb':round(b.rss_gb(),3),'events':b.log},indent=1,default=str)+'\n')
    return res


INK,INK2,GRID='#0b0b0b','#52514e','#e4e3df'
SERIES,SERIES2='#2a78d6','#eb6834'


def plot(R,res,m,n0):
    import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
    from matplotlib.tri import Triangulation
    plt.rcParams.update({'font.size':9,'text.color':INK,'axes.labelcolor':INK,'xtick.color':INK2,'ytick.color':INK2,
                         'axes.edgecolor':INK2,'axes.linewidth':0.6})
    fig,(a,bx)=plt.subplots(1,2,figsize=(12.5,4.8),gridspec_kw={'width_ratios':[1,1.2]})
    c=config();field,_,basis,phi=R.potential(c,m,1.0);mesh=basis.mesh
    ex=element_ex(mesh,field)
    tri=Triangulation(mesh.p[0],mesh.p[1],mesh.t.T)
    pc=a.tripcolor(tri,facecolors=np.abs(ex),cmap='Blues',shading='flat',vmin=0,vmax=np.percentile(np.abs(ex),99.5),rasterized=True)
    fig.colorbar(pc,ax=a,label='|E_x| (V/um) at 1 V gap voltage',shrink=0.8)
    from tfln_mzm.optical_fem import geometry
    for name,shp in geometry(c[0],half=True,extent=extent(1.0)).items():
        if name not in ('ridge','clad','rail_r'):continue
        for poly in getattr(shp,'geoms',[shp]):
            x,y=poly.exterior.xy;a.plot(x,y,color=INK2,lw=0.6)
    for yy in (0,-2):a.axhline(yy,color=INK2,lw=0.6)
    a.set_xlim(0,4);a.set_ylim(-1.0,1.4);a.set_aspect('equal')
    a.set_xlabel('x (um, crystal Z); x = 0 is the symmetry plane');a.set_ylabel('y (um)')
    a.set_title(f"(a) N: quasi-static |E_x| (clamped LN), half domain\nGamma_A = {res['N']['gamma']:.4f}",loc='left',fontsize=9)
    rows=[('N (nominal)',res['N'])]+[(f"{x['factor']} = {x['value']}",x) for x in res['one_factor'] if x['ok']]
    rows+=[(f'corner {s}',res['corners'][s]['record']) for s in ('lo','hi')]
    y=np.arange(len(rows))[::-1]
    bx.axvspan(*WIDE,color='#f1f0ec',lw=0);bx.axvspan(*BAND,color='#d9d8d2',lw=0)
    bx.plot([r[1]['vpil_a_device'] for r in rows],y+0.12,'o',ms=6,color=SERIES,label='convention A (ref [10] Eq.(1)(2))')
    bx.plot([r[1]['vpil_b_device'] for r in rows],y-0.12,'s',ms=5,color=SERIES2,label='convention B (direct perturbation)')
    bx.set_yticks(y);bx.set_yticklabels([r[0] for r in rows],fontsize=8)
    bx.grid(axis='x',color=GRID,lw=0.6);bx.set_axisbelow(True)
    for s in ('top','right'):bx.spines[s].set_visible(False)
    bx.legend(loc='lower right',fontsize=8,frameon=False)
    bx.set_xlabel('VpiL of the device = rail-section VpiL / 0.9 (V*cm); dark band 1.59-1.65 (measured), light 1.32-1.92')
    v=res['verdicts'];bx.set_title(f"(b) VpiL per configuration, duty 0.9\nverdict A: {v['A']['verdict']}, B: {v['B']['verdict']}",loc='left',fontsize=9)
    fig.tight_layout();fig.savefig(OUT/'o02_vpil.png',dpi=150);plt.close(fig)


if __name__=='__main__':
    r=run()
    print(json.dumps({'status':r['status'],'verdicts':r.get('verdicts'),'N':r.get('N')},indent=1))
    sys.exit(0 if r['status']=='COMPLETE' else 2)
