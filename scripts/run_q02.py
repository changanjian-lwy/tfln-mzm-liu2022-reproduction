"""Q02: loaded-line capacitance vs BCB thickness and the Fig.2(b) nm curve (SC-06, Track Q). See Q02 BOUNDARY.

C_avg(t) = (D - f_n) C_H + f_n C_N + (1 - D) C_U from three 2D sections (through T-heads, through necks, unloaded);
r(t) = sqrt(C_avg(t)/C_avg(1.5)) is independent of L, so Delta(t) = 2.25 (r(t) - 1) is compared with the Fig.2(b)
readings at t = 0.5 and 3.5 um; nm(1.5) = c sqrt(L C_avg(1.5)) with L in [L_ext, L_ext + L_int,max] (Q01 method).
G1 (Q01 gates, parent L_ext) -> G2 mesh ladder and domain check on N -> N and 16 one-factor configurations ->
corners per quantity -> G2 mesh check of the corners -> verdicts. Meshes are cached per geometry; permittivity and
bottom-boundary branches reuse them. Outer edges natural (Q01 BOUNDARY section 5). BUDGET.json caps are enforced.
"""
from collections import OrderedDict
from pathlib import Path
import json,sys,time
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parent))
from run_q01 import Budget,BudgetStop,Runner as Q01Runner,readings,rel,sig        # noqa: E402
from tfln_mzm.rf_quasistatic import (C_LIGHT,RFSection,build_mesh,inductance,liu_eps,liu_shapes,loaded_c,  # noqa: E402
                                     solve_capacitance)

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'experiments/track_Q_quasistatic_rf/Q02_loaded_line_bcb_nm'
Q01_RESULTS=ROOT/'experiments/track_Q_quasistatic_rf/Q01_main_electrode_inductance_z0/results.json'
A06=ROOT/'data/digitized/liu2022/a06_curves.csv'
T_POINTS=(0.5,1.5,3.5)
T_CURVE=(0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0)
DUTY,NM_STAR,NM_STAR_U=0.9,2.2500,0.0005
M_START,M_MIN=2.0,0.25
PARENT_TOL,CONV_C,CONV_R=1e-3,2e-3,2e-4
NM_BAND,NM_WIDE=(2.2495,2.2505),(2.20,2.30)
FACTORS=OrderedDict([('w_t',(2.0,[1.0,3.0])),('h_t',(0.5,[0.2,0.8])),('fn',(0.05,[0.0,0.1])),('window',('aligned',['min'])),
                     ('etch',('full',['gap'])),('eps_bcb',(2.65,[2.50])),('ln_scale',(1.0,[0.95])),
                     ('eps_sio2',(3.9,[3.8,4.5])),('substrate',('quartz',['fused'])),('bottom',('air',['pec'])),
                     ('wg',(150.0,[50.0,500.0]))])
FAR=50.0


def nominal():return {k:v[0] for k,v in FACTORS.items()}


def cfg(**ch):
    c=nominal();c.update(ch);return c


def key(c):return tuple(c[k] for k in FACTORS)


def label(c):return {k:v for k,v in c.items() if v!=FACTORS[k][0]}


def section(c,t,kind):
    win=None if c['window']=='aligned' else 1.5+c['w_t']+1.0
    return RFSection(t_bcb=t,kind=kind,w_t=c['w_t'],h_t=c['h_t'],window=win,etch=c['etch'],wg=c['wg'])


def domain_x(wg,d):return max(2000.0,4*(40.0+20.0+wg))*d


class Runner:
    def __init__(s,b):s.b=b;s.meshes={};s.caps={};s.L={}

    def mesh(s,sec,m,d):
        k=(sec,m,d)
        if k not in s.meshes:
            X=domain_x(sec.wg,d);sh,cond,thin=liu_shapes(sec,X)
            s.b.before_solve();t=time.perf_counter()
            mesh,_=build_mesh(sh,cond,m,X,thin=thin,far=FAR)
            s.b.mesh(mesh.t.shape[1],f'{sec.kind} t={sec.t_bcb} m={m} d={d}',time.perf_counter()-t)
            bot=mesh.facets_satisfying(lambda p:np.abs(p[1]+sec.t_box+sec.t_sub)<1e-9)
            s.meshes[k]=(mesh,cond,bot)
        return s.meshes[k]

    def cap(s,c,t,kind,m,d=1.0,vacuum=False):
        sec=section(c,t,kind)
        eps_key=None if vacuum else (c['eps_bcb'],c['ln_scale'],c['eps_sio2'],c['substrate'])
        k=(sec,m,d,eps_key,c['bottom'],vacuum)
        if k not in s.caps:
            mesh,cond,bot=s.mesh(sec,m,d)
            pot={'sig':1.0,'gnd':0.0,'head_s':1.0,'head_g':0.0}
            eps=liu_eps(vacuum=True) if vacuum else liu_eps(c['eps_bcb'],c['ln_scale'],c['eps_sio2'],c['substrate'])
            s.b.before_solve();t0=time.perf_counter()
            val,_,_=solve_capacitance(mesh,eps,[(n,pot[n]) for n in cond],
                                      extra_zero_facets=bot if c['bottom']=='pec' else None)
            s.b.solved(time.perf_counter()-t0,f'{kind} t={t} m={m} d={d} {"vac" if vacuum else eps_key} {c["bottom"]}')
            val*=2
            if not np.isfinite(val) or val<=0:raise ValueError('non-finite capacitance')
            s.caps[k]=val
        return s.caps[k]

    def lext(s,c,m,d=1.0):
        """External inductance (H/m): vacuum, unloaded section at t = 1.5 (position of the electrodes is irrelevant)."""
        return inductance(s.cap(c,1.5,'U',m,d,vacuum=True))

    def evaluate(s,c,m,d=1.0,ts=T_POINTS):
        rec={'C':{}}
        for t in ts:
            cu=s.cap(c,t,'U',m,d);ch=s.cap(c,t,'H',m,d);cn=s.cap(c,t,'N',m,d) if c['fn']>0 else ch
            rec['C'][t]={'U':cu,'H':ch,'N':cn,'avg':loaded_c(cu,ch,cn,DUTY,c['fn'])}
        ca=lambda t:rec['C'][t]['avg']
        rec['r']={t:np.sqrt(ca(t)/ca(1.5)) for t in ts}
        rec['delta']={t:NM_STAR*(rec['r'][t]-1) for t in ts}
        L=s.lext(c,m,d);rec['L_ext']=L
        rec['nm15']=C_LIGHT*np.sqrt(L*ca(1.5))
        rec['z0_qs']=np.sqrt(L/ca(1.5))
        rec['C_needed']=NM_STAR**2/(C_LIGHT**2*L)
        return rec


def pub(c,r):
    out={'config':label(c),'L_ext_nH_per_m':sig(r['L_ext']*1e9,9),'nm_1.5_at_L_ext':sig(r['nm15'],9),
         'Z0_qs_ohm_at_L_ext':sig(r['z0_qs'],9),'C_avg_needed_pF_per_m':sig(r['C_needed']*1e12,9)}
    out['C_pF_per_m']={str(t):{k:sig(v*1e12,9) for k,v in d.items()} for t,d in r['C'].items()}
    out['r']={str(t):sig(v,9) for t,v in r['r'].items()};out['delta']={str(t):sig(v,9) for t,v in r['delta'].items()}
    return out


def fig2b_readings():
    b=pd.read_csv(A06);nb=b[b.series=='nm_at_100GHz'];out={}
    for t in (0.5,3.5):
        w=nb[(nb.x_value-t).abs()<=0.05]
        v=float(w.value.median());u=float(w.uncertainty_value.median())
        out[t]={'n':len(w),'nm_read':v,'u_read':u,'delta':v-NM_STAR,'u':float(np.hypot(u,NM_STAR_U))}
    return out,nb


def delta_verdict(lo,hi,d,u):
    if hi>=d-u and lo<=d+u:return 'compatible'
    a,b=sorted((d/2,2*d))
    if hi<a or lo>b:return 'incompatible'
    return 'undetermined'


def run():
    caps=json.loads((OUT/'BUDGET.json').read_text());b=Budget(caps);R=Runner(b)
    res={'experiment':'Q02','layer':'Q-L2 (quasi-static RF cross-section, perfect conductors, 2D sections averaged)',
         'scope':'BCB sensitivity of the loaded line (t = 0.5, 3.5 vs 1.5 um) vs Fig.2(b), and |nm| at t = 1.5 um',
         'duty':DUTY,'factors':{k:[v[0],v[1]] for k,v in FACTORS.items()},'outer_boundary':'natural (Q01 BOUNDARY section 5)',
         'budget':caps}
    status='COMPLETE'
    try:
        rd=readings();fr,_=fig2b_readings()
        res['readings']={'L_int_max_nH_per_m':sig(rd['L_int_max_nH_per_m']),
                         'fig2b':{str(t):{k:sig(v) for k,v in d.items()} for t,d in fr.items()}}
        g=Q01Runner(b).gates(M_START)
        q01=json.loads(Q01_RESULTS.read_text());lq=q01['configs']['150.0']['L_ext_nH_per_m']
        lp=R.lext(nominal(),M_START)*1e9
        res['G1']={'gates':g,'parent_L_ext_nH_per_m':sig(lp,9),'Q01_L_ext_nH_per_m':lq,'parent_rel':sig(rel(lp,lq)),
                   'pass':bool(g['pass'] and rel(lp,lq)<=PARENT_TOL)}
        if not res['G1']['pass']:status='G1_FAIL';raise StopIteration
        def mesh_change(c,m):
            a,f=R.evaluate(c,m),R.evaluate(c,m/2)
            dc=max(rel(a['C'][t][k],f['C'][t][k]) for t in T_POINTS for k in ('U','H','N','avg'))
            dr=max(abs(a['r'][t]-f['r'][t]) for t in (0.5,3.5))
            dn=rel(a['nm15'],f['nm15'])
            return a,f,{'C':dc,'r':dr,'nm':dn,'pass':bool(dc<=CONV_C and dr<=CONV_R)}
        m=M_START;rounds=[]
        while True:
            a,f,ch=mesh_change(nominal(),m)
            rounds.append({'m':m,'fine':m/2,'change':{k:(sig(v) if k!='pass' else v) for k,v in ch.items()}})
            if ch['pass']:break
            if m/2<=M_MIN+1e-12:res['G2_mesh']={'rounds':rounds,'accepted':None};status='NUMERICAL_FAIL_MESH';raise StopIteration
            m/=2
        res['G2_mesh']={'rounds':rounds,'accepted':m};mesh_n=ch
        n0=R.evaluate(nominal(),m);dm=R.evaluate(nominal(),m,2.0)
        dc=max(rel(n0['C'][t][k],dm['C'][t][k]) for t in T_POINTS for k in ('U','H','N','avg'))
        dr=max(abs(n0['r'][t]-dm['r'][t]) for t in (0.5,3.5));dn=rel(n0['nm15'],dm['nm15'])
        res['G2_domain']={'C':sig(dc),'r':sig(dr),'nm':sig(dn),'pass':bool(dc<=CONV_C and dr<=CONV_R)}
        if not res['G2_domain']['pass']:status='NUMERICAL_FAIL_DOMAIN';raise StopIteration
        oat=OrderedDict()
        for fac,(nom,ends) in FACTORS.items():
            for v in ends:oat[(fac,v)]=(cfg(**{fac:v}),R.evaluate(cfg(**{fac:v}),m))
        res['N']=pub(nominal(),n0);res['one_factor']=[{'factor':k[0],'value':k[1],**pub(*cr)} for k,cr in oat.items()]
        metric={'delta_0.5':lambda r:r['delta'][0.5],'delta_3.5':lambda r:r['delta'][3.5],'nm_1.5':lambda r:r['nm15']}
        def value(fac,v,q):
            return metric[q](n0) if v==FACTORS[fac][0] else metric[q](oat[(fac,v)][1])
        corners=OrderedDict()
        for q in metric:
            for side,pick in (('lo',min),('hi',max)):
                choice={fac:pick([nom]+ends,key=lambda v,fac=fac:value(fac,v,q)) for fac,(nom,ends) in FACTORS.items()}
                c=cfg(**choice);ck=key(c)
                if ck not in {key(x['cfg']) for x in corners.values()}:
                    a,f,mc=mesh_change(c,m);corners[f'{q}_{side}']={'cfg':c,'rec':a,'fine':f,'mesh':mc}
                else:
                    corners[f'{q}_{side}']={'same_as':[k for k,x in corners.items() if 'cfg' in x and key(x['cfg'])==ck][0]}
        res['corners']={k:({'same_as':x['same_as']} if 'same_as' in x else
                           {'config':label(x['cfg']),'record':pub(x['cfg'],x['rec']),
                            'mesh_check':{kk:(sig(vv) if kk!='pass' else vv) for kk,vv in x['mesh'].items()}})
                        for k,x in corners.items()}
        if not all(x['mesh']['pass'] for x in corners.values() if 'mesh' in x):status='NUMERICAL_FAIL_CORNER';raise StopIteration
        recs=[n0]+[r for _,r in oat.values()]+[x['rec'] for x in corners.values() if 'rec' in x]
        u_r=max([mesh_n['r']]+[x['mesh']['r'] for x in corners.values() if 'mesh' in x])+dr
        u_nm=max([mesh_n['nm']]+[x['mesh']['nm'] for x in corners.values() if 'mesh' in x])+dn
        verdicts={}
        for t in (0.5,3.5):
            vals=[r['delta'][t] for r in recs];lo,hi=min(vals)-NM_STAR*u_r,max(vals)+NM_STAR*u_r
            d,u=fr[t]['delta'],fr[t]['u']
            verdicts[f'delta_{t}']={'interval':[sig(lo,6),sig(hi,6)],'read':[sig(d,4),sig(u,2)],
                                    'factor2_band':sorted([sig(d/2,4),sig(2*d,4)]),'verdict':delta_verdict(lo,hi,d,u)}
        lint=rd['L_int_max_nH_per_m']*1e-9
        nlo=min(C_LIGHT*np.sqrt(r['L_ext']*r['C'][1.5]['avg']) for r in recs)*(1-u_nm)
        nhi=max(C_LIGHT*np.sqrt((r['L_ext']+lint)*r['C'][1.5]['avg']) for r in recs)*(1+u_nm)
        v='compatible' if (nhi>=NM_BAND[0] and nlo<=NM_BAND[1]) else ('incompatible' if (nhi<NM_WIDE[0] or nlo>NM_WIDE[1]) else 'undetermined')
        verdicts['nm_1.5']={'interval':[sig(nlo,6),sig(nhi,6)],'band':list(NM_BAND),'wide':list(NM_WIDE),'verdict':v}
        res['verdicts']=verdicts;res['u']={'r':sig(u_r),'nm_rel':sig(u_nm),'mesh_N':{k:(sig(vv) if k!='pass' else vv) for k,vv in mesh_n.items()},
                                           'domain':{'C':sig(dc),'r':sig(dr),'nm':sig(dn)}}
        curve=R.evaluate(nominal(),m,ts=T_CURVE)
        res['N_curve']={str(t):{'r':sig(curve['r'][t],9),'nm_anchored':sig(NM_STAR*curve['r'][t],9),
                                'nm_at_L_ext':sig(C_LIGHT*np.sqrt(curve['L_ext']*curve['C'][t]['avg']),9)} for t in T_CURVE}
        plot(res)
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


def plot(res):
    import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':9,'text.color':INK,'axes.labelcolor':INK,'xtick.color':INK2,'ytick.color':INK2,
                         'axes.edgecolor':INK2,'axes.linewidth':0.6})
    fig,(a,bx)=plt.subplots(1,2,figsize=(12.5,4.8),gridspec_kw={'width_ratios':[1,1.15]})
    _,nb=fig2b_readings()
    a.plot(nb.x_value,nb.value,'.',ms=2,color=INK2,label='Fig.2(b) readings (nm at 100 GHz)')
    a.errorbar([1.5],[NM_STAR],yerr=[NM_STAR_U],fmt='*',ms=9,color=INK,label='star 2.2500 (anchor)')
    ts=[float(t) for t in res['N_curve']]
    a.plot(ts,[res['N_curve'][str(t)]['nm_anchored'] for t in ts],'o-',ms=4,color=SERIES,label='N: 2.25 r(t) (L-free)')
    for t in (0.5,3.5):
        lo,hi=res['verdicts'][f'delta_{t}']['interval'];a.plot([t,t],[NM_STAR+lo,NM_STAR+hi],color=SERIES,lw=5,alpha=0.35,solid_capstyle='butt')
    a.set_xlabel('BCB thickness t (um)');a.set_ylabel('microwave index nm');a.grid(color=GRID,lw=0.6);a.set_axisbelow(True)
    a.legend(fontsize=8,frameon=False)
    v=res['verdicts']
    a.set_title(f"(a) BCB sensitivity, anchored at t = 1.5 um (bars: envelope)\nDelta(0.5): {v['delta_0.5']['verdict']}, Delta(3.5): {v['delta_3.5']['verdict']}",loc='left',fontsize=9)
    rows=[('N (nominal)',res['N'])]+[(f"{x['factor']} = {x['value']}",x) for x in res['one_factor']]
    rows+=[(f'corner {k}',x['record']) for k,x in res['corners'].items() if 'record' in x]
    y=np.arange(len(rows))[::-1]
    bx.axvspan(*NM_WIDE,color='#f1f0ec',lw=0);bx.axvline(NM_STAR,color='#b9b8b2',lw=2)
    bx.plot([r[1]['nm_1.5_at_L_ext'] for r in rows],y,'o',ms=5,color=SERIES,label='nm(1.5) at L_ext (no internal inductance)')
    bx.set_yticks(y);bx.set_yticklabels([r[0] for r in rows],fontsize=7)
    bx.grid(axis='x',color=GRID,lw=0.6);bx.set_axisbelow(True);bx.legend(loc='lower right',fontsize=8,frameon=False)
    bx.set_xlabel('loaded-line nm at t = 1.5 um\nline: 2.25 (Fig.2(b) star); light band 2.20-2.30')
    bx.set_title(f"(b) absolute nm per configuration (duty 0.9)\nverdict: {v['nm_1.5']['verdict']} (interval {v['nm_1.5']['interval'][0]:.3f}-{v['nm_1.5']['interval'][1]:.3f})",loc='left',fontsize=9)
    for ax in (a,bx):
        for s_ in ('top','right'):ax.spines[s_].set_visible(False)
    fig.tight_layout();fig.savefig(OUT/'q02_bcb_nm.png',dpi=150);plt.close(fig)


if __name__=='__main__':
    r=run()
    print(json.dumps({'status':r['status'],'verdicts':r.get('verdicts'),'u':r.get('u')},indent=1))
    sys.exit(0 if r['status']=='COMPLETE' else 2)
