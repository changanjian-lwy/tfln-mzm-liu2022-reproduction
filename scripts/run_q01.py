"""Q01: external inductance of the Liu 2022 main electrodes and the Fig.2(c) Z0 (SC-06, Track Q). See Q01 BOUNDARY.

Quasi-static, perfect conductors, vacuum: L = 1/(c^2 C0) of the signal and grounds only; Z0_pred = c (L_ext + L_int)
/ nm with nm = 2.2500 +- 0.0005 (A06 star) and 0 <= L_int <= 2 Z0 alpha / omega (A03 Fig.2(a)(c) readings).
G1 conformal-mapping gates at m = 2 -> G2 mesh ladder (N and both ground-width endpoints) and domain check (N) ->
G1 gates repeated at the accepted m -> ground-width scan -> verdict. Outer edges are natural (BOUNDARY section 5).
Rounding happens only in results.json. BUDGET.json caps are enforced.
"""
from pathlib import Path
import json,resource,sys,time
import numpy as np
import pandas as pd
from tfln_mzm.rf_quasistatic import (C_LIGHT,build_mesh,c_cpw_exact,half_space_shapes,inductance,main_cpw_shapes,
                                     solve_capacitance)

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'experiments/track_Q_quasistatic_rf/Q01_main_electrode_inductance_z0'
CURVES=ROOT/'data/digitized/liu2022/curves.csv'
S,G,T_AU=80.0,20.0,4.0
WG_N,WG_ENDS=150.0,(50.0,500.0)
WG_SCAN=tuple(float(f'{x:.6g}') for x in np.geomspace(20.0,1000.0,10))
NM,NM_U=2.2500,0.0005
F_REF=100e9
M_START,M_MIN=2.0,0.25
G1_TOL,CONV,DOMAIN=2e-3,1e-3,1e-3
FAR=50.0                                   # far-field size cap (um) before scaling by m


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


def sig(x,n=6):return None if x is None else float(f'{x:.{n}g}')


def rel(a,b):return abs(a-b)/abs(b)


def domain_x(wg,d):return max(2000.0,4*(S/2+G+wg))*d


class Runner:
    def __init__(s,b):s.b=b;s.cache={}

    def _solve(s,label,shapes,conductors,m,X,eps,regions,lines_fn):
        s.b.before_solve()
        t=time.perf_counter();mesh,owned=build_mesh(shapes,conductors,m,X,far=FAR);s.b.mesh(mesh.t.shape[1],label,time.perf_counter()-t)
        out=[]
        for e in eps:
            s.b.before_solve();t=time.perf_counter()
            c,_,_=solve_capacitance(mesh,e,regions,lines_fn(owned))
            s.b.solved(time.perf_counter()-t,label);out.append(2*c)
        return out,mesh.t.shape[1]

    def lext(s,wg,m,d=1.0):
        key=('L',wg,m,d)
        if key not in s.cache:
            X=domain_x(wg,d)
            (c0,),n=s._solve(f'L wg={wg} m={m} d={d}',main_cpw_shapes(S,G,wg,T_AU,X),['sig','gnd'],m,X,[{'air':(1.0,1.0)}],
                             [('sig',1.0),('gnd',0.0)],lambda o:[])
            if not np.isfinite(c0) or c0<=0:raise ValueError('non-finite capacitance')
            s.cache[key]={'C0_pF_per_m':c0*1e12,'L_nH_per_m':inductance(c0)*1e9,'elements':n}
        return s.cache[key]

    def gates(s,m,d=1.0):
        X=domain_x(WG_N,d)
        cases=[('vacuum',{'above':(1.0,1.0),'below':(1.0,1.0)},c_cpw_exact(S,G,WG_N)),
               ('isotropic_4.5',{'above':(1.0,1.0),'below':(4.5,4.5)},c_cpw_exact(S,G,WG_N,4.5)),
               ('anisotropic_27.9_44',{'above':(1.0,1.0),'below':(27.9,44.0)},c_cpw_exact(S,G,WG_N,np.sqrt(27.9*44.0)))]
        lines=lambda o:[('sig',1.0)]+[(k,1.0) for k in o['sig']]+[('gnd',0.0)]+[(k,0.0) for k in o['gnd']]
        vals,n=s._solve(f'gates m={m}',half_space_shapes(S,G,WG_N,X),['sig','gnd'],m,X,[c[1] for c in cases],[],lines)
        rows=[{'case':name,'C_exact_pF_per_m':sig(ex*1e12,9),'C_fem_pF_per_m':sig(v*1e12,9),'rel_err':sig(v/ex-1)}
              for (name,_,ex),v in zip(cases,vals)]
        return {'m':m,'d':d,'elements':n,'cases':rows,'pass':bool(all(abs(r['rel_err'])<=G1_TOL for r in rows))}


def readings():
    d=pd.read_csv(CURVES)
    z=d[d.panel=='fig2c'];a=d[(d.panel=='fig2a')&(d.series=='this_work')]
    near=z[(z.frequency_GHz>=90)&(z.frequency_GHz<=110)];wide=z[(z.frequency_GHz>=10)&(z.frequency_GHz<=200)]
    an=a[(a.frequency_GHz>=95)&(a.frequency_GHz<=105)]
    b_near=[float((near.value-near.uncertainty_value).min()),float((near.value+near.uncertainty_value).max())]
    b_wide=[float((wide.value-wide.uncertainty_value).min()),float((wide.value+wide.uncertainty_value).max())]
    alpha=float((an.value+an.uncertainty_value).max())*1000/(20/np.log(10))       # dB/mm -> Np/m
    lint=2*b_near[1]*alpha/(2*np.pi*F_REF)
    return {'B_near_ohm':b_near,'B_wide_ohm':b_wide,'n_near':len(near),'n_wide':len(wide),'alpha_max_Np_per_m':alpha,
            'n_alpha':len(an),'L_int_max_nH_per_m':lint*1e9}


def z0_interval(l_nh,lint_nh,u):
    lo=C_LIGHT*l_nh*1e-9/(NM+NM_U)*(1-u);hi=C_LIGHT*(l_nh+lint_nh)*1e-9/(NM-NM_U)*(1+u)
    return lo,hi


def verdict(lo,hi,near,wide):
    if hi>=near[0] and lo<=near[1]:return 'compatible'
    if hi<wide[0] or lo>wide[1]:return 'incompatible'
    return 'undetermined'


def run():
    caps=json.loads((OUT/'BUDGET.json').read_text());b=Budget(caps);R=Runner(b)
    res={'experiment':'Q01','layer':'Q-L2 (quasi-static RF cross-section, perfect conductors, 2D, high-frequency limit)',
         'scope':'external inductance of the main CPW (s 80, g 20, t 4 um) and Z0(100 GHz) = c L / nm vs Fig.2(c); no T-rails',
         'geometry_um':{'s':S,'g':G,'t_au':T_AU,'wg_nominal':WG_N,'wg_endpoints':list(WG_ENDS),'wg_scan':list(WG_SCAN)},
         'nm':[NM,NM_U],'outer_boundary':'natural (BOUNDARY section 5)','budget':caps}
    status='COMPLETE'
    try:
        rd=readings();res['readings']={k:(sig(v) if isinstance(v,float) else v) for k,v in rd.items()}
        g=R.gates(M_START);res['G1_start']=g
        if not g['pass']:status='G1_FAIL';raise StopIteration
        m=M_START;rounds=[]
        while True:
            ch={wg:rel(R.lext(wg,m)['L_nH_per_m'],R.lext(wg,m/2)['L_nH_per_m']) for wg in (WG_N,)+WG_ENDS}
            ok=max(ch.values())<=CONV
            rounds.append({'m':m,'fine':m/2,'change':{str(k):sig(v) for k,v in ch.items()},'pass':bool(ok)})
            if ok:break
            if m/2<=M_MIN+1e-12:res['G2_mesh']={'rounds':rounds,'accepted':None};status='NUMERICAL_FAIL_MESH';raise StopIteration
            m/=2
        u_mesh=max(ch.values());res['G2_mesh']={'rounds':rounds,'accepted':m}
        dom=rel(R.lext(WG_N,m)['L_nH_per_m'],R.lext(WG_N,m,2.0)['L_nH_per_m'])
        res['G2_domain']={'change':sig(dom),'pass':bool(dom<=DOMAIN)}
        if not res['G2_domain']['pass']:status='NUMERICAL_FAIL_DOMAIN';raise StopIteration
        g=R.gates(m) if m!=M_START else g;res['G1_accepted']=g
        if not g['pass']:status='G1_FAIL';raise StopIteration
        u=u_mesh+dom;lint=rd['L_int_max_nH_per_m']
        cfg={}
        for wg in (WG_N,)+WG_ENDS:
            r=R.lext(wg,m);lth=inductance(c_cpw_exact(S,G,wg))*1e9;lo,hi=z0_interval(r['L_nH_per_m'],lint,u)
            cfg[wg]={'L_ext_nH_per_m':sig(r['L_nH_per_m'],9),'L_thin_exact_nH_per_m':sig(lth,9),
                     'thickness_reduction':sig(1-r['L_nH_per_m']/lth),'Z0_pred_ohm':[sig(lo,9),sig(hi,9)],'elements':r['elements']}
        res['configs']={str(k):v for k,v in cfg.items()}
        lo=min(v['Z0_pred_ohm'][0] for v in cfg.values());hi=max(v['Z0_pred_ohm'][1] for v in cfg.values())
        scan=[]
        for wg in WG_SCAN:
            r=R.lext(wg,m);zl,zh=z0_interval(r['L_nH_per_m'],lint,u)
            scan.append({'wg_um':wg,'L_ext_nH_per_m':sig(r['L_nH_per_m'],9),'Z0_pred_ohm':[sig(zl,9),sig(zh,9)],
                         'reaches_B_near':bool(zh>=rd['B_near_ohm'][0] and zl<=rd['B_near_ohm'][1])})
        res['wg_scan']=scan
        ok_wg=[s_['wg_um'] for s_ in scan if s_['reaches_B_near']]
        res['verdict']={'Z0_pred_ohm':[sig(lo,9),sig(hi,9)],'u_rel':sig(u),'u_parts':{'mesh':sig(u_mesh),'domain':sig(dom)},
                        'B_near_ohm':[sig(x) for x in rd['B_near_ohm']],'B_wide_ohm':[sig(x) for x in rd['B_wide_ohm']],
                        'verdict':verdict(lo,hi,rd['B_near_ohm'],rd['B_wide_ohm']),
                        'L_imp_nH_per_m':[sig(NM*z/C_LIGHT*1e9) for z in rd['B_near_ohm']],
                        'wg_scanned_reaching_B_near_um':ok_wg}
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
    fig,(a,bx)=plt.subplots(1,2,figsize=(12.0,4.4))
    sc=res['wg_scan'];w=[r['wg_um'] for r in sc]
    a.plot(w,[r['L_ext_nH_per_m'] for r in sc],'o-',ms=4,color=SERIES,label='L_ext, 4 um gold (FEM)')
    a.plot(w,[inductance(c_cpw_exact(S,G,x))*1e9 for x in w],'--',color=INK2,lw=1,label='zero thickness (exact)')
    li=res['verdict']['L_imp_nH_per_m'];a.axhspan(*li,color='#d9d8d2',lw=0,label='nm Z0 / c from Fig.2(c) 90-110 GHz')
    a.set_xscale('log');a.set_xlabel('ground width w_g (um)');a.set_ylabel('inductance (nH/m)')
    a.legend(fontsize=8,frameon=False);a.grid(color=GRID,lw=0.6);a.set_axisbelow(True)
    a.set_title('(a) main-electrode external inductance (s 80, g 20 um)',loc='left',fontsize=9)
    v=res['verdict']
    bx.axhspan(*v['B_wide_ohm'],color='#f1f0ec',lw=0);bx.axhspan(*v['B_near_ohm'],color='#d9d8d2',lw=0)
    for r in sc:bx.plot([r['wg_um']]*2,r['Z0_pred_ohm'],color=SERIES,lw=3,solid_capstyle='butt')
    bx.set_xscale('log');bx.set_xlabel('ground width w_g (um)');bx.set_ylabel('Z0(100 GHz) = c L / nm (ohm)')
    bx.grid(color=GRID,lw=0.6);bx.set_axisbelow(True)
    bx.set_title(f"(b) predicted Z0 (bars: L_int 0..max, nm +-0.0005, +-u)\nbands: Fig.2(c) 90-110 GHz (dark), 10-200 GHz (light); verdict: {v['verdict']}",
                 loc='left',fontsize=9)
    for ax in (a,bx):
        for s_ in ('top','right'):ax.spines[s_].set_visible(False)
    fig.tight_layout();fig.savefig(OUT/'q01_inductance_z0.png',dpi=150);plt.close(fig)


if __name__=='__main__':
    r=run()
    print(json.dumps({'status':r['status'],'verdict':r.get('verdict'),'configs':r.get('configs')},indent=1))
    sys.exit(0 if r['status']=='COMPLETE' else 2)
