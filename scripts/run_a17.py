"""A17: femwell method gate on the Tuncer 1994 CPW benchmark (SC-01 stage 1). See A17 BOUNDARY.

Convergence loop on (mesh scale m, domain scale d) at three frequencies, then one solve per reference
frequency at the accepted (m, d). BUDGET.json caps are enforced; on a cap the partial record is written
and the script exits 2. results.json holds 6-significant-digit values so ARPACK's random start vector and
femwell's single-precision H-field assembly cannot change replayed bytes; timings go to run_log.json.
"""
from pathlib import Path
import json,resource,sys,time
from importlib import metadata
import numpy as np
import pandas as pd
from scipy.special import ellipk
from tfln_mzm.cpw_fem import CPW,build_mesh,solve

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'experiments/track_A_reproduction/A17_femwell_cpw_method_gate'
REF=ROOT/'data/external/tuncer1994_via_femwell/reference_points.csv'
BENCH=CPW(w_sig=7.0,gap=10.0,w_gnd=100.0,t_metal=0.8,sigma=6e7,eps_sub=13.0)  # femwell tutorial values
CONV={'nm':0.005,'alpha':0.02,'absZ0':0.01}          # criterion 2
M_MIN,D_MAX=0.25,8
NM_TOL,NM_ABS,AL_TOL,AL_ABS=0.05,0.10,0.10,0.10      # criterion 3
MIN_PASS,Z0_TOL,CONSIST=12,0.10,0.05                 # criteria 4 and 2b


class BudgetStop(Exception):pass


class Budget:
    def __init__(s,caps):s.caps=caps;s.t0=time.perf_counter();s.solves=0;s.log=[]
    def elapsed(s):return time.perf_counter()-s.t0
    def rss_gb(s):return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e9  # bytes on macOS
    def before_solve(s):
        if s.elapsed()>s.caps['wall_seconds']:raise BudgetStop(f"wall {s.elapsed():.0f} s > {s.caps['wall_seconds']}")
        if s.solves>=s.caps['max_solves']:raise BudgetStop(f"solves {s.solves} reached cap")
    def mesh(s,n,m,d,sec):
        s.log.append({'mesh':[m,d],'elements':n,'seconds':round(sec,2)})
        if n>s.caps['max_mesh_elements']:raise BudgetStop(f'mesh (m={m}, d={d}) has {n} elements')
    def solved(s,r,m,d):
        s.solves+=1;s.log.append({'solve':[m,d,r['f_Hz']],'seconds':round(r['solve_s'],2),'rss_gb':round(s.rss_gb(),3)})
        if r['solve_s']>s.caps['per_solve_seconds']:raise BudgetStop(f"solve took {r['solve_s']:.0f} s")
        if s.rss_gb()>s.caps['max_rss_gb']:raise BudgetStop(f'peak memory {s.rss_gb():.1f} GB')


def sig(x,n=6):
    if isinstance(x,complex):return [sig(x.real,n),sig(x.imag,n)]
    return float(f'{x:.{n}g}')


def z0_conformal(c):
    """Zero-thickness CPW on an infinitely thick substrate with infinite grounds (cheaper check (a))."""
    k=c.w_sig/(c.w_sig+2*c.gap);ee=(c.eps_sub+1)/2
    return 30*np.pi/np.sqrt(ee)*ellipk(1-k*k)/ellipk(k*k)


class Runner:
    def __init__(s,budget):s.b=budget;s.meshes={};s.cache={}
    def at(s,m,d,f):
        key=(m,d,f)
        if key not in s.cache:
            if (m,d) not in s.meshes:
                t=time.perf_counter();mesh=build_mesh(BENCH,m,d);s.b.mesh(int(mesh.nelements),m,d,time.perf_counter()-t);s.meshes[(m,d)]=mesh
            s.b.before_solve();r=solve(BENCH,s.meshes[(m,d)],f);s.cache[key]=r;s.b.solved(r,m,d)
        return s.cache[key]


def checks(r):
    """Criterion 2b for one solved point."""
    vals=[r['neff'],r['Z0_pi'],r['Z0_rlgc'],r['gamma_fem'],r['gamma_rlgc'],r['R'],r['L'],r['G'],r['C']]
    nonfinite=sum(not np.isfinite(complex(v)) for v in vals)
    g=abs(r['gamma_rlgc']/r['gamma_fem']-1);z=abs(r['Z0_rlgc']/r['Z0_pi']-1)
    ok=nonfinite==0 and r['nm']>1 and r['alpha_Np_per_m']>0 and r['Z0_pi'].real>0 and g<=CONSIST and z<=CONSIST
    return {'nonfinite':nonfinite,'gamma_consistency':sig(g),'z0_consistency':sig(z),'pass':bool(ok)}


def rel(a,b):return abs(a-b)/abs(b)


def point(r):
    return {'f_GHz':sig(r['f_Hz']/1e9),'neff':sig(r['neff']),'alpha_dB_per_cm':sig(r['alpha_dB_per_cm']),
            'Z0_pi':sig(r['Z0_pi']),'Z0_rlgc':sig(r['Z0_rlgc']),'R_ohm_per_m':sig(r['R']),'L_H_per_m':sig(r['L']),
            'G_S_per_m':sig(r['G']),'C_F_per_m':sig(r['C']),'contour_to_conduction_current':sig(r['contour_to_conduction_current']),
            'n_elements':r['n_elements'],'checks':checks(r)}


def convergence(run,freqs):
    m,d=1.0,1;rounds=[]
    while True:
        base={f:run.at(m,d,f) for f in freqs};fine={f:run.at(m/2,d,f) for f in freqs};big={f:run.at(m,2*d,f) for f in freqs}
        def change(other):
            return {f'{f/1e9:.4g}GHz':{'nm':sig(rel(base[f]['nm'],other[f]['nm'])),'alpha':sig(rel(base[f]['alpha_Np_per_m'],other[f]['alpha_Np_per_m'])),
                    'absZ0':sig(rel(abs(base[f]['Z0_pi']),abs(other[f]['Z0_pi'])))} for f in freqs}
        cm,cd=change(fine),change(big)
        ok=lambda c:all(v[q]<=CONV[q] for v in c.values() for q in CONV)
        rounds.append({'m':m,'d':d,'mesh_change':cm,'domain_change':cd,'mesh_ok':ok(cm),'domain_ok':ok(cd)})
        if ok(cm) and ok(cd):return (m,d),rounds
        if not ok(cm):
            m/=2
            if m/2<M_MIN-1e-12:return None,rounds
        if not ok(cd):
            d*=2
            if 2*d>D_MAX:return None,rounds


def run():
    caps=json.loads((OUT/'BUDGET.json').read_text());b=Budget(caps);runner=Runner(b)
    ref=pd.read_csv(REF);ref['f_Hz']=10**ref.log10_frequency_GHz*1e9
    freqs=[ref.f_Hz.min(),1e9,ref.f_Hz.max()]
    dist=metadata.distribution('femwell')
    result={'experiment':'A17','layer':'L2','scope':'femwell method gate on the Tuncer 1994 CPW benchmark; not Liu 2022 data',
            'femwell':{'version':dist.version,'source':json.loads(dist.read_text('direct_url.json') or '{}')},
            'structure':{'w_sig_um':BENCH.w_sig,'gap_um':BENCH.gap,'w_gnd_um':BENCH.w_gnd,'t_metal_um':BENCH.t_metal,
                         'sigma_S_per_m':BENCH.sigma,'eps_sub':BENCH.eps_sub},
            'convergence_frequencies_GHz':[sig(f/1e9) for f in freqs],'budget':caps}
    status='COMPLETE'
    try:
        accepted,rounds=convergence(runner,freqs)
        result['convergence']={'rounds':rounds,'accepted':accepted}
        if accepted is None:
            status='NUMERICAL_FAIL_NOT_CONVERGED'
        else:
            m,d=accepted;rows=[]
            for q in ('microwave_index','attenuation'):
                for _,p in ref[ref.quantity==q].iterrows():
                    r=runner.at(m,d,p.f_Hz)
                    model=r['nm'] if q=='microwave_index' else r['alpha_dB_per_cm']
                    tol,ab=(NM_TOL,NM_ABS) if q=='microwave_index' else (AL_TOL,AL_ABS)
                    dv=model-p.value
                    rows.append({'quantity':q,'point':int(p.point),'f_GHz':sig(p.f_Hz/1e9),'reference':float(p.value),
                                 'model':sig(model),'rel_dev':sig(dv/p.value),'within':bool(abs(dv)<=max(tol*abs(p.value),ab)),
                                 'solve':point(r)})
            result['comparison']=rows
            hi=runner.at(m,d,freqs[-1]);z0a=z0_conformal(BENCH)
            n_nm=sum(x['within'] for x in rows if x['quantity']=='microwave_index')
            n_al=sum(x['within'] for x in rows if x['quantity']=='attenuation')
            numfail=sum(not x['solve']['checks']['pass'] for x in rows)
            zdev=rel(hi['Z0_pi'].real,z0a)
            result['criteria']={'2b_numerical_failures_at_comparison_points':numfail,
                '3_nm_within':n_nm,'3_alpha_within':n_al,
                '4_Z0_real_at_fmax_ohm':sig(hi['Z0_pi'].real),'4_Z0_conformal_ohm':sig(z0a),'4_Z0_rel_dev':sig(zdev),
                '4_pass':bool(numfail==0 and n_nm>=MIN_PASS and n_al>=MIN_PASS and zdev<=Z0_TOL)}
            if numfail:status='NUMERICAL_FAIL_CHECKS'
    except BudgetStop as e:
        status=f'BUDGET_STOP: {e}'
    result['status']=status;result['solves']=b.solves
    result['all_solves']=[{'m':k[0],'d':k[1],**point(v)} for k,v in sorted(runner.cache.items())]
    (OUT/'results.json').write_text(json.dumps(result,indent=1)+'\n')
    (OUT/'run_log.json').write_text(json.dumps({'elapsed_s':round(b.elapsed(),1),'peak_rss_gb':round(b.rss_gb(),3),'events':b.log},indent=1)+'\n')
    if 'comparison' in result:plot(result)
    return result


def plot(r):
    import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
    fig,ax=plt.subplots(1,3,figsize=(13,3.8))
    for a,q,lab in ((ax[0],'microwave_index','microwave index nm'),(ax[1],'attenuation','attenuation (dB/cm)')):
        rows=[x for x in r['comparison'] if x['quantity']==q]
        f=[x['f_GHz'] for x in rows]
        a.plot(f,[x['reference'] for x in rows],'k^--',ms=5,label='Tuncer 1994 (via femwell tutorial)')
        a.plot(f,[x['model'] for x in rows],'o-',color='C0',ms=4,label='this project (femwell)')
        for x in rows:
            if not x['within']:a.plot(x['f_GHz'],x['model'],'x',color='C3',ms=9)
        a.set_xscale('log');a.set_xlabel('frequency (GHz)');a.set_ylabel(lab);a.grid(alpha=.3)
    ax[0].legend(fontsize=8)
    pts=sorted((x['solve']['f_GHz'],x['solve']['Z0_pi']) for x in r['comparison'])
    ax[2].plot([p[0] for p in pts],[p[1][0] for p in pts],'o-',ms=3,label='Re Z0')
    ax[2].plot([p[0] for p in pts],[p[1][1] for p in pts],'s-',ms=3,label='Im Z0')
    ax[2].axhline(r['criteria']['4_Z0_conformal_ohm'],color='k',ls=':',label='conformal, zero thickness')
    ax[2].set_xscale('log');ax[2].set_xlabel('frequency (GHz)');ax[2].set_ylabel('Z0 (ohm)');ax[2].legend(fontsize=8);ax[2].grid(alpha=.3)
    m,d=r['convergence']['accepted']
    fig.suptitle(f'A17 femwell CPW method gate (m={m}, d={d}); red x = outside tolerance',fontsize=10)
    fig.tight_layout();fig.savefig(OUT/'tuncer_benchmark.png',dpi=110);plt.close(fig)


if __name__=='__main__':
    r=run();print('status',r['status'],'solves',r['solves'])
    for x in r.get('convergence',{}).get('rounds',[]):print('round m',x['m'],'d',x['d'],'mesh_ok',x['mesh_ok'],'domain_ok',x['domain_ok'])
    if 'criteria' in r:print(json.dumps(r['criteria'],indent=1))
    sys.exit(0 if r['status']=='COMPLETE' else 2)
