"""A17b: the A17 femwell method gate with second-order elements and mode tracking. See A17b BOUNDARY.

NOT RUN: G1 failed before the run (A17b BOUNDARY section 8, low-frequency breakdown at 48 MHz). Kept as the
starting point for any later attempt, which needs its own boundary.

Criteria, structure, reference points and the convergence loop are A17's (run_a17.convergence is reused).
Changes: every solve uses order=2; before solving (m, d, f) the next-coarser level pred(m, d) is solved at
the same frequency and its complex neff seeds the shift (n_guess); BUDGET.json caps are also enforced by a
watchdog thread during a solve, which writes the partial record and exits 2. results.json holds
6-significant-digit values; solves differ between runs at ~1e-4, so A17b is registered by static hash.
"""
from pathlib import Path
import json,os,sys,threading,time
from importlib import metadata
import pandas as pd
import psutil
from run_a17 import (BENCH,REF,MIN_PASS,NM_TOL,NM_ABS,AL_TOL,AL_ABS,Z0_TOL,Budget,BudgetStop,
                     checks,convergence,point,rel,sig,z0_conformal)
from tfln_mzm.cpw_fem import build_mesh,solve

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'experiments/track_A_reproduction/A17b_femwell_cpw_order2_tracking'
ORDER=2
WATCH_S=2


def pred(m,d):
    """Level whose neff seeds the shift at (m, d): next-coarser mesh, else next-smaller domain; (1,1) has none."""
    if m<1:return (2*m,d)
    if d>1:return (m,d//2)
    return None


class Runner:
    def __init__(s,budget):s.b=budget;s.meshes={};s.cache={};s.current=None
    def at(s,m,d,f):
        key=(m,d,f)
        if key not in s.cache:
            p=pred(m,d);guess=None if p is None else s.at(*p,f)['neff']
            if (m,d) not in s.meshes:
                t=time.perf_counter();mesh=build_mesh(BENCH,m,d);s.b.mesh(int(mesh.nelements),m,d,time.perf_counter()-t);s.meshes[(m,d)]=mesh
            s.b.before_solve();s.current=(key,time.perf_counter())
            r=solve(BENCH,s.meshes[(m,d)],f,order=ORDER,n_guess=guess)
            s.current=None;r['n_guess']=guess;r['seeded_by']=p;s.cache[key]=r;s.b.solved(r,m,d)
        return s.cache[key]


class Watchdog(threading.Thread):
    """Checks the BUDGET caps every WATCH_S seconds, also while the main thread is inside a solve."""
    def __init__(s,runner,stop):super().__init__(daemon=True);s.r=runner;s.stop=stop;s.proc=psutil.Process()
    def run(s):
        caps=s.r.b.caps
        while True:
            time.sleep(WATCH_S)
            rss=s.proc.memory_info().rss/1e9;cur=s.r.current;wall=s.r.b.elapsed();reason=None
            if rss>caps['max_rss_gb']:reason=f"resident memory {rss:.1f} GB > {caps['max_rss_gb']}"
            elif wall>caps['wall_seconds']:reason=f"wall {wall:.0f} s > {caps['wall_seconds']}"
            elif cur and time.perf_counter()-cur[1]>caps['per_solve_seconds']:
                reason=f"solve {list(cur[0])} running {time.perf_counter()-cur[1]:.0f} s > {caps['per_solve_seconds']}"
            if reason and s.stop(f'BUDGET_STOP (watchdog): {reason}',cur[0] if cur else None):os._exit(2)


def solve_record(k,v):
    p=point(v)
    return {'m':k[0],'d':k[1],**p,'n_dofs':v['n_dofs'],'seeded_by':list(v['seeded_by']) if v['seeded_by'] else None,
            'n_guess':sig(v['n_guess']) if v['n_guess'] is not None else None}


def run():
    caps=json.loads((OUT/'BUDGET.json').read_text());b=Budget(caps);runner=Runner(b)
    ref=pd.read_csv(REF);ref['f_Hz']=10**ref.log10_frequency_GHz*1e9
    freqs=[ref.f_Hz.min(),1e9,ref.f_Hz.max()]
    dist=metadata.distribution('femwell')
    result={'experiment':'A17b','parent':'A17','layer':'L2',
            'scope':'femwell method gate on the Tuncer 1994 CPW benchmark; not Liu 2022 data',
            'element_order':ORDER,
            'mode_tracking':'n_guess = complex neff of pred(m,d) at the same frequency; pred=(2m,d) if m<1, '
                            'else (m,d/2) if d>1; (1,1) keeps the femwell default shift',
            'femwell':{'version':dist.version,'source':json.loads(dist.read_text('direct_url.json') or '{}')},
            'structure':{'w_sig_um':BENCH.w_sig,'gap_um':BENCH.gap,'w_gnd_um':BENCH.w_gnd,'t_metal_um':BENCH.t_metal,
                         'sigma_S_per_m':BENCH.sigma,'eps_sub':BENCH.eps_sub},
            'convergence_frequencies_GHz':[sig(f/1e9) for f in freqs],'budget':caps}
    lock=threading.Lock();written=[]

    def write(status,in_progress=None):
        """Write results.json and run_log.json once; the watchdog and the main thread race for it."""
        with lock:
            if written:return False
            written.append(status)
            out={**result,'status':status,'solves':b.solves}
            if in_progress:out['in_progress']={'m':in_progress[0],'d':in_progress[1],'f_GHz':sig(in_progress[2]/1e9)}
            out['all_solves']=[solve_record(k,v) for k,v in sorted(dict(runner.cache).items())]
            (OUT/'results.json').write_text(json.dumps(out,indent=1)+'\n')
            (OUT/'run_log.json').write_text(json.dumps({'elapsed_s':round(b.elapsed(),1),'peak_rss_gb':round(b.rss_gb(),3),
                                                        'events':list(b.log)},indent=1)+'\n')
            return True

    Watchdog(runner,write).start()
    status='COMPLETE';error=None
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
    except Exception as e:  # keep the partial record of a long run, then fail loudly
        status=f'ERROR: {type(e).__name__}: {e}';error=e
    if not write(status):  # the watchdog stopped the run first and is about to exit
        time.sleep(10*WATCH_S)
    if error:raise error
    result['status']=status;result['solves']=b.solves
    if 'comparison' in result:plot(result)
    return result


def plot(r):
    import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
    fig,ax=plt.subplots(1,3,figsize=(13,3.8))
    for a,q,lab in ((ax[0],'microwave_index','microwave index nm'),(ax[1],'attenuation','attenuation (dB/cm)')):
        rows=[x for x in r['comparison'] if x['quantity']==q]
        f=[x['f_GHz'] for x in rows]
        a.plot(f,[x['reference'] for x in rows],'k^--',ms=5,label='Tuncer 1994 (via femwell tutorial)')
        a.plot(f,[x['model'] for x in rows],'o-',color='C0',ms=4,label='this project (femwell, order 2)')
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
    fig.suptitle(f'A17b femwell CPW method gate, order 2 + mode tracking (m={m}, d={d}); red x = outside tolerance',fontsize=10)
    fig.tight_layout();fig.savefig(OUT/'tuncer_benchmark.png',dpi=110);plt.close(fig)


if __name__=='__main__':
    r=run();print('status',r['status'],'solves',r['solves'])
    for x in r.get('convergence',{}).get('rounds',[]):print('round m',x['m'],'d',x['d'],'mesh_ok',x['mesh_ok'],'domain_ok',x['domain_ok'])
    if 'criteria' in r:print(json.dumps(r['criteria'],indent=1))
    sys.exit(0 if r['status']=='COMPLETE' else 2)
