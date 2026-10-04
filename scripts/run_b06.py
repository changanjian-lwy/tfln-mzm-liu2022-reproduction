"""B06: (ZL, L) design space with a ripple cap and the data-coverage length limit; grid truth vs GP BO.

See B06 BOUNDARY. Same provider, efficiency J, S11 constraint, BO settings and seeds as B05; the
-3 dB floor is replaced by a peak-to-peak ripple cap P(B) <= R (which implies min M >= -R), and
L <= 14.9 mm keeps the first quarter-wave feature at or above the first supported Z0 point.
"""
from pathlib import Path
import json,sys,time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parent))
from run_b01 import inputs
from run_b05 import line,efficiency_db,label,trace,ZL,BANDS,BO_BANDS,SEEDS,N_INIT,N_EVAL,TOL,C2_MAX,PAPER,INK,INK2,SERIES
from tfln_mzm.surrogate import GP,choose,log_pof

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'experiments/track_B_extensions/B06_flat_covered_design_bo'
LMM=np.round(np.arange(1,14.9+1e-9,.1),10);RIPPLES=[1.0,1.5,2.0,3.0];BO_R=2.0
J=efficiency_db(ZL[None,:],LMM[:,None])


def g1_parent(d):
    f,a,z0=inputs(d,np.arange(.025,200+1e-9,.025),False)
    b01=pd.read_csv(ROOT/'experiments/track_B_extensions/B01_termination_sweep/sweep.csv')
    b01=b01[b01.ZL_ohm!='open'];z=b01.ZL_ohm.astype(float).to_numpy()
    g=grid_rows(f,a,z0,z,6.0,[int(np.searchsorted(f,50.,side='right'))])
    dp=float(abs(g['P'][0]-b01.P50_dB.to_numpy()).max());ds=float(abs(g['c2'][0]-b01.maxS11_0_50GHz_dB.to_numpy()).max())
    return {'B01_P50_max_diff_dB':dp,'B01_maxS11_max_diff_dB':ds,'B01_loads':len(z),'pass':bool(dp<1e-9 and ds<1e-9 and len(z)==61)}


def grid_rows(f,a,z0,zl,lmm,ks):
    m,s=line(f,a,z0,zl[:,None],lmm);mag=np.abs(s)
    with np.errstate(divide='ignore'):sd=20*np.log10(mag)
    lo=np.minimum.accumulate(m,1);hi=np.maximum.accumulate(m,1);cx=np.maximum.accumulate(sd,1)
    return {'P':[np.maximum(0,hi[:,k-1])-np.minimum(0,lo[:,k-1]) for k in ks],'c2':[cx[:,k-1] for k in ks],
            'nonfinite':int((~np.isfinite(m)).sum()+(~np.isfinite(mag)).sum()),'zeros':int((mag==0).sum()),'smax':float(mag.max())}


def grid(d,step):
    f,a,z0=inputs(d,np.arange(step,200+1e-9,step),False)
    ks=[len(f) if b is None else int(np.searchsorted(f,b,side='right')) for b in BANDS]
    P=np.empty((len(BANDS),len(LMM),len(ZL)));c2=np.empty_like(P);bad=zeros=0;smax=0.0
    for i,l in enumerate(LMM):
        g=grid_rows(f,a,z0,ZL,l,ks);bad+=g['nonfinite'];zeros+=g['zeros'];smax=max(smax,g['smax'])
        for b in range(len(BANDS)):P[b,i]=g['P'][b];c2[b,i]=g['c2'][b]
    return {'f':f,'ks':ks,'P':P,'c2':c2,'nonfinite':bad,'S11_exact_zeros':zeros,'max_abs_S11':smax}


def feasible(P,c2,r):
    return (P<=r)&(c2<=C2_MAX)


def optimum(P,c2,r,mask=None):
    ok=feasible(P,c2,r)
    if mask is not None:ok&=mask
    if not ok.any():return None
    i,k=np.unravel_index(np.argmax(np.where(ok,J,-np.inf)),J.shape)
    return {'J_dB':float(J[i,k]),'ZL_ohm':float(ZL[k]),'L_mm':float(LMM[i]),'P_dB':float(P[i,k]),'c2_dB':float(c2[i,k]),
            'feasible_fraction':float(ok.sum()/(ok.size if mask is None else mask.sum())),
            'L_at_14.9mm':bool(i==len(LMM)-1),'ZL_at_bound':bool(k in (0,len(ZL)-1))}


def bo_run(seed,jf,pf,c2f,ok,xs):
    perm=np.random.default_rng(seed).choice(len(jf),N_EVAL,replace=False)
    ev=list(perm[:N_INIT]);done=np.zeros(len(jf),bool);done[ev]=True
    while len(ev)<N_EVAL:
        g1=GP().fit(xs[ev],pf[ev]);g2=GP().fit(xs[ev],c2f[ev])
        m1,s1=g1.predict(xs);m2,s2=g2.predict(xs)
        lp=log_pof(m1,s1,upper=BO_R)+log_pof(m2,s2,upper=C2_MAX)
        fe=[i for i in ev if ok[i]]
        nxt=choose(jf,lp,done,float(jf[fe].max()) if fe else None)
        if nxt is None:break
        ev.append(nxt);done[nxt]=True
    return np.array(ev),perm,{'P_lengthscales':g1.ls.tolist(),'c2_lengthscales':g2.ls.tolist()}


def run():
    t0=time.time()
    d=pd.read_csv(ROOT/'data/digitized/liu2022/curves.csv',dtype={'series':str})
    par=g1_parent(d)
    if not par['pass']:
        raise SystemExit('G1 failed: parent numbers not reproduced; no new numbers reported. '+json.dumps(par,default=float))
    fine=grid(d,.025);coarse=grid(d,.05);print(f'grids done {time.time()-t0:.0f} s',flush=True)
    sub=np.zeros(J.shape,bool);sub[::2,::2]=True
    rows2=[];rows3=[]
    for b,band in enumerate(BANDS):
        for r in RIPPLES:
            of=optimum(fine['P'][b],fine['c2'][b],r);oc=optimum(coarse['P'][b],coarse['c2'][b],r);op=optimum(fine['P'][b],fine['c2'][b],r,sub)
            dj=lambda o:None if (of is None)!=(o is None) else (0.0 if of is None else abs(of['J_dB']-o['J_dB']))
            flip=float((feasible(fine['P'][b],fine['c2'][b],r)!=feasible(coarse['P'][b],coarse['c2'][b],r)).mean())
            rows2.append({'B_GHz':label(band),'R_dB':r,'dJ_freq_dB':dj(oc),'label_flip_fraction':flip,'empty':of is None})
            rows3.append({'B_GHz':label(band),'R_dB':r,'dJ_param_dB':dj(op)})
    c4={'nonfinite':fine['nonfinite']+coarse['nonfinite'],'max_abs_S11':max(fine['max_abs_S11'],coarse['max_abs_S11']),
        'S11_exact_zero_points':fine['S11_exact_zeros']+coarse['S11_exact_zeros'],'supported_points':{'0.025':len(fine['f']),'0.05':len(coarse['f'])},
        'B_full_GHz':float(fine['f'][-1])}
    c4['pass']=bool(c4['nonfinite']==0 and c4['max_abs_S11']<=1)
    criteria={'1_G1_parent':par,
              '2_freq_grid':{'rows':rows2,'pass':all(x['dJ_freq_dB'] is not None and x['dJ_freq_dB']<=.05 and x['label_flip_fraction']<=.005 for x in rows2)},
              '3_param_grid':{'rows':rows3,'pass':all(x['dJ_param_dB'] is not None and x['dJ_param_dB']<=.1 for x in rows3)},'4_finite_passive':c4}
    base={'track':'B','parent':'B05 (A04 L1 provider)','authorization':'D9 / SC-04','changed':'ripple cap R replaces the -3 dB floor; L limited to 1-14.9 mm'}
    if not (criteria['2_freq_grid']['pass'] and criteria['3_param_grid']['pass'] and c4['pass']):
        (OUT/'results.json').write_text(json.dumps({'grade':'NUMERICAL_FAIL',**base,'criteria':criteria},indent=2,default=float)+'\n')
        print('G2 failed; front and BO not reported');return None
    iz,il=int(round((PAPER[0]-10)/.5)),int(round((PAPER[1]-1)/.1))
    front=[];paper=[]
    for b,band in enumerate(BANDS):
        pp,cp=float(fine['P'][b,il,iz]),float(fine['c2'][b,il,iz])
        for r in RIPPLES:
            o=optimum(fine['P'][b],fine['c2'][b],r)
            front.append({'B_GHz':label(band),'R_dB':r,**(o or {'J_dB':None})})
            paper.append({'B_GHz':label(band),'R_dB':r,'P_dB':pp,'c2_dB':cp,'feasible':bool(pp<=r and cp<=C2_MAX),
                          'gap_to_Jstar_dB':None if o is None else o['J_dB']-float(J[il,iz])})
    pd.DataFrame(front).to_csv(OUT/'front.csv',index=False)
    xs=np.stack(np.meshgrid((LMM-1)/13.9,(ZL-10)/110,indexing='ij'),-1).reshape(-1,2)[:,::-1]
    jf=J.ravel();bo={};rows=[];curves={};seed0={}
    for band in BO_BANDS:
        b=BANDS.index(band);pf=fine['P'][b].ravel();c2f=fine['c2'][b].ravel();ok=feasible(pf,c2f,BO_R)
        o=optimum(fine['P'][b],fine['c2'][b],BO_R)
        if o is None:
            bo[f'{band:g}']={'empty_feasible_set':True};continue
        jstar=o['J_dB'];nb=[];nr=[];cb=[];cr=[]
        for s in SEEDS:
            seq,perm,hyp=bo_run(s,jf,pf,c2f,ok,xs)
            n1,k1=trace(seq,jf,ok,jstar);n2,k2=trace(perm,jf,ok,jstar);nb.append(n1);nr.append(n2);cb.append(k1);cr.append(k2)
            for meth,n,k,used in (('BO',n1,k1,len(seq)),('random',n2,k2,len(perm))):
                rows.append({'B_GHz':band,'seed':s,'method':meth,'n_to_success':n,'best_feasible_J_dB':float(k[-1]),'evaluations_used':used})
            if s==0:seed0[band]={'sequence':seq.tolist(),'hyperparameters_final':hyp}
        mb,mr=float(np.median(nb)),float(np.median(nr))
        bo[f'{band:g}']={'empty_feasible_set':False,'J_star_dB':jstar,'BO_successes':int(sum(n<=N_EVAL for n in nb)),'random_successes':int(sum(n<=N_EVAL for n in nr)),
                         'BO_median_n':mb,'random_median_n':mr,'BO_n':nb,'random_n':nr,'pass_5':bool(sum(n<=N_EVAL for n in nb)>=27),'pass_6':bool(mb<=.5*mr),
                         'seed0_hyperparameters_final':seed0[band]['hyperparameters_final']}
        curves[band]=(np.array(cb),np.array(cr),jstar)
        print(f'B={band:g} GHz BO done {time.time()-t0:.0f} s',flush=True)
    pd.DataFrame(rows).to_csv(OUT/'bo_runs.csv',index=False)
    live=[v for v in bo.values() if not v['empty_feasible_set']]
    criteria['5_BO_reliability']={'applicable_targets':len(live),'pass':all(v['pass_5'] for v in live)}
    criteria['6_BO_vs_random']={'applicable_targets':len(live),'pass':all(v['pass_6'] for v in live)}
    result={'grade':'SENSITIVITY_ONLY',**base,'criteria':criteria,'front':front,
            'paper_design':{'ZL_ohm':PAPER[0],'L_mm':PAPER[1],'J_dB':float(J[il,iz]),'by_band':paper},'bo':bo,
            'bo_settings':{'R_dB':BO_R,'n_init':N_INIT,'n_eval':N_EVAL,'seeds':len(SEEDS),'tolerance_dB':TOL,'grid_points':int(J.size)}}
    (OUT/'results.json').write_text(json.dumps(result,indent=2,default=float)+'\n')
    figures(fine,front,seed0,curves,float(J[il,iz]),paper)
    print(f'total {time.time()-t0:.0f} s')
    return result


def figures(fine,front,seed0,curves,jp,paper):
    plt.rcParams.update({'font.size':9,'axes.edgecolor':INK2,'axes.labelcolor':INK,'xtick.color':INK2,'ytick.color':INK2,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(1,3,figsize=(13,4.4),layout='constrained',sharey=True)
    ext=[ZL[0]-.25,ZL[-1]+.25,LMM[0]-.05,LMM[-1]+.05]
    for ax,band,col in zip(axs,BO_BANDS,SERIES):
        b=BANDS.index(band);ok=feasible(fine['P'][b],fine['c2'][b],BO_R)
        ax.imshow(np.where(ok,1.,np.nan),origin='lower',extent=ext,aspect='auto',cmap=matplotlib.colors.ListedColormap([col]),alpha=.22,interpolation='nearest')
        cs=ax.contour(ZL,LMM,J,levels=[-6,-3,0,3,6],colors=INK2,linewidths=.6);ax.clabel(cs,fmt='%+g dB',fontsize=7)
        if band in seed0:
            seq=np.array(seed0[band]['sequence']);zl=ZL[seq%len(ZL)];lm=LMM[seq//len(ZL)]
            ax.scatter(zl[:N_INIT],lm[:N_INIT],s=26,facecolors='none',edgecolors=INK2,lw=1,label='BO seed 0: random start (6)')
            ax.scatter(zl[N_INIT:],lm[N_INIT:],s=14,color=INK,lw=0,label='BO seed 0: chosen points (34)')
        o=next(x for x in front if x['B_GHz']==label(band) and x['R_dB']==BO_R)
        if o['J_dB'] is not None:
            ax.scatter([o['ZL_ohm']],[o['L_mm']],marker='*',s=160,color=col,edgecolors=INK,lw=.8,zorder=5,label=f"grid optimum J* = {o['J_dB']:+.2f} dB")
        ax.scatter([PAPER[0]],[PAPER[1]],marker='s',s=40,facecolors='none',edgecolors=INK,lw=1.2,zorder=5,label='paper design 40 ohm, 6 mm')
        ax.set_title(f'B = {band:g} GHz (shaded: ripple <= {BO_R:g} dB and max S11 <= -10 dB)',fontsize=9,color=INK)
        ax.set_xlabel('termination ZL (ohm)');ax.legend(fontsize=7,loc='upper right',frameon=False)
    axs[0].set_ylabel('electrode length L (mm), limited to 14.9 mm by Z0 data coverage')
    fig.suptitle('B06 | flat, data-covered design space on the A04 L1 provider; contours = low-frequency efficiency J; model sensitivity, not measurement',fontsize=10,color=INK)
    fig.savefig(OUT/'b06_design_space.png',dpi=140);plt.close(fig)
    fig,(a1,a2)=plt.subplots(1,2,figsize=(12,4.4),layout='constrained')
    fb=lambda key:float(fine['f'][-1]) if key=='B_full' else float(key)
    for r,col in zip(RIPPLES,['#2a78d6','#eb6834','#1baf7a','#4a3aa7']):
        pts=[(fb(x['B_GHz']),x['J_dB'],x['L_at_14.9mm']) for x in front if x['R_dB']==r and x['J_dB'] is not None]
        if not pts:continue
        bx,jy,cen=zip(*pts)
        a1.plot(bx,jy,'-',color=col,lw=2,label=f'J*(B), ripple <= {r:g} dB')
        a1.scatter([x for x,c in zip(bx,cen) if not c],[y for y,c in zip(jy,cen) if not c],color=col,s=30,zorder=4)
        a1.scatter([x for x,c in zip(bx,cen) if c],[y for y,c in zip(jy,cen) if c],facecolors='white',edgecolors=col,s=30,zorder=4)
    pf=[fb(p['B_GHz']) for p in paper if p['R_dB']==BO_R and p['feasible']]
    a1.axhline(jp,color=INK2,ls='--',lw=1,label=f'paper design J = {jp:+.2f} dB')
    if pf:a1.scatter(pf,[jp]*len(pf),marker='s',s=36,facecolors='none',edgecolors=INK,zorder=5,label=f'paper design feasible at ripple <= {BO_R:g} dB')
    a1.set_xlabel('target bandwidth B (GHz); hollow = optimum at the 14.9 mm limit');a1.set_ylabel('efficiency J (dB re 6 mm, 50 ohm)')
    a1.grid(alpha=.25);a1.legend(fontsize=7,frameon=False);a1.set_title('efficiency-bandwidth front with a ripple cap (L1 model)',fontsize=9,color=INK)
    k=np.arange(1,N_EVAL+1)
    for (band,(cb,cr,js)),col in zip(curves.items(),SERIES):
        a2.plot(k,(cb>=js-TOL).mean(0),color=col,lw=2,label=f'BO, B = {band:g} GHz')
        a2.plot(k,(cr>=js-TOL).mean(0),color=col,lw=1.5,ls='--',label=f'random, B = {band:g} GHz')
    a2.axvline(N_INIT,color=INK2,lw=.8,ls=':')
    a2.set_xlabel('model evaluations');a2.set_ylabel(f'fraction of 30 seeds within {TOL} dB of J*');a2.set_ylim(-.02,1.02)
    a2.grid(alpha=.25);a2.legend(fontsize=7,frameon=False,ncol=2);a2.set_title(f'search efficiency, ripple <= {BO_R:g} dB, {J.size}-point grid',fontsize=9,color=INK)
    fig.suptitle('B06 | front and Bayesian-optimization efficiency; model sensitivity, not measurement',fontsize=10,color=INK)
    fig.savefig(OUT/'b06_front_bo.png',dpi=140);plt.close(fig)


if __name__=='__main__':
    r=run()
    if r:print(json.dumps({'criteria':{k:v.get('pass') for k,v in r['criteria'].items()},
                           'front':[x for x in r['front'] if x['R_dB']==BO_R],'paper':r['paper_design']['by_band'][2::4],
                           'bo':{k:{x:v.get(x) for x in ('J_star_dB','BO_successes','random_successes','BO_median_n','random_median_n')} for k,v in r['bo'].items()}},indent=1,default=float))
