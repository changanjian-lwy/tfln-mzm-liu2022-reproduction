"""B05: joint (ZL, L) design space on the frozen A04 L1 provider; grid truth vs GP Bayesian optimization.

See B05 BOUNDARY. The efficiency J is closed-form; the constraints c1 (minimum normalized EO
response up to B) and c2 (maximum S11 up to B) come from the model on supported frequencies.
A BO "evaluation" looks up the precomputed grid value, which equals calling the deterministic model.
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
from tfln_mzm.interaction import average_voltage
from tfln_mzm.microwave import propagation
from tfln_mzm.response import eo_response_db,magnitude_db,bandwidth_3db
from tfln_mzm.termination import s11
from tfln_mzm.surrogate import GP,choose,log_pof

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'experiments/track_B_extensions/B05_joint_zl_length_bo'
ZL=np.round(np.arange(10,120+1e-9,.5),10);LMM=np.round(np.arange(1,30+1e-9,.1),10)
BANDS=[50.,75.,100.,125.,150.,175.,None]  # None: last supported frequency (B_full)
BO_BANDS=[50.,100.,150.];SEEDS=range(30);N_INIT=6;N_EVAL=40;TOL=.1;C1_MIN=-3.;C2_MAX=-10.
NG=2.25;F0=1e-3;PAPER=(40.,6.)
INK,INK2='#0b0b0b','#52514e';SERIES=['#2a78d6','#eb6834','#1baf7a']


def efficiency_db(zl,lmm):
    return 20*np.log10((lmm/6)*(zl/(50+zl))/.5)


def line(f,a,z0,zl,lmm):
    """M (dB re DC divider) and complex S11 on supported f; zl is a column of loads, lmm in mm."""
    v=average_voltage(f*1e9,length_m=lmm/1000,z0=z0,zg=50,zl=zl,n_m=NG,n_g=NG,alpha_np_m=a)
    return eo_response_db(v,zl/(50+zl)),s11(propagation(f*1e9,NG,a),lmm/1000,z0,zl)


def label(b):
    return 'B_full' if b is None else f'{b:g}'


def g1_parents(d):
    f,a,z0=inputs(d,np.arange(.025,200+1e-9,.025),False);low=f<=50
    b01=pd.read_csv(ROOT/'experiments/track_B_extensions/B01_termination_sweep/sweep.csv')
    b01=b01[b01.ZL_ohm!='open'];z=b01.ZL_ohm.astype(float).to_numpy()
    m,s=line(f,a,z0,z[:,None],6.0)
    p50=np.maximum(0,m[:,low].max(1))-np.minimum(0,m[:,low].min(1));sd=magnitude_db(s)[:,low].max(1)
    dp=float(abs(p50-b01.P50_dB.to_numpy()).max());ds=float(abs(sd-b01.maxS11_0_50GHz_dB.to_numpy()).max())
    b02=json.loads((ROOT/'experiments/track_B_extensions/B02_length_scan/results.json').read_text())['results']
    dd=[]
    for r in b02:
        if r['value'] not in (9.0,12.0):continue
        m,_=line(f,a,z0,np.array([[r['ZL_ohm']]]),r['value'])
        b=bandwidth_3db(np.concatenate([[F0],f]),np.concatenate([[0.0],m[0]]))
        same=b['censored']==r['D1']['censored'] and (b['censored'] or abs(b['bandwidth_hz']-r['D1']['GHz'])<1e-9)
        dd.append({'L_mm':r['value'],'ZL_ohm':r['ZL_ohm'],'D1_GHz':b['bandwidth_hz'],'B02_D1_GHz':r['D1']['GHz'],'same':bool(same)})
    return {'B01_P50_max_diff_dB':dp,'B01_maxS11_max_diff_dB':ds,'B01_loads':len(z),'B02_D1':dd,
            'pass':bool(dp<1e-9 and ds<1e-9 and len(z)==61 and len(dd)==4 and all(x['same'] for x in dd))}


def grid(d,step):
    f,a,z0=inputs(d,np.arange(step,200+1e-9,step),False)
    ks=[len(f) if b is None else int(np.searchsorted(f,b,side='right')) for b in BANDS]
    c1=np.empty((len(BANDS),len(LMM),len(ZL)));c2=np.empty_like(c1)
    bad=zeros=0;smax=0.0;zl=ZL[:,None]
    for i,l in enumerate(LMM):
        m,s=line(f,a,z0,zl,l);mag=np.abs(s)
        bad+=int((~np.isfinite(m)).sum()+(~np.isfinite(mag)).sum());zeros+=int((mag==0).sum());smax=max(smax,float(mag.max()))
        with np.errstate(divide='ignore'):sd=20*np.log10(mag)
        cm=np.minimum.accumulate(m,1);cx=np.maximum.accumulate(sd,1)
        for b,k in enumerate(ks):c1[b,i]=cm[:,k-1];c2[b,i]=cx[:,k-1]
    return {'f':f,'ks':ks,'c1':c1,'c2':c2,'nonfinite':bad,'S11_exact_zeros':zeros,'max_abs_S11':smax}


J=efficiency_db(ZL[None,:],LMM[:,None])


def feasible(c1,c2):
    return (c1>=C1_MIN)&(c2<=C2_MAX)


def optimum(c1,c2,mask=None):
    ok=feasible(c1,c2)
    if mask is not None:ok&=mask
    if not ok.any():return None
    i,k=np.unravel_index(np.argmax(np.where(ok,J,-np.inf)),J.shape)
    return {'J_dB':float(J[i,k]),'ZL_ohm':float(ZL[k]),'L_mm':float(LMM[i]),'c1_dB':float(c1[i,k]),'c2_dB':float(c2[i,k]),
            'feasible_fraction':float(ok.sum()/(ok.size if mask is None else mask.sum())),
            'censored':{'L_at_30mm':bool(i==len(LMM)-1),'ZL_at_bound':bool(k in (0,len(ZL)-1))}}


def trace(seq,jf,ok,jstar):
    best=-np.inf;n=N_EVAL+1;curve=[]
    for k,i in enumerate(seq,1):
        if ok[i]:best=max(best,jf[i])
        if n>N_EVAL and best>=jstar-TOL:n=k
        curve.append(best)
    curve+= [best]*(N_EVAL-len(curve))
    return n,np.array(curve)


def bo_run(seed,jf,c1f,c2f,ok,xs):
    perm=np.random.default_rng(seed).choice(len(jf),N_EVAL,replace=False)
    ev=list(perm[:N_INIT]);done=np.zeros(len(jf),bool);done[ev]=True
    while len(ev)<N_EVAL:
        g1=GP().fit(xs[ev],c1f[ev]);g2=GP().fit(xs[ev],c2f[ev])
        m1,s1=g1.predict(xs);m2,s2=g2.predict(xs)
        lp=log_pof(m1,s1,lower=C1_MIN)+log_pof(m2,s2,upper=C2_MAX)
        fe=[i for i in ev if ok[i]]
        nxt=choose(jf,lp,done,float(jf[fe].max()) if fe else None)
        if nxt is None:break
        ev.append(nxt);done[nxt]=True
    return np.array(ev),perm,{'c1_lengthscales':g1.ls.tolist(),'c2_lengthscales':g2.ls.tolist()}


def run():
    t0=time.time()
    d=pd.read_csv(ROOT/'data/digitized/liu2022/curves.csv',dtype={'series':str})
    par=g1_parents(d)
    if not par['pass']:
        raise SystemExit('G1 failed: parent numbers not reproduced; no new numbers reported. '+json.dumps(par,default=float))
    fine=grid(d,.025);coarse=grid(d,.05);print(f'grids done {time.time()-t0:.0f} s',flush=True)
    sub=np.zeros(J.shape,bool);sub[::2,::2]=True
    g2=[]
    for b,band in enumerate(BANDS):
        of=optimum(fine['c1'][b],fine['c2'][b]);oc=optimum(coarse['c1'][b],coarse['c2'][b]);op=optimum(fine['c1'][b],fine['c2'][b],sub)
        flip=float((feasible(fine['c1'][b],fine['c2'][b])!=feasible(coarse['c1'][b],coarse['c2'][b])).mean())
        dj=lambda o:None if (of is None)!=(o is None) else (0.0 if of is None else abs(of['J_dB']-o['J_dB']))
        g2.append({'B_GHz':label(band),'dJ_freq_dB':dj(oc),'label_flip_fraction':flip,'dJ_param_dB':dj(op)})
    criteria={'1_G1_parents':par,
              '2_freq_grid':{'rows':[{k:r[k] for k in ('B_GHz','dJ_freq_dB','label_flip_fraction')} for r in g2],
                             'pass':all(r['dJ_freq_dB'] is not None and r['dJ_freq_dB']<=.05 and r['label_flip_fraction']<=.005 for r in g2)},
              '3_param_grid':{'rows':[{k:r[k] for k in ('B_GHz','dJ_param_dB')} for r in g2],
                              'pass':all(r['dJ_param_dB'] is not None and r['dJ_param_dB']<=.1 for r in g2)},
              '4_finite_passive':{'nonfinite':fine['nonfinite']+coarse['nonfinite'],'max_abs_S11':max(fine['max_abs_S11'],coarse['max_abs_S11']),
                                  'S11_exact_zero_points':fine['S11_exact_zeros']+coarse['S11_exact_zeros'],
                                  'supported_points':{'0.025':len(fine['f']),'0.05':len(coarse['f'])},
                                  'supported_points_le_B':{label(b):k for b,k in zip(BANDS,fine['ks'])},'B_full_GHz':float(fine['f'][-1])}}
    c4=criteria['4_finite_passive'];c4['pass']=bool(c4['nonfinite']==0 and c4['max_abs_S11']<=1)
    base={'track':'B','parent':'A04 L1 provider','authorization':'D9 / SC-04','changed':'ZL 10-120 ohm (0.5) x L 1-30 mm (0.1), jointly'}
    if not (criteria['2_freq_grid']['pass'] and criteria['3_param_grid']['pass'] and c4['pass']):
        (OUT/'results.json').write_text(json.dumps({'grade':'NUMERICAL_FAIL',**base,'criteria':criteria},indent=2,default=float)+'\n')
        print('G2 failed; front and BO not reported');return None
    iz,il=int(round((PAPER[0]-10)/.5)),int(round((PAPER[1]-1)/.1))
    front=[];paper=[]
    for b,band in enumerate(BANDS):
        o=optimum(fine['c1'][b],fine['c2'][b]);front.append({'B_GHz':label(band),**{k:v for k,v in o.items() if k!='censored'},**o['censored']})
        c1p,c2p=float(fine['c1'][b,il,iz]),float(fine['c2'][b,il,iz])
        paper.append({'B_GHz':label(band),'c1_dB':c1p,'c2_dB':c2p,'feasible':bool(c1p>=C1_MIN and c2p<=C2_MAX),'gap_to_Jstar_dB':o['J_dB']-float(J[il,iz])})
    pd.DataFrame(front).to_csv(OUT/'front.csv',index=False)
    xs=np.stack(np.meshgrid((LMM-1)/29,(ZL-10)/110,indexing='ij'),-1).reshape(-1,2)[:,::-1]
    jf=J.ravel();bo={};rows=[];curves={};seed0={}
    for band in BO_BANDS:
        b=BANDS.index(band);c1f=fine['c1'][b].ravel();c2f=fine['c2'][b].ravel();ok=feasible(c1f,c2f)
        jstar=front[b]['J_dB'];nb=[];nr=[];cb=[];cr=[]
        for s in SEEDS:
            seq,perm,hyp=bo_run(s,jf,c1f,c2f,ok,xs)
            n1,k1=trace(seq,jf,ok,jstar);n2,k2=trace(perm,jf,ok,jstar);nb.append(n1);nr.append(n2);cb.append(k1);cr.append(k2)
            for meth,n,k,used in (('BO',n1,k1,len(seq)),('random',n2,k2,len(perm))):
                rows.append({'B_GHz':band,'seed':s,'method':meth,'n_to_success':n,'best_feasible_J_dB':float(k[-1]),'evaluations_used':used})
            if s==0:seed0[band]={'sequence':seq.tolist(),'hyperparameters_final':hyp}
        mb,mr=float(np.median(nb)),float(np.median(nr))
        bo[f'{band:g}']={'J_star_dB':jstar,'BO_successes':int(sum(n<=N_EVAL for n in nb)),'random_successes':int(sum(n<=N_EVAL for n in nr)),
                         'BO_median_n':mb,'random_median_n':mr,'BO_n':nb,'random_n':nr,
                         'pass_5':bool(sum(n<=N_EVAL for n in nb)>=27),'pass_6':bool(mb<=.5*mr),'seed0_hyperparameters_final':seed0[band]['hyperparameters_final']}
        curves[band]=(np.array(cb),np.array(cr),jstar)
        print(f'B={band:g} GHz BO done {time.time()-t0:.0f} s',flush=True)
    pd.DataFrame(rows).to_csv(OUT/'bo_runs.csv',index=False)
    criteria['5_BO_reliability']={'pass':all(v['pass_5'] for v in bo.values())}
    criteria['6_BO_vs_random']={'pass':all(v['pass_6'] for v in bo.values())}
    result={'grade':'SENSITIVITY_ONLY',**base,'criteria':criteria,'front':front,
            'paper_design':{'ZL_ohm':PAPER[0],'L_mm':PAPER[1],'J_dB':float(J[il,iz]),'by_band':paper},'bo':bo,
            'bo_settings':{'n_init':N_INIT,'n_eval':N_EVAL,'seeds':len(SEEDS),'tolerance_dB':TOL,'grid_points':int(J.size)}}
    (OUT/'results.json').write_text(json.dumps(result,indent=2,default=float)+'\n')
    figures(fine,front,seed0,curves,J[il,iz],paper)
    print(f'total {time.time()-t0:.0f} s')
    return result


def figures(fine,front,seed0,curves,jp,paper):
    plt.rcParams.update({'font.size':9,'axes.edgecolor':INK2,'axes.labelcolor':INK,'xtick.color':INK2,'ytick.color':INK2,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(1,3,figsize=(13,4.4),layout='constrained',sharey=True)
    ext=[ZL[0]-.25,ZL[-1]+.25,LMM[0]-.05,LMM[-1]+.05]
    for ax,band,col in zip(axs,BO_BANDS,SERIES):
        b=BANDS.index(band);ok=feasible(fine['c1'][b],fine['c2'][b])
        ax.imshow(np.where(ok,1.,np.nan),origin='lower',extent=ext,aspect='auto',cmap=matplotlib.colors.ListedColormap([col]),alpha=.22,interpolation='nearest')
        cs=ax.contour(ZL,LMM,J,levels=[-3,0,3,6,9],colors=INK2,linewidths=.6);ax.clabel(cs,fmt='%+g dB',fontsize=7)
        seq=np.array(seed0[band]['sequence']);zl=ZL[seq%len(ZL)];lm=LMM[seq//len(ZL)]
        ax.scatter(zl[:N_INIT],lm[:N_INIT],s=26,facecolors='none',edgecolors=INK2,lw=1,label='BO seed 0: random start (6)')
        ax.scatter(zl[N_INIT:],lm[N_INIT:],s=14,color=INK,lw=0,label='BO seed 0: chosen points (34)')
        o=front[b];ax.scatter([o['ZL_ohm']],[o['L_mm']],marker='*',s=160,color=col,edgecolors=INK,lw=.8,zorder=5,label=f"grid optimum J* = {o['J_dB']:+.2f} dB")
        ax.scatter([PAPER[0]],[PAPER[1]],marker='s',s=40,facecolors='none',edgecolors=INK,lw=1.2,zorder=5,label='paper design 40 ohm, 6 mm')
        ax.set_title(f'B = {band:g} GHz (shaded: c1 >= -3 dB and max S11 <= -10 dB)',fontsize=9,color=INK)
        ax.set_xlabel('termination ZL (ohm)');ax.legend(fontsize=7,loc='upper right',frameon=False)
    axs[0].set_ylabel('electrode length L (mm)')
    fig.suptitle('B05 | feasible design space on the A04 L1 provider; contours = low-frequency efficiency J; model sensitivity, not measurement',fontsize=10,color=INK)
    fig.savefig(OUT/'b05_design_space.png',dpi=140);plt.close(fig)
    fig,(a1,a2)=plt.subplots(1,2,figsize=(12,4.4),layout='constrained')
    bx=[float(fine['f'][-1]) if r['B_GHz']=='B_full' else float(r['B_GHz']) for r in front]
    a1.plot(bx,[r['J_dB'] for r in front],'-o',color=SERIES[0],lw=2,ms=6,label='grid optimum J*(B)')
    pf=[x for x,p in zip(bx,paper) if p['feasible']]
    a1.axhline(jp,color=INK2,ls='--',lw=1,label=f'paper design J = {jp:+.2f} dB')
    if pf:a1.scatter(pf,[jp]*len(pf),marker='s',s=36,facecolors='none',edgecolors=INK,zorder=5,label='paper design feasible at this B')
    a1.set_xlabel('target bandwidth B (GHz): c1 >= -3 dB and max S11 <= -10 dB up to B');a1.set_ylabel('efficiency J (dB re 6 mm, 50 ohm)')
    a1.grid(alpha=.25);a1.legend(fontsize=8,frameon=False);a1.set_title('efficiency-bandwidth front (L1 model)',fontsize=9,color=INK)
    k=np.arange(1,N_EVAL+1)
    for (band,(cb,cr,js)),col in zip(curves.items(),SERIES):
        a2.plot(k,(cb>=js-TOL).mean(0),color=col,lw=2,label=f'BO, B = {band:g} GHz')
        a2.plot(k,(cr>=js-TOL).mean(0),color=col,lw=1.5,ls='--',label=f'random, B = {band:g} GHz')
    a2.axvline(N_INIT,color=INK2,lw=.8,ls=':');a2.text(N_INIT+.4,.02,'end of random start',fontsize=7,color=INK2)
    a2.set_xlabel('model evaluations');a2.set_ylabel(f'fraction of 30 seeds within {TOL} dB of J*');a2.set_ylim(-.02,1.02)
    a2.grid(alpha=.25);a2.legend(fontsize=7,frameon=False,ncol=2);a2.set_title(f'search efficiency on a {J.size}-point grid',fontsize=9,color=INK)
    fig.suptitle('B05 | front and Bayesian-optimization efficiency; model sensitivity, not measurement',fontsize=10,color=INK)
    fig.savefig(OUT/'b05_front_bo.png',dpi=140);plt.close(fig)


if __name__=='__main__':
    r=run()
    if r:print(json.dumps({'criteria':{k:v.get('pass') for k,v in r['criteria'].items()},'front':r['front'],
                           'bo':{k:{x:v[x] for x in ('J_star_dB','BO_successes','random_successes','BO_median_n','random_median_n')} for k,v in r['bo'].items()}},indent=1,default=float))
