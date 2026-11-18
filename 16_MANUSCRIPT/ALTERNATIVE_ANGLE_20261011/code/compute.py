"""One-pass fixed research components; all writes confined to new runs directory."""
import csv, gzip, hashlib, json, math, random, sys, time
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import numpy as np
import scipy, os
from methods import certificate, float_membership, intervals, count_tests, bh, ranks, unit_contrasts, freeze_training, target_score
ROOT=Path(__file__).resolve().parents[1];AUTH=ROOT.parents[1]
D3=AUTH/'15_FINAL_AUDIT/PHASE5_D3_RANK_TRANSFER_20261010_V1'
P54=AUTH/'15_FINAL_AUDIT/PHASE5_P54_SYNTHESIS_20261010_V1/actual'
FOLDS=[('A_to_B',['A2018','A2019'],['B']),('B_to_A',['B'],['A2018','A2019']),('A2018_to_A2019',['A2018'],['A2019']),('A2019_to_A2018',['A2019'],['A2018']),('C2017_to_C2018',['C2017'],['C2018']),('C2018_to_C2017',['C2018'],['C2017']),('B_Inbred_to_B_Hybrid',['B_Inbred'],['B_Hybrid']),('B_Hybrid_to_B_Inbred',['B_Hybrid'],['B_Inbred'])]
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def save(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x',encoding='utf-8') as f:json.dump(x,f,ensure_ascii=False,indent=2,allow_nan=False)
def sha(p):
    with open(p,'rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def npz(p,**kw):
    with Path(p).open('xb') as f:np.savez_compressed(f,**kw)
def check_sources():
    for rel, expected in read(ROOT/'evidence/PUBLIC_INPUT_HASHES.json').items():
        assert sha(AUTH/rel)==expected,rel

def component_a(out,pilot):
    rng=random.Random(26101101 if pilot else 26101102);cases=[]
    for i in range(100 if pilot else 10000):
        s1=F(rng.randint(1,20),10);s2=F(rng.randint(1,20),10);sb=F(rng.randint(1,20),10)
        a=F(rng.randint(-40,40),10);b=F(rng.randint(-40,40),10);q=rng.choice([F(1,2),F(3,5),F(7,10),F(4,5),F(9,10)])
        cases.append({'id':f'random_{i}','a':str(a),'b':str(b),'alpha':str((s1*s1+s2*s2)/4),'beta':str(s1*s2/2),'z':str(sb*sb),'q':str(q)})
    edges=[('empty',1,1,0,0,1,F(7,10)),('constant_pass',1,1,1,0,1,F(7,10)),('constant_fail',1,1,10,0,1,F(7,10)),('point_qhalf',1,1,2,1,2,F(1,2)),('zero_endpoint',1,1,1,1,1,F(7,10)),('negative',-1,-2,2,1,1,F(7,10)),('opposite',1,-1,2,1,1,F(7,10)),('zero_effect',0,1,2,1,1,F(7,10)),('huge_t',1,10,2,1,1,F(7,10)),('invalid',1,1,1,2,1,F(7,10)),('near_zero',1,1,F(1,10**12),0,F(1,10**12),F(7,10))]
    for name,*v in edges:cases.append(dict(id=name,**{k:str(x) for k,x in zip(['a','b','alpha','beta','z','q'],v)}))
    # Explicit grid-resolution example, distinguished from sampled cases.
    cases.append(dict(id='constructed_isolated_point',a='1',b='1',alpha='2',beta='1',z='20123/10000',q='1/2'))
    for c in cases:
        a,b,A,B,z,q=[F(c[k]) for k in ('a','b','alpha','beta','z','q')]
        c['result']=certificate(a,b,A,B,z,q)
        if c['result']['class'] in ('INVALID_INPUT','NOT_ASSESSABLE'):continue
        c['grid_comparison']={}
        for n in (11,21,101):
            states=[]
            for j in range(n):
                rho=F(-1)+F(2*j,n-1);x=A+B*rho
                if x<=0:continue
                tau=max(F(0),((a-b)**2-x-z)/2)
                states.append(a*b>0 and max(x+tau,z+tau)<=q*(x+z+2*tau))
            c['grid_comparison'][str(n)]={'regular_points':len(states),'any_pass':any(states),'all_regular_points_pass':all(states),'missed_some_domain':c['result']['class']=='SOME_RHO' and not any(states)}
    save(out/'certificate_cases.json',cases)
    rows=[]
    for ix in read(P54/'AB_SYNTHESIS_ALL_FAMILIES.json')['rows']:
        p=P54/ix['path'];assert sha(p)==ix['sha256'];old=read(p);x=old['input']
        for q in ([F(7,10)] if pilot else [F(3,5),F(7,10),F(4,5),F(9,10)]):
            support=True if old['classes']['G4']=='ALL_RHO' else False if old['classes']['G4']=='NONE' else None
            r=certificate(*[F(x[k]) if x[k] is not None else None for k in ('a','b','alpha','beta','z')],q=q,support=support)
            rows.append({'family':old['family'],'q':str(q),'input':x,'support':support,'source_sha256':ix['sha256'],'source_path':str(p.relative_to(AUTH)),'historical_class':old['classes']['joint'],'diagnostics':old['source_diagnostics'],'result':r})
    save(out/'family_certificates.json',rows)
    print('A certificate cases and all historical labels saved',flush=True)

def component_b(out,pilot):
    R=50 if pilot else 5000;rng=np.random.default_rng(26101103 if pilot else 26101104);rows=[]
    for i,(s2,sb,rhot,means) in enumerate(product([.2,.5],[.2,.5],[-.9,0.,.9],[(0.,0.),(.3,.3),(.3,-.3),(.3,.6)])):
        s1=.2;alpha=(s1*s1+s2*s2)/4;beta=s1*s2/2;z=sb*sb;thetaA,thetaB=means
        initial=rng.bit_generator.state;Z=rng.normal(size=(3,R))
        a1=thetaA+s1*Z[0];a2=thetaA+s2*(rhot*Z[0]+math.sqrt(1-rhot*rhot)*Z[1]);a=(a1+a2)/2;b=thetaB+sb*Z[2]
        xtrue=alpha+beta*rhot;d2=(a-b)**2;c=.4
        L=np.maximum(0,np.minimum(z*.3/.7,z-c*d2));U=np.maximum(z*.7/.3,z+c*d2)
        full=(a*b>0)&(alpha-beta>=L)&(alpha+beta<=U)
        selected0=float_membership(a,b,alpha,z);oracle=float_membership(a,b,xtrue,z)
        mu0,hk0,mhk0=intervals(a,b,alpha,z);mut,hkt,mhkt=intervals(a,b,xtrue,z)
        arrays=dict(Z=Z,a=a,b=b,full_regular=full,point0=selected0,point_true=oracle,mu0=mu0,hk0=hk0,mhk0=mhk0,mut=mut,hkt=hkt,mhkt=mhkt)
        npz(out/f'B_{i:02d}.npz',**arrays)
        row=dict(id=i,R=R,s1=s1,s2=s2,sb=sb,true_rho=rhot,thetaA=thetaA,thetaB=thetaB,alpha=alpha,beta=beta,z=z,full_domain_has_nonregular_endpoint=(alpha==beta),rng_initial=initial,rng_final=rng.bit_generator.state)
        metrics={}
        for name,x in [('full_regular_selection',full),('point0_selection',selected0),('true_rho_selection',oracle),('point0_vs_oracle_disagreement',selected0!=oracle),('full_with_mHK0_zero_including',full&(abs(mu0)<=mhk0)),('point0_with_mHK0_zero_including',selected0&(abs(mu0)<=mhk0))]:
            mean=float(np.mean(x));metrics[name]={'estimate':mean,'mcse':math.sqrt(mean*(1-mean)/R),'denominator':R}
        for prefix,mu,hk,mhk in [('rho0',mu0,hk0,mhk0),('true_rho',mut,hkt,mhkt)]:
            for method,width in [('HK',hk),('mHK',mhk)]:
                cov=(abs(mu-thetaA)<=width) if thetaA==thetaB else None
                covp=float(cov.mean()) if cov is not None else None
                metrics[f'{prefix}_{method}']={'coverage':covp,'coverage_mcse':math.sqrt(covp*(1-covp)/R) if covp is not None else None,'coverage_denominator':R if cov is not None else None,'mean_length':float((2*width).mean()),'mean_length_mcse':float((2*width).std(ddof=1)/math.sqrt(R)),'median_length':float(np.median(2*width)),'zero_exclusion':float((abs(mu)>width).mean())}
        row['metrics']=metrics;rows.append(row)
    save(out/'summary_benchmark.json',rows);print('B summary-estimate benchmark saved',flush=True)

def component_c(out,pilot):
    R=20 if pilot else 1000;rng=np.random.default_rng(26101105 if pilot else 26101106);rows=[]
    for idx,(P,m,icc,pattern) in enumerate(product([2,4,8],[5,25],[0.,.5,.9],['null','balanced','unbalanced'])):
        J=24;intercept=np.linspace(-2,2,J);effect=np.zeros(J)
        if pattern=='balanced':effect[:3]=math.log(1.5);effect[3:6]=-math.log(1.5)
        if pattern=='unbalanced':effect[:6]=math.log(1.5)
        truth=effect-effect.mean();nonnull=abs(truth)>1e-12
        pvals=np.empty((R,3,J));latentcover=np.empty((R,J),dtype=bool);fdp=np.empty((R,3));power=np.full((R,3),np.nan);fwer=np.empty((R,3),bool);rawreject=np.empty((R,3));calls=np.empty((R,3,J),bool)
        counts_saved=[];latent_saved=[];rng_initial=rng.bit_generator.state
        sufficient={key:np.empty(shape) for key,shape in [('welch_means',(R,2,J)),('welch_centered_ss',(R,2,J)),('observed_pair_contrasts',(R,P,J)),('latent_pair_contrasts',(R,P,J))]}
        for rep in range(R):
            unit=rng.normal(size=(P,2,1,J))*math.sqrt(icc)
            latent=intercept+unit+rng.normal(size=(P,2,m,J))*math.sqrt(1-icc)
            latent[:,1,:,:]+=effect
            prob=np.exp(latent-latent.max(axis=-1,keepdims=True));prob/=prob.sum(axis=-1,keepdims=True)
            counts=rng.multinomial(10000,prob.reshape(-1,J)).reshape(P,2,m,J)
            latentclr=latent-latent.mean(axis=-1,keepdims=True)
            pv,lm,half,stats_saved=count_tests(counts,latentclr);pvals[rep]=pv
            for key in sufficient:sufficient[key][rep]=stats_saved[key]
            latentcover[rep]=abs(lm-truth)<=half
            for j in range(3):
                selected=bh(pv[j]);calls[rep,j]=selected;false=int(np.sum(selected&~nonnull));total=int(selected.sum())
                fdp[rep,j]=false/max(1,total);fwer[rep,j]=false>0;rawreject[rep,j]=np.mean(pv[j]<=.05)
                if nonnull.any():power[rep,j]=np.sum(selected&nonnull)/int(nonnull.sum())
            if pilot:counts_saved.append(counts.astype(np.uint32));latent_saved.append(latentclr)
        arrays=dict(pvalues=pvals,latent_coverage=latentcover,fdp=fdp,power=power,fwer=fwer,raw_rejection=rawreject,calls=calls,truth=truth)
        arrays.update(sufficient);arrays['welch_n']=np.array([P*m,P*m],dtype=np.int64)
        if pilot:arrays.update(counts=np.stack(counts_saved),latent_clr=np.stack(latent_saved))
        npz(out/f'C_{idx:02d}.npz',**arrays)
        methods=[]
        for j,name in enumerate(['sample_Welch','unit_paired_t','exact_pair_randomization']):
            stats={}
            for key,array in [('FDR',fdp),('FWER',fwer),('unadjusted_rejection_fraction',rawreject),('power',power)]:
                v=array[:,j].astype(float)
                stats[key]={'mean':float(v.mean()),'mcse':float(v.std(ddof=1)/math.sqrt(R)),'denominator':R} if np.isfinite(v).all() else {'mean':None,'mcse':None,'denominator':None}
            methods.append({'method':name,**stats})
        coverage=latentcover.mean(axis=1)
        rows.append({'id':idx,'R':R,'pairs':P,'subsamples':m,'latent_ICC':icc,'pattern':pattern,'taxa':J,'read_depth':10000,'nonnull_CLR_targets':int(nonnull.sum()),'truth':truth.tolist(),'methods':methods,'latent_unit_t_coverage':{'mean':float(coverage.mean()),'mcse':float(coverage.std(ddof=1)/math.sqrt(R)),'denominator':R},'rng_initial':rng_initial,'rng_final':rng.bit_generator.state})
        print(f'C cell {idx+1}/54 {pattern} P={P} m={m} ICC={icc}',flush=True)
    save(out/'count_benchmark.json',rows)

def load_cohort(name):
    meta=read(D3/'inputs'/f'{name}_TRAIN_METADATA.json');runs=[r['run'] for r in meta];positions={r:i for i,r in enumerate(runs)}
    dataset=name[0]
    with gzip.open(D3/'inputs'/f'{dataset}_COUNTS.tsv.gz','rt',encoding='utf-8',newline='') as f:
        reader=csv.reader(f,delimiter='\t');labels=next(reader)[1:];counts=np.zeros((len(meta),len(labels)),np.int64);seen=set()
        for row in reader:
            if row[0] not in positions:continue
            if row[0] in seen:raise ValueError('Duplicate run')
            seen.add(row[0]);counts[positions[row[0]]]=np.asarray(row[1:],dtype=np.int64)
    if seen!=set(runs):raise ValueError('Missing runs')
    names,contrasts=unit_contrasts(ranks(counts),meta)
    return {'id':name,'counts':counts,'metadata':meta,'labels':labels,'pair_names':names,'pair_contrasts':contrasts,'contrast':contrasts.mean(0)}

def component_d(out,pilot):
    cohorts={n:load_cohort(n) for n in ['A2018','A2019','B','B_Inbred','B_Hybrid','C2017','C2018']}
    settings=[(.10,10)] if pilot else list(product([.05,.10,.20],[5,10,20]))
    rows=[];trainings=[];bindings=[]
    for p in sorted((D3/'inputs').iterdir()):
        if p.name=='LABEL_RULE_SOURCE.json':continue
        bindings.append({'path':str(p.relative_to(AUTH)),'sha256':sha(p)})
    for fold,src,targets in FOLDS:
        hist=read(D3/'actual'/fold/'TRAINING_FREEZE.json');allowed={r['taxon']:r['label_decision'] for r in hist['candidates']}
        for prevalence,k in settings:
            frozen=freeze_training([cohorts[n] for n in src],allowed,prevalence,k)
            trainings.append({'fold':fold,'sources':src,'prevalence':prevalence,'k':k,**frozen})
            for target in targets:
                c=cohorts[target];score=target_score(frozen,c['counts'],c['labels'])
                row={'fold':fold,'sources':src,'target':target,'prevalence':prevalence,'k':k,'transductive_target_scope':target,'status':'SCORED' if score is not None else 'UNAVAILABLE','selected':frozen['selected']}
                if score is not None:
                    names,contrasts=unit_contrasts(score,c['metadata'])
                    row.update(equal_pair_mean=float(contrasts.mean()),pairs=[{'pair':n,'contrast':float(v)} for n,v in zip(names,contrasts)])
                else:row.update(equal_pair_mean=None,pairs=None)
                rows.append(row)
    save(out/'rank_ablations.json',rows);save(out/'training_definitions.json',trainings);save(out/'ABC_INPUT_BINDINGS.json',bindings)
    if pilot:
        # Complete input projection for independent pilot calculation without shared loader.
        for name,c in cohorts.items():npz(out/f'D_{name}.npz',counts=c['counts'],contrast=c['contrast'],pair_contrasts=c['pair_contrasts'])
        save(out/'D_cohort_metadata.json',{n:{k:c[k] for k in ('metadata','labels','pair_names')} for n,c in cohorts.items()})
    print('D all fixed training/target records saved',flush=True)

def main():
    stage=sys.argv[1];out=Path(sys.argv[2]);pilot=stage=='pilot';assert stage in ('pilot','formal');check_sources()
    save(out/'ENVIRONMENT.json',{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'threads':{k:os.environ.get(k) for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')},'stage':stage,'protocol_sha256':sha(ROOT/'RESEARCH_PROTOCOL_V1.txt'),'clarifications_sha256':sha(ROOT/'PROTOCOL_CLARIFICATIONS_V1.txt')})
    start=time.perf_counter()
    component_a(out,pilot);component_b(out,pilot);component_c(out,pilot);component_d(out,pilot)
    check_sources();save(out/'RESULT.json',{'status':'COMPUTED_AUDIT_REQUIRED','elapsed_seconds':time.perf_counter()-start,'pilot':pilot,'components':['A','B','C','D'],'science_failure':False})
if __name__=='__main__':main()
