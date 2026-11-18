"""Independent saved-input audit. Does not import producer or certificate modules."""
import gzip,hashlib,json,math,sys
from pathlib import Path
import numpy as np
from scipy import stats,optimize
ROOT=Path(__file__).resolve().parents[1]
TCRIT=stats.t.ppf(.975,1)
class Audit:
    def __init__(self):self.n=0;self.bad=[];self.maxabs=0;self.mode_maxabs=0
    def eq(self,a,b,label):
        self.n+=1
        if a!=b and len(self.bad)<50:self.bad.append(dict(label=label,observed=str(a),expected=str(b)))
    def close(self,a,b,label,tol=1e-10):
        a=np.asarray(a,dtype=float);b=np.asarray(b,dtype=float);self.n+=a.size
        if a.shape!=b.shape:self.bad.append(dict(label=label,shape=[str(a.shape),str(b.shape)]));return
        delta=np.abs(a-b);maximum=float(delta.max(initial=0))
        if tol==1e-7:self.mode_maxabs=max(self.mode_maxabs,maximum)
        else:self.maxabs=max(self.maxabs,maximum)
        if not np.all(np.isfinite(a)&np.isfinite(b)&(delta<=tol*(1+np.abs(b)))) and len(self.bad)<50:self.bad.append(dict(label=label,maxabs=maximum,tolerance=tol))
def readgz(p):
    with gzip.open(p,'rt') as f:return json.load(f)
def mode_definition(x):
    # Independent Python definition of upstream modeest mean-shift, not a second R call.
    x=np.asarray(x);n=len(x);y=np.sort(x);k=math.ceil(n/2)-1;spans=y[k:]-y[:n-k];starts=np.flatnonzero(spans==spans.min())
    # Continuous simulated coefficients have unique shortest intervals. A tie is not silently resolved.
    if len(starts)!=1:raise ValueError('Mean-shift initial interval tie requires explicit review')
    start=starts[0];par=float(y[start:start+k+1].mean());sd=x.std(ddof=1);lo=min(sd,float(np.subtract(*np.quantile(x,[.75,.25])))/1.34)
    if lo==0:lo=sd or abs(x[0]) or 1
    bw=.9*lo*n**(-.2)
    for _ in range(1000):
        w=np.exp(-.5*((x-par)/bw)**2);nxt=float((x*w).sum()/w.sum())
        if abs(nxt/par-1)<math.sqrt(np.finfo(float).eps):break
        par=nxt
    return par
def OLS(Y,X):
    # Full design-matrix least squares independently of paired-contrast producer.
    coef=np.linalg.lstsq(X,Y,rcond=None)[0];res=Y-X@coef;df=X.shape[0]-X.shape[1];se=np.sqrt((res*res).sum(axis=0)/df*np.linalg.inv(X.T@X)[-1,-1]);return coef[-1],se,df
def weights(x,z,d):
    tau=(d*d-x-z)/2
    if tau<0:tau=0
    if x==0 and tau==0:return 1.0
    prec=np.array([1/(x+tau),1/(z+tau)]);return float(max(prec)/sum(prec))
def admissible_pieces(z,d,q,xmin,xmax):
    D=d*d;c=2*q-1;pieces=[]
    lo=max(xmin,0,D-z,z*(1-q)/q);hi=min(xmax,z*q/(1-q))
    if lo<=hi and hi>0:pieces.append((lo,hi))
    if D>0:
        lo=max(xmin,0,z-c*D);hi=min(xmax,D-z,z+c*D)
        if lo<=hi and hi>0:pieces.append((lo,hi))
    return pieces
def statmetrics(b,se,df,tr):
    p=2*stats.t.cdf(-np.abs(b/se),df);order=np.argsort(p);adjust=np.empty(len(p));v=p[order]*len(p)/np.arange(1,len(p)+1)
    for i in range(len(p)-2,-1,-1):v[i]=min(v[i],v[i+1])
    adjust[order]=np.minimum(1,v);rej=adjust<=.05;null=np.abs(tr)<1e-14;alt=~null;width=2*stats.t.ppf(.975,df)*se
    return dict(bias=float(np.average(b-tr)),mse=float(np.average((b-tr)**2)),coverage=float(np.average(np.abs(b-tr)<=width/2)),mean_width=float(np.average(width)),type1=float(np.average(p[null]<=.05)) if null.any() else None,power=float(np.average(p[alt]<=.05)) if alt.any() else None,BH_FDP=float(np.count_nonzero(rej&null)/max(1,np.count_nonzero(rej))),BH_power=float(np.average(rej[alt])) if alt.any() else None,null_count=int(null.sum()),alternative_count=int(alt.sum()),reject_count=int(rej.sum()))
def main():
    scale=sys.argv[1];out=ROOT/'analysis'/scale;A=Audit();manifest=json.loads((ROOT/'data'/scale/'MANIFEST.json').read_text());summary=[];errs={}
    for ri,rec in enumerate(manifest['datasets']):
        p=ROOT/rec['path'];A.eq(hashlib.sha256(p.read_bytes()).hexdigest(),rec['sha256'],'raw hash');dat=readgz(p);fit=json.loads((ROOT/'estimates'/scale/Path(rec['path']).name[:-3]).read_text());report=readgz(out/p.name)
        A.eq((dat['cell'],dat['replicate']),(rec['cell'],rec['replicate']),'keys');A.eq(fit['input_sha256'],rec['sha256'],'native input identity');g=dat['parameters']['pairs'];s=dat['subsamples'];m=100;components={};diffs={}
        pair=np.repeat(np.arange(g),2);arm=np.tile([0.,1.],g);X=np.column_stack([np.eye(g)[pair],arm]);XS=np.repeat(X,s,axis=0)
        for y in ['A1','A2','B']:
            C=np.array(dat[y]['counts']);A.eq(C.shape,(2*g*s,m),'count shape');A.eq(bool(np.issubdtype(C.dtype,np.integer) and (C>=0).all()),True,'integer/nonnegative');A.close(C.sum(axis=1),np.full(2*g*s,10000),'depth',tol=0)
            L=np.log(C+.5);L-=L.mean(axis=1)[:,None];U=np.stack([L[i*s:(i+1)*s].sum(axis=0)/s for i in range(2*g)]);D=U[1::2]-U[::2];diffs[y]=D
            n=fit['years'][y];A.close(n['unit_clr'],U,'unit CLR');P=np.array(n['unit_profile']);A.close(P.sum(axis=0),np.ones(2*g),'closure');A.eq(bool((P>0).all()),True,'positive profile');LP=np.log(P);NP=(LP-LP.mean(axis=0)[None,:]).T;A.close(NP,U,'native measurement');A.close(n['native_clr'],U,'saved native CLR');A.close(n['profile_CLR_maxabs'],np.abs(NP-U).max(),'identity maximum')
            bu,su,dfu=OLS(U,X);bs,ss,dfs=OLS(L,XS);A.close(n['raw_coefficient_log2'],bu/np.log(2),'native original coefficient');bias=mode_definition(np.sqrt(2*g)*bu/np.log(2))/np.sqrt(2*g);A.close(n['bias_log2'],bias,'mode correction',tol=1e-7)
            native=n['native'];bn=np.array([r['log2FoldChange'] for r in native])*np.log(2);sn=np.array([r['lfcSE'] for r in native])*np.log(2)
            A.close(bn,bu-float(n['bias_log2'])*np.log(2),'bias arithmetic');A.close(sn,su,'native SE');A.close([r['df'] for r in native],np.full(m,dfu),'native df');A.close([r['pvalue'] for r in native],2*stats.t.sf(abs(bn/sn),dfu),'native p')
            components[y]={'LinDA_unit':(bn,sn,dfu),'CLR_unit':(bu,su,dfu),'CLR_sample_independence_control':(bs,ss,dfs)}
        for method,obj in report['methods'].items():
            for y in ['A1','A2','B']:
                b,se,df=components[y][method];ann=obj['annual'][y];A.close(ann['effect'],b,'reported effect');A.close(ann['se'],se,'reported SE');A.eq(ann['df'],df,'reported df')
                tr=np.array(dat['beta_B' if y=='B' else 'beta_A']);tr=tr if method=='LinDA_unit' else tr-tr.mean();A.close(ann['truth'],tr,'truth');A.eq(ann['effect_scale'],'absolute_log_mode_zero' if method=='LinDA_unit' else 'centered_log_ratio','scale');met=statmetrics(b,se,df,tr)
                for k,v in met.items():
                    if v is None:A.eq(ann['metrics'][k],None,'metric '+k)
                    else:A.close(ann['metrics'][k],v,'metric '+k)
                summary.append(dict(cell=dat['cell'],replicate=dat['replicate'],method=method,effect_scale=ann['effect_scale'],year=y,**met));errs.setdefault((dat['cell'],method,y),[]).append((b,se))
            b1,s1,_=components['A1'][method];b2,s2,_=components['A2'][method];bb,sb,_=components['B'][method];aa=(b1+b2)/2;av=(s1*s1+s2*s2)/4;cov=s1*s2/2;xmin=(s1-s2)**2/4;xmax=(s1+s2)**2/4
            linked=((diffs['A1']+diffs['A2'])/2).var(axis=0,ddof=1)/g
            for j,r in enumerate(obj['taxa']):
                z=sb[j]**2;d=aa[j]-bb[j];same=aa[j]*bb[j]>0;cap=max(weights(xmin[j],z,d),weights(xmax[j],z,d));capiv=max(weights(xmin[j],z,0),weights(xmax[j],z,0));A.eq(r['taxon'],j,'taxon key');A.close(r['linked_CLR_year_variance'],linked[j],'linked covariance')
                ta=np.array(dat['beta_A']);tb=np.array(dat['beta_B'])
                if method!='LinDA_unit':ta=ta-ta.mean();tb=tb-tb.mean()
                A.close([r['truth_A'],r['truth_B']],[ta[j],tb[j]],'decision truths');A.eq(r['known_opposite'],bool(ta[j]*tb[j]<0),'opposite truth')
                gx=np.maximum(0,av[j]+cov[j]*np.linspace(-1,1,101));gt=np.maximum(0,(d*d-gx-z)/2);gw=np.maximum(gx+gt,z+gt)/(gx+z+2*gt)
                for qi,q in enumerate([.6,.7,.8,.9]):
                    c=r['caps'][qi];pieces=admissible_pieces(z,d,q,xmin[j],xmax[j]);allw=cap<=q;anyw=bool(pieces);A.close(c['required_cap'],cap,'minimal cap');A.close(c['required_cap_iv'],capiv,'IV cap');A.eq(c['weight_class'],'ALL' if allw else 'SOME' if anyw else 'NONE','weight class');A.eq(c['decision_class'],'ALL' if same and allw else 'SOME' if same and anyw else 'NONE','decision class');A.eq(c['all_pass'],bool(same and allw),'ALL');A.eq(c['rho0_pass'],bool(same and weights(av[j],z,d)<=q),'point');A.eq(c['heterogeneity_assisted'],bool(allw and capiv>q),'heterogeneity assistance')
                    A.close([c[k] for k in ['a','b','s1','s2','sb','alpha','beta','z','xmin','xmax','q']],[aa[j],bb[j],s1[j],s2[j],sb[j],av[j],cov[j],z,xmin[j],xmax[j],q],'decision inputs')
                    A.eq(r['grid101_any'][qi],bool((gw<=q).any()),'grid SOME')
                    if pieces:
                        A.close(c['rho_low'],max(-1,(min(p[0] for p in pieces)-av[j])/cov[j]),'rho low');A.close(c['rho_high'],min(1,(max(p[1] for p in pieces)-av[j])/cov[j]),'rho high')
                    else:A.eq(c['rho_low'],None,'no domain')
                    # ALL grid identity is exact because endpoints are included; do not claim superiority.
                    A.eq(r['grid101_all'][qi],allw,'grid all identity')
                tr=(r['truth_A']+r['truth_B'])/2
                for iv in r['intervals']:
                    rho=iv['rho'];x=xmin[j] if rho==-1 else xmax[j] if rho==1 else av[j]
                    if x<=0:A.eq(iv['status'],'NONREGULAR','nonregular');A.eq(iv['width'],None,'nonregular width');continue
                    tt=max(0,((aa[j]-bb[j])**2-x-z)/2);wa=1/(x+tt);wb=1/(z+tt);mu=(wa*aa[j]+wb*bb[j])/(wa+wb);Q=wa*(aa[j]-mu)**2+wb*(bb[j]-mu)**2;h=Q if iv['method']=='HK' else max(1,Q);width=2*TCRIT*math.sqrt(h/(wa+wb))
                    A.close(iv['mean'],mu,'meta mean');A.close(iv['tau2'],tt,'tau');A.close(iv['width'],width,'interval width');A.eq(iv['covered'],bool(abs(mu-tr)<=width/2),'coverage')
        if (ri+1)%10==0:print(json.dumps({'audited':ri+1,'comparisons':A.n,'differences':len(A.bad)}),flush=True)
        if A.bad:break
    if not A.bad:
        saved=json.loads((out/'DATASET_ESTIMATOR_METRICS.json').read_text());A.eq(len(saved),len(summary),'summary count')
        for obs,exp in zip(saved,summary):
            for k,v in exp.items():
                if isinstance(v,(float,int)):A.close(obs[k],v,'summary '+k)
                else:A.eq(obs[k],v,'summary '+k)
        aggregates=json.loads((out/'ESTIMATOR_SUMMARY.json').read_text())
        for row in aggregates:
            key=(row['cell'],row['method'],row['year']);rr=[r for r in summary if (r['cell'],r['method'],r['year'])==key]
            for name in ['bias','mse','coverage','mean_width','type1','power','BH_FDP','BH_power']:
                v=np.array([r[name] for r in rr if r[name] is not None]);o=row['metrics'][name]
                if len(v)>1:A.close(o['mean'],v.mean(),'aggregate mean');A.close(o['MCSE'],v.std(ddof=1)/math.sqrt(len(v)),'dataset MCSE')
                else:A.eq(o,None,'not applicable aggregate')
            e=np.array(errs[key]);A.close(row['metrics']['mean_SE_to_empirical_SD'],np.mean(e[:,1].mean(axis=0)/e[:,0].std(axis=0,ddof=1)),'SE calibration')
        for row in json.loads((out/'EMPIRICAL_YEAR_ERROR_CORRELATION.json').read_text()):
            k=(row['cell'],row['method']);e1=np.array(errs[(*k,'A1')])[:,0];e2=np.array(errs[(*k,'A2')])[:,0];values=[np.corrcoef(e1[:,j],e2[:,j])[0,1] for j in range(100)];A.close(row['empirical_effect_error_rho'],values,'empirical error rho')
        boxes=json.loads((out/'BOX_CASES.json').read_text())
        for i,case in enumerate(boxes):
            x0,x1=case['x'];z0,z1=case['z'];d0,d1=case['d'];d=0 if d0<=0<=d1 else sorted([d0,d1],key=abs)[0];expected=max(weights(x,z,d) for x in [x0,x1] for z in [z0,z1]);A.close(case['cap'],expected,'box corners')
            if i<(64 if scale=='pilot' else 512):
                # Independent one-dimensional restricted likelihood optimization at all four corners.
                for x in [x0,x1]:
                  for z in [z0,z1]:
                    S=x+z;D=d*d;opt=optimize.minimize_scalar(lambda v:math.log(v)+D/v,bounds=(S,S+D+1),method='bounded',options={'xatol':1e-12});v=S if math.log(S)+D/S<=opt.fun else opt.x;t=(v-S)/2;pp=np.array([1/(x+t),1/(z+t)]);A.close(max(pp)/sum(pp),weights(x,z,d),'independent likelihood optimizer',tol=1e-6)
                grid=[weights(x,z,h) for x in np.linspace(x0,x1,7) for z in np.linspace(z0,z1,7) for h in [d0,d,d1]];A.eq(max(grid)<=case['cap']+1e-12,True,'box interior dominance')
    result=dict(scale=scale,status='PASS' if not A.bad else 'FAIL',differences=len(A.bad),examples=A.bad,numeric_comparisons=A.n,maxabs=A.maxabs,mode_maxabs=A.mode_maxabs,raw_counts_reconstructed=True,producer_modules_imported=False,native_LinDA_reexecution=False,rng_resampling=False,independent_linear_component_reconstruction=True)
    with (out/'AUDIT_RESULT.json').open('x') as f:json.dump(result,f,indent=2)
    print(json.dumps(result));return bool(A.bad)
if __name__=='__main__':raise SystemExit(main())
