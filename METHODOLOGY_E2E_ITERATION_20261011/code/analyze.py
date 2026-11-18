"""Saved count/estimator inputs -> decision records. Streaming by dataset."""
import gzip,json,math,sys
from pathlib import Path
import numpy as np
from scipy import stats
from certificate import classify,max_weight,pooled,box_cap
ROOT=Path(__file__).resolve().parents[1]
TCRIT=stats.t.ppf(.975,1)
def readgz(p):
    with gzip.open(p,'rt') as f:return json.load(f)
def bh(p):
    idx=np.argsort(p);v=np.minimum.accumulate((p[idx]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];ans=np.empty(len(p));ans[idx]=np.minimum(1,v);return ans
def metrics(b,se,df,truth):
    p=2*stats.t.sf(np.abs(b/se),df);reject=bh(p)<=.05;null=np.isclose(truth,0,atol=1e-14);alt=~null
    rad=stats.t.ppf(.975,df)*se
    return dict(bias=float(np.mean(b-truth)),mse=float(np.mean((b-truth)**2)),coverage=float(np.mean(np.abs(b-truth)<=rad)),mean_width=float(np.mean(2*rad)),type1=float(np.mean(p[null]<=.05)) if null.any() else None,power=float(np.mean(p[alt]<=.05)) if alt.any() else None,BH_FDP=float((reject&null).sum()/max(1,reject.sum())),BH_power=float(np.mean(reject[alt])) if alt.any() else None,null_count=int(null.sum()),alternative_count=int(alt.sum()),reject_count=int(reject.sum()))
def estimates(dat,fit,method,year):
    g=dat['parameters']['pairs'];s=dat['subsamples'];n=2*g*s
    U=np.array(fit['years'][year]['unit_clr']);D=U[1::2]-U[::2]
    if method=='LinDA_unit':
        native=fit['years'][year]['native'];b=np.array([x['log2FoldChange'] for x in native])*np.log(2);se=np.array([x['lfcSE'] for x in native])*np.log(2);df=g-1
    elif method=='CLR_unit':b=D.mean(axis=0);se=D.std(axis=0,ddof=1)/np.sqrt(g);df=g-1
    else:
        W=np.log(np.array(dat[year]['counts'])+.5);W-=W.mean(axis=1,keepdims=True);C=W.reshape(g,2,s,100);b=(C[:,1].mean(axis=(0,1))-C[:,0].mean(axis=(0,1)))
        resid=C-C.mean(axis=(1,2),keepdims=True)-(np.arange(2)[None,:,None,None]-.5)*b[None,None,None,:];df=n-g-1;se=np.sqrt((resid**2).sum(axis=(0,1,2))/df/(n/4))
    return b,se,df,D
def main():
    scale=sys.argv[1];manifest=json.loads((ROOT/'data'/scale/'MANIFEST.json').read_text());out=ROOT/'analysis'/scale;out.mkdir(parents=True,exist_ok=True)
    summary=[];empirical={};methods=['LinDA_unit','CLR_unit','CLR_sample_independence_control'];qs=[.6,.7,.8,.9]
    for ri,rec in enumerate(manifest['datasets']):
        dat=readgz(ROOT/rec['path']);fp=ROOT/'estimates'/scale/(Path(rec['path']).name[:-3]);fit=json.loads(fp.read_text());assert fit['input_sha256']==rec['sha256']
        report={'cell':dat['cell'],'replicate':dat['replicate'],'input_sha256':rec['sha256'],'methods':{}}
        for method in methods:
            ann={};D={};truth={}
            for y in ['A1','A2','B']:
                b,se,df,D[y]=estimates(dat,fit,method,y);tr=np.array(dat['beta_B' if y=='B' else 'beta_A']);tr=tr if method=='LinDA_unit' else tr-tr.mean();truth[y]=tr
                assert np.isfinite(b).all() and np.isfinite(se).all() and (se>0).all()
                effect_scale='absolute_log_mode_zero' if method=='LinDA_unit' else 'centered_log_ratio'
                ann[y]={'effect_scale':effect_scale,'truth':tr.tolist(),'effect':b.tolist(),'se':se.tolist(),'df':df,'metrics':metrics(b,se,df,tr)}
                summary.append(dict(cell=dat['cell'],replicate=dat['replicate'],method=method,effect_scale=effect_scale,year=y,**ann[y]['metrics']))
                key=(dat['cell'],method,y);empirical.setdefault(key,[]).append((b,se,tr))
            a=(np.array(ann['A1']['effect'])+np.array(ann['A2']['effect']))/2;b=np.array(ann['B']['effect']);s1=np.array(ann['A1']['se']);s2=np.array(ann['A2']['se']);sb=np.array(ann['B']['se'])
            linked_var=((D['A1']+D['A2'])/2).var(axis=0,ddof=1)/dat['parameters']['pairs']
            true_mean=(truth['A1']+truth['B'])/2
            decisions=[]
            for j in range(100):
                cs=[classify(float(a[j]),float(b[j]),float(s1[j]),float(s2[j]),float(sb[j]),q) for q in qs]
                alpha=(s1[j]**2+s2[j]**2)/4;beta=s1[j]*s2[j]/2;z=sb[j]**2;d=a[j]-b[j]
                gridx=np.maximum(0,alpha+beta*np.linspace(-1,1,101));gridw=.5+np.abs(gridx-z)/(2*np.maximum(gridx+z,d*d))
                intervals=[]
                for rho in [-1,0,1]:
                    x=((s1[j]-s2[j])**2/4 if rho==-1 else (s1[j]+s2[j])**2/4 if rho==1 else alpha);p=pooled(a[j],b[j],x,z)
                    if p:
                        for label in ['HK','mHK']:
                            width=2*TCRIT*p[label+'_se'];intervals.append(dict(status='REGULAR',rho=rho,method=label,width=float(width),covered=bool(abs(p['mean']-true_mean[j])<=width/2),mean=float(p['mean']),tau2=float(p['tau2'])))
                    else:
                        for label in ['HK','mHK']:intervals.append(dict(status='NONREGULAR',rho=rho,method=label,width=None,covered=None,mean=None,tau2=None))
                decisions.append(dict(taxon=j,truth_A=float(truth['A1'][j]),truth_B=float(truth['B'][j]),known_opposite=bool(truth['A1'][j]*truth['B'][j]<0),linked_CLR_year_variance=float(linked_var[j]),caps=cs,grid101_all=[bool((gridw<=q).all()) for q in qs],grid101_any=[bool((gridw<=q).any()) for q in qs],intervals=intervals))
            report['methods'][method]={'annual':ann,'taxa':decisions,'warnings':{y:fit['years'][y]['warnings'] for y in ['A1','A2','B']} if method=='LinDA_unit' else {}}
        target=out/Path(rec['path']).name
        with target.open('xb') as raw:
            with gzip.GzipFile(fileobj=raw,mode='wb',mtime=0) as f:f.write(json.dumps(report,separators=(',',':'),allow_nan=False).encode())
        if (ri+1)%10==0:print(json.dumps({'analyzed':ri+1,'total':len(manifest['datasets'])}),flush=True)
    # Dataset-averaged outcomes: MCSE uses independent datasets, not individual taxa.
    agg=[]
    for ci in range(len(manifest['cells'])):
      for method in methods:
       for y in ['A1','A2','B']:
        rows=[r for r in summary if r['cell']==ci and r['method']==method and r['year']==y];statsout={}
        for name in ['bias','mse','coverage','mean_width','type1','power','BH_FDP','BH_power']:
            v=np.array([r[name] for r in rows if r[name] is not None]);statsout[name]={'mean':float(v.mean()),'MCSE':float(v.std(ddof=1)/np.sqrt(len(v))),'n_datasets':len(v)} if len(v)>1 else None
        a=np.array(empirical[(ci,method,y)]);emp_sd=a[:,0,:].std(axis=0,ddof=1)
        statsout['mean_SE_to_empirical_SD']=float(np.mean(a[:,1,:].mean(axis=0)/emp_sd))
        agg.append(dict(cell=ci,parameters=manifest['cells'][ci],method=method,effect_scale='absolute_log_mode_zero' if method=='LinDA_unit' else 'centered_log_ratio',year=y,metrics=statsout))
    with (out/'ESTIMATOR_SUMMARY.json').open('x') as f:json.dump(agg,f,indent=2,allow_nan=False)
    with (out/'DATASET_ESTIMATOR_METRICS.json').open('x') as f:json.dump(summary,f,separators=(',',':'),allow_nan=False)
    corr=[]
    for ci in range(len(manifest['cells'])):
      for method in methods:
        e1=np.array(empirical[(ci,method,'A1')])[:,0,:];e2=np.array(empirical[(ci,method,'A2')])[:,0,:]
        v1=e1-e1.mean(axis=0);v2=e2-e2.mean(axis=0);r=(v1*v2).sum(axis=0)/np.sqrt((v1*v1).sum(axis=0)*(v2*v2).sum(axis=0))
        corr.append(dict(cell=ci,method=method,latent_unit_rho=manifest['cells'][ci]['year_rho'],empirical_effect_error_rho=r.tolist(),n_datasets=len(e1)))
    with (out/'EMPIRICAL_YEAR_ERROR_CORRELATION.json').open('x') as f:json.dump(corr,f,indent=2,allow_nan=False)
    rng=np.random.default_rng(202610110 if scale=='pilot' else 202610111);boxcases=[]
    for i in range(256 if scale=='pilot' else 4096):
        xx=sorted(np.exp(rng.uniform(-4,3,2)));zz=sorted(np.exp(rng.uniform(-4,3,2)));dd=sorted(rng.uniform(-5,5,2));boxcases.append(dict(x=xx,z=zz,d=dd,cap=box_cap(xx,zz,dd)))
    with (out/'BOX_CASES.json').open('x') as f:json.dump(boxcases,f,indent=2)
    with (out/'COMPLETION.json').open('x') as f:json.dump({'scale':scale,'datasets':len(manifest['datasets']),'methods':methods,'taxa_per_method_dataset':100,'caps':qs,'science_status':'PENDING_INDEPENDENT_AUDIT'},f,indent=2)
if __name__=='__main__':main()
