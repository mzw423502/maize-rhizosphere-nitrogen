"""New pure methods. No imports from historical workers or auditors."""
from fractions import Fraction as F
import numpy as np
from scipy.stats import t as student_t, rankdata
from itertools import product

def certificate(a,b,alpha,beta,z,q=F(7,10),support=True):
    vals=[a,b,alpha,beta,z]
    if any(v is None for v in vals):return {'class':'NOT_ASSESSABLE','reason':'MISSING_PARAMETERS'}
    a,b,alpha,beta,z=map(F,vals);q=F(q)
    if not F(1,2)<=q<1 or z<=0 or beta<0 or alpha<beta:return {'class':'INVALID_INPUT'}
    if alpha==0:return {'class':'NOT_ASSESSABLE','reason':'EMPTY_REGULAR_DOMAIN'}
    degen=alpha==beta
    regular_lo=F(-1);regular_hi=F(1)
    if support is None:return {'class':'NOT_ASSESSABLE','reason':'UNKNOWN_SUPPORT'}
    c=2*q-1;d2=(a-b)**2
    L=max(F(0),min(z*(1-q)/q,z-c*d2));U=max(z*q/(1-q),z+c*d2)
    out={'variance_window':[str(L),str(U)],'degenerate_minus_one':degen}
    if a*b<=0 or not support:return {**out,'class':'NONE','rho_window':None,'reason':'SIGN_OR_SUPPORT'}
    if beta==0:
        yes=L<=alpha<=U
        return {**out,'class':'ALL_RHO' if yes else 'NONE','rho_window':['-1','1'] if yes else None,'closed':[True,True] if yes else None}
    lo=max(regular_lo,(L-alpha)/beta);hi=min(regular_hi,(U-alpha)/beta)
    if lo>hi or hi==F(-1) and degen:return {**out,'class':'NONE','rho_window':None}
    all_regular=lo==F(-1) and hi==F(1)
    return {**out,'class':('ALL_REGULAR_RHO_WITH_DEGENERATE_ENDPOINT' if degen else 'ALL_RHO') if all_regular else 'SOME_RHO','rho_window':[str(lo),str(hi)],'closed':[not(degen and lo==-1),True]}

def float_membership(a,b,x,z,q=.7):
    tau=np.maximum(0,((a-b)**2-x-z)/2)
    wa=(z+tau)/(x+z+2*tau)
    return (a*b>0)&(np.maximum(wa,1-wa)<=q)

def intervals(a,b,x,z):
    tau=np.maximum(0,((a-b)**2-x-z)/2)
    wa=1/(x+tau);wb=1/(z+tau);W=wa+wb
    mu=(wa*a+wb*b)/W
    h=wa*(a-mu)**2+wb*(b-mu)**2
    critical=np.tan(np.pi*.475)
    return mu,critical*np.sqrt(h/W),critical*np.sqrt(np.maximum(1,h)/W)

def bh(p):
    p=np.asarray(p);order=np.argsort(p,kind='stable');threshold=.05*np.arange(1,len(p)+1)/len(p)
    ok=p[order]<=threshold
    chosen=np.zeros(len(p),dtype=bool)
    if ok.any():chosen[order[:np.where(ok)[0][-1]+1]]=True
    return chosen

def paired_p(diff):
    n=diff.shape[0];mean=diff.mean(0);se=diff.std(0,ddof=1)/np.sqrt(n)
    stat=np.divide(np.abs(mean),se,out=np.zeros_like(mean),where=se>0)
    p=2*student_t.sf(stat,n-1)
    p=np.where(se>0,p,np.where(mean==0,1.,0.))
    return p,mean,se

def count_tests(counts,latent):
    P,_,m,J=counts.shape
    logc=np.log(counts+.5);clr=logc-logc.mean(axis=-1,keepdims=True)
    ln=clr[:,1,:,:].reshape(P*m,J);hn=clr[:,0,:,:].reshape(P*m,J)
    mean=ln.mean(0)-hn.mean(0);vl=ln.var(0,ddof=1)/len(ln);vh=hn.var(0,ddof=1)/len(hn);vv=vl+vh
    stat=np.divide(abs(mean),np.sqrt(vv),out=np.zeros(J),where=vv>0)
    df=np.divide(vv**2,vl**2/(len(ln)-1)+vh**2/(len(hn)-1),out=np.ones(J),where=vv>0)
    naive=np.where(vv>0,2*student_t.sf(stat,df),np.where(mean==0,1.,0.))
    units=clr.mean(2);diff=units[:,1]-units[:,0]
    paired,_,_=paired_p(diff)
    signs=np.asarray(list(product((-1,1),repeat=P)))
    obs=abs(diff.sum(0));tol=1e-12*np.maximum(1,obs)
    perm=(abs(signs@diff)>=obs-tol).mean(0)
    lp=latent.mean(2);ld=lp[:,1]-lp[:,0]
    _,lm,ls=paired_p(ld);half=student_t.ppf(.975,P-1)*ls
    sufficient={'welch_means':np.stack([hn.mean(0),ln.mean(0)]),
                'welch_centered_ss':np.stack([((hn-hn.mean(0))**2).sum(0),((ln-ln.mean(0))**2).sum(0)]),
                'observed_pair_contrasts':diff,'latent_pair_contrasts':ld}
    return np.stack([naive,paired,perm]),lm,half,sufficient

def ranks(counts):
    raw=(rankdata(counts,axis=0,method='average')-1)/(len(counts)-1)
    sd=raw.std(0,ddof=1)
    return np.divide(raw-raw.mean(0),sd,out=np.zeros_like(raw),where=sd>0)

def unit_contrasts(values,metadata):
    """values can be samples x taxa or a sample score; metadata only downstream."""
    units={}
    for i,r in enumerate(metadata):
        u=units.setdefault(r['unit'],{'pair':r['pair'],'arm':r['arm'],'cells':{}})
        if (u['pair'],u['arm'])!=(r['pair'],r['arm']):raise ValueError('Conflicting unit')
        u['cells'].setdefault(r['cell'],[]).append(i)
    pairs={}
    for name,u in sorted(units.items()):
        v=np.mean([np.mean(values[ix],axis=0) for _,ix in sorted(u['cells'].items())],axis=0)
        if u['arm'] in pairs.setdefault(u['pair'],{}):raise ValueError('Duplicate pair arm')
        pairs[u['pair']][u['arm']]=v
    names=sorted(pairs)
    if any(set(pairs[n])!={'LN','HN'} for n in names):raise ValueError('Missing pair arm')
    return names,np.stack([pairs[n]['LN']-pairs[n]['HN'] for n in names])

def freeze_training(cohorts,allowed,prevalence,k):
    """No target data or target labels in API."""
    labels=cohorts[0]['labels']
    assert all(c['labels']==labels for c in cohorts)
    n=sum(len(c['metadata']) for c in cohorts)
    ge=sum((c['counts']>=10).sum(0) for c in cohorts)
    contrasts=np.mean([c['contrast'] for c in cohorts],axis=0)
    candidate=[i for i,l in enumerate(labels) if allowed[l]=='INCLUDED_BY_REVIEW_RULE' and ge[i]>=prevalence*n-1e-12 and np.isfinite(contrasts[i]) and contrasts[i]!=0]
    candidate.sort(key=lambda i:(-abs(contrasts[i]),labels[i].encode('utf-8')))
    idx=candidate[:k]
    if len(idx)<3:return {'status':'UNAVAILABLE_FEWER_THAN_3','selected':[]}
    denominator=sum(abs(contrasts[i]) for i in idx)
    return {'status':'FROZEN','selected':[{'taxon':labels[i],'contrast':float(contrasts[i]),'weight':float(contrasts[i]/denominator)} for i in idx]}

def target_score(freeze,counts,labels):
    """Only frozen source definition and unlabeled complete named target counts."""
    if freeze['status']!='FROZEN':return None
    if any(r['taxon'] not in labels for r in freeze['selected']):return None
    ix=[labels.index(r['taxon']) for r in freeze['selected']]
    z=ranks(counts[:,ix]);w=np.asarray([r['weight'] for r in freeze['selected']])
    return z@w
