"""Prospectively seeded raw hierarchical counts; no estimator or selection here."""
import gzip,hashlib,itertools,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def cells():
    r=[dict(pairs=g,icc=c,year_rho=r,profile=p) for g,c,r,p in itertools.product([4,12],[.2,.8],[-.5,.8],['global_null','sparse_shared','sparse_heterogeneous'])]
    return r+[dict(pairs=g,icc=.8,year_rho=.8,profile='dense_shared') for g in [4,12]]
def main():
    scale=sys.argv[1];assert scale in ('pilot','formal');out=ROOT/'data'/scale;out.mkdir(parents=True,exist_ok=True)
    reps=2 if scale=='pilot' else 200;m=100;s=4
    mu=np.random.default_rng(117).permutation(np.linspace(-3,3,m))
    records=[]
    for ci,c in enumerate(cells()):
        g=c['pairs'];betaA=np.zeros(m);betaB=np.zeros(m)
        if c['profile']!='global_null':betaA[:(60 if c['profile']=='dense_shared' else 20)]=.8;betaB=betaA.copy()
        if c['profile']=='sparse_heterogeneous':betaB[:10]=-.8;betaB[10:20]=1.6
        for rep in range(reps):
            seed=[20261011,int(scale=='formal'),ci,rep];rng=np.random.default_rng(np.random.SeedSequence(seed))
            pairA=rng.normal(0,.3,(g,m));pairB=rng.normal(0,.3,(g,m))
            u1=rng.normal(size=(g,2,m));u2=c['year_rho']*u1+np.sqrt(1-c['year_rho']**2)*rng.normal(size=u1.shape);uB=rng.normal(size=u1.shape)
            obj={'cell':ci,'replicate':rep,'seed':seed,'parameters':c,'subsamples':s,'taxa':m,'depth':10000,'beta_A':betaA.tolist(),'beta_B':betaB.tolist()}
            for name,u,pair,beta in [('A1',u1,pairA,betaA),('A2',u2,pairA,betaA),('B',uB,pairB,betaB)]:
                eta=mu[None,None,None,:]+pair[:,None,None,:]+np.arange(2)[None,:,None,None]*beta[None,None,None,:]+np.sqrt(c['icc'])*u[:,:,None,:]+rng.normal(0,np.sqrt(1-c['icc']),(g,2,s,m))
                E=eta.reshape(2*g*s,m);P=np.exp(E-E.max(axis=1,keepdims=True));P/=P.sum(axis=1,keepdims=True)
                counts=np.array([rng.multinomial(10000,p) for p in P],dtype=np.int32)
                latent=(E-E.mean(axis=1,keepdims=True)).reshape(2*g,s,m).mean(axis=1)
                obj[name]={'counts':counts.tolist(),'latent_unit_clr':latent.tolist()}
            p=out/f'c{ci:02d}_r{rep:04d}.json.gz'
            with p.open('xb') as raw:
                with gzip.GzipFile(fileobj=raw,mode='wb',mtime=0) as f:f.write(json.dumps(obj,separators=(',',':'),allow_nan=False).encode())
            records.append({'path':p.relative_to(ROOT).as_posix(),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'cell':ci,'replicate':rep})
        print(json.dumps({'completed_cell':ci,'datasets':len(records)}),flush=True)
    with (out/'MANIFEST.json').open('x') as f:json.dump({'scale':scale,'cells':cells(),'replicates_per_cell':reps,'datasets':records,'rng':'NumPy PCG64/SeedSequence','numpy':np.__version__},f,indent=2)
if __name__=='__main__':main()
