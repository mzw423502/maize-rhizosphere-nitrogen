"""New fixed-split utility analysis from immutable saved annual components only."""
import gzip,hashlib,json,math,sys
from collections import defaultdict
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
PREV=ROOT.parent/'METHODOLOGY_E2E_ITERATION_20261011'
METHODS=['LinDA_unit','CLR_unit'];CAPS=[.6,.7,.8,.9]
DOMAINS=['zero','calibrated_point','calibrated_interval','full']
HALFWIDTHS=[0.,.1,.25,.5,1.,2.]
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def put(p,x):
    with p.open('x',encoding='utf8') as f:json.dump(x,f,ensure_ascii=False,allow_nan=False,indent=2)
def read(ci,rep,hashes):
    p=PREV/'analysis/formal'/f'c{ci:02d}_r{rep:04d}.json.gz'
    h=sha(p);hashes[p.name]=h
    with gzip.open(p,'rt') as f:r=json.load(f)
    assert r['cell']==ci and r['replicate']==rep
    return r
def vec(obj,key,year):return np.asarray(obj['annual'][year][key],dtype=float)
def w(x,z,d):return .5+np.abs(x-z)/(2*np.maximum(x+z,d*d))
def vmax(lo,hi,s1,s2,z,d):
    x0=((s1-s2)**2+2*(1+lo)*s1*s2)/4
    x1=((s1-s2)**2+2*(1+hi)*s1*s2)/4
    return np.maximum(w(x0,z,d),w(x1,z,d))
class Moments:
    def __init__(self):self.n=0;self.mean=0.;self.M2=0.
    def add(self,x):
        if x is None:return
        self.n+=1;d=x-self.mean;self.mean+=d/self.n;self.M2+=d*(x-self.mean)
    def result(self):return {'mean':self.mean if self.n else None,'MCSE':math.sqrt(max(0,self.M2)/(self.n-1)/self.n) if self.n>1 else None,'n_datasets':self.n}
def fraction(a,mask):return float(np.mean(a[mask])) if mask.any() else None
def metrics(sel,truth_a,truth_b,a,base,full,snr,iv,available):
    same=truth_a*truth_b>1e-20;opp=truth_a*truth_b < -1e-20
    null=(np.abs(truth_a)<1e-14)&(np.abs(truth_b)<1e-14);one=~(same|opp|null)
    wrong=same&(a*truth_a<0);removed=base&~full&available
    m={'membership':fraction(sel,available),'true_same_retention':fraction(sel,same&available),'true_same_wrong_sign_selected_rate':fraction(sel&wrong,same&available),'true_opposite_membership':fraction(sel,opp&available),'null_membership':fraction(sel,null&available),'one_zero_membership':fraction(sel,one&available),'opposite_fraction_among_selected':fraction(opp,sel),'null_fraction_among_selected':fraction(null,sel),'paired_difference_from_zero':fraction(sel.astype(int)-base.astype(int),available),'full_loss_fraction_of_this_strategy':fraction(~full,sel),'true_same_paired_difference_from_zero':fraction(sel.astype(int)-base.astype(int),same&available),'true_opposite_paired_difference_from_zero':fraction(sel.astype(int)-base.astype(int),opp&available),'removed_from_zero_true_same_fraction':fraction(same,removed),'removed_from_zero_true_opposite_fraction':fraction(opp,removed),'removed_from_zero_null_fraction':fraction(null,removed),'selected_true_same_correct_sign_retention':fraction(sel&~wrong,same&available)}
    for key,mask in [('low',snr<1),('middle',(snr>=1)&(snr<2)),('high',snr>=2)]:m['true_same_retention_snr_'+key]=fraction(sel,same&mask&available)
    for key,mask in [('balanced',iv<=.7),('moderate',(iv>.7)&(iv<=.9)),('dominated',iv>.9)]:m['true_same_retention_IV_'+key]=fraction(sel,same&mask&available)
    magnitudes=np.round(np.minimum(np.abs(truth_a),np.abs(truth_b)),12)
    for mag in sorted(set(magnitudes[same])):m['true_same_retention_magnitude_'+str(float(mag))]=fraction(sel,same&(magnitudes==mag)&available)
    den={'all_taxa':len(a),'available':int(available.sum()),'unavailable':int((~available).sum()),'true_same':int((same&available).sum()),'true_opposite':int((opp&available).sum()),'null':int((null&available).sum()),'one_zero':int((one&available).sum()),'selected':int(sel.sum()),'rho0_selected':int((base&available).sum()),'full_selected':int((full&available).sum()),'rho0_removed_by_full':int(removed.sum())}
    m['full_minus_strategy_membership']=fraction(full.astype(int)-sel.astype(int),available)
    m['full_minus_strategy_true_same_retention']=fraction(full.astype(int)-sel.astype(int),same&available)
    m['full_minus_strategy_true_opposite_membership']=fraction(full.astype(int)-sel.astype(int),opp&available)
    for key,mask in [('low',snr<1),('middle',(snr>=1)&(snr<2)),('high',snr>=2)]:den['true_same_snr_'+key]=int((same&mask&available).sum())
    for key,mask in [('balanced',iv<=.7),('moderate',(iv>.7)&(iv<=.9)),('dominated',iv>.9)]:den['true_same_IV_'+key]=int((same&mask&available).sum())
    for mag in sorted(set(magnitudes[same])):den['true_same_magnitude_'+str(float(mag))]=int((same&(magnitudes==mag)&available).sum())
    return m,den
def main():
    scale=sys.argv[1];assert scale in ['pilot','formal']
    out=ROOT/'results'/scale;out.mkdir(parents=True,exist_ok=False)
    cells=[0,1,2,25] if scale=='pilot' else range(26);eval_ids=range(100,102) if scale=='pilot' else range(100,200)
    manifest=json.loads((PREV/'data/formal/MANIFEST.json').read_text())
    priorseal=json.loads((PREV/'provenance/formal_POST_AUDIT_SEAL.json').read_text())['files']
    hashes={};calrows=[];illustrations={};groups=defaultdict(lambda:defaultdict(Moments));denhist=defaultdict(lambda:defaultdict(lambda:defaultdict(int)))
    records=out/'DATASET_RECORDS.jsonl';casefile=out/'TAXON_DECISIONS.jsonl.gz'
    with records.open('x') as f,gzip.open(casefile,'xt',encoding='utf8') as tf:
      for ci in cells:
        cal={m:{y:[] for y in ['A1','A2']} for m in METHODS}
        for rep in range(100):
            r=read(ci,rep,hashes)
            for m in METHODS:
                for y in ['A1','A2']:cal[m][y].append(vec(r['methods'][m],'effect',y)-vec(r['methods'][m],'truth',y))
        rh={};bounds={}
        for m in METHODS:
            e1=np.array(cal[m]['A1']);e2=np.array(cal[m]['A2']);e1-=e1.mean(axis=0);e2-=e2.mean(axis=0)
            den=np.sqrt((e1*e1).sum(axis=0)*(e2*e2).sum(axis=0));rr=np.divide((e1*e2).sum(axis=0),den,out=np.full(100,np.nan),where=den>0)
            valid=np.isfinite(rr)&(np.abs(rr)<1);lo=np.full(100,np.nan);hi=lo.copy()
            lo[valid]=np.maximum(-1,np.tanh(np.arctanh(rr[valid])-1.96/math.sqrt(97)));hi[valid]=np.minimum(1,np.tanh(np.arctanh(rr[valid])+1.96/math.sqrt(97)));valid&=np.isfinite(lo)&np.isfinite(hi)&(lo<=hi);rr[~valid]=np.nan;lo[~valid]=np.nan;hi[~valid]=np.nan;rh[m]=rr;bounds[m]=(lo,hi)
            clean=lambda v:[float(x) if np.isfinite(x) else None for x in v]
            calrows.append({'cell':ci,'parameters':manifest['cells'][ci],'method':m,'n_calibration':100,'rho':clean(rr),'lower':clean(lo),'upper':clean(hi),'availability':valid.tolist(),'halfwidth_on_z':1.96/math.sqrt(97),'reference_status':'simulation-only truth-referenced independent calibration, not known true rho or real-data estimator','approximation':'Fisher z; not simultaneous or guaranteed under nonnormal errors'})
        del cal
        for rep in eval_ids:
          r=read(ci,rep,hashes)
          for m in METHODS:
            obj=r['methods'][m];a=(vec(obj,'effect','A1')+vec(obj,'effect','A2'))/2;b=vec(obj,'effect','B');s1=vec(obj,'se','A1');s2=vec(obj,'se','A2');sb=vec(obj,'se','B');ta=vec(obj,'truth','A1');tb=vec(obj,'truth','B');z=sb**2;d=a-b
            assert all(np.isfinite(v).all() for v in [a,b,s1,s2,sb,ta,tb]) and (np.minimum(np.minimum(s1,s2),sb)>0).all()
            sign=a*b>0;alpha=(s1*s1+s2*s2)/4;snr=np.minimum(np.abs(a)/np.sqrt(alpha),np.abs(b)/sb)
            req={k:vmax(lo,hi,s1,s2,z,d) for k,(lo,hi) in zip(DOMAINS,[(0,0),(rh[m],rh[m]),bounds[m],(-1,1)])}
            reqiv=vmax(-1,1,s1,s2,z,np.zeros(100))
            masks={k:[(sign&(v<=q)) for q in CAPS] for k,v in req.items()}
            tf.write(json.dumps({'cell':ci,'replicate':rep,'method':m,'required_caps':{k:clean(v) for k,v in req.items()},'required_cap_IV_full':reqiv.tolist(),'same_estimated_sign':sign.tolist(),'nonregular_full_endpoint':((s1-s2)==0).tolist(),'truth_A':ta.tolist(),'truth_B':tb.tolist(),'selected':{k:[x.tolist() for x in v] for k,v in masks.items()}},allow_nan=False)+'\n')
            for qi,q in enumerate(CAPS):
              base=masks['zero'][qi];full=masks['full'][qi]
              assert np.all(~full|base)
              for dom in DOMAINS:
                sel=masks[dom][qi];valid=np.isfinite(req[dom]);met,den=metrics(sel,ta,tb,a,base,full,snr,reqiv,valid)
                met.update(mean_required_cap=fraction(req[dom],valid),mean_required_cap_IV_full=float(reqiv.mean()),heterogeneity_assisted_membership=float((sel&(reqiv>q)).mean()) if dom=='full' else None)
                key=(ci,m,q,dom);row={'cell':ci,'replicate':rep,'method':m,'q':q,'domain':dom,'metrics':met,'denominators':den}
                f.write(json.dumps(row)+'\n')
                for k,v in met.items():groups[key][k].add(v)
                for k,v in den.items():denhist[key][k][str(v)]+=1
              for width in HALFWIDTHS:
                domain='width_'+str(width);lo=np.maximum(-1,rh[m]-width);hi=np.minimum(1,rh[m]+width);v=vmax(lo,hi,s1,s2,z,d);sel=sign&(v<=q)
                valid=np.isfinite(v);met={'membership':fraction(sel,valid),'true_same_retention':fraction(sel,(ta*tb>1e-20)&valid),'true_opposite_membership':fraction(sel,(ta*tb < -1e-20)&valid),'mean_required_cap':fraction(v,valid),'mean_domain_width':fraction(hi-lo,valid)};key=(ci,m,q,domain)
                f.write(json.dumps({'cell':ci,'replicate':rep,'method':m,'q':q,'domain':domain,'metrics':met,'denominators':{'all_taxa':100,'available':int(valid.sum()),'unavailable':int((~valid).sum())}})+'\n')
                for k,vv in met.items():groups[key][k].add(vv)
            for typ,test in [('heterogeneity_assisted',(reqiv>.8)),('inverse_variance_sufficient',(reqiv<=.8))]:
                valid=(ta*tb>1e-20)&sign&(req['full']>.6)&(req['full']<=.8)&test&np.isfinite(rh[m])&np.isfinite(bounds[m][0])&np.isfinite(bounds[m][1])
                key=m+'_'+typ
                if key not in illustrations and valid.any():
                    j=int(np.flatnonzero(valid)[0]);illustrations[key]={'cell':ci,'replicate':rep,'method':m,'taxon':j,'parameters':manifest['cells'][ci],'a':float(a[j]),'b':float(b[j]),'s1':float(s1[j]),'s2':float(s2[j]),'sb':float(sb[j]),'truth_A':float(ta[j]),'truth_B':float(tb[j]),'qstar_REML':float(req['full'][j]),'qstar_IV':float(reqiv[j]),'rho_hat':float(rh[m][j]),'calibration_interval':[float(bounds[m][0][j]),float(bounds[m][1][j])],'selection':'first record under fixed ordering satisfying predeclared illustrative condition'}
        print(json.dumps({'completed_cell':ci,'scale':scale}),flush=True)
    for name,h in hashes.items():assert h==priorseal['analysis/formal/'+name],name
    summaries=[{'cell':ci,'parameters':manifest['cells'][ci],'method':m,'q':q,'domain':dom,'metrics':{k:v.result() for k,v in vv.items()},'denominator_histograms':{k:dict(v) for k,v in denhist[(ci,m,q,dom)].items()},'uncertainty':'dataset MCSE conditional on fixed calibration sample'} for (ci,m,q,dom),vv in groups.items()]
    put(out/'SUMMARY.json',summaries);put(out/'CALIBRATION.json',calrows);put(out/'ILLUSTRATIVE_CASES.json',illustrations)
    put(out/'RESULT.json',{'status':'AWAITING_INDEPENDENT_AUDIT','source_root':str(PREV),'source_hashes':hashes,'prior_seal_sha256':sha(PREV/'provenance/formal_POST_AUDIT_SEAL.json'),'calibration_replicates':[0,99],'evaluation_replicates':[100,max(eval_ids)],'scenario_count':len(list(cells)),'summary_rows':len(summaries),'new_fits':False,'new_simulated_datasets':False,'original_native_exit_code':None,'tolerance_classification':0,'truth_ledger_separate':True,'FDR_claim':False})
if __name__=='__main__':main()
