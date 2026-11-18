"""P54 producer. All new arithmetic; no import/call of a sealed analysis."""
from fractions import Fraction as F
from exact_math import rat,s,dec,add,scale,ev,roots,polys,numeric,constants

def val(x):
    try:return rat(x) if x is not None else None
    except (ValueError,TypeError,ZeroDivisionError):return None
def prevalence(records):
    if not records:return 'UNKNOWN'
    if any(r!=records[0] for r in records):return 'UNKNOWN'
    r=records[0]
    try:
        n=int(r['retained_samples']);k=int(r['samples_ge_10_reads'])
        if n<=0 or k<0 or k>n:return 'UNKNOWN'
        p=k*10>=n
        literal=r['passes_frozen_filter']
        if literal not in ('TRUE','FALSE'):return 'UNKNOWN'
        if (literal=='TRUE')!=p:return 'UNKNOWN'
        return 'PASS' if p else 'FAIL'
    except (KeyError,ValueError,TypeError):return 'UNKNOWN'
def conj(items):
    if 'FAIL' in items:return 'FAIL'
    if any(x!='PASS' for x in items):return 'UNKNOWN'
    return 'PASS'
def base_gates(a,b,p):
    if a==0 or b==0:g1='FAIL'
    elif a is None or b is None:g1='UNKNOWN'
    else:g1='PASS' if a*b>0 else 'FAIL'
    return g1,conj([prevalence(p.get('A',[])),prevalence(p.get('B',[]))])
def extract(entry):
    ar=entry['experiments'].get('A');br=entry['experiments'].get('B')
    a=val(ar.get('theta')) if ar else None;b=val(br.get('theta')) if br else None
    alpha=val(ar.get('variance_intercept')) if ar else None;beta=val(ar.get('variance_slope')) if ar else None
    z=None
    if br and len(br['years'])==1:
        se=val(br['years'][0].get('se_literal'))
        if se is not None and se>0:z=se*se
    return a,b,alpha,beta,z
def outclass(cells,key,degenerate):
    values=[c['gates'][key] if key!='joint' else c['joint'] for c in cells if c['state']=='REGULAR']
    if not values or any(x not in ('PASS','FAIL') for x in values):return 'NOT_ASSESSABLE'
    if all(x=='PASS' for x in values):return 'ALL_REGULAR_RHO_WITH_DEGENERATE_ENDPOINT' if degenerate else 'ALL_RHO'
    if all(x=='FAIL' for x in values):return 'NONE'
    return 'SOME_RHO'
def domains(cells):
    return {key:{state:[dict(kind=c['kind'],lo=c['lo'],hi=c['hi'],closed=c['closed']) for c in cells
               if c['state']=='REGULAR' and (c['joint'] if key=='joint' else c['gates'][key])==state]
          for state in ['PASS','FAIL','UNKNOWN','NOT_APPLICABLE_NO_COMMON_DIRECTION']}
          for key in ['G1','G2','G3','G4','G5','joint']}
def limit_zero(a,b,z):
    return dict(meaning='ONE_SIDED_LIMIT_NOT_REGULAR_ESTIMATE',
                numeric=numeric(a,b,z,F(0),endpoint_limit=True))
def evaluate(entry):
    a,b,alpha,beta,z=extract(entry);g1,g4=base_gates(a,b,entry['prevalence_records'])
    result=dict(family=entry['family'],input_source='BOUND_INPUT_EXACT_FAMILY_KEY',
        primary_experiments=['A','B'],source_diagnostics={e:(r['diagnostic_stratum'] if r else 'ABSENT')
                for e,r in entry['experiments'].items()},
        formal_P1_status='NOT_RELEASED',P1_count=None,signature_weight=None,Bonn_release=False,
        P=None,FDR=None,estimand='MODEL_CONDITIONAL_AB_NATURAL_LOG_ASSOCIATION')
    result['C_auxiliary']=dict(role='AUXILIARY_NOT_POOLED_NOT_AN_AB_GATE',
        presence='PRESENT' if entry['experiments'].get('C') else 'ABSENT',
        source_record=entry['experiments'].get('C'),
        provenance='Literal bound C record, including full affine rho_C precision, year estimates, diagnostics and missing statuses; no new C computation')
    available=all(v is not None for v in [a,b,alpha,beta,z])
    if available and (z<=0 or beta<0 or alpha-beta<0):raise ValueError('Invalid saved variance geometry; no repair/floor')
    if not available or (alpha==0 and beta==0):
        gs={'G1':g1,'G2':'UNKNOWN','G3':('PASS' if g1=='PASS' else 'NOT_APPLICABLE_NO_COMMON_DIRECTION'),'G4':g4,'G5':'UNKNOWN'}
        cell=dict(kind='INTERVAL',lo='-1/1',hi='1/1',closed=[True,True],state='REGULAR',gates=gs,joint=conj(list(gs.values())),numeric=None,
                  missing_reason='SAVED_EFFECT_OR_POSITIVE_PRECISION_UNAVAILABLE')
        result.update(precision_supported=False,cells=[cell],certificates=[],domain=domains([cell]),
            classes={k:outclass([cell],k,False) for k in [*gs,'joint']},degenerate_endpoints=[],
            input={k:(s(v) if v is not None else None) for k,v in zip(['a','b','alpha','beta','z'],[a,b,alpha,beta,z])})
        return result
    result['precision_supported']=True
    D=a-b;cuts={F(-1),F(1)}
    if beta:
        for r in [(D*D-z-alpha)/beta,-alpha/beta]:
            if -1<=r<=1:cuts.add(r)
    initial=sorted(cuts);certs=[]
    for lo,hi in zip(initial,initial[1:]):
        mid=(lo+hi)/2;x=alpha+beta*mid
        if x<=0:continue
        branch='POSITIVE_TAU' if D*D>x+z else 'ZERO_TAU'
        for name,pair in polys(a,b,z,alpha,beta,branch).items():
            for bound,p in enumerate(pair):
                rb,proof=roots(p,lo,hi)
                for rr in rb:
                    cuts.update([F(rr['lo']),F(rr['hi'])])
                certs.append(dict(name=name,bound=bound,branch=branch,domain=[s(lo),s(hi)],
                    coefficients=[s(v) for v in p],roots=rb,proof=proof))
    cuts=sorted(cuts);atoms=[]
    for i,r in enumerate(cuts):
        atoms.append(('POINT',r,r))
        if i+1<len(cuts):atoms.append(('INTERVAL',r,cuts[i+1]))
    cells=[];degen=[]
    for kind,lo,hi in atoms:
        rho=(lo+hi)/2;x=alpha+beta*rho
        cell=dict(kind=kind,lo=s(lo),hi=s(hi),closed=[kind=='POINT',kind=='POINT'],sample_rho=s(rho))
        if x==0:
            assert kind=='POINT'
            cell.update(state='DEGENERATE_ENDPOINT',gates=None,joint=None,numeric=None,limit=limit_zero(a,b,z));degen.append(s(rho))
            cells.append(cell);continue
        t=max(F(0),(D*D-x-z)/2);wa=1/(x+t);wb=1/(z+t);p=wa/(wa+wb);mu=p*a+(1-p)*b
        gs={'G1':g1,'G2':('PASS' if g1=='PASS' and mu*a>0 else ('NOT_APPLICABLE_NO_COMMON_DIRECTION' if g1=='FAIL' else 'UNKNOWN')),
            'G3':('PASS' if g1=='PASS' else 'NOT_APPLICABLE_NO_COMMON_DIRECTION'),
            'G4':g4,'G5':('PASS' if F(3,10)<=p<=F(7,10) else 'FAIL')}
        branch='POSITIVE_TAU' if t>0 else 'ZERO_TAU'
        interval_states={}
        for name in ['mHK_zero','HK_zero','PI_zero','LOO_A_zero','LOO_B_zero']:
            pp=polys(a,b,z,alpha,beta,branch)[name];values=[ev(q,rho) for q in pp]
            uncertain=False
            for cert in certs:
                cl,ch=map(F,cert['domain'])
                if cert['name']!=name or cert['branch']!=branch or not cl<=rho<=ch:continue
                for rr in cert['roots']:
                    rl,rh=F(rr['lo']),F(rr['hi'])
                    if rl<rh and rl<=rho<=rh:uncertain=True
            center=a if name=='LOO_A_zero' else b if name=='LOO_B_zero' else mu
            if not uncertain and min(values)>0:st='POSITIVE_EXCLUDES_ZERO' if center>0 else 'NEGATIVE_EXCLUDES_ZERO'
            elif not uncertain and max(values)<0:st='CONTAINS_ZERO'
            elif values[0]==values[1]==0:st='TOUCHES_ZERO'
            else:st='UNKNOWN_ZERO_BOUNDARY'
            interval_states[name]=st
        reverse=[]
        for name in ['LOO_A_zero','LOO_B_zero']:
            st=interval_states[name]
            reverse.append('UNDEFINED_POOLED_ZERO' if mu==0 else ('UNKNOWN' if st=='UNKNOWN_ZERO_BOUNDARY'
                else bool((mu>0 and st=='NEGATIVE_EXCLUDES_ZERO') or (mu<0 and st=='POSITIVE_EXCLUDES_ZERO'))))
        q=D*D/(x+z);i2=F(0) if q==0 else max(F(0),100*(q-1)/q)
        csrc=entry['experiments'].get('C');ctheta=val(csrc.get('theta')) if csrc else None
        cdir='UNAVAILABLE' if ctheta is None else 'ZERO_REFERENCE' if mu==0 or ctheta==0 else 'SAME_SIGN' if ctheta*mu>0 else 'OPPOSITE_SIGN'
        cell.update(state='REGULAR',branch=branch,gates=gs,joint=conj(list(gs.values())),numeric=numeric(a,b,z,x),
            interval_zero_states=interval_states,CR_components=dict(I2_75=i2>=75,opposing_experiment_interval=reverse,
                PI_reverse_extent_status=interval_states['PI_zero'],complete_C_R='NOT_RELEASED'),
            C_auxiliary_sign_relation=cdir)
        cells.append(cell)
    result.update(input={k:s(v) for k,v in zip(['a','b','alpha','beta','z'],[a,b,alpha,beta,z])},
        cells=cells,certificates=certs,domain=domains(cells),
        classes={k:outclass(cells,k,bool(degen)) for k in ['G1','G2','G3','G4','G5','joint']},degenerate_endpoints=degen,
        C_auxiliary_source='C_auxiliary.source_record is the complete literal bound C record')
    return result
