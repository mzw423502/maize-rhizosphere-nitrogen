"""Pure D3 algorithms. Training API cannot receive test data or labels."""
import math
from collections import defaultdict

def rank_z(values):
    n=len(values)
    if n<2:raise ValueError('N_LT_2')
    order=sorted(range(n),key=values.__getitem__)
    if values[order[0]]==values[order[-1]]:return [0.5]*n,[0.0]*n
    ranks=[0.0]*n;i=0
    while i<n:
        j=i+1
        while j<n and values[order[j]]==values[order[i]]:j+=1
        scaled=((i+j-1)/2)/(n-1)
        for k in order[i:j]:ranks[k]=scaled
        i=j
    mean=math.fsum(ranks)/n
    sd=math.sqrt(math.fsum((x-mean)**2 for x in ranks)/(n-1))
    return ranks,[(x-mean)/sd for x in ranks]

def structure(meta):
    if len(meta)<2:raise ValueError('N_LT_2')
    units={};pairs=defaultdict(dict)
    for i,r in enumerate(meta):
        if r['arm'] not in ('LN','HN'):raise ValueError('UNKNOWN_N_ARM')
        u=units.setdefault(r['unit'],{'pair':r['pair'],'arm':r['arm'],'cells':defaultdict(list)})
        if (u['pair'],u['arm'])!=(r['pair'],r['arm']):raise ValueError('UNIT_METADATA_CONFLICT')
        u['cells'][r['cell']].append(i)
    for name,u in units.items():
        if u['arm'] in pairs[u['pair']]:raise ValueError('MULTIPLE_UNITS_PER_PAIR_ARM')
        pairs[u['pair']][u['arm']]=name
    if any(set(p)!={'LN','HN'} for p in pairs.values()):raise ValueError('STRUCTURALLY_UNESTIMABLE_MISSING_ARM')
    return units,dict(pairs)

def summarize(values,meta):
    if len(values)!=len(meta):raise ValueError('IDENTITY_LENGTH_MISMATCH')
    units,pairs=structure(meta);unit_rows=[];means={}
    for name,u in sorted(units.items()):
        cellmeans=[math.fsum(values[i] for i in ix)/len(ix) for _,ix in sorted(u['cells'].items())]
        m=math.fsum(cellmeans)/len(cellmeans);means[name]=m
        unit_rows.append({'unit':name,'pair':u['pair'],'arm':u['arm'],'unit_mean':m,'n_samples':sum(map(len,u['cells'].values())),'n_cells':len(cellmeans)})
    pair_rows=[{'pair':p,'LN_unit':u['LN'],'HN_unit':u['HN'],'LN_mean':means[u['LN']],'HN_mean':means[u['HN']],
                'contrast_LN_minus_HN':means[u['LN']]-means[u['HN']]} for p,u in sorted(pairs.items())]
    cc=[r['contrast_LN_minus_HN'] for r in pair_rows];mean=math.fsum(cc)/len(cc)
    return {'units':unit_rows,'pairs':pair_rows,'equal_pair_mean':mean,'pair_min':min(cc),'pair_max':max(cc),
            'pair_sample_sd':math.sqrt(math.fsum((x-mean)**2 for x in cc)/(len(cc)-1)) if len(cc)>1 else None,
            'biological_SE':None,'biological_CI':None,'biological_P':None,'biological_FDR':None}

def train(train_cohorts,allowed):
    """train_cohorts=list of {id,counts,metadata}; no test arguments."""
    for c in train_cohorts:structure(c['metadata'])
    labels=list(train_cohorts[0]['counts'])
    if any(list(c['counts'])!=labels for c in train_cohorts):raise ValueError('TRAIN_VOCABULARY_MISMATCH')
    n=sum(len(c['metadata']) for c in train_cohorts);rows=[]
    for label in labels:
        decision=allowed[label]
        ge=sum(sum(x>=10 for x in c['counts'][label]) for c in train_cohorts)
        eligible=decision=='INCLUDED_BY_REVIEW_RULE' and 10*ge>=n
        row={'taxon':label,'label_decision':decision,'train_samples_ge10':ge,'train_n':n,'eligible':eligible,
             'cohort_contrasts':{},'training_contrast':None,'weight':None,'selected':False}
        if eligible:
            for c in train_cohorts:
                _,z=rank_z(c['counts'][label]);row['cohort_contrasts'][c['id']]=summarize(z,c['metadata'])['equal_pair_mean']
            row['training_contrast']=math.fsum(row['cohort_contrasts'].values())/len(train_cohorts)
        rows.append(row)
    candidates=[r for r in rows if r['eligible'] and math.isfinite(r['training_contrast']) and r['training_contrast']!=0]
    candidates.sort(key=lambda r:(-abs(r['training_contrast']),r['taxon'].encode('utf-8')))
    selected=candidates[:10]
    if len(selected)<3:return {'status':'NO_SCORE_FEWER_THAN_3','candidates':rows,'selected':[],'train_n':n}
    denom=math.fsum(abs(r['training_contrast']) for r in selected)
    for r in selected:r['selected']=True;r['weight']=r['training_contrast']/denom
    return {'status':'FROZEN_SCORE','candidates':rows,'selected':[{'taxon':r['taxon'],'contrast':r['training_contrast'],'weight':r['weight']} for r in selected],'train_n':n}

def score(freeze,counts):
    """Only frozen weights and unlabeled test counts. No labels/units/target filters."""
    if freeze['status']!='FROZEN_SCORE':return {'status':freeze['status'],'scores':None,'components':None}
    missing=[r['taxon'] for r in freeze['selected'] if r['taxon'] not in counts]
    if missing:return {'status':'INPUT_GATE_UNAVAILABLE_MISSING_SELECTED_COLUMN','missing':missing,'scores':None,'components':None}
    zz={};constant=[]
    for r in freeze['selected']:
        label=r['taxon'];_,z=rank_z(counts[label]);zz[label]=z
        if not any(z):constant.append(label)
    n=len(next(iter(zz.values())))
    assert all(len(z)==n for z in zz.values())
    result=[math.fsum(r['weight']*zz[r['taxon']][i] for r in freeze['selected']) for i in range(n)]
    return {'status':'SCORED_WITHOUT_TEST_LABELS','scores':result,'components':zz,'constant_selected_test_families':constant,'weights_redistributed':False}
