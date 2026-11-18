"""D2 saved-component arithmetic. No model, RNG, old helper or network imports."""
import csv
import ctypes
import hashlib
import json
import math
import os
import platform
import shutil
import subprocess
import sys
import time
import traceback
from collections import Counter
from datetime import datetime, timezone
from decimal import Decimal, localcontext
from pathlib import Path

ROOT = Path('.')
OUT = Path(__file__).resolve().parent
P51 = ROOT / '15_FINAL_AUDIT/PHASE5_P51_ENTRY_CONTRACT_20261010_V1'
COHORTS = {'A':['A_Y2018','A_Y2019'], 'B':['B_2022'], 'C':['C_2017','C_2018']}
GIB = 1024**3

def now(): return datetime.now(timezone.utc).isoformat()
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p): return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def save(p, obj):
    with Path(p).open('x', encoding='utf-8', newline='\n') as f:
        json.dump(obj, f, ensure_ascii=False, indent=2, allow_nan=False)
        f.write('\n')
def git(*args):
    return subprocess.check_output(['git','-C',str(ROOT),*args], creationflags=0x08000000)
def value(s):
    if s is None: return None
    try: d=Decimal(str(s))
    except Exception: return None
    return d if d.is_finite() else None
def num(d):
    if d is None: return None
    f=float(d)
    if not math.isfinite(f): raise ArithmeticError('Nonfinite derived binary64 output')
    return f
def memory():
    class PMC(ctypes.Structure):
        _fields_=[('cb',ctypes.c_ulong),('PageFaultCount',ctypes.c_ulong)]+[(n,ctypes.c_size_t) for n in ['PeakWorkingSetSize','WorkingSetSize','QuotaPeakPagedPoolUsage','QuotaPagedPoolUsage','QuotaPeakNonPagedPoolUsage','QuotaNonPagedPoolUsage','PagefileUsage','PeakPagefileUsage','PrivateUsage']]
    x=PMC();x.cb=ctypes.sizeof(x)
    k=ctypes.WinDLL('kernel32');k.GetCurrentProcess.restype=ctypes.c_void_p
    p=ctypes.WinDLL('psapi');p.GetProcessMemoryInfo.argtypes=[ctypes.c_void_p,ctypes.c_void_p,ctypes.c_ulong]
    if not p.GetProcessMemoryInfo(k.GetCurrentProcess(),ctypes.byref(x),x.cb): raise ctypes.WinError()
    return {'working_set_bytes':x.WorkingSetSize,'peak_working_set_bytes':x.PeakWorkingSetSize,'private_bytes':x.PrivateUsage,'peak_pagefile_bytes':x.PeakPagefileUsage}
def resources():
    class MS(ctypes.Structure):
        _fields_=[('length',ctypes.c_ulong),('load',ctypes.c_ulong)]+[(n,ctypes.c_ulonglong) for n in ['total_physical','available_physical','total_pagefile','available_pagefile','total_virtual','available_virtual','extended']]
    x=MS();x.length=ctypes.sizeof(x)
    if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(x)): raise ctypes.WinError()
    return {'utc':now(),'available_physical_bytes':x.available_physical,'commit_headroom_bytes':x.available_pagefile,'C_free_bytes':shutil.disk_usage('C:/').free}

def derive(years):
    """Decimal80 computations; strict threshold decisions use decimal source literals."""
    with localcontext() as ctx:
        ctx.prec=80
        b=[value(y.get('lfc_literal')) for y in years]
        s=[value(y.get('se_literal')) for y in years]
        out={'theta':None,'point_status':'UNAVAILABLE_MISSING_OR_NONFINITE_YEAR_LFC',
             'precision_status':'UNAVAILABLE_POINT','point_sign':None,'variance_intercept':None,
             'variance_slope':None,'se_min':None,'se_max':None,'critical_rho':None,
             'rho_zero_exclusion_class':None,'rho_interval_low':None,'rho_interval_high':None,
             'rho_low_inclusive':None,'rho_high_inclusive':None,
             'degenerate_minus_one':None,'contour_min_low':None,'contour_min_high':None,
             'contour_max_low':None,'contour_max_high':None,'precision_ratio_min_se':None,
             'precision_ratio_max_se':None,'B_two_SE_excludes_zero':None}
        if not b or any(x is None for x in b): return out
        if len(b) not in (1,2): raise ValueError('Only one or two saved years')
        theta=sum(b)/len(b)
        out.update(theta=num(theta),point_status='AVAILABLE',point_sign='POSITIVE' if theta>0 else 'NEGATIVE' if theta<0 else 'ZERO')
        if any(x is None or x<=0 for x in s):
            out['precision_status']='UNAVAILABLE_MISSING_NONFINITE_OR_NONPOSITIVE_NATIVE_SE'
            return out
        if len(b)==1:
            lo=hi=s[0]
            out.update(precision_status='SINGLE_YEAR_NATIVE',rho_zero_exclusion_class='NOT_APPLICABLE_SINGLE_YEAR',B_two_SE_excludes_zero=abs(theta)>2*lo)
        else:
            a=(s[0]**2+s[1]**2)/4
            k=s[0]*s[1]/2
            lo=abs(s[0]-s[1])/2;hi=(s[0]+s[1])/2
            critical=(theta**2-s[0]**2-s[1]**2)/(2*s[0]*s[1])
            # Endpoint comparisons avoid rounding the critical value for boundary logic.
            min_excl=theta!=0 and theta**2>(s[0]-s[1])**2
            max_excl=theta!=0 and theta**2>(s[0]+s[1])**2
            if not min_excl:
                cls='NO_RHO'; il=ih=None;lin=hin=False
            elif max_excl:
                cls='ALL_RHO';il=Decimal(-1);ih=Decimal(1);lin=hin=True
            elif theta**2==(s[0]+s[1])**2:
                cls='ALL_EXCEPT_PLUS_ONE';il=Decimal(-1);ih=Decimal(1);lin=True;hin=False
            else:
                cls='PARTIAL_RHO';il=Decimal(-1);ih=critical;lin=True;hin=False
            out.update(precision_status='FULL_RHO_CONDITIONAL',variance_intercept=num(a),variance_slope=num(k),
                       critical_rho=num(critical),rho_zero_exclusion_class=cls,
                       rho_interval_low=num(il),rho_interval_high=num(ih),rho_low_inclusive=lin,rho_high_inclusive=hin,
                       degenerate_minus_one=lo==0)
        out.update(se_min=num(lo),se_max=num(hi),contour_min_low=num(theta-2*lo),contour_min_high=num(theta+2*lo),
                   contour_max_low=num(theta-2*hi),contour_max_high=num(theta+2*hi),
                   precision_ratio_min_se=num(abs(theta)/lo) if lo>0 else None,precision_ratio_max_se=num(abs(theta)/hi))
        return out

def assemble(experiment, family, years):
    result=derive(years)
    clean=all(y['diagnostic_state']=='NO_CAPTURED_FINAL_FLAG' for y in years)
    return {'experiment':experiment,'family':family,'role':'AUXILIARY' if experiment=='C' else 'PRIMARY_DISCOVERY',
            'year_count':len(years),'diagnostic_stratum':'ALL_YEARS_NO_CAPTURED_FLAG' if clean else 'ORIGINAL_DIAGNOSTICS_RETAINED',
            'years':years,**result,'bio_CI':None,'bio_P':None,'bio_FDR':None,'meta_weight':None,'P1':None,
            'contour_definition':'theta plus/minus two native model standard errors; conditional, not biological95CI'}

def source_records():
    bindings=read(P51/'SAVED_COMPONENT_BINDINGS.json')['rows']
    p51src=read(P51/'SOURCE_BINDINGS.json')['bindings']
    joint=[b for b in p51src if b['path'].endswith('ALL_FAMILY_DIRECTION_COMPARISON.json')][0]
    checks=[]
    for item in bindings+[joint]:
        p=ROOT/item['path']
        if sha(p)!=item['sha256']: raise ValueError('Source SHA mismatch '+str(p))
        blob=git('hash-object','--no-filters',str(p)).decode().strip()
        if blob!=git('rev-parse','HEAD:'+item['path']).decode().strip(): raise ValueError('Source raw/HEAD mismatch')
        checks.append({'path':item['path'],'sha256':sha(p),'blob':blob})
    sealed=read(ROOT/joint['path'])['rows']
    if len(sealed)!=835: raise ValueError('835 universe changed')
    table={}
    for item in bindings:
        with (ROOT/item['path']).open(encoding='utf-8-sig',newline='') as f:
            rows=list(csv.DictReader(f,delimiter='\t'))
        m={}
        for i,r in enumerate(rows,2):
            if r['coefficient']!='NLN' or r['family'] in m: raise ValueError('Source row identity')
            m[r['family']]={'lfc_literal':r['lfc'],'se_literal':r['se'],'dof_literal':r['dof'],
                              'source_path':item['path'],'source_line':i,'source_sha256':item['sha256']}
        table[item['cohort']]=m
    with (P51/'EXPERIMENT_ENTRY.tsv').open(encoding='utf-8-sig',newline='') as f: entries=list(csv.DictReader(f,delimiter='\t'))
    if len(entries)!=835: raise ValueError('Entry universe changed')
    entry={(r['cohort'],r['family_label']):r for r in entries}
    if len(entry)!=835: raise ValueError('Duplicate entry key')
    grouped={}
    for j in sealed:
        c,fam=j['cohort'],j['family'];e=c[0];r=entry[(c,fam)]
        if not r['effect_scale'].startswith('ANCOM native natural-log coefficient NLN'): raise ValueError('Scale')
        native=table[c].get(fam)
        expected=j['ANCOM_estimate_status']=='ESTIMATE_PRESENT_NOT_REQUALIFIED_FOR_INFERENCE'
        if (native is not None)!=expected: raise ValueError('Native/835 existence mismatch')
        year={'cohort':c,'family':fam,'diagnostic_state':r['diagnostic_state'],
              'diagnostics_literal':j['ANCOM_final_diagnostic_literal'],
              'input_support_literal':j['ANCOM_input_support_literal'],
              'source_status':j['ANCOM_estimate_status'],'N_contrast':r['N_contrast'],
              'reference_context':r['reference_compatibility'],'unit_count':int(r['unit_count']),'pair_count':int(r['pair_count']),
              'lfc_literal':None,'se_literal':None,'dof_literal':None}
        if native: year.update(native)
        grouped.setdefault((e,fam),{})[c]=year
    keys=sorted(grouped,key=lambda k:(k[0],k[1].encode('utf-8')))
    if Counter(k[0] for k in keys)!=Counter(A=197,B=119,C=161): raise ValueError('477 experiment universe')
    records=[assemble(e,f,[grouped[(e,f)][c] for c in COHORTS[e]]) for e,f in keys]
    return records,checks

def actual():
    from runtime_guard import owner, check
    with owner():
        snapshot=check(ignore_own_dispatch_parent=True)
        if not snapshot['effective_resource_gate_passed']:
            save(OUT/'CALCULATOR_PRESTART_HOLD.json',{'utc':now(),'snapshot':snapshot,'actual_started':False})
            raise SystemExit(75)
        _actual_after_guard(snapshot)

def _actual_after_guard(snapshot):
    folder=OUT/'actual'
    if folder.exists(): raise RuntimeError('Actual namespace exists; never restart')
    auth=read(OUT/'AUTHORIZATION.json')
    if sha(P51/'NEXT_STATISTICAL_CONTRACT.json')!=auth['contract_sha256']: raise ValueError('Contract changed')
    ready=read(OUT/'PREPARATION_REVIEW.json')
    if ready['actual_implementation_review']!='PASS': raise ValueError('No implementation release')
    if not (OUT/'benchmark'/'INDEPENDENT_AUDIT.json').exists(): raise ValueError('Benchmark audit missing')
    if read(OUT/'benchmark'/'INDEPENDENT_AUDIT.json')['differences']: raise ValueError('Benchmark audit failed')
    r=resources();gate=snapshot['effective_resource_gate_bytes']
    if r['available_physical_bytes']<gate['available_physical'] or r['commit_headroom_bytes']<gate['commit_headroom'] or r['C_free_bytes']<gate['C_disk']:
        save(OUT/'CALCULATOR_PRESTART_HOLD.json',{'utc':now(),'resources':r,'actual_started':False})
        raise SystemExit(75)
    folder.mkdir()
    save(folder/'START.json',{'utc':now(),'pid':os.getpid(),'resources':r,'contract_sha256':auth['contract_sha256'],
                             'source_code_sha256':sha(__file__),'runtime':sys.version,'threads':1,'current_checks':snapshot})
    started=time.perf_counter()
    try:
        rows,checks=source_records()
        save(folder/'RESULT.json',{'rows':rows,'input_checks':checks,'scope':'D2 conditional saved-component calculation'})
        summary={e:{'universe':sum(x['experiment']==e for x in rows),'point_status':dict(Counter(x['point_status'] for x in rows if x['experiment']==e)),
                    'precision_status':dict(Counter(x['precision_status'] for x in rows if x['experiment']==e)),
                    'contour_class':dict(Counter(x['rho_zero_exclusion_class'] for x in rows if x['experiment']==e))} for e in COHORTS}
        save(folder/'COMPLETION.json',{'utc':now(),'exit_code':0,'elapsed_seconds':time.perf_counter()-started,'memory':memory(),
                                     'records':len(rows),'summary':summary,'independent_audit':'NOT_YET_RUN','review':'NOT_YET_RUN'})
    except Exception:
        save(folder/'FAILURE.json',{'utc':now(),'traceback':traceback.format_exc(),'memory':memory(),'automatic_retry':False})
        raise

if __name__=='__main__':
    if sys.argv[1:]!=['actual']: raise SystemExit('Only explicit actual entry allowed')
    actual()
