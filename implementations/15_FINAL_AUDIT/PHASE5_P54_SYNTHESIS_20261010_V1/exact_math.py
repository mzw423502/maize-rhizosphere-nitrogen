"""Producer: exact rational coefficients/root enclosures; 60+ digit constants."""
from fractions import Fraction as F
from math import isqrt
from decimal import Decimal, localcontext
from functools import lru_cache

def rat(x): return F(str(x))
def s(q): return str(q.numerator)+'/'+str(q.denominator)
def dec(q,d=55):
    with localcontext() as c:
        c.prec=d
        return str(Decimal(q.numerator)/Decimal(q.denominator))
def bounds(q,d=75):
    k=10**d;n=q.numerator*k//q.denominator
    return F(n,k),F(n+1,k)
def sqrt_bounds(q,d=85):
    assert q>=0
    if q==0:return F(0),F(0)
    an=isqrt(q.numerator);ad=isqrt(q.denominator)
    if an*an==q.numerator and ad*ad==q.denominator:return F(an,ad),F(an,ad)
    # Scale the grid with the operand, so small variances retain d significant digits.
    exponent=len(str(q.numerator))-len(str(q.denominator))
    power=F(10**exponent) if exponent>=0 else F(1,10**(-exponent))
    if q<power:exponent-=1
    grid_digits=max(0,d-exponent//2)
    k=10**grid_digits;n=isqrt(q.numerator*k*k//q.denominator)
    return F(n,k),F(n+1,k)
def atan_bounds(inv):
    x=F(1,inv);total=F(0);power=x;n=0
    while True:
        term=power/(2*n+1);total+=term if n%2==0 else -term
        nxt=power*x*x/(2*n+3)
        if nxt<F(1,10**95):
            other=total+(-nxt if n%2==0 else nxt)
            return min(total,other),max(total,other)
        power*=x*x;n+=1
def sincos_bounds(x,cos=False):
    total=F(0);term=F(1) if cos else x;n=0
    while True:
        total+=term if n%2==0 else -term
        den=(2*n+1)*(2*n+2) if cos else (2*n+2)*(2*n+3)
        nxt=term*x*x/den
        if nxt<F(1,10**90):
            other=total+(-nxt if n%2==0 else nxt)
            return min(total,other),max(total,other)
        term=nxt;n+=1
def normal_integral_bounds(z):
    # Integral exp(-t^2/2) dt from 0 to z; z near 1.96, decreasing alternating terms.
    assert F(19,10)<z<F(2)
    total=F(0);power=z;fact=1;twop=1;n=0
    while True:
        term=power/(twop*fact*(2*n+1));total+=term if n%2==0 else -term
        nxt=power*z*z/(twop*2*fact*(n+1)*(2*n+3))
        if nxt<F(1,10**85):
            other=total+(-nxt if n%2==0 else nxt)
            return min(total,other),max(total,other)
        power*=z*z;fact*=n+1;twop*=2;n+=1
@lru_cache(None)
def constants():
    a,b=atan_bounds(5);c,d=atan_bounds(239)
    pl,pu=16*a-4*d,16*b-4*c
    pl=bounds(pl,80)[0];pu=bounds(pu,80)[1]
    xl,xu=pl/40,pu/40
    sl,su=sincos_bounds(xl);cl,cu=sincos_bounds(xl,True)
    tl,tu=sincos_bounds(xu);dl,du=sincos_bounds(xu,True)
    Tlo=bounds((dl/tu)**2,65)[0];Thi=bounds((cu/sl)**2,65)[1]
    # Approximate z by Decimal bisection; certify enclosing rationals below.
    with localcontext() as ctx:
        ctx.prec=95;pi=Decimal(pl.numerator)/Decimal(pl.denominator);scale=(2*pi).sqrt()
        low=Decimal('1.95');high=Decimal('1.97')
        for _ in range(290):
            z=(low+high)/2;term=z;total=z
            for n in range(1,180):
                term*= -z*z/Decimal(2*n)
                total+=term/Decimal(2*n+1)
            if Decimal('.5')+total/scale<Decimal('.975'):low=z
            else:high=z
        zl=bounds(F(str(low)),65)[0];zu=bounds(F(str(high)),65)[1]
    kl,ku=sqrt_bounds(2*pl)[0],sqrt_bounds(2*pu)[1]
    il,iu=normal_integral_bounds(zl)
    jl,ju=normal_integral_bounds(zu)
    assert F(1,2)+iu/kl<F(39,40)<F(1,2)+jl/ku
    return dict(T2=(Tlo,Thi),Z2=(zl*zl,zu*zu),pi=(pl,pu),z=(zl,zu),
        proof='Machin alternating arctan; alternating sin/cos and Gaussian integral; integer square-root brackets; rational CDF endpoint signs')
def add(a,b):
    return [ (a[i] if i<len(a) else F(0))+(b[i] if i<len(b) else F(0)) for i in range(max(len(a),len(b))) ]
def mul(a,b):
    o=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):o[i+j]+=x*y
    return o
def scale(a,k):return [x*k for x in a]
def ev(a,x):
    v=F(0)
    for c in reversed(a):v=v*x+c
    return v
def roots(poly,lo,hi):
    p=list(poly)
    while len(p)>1 and p[-1]==0:p.pop()
    assert len(p)<=3
    if len(p)==1:return [],('IDENTICALLY_ZERO' if p[0]==0 else 'NONZERO_CONSTANT')
    if len(p)==2:
        x=-p[0]/p[1]
        return ([dict(lo=s(x),hi=s(x),multiplicity=1)] if lo<=x<=hi else []),'LINEAR_EXACT'
    c,b,a=p;dis=b*b-4*a*c
    if dis<0:return [],'NEGATIVE_DISCRIMINANT'
    if dis==0:
        x=-b/(2*a)
        return ([dict(lo=s(x),hi=s(x),multiplicity=2)] if lo<=x<=hi else []),'DOUBLE_ROOT_EXACT'
    dl,du=sqrt_bounds(dis,110)
    out=[]
    for sign in [-1,1]:
        u=(-b+sign*dl)/(2*a);v=(-b+sign*du)/(2*a);left,right=min(u,v),max(u,v)
        if right<lo or left>hi:continue
        out.append(dict(lo=s(max(left,lo)),hi=s(min(right,hi)),multiplicity=1,
             sqrt_discriminant_bracket=[s(dl),s(du)],unclipped_root_bracket=[s(left),s(right)]))
    return out,'QUADRATIC_DISCRIMINANT_CERTIFIED'
def polys(a,b,z,alpha,beta,branch):
    X=[alpha,beta];S=add(X,[z]);D=a-b;cst=constants()
    if branch=='ZERO_TAU':
        U=add(scale(X,b),[z*a]);mu=U
        common=mul(mul(X,[z]),S);hk=scale(X,D*D*z)
        cn=mul(U,U);pin=common
    else:
        U=add(scale(X,-1),[D*(a+b)+z]);mu=U
        common=add([D**4],scale(mul(add(X,[-z]),add(X,[-z])),-1))
        hk=common;cn=mul(U,U)
        pin=add(common,scale(add(scale(X,-1),[D*D-z]),2*D*D))
    out={'mu_zero':(mu,mu),'I2_75':(add(X,[z-D*D/4]),)*2}
    for name,v in [('mHK_zero',common),('HK_zero',hk),('PI_zero',pin)]:
        out[name]=(add(cn,scale(v,-cst['T2'][0])),add(cn,scale(v,-cst['T2'][1])))
    out['LOO_A_zero']=(add([a*a],scale(X,-cst['Z2'][0])),add([a*a],scale(X,-cst['Z2'][1])))
    out['LOO_B_zero']=([b*b-cst['Z2'][0]*z],[b*b-cst['Z2'][1]*z])
    if branch=='ZERO_TAU':
        out['weight_A_30']=(add(scale(S,-F(3,10)),[z]),)*2
        out['weight_A_70']=(add(scale(S,-F(7,10)),[z]),)*2
    else:
        out['weight_A_30']=(add(scale(X,-1),[z+F(2,5)*D*D]),)*2
        out['weight_A_70']=(add(scale(X,-1),[z-F(2,5)*D*D]),)*2
    return out
def numeric(a,b,z,x,endpoint_limit=False):
    assert z>0 and (x>0 or (x==0 and endpoint_limit))
    D=a-b;t=max(F(0),(D*D-x-z)/2)
    if x==0 and t==0:
        # Analytic one-sided limits, never evaluate an infinite inverse variance.
        pa=F(1);pb=F(0);mu=a;v=F(0);q=D*D/z
    else:
        wa=1/(x+t);wb=1/(z+t);pa=wa/(wa+wb);pb=wb/(wa+wb)
        mu=pa*a+pb*b;v=1/(wa+wb);q=wa*(a-mu)**2+wb*(b-mu)**2
    q0=D*D/(x+z);i=F(0) if q0==0 else max(F(0),100*(q0-1)/q0)
    vals={'rho_variance_A':x,'tau2':t,'weight_A':pa,'weight_B':pb,
        'mu':mu,'Vw':v,'q_HK':q,'V_HK':q*v,'V_mHK':max(F(1),q)*v,'Q':q0,'I2':i,'within_weight_A':z/(x+z)}
    def dd(v):
        with localcontext() as c:
            c.prec=55;return Decimal(v.numerator)/Decimal(v.denominator)
    with localcontext() as c:
        c.prec=55
        T=dd(sum(constants()['T2'])/2).sqrt();Z=dd(sum(constants()['Z2'])/2).sqrt()
        for name,var in [('CI_mHK',vals['V_mHK']),('CI_HK',vals['V_HK']),('PI_mHK',t+vals['V_mHK'])]:
            w=T*dd(var).sqrt();vals[name+'_low']=str(dd(mu)-w);vals[name+'_high']=str(dd(mu)+w)
        for name,y,var in [('LOO_A',a,x),('LOO_B',b,z)]:
            w=Z*dd(var).sqrt();vals[name+'_low']=str(dd(y)-w);vals[name+'_high']=str(dd(y)+w)
    return {k:(dec(v) if isinstance(v,F) else v) for k,v in vals.items()}
