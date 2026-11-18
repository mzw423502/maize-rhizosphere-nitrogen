"""Reusable deterministic two-study balance certificate. No inferential p values."""
import math
def max_weight(x,z,d):
    if x<0 or z<=0:raise ValueError('Require nonnegative x and positive z')
    return .5+abs(x-z)/(2*max(x+z,d*d))
def cap_required(alpha,beta,z,d):
    if beta<0 or alpha<beta or z<=0:raise ValueError('Invalid variance domain')
    return max(max_weight(alpha-beta,z,d),max_weight(alpha+beta,z,d))
def window(z,d,q):
    if not (.5<=q<1 and z>0):raise ValueError('Invalid cap or variance')
    c=2*q-1
    return max(0,min(z*(1-q)/q,z-c*d*d)),max(z*q/(1-q),z+c*d*d)
def classify(a,b,s1,s2,sb,q=.7):
    alpha=(s1*s1+s2*s2)/4;beta=s1*s2/2;z=sb*sb;d=a-b
    if min(s1,s2,sb)<=0 or not all(map(math.isfinite,[a,b,s1,s2,sb])):return {'status':'UNCLASSIFIABLE'}
    # alpha-beta written as a square avoids a negative cancellation artifact.
    xmin=(s1-s2)**2/4;xmax=(s1+s2)**2/4
    lo,hi=window(z,d,q);overlo=max(xmin,lo);overhi=min(xmax,hi)
    regular_nonempty=(overlo<=overhi and overhi>0)
    required=max(max_weight(xmin,z,d),max_weight(xmax,z,d))
    reqiv=max(max_weight(xmin,z,0),max_weight(xmax,z,0))
    allw=required<=q;anyw=regular_nonempty
    same=(a*b)>0
    return dict(status='REGULAR' if xmin>0 else 'ZERO_ENDPOINT_LIMIT',a=a,b=b,s1=s1,s2=s2,sb=sb,alpha=alpha,beta=beta,z=z,xmin=xmin,xmax=xmax,q=q,required_cap=required,required_cap_iv=reqiv,sign_agree=same,weight_class=('ALL' if allw else 'SOME' if anyw else 'NONE'),decision_class=('ALL' if same and allw else 'SOME' if same and anyw else 'NONE'),rho0_pass=bool(same and max_weight(alpha,z,d)<=q),all_pass=bool(same and allw),heterogeneity_assisted=bool(allw and reqiv>q),rho_low=(max(-1,(lo-alpha)/beta) if beta>0 and anyw else None),rho_high=(min(1,(hi-alpha)/beta) if beta>0 and anyw else None),zero_endpoint_nonregular=bool(xmin==0))
def box_cap(x_range,z_range,d_range):
    x0,x1=x_range;z0,z1=z_range;d0,d1=d_range
    if not (0<x0<=x1 and 0<z0<=z1 and d0<=d1):raise ValueError('Invalid positive-variance box')
    d=0 if d0<=0<=d1 else min(abs(d0),abs(d1))
    return max(max_weight(x,z,d) for x in [x0,x1] for z in [z0,z1])
def pooled(a,b,x,z):
    if x<=0 or z<=0:return None
    tau=max(0,((a-b)**2-x-z)/2);wa=1/(x+tau);wb=1/(z+tau);w=wa+wb
    mean=(wa*a+wb*b)/w;Q=wa*(a-mean)**2+wb*(b-mean)**2
    return dict(mean=mean,tau2=tau,HK_se=math.sqrt(Q/w),mHK_se=math.sqrt(max(1,Q)/w),max_weight=max(wa,wb)/w)
