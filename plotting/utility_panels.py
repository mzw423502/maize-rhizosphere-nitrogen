"""Render only independently audited saved decision summaries; never fit or audit."""
import hashlib,importlib.util,json,math,shutil
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties,findfont
from decimal import Decimal,ROUND_HALF_UP
import numpy as np
from PIL import Image
DOM=['zero','calibrated_point','calibrated_interval','full'];DL=['Zero','Cal.\npoint','Cal.\ninterval','Full']
WIDTHS=[0.,.1,.25,.5,1.,2.]
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def put(p,obj):p.write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False),encoding='utf8')
def read(p):return json.loads(p.read_text(encoding='utf8'))
def metric(row,name):return row['metrics'][name]['mean']
def short(row,profile=True):
    p=row['parameters'];pr={'sparse_shared':'S','sparse_heterogeneous':'H','dense_shared':'D'}[p['profile']]
    return f"{pr} | {p['pairs']} / {p['icc']:.1f} / {p['year_rho']:+.1f}"
def axe(fig,box,rect):
    x,y,w,h=box;l,b,aw,ah=rect;return fig.add_axes([x+w*l,y+h*b,w*aw,h*ah])
def matrix(ax,values,rows,cols,vmax,title,cmap='Blues'):
    im=ax.imshow(values,aspect='auto',vmin=0,vmax=vmax,cmap=cmap,interpolation='none')
    ax.set_yticks(range(len(rows)),rows,fontsize=6);ax.set_xticks(range(len(cols)),cols,fontsize=6)
    ax.tick_params(length=0);ax.set_title(title,loc='left',fontsize=7,fontweight='bold',pad=8)
    for i in range(len(rows)):
        for j in range(len(cols)):
            x=values[i,j]
            if np.isfinite(x):ax.text(j,i,str(Decimal(str(round(100*x,10))).quantize(Decimal('.1'),rounding=ROUND_HALF_UP)),ha='center',va='center',fontsize=6,color='white' if x/vmax>.62 else '#24303A')
    return im
def render_utility(current, data_root, output_root):
    bundle=read(data_root/current/'data.json')
    rows=bundle['rows'];lookup={(r['cell'],r['method'],r['q'],r['domain']):r for r in rows}
    cases=bundle['illustrative_cases']
    alt=[r for r in rows if r['method']=='LinDA_unit' and r['q']==.7 and r['domain']=='zero' and r['parameters']['profile']!='global_null']
    alt.sort(key=lambda r:({'sparse_shared':0,'sparse_heterogeneous':1,'dense_shared':2}[r['parameters']['profile']],r['cell']))
    hetero=[r for r in alt if r['parameters']['profile']=='sparse_heterogeneous']
    outlines={
      'Figure_3':{'claim':'The complete domain reduces discordant membership at a substantial true-signal retention cost; minimum caps identify why a decision changes.','width_mm':180,'ncol':12,'row_heights_mm':[66,62],'panels':[
        {'letter':'a','role':'hero','message':'Complete domains retain fewer true same-direction effects','chart_family':'matrix','row':0,'col':0,'colspan':6,'rowspan':2,'label_budget':0,'ask':'All 18 nonnull LinDA cells, four domains, q=0.7; true-same denominator, percentages.'},
        {'letter':'b','role':'primary','message':'Opposite-effect membership also falls','chart_family':'matrix','row':0,'col':6,'colspan':6,'label_budget':0,'ask':'All eight heterogeneous LinDA cells, four domains, opposite-truth denominator, percentages; do not call it FDR.'},
        {'letter':'c','role':'supporting','message':'A permissive cap can rely on heterogeneity','chart_family':'interval bars','row':1,'col':6,'colspan':6,'label_budget':4,'ask':'Two deterministic LinDA illustrative records; REML and inverse-variance minimal caps; highlight threshold 0.8 without interpreting agreement.'}]},
      'Supplementary_Figure_S3':{'claim':'Wider domains and stricter caps jointly reduce same-direction membership across the complete nonnull scenario set.','width_mm':180,'ncol':12,'row_heights_mm':[99,99],'panels':[{'letter':chr(97+i),'role':'primary','message':f'Maximum weight cap {q:.1f}','chart_family':'matrix','row':i//2,'col':(i%2)*6,'colspan':6,'label_budget':0,'ask':f'18 nonnull LinDA cells and six fixed calibrated half-widths, q={q}; all true-same retention fractions.'} for i,q in enumerate([.6,.7,.8,.9])]}}
    def draw(name,letter,fig,box):
        if name=='Figure_3' and letter in ['a','b']:
            is_a=letter=='a';selected=alt if is_a else hetero
            values=np.array([[metric(lookup[(r['cell'],'LinDA_unit',.7,d)],'true_same_retention' if is_a else 'true_opposite_membership') for d in DOM] for r in selected])
            ax=axe(fig,box,[.34,.17 if is_a else .24,.60,.73 if is_a else .59])
            im=matrix(ax,values,[short(r) for r in selected],DL,1 if is_a else .16,'True same-direction retained (%)' if is_a else 'Opposite-effect members (%)','Blues' if is_a else 'Oranges')
            if is_a:
                for y in [7.5,15.5]:ax.axhline(y,color='white',lw=2)
            cbax=axe(fig,box,[.34,.060 if is_a else .065,.60,.016 if is_a else .027]);cb=fig.colorbar(im,cax=cbax,orientation='horizontal');cb.ax.tick_params(labelsize=6,length=2)
            cb.set_ticks([0,.5,1] if is_a else [0,.08,.16],labels=['0','50','100'] if is_a else ['0','8','16'])
            return
        if name=='Figure_3':
            ax=axe(fig,box,[.29,.30,.65,.49]);keys=['LinDA_unit_heterogeneity_assisted','LinDA_unit_inverse_variance_sufficient']
            for i,k in enumerate(keys):
                r=cases[k];a=r['qstar_REML'];b=r['qstar_IV'];y=1-i
                ax.barh(y,.5,left=.5,height=.34,color='#EDF0F2')
                ax.barh(y,b-a,left=a,height=.34,color='#CFAF8E')
                ax.plot(a,y,'|',ms=16,mew=2.8,color='#386A8C');ax.plot(b,y,'|',ms=16,mew=2.8,color='#A65E34')
                ax.text(a,y+.28,f'{a:.3f}',color='#386A8C',ha='center',fontsize=6)
                ax.text(b,y-.38,f'{b:.3f}',color='#A65E34',ha='center',fontsize=6)
            ax.axvline(.8,color='#77828A',ls=':',lw=.9)
            ax.set(xlim=(.5,1),ylim=(-.6,1.55),yticks=[0,1],yticklabels=['IV sufficient','Heterogeneity\nassisted'],xlabel='Minimum complete-domain cap')
            ax.tick_params(labelsize=6);ax.set_title('Similar caps, different mechanisms',loc='left',fontsize=7,fontweight='bold',pad=10)
            ax.plot([],[],'|',ms=10,mew=2.8,color='#386A8C',label='REML');ax.plot([],[],'|',ms=10,mew=2.8,color='#A65E34',label='Inverse variance')
            ax.legend(loc='upper left',bbox_to_anchor=(-.03,-.42),ncol=2,fontsize=6,handlelength=.8,columnspacing=1)
            return
        qi=ord(letter)-97;q=[.6,.7,.8,.9][qi]
        values=np.array([[metric(lookup[(r['cell'],'LinDA_unit',q,'width_'+str(w))],'true_same_retention') for w in WIDTHS] for r in alt])
        ax=axe(fig,box,[.33,.19,.63,.70]);im=matrix(ax,values,[short(r) for r in alt],[str(w).rstrip('0').rstrip('.') for w in WIDTHS],1,f'Cap q = {q:.1f}')
        ax.set_xlabel('Calibrated-domain half-width',fontsize=6)
        for y in [7.5,15.5]:ax.axhline(y,color='white',lw=2)
        cbax=axe(fig,box,[.33,.055,.63,.021]);cb=fig.colorbar(im,cax=cbax,orientation='horizontal');cb.set_ticks([0,.5,1],labels=['0%','50%','100%']);cb.ax.tick_params(labelsize=6,length=2)
    name='Figure_3' if current=='Figure_4' else current
    outline=outlines[name]; heights=outline['row_heights_mm']; gap=0.016
    height=sum(heights);fig=plt.figure(figsize=(180/25.4,height/25.4))
    for p in outline['panels']:
        row=p['row']; x=p['col']/12; top=sum(heights[:row])/height
        w=p['colspan']/12; hgt=sum(heights[row:row+p.get('rowspan',1)])/height
        box=(x+gap,1-top-hgt+gap,w-2*gap,hgt-2*gap)
        draw(name,p['letter'],fig,box);fig.text(x+gap,1-top-gap,p['letter'],fontweight='bold',va='top')
    target=output_root/'figures'/current
    target.mkdir(parents=True,exist_ok=True)
    for ext in ('pdf','svg','png'):
        fig.savefig(target/(current+'.'+ext),dpi=300,bbox_inches='tight')
    plt.close(fig)
