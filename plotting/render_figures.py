from pathlib import Path
import json,runpy,sys,shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
import panels
from utility_panels import render_utility
h=vars(panels);panels.setup()
frozen=Path(__file__).resolve().parents[1]/'figures'
out=Path(sys.argv[1]) if len(sys.argv)>1 else Path('rendered')
out.mkdir(parents=True,exist_ok=True)
def build_s1(d):
    # Reuse the verified artist arrays, normalizations and colormaps.
    old=h['figS1'](d); data=[ax.images[0] for ax in old.axes if ax.images]
    vals=[np.array(img.get_array()) for img in data]; rows=d['homogeneous_cells']; null=d['null_cells']
    assert len(rows)==24 and [x['id'] for x in rows[::2]]==[x['id'] for x in null]
    fig=plt.figure(figsize=(7.15,5.8)); grid=fig.add_gridspec(2,4,left=.02,right=.99,top=.85,bottom=.135,width_ratios=[2.5,2.65,3.15,3.15],height_ratios=[1,.025],wspace=.16,hspace=.25)
    t=fig.add_subplot(grid[0,0]);t.set(xlim=(0,4),ylim=(23.5,-.5));t.axis('off')
    t.set_title('Field design',loc='left',size=8,pad=27)
    col=['$s_2$','$s_B$','ρ','θ']; keys=['s2','sb','true_rho','thetaA']
    for j,(label,key,span) in enumerate(zip(col,keys,[12,6,2,1])):
        t.text(j+.5,-1.2,label,ha='center',size=8)
        for i in range(0,24,span):
            assert all(r[key]==rows[i][key] for r in rows[i:i+span])
            t.text(j+.5,i+(span-1)/2,('%g'%rows[i][key]).replace('-','−'),va='center',ha='center',size=7)
    for i in (6,12,18):t.axhline(i-.5,xmin=0 if i==12 else .25,color='#bbc5cc',lw=.65)
    for j,(img,v,title,labels) in enumerate(zip(data,vals,['a   Null selection','b   Interval coverage','c   Mean interval width'],[['Point','Full','Oracle'],['HK','mHK','HK','mHK'],['HK','mHK','HK','mHK']])):
        ax=fig.add_subplot(grid[0,j+1]);ax.set_title(title,loc='left',size=8,pad=27)
        if j==0:
            arr=np.full((24,3),np.nan);arr[::2]=v
        else:arr=v
        im=ax.imshow(np.ma.masked_invalid(arr),interpolation='nearest',aspect='auto',norm=img.norm,cmap=img.cmap.with_extremes(bad='#f1f2f3'))
        np.testing.assert_array_equal(arr[::2] if j==0 else arr,v)
        ax.set_yticks([]);ax.set_xticks(range(len(labels)),labels,fontsize=6.7);ax.xaxis.tick_top();ax.tick_params(length=0,pad=3)
        if j:
            for u,name in [(.25,'ρ = 0'),(.75,'True ρ')]:ax.text(u,1.046,name,transform=ax.transAxes,ha='center',size=6.7)
        for idx in (6,12,18):ax.axhline(idx-.5,color='#bbc5cc',lw=.65)
        for (iy,ix),val in np.ndenumerate(arr):
            if np.isfinite(val):
                ax.text(ix,iy,format(val,'.2f' if j==2 else '.3f'),ha='center',va='center',size=6.3,color='white' if img.norm(val)>.60 else '#20262b')
        ca=fig.add_subplot(grid[1,j+1]);fig.colorbar(im,cax=ca,orientation='horizontal');ca.tick_params(labelsize=6,length=2)
    fig.suptitle('Taxon selection and interval precision across simulated field designs',x=.02,y=.988,ha='left',weight='bold',fontsize=9)
    fig.text(.02,.045,'Gray rows in a: null selection applies only at θ = 0.     Point: ρ = 0; Full: complete domain; Oracle: true ρ.',size=7)
    plt.close(old);return fig
mapping={'Supplementary_Figure_S2':'Supplementary_Figure_S2','Figure_1':'Figure_5','Figure_2':'Figure_2','Figure_3':'Figure_1','Figure_5':'Figure_4','Supplementary_Figure_S1':'Supplementary_Figure_S1'}
for name,old in mapping.items():
    path=out/'figures'/name;path.mkdir(parents=True,exist_ok=True)
    f=frozen/name/'data.json';d=json.loads(f.read_bytes())
    if name=='Supplementary_Figure_S1':fig=build_s1(d)
    else:
        fig=h['BUILD'][old](d);title=h['SPECS'][old][2]
        if name=='Figure_1':
            title='Low-nitrogen plots score higher on average despite a local reversal'
            fig.axes[-1].set_ylabel('Standardized rank score (LN − HN)')
            for ax in fig.axes[:2]:
                im=ax.images[0]
                for tx in ax.texts:
                    col,row=tx.get_position()
                    rgba=im.cmap(im.norm(im.get_array()[int(row),int(col)]))
                    lum=sum(w*((v/12.92) if v<=.04045 else ((v+.055)/1.055)**2.4) for w,v in zip([.2126,.7152,.0722],rgba))
                    tx.set_color('white' if lum<.179 else '#111111')
        if name=='Figure_2':fig.axes[1].get_legend().get_texts()[-1].set_text('Unflagged in A/B & full-domain stable')
        if name=='Figure_3':title='Cross-year error correlation governs taxon qualification across field-design simulations'
        if name=='Figure_5':
            fig.axes[1].set_title('b  Dense shift scenario (60/100 taxa, 12 pairs)',loc='left',fontsize=8,fontweight='bold',pad=8)
        fig.suptitle(title,x=.01,ha='left',fontsize=9,fontweight='bold',y=.97 if name=='Figure_3' else 1.055)
    for ext in ['png','svg','pdf']:fig.savefig(path/(name+'.'+ext),bbox_inches='tight',dpi=300)
    plt.close(fig)

for name in ('Figure_4','Supplementary_Figure_S3'):
    render_utility(name,frozen,out)
