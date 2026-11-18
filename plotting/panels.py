"""Report-layer redesign of six frozen historical figures.

This script never fits a model, resamples, or recalculates scientific outputs.  It
copies the supplied frozen JSON into a new isolated asset tree, hashes it, and
renders only deterministic display transforms of the saved estimates.  MCSE is
retained verbatim in each copied JSON and is described in the captions; it is
not silently converted into graphical error bars.
"""
from __future__ import annotations

import argparse, hashlib, importlib.util, json, shutil, textwrap
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
from PIL import Image

MM = 1 / 25.4
BLUE = "#386A8C"; BLUE2 = "#84AEC4"; ORANGE = "#B36A3C"; GOLD = "#D8B56A"; INK = "#20262B"; GREY = "#E7EAED"

SPECS = {
    "Figure_1": ("Figure_1", "Figure_1", "Covariance-domain qualification across 26 end-to-end scenarios", "A common scenario matrix makes visible how pointwise and full-domain membership, heterogeneity, and weight caps jointly define a qualification decision."),
    "Figure_2": ("Figure_2", "Figure_2", "Exact variance window and its classification consequences", "The analytic window separates full-, partial-, and non-qualification while retaining the diagnostic subset."),
    "Figure_4": ("Figure_3", "Figure_4", "Unit structure, reference assumptions and interval precision", "Global-null error, annual coverage, and interval behavior answer distinct estimands and must not be conflated."),
    "Figure_5": ("Figure_6", "Figure_5", "Transported directional scores retain pairwise reversals", "A shared zero-centered score scale exposes both positive average transfer and directionally discordant pairs."),
    "Supplementary_Figure_S1": ("Figure_4", "Supplementary_Figure_S1", "Gaussian qualification and calibration are distinct properties", "Selection under a full covariance domain and interval coverage address different questions."),
    "Supplementary_Figure_S2": ("Figure_5", "Supplementary_Figure_S2", "Allocation-unit structure governs calibration precision", "Global-null error and latent-unit coverage change systematically with paired units, subsampling, and intraclass correlation."),
}

def sha(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(1<<20), b""): h.update(b)
    return h.hexdigest()

def setup():
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8,"axes.spines.right":False,"axes.spines.top":False,"axes.titlesize":8,"savefig.dpi":300})
    plt.rcParams.update({"figure.facecolor":"white","axes.facecolor":"white","savefig.facecolor":"white","axes.labelcolor":INK,"text.color":INK,"xtick.color":INK,"ytick.color":INK})

def metric(x, name, key="mean"):
    v=x["metrics"].get(name)
    return np.nan if v is None else v.get(key, np.nan)

def lab(p):
    return f"g={p['pairs']}  c={p['icc']:.1f}  r={p['year_rho']:+.1f}  |  {p['profile'].replace('_',' ')}"

def heat(ax, a, rows, cols, title, cmap="Blues", norm=None, fmt=".2f", group_breaks=()):
    im=ax.imshow(a, aspect="auto", interpolation="nearest", cmap=cmap, norm=norm)
    ax.set_title(title, loc="left", fontweight="bold", pad=8)
    ax.set_xticks(np.arange(len(cols)), cols, rotation=0)
    ax.set_yticks(np.arange(len(rows)), rows, fontsize=6)
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            val=a[i,j]
            if np.isfinite(val): ax.text(j,i,format(val,fmt),ha="center",va="center",fontsize=6,color="white" if (norm(val) if norm else val)>0.60 else INK)
    for cut in group_breaks: ax.axhline(cut-.5,color="white",lw=1.4)
    ax.tick_params(length=0); return im

def fig1(d):
    rows=[]; values=[]
    # Exactly 26 scenario cells at cap=0.7 and unit-aware LinDA; all saved records remain copied.
    rec=[x for x in d["decision"] if x["method"]=="LinDA_unit" and abs(x["cap"]-.7)<1e-10]
    profile_order={"global_null":0,"sparse_shared":1,"sparse_heterogeneous":2,"dense_shared":3}
    rec.sort(key=lambda x:(profile_order.get(x["parameters"]["profile"],99),x["parameters"]["pairs"],x["parameters"]["icc"],x["parameters"]["year_rho"]))
    for x in rec:
        m=x["metrics"]; rows.append(lab(x["parameters"])); values.append([metric(x,"rho0_membership"),metric(x,"all_domain_membership"),metric(x,"rho0_qualified_not_full_fraction"),metric(x,"heterogeneity_assisted_ALL"),metric(x,"mean_minimal_cap_REML"),metric(x,"mean_minimal_cap_IV")])
    a=np.array(values); fig,ax=plt.subplots(figsize=(180*MM,112*MM))
    fig.subplots_adjust(left=.245,right=.90,top=.89,bottom=.09)
    breaks=[]; cursor=0
    for profile in ["global_null","sparse_shared","sparse_heterogeneous","dense_shared"]:
        n=sum(x["parameters"]["profile"]==profile for x in rec); cursor+=n
        if n and cursor<len(rec): breaks.append(cursor)
    im=heat(ax,a,rows,["at ρ=0","all ρ","lost after\nfull domain","heterogeneity-\nassisted","REML q*","IV q*"],"",cmap="YlGnBu",norm=plt.Normalize(0,1),group_breaks=breaks)
    ax.axvline(1.5,color=INK,lw=.8); ax.axvline(3.5,color=INK,lw=.8)
    # group headers sit inside the axes margin, not over row labels or suptitle.
    for xpos,label in [(1/12,"Membership"),(3/6,"Domain consequence"),(5/6,"Required cap")]:
        ax.text(xpos,1.014,label,transform=ax.transAxes,ha="center",va="bottom",fontsize=7,fontweight="bold")
    c=fig.colorbar(im,ax=ax,pad=.012,shrink=.42); c.set_label("proportion or weight cap")
    ax.set_ylabel("end-to-end scenario (g, c, year correlation r, profile)")
    return fig

def fig2(d):
    x=np.asarray(d["d_squared_over_z"]); lo=np.asarray(d["lower_x_over_z"]); hi=np.asarray(d["upper_x_over_z"])
    fig,axs=plt.subplots(1,2,figsize=(180*MM,68*MM),gridspec_kw={"width_ratios":[1.22,1]},layout="constrained")
    ax=axs[0]; ax.plot(x,lo,color=BLUE,lw=1.3,label="lower boundary"); ax.plot(x,hi,color=ORANGE,lw=1.3,label="upper boundary"); ax.fill_between(x,lo,hi,color=BLUE2,alpha=.24,label="admissible window"); ax.axhline(1,color=INK,lw=.8,ls="--",label="equal variance (x/z=1)")
    ax.set(xlabel="d² / z",ylabel="x / z",title="a  Exact analytic window")
    # Keep all analytic curves unobscured: the key sits outside the left axes.
    ax.legend(frameon=False,ncol=2,fontsize=6,loc="upper center",bbox_to_anchor=(.50,-.34))
    # Horizontal q rows preserve the literal four caps, every 214-label category
    # count, and the diagnostic intersection without placing a legend over data.
    cats=["ALL_RHO","SOME_RHO","NONE","NOT_ASSESSABLE"]; catlabs=["Full","Partial","None","Not assessable"]; colors=[BLUE,GOLD,ORANGE,"#B9BEC4"]
    qlabs=["q=0.60","q=0.70","q=0.80","q=0.90"]
    for yi,r in enumerate(d["classes"]):
        c=r["classes"]; left=0
        for cat,col in zip(cats,colors):
            width=c[cat]; axs[1].barh(yi,width,left=left,color=col,height=.62)
            if width>=22 or cat=="ALL_RHO": axs[1].text(left+width/2,yi,str(width),ha="center",va="center",fontsize=5.6 if width<22 else 6,color="white" if cat in ("ALL_RHO","NONE") else INK)
            left+=width
        # Exact subset is a narrow outlined segment within Full, labelled explicitly.
        diag=r["classes_with_original_AB_no_captured_flag"].get("ALL_RHO",0)
        axs[1].add_patch(plt.Rectangle((0,yi-.31),diag,.62,fill=False,ec=INK,lw=1.05,zorder=4))
        axs[1].text(216,yi,f"subset={diag}",ha="left",va="center",fontsize=5.8,color=INK,clip_on=False)
    assert [r["classes_with_original_AB_no_captured_flag"].get("ALL_RHO",0) for r in d["classes"]]==[1,2,3,3]
    axs[1].set(title="b  Classification of all 214 labels",xlabel="labels",xlim=(0,244),yticks=range(4),yticklabels=qlabs)
    axs[1].invert_yaxis(); axs[1].grid(axis="x",lw=.35,color="#D8DDE1")
    from matplotlib.patches import Patch
    handles=[Patch(facecolor=c,label=l) for c,l in zip(colors,catlabs)]+[Patch(facecolor="none",edgecolor=INK,label="outlined: original A/B no-flag ∩ Full")]
    axs[1].legend(handles=handles,frameon=False,fontsize=5.5,ncol=2,loc="upper left",bbox_to_anchor=(0,-.24))
    return fig

def fig4(d):
    est=d["estimator"]; inter=d["interval"]; years=["A1","A2","B"]
    methods=["LinDA_unit","CLR_unit","CLR_sample_independence_control"]; method_labels=["LinDA unit","CLR unit","CLR sample"]
    # Exact global-null projection: eight cells × three annual estimands = 24 rows.
    null_cells=sorted({x["cell"] for x in est if x["parameters"]["profile"]=="global_null"})
    rows=[(c,y) for c in null_cells for y in years]
    assert len(rows)==24 and all(any(x["cell"]==c and x["year"]==y and x["method"]==m for x in est) for c,y in rows for m in methods)
    a=np.array([[next(metric(x,"BH_FDP") for x in est if x["cell"]==c and x["year"]==y and x["method"]==m) for m in methods] for c,y in rows])
    # Explicit dense cell 25 projection; never substitute a first-N proxy.
    b=np.array([[next(metric(x,"coverage") for x in est if x["cell"]==25 and x["year"]==y and x["method"]==m) for m in ["LinDA_unit","CLR_unit"]] for y in years])
    assert np.allclose(b[:,0],[.65355,.66805,.6544]) and np.allclose(b[:,1],[.94625,.9438,.9473])
    rhos=[-1,0,1]
    cov=np.array([[next(metric(x,"coverage") for x in inter if x["cell"]==25 and x["method"]=="LinDA_unit" and x["interval"]=="mHK" and x["rho"]==rho) for rho in rhos]])
    wid=np.array([[next(metric(x,"mean_width") for x in inter if x["cell"]==25 and x["method"]=="LinDA_unit" and x["interval"]=="mHK" and x["rho"]==rho) for rho in rhos]])
    assert np.allclose(cov[0],[.64105,.99995,1.0]) and np.allclose(wid[0],[4.296695277239986,6.767796355600095,7.47383038646405])
    fig=plt.figure(figsize=(180*MM,112*MM),layout="constrained"); gs=fig.add_gridspec(2,2,width_ratios=[1,1.15]); ax1=fig.add_subplot(gs[:,0]); ax2=fig.add_subplot(gs[0,1])
    im1=heat(ax1,a,[f"null {c+1} · {y}" for c,y in rows],method_labels,"a  Global-null BH FDP (24 scenario-year rows)",cmap="Oranges",norm=plt.Normalize(0,1),fmt=".3f")
    im2=heat(ax2,b,years,["LinDA unit","CLR unit"],"b  Coverage: dense cell 25",cmap="Blues",norm=plt.Normalize(.60,1),fmt=".3f")
    # Coverage and width are separate, aligned horizontal bars in the same bottom-right panel.
    low=gs[1,1].subgridspec(1,2,wspace=.24); ca=fig.add_subplot(low[0,0]); wa=fig.add_subplot(low[0,1])
    ypos=np.arange(3); rho_labels=["ρ=-1","ρ=0","ρ=1"]
    ca.barh(ypos,cov[0],color=BLUE2,height=.58); ca.axvline(.95,color=INK,lw=.8,ls="--",label="nominal 0.95")
    ca.set(yticks=ypos,yticklabels=rho_labels,xlim=(.55,1.03),xlabel="coverage",title="c  LinDA mHK coverage"); ca.invert_yaxis(); ca.legend(frameon=False,fontsize=5.6,loc="lower left")
    for y,v in zip(ypos,cov[0]): ca.text(min(v+.012,1.018),y,f"{v:.5f}",va="center",fontsize=5.8)
    wa.barh(ypos,wid[0],color=GOLD,height=.58); wa.set(yticks=ypos,yticklabels=[],xlim=(0,max(wid[0])*1.20),xlabel="mean interval width",title="LinDA mHK width"); wa.invert_yaxis()
    for y,v in zip(ypos,wid[0]): wa.text(v+max(wid[0])*.025,y,f"{v:.3f}",va="center",fontsize=5.8)
    for im,ax in [(im1,ax1),(im2,ax2)]: fig.colorbar(im,ax=ax,pad=.01,shrink=.58)
    return fig

def fig5(d):
    a=np.asarray(d["mean_matrix"],float); p=np.asarray(d["primary_matrix"],float); labels=d["target_labels"]
    settings=[f"p={pr:.2f}\nk={int(k)}" for pr,k in d["settings"]]
    vmax=max(abs(np.nanmin(a)),abs(np.nanmax(a)),abs(np.nanmin(p)),abs(np.nanmax(p))); norm=TwoSlopeNorm(vmin=-vmax,vcenter=0,vmax=vmax)
    fig,axs=plt.subplots(2,1,figsize=(180*MM,105*MM),gridspec_kw={"height_ratios":[1.25,1]},layout="constrained")
    shorthand={
        "B|2022|fieldRep|1":"A→B pair 1", "B|2022|fieldRep|2":"A→B pair 2",
        "A|Y2018|NS|N":"B→A2018 North", "A|Y2018|NS|S":"B→A2018 South",
        "A|Y2019|NS|N":"B→A2019 North", "A|Y2019|NS|S":"B→A2019 South"}
    pair_labels=[shorthand[x[2]] for x in d["primary_pair_keys"]]
    for ax,mat,rows,title in [(axs[0],a,labels,"a  Mean transported directional score"),(axs[1],p,pair_labels,"b  Primary-pair score")]:
        im=heat(ax,mat,rows,settings,title,cmap="RdBu_r",norm=norm,fmt="+.2f")
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                if mat[i,j]<0: ax.add_patch(plt.Rectangle((j-.5,i-.5),1,1,fill=False,ec=INK,lw=.8))
        for x in [2.5,5.5]: ax.axvline(x,color=INK,lw=.7)
    c=fig.colorbar(im,ax=axs,pad=.01,shrink=.55); c.set_label("saved directional score; outlined cells are negative")
    return fig

def figS1(d):
    nul=d["null_cells"]; hom=d["homogeneous_cells"]
    def keys(rows): return [f"id={x['id']}; s2={x['s2']:.1f}; sB={x['sb']:.1f}; ρ={x['true_rho']:+.1f}; θ={x['thetaA']:.1f}" for x in rows]
    a=np.array([[x['metrics']['point0_selection']['estimate'],x['metrics']['full_regular_selection']['estimate'],x['metrics']['true_rho_selection']['estimate']] for x in nul])
    b=np.array([[x['metrics']['rho0_HK']['coverage'],x['metrics']['rho0_mHK']['coverage'],x['metrics']['true_rho_HK']['coverage'],x['metrics']['true_rho_mHK']['coverage']] for x in hom])
    w=np.array([[x['metrics'][z]['mean_length'] for z in ['rho0_HK','rho0_mHK','true_rho_HK','true_rho_mHK']] for x in hom])
    fig=plt.figure(figsize=(180*MM,150*MM),layout="constrained"); gs=fig.add_gridspec(2,2,width_ratios=[.82,1.4]); axa=fig.add_subplot(gs[:,0]); axb=fig.add_subplot(gs[0,1]); axc=fig.add_subplot(gs[1,1])
    colsel=["ρ=0\npoint","full\ndomain","true ρ\noracle"]
    colint=["ρ=0\nHK","ρ=0\nmHK","true ρ\nHK","true ρ\nmHK"]
    im1=heat(axa,a,keys(nul),colsel,"a  Selection frequency under the global null",cmap="YlOrBr",norm=plt.Normalize(0,.55),fmt=".3f")
    im2=heat(axb,b,keys(hom),colint,"b  Interval coverage in homogeneous settings",cmap="Blues",norm=plt.Normalize(.75,1),fmt=".3f")
    im3=heat(axc,w,["" for _ in hom],colint,"c  Mean interval width in the same settings",cmap="YlOrBr",norm=plt.Normalize(np.nanmin(w),np.nanmax(w)),fmt=".2f")
    fig.colorbar(im1,ax=axa,pad=.01,shrink=.62); fig.colorbar(im2,ax=axb,pad=.01,shrink=.62); fig.colorbar(im3,ax=axc,pad=.01,shrink=.62)
    return fig

def figS2(d):
    gn=d["global_null_cells"]; cov=d["all_cells_latent_coverage"]
    methods=["sample_Welch","unit_paired_t","exact_pair_randomization"]
    ga=np.array([[next(y['FDR']['mean'] for y in x['methods'] if y['method']==m) for m in methods] for x in gn])
    gl=[f"g={x['pairs']}; n={x['subsamples']}; ICC={x['latent_ICC']:.1f}" for x in gn]
    groups=["null","balanced","unbalanced"]
    ids=sorted({(x['pairs'],x['subsamples'],x['latent_ICC']) for x in cov})
    assert len(ids)==18 and all(sum(x['pattern']==g and (x['pairs'],x['subsamples'],x['latent_ICC'])==i for x in cov)==1 for g in groups for i in ids)
    cover=np.array([[next(x['latent_unit_t_coverage']['mean'] for x in cov if x['pattern']==g and (x['pairs'],x['subsamples'],x['latent_ICC'])==i) for g in groups] for i in ids])
    labs=[f"g={p}; n={s}; ICC={icc:.1f}" for p,s,icc in ids]
    fig=plt.figure(figsize=(180*MM,126*MM),layout="constrained"); gs=fig.add_gridspec(1,2,width_ratios=[1.15,1.35]); a=fig.add_subplot(gs[0]); b=fig.add_subplot(gs[1]); im=heat(a,ga,gl,["Welch","paired t","Pair sign"],"a  Global-null FDR",cmap="Oranges",norm=plt.Normalize(0,1),fmt=".3f")
    im2=heat(b,cover,labs,groups,"b  Latent-unit coverage (all 54 scenario-profile cells)",cmap="Blues",norm=plt.Normalize(.80,1),fmt=".3f")
    fig.colorbar(im,ax=a,pad=.01,shrink=.45); fig.colorbar(im2,ax=b,pad=.01,shrink=.45)
    return fig

BUILD={"Figure_1":fig1,"Figure_2":fig2,"Figure_4":fig4,"Figure_5":fig5,"Supplementary_Figure_S1":figS1,"Supplementary_Figure_S2":figS2}