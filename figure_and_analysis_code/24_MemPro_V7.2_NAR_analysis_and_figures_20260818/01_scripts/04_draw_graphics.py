from pathlib import Path
import json
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np


ROOT=Path(r"D:\finale\24_MemPro_V7.2_NAR_analysis_and_figures_20260818")
OUT=ROOT/"05_graphics"
PALETTE=["#7b95c6","#49c2d9","#a1d8e8","#67a583","#a2c986","#d0e2c0","#fded95","#ffc1a6","#f59c7c","#f47254","#c85e62"]
BLUE,CYAN,LTBLUE,GREEN,LTGREEN,PALEGREEN,YELLOW,PEACH,SALMON,ORANGE,RED=PALETTE
DARK="#26323a";MID="#69757d";LIGHT="#edf1f3"
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],"font.size":7,"pdf.fonttype":42,"svg.fonttype":"none","text.color":DARK})


def save(fig,name):
    fig.savefig(OUT/f"{name}.png",dpi=450,bbox_inches="tight",facecolor="white")
    fig.savefig(OUT/f"{name}.svg",bbox_inches="tight",facecolor="white")
    fig.savefig(OUT/f"{name}.pdf",bbox_inches="tight",facecolor="white")
    plt.close(fig)


def rounded(ax,x,y,w,h,text,color,fs=7,weight="normal"):
    ax.add_patch(patches.FancyBboxPatch((x,y),w,h,boxstyle="round,pad=.012,rounding_size=.018",fc=color,ec="white",lw=.8))
    ax.text(x+w/2,y+h/2,text,ha="center",va="center",fontsize=fs,fontweight=weight)


def arrow(ax,x0,y0,x1,y1):
    ax.annotate("",xy=(x1,y1),xytext=(x0,y0),arrowprops=dict(arrowstyle="-|>",lw=1.1,color=MID,shrinkA=2,shrinkB=2))


def bilayer(ax,x,y,w,h=.08):
    xs=np.linspace(x,x+w,24)
    for yy in [y,y+h]:
        ax.scatter(xs,[yy]*len(xs),s=18,c=LTBLUE,edgecolors=GREEN,linewidths=.25,zorder=2)
    for xx in xs[::2]:
        ax.plot([xx,xx-.006],[y+.006,y+h/2],color=CYAN,lw=.6);ax.plot([xx,xx+.006],[y+h-.006,y+h/2],color=CYAN,lw=.6)


def protein_helix(ax,x,y,scale=1,color=BLUE):
    t=np.linspace(0,4*np.pi,100);xx=x+.018*scale*np.sin(t);yy=y+.10*scale*(t/(4*np.pi)-.5);ax.plot(xx,yy,color=color,lw=3,solid_capstyle="round")


def molecule_icon(ax,x,y,scale=.04,color=ORANGE):
    theta=np.linspace(0,2*np.pi,7)+np.pi/6;ax.plot(x+scale*np.cos(theta),y+scale*np.sin(theta),color=color,lw=1.5)
    ax.plot([x+scale*.87,x+scale*1.6],[y+scale*.5,y+scale*.85],color=color,lw=1.5);ax.scatter([x+scale*1.7],[y+scale*.9],s=18,color=RED)


def graphical_abstract():
    fig,ax=plt.subplots(figsize=(12,6.4));ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis("off")
    ax.text(.02,.96,"MemPro",fontsize=24,fontweight="bold",color=BLUE,va="top");ax.text(.02,.88,"Evidence-aware membrane-target knowledge for drug-discovery hypothesis generation",fontsize=10,color=MID,va="top")
    # Sources
    ax.text(.11,.79,"PUBLIC DATA LAYERS",ha="center",fontweight="bold",fontsize=8)
    groups=[("Protein", "UniProt · HPA\nHTP · Membranome",LTBLUE),("Interaction","ChEMBL · BindingDB\nPubChem · BRENDA",PALEGREEN),("Structure","PDBe · OPM · PDBTM\nPDBbind",YELLOW),("Disease","Open Targets · UniProt\nMONDO · DO · Uberon",PEACH)]
    for i,(lab,txt,col) in enumerate(groups):
        y=.66-i*.14;rounded(ax,.02,y,.18,.105,f"{lab}\n{txt}",col,fs=6.5,weight="bold" if i==0 else "normal")
    arrow(ax,.21,.47,.28,.47)
    # Standardization gate
    rounded(ax,.28,.25,.18,.47,"CANONICALISATION\n\nProtein · isoform · form\nCompound parent · form\nMONDO disease identity\n\nDEDUPLICATION\nLINEAGE · E / BE TIERS\nMISSINGNESS SEMANTICS",LIGHT,fs=7,weight="bold")
    arrow(ax,.47,.47,.53,.47)
    # Core knowledge graph
    ax.add_patch(patches.Circle((.65,.49),.18,fc="#f8fafb",ec=BLUE,lw=1.4))
    bilayer(ax,.55,.57,.20,.075);protein_helix(ax,.61,.61,1.15,BLUE);molecule_icon(ax,.70,.58,.035,ORANGE)
    ax.add_patch(patches.Circle((.64,.41),.045,fc=PEACH,ec=RED,lw=.8));ax.text(.64,.41,"D",ha="center",va="center",fontweight="bold",color=RED)
    ax.plot([.62,.57],[.45,.55],color=MID,lw=.8);ax.plot([.67,.58],[.45,.55],color=MID,lw=.8);ax.plot([.67,.45],[.66,.44],color=MID,lw=.8)
    ax.text(.65,.75,"MEMPRO MULTILAYER CORE",ha="center",fontweight="bold",fontsize=8)
    satellites=[(.53,.36,"Five-axis\nclassification",LTBLUE),(.74,.35,"Disease · anatomy\nexpression",PEACH),(.53,.78,"Evidence\nprovenance",PALEGREEN),(.76,.75,"Sites · structures\ncomplexes",YELLOW)]
    for x,y,txt,col in satellites:rounded(ax,x,y,.13,.075,txt,col,fs=5.8)
    arrow(ax,.84,.49,.89,.49)
    # Outcomes
    ax.text(.94,.79,"HYPOTHESIS OUTPUTS",ha="center",fontweight="bold",fontsize=8)
    outputs=[("Protein–compound\npreference",LTBLUE),("Evidence-aware\ninterpretation",PALEGREEN),("Disease-relevant,\nchemically underexplored\ntargets",PEACH)]
    for i,(txt,col) in enumerate(outputs):rounded(ax,.88,.63-i*.18,.115,.12,txt,col,fs=6.2,weight="bold")
    ax.text(.50,.05,"Internal design draft — statistical associations support prioritisation, not causal or experimental validation.",ha="center",fontsize=6,color=RED)
    save(fig,"Graphical_Abstract_MemPro_V72_internal_draft")


def section_protein():
    fig,ax=plt.subplots(figsize=(8,4.5));ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis("off");ax.text(.03,.94,"Protein module | membrane attachment and five-axis identity",fontweight="bold",fontsize=12)
    bilayer(ax,.04,.47,.48,.12)
    protein_helix(ax,.14,.55,1.5,BLUE);ax.text(.14,.34,"A · transmembrane",ha="center",fontweight="bold")
    protein_helix(ax,.30,.63,.45,GREEN);ax.plot([.30,.30],[.55,.59],color=ORANGE,lw=1.2);ax.scatter([.30],[.55],s=35,color=ORANGE);ax.text(.30,.34,"B · lipid/monolayer",ha="center",fontweight="bold")
    ax.add_patch(patches.Ellipse((.45,.67),.11,.08,fc=PEACH,ec=RED,lw=.8));ax.plot([.45,.45],[.59,.63],color=MID,lw=1);ax.text(.45,.34,"C · peripheral",ha="center",fontweight="bold")
    labels=["Structural family","Molecular function","Biological process","Membrane role","Specialist class"]
    center=(.77,.55);ax.add_patch(patches.Circle(center,.10,fc=LIGHT,ec=BLUE,lw=1.2));ax.text(*center,"One protein\nmultiple valid views",ha="center",va="center",fontweight="bold")
    for i,label in enumerate(labels):
        a=2*np.pi*i/5+np.pi/2;x=center[0]+.19*np.cos(a);y=center[1]+.29*np.sin(a);rounded(ax,x-.09,y-.045,.18,.09,label,PALETTE[i*2],fs=6);ax.plot([center[0]+.10*np.cos(a),x-.09*np.cos(a)],[center[1]+.10*np.sin(a),y-.045*np.sin(a)],color=MID,lw=.6)
    save(fig,"Section_Graphic_Protein")


def section_compound():
    fig,ax=plt.subplots(figsize=(8,4.5));ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis("off");ax.text(.03,.94,"Compound module | identity, scaffold and target breadth",fontweight="bold",fontsize=12)
    molecule_icon(ax,.15,.58,.09,BLUE);ax.text(.15,.38,"Canonical parent",ha="center",fontweight="bold")
    arrow(ax,.25,.58,.34,.58);molecule_icon(ax,.42,.68,.055,GREEN);molecule_icon(ax,.42,.48,.055,ORANGE);ax.text(.42,.32,"Salt · stereoisomer · charge form",ha="center")
    arrow(ax,.51,.58,.60,.58);rounded(ax,.61,.48,.14,.20,"Murcko scaffold\nMorgan fingerprint",LTBLUE,fs=7,weight="bold")
    arrow(ax,.76,.58,.82,.58);ax.scatter([.90],[.58],s=300,color=PEACH,edgecolors=RED);ax.scatter([.85,.92,.95],[.72,.78,.40],s=[70,90,60],color=[BLUE,GREEN,ORANGE]);
    for x,y in [(.85,.72),(.92,.78),(.95,.40)]:ax.plot([.90,x],[.58,y],color=MID,lw=.8)
    ax.text(.90,.28,"Selective ↔ polypharmacological",ha="center",fontweight="bold")
    save(fig,"Section_Graphic_Compound")


def section_evidence():
    fig,ax=plt.subplots(figsize=(8,4.5));ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis("off");ax.text(.03,.94,"Interaction module | record, experiment and structure are different units",fontweight="bold",fontsize=12)
    rounded(ax,.04,.53,.16,.19,"Database records\nChEMBL · PubChem\nBindingDB · PDBe",LTBLUE,fs=7);arrow(ax,.21,.62,.30,.62);rounded(ax,.31,.53,.16,.19,"Mirror control\nsource record IDs\npublication/assay keys",PALEGREEN,fs=7);arrow(ax,.48,.62,.57,.62);rounded(ax,.58,.53,.16,.19,"Evidence tiers\nBE1 · BE2 · BE3\ndirectness + context",YELLOW,fs=7);arrow(ax,.75,.62,.82,.62);rounded(ax,.83,.53,.14,.19,"Pair summary\nprovenance retained",PEACH,fs=7)
    ax.text(.50,.30,"Contributing database count  ≠  putative independent experiment count",ha="center",fontsize=9,fontweight="bold",color=RED)
    save(fig,"Section_Graphic_Evidence")


def section_disease():
    fig,ax=plt.subplots(figsize=(8,4.5));ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis("off");ax.text(.03,.94,"Disease module | canonical identity, hierarchy and multi-label anatomy",fontweight="bold",fontsize=12)
    rounded(ax,.04,.50,.16,.20,"Source diseases\nOMIM · Orphanet\nEFO · UniProt",LTBLUE,fs=7);arrow(ax,.21,.60,.31,.60);rounded(ax,.32,.50,.16,.20,"MONDO exact-only\ncanonical identity",PALEGREEN,fs=7,weight="bold");arrow(ax,.49,.60,.57,.60)
    rounded(ax,.58,.66,.17,.13,"Therapeutic areas",YELLOW,fs=7);rounded(ax,.58,.47,.17,.13,"Disease hierarchy",PEACH,fs=7);rounded(ax,.58,.28,.17,.13,"DO/Uberon anatomy",SALMON,fs=7)
    for y in [.725,.535,.345]:arrow(ax,.75,y,.83,.55)
    ax.add_patch(patches.Circle((.90,.55),.09,fc=LIGHT,ec=RED,lw=1.2));ax.text(.90,.55,"One disease\nmay map to\nmultiple systems",ha="center",va="center",fontsize=6.5,fontweight="bold")
    save(fig,"Section_Graphic_Disease")


def section_site_complex():
    fig,ax=plt.subplots(figsize=(8,4.5));ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis("off");ax.text(.03,.94,"Site and complex module | coordinate usability and biological context are separate",fontweight="bold",fontsize=12)
    bilayer(ax,.04,.48,.36,.10);protein_helix(ax,.18,.57,1.5,BLUE);molecule_icon(ax,.23,.60,.035,ORANGE);ax.add_patch(patches.Circle((.23,.60),.075,fill=False,ec=RED,lw=1,ls="--"));ax.text(.20,.33,"Residues → chain → coordinates\nS1 / S2 / S3",ha="center",fontweight="bold")
    arrow(ax,.41,.56,.49,.56)
    for i,(x,c) in enumerate([(.58,BLUE),(.68,GREEN),(.63,ORANGE)]):ax.add_patch(patches.Ellipse((x,.58+(i-1)*.05),.14,.10,fc=c,alpha=.65,ec="white"))
    molecule_icon(ax,.64,.58,.03,RED);ax.text(.63,.33,"Complex entity\ncomponents · stoichiometry · state",ha="center",fontweight="bold")
    arrow(ax,.73,.56,.80,.56);rounded(ax,.81,.43,.16,.25,"complex-specific\nsource-asserted context\ncomponent-context-only",PEACH,fs=7,weight="bold")
    save(fig,"Section_Graphic_Site_Complex")


def section_integrated():
    fig,ax=plt.subplots(figsize=(8,4.5));ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis("off");ax.text(.03,.94,"Integrated analysis | transparent dimensions support follow-up",fontweight="bold",fontsize=12)
    inputs=[("Disease\nsupport",LTBLUE),("Anatomical\nconcordance",PALEGREEN),("Chemical\ncoverage",YELLOW),("Structural\nreadiness",PEACH)]
    for i,(txt,col) in enumerate(inputs):rounded(ax,.04,.73-i*.17,.18,.11,txt,col,fs=7,weight="bold");arrow(ax,.23,.785-i*.17,.39,.57)
    ax.add_patch(patches.FancyBboxPatch((.40,.38),.20,.35,boxstyle="round,pad=.018",fc=LIGHT,ec=BLUE,lw=1.2));ax.text(.50,.555,"Transparent\n2D/Pareto landscape",ha="center",va="center",fontweight="bold",fontsize=9)
    arrow(ax,.61,.56,.70,.56);rounded(ax,.71,.46,.24,.20,"Hypothesis-generation subset\n→ literature\n→ assay\n→ docking",SALMON,fs=8,weight="bold")
    ax.text(.50,.16,"Prioritisation is not causal inference or binding validation",ha="center",color=RED,fontweight="bold")
    save(fig,"Section_Graphic_Integrated")


def main():
    OUT.mkdir(parents=True,exist_ok=True);graphical_abstract();section_protein();section_compound();section_evidence();section_disease();section_site_complex();section_integrated();print(json.dumps({"graphical_abstract":1,"section_graphics":6,"formats":["png","svg","pdf"]}))


if __name__=="__main__":main()
