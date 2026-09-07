"""Visual-only v2 redraw of the frozen C2 compound-source overlap tables."""
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, LogLocator
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results'/'compound_figure_pool'
I=OUT/'C2_compound_source_intersections.tsv';M=OUT/'C2_compound_source_membership.tsv';S=OUT/'C2_source_overlap_summary.tsv'
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans','sans-serif'],'svg.fonttype':'none','pdf.fonttype':42,'figure.facecolor':'white','savefig.facecolor':'white','axes.linewidth':.8})
INK,SLATE,NAVY,TEAL,LIGHT='#24364B','#5E7182','#173B6C','#168C8C','#DCE4E9'
SOURCES=['PubChem BioAssay','ChEMBL','BindingDB','PDBe','IUPHAR/BPS Guide to PHARMACOLOGY','PDBbind','BRENDA','sc-PDB']
LABEL={'IUPHAR/BPS Guide to PHARMACOLOGY':'GtoPdb (IUPHAR/BPS)'}
def fmt(n):return f'{n/1000:.1f}k' if n>=1000 else str(int(n))
def main():
 inter=pd.read_csv(I,sep='\t');membership=pd.read_csv(M,sep='\t');summary=pd.read_csv(S,sep='\t');N=len(membership)
 # Exclude the deliberately aggregated low-frequency node from the display,
 # not from the analysis tables or denominator.
 eligible=inter.loc[~inter.source_combination.str.contains('Other sources',regex=False)].copy().sort_values('canonical_compound_count',ascending=False)
 top=eligible.head(12).copy(); single=eligible.loc[eligible.source_count_display.eq(1)&eligible.source_combination.isin(SOURCES)]
 extra=single.loc[~single.source_combination.isin(top.source_combination)].head(2)
 shown=pd.concat([top,extra]).drop_duplicates('source_combination').sort_values('canonical_compound_count',ascending=False).reset_index(drop=True)
 source_counts=membership.drop(columns='compound_internal_id').sum().reindex(SOURCES).fillna(0).astype(int)
 fig=plt.figure(figsize=(1800/150,1200/150),dpi=600);axbar=fig.add_axes([.38,.63,.53,.23]);axmat=fig.add_axes([.38,.19,.53,.34]);axset=fig.add_axes([.07,.19,.17,.34]);x=np.arange(len(shown));singlemask=shown.source_count_display.eq(1)
 axbar.bar(x,shown.canonical_compound_count,color=np.where(singlemask,TEAL,NAVY),width=.66);axbar.set_yscale('log');axbar.set(xticks=[],ylabel='Intersection size\n(log scale)',ylim=(max(1,shown.canonical_compound_count.min()*.55),shown.canonical_compound_count.max()*2.5));axbar.yaxis.set_major_locator(LogLocator(base=10));axbar.yaxis.set_major_formatter(FuncFormatter(lambda v,_:fmt(v)));axbar.tick_params(labelsize=8);axbar.spines[['top','right','bottom']].set_visible(False);axbar.grid(axis='y',which='major',color='#E5EAED',lw=.65);axbar.set_axisbelow(True)
 for i,r in shown.head(6).iterrows():axbar.text(i,r.canonical_compound_count*1.23,fmt(r.canonical_compound_count),ha='center',va='bottom',fontsize=8,color=INK)
 y=np.arange(len(SOURCES));axmat.set(xticks=x,yticks=y,yticklabels=[LABEL.get(s,s) for s in SOURCES],xlim=(-.6,len(shown)-.4),ylim=(len(y)-.5,-.5));axmat.tick_params(axis='y',labelsize=8,length=0);axmat.tick_params(axis='x',bottom=False,labelbottom=False);axmat.spines[:].set_visible(False)
 for i,r in shown.iterrows():
  inc=r.source_combination.split('|');yy=[SOURCES.index(s) for s in inc if s in SOURCES];axmat.scatter([i]*len(y),y,s=30,color=LIGHT,zorder=1);axmat.scatter([i]*len(yy),yy,s=62,color=TEAL if len(yy)==1 else NAVY,zorder=3)
  if len(yy)>1:axmat.plot([i,i],[min(yy),max(yy)],color=NAVY,lw=1.65,zorder=2)
 axset.barh(y,source_counts.values,color=NAVY,height=.54);axset.set_xscale('log');axset.set(yticks=y,yticklabels=['']*len(y),xlabel='Source coverage\n(log scale)');axset.xaxis.set_major_locator(LogLocator(base=10));axset.xaxis.set_major_formatter(FuncFormatter(lambda v,_:fmt(v)));axset.invert_yaxis();axset.tick_params(axis='x',labelsize=7);axset.spines[['top','right','left']].set_visible(False)
 for i,v in enumerate(source_counts.values):axset.text(v*1.12,i,fmt(v),va='center',fontsize=7.2,color=INK)
 axbar.text(0,1.28,'Compound-source overlap and complementarity',transform=axbar.transAxes,fontsize=13,weight='bold',color=INK);axbar.text(0,1.12,'Canonical compounds counted once per source; top source intersections shown.',transform=axbar.transAxes,fontsize=8.5,color=SLATE);axbar.text(0,1.00,f'Analysis N = {N:,} canonical compounds with direct source provenance',transform=axbar.transAxes,fontsize=8,color=SLATE);fig.text(.38,.095,'Teal = source-only intersection; complete intersections are provided in source data.',fontsize=7.5,color=SLATE)
 for s,k in {'.svg':{},'.pdf':{},'.png':{'dpi':600}}.items():fig.savefig(OUT/f'C2_compound_source_upset_v2{s}',bbox_inches='tight',**k)
 plt.close(fig)
if __name__=='__main__':main()
