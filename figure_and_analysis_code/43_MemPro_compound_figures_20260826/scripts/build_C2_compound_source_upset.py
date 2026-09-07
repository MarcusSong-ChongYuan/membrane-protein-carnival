"""Canonical-compound provenance overlap, using direct formal source membership only."""
from pathlib import Path
import json,re,itertools
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results'/'compound_figure_pool';M=Path(r'D:\finale\FORMAL\01_core_v72\01_release_tables\compound_master_v72.tsv.gz')
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans','sans-serif'],'svg.fonttype':'none','pdf.fonttype':42,'figure.facecolor':'white','savefig.facecolor':'white','axes.linewidth':.8})
INK,SLATE,NAVY,TEAL='#24364B','#5E7182','#173B6C','#168C8C'
ALIASES={'CHEMBL':'ChEMBL','ChEMBL database':'ChEMBL','PubChem':'PubChem BioAssay','IUPHAR/BPS Guide to Pharmacology':'IUPHAR/BPS Guide to PHARMACOLOGY'}
def split(v):return [ALIASES.get(x.strip(),x.strip()) for x in re.split(r'[;|]',str(v)) if x.strip() and x.strip().lower()!='nan']
def main():
 OUT.mkdir(parents=True,exist_ok=True);d=pd.read_csv(M,sep='\t',compression='gzip',usecols=['compound_internal_id','source_databases'],low_memory=False).drop_duplicates('compound_internal_id');mem=[]
 for r in d.itertuples(index=False):
  for s in set(split(r.source_databases)):mem.append((r.compound_internal_id,s,1))
 long=pd.DataFrame(mem,columns=['compound_internal_id','source_database','direct_source_membership']).drop_duplicates();universe=long.compound_internal_id.nunique();source_counts=long.groupby('source_database').compound_internal_id.nunique().sort_values(ascending=False);top=source_counts.head(8).index.tolist();other=set(source_counts.index)-set(top)
 def combo(g):
  ss=set(g.source_database);show=[s for s in top if s in ss];
  if ss&other:show.append('Other sources')
  return '|'.join(show)
 combos=long.groupby('compound_internal_id').apply(combo,include_groups=False).rename('source_combination').reset_index(); inter=combos.source_combination.value_counts().rename_axis('source_combination').reset_index(name='canonical_compound_count');inter['fraction']=inter.canonical_compound_count/universe;inter['source_count_display']=inter.source_combination.str.count(r'\|')+1;inter.to_csv(OUT/'C2_compound_source_intersections.tsv',sep='\t',index=False)
 mat=pd.crosstab(long.compound_internal_id,long.source_database).astype(int);mat.insert(0,'compound_internal_id',mat.index);mat.to_csv(OUT/'C2_compound_source_membership.tsv',sep='\t',index=False)
 pairs=[]
 for a,b in itertools.combinations(source_counts.index,2):
  ia=((mat[a].to_numpy()>0)&(mat[b].to_numpy()>0)).sum();u=((mat[a].to_numpy()>0)|(mat[b].to_numpy()>0)).sum();pairs.append({'source_a':a,'source_b':b,'intersection_count':int(ia),'jaccard_index':float(ia/u) if u else np.nan})
 pd.DataFrame(pairs).to_csv(OUT/'C2_source_pairwise_overlap.tsv',sep='\t',index=False)
 nsrc=long.groupby('compound_internal_id').source_database.nunique();summ=pd.DataFrame({'source_support_class':['exactly 1','exactly 2','exactly 3','≥4'],'canonical_compound_count':[(nsrc==1).sum(),(nsrc==2).sum(),(nsrc==3).sum(),(nsrc>=4).sum()]});summ['fraction']=summ.canonical_compound_count/universe;summ.to_csv(OUT/'C2_source_overlap_summary.tsv',sep='\t',index=False)
 shown=inter.head(20).copy();fig=plt.figure(figsize=(178/25.4,112/25.4),dpi=600);axbar=fig.add_axes([.46,.62,.45,.25]);axmat=fig.add_axes([.46,.22,.45,.32]);axset=fig.add_axes([.04,.22,.12,.32]);x=np.arange(len(shown));singleton=shown.source_count_display.eq(1);axbar.bar(x,shown.canonical_compound_count,color=np.where(singleton,TEAL,NAVY),width=.72);axbar.set(xticks=[],ylabel='Intersection\nsize');axbar.tick_params(labelsize=7);axbar.spines[['top','right','bottom']].set_visible(False)
 for i,r in shown.iterrows():axbar.text(i, r.canonical_compound_count+max(1,shown.canonical_compound_count.max()*.015),f'{r.canonical_compound_count/1000:.1f}k' if r.canonical_compound_count>=1000 else str(r.canonical_compound_count),ha='center',fontsize=6.2,color=INK)
 display_sources=top+['Other sources'];y=np.arange(len(display_sources));axmat.set(xticks=x,yticks=y,yticklabels=display_sources,xlim=(-.6,len(shown)-.4),ylim=(len(y)-.5,-.5));axmat.tick_params(axis='y',labelsize=7,length=0);axmat.tick_params(axis='x',bottom=False,labelbottom=False);axmat.spines[:].set_visible(False)
 for i,r in shown.iterrows():
  ss=r.source_combination.split('|');yy=[display_sources.index(s) for s in ss if s in display_sources];axmat.scatter([i]*len(y),y,s=12,color='#D9E0E4',zorder=1);axmat.scatter([i]*len(yy),yy,s=32,color=TEAL if len(yy)==1 else NAVY,zorder=2); 
  if len(yy)>1:axmat.plot([i,i],[min(yy),max(yy)],color=NAVY,lw=1.1,zorder=1)
 sc=source_counts.reindex(top);axset.barh(np.arange(len(top)),sc.values,color=NAVY,height=.55);axset.set(yticks=np.arange(len(top)),yticklabels=['']*len(top),xlabel='Set size');axset.invert_yaxis();axset.tick_params(axis='x',labelsize=6.5);axset.spines[['top','right','left']].set_visible(False)
 for i,v in enumerate(sc.values):axset.text(v,i,f'{v/1000:.1f}k\n({v/universe:.1%})',va='center',fontsize=6.2,color=INK)
 axbar.text(0,1.2,'Compound-source overlap and complementarity',transform=axbar.transAxes,fontsize=9,weight='bold',color=INK);axbar.text(0,1.06,f'Analysis N = {universe:,} canonical compounds with direct formal source provenance; top 20 intersections shown',transform=axbar.transAxes,fontsize=6.8,color=SLATE);fig.text(.29,.12,'Filled dot: source included; teal bars/dots: singleton-source intersections. Canonical compounds are counted once per source.',fontsize=6.6,color=SLATE)
 for s,k in {'.svg':{},'.pdf':{},'.png':{'dpi':600}}.items():fig.savefig(OUT/f'C2_compound_source_upset{s}',bbox_inches='tight',**k)
 plt.close(fig);(OUT/'C2_source_overlap_manifest.json').write_text(json.dumps({'universe':'all formal canonical compounds with direct source provenance','N':int(universe),'sources_top8':top,'source_field':'source_databases; direct membership only'},indent=2),encoding='utf-8')
if __name__=='__main__':main()
