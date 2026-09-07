from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results'/'compound_figure_pool';M=Path(r'D:\finale\FORMAL\01_core_v72\01_release_tables\compound_master_v72.tsv.gz')
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans','sans-serif'],'svg.fonttype':'none','pdf.fonttype':42,'figure.facecolor':'white','savefig.facecolor':'white'})
def main():
 d=pd.read_csv(M,sep='\t',compression='gzip',usecols=['compound_internal_id','compound_scope_class','formal_charge']);d=d[d.compound_scope_class.isin(['core_small_molecule','structure_resolved_small_molecule'])];q=pd.to_numeric(d.formal_charge,errors='coerce');lab=pd.Series('Unresolved',index=d.index);lab[q.eq(0)]='Neutral recorded structure';lab[q.gt(0)]='Positively charged recorded structure';lab[q.lt(0)]='Negatively charged recorded structure';tab=lab.value_counts().rename_axis('charge_class').reset_index(name='canonical_compound_count');tab['fraction']=tab.canonical_compound_count/len(d);tab.to_csv(OUT/'C6_charge_class_distribution.tsv',sep='\t',index=False);cols=['#7B95C6','#F59C7C','#49C2D9','#A8B1B8'];fig=plt.figure(figsize=(178/25.4,78/25.4),dpi=600);ax=fig.add_axes([.24,.14,.52,.72]);wedges,_=ax.pie(tab.canonical_compound_count,colors=cols[:len(tab)],startangle=90,wedgeprops={'width':.38,'edgecolor':'white','linewidth':1});ax.text(0,0,'Canonical\nsmall molecules\n' + f'n = {len(d):,}',ha='center',va='center',fontsize=8,color='#24364B');
 for w,r in zip(wedges,tab.itertuples(index=False)):
  ang=(w.theta1+w.theta2)/2; x,y=__import__('numpy').cos(__import__('numpy').deg2rad(ang)),__import__('numpy').sin(__import__('numpy').deg2rad(ang));ax.annotate(f'{r.charge_class}\n{r.canonical_compound_count:,} ({r.fraction:.1%})',xy=(x*.82,y*.82),xytext=(x*1.27,y*1.27),ha='left' if x>0 else 'right',va='center',fontsize=6.7,arrowprops={'arrowstyle':'-','lw':.6,'color':'#5E7182'})
 ax.text(.5,-.12,'Formal charge of recorded standardized structure; not physiological charge at pH 7.4.',transform=ax.transAxes,ha='center',fontsize=6.5,color='#5E7182');
 for s,k in {'.svg':{},'.pdf':{},'.png':{'dpi':600}}.items():fig.savefig(OUT/f'C6_charge_class_donut{s}',bbox_inches='tight',**k)
 plt.close(fig)
if __name__=='__main__':main()
