from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results'/'compound_figure_pool';SRC=OUT/'C3_scaffold_rank_frequency.tsv'
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans','sans-serif'],'svg.fonttype':'none','pdf.fonttype':42,'figure.facecolor':'white','savefig.facecolor':'white'})
def main():
 d=pd.read_csv(SRC,sep='\t');d[['scaffold_rank','cumulative_compound_coverage_fraction']].to_csv(OUT/'C5_scaffold_cumulative_coverage.tsv',sep='\t',index=False);fig=plt.figure(figsize=(178/25.4,82/25.4),dpi=600);ax=fig.add_axes([.15,.2,.7,.65]);ax.plot(d.scaffold_rank,d.cumulative_compound_coverage_fraction*100,color='#168C8C',lw=1.5);ax.set_xscale('log');ax.set(xlabel='Number of exact Bemis–Murcko scaffolds',ylabel='Cumulative fraction of compounds (%)',ylim=(0,102),xlim=(1,len(d)));ax.tick_params(axis='both',labelsize=8);ax.grid(color='#E4EAED',lw=.65);ax.set_axisbelow(True);ax.spines[['top','right']].set_visible(False)
 for pct in [25,50,75,80,90]:
  r=d.loc[d.cumulative_compound_coverage_fraction.ge(pct/100)].iloc[0];ax.scatter(r.scaffold_rank,pct,s=27,color='#D6A23A',edgecolor='white',lw=.7,zorder=3);ax.annotate(f'{pct}%: {int(r.scaffold_rank):,}',xy=(r.scaffold_rank,pct),xytext=(r.scaffold_rank*1.18,pct-8),fontsize=7,arrowprops={'arrowstyle':'-','lw':.6,'color':'#5E7182'})
 ax.text(0,1.04,f'Analysis N = {int(d.cumulative_unique_compounds.iloc[-1]):,} scaffold-bearing canonical small molecules',transform=ax.transAxes,fontsize=7.2,color='#5E7182')
 for s,k in {'.svg':{},'.pdf':{},'.png':{'dpi':600}}.items():fig.savefig(OUT/f'C5_scaffold_cumulative_coverage{s}',bbox_inches='tight',**k)
 plt.close(fig)
if __name__=='__main__':main()
