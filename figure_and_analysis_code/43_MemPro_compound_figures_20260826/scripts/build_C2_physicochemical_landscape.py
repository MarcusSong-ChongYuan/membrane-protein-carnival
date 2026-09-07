"""C2: Compact robust-scaled physicochemical distribution landscape."""
from __future__ import annotations
from pathlib import Path
import json
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results'/'compound_figure_pool';MASTER=Path(r'D:\finale\FORMAL\01_core_v72\01_release_tables\compound_master_v72.tsv.gz')
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans','sans-serif'],'svg.fonttype':'none','pdf.fonttype':42,'figure.facecolor':'white','savefig.facecolor':'white','axes.linewidth':.8})
INK,SLATE='#24364B','#5E7182';COLS=['#173B6C','#168C8C','#6F9B6D','#887AA9','#D6A23A','#B77B8C']
PROPS=[('molecular_weight','Molecular weight (Da)'),('xlogp','XlogP'),('tpsa','TPSA (Å²)'),('hbond_donor_count','H-bond donors'),('hbond_acceptor_count','H-bond acceptors'),('rotatable_bond_count','Rotatable bonds')]
def main():
 OUT.mkdir(parents=True,exist_ok=True);use=['compound_internal_id','compound_scope_class']+[p[0] for p in PROPS];d=pd.read_csv(MASTER,sep='\t',compression='gzip',usecols=use,low_memory=False);d=d[d.compound_scope_class.isin(['core_small_molecule','structure_resolved_small_molecule'])].copy();d.to_csv(OUT/'C2_physicochemical_raw_source.tsv',sep='\t',index=False)
 summary=[];plot=[]
 for key,label in PROPS:
  v=pd.to_numeric(d[key],errors='coerce').replace([np.inf,-np.inf],np.nan).dropna().to_numpy();q5,q25,med,q75,q95=np.quantile(v,[.05,.25,.5,.75,.95]);iqr=q75-q25
  summary.append({'property':key,'display_label':label,'N':len(v),'missing_count':len(d)-len(v),'missing_fraction':1-len(v)/len(d),'median':med,'IQR':iqr,'mean':v.mean(),'p5':q5,'p95':q95,'minimum':v.min(),'maximum':v.max(),'plot_transform':'robust z=(x-median)/IQR; values outside [-4,4] binned into edge bins only for visualization'})
  z=(v-med)/(iqr if iqr>0 else 1);z=np.clip(z,-4,4);h,e=np.histogram(z,bins=120,range=(-4,4),density=True);plot.append((label,(e[:-1]+e[1:])/2,h))
 s=pd.DataFrame(summary);s.to_csv(OUT/'C2_physicochemical_summary.tsv',sep='\t',index=False)
 fig=plt.figure(figsize=(178/25.4,100/25.4),dpi=600);ax=fig.add_axes([.16,.17,.70,.70])
 for i,(label,x,h) in enumerate(plot[::-1]):
  y=i+h/h.max()*.72;c=COLS[len(plot)-1-i];ax.fill_between(x,i,y,color=c,alpha=.52,lw=0);ax.plot(x,y,color=c,lw=1);ax.text(-4.25,i+.23,label,ha='right',va='center',fontsize=7,color=INK)
 ax.axvline(0,color='#9AA8B1',lw=.8,ls='--');ax.set(yticks=[],xlim=(-4.4,4.1),xlabel='Robust-scaled property value: (value − median) / IQR');ax.tick_params(axis='x',labelsize=8,colors=INK);ax.spines[['top','right','left']].set_visible(False);ax.grid(axis='x',color='#E5EAED',lw=.6);ax.set_axisbelow(True)
 ax.text(0,1.03,f'Analysis N = {len(d):,} formal small molecules; distributions use all valid values',transform=ax.transAxes,fontsize=7.2,color=SLATE);ax.text(0,-.26,'Rows are aligned by robust scale only to compare distribution shape. Original units and untrimmed values are retained in Source Data.',transform=ax.transAxes,fontsize=6.8,color=SLATE)
 for sfx,k in {'.svg':{},'.pdf':{},'.png':{'dpi':600}}.items():fig.savefig(OUT/f'C2_physicochemical_landscape{sfx}',bbox_inches='tight',**k)
 plt.close(fig);(OUT/'C2_physchem_qc.txt').write_text(f'Population\t{len(d):,} existing formal small molecules\nInvalid-value rule\tNaN and ±Inf excluded only property-wise; no numerical outliers removed\n',encoding='utf-8');(OUT/'C2_reproducibility_manifest.json').write_text(json.dumps({'figure':'C2_physicochemical_landscape','input':str(MASTER),'analysis_unit':'one canonical compound'},indent=2),encoding='utf-8')
if __name__=='__main__':main()
