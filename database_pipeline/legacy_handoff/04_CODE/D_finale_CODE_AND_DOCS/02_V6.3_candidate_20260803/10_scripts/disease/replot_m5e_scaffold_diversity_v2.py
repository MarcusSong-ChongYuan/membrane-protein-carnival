#!/usr/bin/env python3
"""Replot the completed full-release scaffold analysis with ASCII-safe labels."""

from __future__ import annotations

import argparse, json
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

C={"navy":"#153B55","blue":"#3B73A3","teal":"#2A9D8F","orange":"#F4A261","gray":"#AEB7C1","light":"#E8EDF1","dark":"#5F6B75","purple":"#7A5195"}
def ar():
 p=argparse.ArgumentParser(); p.add_argument('--frequency',type=Path,required=True); p.add_argument('--summary',type=Path,required=True); p.add_argument('--figure-dir',type=Path,required=True); p.add_argument('--qa-dir',type=Path,required=True); return p.parse_args()
def main():
 a=ar(); s=pd.read_csv(a.summary,sep='\t').iloc[0]; f=pd.read_csv(a.frequency,sep='\t'); mpl.rcParams.update({"font.family":"Arial","font.size":8,"axes.titlesize":9,"axes.titleweight":"bold","axes.labelsize":8,"xtick.labelsize":7,"ytick.labelsize":7,"axes.spines.top":False,"axes.spines.right":False,"figure.facecolor":"white","savefig.facecolor":"white","pdf.fonttype":42,"svg.fonttype":"none"})
 fig=plt.figure(figsize=(11.4,7.4)); gs=fig.add_gridspec(2,2,hspace=.42,wspace=.28); axes=[fig.add_subplot(gs[i]) for i in range(4)]
 for ax,l,t in zip(axes,'ABCD',['Structure-resolution coverage','Scaffold rank-frequency','Cumulative chemical coverage','Most frequent Murcko scaffolds']): ax.text(-.08,1.06,l,transform=ax.transAxes,fontsize=11,fontweight='bold',color=C['navy']); ax.set_title(t,loc='left',pad=7)
 total=int(s.core_qc_compounds); vals=[int(s.scaffold_bearing_compounds),int(s.acyclic_compounds),int(s.missing_smiles),int(s.invalid_smiles)]; labs=['Scaffold-bearing','Acyclic','Missing SMILES','Invalid SMILES']; axes[0].barh(range(4),vals,color=[C['teal'],C['orange'],C['gray'],C['purple']]); axes[0].set_yticks(range(4),labs); axes[0].invert_yaxis(); axes[0].set_xlabel('Canonical compounds (count)'); axes[0].grid(axis='x',color=C['light'],lw=.6)
 for i,v in enumerate(vals): axes[0].text(v,i,f' {v:,} ({v/total:.1%})',va='center',fontsize=6.8)
 counts=f.compound_count.to_numpy(float); ranks=np.arange(1,len(f)+1); axes[1].plot(ranks,counts,color=C['blue'],lw=1.2); axes[1].set_xscale('log'); axes[1].set_yscale('log'); axes[1].set_xlabel('Scaffold rank'); axes[1].set_ylabel('Compounds per scaffold'); axes[1].grid(which='both',color=C['light'],lw=.5); axes[1].text(.98,.96,f"Unique scaffolds: {int(s.unique_murcko_scaffolds):,}\nSingletons: {int(s.singleton_scaffolds):,} ({s.singleton_scaffold_fraction:.1%})",transform=axes[1].transAxes,ha='right',va='top',fontsize=7,color=C['dark'])
 cumulative=np.cumsum(counts)/counts.sum()*100; axes[2].plot(ranks,cumulative,color=C['teal'],lw=1.6); axes[2].set_xscale('log'); axes[2].set_ylim(0,101); axes[2].set_xlabel('Top-ranked scaffolds included'); axes[2].set_ylabel('Scaffold-bearing compounds covered (%)'); axes[2].grid(which='both',color=C['light'],lw=.5); axes[2].scatter([10,100],[s.top10_scaffold_compound_share*100,s.top100_scaffold_compound_share*100],color=[C['orange'],C['purple']],zorder=3); axes[2].text(10,s.top10_scaffold_compound_share*100,f'  Top 10: {s.top10_scaffold_compound_share:.1%}',fontsize=7,va='center'); axes[2].text(100,s.top100_scaffold_compound_share*100,f'  Top 100: {s.top100_scaffold_compound_share:.1%}',fontsize=7,va='center')
 top=f.head(12); mols=[Chem.MolFromSmiles(x) for x in top.murcko_scaffold_smiles]; legends=[f'#{i+1}  n={int(v):,}' for i,v in enumerate(top.compound_count)]; grid=Draw.MolsToGridImage(mols,molsPerRow=4,subImgSize=(180,125),legends=legends,useSVG=False); axes[3].imshow(np.asarray(grid)); axes[3].axis('off'); axes[3].text(0,-.06,'Connectivity frameworks; counts are not activity or approval rankings.',transform=axes[3].transAxes,fontsize=6.7,color=C['dark'],va='top')
 fig.suptitle('M5E - Full-release Bemis-Murcko scaffold diversity',x=.06,ha='left',fontsize=12,fontweight='bold',color=C['navy']); fig.text(.06,.93,f'All core, QC-passed canonical compounds were processed (n={total:,}); no embedding or random subsampling was used.',fontsize=7.5,color=C['dark']); fig.subplots_adjust(left=.10,right=.98,top=.86,bottom=.08)
 paths={}
 for ext in ('png','svg','pdf'):
  d=a.figure_dir/ext; d.mkdir(parents=True,exist_ok=True); p=d/f'M5E_scaffold_diversity.{ext}'; fig.savefig(p,dpi=450 if ext=='png' else None,bbox_inches='tight'); paths[ext]=str(p)
 plt.close(fig)
 (a.figure_dir/'CAPTION_M5E_SCAFFOLD_DIVERSITY.md').write_text("""# M5E caption\n\n**M5E. Full-release Bemis-Murcko scaffold diversity of canonical compounds.** All core, QC-passed canonical compounds in the frozen V6.2 small-molecule master were processed without random sampling. Ring-containing molecules were reduced to Bemis-Murcko frameworks with RDKit; acyclic molecules, missing SMILES and invalid SMILES are reported separately rather than assigned artificial scaffolds. The rank-frequency and cumulative-coverage panels quantify scaffold reuse, while the structure gallery shows the twelve most frequent frameworks. Scaffold frequency describes chemical redundancy and diversity; it is not a measure of bioactivity, approval status or docking suitability.\n""",encoding='utf-8')
 out={'status':'PASS','module_version':'0.2.0-plot-refresh','data_unchanged':True,'figure_paths':paths,'counts':{k:(v.item() if hasattr(v,'item') else v) for k,v in s.to_dict().items()}}; (a.qa_dir/'M5E_SCAFFOLD_DIVERSITY_V0_2_PLOT_VALIDATION.json').write_text(json.dumps(out,indent=2),encoding='utf-8'); print(json.dumps(out,indent=2))
if __name__=='__main__': main()
