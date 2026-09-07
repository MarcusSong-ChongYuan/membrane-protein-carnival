#!/usr/bin/env python3
"""Generate evidence-lineage and docking-priority figures for V6.3."""

from __future__ import annotations

import argparse, json
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

C={"navy":"#153B55","blue":"#3B73A3","teal":"#2A9D8F","orange":"#F4A261","red":"#C95B5B","purple":"#7A5195","gray":"#AEB7C1","light":"#E8EDF1","dark":"#5F6B75","green":"#6A994E"}

def ar():
 p=argparse.ArgumentParser(); p.add_argument('--pair',type=Path,required=True); p.add_argument('--review-sample',type=Path,required=True); p.add_argument('--docking',type=Path,required=True); p.add_argument('--validation',type=Path,required=True); p.add_argument('--figure-dir',type=Path,required=True); p.add_argument('--qa-dir',type=Path,required=True); return p.parse_args()
def setup(): mpl.rcParams.update({"font.family":"Arial","font.size":8,"axes.titlesize":9,"axes.titleweight":"bold","axes.labelsize":8,"xtick.labelsize":7,"ytick.labelsize":7,"axes.spines.top":False,"axes.spines.right":False,"figure.facecolor":"white","savefig.facecolor":"white","pdf.fonttype":42,"svg.fonttype":"none"})
def pan(ax,l,t): ax.text(-.10,1.06,l,transform=ax.transAxes,fontsize=11,fontweight='bold',color=C['navy']); ax.set_title(t,loc='left',pad=7)
def save(fig,base,stem):
 out={}
 for ext in ('png','svg','pdf'):
  d=base/ext; d.mkdir(parents=True,exist_ok=True); p=d/f'{stem}.{ext}'; fig.savefig(p,dpi=450 if ext=='png' else None,bbox_inches='tight'); out[ext]=str(p)
 plt.close(fig); return out

def m7(a,pair,sample,val):
 fig=plt.figure(figsize=(11.6,7.6)); gs=fig.add_gridspec(2,2,hspace=.44,wspace=.32); aa,ab,ac,ad=[fig.add_subplot(gs[i]) for i in range(4)]
 pan(aa,'A','Database-count semantics audit')
 audit=pair.legacy_source_count_semantics_audit_v63.value_counts(); labels=['Legacy count matches','Legacy count differs']; vals=[audit.get('legacy_count_matches_database_count',0),audit.get('legacy_count_requires_review',0)]
 aa.barh(range(2),vals,color=[C['teal'],C['orange']]); aa.set_yticks(range(2),labels); aa.invert_yaxis(); aa.set_xlabel('Protein–compound pairs'); aa.grid(axis='x',color=C['light'],lw=.6)
 for i,v in enumerate(vals): aa.text(v,i,f' {v:,} ({v/len(pair):.1%})',va='center',fontsize=7)
 pan(ab,'B','Distinct lineage keys by evidence domain')
 counts=val['counts']; names=['Structure\n(evidence)','Structure\n(site/chain)','PubChem\nassay','Literature\nproxy']; vals=[counts['structure_lineage_keys'],counts['structure_site_lineage_keys'],counts['pubchem_lineage_keys'],counts['literature_lineage_keys']]
 ab.bar(range(4),vals,color=[C['blue'],C['navy'],C['purple'],C['green']]); ab.set_xticks(range(4),names); ab.set_yscale('log'); ab.set_ylabel('Distinct lineage keys (log scale)'); ab.grid(axis='y',which='both',color=C['light'],lw=.6)
 for i,v in enumerate(vals): ab.text(i,v*1.08,f'{v:,}',ha='center',fontsize=7)
 pan(ac,'C','Pair-level lineage coverage')
 cov={"Structure":(pd.to_numeric(pair.independent_structure_count_v63)>0).mean(),"PubChem assay":(pd.to_numeric(pair.independent_pubchem_assay_count_v63)>0).mean(),"Literature proxy":(pd.to_numeric(pair.independent_literature_experiment_count_v63)>0).mean(),"≥2 modalities":(pd.to_numeric(pair.independent_evidence_modality_count_v63)>=2).mean()}
 ac.barh(range(4),np.array(list(cov.values()))*100,color=[C['blue'],C['purple'],C['green'],C['teal']]); ac.set_yticks(range(4),list(cov)); ac.invert_yaxis(); ac.set_xlim(0,105); ac.set_xlabel('Protein–compound pairs (%)'); ac.grid(axis='x',color=C['light'],lw=.6)
 for i,v in enumerate(cov.values()): ac.text(v*100,i,f' {v:.1%}',va='center',fontsize=7)
 pan(ad,'D','Bounded validation sample, not a claim of full manual review')
 pri=sample.review_priority.value_counts().reindex(['P1','P2','P3']).fillna(0)
 ad.bar(range(3),pri.values,color=[C['red'],C['orange'],C['gray']]); ad.set_xticks(range(3),['P1 conflict','P2 BE1/docking','P3 BE2']); ad.set_ylabel('Sampled evidence rows'); ad.grid(axis='y',color=C['light'],lw=.6)
 for i,v in enumerate(pri.values): ad.text(i,v,f'{int(v):,}',ha='center',va='bottom',fontsize=7)
 ad.text(.02,.96,f"Audit universe: {counts['review_queue_rows']:,} machine-selected rows\nDeterministic stratified sample: {counts['review_sample_rows']:,} rows",transform=ad.transAxes,va='top',fontsize=7,color=C['dark'])
 fig.suptitle('M7 — Evidence provenance is separated from database contribution counts',x=.06,ha='left',fontsize=12,fontweight='bold',color=C['navy']); fig.text(.06,.93,'Structure, PubChem assay, literature and evidence-modality keys are counted independently; database labels are never presented as independent experiments.',ha='left',fontsize=7.5,color=C['dark']); fig.subplots_adjust(left=.13,right=.98,top=.86,bottom=.09)
 paths=save(fig,a.figure_dir,'M7_evidence_lineage_v63')
 (a.figure_dir/'CAPTION_M7_EVIDENCE_LINEAGE.md').write_text("""# M7 caption\n\n**M7. Evidence provenance and source-count semantics.** The legacy field named `independent_source_count` is audited against a literal count of distinct contributing database labels (A); mismatches are retained for review and are not interpreted as biological disagreement. Structure, PubChem assay and literature-proxy lineage keys use different provenance fields and are therefore reported separately (B–C). The complete machine-selected audit universe is not claimed to have been manually reviewed; panel D reports the fixed, source- and tier-stratified validation sample. A literature key is a conservative experiment proxy and may still group or split experiments imperfectly when source metadata are sparse.\n",encoding='utf-8')
 return paths,cov

def m8(a,dock):
 d=dock.copy(); d['old']=pd.to_numeric(d.ranking_score_v2,errors='coerce').fillna(0); d['new']=pd.to_numeric(d.ranking_score_v63,errors='coerce').fillna(0); d['delta']=d['new']-d['old']; d['bonus']=pd.to_numeric(d.lineage_bonus_v63,errors='coerce').fillna(0); d['penalty']=pd.to_numeric(d.conflict_penalty_v63,errors='coerce').fillna(0)
 fig=plt.figure(figsize=(11.6,7.6)); gs=fig.add_gridspec(2,2,hspace=.44,wspace=.32); aa,ab,ac,ad=[fig.add_subplot(gs[i]) for i in range(4)]
 pan(aa,'A','Existing shortlist membership is preserved')
 tier=d.docking_tier.value_counts().sort_index(); aa.bar(range(len(tier)),tier.values,color=C['blue']); aa.set_xticks(range(len(tier)),tier.index,rotation=25,ha='right'); aa.set_ylabel('Candidate pairs'); aa.grid(axis='y',color=C['light'],lw=.6)
 for i,v in enumerate(tier.values): aa.text(i,v,f'{v:,}',ha='center',va='bottom',fontsize=6.8)
 pan(ab,'B','V6.3 score change from lineage and conflict terms')
 bins=np.arange(np.floor(d.delta.min())-.5,np.ceil(d.delta.max())+1.5,1); ab.hist(d.delta,bins=bins,color=C['teal'],edgecolor='white',lw=.3); ab.axvline(0,color=C['dark'],lw=.8); ab.set_xlabel('V6.3 score − V6.2 score'); ab.set_ylabel('Candidate pairs'); ab.set_yscale('log'); ab.grid(axis='y',which='both',color=C['light'],lw=.6)
 pan(ac,'C','Lineage bonus and conflict penalty')
 bonus=d.bonus.value_counts().sort_index(); ac.plot(bonus.index,bonus.values,marker='o',ms=3,color=C['purple'],label='Lineage bonus'); pen=d.penalty.value_counts().sort_index(); ac.scatter(pen.index,pen.values,color=C['red'],s=24,label='Conflict penalty'); ac.set_yscale('log'); ac.set_xlabel('Score term'); ac.set_ylabel('Candidate pairs (log scale)'); ac.grid(color=C['light'],lw=.6); ac.legend(frameon=False,fontsize=7)
 pan(ad,'D','Docking candidates by V6.3 membrane role')
 role=d.membrane_role_primary_v63.replace('',np.nan).fillna('unclassified').value_counts().sort_values().tail(10); ad.barh(range(len(role)),role.values,color=[C['orange'] if x=='family_defined_membrane_role_unresolved' else C['blue'] for x in role.index]); ad.set_yticks(range(len(role)),[x.replace('_',' ') for x in role.index]); ad.set_xlabel('Candidate pairs'); ad.grid(axis='x',color=C['light'],lw=.6)
 for i,v in enumerate(role.values): ad.text(v,i,f' {v:,}',va='center',fontsize=6.5)
 fig.suptitle('M8 — Lineage-aware reranking of the frozen V6.2 docking shortlist',x=.06,ha='left',fontsize=12,fontweight='bold',color=C['navy']); fig.text(.06,.93,'No pair was newly admitted or removed: V6.3 adds capped lineage bonuses and an explicit positive–negative conflict penalty to the existing score.',ha='left',fontsize=7.5,color=C['dark']); fig.subplots_adjust(left=.14,right=.98,top=.86,bottom=.10)
 paths=save(fig,a.figure_dir,'M8_docking_lineage_rerank_v63')
 (a.figure_dir/'CAPTION_M8_DOCKING_LINEAGE_RERANK.md').write_text("""# M8 caption\n\n**M8. Lineage-aware reranking of the frozen V6.2 docking shortlist.** Candidate membership and docking tier are preserved from V6.2 (A). The V6.3 score adds a capped bonus for distinct structure, literature-proxy, PubChem-assay and evidence-modality lineages and subtracts a fixed penalty for a positive–negative evidence conflict (B–C). The five-axis membrane-role annotation is joined only for interpretation (D); it does not create new docking candidates. This is a prioritization score, not a predicted binding affinity.\n",encoding='utf-8')
 return paths,{"rows":len(d),"delta_min":float(d.delta.min()),"delta_median":float(d.delta.median()),"delta_max":float(d.delta.max())}

def main():
 a=ar(); setup(); a.figure_dir.mkdir(parents=True,exist_ok=True); a.qa_dir.mkdir(parents=True,exist_ok=True)
 pair=pd.read_csv(a.pair,sep='\t',compression='gzip',dtype=str,keep_default_na=False,usecols=['legacy_source_count_semantics_audit_v63','independent_structure_count_v63','independent_pubchem_assay_count_v63','independent_literature_experiment_count_v63','independent_evidence_modality_count_v63'])
 sample=pd.read_csv(a.review_sample,sep='\t',dtype=str,keep_default_na=False); dock=pd.read_csv(a.docking,sep='\t',compression='gzip',dtype=str,keep_default_na=False); val=json.loads(a.validation.read_text(encoding='utf-8'))
 p7,c7=m7(a,pair,sample,val); p8,c8=m8(a,dock); out={"status":"PASS","M7":p7,"M8":p8,"counts":{"M7_pair_rows":len(pair),"M8":c8,"M7_coverage":c7}}; (a.qa_dir/'V63_M7_M8_FIGURE_VALIDATION.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8'); print(json.dumps(out,ensure_ascii=False,indent=2))
if __name__=='__main__': main()
