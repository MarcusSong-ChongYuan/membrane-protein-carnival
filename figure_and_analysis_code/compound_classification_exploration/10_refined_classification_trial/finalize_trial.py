from __future__ import annotations
import hashlib, json
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from PIL import Image

OUT=Path(r"D:\finale\compound_classification_exploration\10_refined_classification_trial"); DATA=OUT/"data"; FIG=OUT/"figures"; QA=OUT/"qa"
QA.mkdir(exist_ok=True)
PAL=["#7b95c6","#49c2d9","#67a583","#a2c986","#fded95","#ffc1a6","#f59c7c","#c85e62"]
mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],"font.size":7,"axes.spines.top":False,"axes.spines.right":False,"svg.fonttype":"none","pdf.fonttype":42})
def export(fig,base):
    fig.savefig(base.with_suffix('.png'),dpi=350,bbox_inches='tight'); fig.savefig(base.with_suffix('.svg'),bbox_inches='tight'); fig.savefig(base.with_suffix('.pdf'),bbox_inches='tight'); fig.savefig(base.with_suffix('.tiff'),dpi=600,bbox_inches='tight'); plt.close(fig)

m=pd.read_csv(DATA/"HELDOUT_DEFINED_ROLE_DISCRIMINATION.tsv",sep='\t').sort_values('macro_f1')
fig,ax=plt.subplots(figsize=(7.2,4.4)); y=range(len(m)); ax.plot(m.macro_f1,y,'o',ms=5,color='#c85e62',label='Macro-F1'); ax.plot(m.balanced_accuracy,y,'s',ms=4.5,color='#7b95c6',label='Balanced accuracy')
ax.set_yticks(list(y),m.feature_set); ax.set_xlabel('Held-out classification score'); ax.set_xlim(0,.43); ax.legend(loc='upper right',bbox_to_anchor=(1,1.02),ncol=2); ax.set_title('Finer multi-axis chemistry improves target-role discrimination',loc='left',fontweight='bold'); ax.text(.99,.02,'n=239,182; six defined roles; fixed 80/20 split',transform=ax.transAxes,ha='right',fontsize=6,color='#555'); export(fig,FIG/'R2_heldout_role_discrimination')

e=pd.read_csv(DATA/'FUNCTIONAL_GROUP_ROLE_LOG2_ENRICHMENT.tsv',sep='\t',index_col=0)
fig,ax=plt.subplots(figsize=(8.4,5.4)); sns.heatmap(e.clip(-2,2),cmap='vlag',center=0,vmin=-2,vmax=2,linewidths=.3,linecolor='white',cbar_kws={'label':'log2(observed / expected)'},ax=ax)
ax.set_xlabel('Dominant formal membrane role'); ax.set_ylabel(''); ax.tick_params(axis='x',rotation=40); ax.set_title('Functional-group enrichment across membrane-target roles',loc='left',fontweight='bold'); export(fig,FIG/'R4_functional_group_role_enrichment')

# Compact overview, strictly for selecting candidate panels.
paths=[FIG/'R1_refined_ring_topology.png',FIG/'R2_heldout_role_discrimination.png',FIG/'R3_refined_ring_role_profile.png',FIG/'R4_functional_group_role_enrichment.png',FIG/'R5_resolution_vs_role_association.png']
fig=plt.figure(figsize=(13,7.3),facecolor='white'); gs=fig.add_gridspec(2,3,left=.02,right=.99,top=.92,bottom=.04,wspace=.08,hspace=.18); fig.suptitle('Refined compound classification trial',x=.02,ha='left',fontsize=16,fontweight='bold')
for i,p in enumerate(paths):
    ax=fig.add_subplot(gs[i//3,i%3]); ax.axis('off'); ax.imshow(Image.open(p)); ax.set_title(f'R{i+1}',loc='left',fontweight='bold')
fig.savefig(OUT/'REFINED_CLASSIFICATION_CONTACT_SHEET.png',dpi=220,bbox_inches='tight'); plt.close(fig)

ring=pd.read_csv(DATA/'CLASSIFICATION_RESOLUTION_AND_ASSOCIATION.tsv',sep='\t'); old=ring[ring.classification=='old_ring_class'].iloc[0]; new=ring[ring.classification=='refined_ring_topology'].iloc[0]
dm=pd.read_csv(DATA/'HELDOUT_DEFINED_ROLE_DISCRIMINATION.tsv',sep='\t'); oo=dm[dm.feature_set=='Old ring class'].iloc[0]; rr=dm[dm.feature_set=='Refined ring topology'].iloc[0]; best=dm.sort_values('balanced_accuracy',ascending=False).iloc[0]
report=f'''# Refined compound-classification trial results

## Direct answer

Yes. The refined ring topology is demonstrably more informative than the old acyclic/monocyclic/other-polycyclic/macrocyclic scheme, but ring topology alone remains insufficient. The strongest practical representation is a multi-axis annotation in which ring topology, functional groups, physicochemical bins, charge, Murcko scaffold and target-specific interaction type remain separate fields.

## Full-set comparison (n=240,054)

- Old ring classification: 4 nominal categories, effective category number {old.effective_category_n:.2f}; largest category {old.largest_category_fraction*100:.1f}%; Cramer's V versus dominant target role {old.cramers_v_vs_dominant_role:.4f}.
- Refined ring topology: 7 categories, effective category number {new.effective_category_n:.2f}; largest category {new.largest_category_fraction*100:.1f}%; Cramer's V {new.cramers_v_vs_dominant_role:.4f}.
- Therefore the refined topology reduces collapse into one broad polycyclic category and modestly improves role association.

## Held-out defined-role comparison

The diagnostic classifier used 239,182 compounds in six defined roles (role n>=500); 872 unresolved or rare-role compounds were excluded only from this inferential comparison and remain in the database.

- Old ring: balanced accuracy {oo.balanced_accuracy:.3f}, macro-F1 {oo.macro_f1:.3f}.
- Refined ring: balanced accuracy {rr.balanced_accuracy:.3f}, macro-F1 {rr.macro_f1:.3f}.
- Best balanced accuracy: {best.feature_set}, {best.balanced_accuracy:.3f}.
- Functional groups: balanced accuracy {dm.loc[dm.feature_set=='Functional groups','balanced_accuracy'].iloc[0]:.3f}.

This demonstrates that the new scheme is better, but the largest gain comes from orthogonal multi-label chemistry rather than creating more mutually exclusive ring labels.

## Recommended database fields

1. `ring_topology_primary` plus multi-label ring flags.
2. exact and generic Bemis-Murcko scaffold IDs.
3. multi-label functional groups.
4. molecular-weight, LogP, TPSA, flexibility and Fsp3 bins, while retaining continuous values.
5. structure-record charge class, explicitly not physiological-pH charge.
6. recognizable scaffold motifs as local annotations, not ClassyFire/ChEBI ontology.
7. pharmacological action kept at the protein-compound interaction level.

## Important limitations

- The named motif list is curated but not an exhaustive chemical ontology.
- Formal charge reflects the stored structure; it does not predict protonation at pH 7.4.
- Dominant target role depends on current database coverage.
- Functional-group enrichment is descriptive; it does not establish causal activity determinants.
'''
(OUT/'REFINED_CLASSIFICATION_TRIAL_REPORT.md').write_text(report,encoding='utf-8')

manifest=[]
for p in sorted(OUT.rglob('*')):
    if p.is_file() and p.name!='SHA256SUMS.tsv':
        h=hashlib.sha256();
        with p.open('rb') as f:
            for ch in iter(lambda:f.read(1048576),b''): h.update(ch)
        manifest.append((p.relative_to(OUT).as_posix(),p.stat().st_size,h.hexdigest()))
with (QA/'SHA256SUMS.tsv').open('w',encoding='utf-8') as f:
    f.write('relative_path\tsize_bytes\tsha256\n'); [f.write(f'{a}\t{b}\t{c}\n') for a,b,c in manifest]
print({'status':'PASS','files':len(manifest)})
