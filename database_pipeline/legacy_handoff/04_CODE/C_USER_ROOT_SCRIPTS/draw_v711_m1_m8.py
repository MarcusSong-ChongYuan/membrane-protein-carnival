import json,math,textwrap
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch
import seaborn as sns

ROOT=Path(r'D:\finale\21_MemPro_V7.1.1_figures_20260814');ST=ROOT/'02_statistics';FIG=ROOT/'03_figures';CAP=ROOT/'04_captions'
S=json.loads((ST/'V711_FIGURE_STATISTICS.json').read_text(encoding='utf-8'))
for x in ['png','svg','pdf']:(FIG/x).mkdir(parents=True,exist_ok=True)
CAP.mkdir(parents=True,exist_ok=True)
sns.set_theme(style='whitegrid',context='paper')
plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.titlesize':11,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,'figure.titlesize':16,'text.color':'#30343b','axes.labelcolor':'#30343b','axes.titlecolor':'#30343b','axes.edgecolor':'#9aa1a9','grid.color':'#d9dde2','grid.linewidth':0.6,'legend.frameon':False,'pdf.fonttype':42,'svg.fonttype':'none'})
PAIRED=sns.color_palette('Paired',12);ROCKET=sns.color_palette('rocket',8);MUTED=sns.color_palette('muted',10)
def fmt(n):
 n=float(n)
 return f'{n/1e6:.2f}M' if n>=1e6 else (f'{n/1e3:.1f}k' if n>=1e3 else f'{int(n)}')
def title(fig,t,sub):
 fig.suptitle(t,x=.02,y=.99,ha='left',fontweight='bold');fig.text(.02,.955,sub,ha='left',va='top',fontsize=9,color='#59616b')
def panel(ax,label,t):ax.set_title(f'{label}  {t}',loc='left',fontweight='bold',pad=8)
def footer(fig,text):fig.text(.02,.01,'Source: '+text,fontsize=7,color='#6b7280',ha='left')
def clean(ax):ax.spines[['top','right']].set_visible(False)
def save(fig,name):
 for ext in ['png','svg','pdf']:
  fig.savefig(FIG/ext/f'{name}.{ext}',dpi=320 if ext=='png' else None,bbox_inches='tight',facecolor='white')
 plt.close(fig)
def hbar(ax,data,n=10,color=None,xlabel='Count'):
 arr=sorted(data.items(),key=lambda x:x[1],reverse=True)[:n][::-1];labels=[textwrap.shorten(str(x[0]),34,placeholder='…') for x in arr];vals=[x[1] for x in arr]
 ax.barh(range(len(vals)),vals,color=color or PAIRED[1],alpha=.78);ax.set_yticks(range(len(vals)),labels);ax.set_xlabel(xlabel);clean(ax)
 for i,v in enumerate(vals):ax.text(v,i,' '+fmt(v),va='center',fontsize=7)
def donut(ax,data,colors=None,center=''):
 vals=list(data.values());labs=list(data.keys());colors=colors or PAIRED[:len(vals)]
 ax.pie(vals,startangle=90,colors=colors,wedgeprops={'width':.38,'edgecolor':'white'},autopct=lambda p:f'{p:.1f}%' if p>=4 else '',pctdistance=.78,textprops={'fontsize':7})
 ax.text(0,0,center,ha='center',va='center',fontweight='bold');ax.legend(labs,loc='lower center',bbox_to_anchor=(.5,-.18),ncol=min(3,len(labs)),fontsize=7)
def lollipop(ax,data,n=12,color=None,xlabel='Proteins'):
 arr=sorted(data.items(),key=lambda x:x[1],reverse=True)[:n][::-1];y=np.arange(len(arr));v=np.array([x[1] for x in arr]);ax.hlines(y,0,v,color='#bac1c8',lw=1.2);ax.scatter(v,y,s=32,color=color or PAIRED[3],alpha=.8,zorder=3);ax.set_yticks(y,[textwrap.shorten(str(x[0]),30,placeholder='…') for x in arr]);ax.set_xlabel(xlabel);clean(ax)
 for yy,vv in zip(y,v):ax.text(vv,yy,' '+fmt(vv),va='center',fontsize=7)
def add_box(ax,xy,w,h,text,fc,ec='#8b949e'):
 r=patches.FancyBboxPatch(xy,w,h,boxstyle='round,pad=.015,rounding_size=.015',facecolor=fc,edgecolor=ec,lw=.8);ax.add_patch(r);ax.text(xy[0]+w/2,xy[1]+h/2,text,ha='center',va='center',fontsize=8);return r
def ribbon(ax,x0,y0,x1,y1,width,color,alpha=.35):
 verts=[(x0,y0-width/2),(x0+(x1-x0)*.45,y0-width/2),(x0+(x1-x0)*.55,y1-width/2),(x1,y1-width/2),(x1,y1+width/2),(x0+(x1-x0)*.55,y1+width/2),(x0+(x1-x0)*.45,y0+width/2),(x0,y0+width/2),(x0,y0-width/2)]
 codes=[MPath.MOVETO,MPath.CURVE4,MPath.CURVE4,MPath.CURVE4,MPath.LINETO,MPath.CURVE4,MPath.CURVE4,MPath.CURVE4,MPath.CLOSEPOLY];ax.add_patch(PathPatch(MPath(verts,codes),facecolor=color,edgecolor='none',alpha=alpha))

# M1
fig=plt.figure(figsize=(16,10));gs=fig.add_gridspec(2,2,height_ratios=[1.05,1],hspace=.35,wspace=.28);title(fig,'M1 | MemPro V7.1.1 architecture and scale','From heterogeneous records to explicit protein, compound, disease, complex and evidence entities')
ax=fig.add_subplot(gs[0,:]);panel(ax,'A','Integration workflow');ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis('off')
groups=[(.02,['Protein sources','UniProt · HPA · HTP','Membranome · OPM/PDBTM']),(.27,['Identity harmonization','UniProt · InChIKey','MONDO · source IDs']),(.52,['Evidence controls','Lineage · deduplication','positive/contextual negative']),(.77,['V7.1.1 release','24 entity/relation tables','manifest + semantic QA'])]
for i,(x,txt) in enumerate(groups):
 add_box(ax,(x,.35),.19,.35,'\n'.join(txt),PAIRED[(i*2)+1]);
 if i<3:ax.annotate('',xy=(groups[i+1][0]-.015,.525),xytext=(x+.195,.525),arrowprops={'arrowstyle':'->','lw':1.2,'color':'#606872'})
ax.text(.02,.16,'Raw records',fontsize=8,color='#6b7280');ax.text(.27,.16,'Canonical parent/form + protein-form granularity',fontsize=8,color='#6b7280');ax.text(.77,.16,'Frozen, traceable and reproducible',fontsize=8,color='#6b7280')
ax=fig.add_subplot(gs[1,0]);panel(ax,'B','Release entity scale (log10)');entities=S['release']['entities'];arr=sorted(entities.items(),key=lambda x:x[1]);y=np.arange(len(arr));v=np.array([x[1] for x in arr]);ax.barh(y,np.log10(v),color=PAIRED[1],alpha=.78);ax.set_yticks(y,[x[0] for x in arr]);ax.set_xlabel('log10(count)');clean(ax)
for yy,vv in zip(y,v):ax.text(math.log10(vv),yy,' '+fmt(vv),va='center',fontsize=7)
ax=fig.add_subplot(gs[1,1]);panel(ax,'C','Positive evidence contribution by database');hbar(ax,S['evidence']['source'],n=12,color=PAIRED[3],xlabel='Evidence records')
footer(fig,'MemPro V7.1.1 source registry; positive_interaction_evidence_v71; release manifest. Statistical units are labeled per panel.');save(fig,'M1_V711_architecture_scale_sources')

# M2
fig,axs=plt.subplots(2,2,figsize=(16,11),constrained_layout=True);title(fig,'M2 | Human membrane-protein landscape','Membrane mechanism, evidence strength and five-axis functional annotation')
panel(axs[0,0],'A','Membrane classes');donut(axs[0,0],S['protein']['class'],center='n=7,800')
panel(axs[0,1],'B','Membrane evidence levels');hbar(axs[0,1],S['protein']['evidence'],n=8,color=PAIRED[5],xlabel='Unique proteins')
panel(axs[1,0],'C','Five-axis annotation coverage');cov=S['protein']['axis_coverage'];lollipop(axs[1,0],cov,n=8,color=PAIRED[7],xlabel='Unique proteins with ≥1 annotation');axs[1,0].axvline(7800,color='#9aa1a9',ls='--',lw=.8)
panel(axs[1,1],'D','Top membrane roles');roles=dict(S['protein']['axis_top'].get('membrane_role',[]));lollipop(axs[1,1],roles,n=12,color=PAIRED[9],xlabel='Annotation assertions')
footer(fig,'UniProtKB 2026_02; HPA; HTP; Membranome; OPM/PDBTM; GO 2026-06-15; InterPro/Pfam frozen 2026-08-03; Reactome 97; specialist classifications.');save(fig,'M2_V711_membrane_proteome_five_axis')

# M3
fig,axs=plt.subplots(2,2,figsize=(16,11),constrained_layout=True);title(fig,'M3 | Expression and localization atlas','RNA abundance, IHC protein staining and tissue proteomics are separated; missing values are excluded from measured denominators')
panel(axs[0,0],'A','Expression records by corrected layer');hbar(axs[0,0],S['expression']['layer'],n=8,color=PAIRED[1],xlabel='Measurement records')
panel(axs[0,1],'B','Detection composition within each layer');layers=list(S['expression']['detection']);cats=['DETECTED','NOT_DETECTED','NOT_APPLICABLE'];bottom=np.zeros(len(layers))
for i,c in enumerate(cats):
 vals=np.array([S['expression']['detection'][x].get(c,0) for x in layers]);den=np.array([sum(S['expression']['detection'][x].values()) for x in layers]);pct=np.divide(vals,den,where=den>0)*100;axs[0,1].bar(np.arange(len(layers)),pct,bottom=bottom,label=c,color=PAIRED[i*2+1],alpha=.8);bottom+=pct
axs[0,1].set_xticks(range(len(layers)),[x.replace('_','\n') for x in layers]);axs[0,1].set_ylabel('Records (%)');axs[0,1].legend(ncol=3,fontsize=7,loc='upper center');clean(axs[0,1])
panel(axs[1,0],'C','Detected-location breadth per protein');order=['tissue_RNA','cell_type_RNA','normal_tissue_IHC','tissue_protein_MS'];data=[S['expression']['breadth'].get(x,[]) for x in order];parts=axs[1,0].violinplot(data,showextrema=False,showmedians=False)
for i,b in enumerate(parts['bodies']):b.set_facecolor(PAIRED[i*2+1]);b.set_alpha(.35);b.set_edgecolor('none')
axs[1,0].boxplot(data,widths=.15,showfliers=False,medianprops={'color':'#30343b'},boxprops={'color':'#59616b'},whiskerprops={'color':'#59616b'},capprops={'color':'#59616b'});axs[1,0].set_xticks(range(1,5),['Tissue RNA','Cell-type RNA','IHC staining','Tissue MS']);axs[1,0].set_ylabel('Detected tissues/cell types per protein');clean(axs[1,0])
panel(axs[1,1],'D','Protein identity mapping status');donut(axs[1,1],S['expression']['mapping'],center='3.05M\nrecords')
footer(fig,'Human Protein Atlas frozen source records; V7.1.1 semantic correction. Denominator: mapped + measured only; 103,709 empty tissue-MS cells remain missing/not-measured.');save(fig,'M3_V711_expression_localization_corrected')

# M4
fig=plt.figure(figsize=(16,11));gs=fig.add_gridspec(2,2,hspace=.35,wspace=.35);title(fig,'M4 | Protein–disease ontology landscape','Canonical MONDO identity, exact source mappings and multi-label DO/Uberon anatomy')
ax=fig.add_subplot(gs[0,0]);panel(ax,'A','Disease evidence strength');donut(ax,S['disease']['evidence'],center='6,753\nrelations')
ax=fig.add_subplot(gs[0,1]);panel(ax,'B','Contributing disease sources');lollipop(ax,S['disease']['sources'],n=10,color=PAIRED[5],xlabel='Relations carrying source')
ax=fig.add_subplot(gs[1,:]);panel(ax,'C','Top anatomical systems × evidence level')
systems=sorted(S['disease']['system_evidence'].items(),key=lambda x:sum(x[1].values()),reverse=True)[:14];levels=['medium','high','very_high'];mat=np.array([[x[1].get(l,0) for l in levels] for x in systems],float);sns.heatmap(np.log1p(mat),ax=ax,cmap='rocket',annot=mat.astype(int),fmt='d',cbar_kws={'label':'log1p(relations)'},yticklabels=[x[0] for x in systems],xticklabels=[x.replace('_',' ').title() for x in levels],linewidths=.4,linecolor='white');ax.set_xlabel('Evidence level');ax.set_ylabel('Ontology-derived anatomical system')
footer(fig,'Open Targets 26.06; UniProtKB disease annotation; MONDO 2026-07-06; Disease Ontology; Uberon basic. Multi-label assignments, not keyword-only forced categories.');save(fig,'M4_V711_disease_ontology_multilabel')

# M5
fig=plt.figure(figsize=(16,12));gs=fig.add_gridspec(2,2,hspace=.42,wspace=.35);title(fig,'M5 | Chemical identity and biological-status space','Canonical structures, parent/form hierarchy and interpretable compound-status intersections')
ax=fig.add_subplot(gs[0,0]);panel(ax,'A','Biological-status UpSet (top intersections)');inter=S['compound']['bio_intersections'];arr=sorted(inter.items(),key=lambda x:x[1],reverse=True)[:10];sets=['is_approved_drug','is_clinical_candidate','is_endogenous_ligand','is_natural_product','is_chemical_probe'];x=np.arange(len(arr));ax.bar(x,[v for _,v in arr],color=PAIRED[1],alpha=.8);ax.set_yscale('log');ax.set_ylabel('Canonical compounds (log scale)');ax.set_xticks(x,['' for _ in x]);clean(ax)
inset=ax.inset_axes([0,-.45,1,.38]);inset.set_xlim(-.5,len(arr)-.5);inset.set_ylim(-.5,len(sets)-.5);inset.axis('off')
for i,s in enumerate(sets):inset.text(-.65,len(sets)-1-i,s.replace('is_','').replace('_',' '),ha='right',va='center',fontsize=7)
for j,(k,v) in enumerate(arr):
 active=set(k.split('&')) if k!='None' else set()
 ys=[]
 for i,s in enumerate(sets):
  yy=len(sets)-1-i;on=s in active;inset.scatter(j,yy,s=28,color='#30343b' if on else '#d5d9de');
  if on:ys.append(yy)
 if len(ys)>1:inset.plot([j,j],[min(ys),max(ys)],color='#30343b',lw=1)
ax=fig.add_subplot(gs[0,1]);panel(ax,'B','Top computed structural classes');lollipop(ax,S['compound']['structural_class'],n=12,color=PAIRED[3],xlabel='Canonical compounds')
ax=fig.add_subplot(gs[1,0]);panel(ax,'C','Descriptor distributions (deterministic QC sample)');keys=['molecular_weight','xlogp','tpsa','rotatable_bond_count'];vals=[S['compound']['descriptor_sample'][k] for k in keys];# robust 1–99% clipping for readable violins
vals=[list(np.clip(v,np.percentile(v,1),np.percentile(v,99))) for v in vals];parts=ax.violinplot(vals,showextrema=False)
for i,b in enumerate(parts['bodies']):b.set_facecolor(ROCKET[i+2]);b.set_alpha(.35);b.set_edgecolor('none')
ax.boxplot(vals,widths=.15,showfliers=False);ax.set_xticks(range(1,5),['MW (Da)','XlogP','TPSA (Å²)','Rotatable bonds']);ax.set_ylabel('Value (property-specific scale)');clean(ax)
ax=fig.add_subplot(gs[1,1]);panel(ax,'D','Chemical identity confidence');hbar(ax,S['compound']['identity'],n=10,color=PAIRED[7],xlabel='Canonical compounds')
footer(fig,'Standardized compound master integrating ChEMBL, BindingDB, PubChem, BRENDA and additional mapped identifiers. UpSet memberships are biological statuses, not approval totals.');save(fig,'M5_V711_compound_upset_structure_properties')

# M6
fig=plt.figure(figsize=(16,11));gs=fig.add_gridspec(2,2,hspace=.38,wspace=.35);title(fig,'M6 | Binding evidence, sites and coordinate readiness','Interaction inclusion is independent of residue or coordinate availability')
ax=fig.add_subplot(gs[0,0]);panel(ax,'A','Residue → coordinate status (alluvial)');ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis('off')
left=[('Complete residues',49318),('Partial residues',3077),('Source-only residues',1886),('No residues reported',1329)];right=[('Complete coordinates',49097),('Partial coordinates',3077),('Source structure only',3436)];tot=55610
ly=[.82,.57,.34,.13];ry=[.76,.48,.22]
for i,(lab,v) in enumerate(left):add_box(ax,(.02,ly[i]-.055),.25,.11,f'{lab}\n{fmt(v)}',PAIRED[i*2+1])
for i,(lab,v) in enumerate(right):add_box(ax,(.73,ry[i]-.055),.25,.11,f'{lab}\n{fmt(v)}',PAIRED[i*2+2])
links=[(0,0,49097),(0,2,221),(1,1,3077),(2,2,1886),(3,2,1329)]
for a,b,v in links:ribbon(ax,.27,ly[a],.73,ry[b],max(.006,v/tot*.16),PAIRED[a*2+1],.35)
ax=fig.add_subplot(gs[0,1]);panel(ax,'B','Pair-level coordinate funnel');f=S['site']['pair_funnel'];labels=list(f);vals=list(f.values());y=np.arange(len(vals));ax.barh(y,vals,color=sns.color_palette('rocket',len(vals)),alpha=.82);ax.set_yticks(y,[x.replace('_',' ') for x in labels]);ax.invert_yaxis();ax.set_xlabel('Unique protein–compound pairs');clean(ax)
for i,v in enumerate(vals):ax.text(v,i,' '+fmt(v),va='center',fontsize=7)
ax=fig.add_subplot(gs[1,0]);panel(ax,'C','Chain assignment status');donut(ax,S['site']['chain'],center='55,610\nsites')
ax=fig.add_subplot(gs[1,1]);panel(ax,'D','Structure–query ligand relationship');hbar(ax,S['site']['ligand'],n=8,color=PAIRED[9],xlabel='Site assertions')
footer(fig,'Binding-site assertions linked to positive evidence; PDBe/PDB/SIFTS-derived coordinates. Only 663 exact HET-on-target-chain contact pockets are called exact co-crystals.');save(fig,'M6_V711_binding_sites_alluvial_docking')

# M7
fig,axs=plt.subplots(2,2,figsize=(16,11),constrained_layout=True);title(fig,'M7 | Provenance, negative evidence and release QA','Database contribution is separated from experiment and structure independence')
panel(axs[0,0],'A','Experiment-lineage status');hbar(axs[0,0],S['lineage']['experiment_status'],n=8,color=PAIRED[1],xlabel='Evidence records')
panel(axs[0,1],'B','Structure-lineage status');hbar(axs[0,1],S['lineage']['structure_status'],n=8,color=PAIRED[3],xlabel='Evidence records')
panel(axs[1,0],'C','Negative evidence disposition');neg=S['lineage']['negative'];axs[1,0].bar(neg.keys(),neg.values(),color=[PAIRED[5],PAIRED[7]],alpha=.8);axs[1,0].set_yscale('symlog',linthresh=1);axs[1,0].set_ylabel('Records (symlog)');clean(axs[1,0]);
for i,v in enumerate(neg.values()):axs[1,0].text(i,v,fmt(v),ha='center',va='bottom',fontsize=8)
panel(axs[1,1],'D','Blocking QA gates');gates={'Primary keys':0,'Foreign keys':0,'Semantic layer':0,'Manifest hashes':0};axs[1,1].scatter(range(len(gates)),[1]*len(gates),s=120,color=PAIRED[9],alpha=.8);axs[1,1].set_xticks(range(len(gates)),gates.keys(),rotation=20,ha='right');axs[1,1].set_yticks([1],['PASS']);axs[1,1].set_ylim(.7,1.3);axs[1,1].grid(axis='x',visible=False);clean(axs[1,1]);axs[1,1].text(1.5,.78,'0 blocking errors · 66 manifest files · 0 hash mismatches',ha='center',fontsize=8)
footer(fig,'Evidence lineage keys (publication/assay/target/compound/measurement and PDB/chain/ligand); V7.1.1 QA and manifest. Isolated inactive records are archived.');save(fig,'M7_V711_lineage_negative_QA')

# M8
fig=plt.figure(figsize=(16,10));gs=fig.add_gridspec(1,2,wspace=.32);title(fig,'M8 | Release summary and docking-ready subset','A compact handoff view for publication, website and structure-based follow-up')
ax=fig.add_subplot(gs[0,0]);panel(ax,'A','Released entity universe');entities=S['release']['entities'];arr=sorted(entities.items(),key=lambda x:x[1]);y=np.arange(len(arr));v=np.array([x[1] for x in arr]);ax.hlines(y,0,np.log10(v),color='#c7ccd2');ax.scatter(np.log10(v),y,s=48,color=PAIRED[1],alpha=.82);ax.set_yticks(y,[x[0] for x in arr]);ax.set_xlabel('log10(count)');clean(ax)
for yy,vv in zip(y,v):ax.text(math.log10(vv),yy,' '+fmt(vv),va='center',fontsize=7)
ax=fig.add_subplot(gs[0,1]);panel(ax,'B','From positive pairs to coordinate-ready pairs');f=S['release']['docking_funnel'];labels=['All positive pairs','With site assertion','Complete/partial coordinates','Complete coordinates','Coordinate docking eligible'];vals=[f['all_positive_pairs'],f['pairs_with_site_assertion'],f['pairs_with_complete_or_partial_coordinates'],f['pairs_with_complete_coordinates'],f['pairs_coordinate_docking_eligible']];maxv=max(vals)
for i,(lab,v) in enumerate(zip(labels,vals)):
 w=.82*(v/maxv)**.45;x=.5-w/2;y0=.88-i*.17;ax.add_patch(patches.FancyBboxPatch((x,y0),w,.1,boxstyle='round,pad=.005',facecolor=ROCKET[i+2],edgecolor='none',alpha=.75));ax.text(.5,y0+.05,f'{lab}: {fmt(v)}',ha='center',va='center',fontsize=8,color='white' if i>=2 else '#30343b')
ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis('off');ax.text(.5,.04,'Site absence does not remove a valid interaction; it only limits site-directed docking.',ha='center',fontsize=8,color='#59616b')
footer(fig,'MemPro V7.1.1 frozen release. Pair funnel is derived from positive evidence-linked site assertions, not from the separate 14k HPC benchmark package.');save(fig,'M8_V711_release_docking_readiness')

# Captions
caps={
'M1':'M1 summarizes the V7.1.1 integration architecture, logarithmic entity scale and positive-evidence contribution by source database. Counts refer to explicit database entities or evidence records, not independent experiments.',
'M2':'M2 describes the confirmed human membrane-protein layer by A/B/C membrane mechanism, E-level evidence and five independent annotation axes. Coverage is the number of unique proteins with at least one annotation on an axis.',
'M3':'M3 separates tissue RNA, single-cell RNA, IHC protein staining and tissue mass spectrometry. Expression breadth is computed only among mapped and measured records. Missing tissue-MS cells are not interpreted as zero or not detected.',
'M4':'M4 uses MONDO canonical disease identity and ontology-derived multi-label DO/Uberon anatomical systems. A disease may contribute to more than one anatomical system; heatmap cells show unique protein–disease relations by evidence level.',
'M5':'M5 shows biological-status intersections using an UpSet layout, computed structural classes, descriptor distributions from a deterministic QC sample and identity confidence. Biological-status intersections are not equivalent to approved-drug totals.',
'M6':'M6 separates source-reported residues from coordinate mapping. Interactions remain included without residues. Complete/partial coordinates and unique chain assignment determine coordinate-ready status; only strictly verified co-crystal contacts are labelled exact.',
'M7':'M7 separates contributing databases from experiment and structure lineage. Isolated inactive evidence is archived, while only same-canonical-pair contextual negatives enter the public conflict layer. All release QA gates pass.',
'M8':'M8 provides the final release scale and pair-level coordinate-readiness funnel. Coordinate readiness is a subset used for structural follow-up and does not define interaction validity.'}
for k,v in caps.items():(CAP/f'{k}_V711_caption.md').write_text(f'# {k} caption\n\n{v}\n',encoding='utf-8')
print(json.dumps({'figures':8,'formats':['png','svg','pdf'],'captions':8},ensure_ascii=False))
