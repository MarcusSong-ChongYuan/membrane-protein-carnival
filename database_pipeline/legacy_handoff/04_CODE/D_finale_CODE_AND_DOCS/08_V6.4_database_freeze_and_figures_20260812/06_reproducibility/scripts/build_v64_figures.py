from pathlib import Path
import pandas as pd,numpy as np,json,matplotlib as mpl,matplotlib.pyplot as plt,seaborn as sns,textwrap
ROOT=Path(r'D:\7.22\v64_candidate_working');S=ROOT/'stats';F=ROOT/'figures';
for x in ['png','svg','pdf','captions']: (F/x).mkdir(parents=True,exist_ok=True)
C={'navy':'#334E60','blue':'#6F8798','teal':'#4F8F88','orange':'#C98257','purple':'#80739E','gray':'#9AA1A6','red':'#B76E62','light':'#E8ECEF','text':'#2F3437'}
def setup():
 sns.set_theme(style='whitegrid',context='paper');mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Microsoft YaHei','DejaVu Sans'],'font.size':8.5,'axes.titlesize':10,'axes.titleweight':'bold','axes.labelsize':8.5,'xtick.labelsize':7.2,'ytick.labelsize':7.2,'legend.fontsize':7,'figure.titlesize':14,'text.color':C['text'],'axes.labelcolor':C['text'],'xtick.color':'#454A4E','ytick.color':'#454A4E','axes.edgecolor':'#9AA1A6','axes.linewidth':.6,'grid.color':'#D6D9DC','grid.linewidth':.45,'grid.alpha':.5,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
def fig(title):setup();f,a=plt.subplots(2,3,figsize=(14.2,9.0));f.subplots_adjust(left=.08,right=.98,top=.90,bottom=.09,wspace=.36,hspace=.52);f.suptitle(title,x=.01,y=.985,ha='left',color=C['navy'],fontweight='bold');return f,a.ravel()
def panel(ax,l,t,src):ax.text(-.08,1.08,l,transform=ax.transAxes,fontsize=11,fontweight='bold',va='top',color=C['navy']);ax.set_title(t,loc='left',pad=8);ax.text(0,-.19,'Source: '+src,transform=ax.transAxes,fontsize=6.2,color='#626A70',va='top',clip_on=False)
def clean(ax,axis='x'):ax.spines[['top','right']].set_visible(False);ax.grid(False);ax.grid(True,axis=axis,zorder=0)
def barh(ax,ser,color):ser=ser.sort_values();ax.barh(ser.index,ser.values,color=color,alpha=.86);mx=max(ser.values) if len(ser) else 1
# labels optional
def save(f,stem,cap):
 for e in ['png','svg','pdf']:f.savefig((F/e)/(f'{stem}.{e}'),dpi=600 if e=='png' else None,bbox_inches='tight')
 plt.close(f);(F/'captions'/f'{stem}.md').write_text(cap,encoding='utf-8')
def read(n):return pd.read_csv(S/n,sep='\t')
def top(df,key,val,n=10):return df.set_index(key)[val].head(n)
# M1
f,a=fig('M1  MemPro V6.4 candidate: traceable integration and controlled corrections')
panel(a[0],'A','Integrated entity scale','V6.4 candidate validation and frozen manifests');m=pd.Series({'Proteins':10997,'Isoforms':20211,'Complexes':19791,'Diseases':36153,'Compounds':2460533,'Binding evidence':3003306});barh(a[0],m,C['blue']);a[0].set_xscale('log');a[0].set_xlabel('Entities / records (log10)');clean(a[0])
panel(a[1],'B','Default release is distinct from all evidence','V6.4 release-scope accounting');x=pd.DataFrame({'Evidence rows':[991847,3003306-991847],'Unique pairs':[576038,1639146-576038]},index=['Default release','Review / non-default']);x.plot(kind='bar',ax=a[1],color=[C['teal'],C['purple']],alpha=.85);a[1].set_ylabel('Count');a[1].tick_params(axis='x',rotation=0);clean(a[1],'y')
panel(a[2],'C','PDB identity correction was targeted','RCSB current holdings snapshot 2026-08-12');q=pd.Series({'Unchanged':3002032,'Cleaned partial':1272,'All IDs invalid':2});barh(a[2],q,C['orange']);a[2].set_xscale('log');a[2].set_xlabel('Evidence rows (log10)');clean(a[2])
panel(a[3],'D','IHC identity correction preserves distinct biology','HPA 25.1; V6.4 deterministic identity');q=pd.Series({'Retained distinct rows':646118,'Exact duplicates removed':35351,'Legacy-ID conflicts rekeyed':9506});barh(a[3],q,C['teal']);a[3].set_xlabel('IHC records');clean(a[3])
panel(a[4],'E','Primary and foreign-key gates','V6.4 relational-integrity reports');a[4].axis('off');gates=[('Protein / isoform / gene PK','PASS'),('Complex / disease PK','PASS'),('Evidence / site PK','PASS'),('Protein-compound FK','PASS'),('Disease FK','PASS'),('HPA FK','PASS')];
for i,(t,s) in enumerate(gates):y=.9-i*.14;a[4].text(.03,y,t,fontsize=8);a[4].text(.92,y,s,ha='right',fontweight='bold',color=C['teal'])
panel(a[5],'F','Redistribution policy is source-aware','V6.4 source-license registry');a[5].axis('off');items=[('Public / attribution','UniProt, HPA, MONDO, DO, Uberon'),('Share-alike / conditional','ChEMBL, BindingDB-derived ChEMBL'),('Contributor-specific','PubChem BioAssay'),('Internal pending permission','BRENDA, PDBbind raw fields')];
for i,(t,s) in enumerate(items):y=.86-i*.2;a[5].text(.03,y,t,fontweight='bold',color=[C['teal'],C['blue'],C['orange'],C['red']][i]);a[5].text(.03,y-.08,s,fontsize=7)
save(f,'M1_V64_architecture_QA','V6.4 candidate integration, release-scope separation, controlled PDB/IHC corrections, relational QA, and conservative redistribution status. Counts are derived from frozen candidate tables; Docking outputs are not required.')
# M2
f,a=fig('M2  Human membrane-protein universe and five-axis classification')
for ax,l,t,n,k,v,c in [(a[0],'A','Membrane classes','protein_membrane_class_counts_v64.tsv','membrane_class','protein_count',C['blue']),(a[1],'B','Membrane-evidence levels','protein_evidence_level_counts_v64.tsv','evidence_level','protein_count',C['teal']),(a[2],'C','Molecular-function axis','protein_molecular_function_counts_v64.tsv','molecular_function','protein_count',C['orange']),(a[3],'D','Membrane-role axis','protein_membrane_role_counts_v64.tsv','membrane_role','protein_count',C['purple']),(a[4],'E','Most frequent structural families','protein_structural_family_counts_v64.tsv','structural_family','protein_count',C['blue']),(a[5],'F','Structure coverage','protein_structure_coverage_counts_v64.tsv','structure_status','protein_count',C['teal'])]:
 panel(ax,l,t,'UniProtKB 2026_02; GO; InterPro/Pfam; Reactome 97; OPM/PDBTM; AlphaFoldDB');d=read(n);ser=d.set_index(k)[v].head(10);barh(ax,ser,c);ax.set_xlabel('Unique canonical proteins');clean(ax)
save(f,'M2_V64_protein_landscape','Composition of 10,997 canonical membrane proteins. Five-axis annotations remain multi-dimensional: structural family, molecular function, biological process, membrane role and specialist classification.')
# M3
f,a=fig('M3  HPA expression and subcellular-localization coverage')
panel(a[0],'A','Protein-to-HPA mapping','HPA 25.1 and exact Ensembl mapping');d=read('protein_hpa_mapping_counts_v64.tsv');barh(a[0],d.set_index('hpa_mapping_status').protein_count,C['blue']);a[0].set_xlabel('Proteins');clean(a[0])
panel(a[1],'B','Expression breadth by measurement layer','HPA 25.1; RNA nTPM/nCPM and IHC ordinal calls');d=read('expression_breadth_summary_v64.tsv');xx=np.arange(len(d));a[1].errorbar(xx,d['median'],yerr=[d['median']-d['p25'],d['p75']-d['median']],fmt='o',color=C['teal'],capsize=4);a[1].set_xticks(xx,d['layer'],rotation=20,ha='right');a[1].set_ylabel('Detected tissues / cell types');clean(a[1],'y')
panel(a[2],'C','IHC detection calls after identity correction','HPA normal-tissue IHC 25.1');d=read('ihc_detection_counts_v64.tsv');barh(a[2],d.set_index('detection_status').records,C['orange']);a[2].set_xlabel('IHC records');clean(a[2])
panel(a[3],'D','IHC reliability','HPA 25.1');d=read('ihc_reliability_counts_v64.tsv');barh(a[3],d.set_index('reliability').records,C['purple']);a[3].set_xlabel('IHC records');clean(a[3])
panel(a[4],'E','Most frequent subcellular locations','HPA 25.1 localization');d=read('subcellular_location_counts_v64.tsv');barh(a[4],d.set_index('location').protein_count.head(12),C['blue']);a[4].set_xlabel('Protein annotations');clean(a[4])
panel(a[5],'F','Missingness is displayed, not interpreted as absence','HPA 25.1');a[5].axis('off');notes=['Mapped proteins: 10,664','No HPA Ensembl match: 333','RNA zero-detection and unmapped/missing are kept separate','IHC values are ordinal staining categories, not RNA quantities','Cell-type RNA uses nCPM; tissue RNA uses nTPM'];
for i,t in enumerate(notes):a[5].text(.04,.88-i*.16,'• '+t,fontsize=8)
save(f,'M3_V64_expression_localization','HPA expression and localization coverage. RNA and IHC layers are not combined into a single quantitative scale; mapping missingness is separated from biological non-detection.')
# M4
f,a=fig('M4  Canonical protein–disease relations with ontology-derived multi-label classification')
panel(a[0],'A','Disease evidence levels','Open Targets 26.06 and UniProtKB');d=read('disease_evidence_level_counts_v64.tsv');barh(a[0],d.set_index('evidence_level').protein_disease_pairs,C['blue']);a[0].set_xlabel('Canonical protein–disease pairs');clean(a[0])
panel(a[1],'B','Contributing disease sources','Open Targets and UniProtKB source evidence');d=read('disease_source_counts_v64.tsv');barh(a[1],d.set_index('source_database').protein_disease_pairs.head(10),C['teal']);a[1].set_xlabel('Source-supported pairs');clean(a[1])
panel(a[2],'C','Therapeutic-area memberships','Open Targets Platform 26.06, multi-label');d=read('disease_therapeutic_area_counts_v64.tsv');barh(a[2],d.set_index('therapeutic_area').relations.head(12),C['orange']);a[2].set_xlabel('Disease memberships');clean(a[2])
panel(a[3],'D','Anatomical-system memberships','MONDO/DO axioms to Uberon, multi-label');d=read('disease_anatomical_system_counts_v64.tsv');barh(a[3],d.set_index('anatomical_system').relations.head(12),C['purple']);a[3].set_xlabel('Disease–system memberships');clean(a[3])
panel(a[4],'E','Identity policy','MONDO 2026-07-06 official exact/equivalent mappings');a[4].axis('off');notes=['Direct MONDO IDs enter canonical layer','Official 1:1 exact/equivalent xrefs may merge','Broad, narrow, related and 1:N mappings never merge identity','Source IDs and names are permanently retained','Disease subtype and parent remain hierarchical entities'];
for i,t in enumerate(notes):a[4].text(.04,.88-i*.16,'• '+t,fontsize=8)
panel(a[5],'F','Release scale','V6.4 disease candidate');q=pd.Series({'Canonical disease entities':36153,'Protein–disease relations':9096,'Default relations':8822});barh(a[5],q,C['blue']);a[5].set_xlabel('Entities / relations');clean(a[5])
save(f,'M4_V64_disease_ontology','Canonical disease identity and ontology-derived multi-label therapeutic-area/anatomical-system classification. Counts are memberships and therefore may exceed unique disease counts.')
# M5
f,a=fig('M5  Canonical small-molecule identity, status and chemical properties')
panel(a[0],'A','Identity-confidence distribution','MemPro compound master v1.3 with V6.4 QA');d=read('compound_identity_confidence_counts_v64.tsv');barh(a[0],d.set_index('identity_confidence').compound_count,C['blue']);a[0].set_xlabel('Canonical compounds');clean(a[0])
panel(a[1],'B','Biological-status intersections','ChEMBL, DrugCentral, IUPHAR, ChEBI and curated status fields');d=read('compound_biostatus_intersections_v64.tsv');d=d[d.biological_status_intersection!='none'].head(12);barh(a[1],d.set_index('biological_status_intersection').compound_count,C['teal']);a[1].set_xlabel('Annotated compounds');clean(a[1])
panel(a[2],'C','Unified structural classes','RDKit deterministic rule hierarchy');d=read('compound_structural_class_counts_v64.tsv').head(10);barh(a[2],d.set_index('structural_class').compound_count,C['orange']);a[2].set_xlabel('Compounds');clean(a[2])
panel(a[3],'D','Physicochemical descriptor quartiles','RDKit/PubChem-normalized descriptors');d=read('compound_descriptor_summary_v64.tsv');xx=np.arange(len(d));a[3].errorbar(xx,d['median'],yerr=[d['median']-d['p25'],d['p75']-d['median']],fmt='o',color=C['purple'],capsize=4);a[3].set_xticks(xx,d.descriptor,rotation=20,ha='right');a[3].set_ylabel('Descriptor value (mixed units; see table)');clean(a[3],'y')
panel(a[4],'E','Parent/form identity model','InChIKey, canonical SMILES and exact form identifiers');a[4].axis('off');notes=['Canonical parent: standard structure identity','Form layer: salt, stereoisomer, charge/protonation representation','InChIKey supports fixed-length identity lookup','SMILES preserves an explicit graph representation','Names alone never trigger automatic identity merge'];
for i,t in enumerate(notes):a[4].text(.04,.88-i*.16,'• '+t,fontsize=8)
panel(a[5],'F','Release scale and uncertainty','V6.4 compound identity layer');q=pd.Series({'All canonical compounds':2460533,'High confidence':1117262,'Medium confidence':898802,'Low / unresolved structure':444469});barh(a[5],q,C['blue']);a[5].set_xscale('log');a[5].set_xlabel('Compounds (log10)');clean(a[5])
save(f,'M5_V64_compound_landscape','Canonical compounds are separated from exact forms. Biological-status counts are intersections among annotated compounds, not totals of approved drugs. Descriptor panels summarize identity-QC-passed values and are not drug-likeness conclusions.')
# M6
f,a=fig('M6  Protein–small-molecule evidence, sources and binding-site readiness')
panel(a[0],'A','Evidence tiers: all evidence layer','V6.4 cleaned binding evidence');d=read('binding_evidence_tier_counts_v64.tsv');barh(a[0],d.set_index('evidence_tier').evidence_rows,C['blue']);a[0].set_xlabel('Evidence rows');clean(a[0])
panel(a[1],'B','Contributing databases','V6.4 cleaned binding evidence');d=read('binding_source_counts_v64.tsv').head(12);barh(a[1],d.set_index('source_database').evidence_rows,C['teal']);a[1].set_xscale('log');a[1].set_xlabel('Evidence rows (log10)');clean(a[1])
panel(a[2],'C','PDB identity status','RCSB holdings 2026-08-12');d=read('binding_pdb_identity_counts_v64.tsv');barh(a[2],d.set_index('pdb_identity_status').evidence_rows,C['orange']);a[2].set_xscale('log');a[2].set_xlabel('Evidence rows (log10)');clean(a[2])
panel(a[3],'D','Binding-site types','V6.4 cleaned binding-site instances');d=read('binding_site_type_counts_v64.tsv').head(12);barh(a[3],d.set_index('site_type').site_rows,C['purple']);a[3].set_xlabel('Site instances');clean(a[3])
panel(a[4],'E','Binding-site sources','PDBe/PDBbind/BioLiP/BindingDB/UniProt and integrated sources');d=read('binding_site_source_counts_v64.tsv').head(12);barh(a[4],d.set_index('source_database').site_rows,C['blue']);a[4].set_xlabel('Site instances');clean(a[4])
panel(a[5],'F','Release and Docking interpretation','V6.4 release-scope accounting');a[5].axis('off');notes=['Default release: 991,847 evidence rows / 576,038 unique pairs','All evidence: 3,003,306 rows / 1,639,146 unique pairs','95,598 binding-site instances before site-level default filtering','Two structural rows lost all valid PDB IDs and left default structural layer','Docking scores and predicted pockets are not binding evidence'];
for i,t in enumerate(notes):a[5].text(.04,.88-i*.16,'• '+t,fontsize=8)
save(f,'M6_V64_binding_evidence_sites','Cleaned protein–small-molecule evidence and binding-site landscape. All-evidence and default-release scopes are shown separately; Docking readiness is not equivalent to experimental evidence strength.')
print('FIGURES_DONE')
