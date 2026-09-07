from pathlib import Path
import csv,gzip,json,re,hashlib,collections,sys,math
ROOT=Path(r'D:\finale\02_V6.3_candidate_20260803'); V62=Path(r'D:\finale\01_正式数据_V6.2'); OUT=Path(r'D:\7.22\v64_candidate_working')
files=[
 V62/'human_membrane_protein_master_v6_2.tsv',ROOT/'05_protein_annotation'/'human_membrane_protein_master_v6_3_candidate.tsv.gz',
 V62/'binding_evidence_master_v6_2.tsv.gz',V62/'binding_site_instances_v6_2.tsv.gz',V62/'small_molecule_master_v1_3.tsv',V62/'compound_form_hierarchy_v1_3.tsv',V62/'protein_gene_disease_relations_v6_2.tsv',
 ROOT/'02_identity'/'canonical_protein_entity_v0_1.tsv.gz',ROOT/'02_identity'/'protein_isoform_v0_2.tsv.gz',ROOT/'03_complex'/'complex_target_master_v0_2.tsv',ROOT/'03_complex'/'complex_target_components_v0_2.tsv',ROOT/'04_disease'/'disease_entity_master_v0_1.tsv',ROOT/'04_disease'/'protein_disease_relation_v2_1.tsv',
 V62/'protein_tissue_expression_v2.tsv.gz',V62/'protein_cell_type_expression_v2.tsv.gz',V62/'protein_tissue_cell_ihc_expression_v2.tsv.gz',V62/'protein_subcellular_localization_v2.tsv.gz']
UNIPROT=re.compile(r'^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?$');INCHIKEY=re.compile(r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$');PDB=re.compile(r'^[0-9][A-Za-z0-9]{3}$');CID=re.compile(r'^\d+$');CHEMBL=re.compile(r'^CHEMBL\d+$',re.I);MONDO=re.compile(r'^MONDO:\d{7}$'); UNITLIKE=re.compile(r'^\d+(?:\.\d+)?\s*(?:pM|nM|uM|µM|mM|M)$',re.I)
def op(p):return gzip.open(p,'rt',encoding='utf-8-sig',errors='replace',newline='') if p.suffix=='.gz' else open(p,'rt',encoding='utf-8-sig',errors='replace',newline='')
def audit(p):
 c=collections.Counter();samples=collections.defaultdict(list);pk_candidates=[];fk=collections.Counter();header=[]
 with op(p) as f:
  rd=csv.DictReader(f,delimiter='\t');header=[x.lstrip('\ufeff') for x in (rd.fieldnames or [])]; rd.fieldnames=header
  ids=[x for x in header if x.endswith('_id') or x.endswith('_ids') or 'accession' in x or x in ('pdb_ids','pubchem_cids','chembl_ids','standard_inchikey','exact_inchikey','canonical_disease_id','target_uniprot_id','canonical_uniprot_accession')]
  for row in rd:
   c['rows']+=1
   if len(row)!=len(header):c['column_count_anomaly']+=1
   for col in ids:
    v=(row.get(col) or '').strip()
    if not v:continue
    vals=re.split(r'[;,|]',v)
    for z in vals:
     z=z.strip()
     if not z:continue
     reason=None
     lc=col.lower()
     if 'uniprot' in lc or lc=='target_uniprot_id': reason=None if UNIPROT.fullmatch(z) else 'invalid_uniprot'
     elif lc in ('standard_inchikey','exact_inchikey'):reason=None if INCHIKEY.fullmatch(z) else 'invalid_inchikey'
     elif 'pdb' in lc: reason='unit_like_in_pdb' if UNITLIKE.fullmatch(z) else (None if PDB.fullmatch(z) else 'invalid_pdb_format')
     elif lc=='canonical_disease_id' and z.startswith('MONDO:'):reason=None if MONDO.fullmatch(z) else 'invalid_mondo'
     elif lc=='pubchem_cids':reason=None if CID.fullmatch(z) else 'invalid_pubchem_cid'
     elif lc=='chembl_ids':reason=None if CHEMBL.fullmatch(z) else 'invalid_chembl_id'
     if reason:
      c[reason]+=1
      if len(samples[reason])<20:samples[reason].append({'column':col,'value':z,'row_id':next((row.get(x,'') for x in header if x.endswith('_id') and row.get(x)), '')})
   for col,v0 in row.items():
    v=(v0 or '').strip()
    if not v:continue
    if (col.endswith('_id') or col.endswith('_ids')) and UNITLIKE.fullmatch(v):
     c['unit_like_in_id_field']+=1
     if len(samples['unit_like_in_id_field'])<20:samples['unit_like_in_id_field'].append({'column':col,'value':v})
    if col.endswith('_flag') or col.startswith('is_') or col in ('reviewed','website_default','default_release_inclusion'):
     if v.lower() not in {'0','1','true','false','yes','no','y','n','na','n/a','unknown','review',''}:c['noncanonical_boolean']+=1
 return {'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'columns':header,'counts':dict(c),'samples':dict(samples)}
res=[]
for i,p in enumerate(files,1):
 print(f'[{i}/{len(files)}] {p.name}',flush=True);res.append(audit(p))
(OUT/'qa'/'V64_FIELD_AUDIT_RAW.json').write_text(json.dumps(res,ensure_ascii=False,indent=2),encoding='utf-8')
with open(OUT/'qa'/'V64_FIELD_AUDIT_SUMMARY.tsv','w',encoding='utf-8',newline='') as f:
 fields=['file','rows','bytes','invalid_uniprot','invalid_inchikey','invalid_pdb_format','unit_like_in_pdb','invalid_mondo','invalid_pubchem_cid','invalid_chembl_id','unit_like_in_id_field','column_count_anomaly','noncanonical_boolean'];w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
 for x in res:
  c=x['counts'];w.writerow({'file':x['path'],'bytes':x['bytes'],**{k:c.get(k,0) for k in fields if k not in ('file','bytes')}})
print('DONE')
