import csv,gzip,hashlib,json
from pathlib import Path
root=Path(r"D:\finale\09_V7_data_freeze_working_20260813");out=root/'05_v7_automatic_release_candidate'/'protein_master_V7.tsv.gz';old=root/'03_v7_relation_rebuild'/'canonical_protein_entity_V7.tsv.gz';jsondir=root/'02_protein_audit'/'uniprot_current_jsonl'
legacy={}
with gzip.open(old,'rt',encoding='utf-8',newline='') as f:
 for r in csv.DictReader(f,delimiter='\t'):legacy[r['canonical_uniprot_accession']]=r
xrefs={}
for p in sorted(jsondir.glob('batch_*.json')):
 d=json.loads(p.read_text(encoding='utf-8'))
 for r in d.get('results',[]):
  x={'HGNC':set(),'EnsemblGene':set(),'NCBIGene':set()}
  for q in r.get('uniProtKBCrossReferences',[]):
   db=q.get('database');qid=q.get('id','')
   if db=='HGNC':x['HGNC'].add(qid)
   elif db=='GeneID':x['NCBIGene'].add(qid)
   elif db=='Ensembl':
    for z in q.get('properties',[]):
     if z.get('key')=='GeneId':x['EnsemblGene'].add(z.get('value',''))
  xrefs[r['primaryAccession']]={k:';'.join(sorted(v)) for k,v in x.items()}
rows=[]
with gzip.open(out,'rt',encoding='utf-8',newline='') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  acc=r['canonical_uniprot_accession'];l=legacy.get(acc,{});x=xrefs.get(acc,{})
  r.update({'gene_entity_id':l.get('gene_entity_id',''),'hgnc_ids':x.get('HGNC',''),'ensembl_gene_ids':x.get('EnsemblGene',''),'ncbi_gene_ids':x.get('NCBIGene',''),'canonical_sequence':l.get('canonical_sequence',''),'canonical_sequence_sha256':l.get('canonical_sequence_sha256',''),'accession_history_status':'current_canonical_accession_validated_2026_08_13'})
  rows.append(r)
fields=list(rows[0])
with gzip.open(out,'wt',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
report={'rows':len(rows),'with_sequence':sum(bool(r['canonical_sequence']) for r in rows),'with_hgnc':sum(bool(r['hgnc_ids']) for r in rows),'with_ensembl_gene':sum(bool(r['ensembl_gene_ids']) for r in rows),'with_ncbi_gene':sum(bool(r['ncbi_gene_ids']) for r in rows),'unique_accessions':len({r['canonical_uniprot_accession'] for r in rows})}
(out.parent/'V7_T11_CANONICAL_PROTEIN_REPORT.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(report,ensure_ascii=False,indent=2))
