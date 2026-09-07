#!/usr/bin/env python3
import csv, json, re, time, urllib.request, urllib.error
from pathlib import Path
from collections import Counter

ROOT=Path(__file__).resolve().parent
MANIFEST=ROOT/'docking_manifest_D1_AE1_candidate.tsv'
CACHE=ROOT/'source_sifts_api'; CACHE.mkdir(exist_ok=True)
META=ROOT/'metadata'; META.mkdir(exist_ok=True)

AA={'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY','H':'HIS','I':'ILE','K':'LYS','L':'LEU','M':'MET','N':'ASN','P':'PRO','Q':'GLN','R':'ARG','S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR'}
def residues(s):
    out=[]
    for tok in re.split(r'[;\s]+',s or ''):
        tok=tok.split(':')[-1]
        m=re.search(r'([A-Za-z]{1,3})(-?\d+)([A-Za-z]?)',tok)
        if m:
            aa=m.group(1).upper(); out.append((AA.get(aa,aa),int(m.group(2)),m.group(3)))
    return out

with MANIFEST.open(encoding='utf-8',newline='') as f:
    reader=csv.DictReader(f,delimiter='\t'); fields=reader.fieldnames; tasks=list(reader)
pdbs=sorted({r['pdb_id'].lower() for r in tasks})
download=[]
for i,pdb in enumerate(pdbs,1):
    path=CACHE/f'{pdb}.json'
    if not path.exists() or path.stat().st_size<10:
        url=f'https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb}'
        ok=False; err=''
        for n in range(4):
            try:
                with urllib.request.urlopen(url,timeout=60) as r: path.write_bytes(r.read())
                ok=True; break
            except Exception as e: err=str(e); time.sleep(2**n)
        download.append([pdb,'PASS' if ok else 'FAIL',err])
    if i%100==0: print(f'sifts {i}/{len(pdbs)}',flush=True)
with (META/'sifts_download_status.tsv').open('w',encoding='utf-8',newline='') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['pdb_id','status','error']);w.writerows(download)

extra=['sifts_target_present','sifts_author_chains','sifts_chain_status','requested_residue_count','sifts_residue_range_covered_count','structure_qc_stage']
out=[]
for r in tasks:
    pdb=r['pdb_id'].lower(); up=r['target_uniprot']; req=residues(r.get('evidence_residues',''))
    maps=[]
    try:
        data=json.loads((CACHE/f'{pdb}.json').read_text())
        rec=data.get(pdb,{}).get('UniProt',{}).get(up,{})
        maps=rec.get('mappings',[])
    except Exception: pass
    chains=sorted({str(m.get('chain_id') or '') for m in maps if m.get('chain_id')})
    covered=0
    for aa,n,ins in req:
        if any((m.get('unp_start') or 10**12)<=n<=(m.get('unp_end') or -10**12) for m in maps): covered+=1
    x=dict(r); x.update({'sifts_target_present':'1' if maps else '0','sifts_author_chains':';'.join(chains),
        'sifts_chain_status':('UNIQUE' if len(chains)==1 else 'MULTIPLE_EQUIVALENT' if len(chains)>1 else 'NO_TARGET_MAPPING'),
        'requested_residue_count':str(len(req)),'sifts_residue_range_covered_count':str(covered),
        'structure_qc_stage':'SIFTS_PASS_RESIDUE_COORDINATES_PENDING' if maps and req and covered==len(req) else 'SIFTS_REVIEW'})
    out.append(x)
with (ROOT/'docking_manifest_D1_sifts_qc.tsv').open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields+extra,delimiter='\t');w.writeheader();w.writerows(out)
summary={'tasks':len(out),'pdbs':len(pdbs),'status':dict(Counter(r['structure_qc_stage'] for r in out)),'chain_status':dict(Counter(r['sifts_chain_status'] for r in out)),'download_attempts':len(download),'download_failures':sum(r[1]=='FAIL' for r in download)}
(META/'D1_SIFTS_QC_SUMMARY.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))
