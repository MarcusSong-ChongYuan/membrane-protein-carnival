#!/usr/bin/env python3
import csv,json,re,os
from pathlib import Path
from collections import Counter

ROOT=Path(__file__).resolve().parent
SRC=ROOT/'docking_manifest_D1_G1_EXPERIMENTAL_RESIDUES.tsv'
SIFTS=ROOT/'source_sifts_api'
PDB_ROOT=Path(os.environ.get('D1_PDB_CACHE','/home/csong/docking/handoff/MemPro_Docking_Dscope_RCSBclean_14333_20260811/pdb_structures'))
META=ROOT/'metadata'; META.mkdir(exist_ok=True)
AA={'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY','H':'HIS','I':'ILE','K':'LYS','L':'LEU','M':'MET','N':'ASN','P':'PRO','Q':'GLN','R':'ARG','S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR'}

def parse_req(s):
    out=[]
    for tok in re.split(r'[;\s]+',s or ''):
        tok=tok.strip()
        if not tok: continue
        chain_hint=''
        if ':' in tok:
            chain_hint,tok=tok.rsplit(':',1)
        m=re.search(r'([A-Za-z]{1,3})(-?\d+)([A-Za-z]?)',tok)
        if m:
            aa=m.group(1).upper(); out.append((AA.get(aa,aa),int(m.group(2)),m.group(3),chain_hint))
    return out

def atoms(path):
    idx={}
    with path.open(encoding='utf-8',errors='ignore') as f:
        for line in f:
            if not line.startswith(('ATOM  ','HETATM')) or len(line)<54: continue
            chain=line[21:22].strip(); rn=line[17:20].strip().upper(); ins=line[26:27].strip()
            try:num=int(line[22:26]); xyz=(float(line[30:38]),float(line[38:46]),float(line[46:54]))
            except ValueError:continue
            elem=(line[76:78].strip() if len(line)>=78 else '') or line[12:16].strip()[0]
            if elem.upper()!='H': idx.setdefault((chain,rn,num,ins),[]).append(xyz)
    return idx

def sifts_chains(pdb,up):
    try:
        d=json.loads((SIFTS/f'{pdb}.json').read_text()); rec=d.get(pdb,{}).get('UniProt',{}).get(up,{})
        return sorted({str(m.get('chain_id')) for m in rec.get('mappings',[]) if m.get('chain_id')})
    except Exception:return []

with SRC.open(encoding='utf-8',newline='') as f: rd=csv.DictReader(f,delimiter='\t'); fields=rd.fieldnames; rows=list(rd)
extra=['coordinate_qc_status','coordinate_qc_reason','actual_chain_used','candidate_sifts_chains','requested_residue_count_coord','matched_residue_count_coord','matched_residue_fraction_coord','grid_center_x','grid_center_y','grid_center_z','grid_size_x','grid_size_y','grid_size_z']
out=[]; cache={}
for i,r in enumerate(rows,1):
    pdb=r['pdb_id'].lower(); up=r['target_uniprot']; req=parse_req(r.get('evidence_residues','')); path=PDB_ROOT/f'{pdb}.pdb'; reasons=[]
    candidates=sifts_chains(pdb,up)
    if not path.exists() or path.stat().st_size==0: reasons.append('PDB_FILE_MISSING')
    if not candidates: reasons.append('NO_SIFTS_TARGET_CHAIN')
    if not req: reasons.append('NO_PARSEABLE_RESIDUES')
    best_chain=''; best=[]; tied=[]
    if not reasons:
        if pdb not in cache: cache[pdb]=atoms(path)
        idx=cache[pdb]; scores=[]
        for ch in candidates:
            matched=[]
            for aa,n,ins,hint in req:
                xyz=idx.get((ch,aa,n,ins),[])
                if not xyz and not ins:
                    xyz=[q for (c,rn,num,ic),vals in idx.items() if c==ch and rn==aa and num==n for q in vals]
                if xyz: matched.append((aa,n,ins,xyz))
            scores.append((len(matched),ch,matched))
        scores.sort(reverse=True,key=lambda x:x[0]); top=scores[0][0] if scores else 0; tied=[x for x in scores if x[0]==top]
        # Prefer manifest chain only when it is a SIFTS-confirmed top-scoring chain.
        manifest=(r.get('chain') or '').strip()
        chosen=next((x for x in tied if x[1]==manifest),None) if manifest else None
        if chosen is None and len(tied)==1: chosen=tied[0]
        if chosen is None: reasons.append('AMBIGUOUS_EQUIVALENT_CHAINS')
        else: best_chain=chosen[1]; best=chosen[2]
    frac=len(best)/len(req) if req else 0
    if not reasons and frac<0.8: reasons.append('RESIDUE_COORDINATE_COVERAGE_LT_80PCT')
    center=['','','']; size=['','','']
    if not reasons:
        xyz=[q for _,_,_,vals in best for q in vals]
        mins=[min(q[j] for q in xyz) for j in range(3)]; maxs=[max(q[j] for q in xyz) for j in range(3)]
        center=[f'{(mins[j]+maxs[j])/2:.3f}' for j in range(3)]
        dims=[max(20.0,maxs[j]-mins[j]+12.0) for j in range(3)]
        if any(v>30.0 for v in dims): reasons.append('BOX_AXIS_GT_30A'); center=['','','']
        else:size=[f'{v:.1f}' for v in dims]
    x=dict(r); x.update({'coordinate_qc_status':'PASS' if not reasons else 'HOLD','coordinate_qc_reason':';'.join(reasons),'actual_chain_used':best_chain,'candidate_sifts_chains':';'.join(candidates),'requested_residue_count_coord':str(len(req)),'matched_residue_count_coord':str(len(best)),'matched_residue_fraction_coord':f'{frac:.4f}','grid_center_x':center[0],'grid_center_y':center[1],'grid_center_z':center[2],'grid_size_x':size[0],'grid_size_y':size[1],'grid_size_z':size[2]})
    out.append(x)
    if i%200==0: print(f'coordinate {i}/{len(rows)}',flush=True)
for status in ['PASS','HOLD']:
    with (ROOT/f'docking_manifest_D1_G1_coordinate_{status.lower()}.tsv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields+extra,delimiter='\t');w.writeheader();w.writerows(r for r in out if r['coordinate_qc_status']==status)
summary={'input':len(out),'pass':sum(r['coordinate_qc_status']=='PASS' for r in out),'hold':sum(r['coordinate_qc_status']=='HOLD' for r in out),'hold_reasons':dict(Counter(x for r in out for x in r['coordinate_qc_reason'].split(';') if x)),'pdb_cache':str(PDB_ROOT)}
(META/'D1_G1_COORDINATE_QC_SUMMARY.json').write_text(json.dumps(summary,indent=2),encoding='utf-8'); print(json.dumps(summary,indent=2))
