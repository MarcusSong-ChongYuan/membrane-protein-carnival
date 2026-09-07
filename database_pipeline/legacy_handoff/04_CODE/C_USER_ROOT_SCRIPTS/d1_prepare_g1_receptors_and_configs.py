#!/usr/bin/env python3
import csv,hashlib,json,os,shutil,subprocess
from pathlib import Path
from collections import Counter
ROOT=Path(__file__).resolve().parent
SRC=ROOT/'docking_manifest_D1_G1_prefilter_pass.tsv'
PDB_CACHE=Path('/home/csong/docking/handoff/MemPro_Docking_Dscope_RCSBclean_14333_20260811/pdb_structures')
LIG_CACHE=Path('/home/csong/docking/handoff/MemPro_Docking_Dscope_RCSBclean_14333_20260811/ligands_pdbqt')
RAW=ROOT/'receptors_pdbqt_raw'; CLEAN=ROOT/'receptors_pdbqt_clean'; LIG=ROOT/'ligands_pdbqt'; CONF=ROOT/'vina_configs'; META=ROOT/'metadata'; LOG=ROOT/'logs'
for d in [RAW,CLEAN,LIG,CONF,META,LOG]:d.mkdir(exist_ok=True)
with SRC.open(encoding='utf-8',newline='') as f: rd=csv.DictReader(f,delimiter='\t');fields=rd.fieldnames;rows=list(rd)

def unsupported(line):
    if not line.startswith(('ATOM  ','HETATM')):return None
    typ=line[77:79].strip() if len(line)>=79 else ''
    if typ:return None
    try:return (float(line[30:38]),float(line[38:46]),float(line[46:54]),line[17:20].strip(),line[21:22].strip(),line[22:26].strip(),line[12:16].strip())
    except:return (None,None,None,'','','','')
def critical(a,r):
    if a[0] is None:return True
    return all(abs(a[j]-float(r[f'grid_center_{"xyz"[j]}']))<=float(r[f'grid_size_{"xyz"[j]}'])/2+5 for j in range(3))
def sha(path):
    h=hashlib.sha256();
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

# prepare each unique receptor with Open Babel, preserving unsupported atoms until task-aware QC
prep={}
for pdb in sorted({r['pdb_id'].lower() for r in rows}):
    src=PDB_CACHE/f'{pdb}.pdb'; out=RAW/f'{pdb}.pdbqt'; tmp=RAW/f'{pdb}.tmp.pdbqt'
    if not src.exists():prep[pdb]=('PDB_MISSING','');continue
    cp=subprocess.run(['timeout','180','obabel',str(src),'-O',str(tmp),'-h'],stdout=subprocess.DEVNULL,stderr=subprocess.PIPE,text=True)
    if cp.returncode!=0 or not tmp.exists():prep[pdb]=('OBABEL_FAILED',cp.stderr[-300:]);continue
    with tmp.open(encoding='utf-8',errors='ignore') as fi,out.open('w',encoding='utf-8',newline='\n') as fo:
        for line in fi:
            if line.startswith(('ATOM','HETATM','TER')):fo.write(line.rstrip('\r\n')+'\n')
    tmp.unlink(missing_ok=True);prep[pdb]=('PASS',sha(out))

extra=['receptor_prep_status','unsupported_atom_total','unsupported_atom_in_grid_or_5A','ligand_cache_status','production_config_status']
outrows=[]; receptor_safe={}
for r in rows:
    pdb=r['pdb_id'].lower(); tid=r['task_id']; status=prep.get(pdb,('PDB_MISSING',''))[0]; atoms=[]
    raw=RAW/f'{pdb}.pdbqt'
    if raw.exists():
        with raw.open(encoding='utf-8',errors='ignore') as f: atoms=[a for line in f if (a:=unsupported(line))]
    bad=[a for a in atoms if critical(a,r)]
    if status=='PASS' and bad:status='PREP_FAILED_UNSUPPORTED_CRITICAL_ATOM'
    ligsrc=LIG_CACHE/f"{r['compound_internal_id']}.pdbqt"; ligstatus='PASS' if ligsrc.exists() and ligsrc.stat().st_size else 'LIGAND_PDBQT_MISSING'
    if ligstatus=='PASS':
        dst=LIG/ligsrc.name
        if not dst.exists():os.link(ligsrc,dst)
    config='HOLD'
    if status=='PASS' and ligstatus=='PASS':
        safe=CLEAN/f'{pdb}.pdbqt'
        if pdb not in receptor_safe:
            with raw.open(encoding='utf-8',errors='ignore') as fi,safe.open('w',encoding='utf-8',newline='\n') as fo:
                for line in fi:
                    if unsupported(line) is None:fo.write(line.rstrip('\r\n')+'\n')
            receptor_safe[pdb]=sha(safe)
        txt=(f'receptor = receptors_pdbqt_clean/{pdb}.pdbqt\nligand = ligands_pdbqt/{r["compound_internal_id"]}.pdbqt\n'
             f'center_x = {r["grid_center_x"]}\ncenter_y = {r["grid_center_y"]}\ncenter_z = {r["grid_center_z"]}\n'
             f'size_x = {r["grid_size_x"]}\nsize_y = {r["grid_size_y"]}\nsize_z = {r["grid_size_z"]}\n'
             'thread = 8000\nnum_modes = 9\nenergy_range = 3\n'
             f'out = results/{tid}_out.pdbqt\n')
        (CONF/f'{tid}.conf').write_text(txt,encoding='ascii');config='READY_FOR_BENCHMARK_NOT_PRODUCTION'
    x=dict(r);x.update({'receptor_prep_status':status,'unsupported_atom_total':str(len(atoms)),'unsupported_atom_in_grid_or_5A':str(len(bad)),'ligand_cache_status':ligstatus,'production_config_status':config});outrows.append(x)
with (ROOT/'docking_manifest_D1_G1_receptor_qc.tsv').open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields+extra,delimiter='\t');w.writeheader();w.writerows(outrows)
summary={'input':len(rows),'configs':sum(r['production_config_status'].startswith('READY') for r in outrows),'receptor_status':dict(Counter(r['receptor_prep_status'] for r in outrows)),'ligand_status':dict(Counter(r['ligand_cache_status'] for r in outrows)),'unique_receptors_clean':len(receptor_safe)}
(META/'D1_G1_RECEPTOR_QC_SUMMARY.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(json.dumps(summary,indent=2))
