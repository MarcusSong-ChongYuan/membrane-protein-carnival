#!/usr/bin/env python3
import csv, json, urllib.request
from pathlib import Path
import gemmi

root = Path(__file__).resolve().parent
manifest = root / 'docking_manifest.tsv'
pdb_dir = root / 'pdb_structures'
cif_dir = root / 'source_mmcif'
meta = root / 'metadata'
pdb_dir.mkdir(exist_ok=True); cif_dir.mkdir(exist_ok=True); meta.mkdir(exist_ok=True)

needed = {}
with manifest.open(encoding='utf-8', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        pdb = (row.get('pdb_id') or '').strip().lower()
        chain = (row.get('chain') or '').strip()
        if pdb and chain:
            needed.setdefault(pdb, set()).add(chain)

rows = []
for pdb, chains in sorted(needed.items()):
    out = pdb_dir / f'{pdb}.pdb'
    valid_existing = out.exists() and out.stat().st_size > 0
    if valid_existing:
        try:
            st0 = gemmi.read_structure(str(out))
            valid_existing = bool(st0) and any(c.name in chains for c in st0[0])
        except Exception:
            valid_existing = False
    if valid_existing:
        continue

    cif = cif_dir / f'{pdb}.cif'
    try:
        if not cif.exists() or cif.stat().st_size == 0:
            urllib.request.urlretrieve(f'https://files.rcsb.org/download/{pdb}.cif', cif)
        st = gemmi.read_structure(str(cif))
        available = {c.name for c in st[0]}
        selected = sorted(set(chains) & available)
        if not selected:
            raise RuntimeError(f'requested chains {sorted(chains)} absent; available sample={sorted(available)[:20]}')
        new = gemmi.Structure()
        new.name = st.name
        new.cell = st.cell
        new.spacegroup_hm = st.spacegroup_hm
        model = gemmi.Model(st[0].name)
        for ch in st[0]:
            if ch.name in selected:
                model.add_chain(ch.clone())
        new.add_model(model)
        new.setup_entities()
        new.assign_label_seq_id()
        new.write_pdb(str(out))
        check = gemmi.read_structure(str(out))
        atom_count = sum(1 for model in check for chain in model for residue in chain for atom in residue)
        if atom_count == 0:
            raise RuntimeError('zero atoms after PDB export')
        rows.append({'pdb_id':pdb,'requested_chains':';'.join(sorted(chains)),
                     'exported_chains':';'.join(selected),'cif_bytes':cif.stat().st_size,
                     'pdb_bytes':out.stat().st_size,'atom_count':atom_count,'status':'PASS','error':''})
    except Exception as exc:
        if out.exists(): out.unlink()
        rows.append({'pdb_id':pdb,'requested_chains':';'.join(sorted(chains)),
                     'exported_chains':'','cif_bytes':cif.stat().st_size if cif.exists() else 0,
                     'pdb_bytes':0,'atom_count':0,'status':'FAIL','error':str(exc)})

fields=['pdb_id','requested_chains','exported_chains','cif_bytes','pdb_bytes','atom_count','status','error']
with (meta/'gemmi_cif_chain_extraction.tsv').open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(rows)
summary={'processed':len(rows),'pass':sum(r['status']=='PASS' for r in rows),
         'fail':sum(r['status']=='FAIL' for r in rows),'gemmi_version':gemmi.__version__}
(meta/'gemmi_cif_chain_extraction_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary))
if summary['fail']:
    raise SystemExit(2)
