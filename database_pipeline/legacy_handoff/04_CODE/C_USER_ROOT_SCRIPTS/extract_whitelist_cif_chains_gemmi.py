#!/usr/bin/env python3
"""Repair only explicitly whitelisted PDB IDs from RCSB mmCIF."""
import csv, json, urllib.request
from pathlib import Path
import gemmi

ROOT=Path(__file__).resolve().parent
ALLOW={'22tu','5lks','7mq8','7o9k','7u9q','8g60','8qyx','9ndp','9p7h','9pa7'}
chains={x:set() for x in ALLOW}
with (ROOT/'docking_manifest.tsv').open(encoding='utf-8',newline='') as f:
    for r in csv.DictReader(f,delimiter='\t'):
        p=(r.get('pdb_id') or '').lower().strip()
        if p in ALLOW and (r.get('chain') or '').strip(): chains[p].add(r['chain'].strip())

rows=[]; cifdir=ROOT/'source_mmcif_whitelist'; cifdir.mkdir(exist_ok=True)
for pdb in sorted(ALLOW):
    out=ROOT/'pdb_structures'/f'{pdb}.pdb'; cif=cifdir/f'{pdb}.cif'
    try:
        if not cif.exists(): urllib.request.urlretrieve(f'https://files.rcsb.org/download/{pdb}.cif',cif)
        st=gemmi.read_structure(str(cif)); available={c.name for c in st[0]}; selected=sorted(chains[pdb]&available)
        if not selected: raise RuntimeError(f'requested={sorted(chains[pdb])}; available={sorted(available)[:30]}')
        new=gemmi.Structure(); new.name=st.name; new.cell=st.cell; new.spacegroup_hm=st.spacegroup_hm
        model=gemmi.Model(st[0].name)
        for ch in st[0]:
            if ch.name in selected: model.add_chain(ch.clone())
        new.add_model(model); new.setup_entities(); new.assign_label_seq_id()
        tmp=out.with_suffix('.pdb.tmp'); new.write_pdb(str(tmp)); check=gemmi.read_structure(str(tmp))
        atoms=sum(1 for m in check for c in m for res in c for a in res)
        if atoms==0: raise RuntimeError('zero atoms')
        tmp.replace(out)
        rows.append([pdb,';'.join(sorted(chains[pdb])),';'.join(selected),cif.stat().st_size,out.stat().st_size,atoms,'PASS',''])
    except Exception as e:
        rows.append([pdb,';'.join(sorted(chains[pdb])),'',cif.stat().st_size if cif.exists() else 0,0,0,'FAIL',str(e)])

path=ROOT/'metadata'/'gemmi_whitelist_extraction.tsv'; path.parent.mkdir(exist_ok=True)
with path.open('w',encoding='utf-8',newline='') as f:
    w=csv.writer(f,delimiter='\t'); w.writerow(['pdb_id','requested_chains','exported_chains','cif_bytes','pdb_bytes','atom_count','status','error']); w.writerows(rows)
summary={'allowlist':len(ALLOW),'pass':sum(r[6]=='PASS' for r in rows),'fail':sum(r[6]=='FAIL' for r in rows),'gemmi':gemmi.__version__}
(ROOT/'metadata'/'gemmi_whitelist_extraction_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary)); raise SystemExit(0 if summary['fail']==0 else 2)
